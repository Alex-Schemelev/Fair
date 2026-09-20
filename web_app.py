#!/usr/bin/env python3
"""
Единый скрипт запуска веб‑приложения:
- отдаёт собранный React‑фронтенд
- предоставляет API /api/run, которое дергает In_silico_PCR.py
"""

import os
import tempfile
import subprocess
import json
import zipfile
import sys

# Ensure project root is on sys.path (portable Python / non-cwd launches).
_ROOT = os.path.dirname(os.path.abspath(__file__))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

import pandas as pd
from flask import Flask, request, jsonify, send_from_directory, send_file
import shutil
import uuid

from pcrlib.paths import (
    REF_FILENAME,
    PRIMERS_FILENAME,
    CONSERVATION_FILENAME,
    LEGACY_REF_NAMES,
    find_existing,
)
from pcrlib.iupac import mismatch_stats_vs_consensus
from pcrlib.align import simple_align
from pcrlib.conservation import (
    compute_conservation_for_primers,
    compute_consensus_qc_for_products,
)
from pcrlib.filters import apply_product_filters
from pcrlib.annotation import (
    rebuild_annotation_products,
    stable_product_id,
    extract_annotation_map,
)
from pcrlib.session_io import (
    load_dimers_matrix_from_run,
    safe_extract_zip,
    slim_session_payload,
    cleanup_old_runs,
)

# Значения по умолчанию — синхронизированы с In_silico_PCR.py
MIN_AMPLICON = 100
MAX_AMPLICON = 10000
MAX_DEGENERATE_VARIANTS = 50

MAX_TM_DIFF = 7.0
MAX_3PRIME_MATCHES = 3

ANNEALING_TEMP = 60.0
MAX_HUMAN_DELTA = 1500

MAX_MISMATCHES_TOTAL = 0
MAX_MISMATCHES_3PRIME = 0
THREE_PRIME_WINDOW = 5

FRONTEND_DIR = os.path.join(os.path.dirname(__file__), "frontend")
FRONTEND_DIST = os.path.join(FRONTEND_DIR, "dist")


def ensure_frontend_built():
    """
    Гарантирует наличие React‑интерфейса в frontend/dist/index.html.

    Канон — файлы в frontend/dist (их правим вручную). Встроенный HTML‑шаблон
    записывается только если index.html ещё нет (первый запуск / чистая копия).
    """
    index_html = os.path.join(FRONTEND_DIST, "index.html")
    os.makedirs(FRONTEND_DIST, exist_ok=True)
    if os.path.isfile(index_html) and os.path.getsize(index_html) > 200:
        return

    html = """<!doctype html>
<html lang="ru">
  <head>
    <meta charset="UTF-8" />
    <title>In silico PCR UI</title>
    <meta name="viewport" content="width=device-width, initial-scale=1.0" />
  </head>
  <body>
    <h1>Нет frontend/dist/index.html</h1>
    <p>Скопируйте интерфейс в <code>frontend/dist/</code> (index.html + annotation_map.js).</p>
  </body>
</html>
"""
    with open(index_html, "w", encoding="utf-8") as f:
        f.write(html)


app = Flask(
    __name__,
    static_folder=FRONTEND_DIST,
    static_url_path="/",
)

# Хранилище результатов последнего запуска (в памяти процесса)
RUNS_DIR = os.path.join(os.path.dirname(__file__), "runs")
os.makedirs(RUNS_DIR, exist_ok=True)
LATEST_RUN = {"id": None, "dir": None}


@app.post("/api/run")
def run_pcr():
    ref_file = request.files.get("reference")
    if ref_file is None:
        return jsonify({"error": "reference file is required"}), 400

    payload_raw = request.form.get("payload", "{}")
    payload = json.loads(payload_raw)

    primers = payload.get("primers", [])
    blast_enabled = bool(payload.get("blast_enabled", True))
    blast_human_enabled = bool(payload.get("blast_human_enabled", False))
    params = payload.get("params", {})
    filters = payload.get("filters", {})

    tmpdir = tempfile.mkdtemp(prefix="in_silico_pcr_")

    ref_path = os.path.join(tmpdir, REF_FILENAME)
    ref_file.save(ref_path)

    primers_csv_path = os.path.join(tmpdir, PRIMERS_FILENAME)

    primers_table_file = request.files.get("primers_table")
    if primers_table_file:
        # Парсим загруженную таблицу: первая колонка — имя, вторая — последовательность
        uploaded_path = os.path.join(tmpdir, "uploaded_primers_table.txt")
        primers_table_file.save(uploaded_path)
        try:
            df = pd.read_csv(uploaded_path, sep=None, engine="python", header=None)
        except Exception as e:
            shutil.rmtree(tmpdir, ignore_errors=True)
            return jsonify({"error": f"Не удалось прочитать таблицу праймеров: {e}"}), 400

        if df.shape[1] < 2:
            shutil.rmtree(tmpdir, ignore_errors=True)
            return jsonify({"error": "В таблице праймеров должно быть минимум две колонки (имя, последовательность)."}), 400

        df = df.iloc[:, :2]
        df.columns = ["Name", "Sequence"]
        df["Type"] = "auto"
        df.to_csv(primers_csv_path, sep=";", index=False)
        primers_dict = {str(r[0]): str(r[1]) for r in df[["Name", "Sequence"]].itertuples(index=False, name=None)}
    else:
        # Используем праймеры из формы (JSON), генерируем CSV в нужном формате
        with open(primers_csv_path, "w", encoding="utf-8") as f:
            f.write("Name;Sequence;Type\n")
            for p in primers:
                name = p.get("name", "")
                seq = p.get("sequence", "")
                if not name or not seq:
                    continue
                f.write(f"{name};{seq};auto\n")
        primers_dict = {p.get("name"): p.get("sequence") for p in primers if p.get("name") and p.get("sequence")}

    # опциональная таблица консервативности
    cons_file = request.files.get("conservation_table")
    cons_path = None
    if cons_file:
        cons_path = os.path.join(tmpdir, CONSERVATION_FILENAME)
        cons_file.save(cons_path)

    script_path = os.path.join(os.path.dirname(__file__), "In_silico_PCR.py")
    env = os.environ.copy()
    env["BLAST_ENABLED"] = "1" if blast_enabled else "0"
    env["BLAST_HUMAN_ENABLED"] = "1" if blast_human_enabled else "0"
    env["PRIMER_CSV"] = PRIMERS_FILENAME
    env["TARGET_GB"] = REF_FILENAME

    # Передаём числовые параметры в In_silico_PCR.py через переменные окружения
    def set_int_env(key):
        if key in params and params[key] is not None:
            try:
                env[key] = str(int(params[key]))
            except Exception:
                pass

    def set_float_env(key):
        if key in params and params[key] is not None:
            try:
                env[key] = str(float(params[key]))
            except Exception:
                pass

    set_int_env("MIN_AMPLICON")
    set_int_env("MAX_AMPLICON")
    set_int_env("MAX_DEGENERATE_VARIANTS")
    # MAX_TM_DIFF опционален: пустое значение = без фильтра (в пайплайн не передаём)
    if params.get("MAX_TM_DIFF") not in (None, ""):
        set_float_env("MAX_TM_DIFF")
    set_int_env("MAX_3PRIME_MATCHES")
    set_float_env("ANNEALING_TEMP")
    set_int_env("MAX_HUMAN_DELTA")
    set_int_env("MAX_MISMATCHES_TOTAL")
    set_int_env("MAX_MISMATCHES_3PRIME")
    set_int_env("THREE_PRIME_WINDOW")

    log_path = os.path.join(tmpdir, "pipeline.log")
    try:
        with open(log_path, "w", encoding="utf-8") as logf:
            logf.write("executable: " + str(sys.executable) + "\n")
            logf.write("script: " + str(script_path) + "\n")
            logf.write("BLAST_ENABLED=" + str(env.get("BLAST_ENABLED")) + "\n")
            logf.write("BLAST_HUMAN_ENABLED=" + str(env.get("BLAST_HUMAN_ENABLED")) + "\n\n")
            logf.flush()
            proc = subprocess.run(
                [sys.executable, script_path],
                cwd=tmpdir,
                check=False,
                env=env,
                stdout=logf,
                stderr=subprocess.STDOUT,
            )
        if proc.returncode != 0:
            tail = ""
            try:
                with open(log_path, "r", encoding="utf-8", errors="replace") as lf:
                    tail = lf.read()[-4000:]
            except Exception:
                pass
            shutil.rmtree(tmpdir, ignore_errors=True)
            return jsonify({
                "error": "In_silico_PCR.py failed (code " + str(proc.returncode) + ")",
                "pipeline_log_tail": tail,
            }), 500
    except Exception as e:
        shutil.rmtree(tmpdir, ignore_errors=True)
        return jsonify({"error": f"In_silico_PCR.py failed: {e}"}), 500

    # Копируем все выходные файлы в постоянную папку запуска, чтобы можно было скачивать
    run_id = uuid.uuid4().hex[:12]
    run_dir = os.path.join(RUNS_DIR, run_id)
    os.makedirs(run_dir, exist_ok=True)
    for fn in os.listdir(tmpdir):
        src = os.path.join(tmpdir, fn)
        if os.path.isfile(src):
            try:
                shutil.copy2(src, os.path.join(run_dir, fn))
            except Exception:
                pass
    LATEST_RUN["id"] = run_id
    LATEST_RUN["dir"] = run_dir

    results_csv = os.path.join(tmpdir, "PCR_results_primer3py.csv")
    if not os.path.exists(results_csv):
        shutil.rmtree(tmpdir, ignore_errors=True)
        return jsonify({"error": "PCR_results_primer3py.csv not found"}), 500

    # Файл может существовать, но быть пустым (нет найденных ампликонов)
    try:
        if os.path.getsize(results_csv) == 0:
            df = pd.DataFrame([])
        else:
            df = pd.read_csv(results_csv)
    except pd.errors.EmptyDataError:
        df = pd.DataFrame([])

    dimers = []
    if not df.empty:
        for _, row in df.iterrows():
            seq1 = str(row.get("Forward_sequence", ""))
            seq2 = str(row.get("Reverse_sequence", ""))
            if not seq1 or not seq2:
                continue
            offset, matches = simple_align(seq1, seq2)
            dimers.append(
                {
                    "primer1": row.get("Forward_primer", ""),
                    "primer2": row.get("Reverse_primer", ""),
                    "seq1": seq1,
                    "seq2": seq2,
                    "offset": int(offset),
                    "matches": matches,
                    "Tm_F": float(row.get("Tm_F", 0)),
                    "Tm_R": float(row.get("Tm_R", 0)),
                    "Heterodimer_TM": float(row.get("Heterodimer_TM", 0)),
                    "Heterodimer_dG": float(row.get("Heterodimer_dG", 0)),
                    "three_prime_matches": int(row.get("3prime_matches", 0)),
                }
            )

    # Читаем отчёт по посадкам праймеров
    bindings = []
    bindings_csv = os.path.join(tmpdir, "primer_binding_sites.csv")
    if os.path.exists(bindings_csv) and os.path.getsize(bindings_csv) > 0:
        try:
            df_bind = pd.read_csv(bindings_csv)
            bindings = df_bind.to_dict(orient="records")
        except Exception:
            bindings = []

    # Читаем матрицу димеров (гомо- и гетеродимеры)
    dimers_matrix = []
    dimers_csv = os.path.join(tmpdir, "primer_dimers_matrix.csv")
    if os.path.exists(dimers_csv) and os.path.getsize(dimers_csv) > 0:
        try:
            df_dim = pd.read_csv(dimers_csv)
            dimers_matrix = df_dim.to_dict(orient="records")
        except Exception:
            dimers_matrix = []

    # QC по гомодимерам
    homo_qc = []
    homo_qc_csv = os.path.join(tmpdir, "primer_homodimer_qc.csv")
    if os.path.exists(homo_qc_csv) and os.path.getsize(homo_qc_csv) > 0:
        try:
            df_hq = pd.read_csv(homo_qc_csv)
            homo_qc = df_hq.to_dict(orient="records")
        except Exception:
            homo_qc = []

    # BLAST-посадки праймеров на референс
    blast_bindings = []
    blast_csv = os.path.join(tmpdir, "primer_blast_ref_hits.csv")
    if os.path.exists(blast_csv) and os.path.getsize(blast_csv) > 0:
        try:
            df_blast = pd.read_csv(blast_csv)
            df_blast = df_blast.where(pd.notnull(df_blast), None)
            tmp_b = df_blast.to_dict(orient="records")
            cleaned_b = []
            for row in tmp_b:
                new_row = {}
                for k, v in row.items():
                    import math
                    if isinstance(v, float) and math.isnan(v):
                        new_row[k] = None
                    else:
                        new_row[k] = v
                cleaned_b.append(new_row)
            blast_bindings = cleaned_b
        except Exception:
            blast_bindings = []

    # Сводка BLAST по праймерам на человеке
    blast_human_summary = []
    blast_human_csv = os.path.join(tmpdir, "primer_human_blast_summary.csv")
    if os.path.exists(blast_human_csv) and os.path.getsize(blast_human_csv) > 0:
        try:
            df_bh = pd.read_csv(blast_human_csv)
            # предварительно заменяем NaN в DataFrame
            df_bh = df_bh.where(pd.notnull(df_bh), None)
            tmp = df_bh.to_dict(orient="records")
            # дополнительная защита: заменяем все float('nan') на None вручную
            cleaned = []
            for row in tmp:
                new_row = {}
                for k, v in row.items():
                    try:
                        import math
                        if isinstance(v, float) and math.isnan(v):
                            new_row[k] = None
                        else:
                            new_row[k] = v
                    except Exception:
                        new_row[k] = v
                cleaned.append(new_row)
            blast_human_summary = cleaned
        except Exception:
            blast_human_summary = []

    # Таблица продуктов (все строки результатов)
    products = df.to_dict(orient="records") if not df.empty else []
    for p in products:
        try:
            p["product_id"] = stable_product_id(
                p.get("Forward_primer", ""),
                p.get("Reverse_primer", ""),
                int(p.get("Amplicon_start")),
                int(p.get("Amplicon_end")),
            )
        except Exception:
            continue

    conservation = []
    if cons_path and bindings:
        try:
            conservation = compute_conservation_for_primers(bindings, products, primers_dict, cons_path)
        except Exception:
            conservation = []

    consensus_qc = []
    if products and conservation:
        try:
            consensus_qc = compute_consensus_qc_for_products(products, conservation, params)
            # сохраняем в run_dir, чтобы можно было скачать
            try:
                out_path = os.path.join(run_dir, "product_consensus_qc.csv")
                pd.DataFrame(consensus_qc).to_csv(out_path, sep=";", index=False)
            except Exception:
                pass
        except Exception:
            consensus_qc = []

    # Таблица гомодимеров/шпилек только для праймеров, которые используются в продуктах,
    # с пометкой, проходит ли КАЖДЫЙ праймер сам по себе по параметрам сравнения с консенсусом
    # (mismatches_total / mismatches_3prime / отсутствие инделов), независимо от продуктов.
    try:
        if homo_qc and products and conservation:
            used_primers = set()
            for p in products:
                fp = p.get("Forward_primer")
                rp = p.get("Reverse_primer")
                if fp:
                    used_primers.add(fp)
                if rp:
                    used_primers.add(rp)

            # считаем статус по консервативности для каждого праймера отдельно,
            # используя те же пороги, что и для QC по продуктам
            max_total = int((params or {}).get("MAX_MISMATCHES_TOTAL", MAX_MISMATCHES_TOTAL))
            max_3p = int((params or {}).get("MAX_MISMATCHES_3PRIME", MAX_MISMATCHES_3PRIME))
            win3 = int((params or {}).get("THREE_PRIME_WINDOW", THREE_PRIME_WINDOW))

            def _opt_int(key):
                v = (params or {}).get(key, "")
                if v is None or v == "":
                    return None
                try:
                    return int(v)
                except Exception:
                    return None

            max_indels = _opt_int("MAX_INDELS_TOTAL")
            max_indels_3p = _opt_int("MAX_INDELS_3PRIME")

            primer_pass = {}
            for r in conservation or []:
                name = r.get("Primer")
                if not name:
                    continue
                stats = mismatch_stats_vs_consensus(
                    r.get("Primer_seq", ""),
                    r.get("Consensus_seq", ""),
                    three_prime_window=win3,
                )
                passed = (stats["Mismatches_total"] <= max_total) and (stats["Mismatches_3prime"] <= max_3p)
                if max_indels is not None and stats["Indels_total"] > max_indels:
                    passed = False
                if max_indels_3p is not None and stats["Indels_3prime"] > max_indels_3p:
                    passed = False
                primer_pass[name] = bool(passed)

            filtered = []
            for row in homo_qc:
                pname = row.get("Primer")
                if pname not in used_primers:
                    continue
                status = "NA"
                if pname in primer_pass:
                    status = "OK" if primer_pass[pname] else "FAIL"
                new_row = dict(row)
                new_row["Conservation_status"] = status
                filtered.append(new_row)

            if filtered:
                out_path = os.path.join(run_dir, "primer_homodimer_qc_with_conservation.csv")
                pd.DataFrame(filtered).to_csv(out_path, sep=";", index=False)
    except Exception:
        pass

    # Таблица продуктов "как PCR_results_primer3py", но с колонкой статуса по консервативности
    # (на основе QC по консенсусу/инделам, который показываем в UI)
    try:
        if products:
            merged_rows = []
            for i, p in enumerate(products):
                status = "NA"
                if consensus_qc and i < len(consensus_qc):
                    status = "OK" if bool(consensus_qc[i].get("Pass")) else "FAIL"
                merged = dict(p)
                merged["Conservation_status"] = status
                merged_rows.append(merged)
            out_path = os.path.join(run_dir, "PCR_results_primer3py_with_conservation.csv")
            pd.DataFrame(merged_rows).to_csv(out_path, sep=";", index=False)
    except Exception:
        pass

    # Данные для интерактивной схемы аннотаций
    annotation_map = {"length": 0, "id": "", "features": []}
    try:
        annotated_gb = find_existing(run_dir, ("annotated_primer3py.gb",))
        source_gb = find_existing(run_dir, LEGACY_REF_NAMES)
        gb_for_map = annotated_gb or source_gb
        if gb_for_map:
            annotation_map = extract_annotation_map(gb_for_map, products=products, bindings=bindings)
            try:
                with open(os.path.join(run_dir, "annotation_map.json"), "w", encoding="utf-8") as jf:
                    json.dump(annotation_map, jf, ensure_ascii=False)
            except Exception:
                pass
    except Exception:
        annotation_map = {"length": 0, "id": "", "features": []}

    # Опциональная фильтрация продуктов по UI-правилам
    products_before = len(products or [])
    products = apply_product_filters(
        products, conservation, homo_qc, dimers_matrix, params, filters
    )
    # consensus_qc пересчитываем для отфильтрованного списка
    if products and conservation:
        try:
            consensus_qc = compute_consensus_qc_for_products(products, conservation, params)
        except Exception:
            consensus_qc = []
    else:
        consensus_qc = []

    annotation_map = rebuild_annotation_products(annotation_map, products)
    try:
        with open(os.path.join(run_dir, "annotation_map.json"), "w", encoding="utf-8") as jf:
            json.dump(annotation_map, jf, ensure_ascii=False)
    except Exception:
        pass

    # dimers для UI — только по оставшимся продуктам
    if products:
        product_keys = {(p.get("Forward_primer"), p.get("Reverse_primer")) for p in products}
        dimers = [d for d in dimers if (d.get("primer1"), d.get("primer2")) in product_keys]
    else:
        dimers = []

    # сохраняем отфильтрованные продукты + QC
    try:
        pd.DataFrame(products).to_csv(
            os.path.join(run_dir, "PCR_results_filtered.csv"), sep=";", index=False
        )
    except Exception:
        pass
    try:
        if consensus_qc:
            pd.DataFrame(consensus_qc).to_csv(
                os.path.join(run_dir, "product_consensus_qc.csv"), sep=";", index=False
            )
    except Exception:
        pass

    response_payload = {
        "status": "ok",
        "blast_enabled": blast_enabled,
        "blast_human_enabled": blast_human_enabled,
        "dimers": dimers,
        "has_results": bool(products),
        "bindings": bindings,
        "dimers_matrix": dimers_matrix,
        "products": products,
        "products_before_filter": products_before,
        "products_after_filter": len(products or []),
        "consensus_qc": consensus_qc,
        "homo_qc": homo_qc,
        "blast_bindings": blast_bindings,
        "blast_human_summary": blast_human_summary,
        "run_id": run_id,
        "available_files": sorted([f for f in os.listdir(run_dir) if os.path.isfile(os.path.join(run_dir, f))]),
        "conservation": conservation,
        "annotation_map": annotation_map,
        "filters_applied": filters,
        "params": params,
        "pcr_schemes": [],
        "lane_maps": {"ref": {}, "primer": {}, "product": {}},
    }

    # Сохраняем снимок сессии без тяжёлой dimers_matrix (она в CSV на диске)
    try:
        session_to_save = slim_session_payload(response_payload)
        with open(os.path.join(run_dir, "session.json"), "w", encoding="utf-8") as sf:
            json.dump(session_to_save, sf, ensure_ascii=False)
        response_payload["available_files"] = sorted(
            [f for f in os.listdir(run_dir) if os.path.isfile(os.path.join(run_dir, f))]
        )
    except Exception:
        pass

    # очистка временного каталога пайплайна
    try:
        shutil.rmtree(tmpdir, ignore_errors=True)
    except Exception:
        pass

    return jsonify(response_payload)


@app.post("/api/run/<run_id>/update_session")
def update_session(run_id):
    """Обновить session.json (например, схемы ПЦР), не пересчитывая пайплайн."""
    run_dir = os.path.join(RUNS_DIR, run_id)
    if not os.path.isdir(run_dir):
        return jsonify({"error": "run not found"}), 404
    session_path = os.path.join(run_dir, "session.json")
    payload = request.get_json(silent=True) or {}
    session = {}
    if os.path.exists(session_path):
        try:
            with open(session_path, "r", encoding="utf-8") as f:
                session = json.load(f)
        except Exception:
            session = {}
    if "pcr_schemes" in payload:
        session["pcr_schemes"] = payload.get("pcr_schemes") or []
    # разрешаем точечно обновлять и другие ключи UI-состояния при необходимости
    for key in ("lane_maps",):
        if key in payload:
            session[key] = payload.get(key)
    session["run_id"] = run_id
    session = slim_session_payload(session)
    try:
        with open(session_path, "w", encoding="utf-8") as f:
            json.dump(session, f, ensure_ascii=False)
    except Exception as e:
        return jsonify({"error": str(e)}), 500
    return jsonify({
        "status": "ok",
        "run_id": run_id,
        "pcr_schemes": session.get("pcr_schemes", []),
        "lane_maps": session.get("lane_maps", {}),
    })


@app.get("/api/export/<run_id>")
def export_run(run_id):
    """Скачать весь прогон как zip (включая session.json)."""
    run_dir = os.path.join(RUNS_DIR, run_id)
    if not os.path.isdir(run_dir):
        return jsonify({"error": "run not found"}), 404

    # убедимся, что session.json есть
    session_path = os.path.join(run_dir, "session.json")
    if not os.path.exists(session_path):
        return jsonify({"error": "session.json not found — этот прогон нельзя восстановить"}), 400

    zip_path = os.path.join(tempfile.gettempdir(), f"in_silico_pcr_run_{run_id}.zip")
    try:
        with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
            for root, _dirs, files in os.walk(run_dir):
                for fn in files:
                    full = os.path.join(root, fn)
                    arc = os.path.relpath(full, run_dir)
                    zf.write(full, arcname=arc)
        return send_file(
            zip_path,
            as_attachment=True,
            download_name=f"in_silico_pcr_run_{run_id}.zip",
        )
    except Exception as e:
        return jsonify({"error": str(e)}), 500


@app.post("/api/load_run")
def load_run():
    """Загрузить ранее сохранённый zip-прогон и вернуть session.json как результат."""
    zf_file = request.files.get("run_archive")
    if zf_file is None:
        return jsonify({"error": "run_archive is required"}), 400

    run_id = uuid.uuid4().hex[:12]
    run_dir = os.path.join(RUNS_DIR, run_id)
    os.makedirs(run_dir, exist_ok=True)

    tmp_zip = os.path.join(tempfile.gettempdir(), f"upload_{run_id}.zip")
    try:
        zf_file.save(tmp_zip)
        with zipfile.ZipFile(tmp_zip, "r") as zf:
            safe_extract_zip(zf, run_dir)
    except Exception as e:
        shutil.rmtree(run_dir, ignore_errors=True)
        return jsonify({"error": f"Не удалось распаковать архив: {e}"}), 400
    finally:
        try:
            os.remove(tmp_zip)
        except Exception:
            pass

    session_path = os.path.join(run_dir, "session.json")
    # иногда session.json лежит во вложенной папке
    if not os.path.exists(session_path):
        for root, _dirs, files in os.walk(run_dir):
            if "session.json" in files:
                session_path = os.path.join(root, "session.json")
                # если распаковалось во вложенную папку — поднимем файлы наверх
                nested = root
                if os.path.abspath(nested) != os.path.abspath(run_dir):
                    for fn in os.listdir(nested):
                        src = os.path.join(nested, fn)
                        dst = os.path.join(run_dir, fn)
                        if not os.path.exists(dst):
                            shutil.move(src, dst)
                    session_path = os.path.join(run_dir, "session.json")
                break

    if not os.path.exists(session_path):
        return jsonify({"error": "В архиве нет session.json"}), 400

    try:
        with open(session_path, "r", encoding="utf-8") as f:
            payload = json.load(f)
    except Exception as e:
        return jsonify({"error": f"Не удалось прочитать session.json: {e}"}), 400

    payload["run_id"] = run_id
    # восстановить dimers_matrix из CSV, если выкинули из session
    if not payload.get("dimers_matrix"):
        payload["dimers_matrix"] = load_dimers_matrix_from_run(run_dir)
    # стабильные product_id для старых прогонов без поля
    for p in payload.get("products") or []:
        try:
            fp = p.get("Forward_primer", "")
            rp = p.get("Reverse_primer", "")
            a_start = int(p.get("Amplicon_start"))
            a_end = int(p.get("Amplicon_end"))
            p["product_id"] = p.get("product_id") or stable_product_id(fp, rp, a_start, a_end)
        except Exception:
            continue
    # пересобрать annotation products со стабильными id, если карта есть
    if payload.get("annotation_map") and payload.get("products") is not None:
        try:
            payload["annotation_map"] = rebuild_annotation_products(
                payload.get("annotation_map"), payload.get("products")
            )
        except Exception:
            pass

    # lane_maps already in session payload if present — leave as-is
    payload["available_files"] = sorted(
        [f for f in os.listdir(run_dir) if os.path.isfile(os.path.join(run_dir, f))]
    )
    # обновим session с новым run_id (без dimers_matrix)
    try:
        with open(os.path.join(run_dir, "session.json"), "w", encoding="utf-8") as f:
            json.dump(slim_session_payload(payload), f, ensure_ascii=False)
    except Exception:
        pass

    LATEST_RUN["id"] = run_id
    LATEST_RUN["dir"] = run_dir
    return jsonify(payload)


@app.get("/")
def index():
    index_path = os.path.join(app.static_folder, "index.html")
    if os.path.exists(index_path):
        return send_from_directory(app.static_folder, "index.html")
    return (
        "<h1>React‑сборка не найдена</h1>"
        "<p>Соберите фронтенд в папке <code>frontend</code> с помощью <code>npm run build</code>, "
        "чтобы web_app.py мог отдать готовый интерфейс.</p>"
    )


@app.get("/api/download/<run_id>/<path:filename>")
def download_file(run_id, filename):
    run_dir = os.path.join(RUNS_DIR, run_id)
    if not os.path.isdir(run_dir):
        return jsonify({"error": "run not found"}), 404
    file_path = os.path.join(run_dir, filename)
    # защита от выхода за каталог
    if not os.path.abspath(file_path).startswith(os.path.abspath(run_dir) + os.sep):
        return jsonify({"error": "invalid path"}), 400
    if not os.path.exists(file_path):
        return jsonify({"error": "file not found"}), 404
    return send_from_directory(run_dir, filename, as_attachment=True)


if __name__ == "__main__":
    ensure_frontend_built()
    cleanup_old_runs(RUNS_DIR)

    host = os.environ.get("HOST", "127.0.0.1")
    port = int(os.environ.get("PORT", "5000"))
    # Open browser only for local interactive runs (not on a VPS).
    open_browser = os.environ.get(
        "OPEN_BROWSER",
        "1" if host in ("127.0.0.1", "localhost") else "0",
    ) == "1"

    if open_browser:
        import webbrowser
        from threading import Timer

        def _open():
            webbrowser.open(f"http://{host}:{port}/")

        Timer(1.0, _open).start()

    # debug=False, чтобы отключить авто-перезапуск (watchdog),
    # который может рвать запросы во время длительного анализа
    app.run(host=host, port=port, debug=False)
