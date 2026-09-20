import pandas as pd

from .iupac import (
    consensus_iupac_from_row,
    iupac_complement,
    mismatch_stats_vs_consensus,
)

# Defaults synchronized with web_app / In_silico_PCR.py
MAX_MISMATCHES_TOTAL = 0
MAX_MISMATCHES_3PRIME = 0
THREE_PRIME_WINDOW = 5


def compute_conservation_for_primers(bindings, products, primers_dict, cons_table_path):
    """
    bindings: list of dicts from primer_binding_sites.csv (filtered bindings)
    primers_dict: name->sequence (iupac/original), 5'->3'
    returns: list of per-primer reports
    """
    # читаем как строки, чтобы не потерять координаты вида "7+1"
    try:
        df = pd.read_csv(cons_table_path, sep="\t", dtype={"Coordinates": str})
    except Exception:
        df = pd.read_csv(cons_table_path, sep=None, engine="python", dtype={"Coordinates": str})

    # карты координат:
    # - int_rows: coordinate(int) -> row dict
    # - ins_rows: base_coordinate(int) -> list of (token, row dict) in appearance order, где token выглядит как "7+1"
    int_rows = {}
    ins_rows = {}

    def parse_coord_token(tok):
        tok = str(tok).strip()
        if "+" in tok:
            base, suf = tok.split("+", 1)
            try:
                return ("ins", int(base), tok)
            except Exception:
                return (None, None, tok)
        try:
            return ("int", int(tok), tok)
        except Exception:
            return (None, None, tok)

    for _, r in df.iterrows():
        kind, base, token = parse_coord_token(r.get("Coordinates"))
        if kind is None:
            continue
        rowd = r.to_dict()
        if kind == "int":
            int_rows[base] = rowd
        else:
            ins_rows.setdefault(base, []).append((token, rowd))

    # Индексируем посадки по праймеру для быстрого поиска координат
    bindings_by_primer = {}
    for b in bindings:
        pname = b.get("Primer")
        if not pname:
            continue
        bindings_by_primer.setdefault(pname, []).append(b)

    # Выбираем посадку, которая реально использовалась в продуктах (PCR),
    # иначе fallback: "лучшая" по mismatches.
    chosen = {}

    # сначала собираем кандидатов из продуктов (берём самый короткий продукт для детерминизма)
    try:
        products_sorted = sorted(products or [], key=lambda r: float(r.get("Product_size", 1e18)))
    except Exception:
        products_sorted = products or []

    for prod in products_sorted:
        fp = prod.get("Forward_primer")
        rp = prod.get("Reverse_primer")
        a_start = prod.get("Amplicon_start")
        a_end = prod.get("Amplicon_end")
        try:
            a_start = int(a_start)
            a_end = int(a_end)
        except Exception:
            continue

        # forward primer: посадка начинается в Amplicon_start на '+' цепи
        if fp and fp not in chosen:
            for b in bindings_by_primer.get(fp, []):
                try:
                    if b.get("Strand") == "+" and int(b.get("Start")) == a_start:
                        chosen[fp] = {"Start": int(b.get("Start")), "End": int(b.get("End")), "Strand": b.get("Strand")}
                        break
                except Exception:
                    continue

        # reverse primer: посадка заканчивается в Amplicon_end на '-' цепи
        if rp and rp not in chosen:
            for b in bindings_by_primer.get(rp, []):
                try:
                    if b.get("Strand") == "-" and int(b.get("End")) == a_end:
                        chosen[rp] = {"Start": int(b.get("Start")), "End": int(b.get("End")), "Strand": b.get("Strand")}
                        break
                except Exception:
                    continue

    # fallback: лучшая по mismatches
    if len(chosen) < len(bindings_by_primer):
        for pname, rows in bindings_by_primer.items():
            if pname in chosen:
                continue
            best = None
            best_key = None
            for b in rows:
                try:
                    mm = int(b.get("Total_mismatches", 0))
                    mm3 = int(b.get("Mismatches_3prime", 0))
                    key = (mm, mm3)
                except Exception:
                    continue
                if best is None or key < best_key:
                    try:
                        best = {"Start": int(b.get("Start")), "End": int(b.get("End")), "Strand": b.get("Strand")}
                        best_key = key
                    except Exception:
                        continue
            if best is not None:
                chosen[pname] = best

    reports = []
    for name, hit in chosen.items():
        primer_seq = str(primers_dict.get(name, "")).upper()
        if not primer_seq:
            continue
        start = hit["Start"]
        end = hit["End"]
        strand = hit.get("Strand", "+")
        L = len(primer_seq)
        # входные координаты bindings теперь 1-based inclusive
        if (end - start + 1) != L:
            # не должно быть, но пропустим
            continue

        # строим выравнивание "праймер vs таблица координат":
        # для каждой базовой координаты (int) добавляем одну букву праймера,
        # а для вставок вида k+1 добавляем '-' (индел) в праймер и консенсус/S_index по строке вставки.
        primer_aln = []
        consensus_aln = []
        s_index_aln = []
        s_index_for_avg = []

        # последовательность опорных координат (int) по направлению посадки
        if strand == "+":
            coords = list(range(start, end + 1))
            primer_iter = list(primer_seq)
        else:
            coords = list(range(end, start - 1, -1))
            primer_iter = list(primer_seq)

        pi = 0
        for base_coord in coords:
            base_row = int_rows.get(base_coord)
            pch = primer_iter[pi] if pi < len(primer_iter) else " "
            pi += 1
            primer_aln.append(pch)
            if base_row is None:
                consensus_aln.append(" ")
                s_index_aln.append(None)
            else:
                cons_ch = consensus_iupac_from_row(base_row, min_freq=0.05)
                if strand == "-":
                    cons_ch = iupac_complement(cons_ch)
                consensus_aln.append(cons_ch)
                try:
                    v = float(base_row.get("S_index"))
                except Exception:
                    v = None
                s_index_aln.append(v)
                if isinstance(v, (int, float)):
                    s_index_for_avg.append(v)

            # вставки после этой координаты
            for _tok, ins_row in ins_rows.get(base_coord, []):
                cons_ch = consensus_iupac_from_row(ins_row, min_freq=0.05)
                # инделы, которые "пустые" (консенсус = пробел),
                # считаем несущественными и полностью пропускаем
                if cons_ch == " ":
                    continue
                primer_aln.append("-")
                if strand == "-":
                    cons_ch = iupac_complement(cons_ch)
                consensus_aln.append(cons_ch)
                try:
                    s_index_aln.append(float(ins_row.get("S_index")))
                except Exception:
                    s_index_aln.append(None)

        avg = None
        if s_index_for_avg:
            avg = sum(s_index_for_avg) / len(s_index_for_avg)

        reports.append({
            "Primer": name,
            "Strand": strand,
            "Start": start,
            "End": end,
            "Primer_seq": "".join(primer_aln),
            "Consensus_seq": "".join(consensus_aln),
            "S_index": s_index_aln,
            "S_index_avg": avg,
        })
    return reports


def compute_consensus_qc_for_products(products, conservation_reports, params):
    cons_by_primer = {r.get("Primer"): r for r in (conservation_reports or []) if r.get("Primer")}

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

    def primer_qc(primer_name: str):
        r = cons_by_primer.get(primer_name)
        if not r:
            return {"Primer": primer_name, "Status": "no_data"}
        stats = mismatch_stats_vs_consensus(r.get("Primer_seq", ""), r.get("Consensus_seq", ""), three_prime_window=win3)
        passed = (stats["Mismatches_total"] <= max_total) and (stats["Mismatches_3prime"] <= max_3p)
        if max_indels is not None and stats["Indels_total"] > max_indels:
            passed = False
        if max_indels_3p is not None and stats["Indels_3prime"] > max_indels_3p:
            passed = False
        return {
            "Primer": primer_name,
            "Matches": stats["Matches"],
            "Mismatches_total": stats["Mismatches_total"],
            "Mismatches_3prime": stats["Mismatches_3prime"],
            "Indels_total": stats["Indels_total"],
            "Indels_3prime": stats["Indels_3prime"],
            "Has_indels": stats["Has_indels"],
            "Pass": bool(passed),
        }

    rows = []
    for p in products or []:
        fp = p.get("Forward_primer")
        rp = p.get("Reverse_primer")
        fqc = primer_qc(fp) if fp else {"Primer": None, "Status": "no_data"}
        rqc = primer_qc(rp) if rp else {"Primer": None, "Status": "no_data"}
        overall = bool(fqc.get("Pass")) and bool(rqc.get("Pass"))
        rows.append({
            "Forward_primer": fp,
            "Reverse_primer": rp,
            "Amplicon_start": p.get("Amplicon_start"),
            "Amplicon_end": p.get("Amplicon_end"),
            "Product_size": p.get("Product_size"),
            "F_Matches": fqc.get("Matches"),
            "F_Mismatches_total": fqc.get("Mismatches_total"),
            "F_Mismatches_3prime": fqc.get("Mismatches_3prime"),
            "F_Indels_total": fqc.get("Indels_total"),
            "F_Indels_3prime": fqc.get("Indels_3prime"),
            "F_Has_indels": fqc.get("Has_indels"),
            "F_Pass": fqc.get("Pass"),
            "R_Matches": rqc.get("Matches"),
            "R_Mismatches_total": rqc.get("Mismatches_total"),
            "R_Mismatches_3prime": rqc.get("Mismatches_3prime"),
            "R_Indels_total": rqc.get("Indels_total"),
            "R_Indels_3prime": rqc.get("Indels_3prime"),
            "R_Has_indels": rqc.get("Has_indels"),
            "R_Pass": rqc.get("Pass"),
            "Pass": overall,
        })
    return rows
