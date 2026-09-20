#!/usr/bin/env python3
"""
FULL IN SILICO PCR PIPELINE (С ПОДРОБНЫМ ОТЧЁТОМ И ПРОГРЕССОМ)

В одном скрипте:
- Primer3 термодинамика
- Hairpin / Homodimer / Heterodimer
- 3'-complementarity
- Поиск ампликонов на целевом геноме
- BLAST праймеров против человека
- QC неспецифических ампликонов
- RISK помечаются, НЕ удаляются
- Подробный отчёт по каждому праймеру
- Подробный вывод прогресса в консоль
"""

import pandas as pd
from itertools import product
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
import primer3
import subprocess, os, sys

# =============================================================
# НАСТРОЙКИ
# =============================================================

PRIMER_CSV = os.environ.get("PRIMER_CSV", "primers.csv")
TARGET_GB = os.environ.get("TARGET_GB", "reference.gb")

BASE_DIR = os.path.dirname(__file__)


def _resolve_blastn_path():
    """Windows bundle, Linux local bin, PATH, or BLASTN_PATH env."""
    override = os.environ.get("BLASTN_PATH")
    if override:
        return override
    candidates = [
        os.path.join(BASE_DIR, "ncbi-blast-2.17.0+", "bin", "blastn.exe"),
        os.path.join(BASE_DIR, "ncbi-blast-2.17.0+", "bin", "blastn"),
    ]
    for path in candidates:
        if os.path.exists(path):
            return path
    from shutil import which
    found = which("blastn")
    return found or candidates[0]


BLASTN_PATH = _resolve_blastn_path()
# для BLAST по человеку используем ОТНОСИТЕЛЬНЫЙ путь к базе,
# а cwd указываем BASE_DIR, чтобы избежать проблем с пробелами/скобками в пути
DB_HUMAN = os.path.join("human_genome", "GCF_000001405.39_top_level")

# Базовые значения настроек (по умолчанию)
MIN_AMPLICON_DEFAULT = 100
MAX_AMPLICON_DEFAULT = 10000
MAX_DEGENERATE_VARIANTS_DEFAULT = 50

MAX_TM_DIFF_DEFAULT = 7.0
MAX_3PRIME_MATCHES_DEFAULT = 3

ANNEALING_TEMP_DEFAULT = 60.0
MAX_HUMAN_DELTA_DEFAULT = 1500

MAX_MISMATCHES_TOTAL_DEFAULT = 0
MAX_MISMATCHES_3PRIME_DEFAULT = 0
THREE_PRIME_WINDOW_DEFAULT = 5

# Текущие значения могут переопределяться через переменные окружения
MIN_AMPLICON = int(os.getenv("MIN_AMPLICON", str(MIN_AMPLICON_DEFAULT)))
MAX_AMPLICON = int(os.getenv("MAX_AMPLICON", str(MAX_AMPLICON_DEFAULT)))
MAX_DEGENERATE_VARIANTS = int(os.getenv("MAX_DEGENERATE_VARIANTS", str(MAX_DEGENERATE_VARIANTS_DEFAULT)))

MAX_TM_DIFF = float(os.getenv("MAX_TM_DIFF", str(MAX_TM_DIFF_DEFAULT)))
MAX_3PRIME_MATCHES = int(os.getenv("MAX_3PRIME_MATCHES", str(MAX_3PRIME_MATCHES_DEFAULT)))

ANNEALING_TEMP = float(os.getenv("ANNEALING_TEMP", str(ANNEALING_TEMP_DEFAULT)))
MAX_HUMAN_DELTA = int(os.getenv("MAX_HUMAN_DELTA", str(MAX_HUMAN_DELTA_DEFAULT)))

MAX_MISMATCHES_TOTAL = int(os.getenv("MAX_MISMATCHES_TOTAL", str(MAX_MISMATCHES_TOTAL_DEFAULT)))
MAX_MISMATCHES_3PRIME = int(os.getenv("MAX_MISMATCHES_3PRIME", str(MAX_MISMATCHES_3PRIME_DEFAULT)))
THREE_PRIME_WINDOW = int(os.getenv("THREE_PRIME_WINDOW", str(THREE_PRIME_WINDOW_DEFAULT)))

OUT_CSV    = "PCR_results_primer3py.csv"
OUT_GB     = "annotated_primer3py.gb"
OUT_REPORT = "primer3_dimer_report.txt"

# =============================================================
# ТЕРМОДИНАМИКА
# =============================================================

PCR_PARAMS_TM = {
    'mv_conc': 50.0,
    'dv_conc': 1.5,
    'dntp_conc': 0.6,
    'dna_conc': 250.0,
}

PCR_PARAMS_STRUCT = {**PCR_PARAMS_TM, 'temp_c': 37.0}

# =============================================================
# ЗАГРУЗКА РЕФЕРЕНСА (GenBank/FASTA)
# =============================================================


def load_reference_record(path):
    """
    Пытается прочитать референс как GenBank, если не получается — как FASTA.
    Бросает RuntimeError с понятным сообщением, если ни один формат не подошёл.
    """
    path = str(path)
    tried = []

    for fmt in ("genbank", "fasta"):
        try:
            record = next(SeqIO.parse(path, fmt))
            # Biopython требует molecule_type для записи в GenBank
            if "molecule_type" not in record.annotations:
                record.annotations["molecule_type"] = "DNA"
            return record
        except StopIteration:
            tried.append(fmt)
            continue

    raise RuntimeError(
        f"Не удалось прочитать файл референса '{path}' ни как GenBank, ни как FASTA "
        f"(проверены форматы: {', '.join(tried)}). Убедитесь, что файл не пустой и имеет корректный формат."
    )

# =============================================================
# IUPAC
# =============================================================

IUPAC_TO_DNA = {
    "A": ["A"], "C": ["C"], "G": ["G"], "T": ["T"],
    "R": ["A", "G"], "Y": ["C", "T"], "S": ["G", "C"],
    "W": ["A", "T"], "K": ["G", "T"], "M": ["A", "C"],
    "B": ["C", "G", "T"], "D": ["A", "G", "T"],
    "H": ["A", "C", "T"], "V": ["A", "C", "G"],
    "N": ["A", "C", "G", "T"]
}

# =============================================================
# УТИЛИТЫ
# =============================================================

def infer_direction_from_bindings(binding_sites, fallback="Forward"):
    """
    Определяет направление праймера по ЛУЧШЕЙ посадке на референс.
    Критерий "лучшести" совпадает с используемым при анализе продуктов:
      1) минимальные mismatches на 3'-конце
      2) затем минимальные mismatches по всему праймеру
      3) затем координаты (для детерминизма)

    Это важно, когда праймер имеет посадки на обе цепи: направление должно
    соответствовать наиболее вероятной (лучшей) посадке, а не "по умолчанию".
    """
    if not binding_sites:
        return fallback

    try:
        # hit: (start, end, strand, used_variant, total_mismatches, mismatches_3prime)
        best = min(binding_sites, key=lambda h: (h[5], h[4], h[0], h[1]))
        return "Forward" if best[2] == "+" else "Reverse"
    except Exception:
        return fallback


def expand_iupac_to_dna(seq, max_variants=50):
    pools = [IUPAC_TO_DNA.get(c.upper(), IUPAC_TO_DNA['N']) for c in seq]
    variants = []
    for i, comb in enumerate(product(*pools)):
        if i >= max_variants:
            break
        variants.append("".join(comb))
    return variants


def normalize_tm(tm):
    if tm is None or tm < 0:
        return 0.0
    return float(tm)

# =============================================================
# PRIMER3 АНАЛИЗ
# =============================================================

def safe_primer_analysis(seq):
    res = {'tm': None, 'hairpin_tm': 0, 'hairpin_dg': 0, 'homodimer_tm': 0, 'homodimer_dg': 0}

    tm = primer3.calc_tm(seq, **PCR_PARAMS_TM)
    res['tm'] = float(tm)

    hp = primer3.calc_hairpin(seq, **PCR_PARAMS_STRUCT)
    if hp.structure_found and hp.tm is not None:
        res['hairpin_tm'] = normalize_tm(hp.tm)
        res['hairpin_dg'] = float(hp.dg) / 1000

    hd = primer3.calc_homodimer(seq, **PCR_PARAMS_STRUCT)
    if hd.structure_found and hd.tm is not None:
        res['homodimer_tm'] = normalize_tm(hd.tm)
        res['homodimer_dg'] = float(hd.dg) / 1000

    return res


def heterodimer_analysis(seq_f, seq_r):
    res = {'heterodimer_tm': 0, 'heterodimer_dg': 0}
    hd = primer3.calc_heterodimer(seq_f, seq_r, **PCR_PARAMS_STRUCT)
    if hd.structure_found and hd.tm is not None:
        res['heterodimer_tm'] = normalize_tm(hd.tm)
        res['heterodimer_dg'] = float(hd.dg) / 1000
    return res

# =============================================================
# 3'-CHECK
# =============================================================

def check_3prime_complementarity(seq1, seq2, max_len=5):
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    end1 = seq1[-max_len:]
    end2 = seq2[-max_len:]
    matches = 0
    for i in range(min(len(end1), len(end2))):
        if comp.get(end1[i], '') == end2[::-1][i]:
            matches += 1
    return matches

# =============================================================
# ПОИСК ПОСАДОК НА ЦЕЛЕВОМ ГЕНОМЕ
# =============================================================

def _mismatch_counts(primer, template_segment, window_3prime=5):
    """
    Возвращает (total_mismatches, mismatches_3prime) между праймером и сегментом шаблона.
    Предполагается, что длины primer и template_segment совпадают и уже приведены к верхнему регистру.
    """
    total = 0
    mismatches_3p = 0

    L = len(primer)
    w = min(window_3prime, L)

    # последние w оснований — это 3'-конец праймера
    for i, (p, t) in enumerate(zip(primer, template_segment)):
        if p != t:
            total += 1
            # позиции, попадающие в 3'-окно
            if i >= L - w:
                mismatches_3p += 1

    return total, mismatches_3p


def find_binding_sites(template, dna_variants, max_variants=5,
                       max_mismatches_total=None,
                       max_mismatches_3prime=None,
                       window_3prime=None):
    """
    Поиск посадок праймера на целевом геноме с учётом допуска по mismatches.

    Возвращает список кортежей:
      (start, end, strand, used_variant, total_mismatches, mismatches_3prime)
    """
    hits = []
    template = str(template).upper()

    if max_mismatches_total is None:
        max_mismatches_total = MAX_MISMATCHES_TOTAL
    if max_mismatches_3prime is None:
        max_mismatches_3prime = MAX_MISMATCHES_3PRIME
    if window_3prime is None:
        window_3prime = THREE_PRIME_WINDOW

    for var in dna_variants[:max_variants]:
        var = var.upper()
        L = len(var)

        # прямой поиск (праймер на + цепи)
        for pos in range(0, len(template) - L + 1):
            segment = template[pos:pos + L]
            total_mm, mm_3p = _mismatch_counts(var, segment, window_3prime=window_3prime)
            if total_mm <= max_mismatches_total and mm_3p <= max_mismatches_3prime:
                hits.append((pos, pos + L, '+', var, total_mm, mm_3p))

        # поиск обратного комплемента (праймер на - цепи)
        rc = str(Seq(var).reverse_complement())
        for pos in range(0, len(template) - L + 1):
            segment = template[pos:pos + L]
            total_mm, mm_3p = _mismatch_counts(rc, segment, window_3prime=window_3prime)
            if total_mm <= max_mismatches_total and mm_3p <= max_mismatches_3prime:
                hits.append((pos, pos + L, '-', var, total_mm, mm_3p))

    uniq, seen = [], set()
    for h in hits:
        key = (h[0], h[1], h[2])
        if key not in seen:
            uniq.append(h)
            seen.add(key)

    return uniq

# =============================================================
# BLAST
# =============================================================

def blast_primer(seq):
    # временные файлы создаём в каталоге проекта, чтобы blastn их видел при cwd=BASE_DIR
    fasta = os.path.join(BASE_DIR, "tmp.fa")
    out   = os.path.join(BASE_DIR, "tmp.out")
    with open(fasta, "w") as f:
        f.write(">q\n" + seq + "\n")

    # если blastn.exe или база не найдены — тихо возвращаем пустой список,
    # чтобы не ронять весь пайплайн
    if not os.path.exists(BLASTN_PATH):
        print(f"[BLAST] Внимание: blastn не найден по пути: {BLASTN_PATH}. BLAST пропущен для этого праймера.")
        return []

    # проверка наличия базы: ищем хотя бы один файл, начинающийся с DB_HUMAN
    db_dir_rel = os.path.dirname(DB_HUMAN) or "."
    db_dir = os.path.join(BASE_DIR, db_dir_rel)
    db_prefix = os.path.basename(DB_HUMAN)
    if (not os.path.isdir(db_dir)) or (not any(fn.startswith(db_prefix) for fn in os.listdir(db_dir))):
        print(f"[BLAST] Внимание: база человека не найдена по префиксу: {os.path.join(db_dir_rel, db_prefix)}. "
              f"Ищем файлы в каталоге: {db_dir}")
        return []

    # Расширенный формат: добавляем mismatch и length для последующего анализа
    cmd = [BLASTN_PATH, "-task", "blastn-short", "-query", fasta, "-db", DB_HUMAN,
           "-word_size", "7", "-perc_identity", "70", "-evalue", "100",
           "-outfmt", "6 sseqid sstart send mismatch length", "-out", out]

    try:
        print(f"[BLAST human] Seq: {seq}")
        print(f"[BLAST human] Running: {' '.join(cmd)}")
        proc = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
            cwd=BASE_DIR,
        )
        print(f"[BLAST human] returncode={proc.returncode}")
        if proc.stderr:
            print("[BLAST human stderr]:")
            # ограничим объём вывода stderr
            print(proc.stderr[:1000])
    except FileNotFoundError:
        print(f"[BLAST] Внимание: не удалось запустить blastn по пути: {BLASTN_PATH}.")
        return []

    hits = []
    if os.path.exists(out):
        with open(out) as f:
            for line in f:
                parts = line.split()
                if len(parts) < 5:
                    continue
                sid, sstart, send, mismatch, length = parts
                hits.append({
                    'chrom': sid,
                    'sstart': int(sstart),
                    'send': int(send),
                    'mismatch': int(mismatch),
                    'length': int(length),
                })

    print(f"[BLAST human] Hits found: {len(hits)}")

    for fn in [fasta, out]:
        try: os.remove(fn)
        except: pass

    return hits


def blast_primer_on_reference(seq, subject_path, template_seq):
    """
    BLAST праймера против референса (GenBank/FASTA) для поиска лучшей посадки.
    Использует blastn-short с ключом -subject (без предварительного makeblastdb).
    Возвращает словарь с полями:
      strand ('+' / '-'), qstart, qend, sstart, send,
      total_mismatches, mismatches_3prime
    или None, если хитов нет или BLAST недоступен.
    """
    # временные файлы для BLAST по референсу тоже создаём в каталоге проекта,
    # но сам BLAST запускаем в текущем cwd (там лежит референс)
    fasta = os.path.join(BASE_DIR, "tmp_ref.fa")
    out   = os.path.join(BASE_DIR, "tmp_ref.out")
    with open(fasta, "w") as f:
        f.write(">q\n" + seq + "\n")

    # BLAST -subject надёжнее работает с FASTA, поэтому генерируем временный subject.fa
    subject_fa = os.path.join(BASE_DIR, "tmp_subject_ref.fa")

    if not os.path.exists(BLASTN_PATH):
        for fn in [fasta]:
            try: os.remove(fn)
            except: pass
        return None

    try:
        with open(subject_fa, "w", encoding="ascii") as sf:
            sf.write(">ref\n")
            s = str(template_seq).upper()
            # пишем в FASTA по 60 символов на строку
            for i in range(0, len(s), 60):
                sf.write(s[i:i+60] + "\n")
    except Exception:
        for fn in [fasta]:
            try: os.remove(fn)
            except: pass
        return None

    # Требуем полное покрытие запроса и запрещаем гэпы (инделы)
    # + берём координаты и базовые метрики для выбора лучшего хита
    outfmt = "6 sstrand qstart qend sstart send mismatch length"
    cmd = [
        BLASTN_PATH,
        "-task", "blastn-short",
        "-query", fasta,
        "-subject", subject_fa,
        "-word_size", "7",
        "-perc_identity", "70",
        # на референсе сейчас жёстко требуем только отсутствие инделов
        "-ungapped",
        "-evalue", "100",
        "-outfmt", outfmt,
        "-out", out,
    ]

    try:
        print(f"[BLAST ref] Seq: {seq}")
        print(f"[BLAST ref] Subject length: {len(str(template_seq))}")
        print(f"[BLAST ref] Running: {' '.join(cmd)}")
        proc = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )
        print(f"[BLAST ref] returncode={proc.returncode}")
        if proc.stderr:
            print("[BLAST ref stderr]:")
            print(proc.stderr[:1000])
    except FileNotFoundError:
        for fn in [fasta, out]:
            try: os.remove(fn)
            except: pass
        return None

    best = None
    if os.path.exists(out):
        try:
            out_size = os.path.getsize(out)
        except Exception:
            out_size = -1
        print(f"[BLAST ref] out file size: {out_size}")
        with open(out) as f:
            for line in f:
                parts = line.strip().split()
                # outfmt: sstrand qstart qend sstart send mismatch length  -> 7 полей
                if len(parts) < 7:
                    continue
                sstrand, qstart, qend, sstart, send, mismatch, length = parts
                length = int(length)
                mismatch = int(mismatch)
                if best is None or length > best["length"] or (length == best["length"] and mismatch < best["blast_mismatches"]):
                    best = {
                        "strand": "+" if sstrand.startswith("plus") else "-",
                        "qstart": int(qstart),
                        "qend": int(qend),
                        "sstart": int(sstart),
                        "send": int(send),
                        "blast_mismatches": mismatch,
                        "length": length,
                    }

    if best is None:
        print(f"[BLAST ref] No hits for seq={seq}")

    for fn in [fasta, out]:
        try: os.remove(fn)
        except: pass
    try:
        os.remove(subject_fa)
    except:
        pass

    if best is None:
        return None

    # Пересчитываем mismatches по референсной последовательности тем же методом,
    # что и в "обычном" поиске посадок, чтобы метрики совпадали.
    template = str(template_seq).upper()
    # BLAST может вернуть неполное покрытие запроса; расширяем координаты до полной длины праймера
    primer_len = len(seq)
    qstart = int(best["qstart"])
    qend = int(best["qend"])
    s0 = int(best["sstart"])
    s1 = int(best["send"])
    step = 1 if s1 >= s0 else -1
    s0_full = s0 - step * (qstart - 1)
    s1_full = s1 + step * (primer_len - qend)

    start = min(s0_full, s1_full) - 1  # BLAST 1-based inclusive
    end = max(s0_full, s1_full)        # python slice end exclusive
    if start < 0 or end > len(template):
        print(f"[BLAST ref] Out of range after extension: s0={s0} s1={s1} q={qstart}-{qend} -> full={s0_full}-{s1_full}, template_len={len(template)}")
        return None

    segment = template[start:end]
    if best["strand"] == "-":
        segment = str(Seq(segment).reverse_complement())

    if len(segment) != len(seq):
        print(f"[BLAST ref] Segment length mismatch after extension: got={len(segment)} expected={len(seq)}")
        return None

    total_mm, mm_3p = _mismatch_counts(seq.upper(), segment.upper(), window_3prime=THREE_PRIME_WINDOW)
    best["total_mismatches"] = total_mm
    best["mismatches_3prime"] = mm_3p
    best["sstart_full"] = s0_full
    best["send_full"] = s1_full
    return best

# =============================================================
# ОСНОВНОЙ КОД
# =============================================================

print("="*60)
print("IN SILICO PCR PIPELINE STARTED")
print("="*60)

primers_df = pd.read_csv(PRIMER_CSV, sep=";")
record = load_reference_record(TARGET_GB)
template = record.seq

report = open(OUT_REPORT, "w", encoding="utf-8")
report.write("=== ОТЧЁТ ПО ПРАЙМЕРАМ (термодинамика + структуры) ===\n")
report.write(f"Ta = {ANNEALING_TEMP} °C\n")
report.write("Критерий: Tm_struct <= min(Tm_primer, Ta) - 15\n\n")

primer_db = {}

# ---------- АНАЛИЗ ПРАЙМЕРОВ ----------

print(f"[STEP 1] Анализ праймеров: {len(primers_df)} шт")

for idx, row in primers_df.iterrows():
    name = row["Name"]
    seq_iupac = row["Sequence"]
    ptype = row["Type"]

    print(f"  [{idx+1}/{len(primers_df)}] Primer {name}")

    dna_variants = expand_iupac_to_dna(seq_iupac, MAX_DEGENERATE_VARIANTS)
    variant = dna_variants[0]

    analysis = safe_primer_analysis(variant)
    binding_sites = find_binding_sites(template, dna_variants)

    primer_tm = analysis['tm']
    limit = min(primer_tm, ANNEALING_TEMP) - 15.0

    # --- ОТЧЁТ ПО ПРАЙМЕРУ ---
    report.write(f"{name} ({ptype})\n")
    report.write(f"  Seq: {variant}\n")
    report.write(f"  Tm primer = {primer_tm:.2f}\n")

    # Hairpin
    hp_tm = analysis['hairpin_tm']
    hp_dg = analysis['hairpin_dg']

    if hp_tm >= primer_tm:
        report.write(f"  ❌ Hairpin Tm={hp_tm:.2f}, ΔG={hp_dg:.2f}  (СТРУКТУРА СТАБИЛЬНЕЕ ПРАЙМЕРА)\n")
    elif hp_tm > limit:
        report.write(f"  ❌ Hairpin Tm={hp_tm:.2f}, ΔG={hp_dg:.2f}  (выше лимита {limit:.2f})\n")
    else:
        report.write(f"  ✔ Hairpin OK  Tm={hp_tm:.2f}, ΔG={hp_dg:.2f}\n")

    # Homodimer
    hd_tm = analysis['homodimer_tm']
    hd_dg = analysis['homodimer_dg']

    if hd_tm >= primer_tm:
        report.write(f"  ❌ Homodimer Tm={hd_tm:.2f}, ΔG={hd_dg:.2f}  (СТРУКТУРА СТАБИЛЬНЕЕ ПРАЙМЕРА)\n")
    elif hd_tm > limit:
        report.write(f"  ❌ Homodimer Tm={hd_tm:.2f}, ΔG={hd_dg:.2f}  (выше лимита {limit:.2f})\n")
    else:
        report.write(f"  ✔ Homodimer OK  Tm={hd_tm:.2f}, ΔG={hd_dg:.2f}\n")

    report.write(f"  Binding sites on target: {len(binding_sites)}\n")
    report.write("\n")

    direction = infer_direction_from_bindings(binding_sites)

    primer_db[name] = {
        'direction': direction,
        'iupac_seq': seq_iupac,
        'dna_variants': dna_variants,
        'tm': primer_tm,
        'hairpin_tm': hp_tm,
        'hairpin_dg': hp_dg,
        'homodimer_tm': hd_tm,
        'homodimer_dg': hd_dg,
        'binding_sites': binding_sites
    }

print("[STEP 1] DONE\n")

# Экспорт посадок праймеров в CSV
bindings_rows = []
for name, pdata in primer_db.items():
    for start, end, strand, var, total_mm, mm_3p in pdata["binding_sites"]:
        # Публичные координаты (CSV/UI): 1-based, inclusive
        start_1 = int(start) + 1
        end_1 = int(end)  # end is exclusive (0-based), so inclusive 1-based equals end
        bindings_rows.append({
            "Primer": name,
            "Variant": var,
            "Start": start_1,
            "End": end_1,
            "Strand": strand,
            "Total_mismatches": total_mm,
            "Mismatches_3prime": mm_3p,
        })

pd.DataFrame(bindings_rows).to_csv("primer_binding_sites.csv", index=False)

# BLAST-посадки праймеров на референс (лучший хит)
blast_rows = []
blast_enabled = os.getenv("BLAST_ENABLED", "1") == "1"
blast_human_enabled = os.getenv("BLAST_HUMAN_ENABLED", "0") == "1"

if blast_enabled:
    for name, pdata in primer_db.items():
        hit = blast_primer_on_reference(pdata['dna_variants'][0], TARGET_GB, template)
        if hit is None:
            blast_rows.append({
                "Primer": name,
                "Strand": "",
                "Q_start": "",
                "Q_end": "",
                "S_start": "",
                "S_end": "",
                "Total_mismatches": "",
                "Mismatches_3prime": "",
            })
        else:
            blast_rows.append({
                "Primer": name,
                "Strand": hit["strand"],
                "Q_start": hit["qstart"],
                "Q_end": hit["qend"],
                "S_start": hit["sstart"],
                "S_end": hit["send"],
                "Total_mismatches": hit["total_mismatches"],
                "Mismatches_3prime": hit["mismatches_3prime"],
            })

    pd.DataFrame(blast_rows).to_csv("primer_blast_ref_hits.csv", index=False)

# Отчёт по гомодимерам (лабораторный стандарт)
homo_qc_rows = []
for name, pdata in primer_db.items():
    tm_primer = pdata["tm"]
    hp_tm = pdata["hairpin_tm"]
    hp_dg = pdata["hairpin_dg"]
    homo_tm = pdata["homodimer_tm"]
    homo_dg = pdata["homodimer_dg"]
    limit = min(tm_primer, ANNEALING_TEMP) - 15.0

    # Hairpin QC
    if hp_tm >= tm_primer:
        hp_status = "BAD"
        hp_reason = "СТРУКТУРА СТАБИЛЬНЕЕ ПРАЙМЕРА"
    elif hp_tm > limit:
        hp_status = "BAD"
        hp_reason = f"выше лимита {limit:.2f}"
    else:
        hp_status = "OK"
        hp_reason = "OK"

    if homo_tm >= tm_primer:
        status = "BAD"
        reason = "СТРУКТУРА СТАБИЛЬНЕЕ ПРАЙМЕРА"
    elif homo_tm > limit:
        status = "BAD"
        reason = f"выше лимита {limit:.2f}"
    else:
        status = "OK"
        reason = "OK"

    homo_qc_rows.append({
        "Primer": name,
        "Tm_primer": tm_primer,
        "Hairpin_Tm": hp_tm,
        "Hairpin_dG": hp_dg,
        "Hairpin_Status": hp_status,
        "Hairpin_Reason": hp_reason,
        "Homodimer_Tm": homo_tm,
        "Homodimer_dG": homo_dg,
        "Status": status,
        "Reason": reason,
    })

pd.DataFrame(homo_qc_rows).to_csv("primer_homodimer_qc.csv", index=False)

# Экспорт полной матрицы димеров (гомо- и гетеродимеры) + QC-статус
dimer_rows = []
primer_names = list(primer_db.keys())

def dimer_status(tm_dimer, tm_p1, tm_p2, ta):
    """
    Критерий стабильности:
      Tm_struct <= min(Tm_primer1, Tm_primer2, Ta) - 15  → нестабильный (OK)
      иначе → стабильный (BAD)
    """
    try:
        limit = min(float(tm_p1), float(tm_p2), float(ta)) - 15.0
        if tm_dimer is None:
            return "OK"
        if float(tm_dimer) <= limit:
            return "OK"
        return "BAD"
    except Exception:
        return "UNKNOWN"

for i, name_i in enumerate(primer_names):
    seq_i = primer_db[name_i]['dna_variants'][0]
    tm_i = primer_db[name_i]['tm']

    # гомодимер (i,i)
    homo_tm = primer_db[name_i]['homodimer_tm']
    homo_dg = primer_db[name_i]['homodimer_dg']
    status_homo = dimer_status(homo_tm, tm_i, tm_i, ANNEALING_TEMP)
    dimer_rows.append({
        "Primer1": name_i,
        "Primer2": name_i,
        "Type": "homodimer",
        "Tm": homo_tm,
        "dG": homo_dg,
        "Tm_P1": tm_i,
        "Tm_P2": tm_i,
        "Status": status_homo,
    })

    # гетеродимеры (i,j), j > i
    for name_j in primer_names[i+1:]:
        seq_j = primer_db[name_j]['dna_variants'][0]
        tm_j = primer_db[name_j]['tm']
        hd = heterodimer_analysis(seq_i, seq_j)
        tm_hd = hd['heterodimer_tm']
        status_het = dimer_status(tm_hd, tm_i, tm_j, ANNEALING_TEMP)
        dimer_rows.append({
            "Primer1": name_i,
            "Primer2": name_j,
            "Type": "heterodimer",
            "Tm": tm_hd,
            "dG": hd['heterodimer_dg'],
            "Tm_P1": tm_i,
            "Tm_P2": tm_j,
            "Status": status_het,
        })

pd.DataFrame(dimer_rows).to_csv("primer_dimers_matrix.csv", index=False)

# ---------- ПАРЫ ПРАЙМЕРОВ ----------

print("[STEP 2] Анализ пар праймеров")

results = []

f_list = [k for k, v in primer_db.items() if v['direction'] == 'Forward']
r_list = [k for k, v in primer_db.items() if v['direction'] == 'Reverse']

pair_total = len(f_list) * len(r_list)
pair_count = 0

for f_name in f_list:
    for r_name in r_list:
        pair_count += 1
        print(f"  [{pair_count}/{pair_total}] Pair {f_name} + {r_name}")

        f_data = primer_db[f_name]
        r_data = primer_db[r_name]

        tm_diff = abs(f_data['tm'] - r_data['tm'])
        hetero = heterodimer_analysis(f_data['dna_variants'][0], r_data['dna_variants'][0])
        matches_3p = check_3prime_complementarity(f_data['dna_variants'][0], r_data['dna_variants'][0])

        # Для анализа продуктов используем "лучшие" посадки праймеров:
        # 1) минимальные mismatches на 3'-конце
        # 2) затем минимальные mismatches по всему праймеру
        # (для детерминизма добавляем координаты как третий критерий)
        f_candidates = [h for h in f_data['binding_sites'] if h[2] == '+']
        r_candidates = [h for h in r_data['binding_sites'] if h[2] == '-']

        if not f_candidates or not r_candidates:
            continue

        f_hit = min(f_candidates, key=lambda h: (h[5], h[4], h[0], h[1]))
        r_hit = min(r_candidates, key=lambda h: (h[5], h[4], h[0], h[1]))

        # Проверяем ориентацию и размер продукта
        if not (r_hit[0] > f_hit[1]):
            continue

        size = r_hit[0] - f_hit[1]
        if not (MIN_AMPLICON <= size <= MAX_AMPLICON):
            continue

        results.append({
            'Forward_primer': f_name,
            'Reverse_primer': r_name,
            # Публичные координаты (CSV/UI): 1-based, inclusive
            'Amplicon_start': int(f_hit[0]) + 1,
            'Amplicon_end': int(r_hit[1]),
            # Внутренние координаты для GenBank (0-based, end-exclusive)
            'Amplicon_start0': int(f_hit[0]),
            'Amplicon_end0': int(r_hit[1]),
            'Product_size': size,
            'F_total_mismatches': f_hit[4],
            'F_3prime_mismatches': f_hit[5],
            'R_total_mismatches': r_hit[4],
            'R_3prime_mismatches': r_hit[5],
            'Tm_F': f_data['tm'],
            'Tm_R': r_data['tm'],
            'Tm_diff': round(tm_diff, 1),
            'Heterodimer_TM': hetero['heterodimer_tm'],
            'Heterodimer_dG': hetero['heterodimer_dg'],
            '3prime_matches': matches_3p,
            'Forward_sequence': f_data['iupac_seq'],
            'Reverse_sequence': r_data['iupac_seq'],
        })

print(f"[STEP 2] DONE — найдено пар с ампликонами: {len(results)}\n")

# ---------- BLAST ----------

print("[STEP 3] BLAST праймеров против человека")

human_hits = {}

if blast_human_enabled:
    for i, (pname, pdata) in enumerate(primer_db.items()):
        print(f"  [{i+1}/{len(primer_db)}] BLAST {pname}")
        human_hits[pname] = blast_primer(pdata['dna_variants'][0])
    print("[STEP 3] DONE\n")

    # Сводка по BLAST праймеров против человека (по каждому праймеру)
    human_blast_rows = []
    for pname, hits in human_hits.items():
        if not hits:
            human_blast_rows.append({
                "Primer": pname,
                "Hits_count": 0,
                "Best_mismatch": "",
                "Best_length": "",
            })
        else:
            best = min(hits, key=lambda h: (h['mismatch'], -h['length']))
            human_blast_rows.append({
                "Primer": pname,
                "Hits_count": len(hits),
                "Best_mismatch": best['mismatch'],
                "Best_length": best['length'],
            })
        print(f"[BLAST human summary] {pname}: {human_blast_rows[-1]['Hits_count']} hits, "
              f"best mismatch={human_blast_rows[-1]['Best_mismatch']}, "
              f"best length={human_blast_rows[-1]['Best_length']}")
    pd.DataFrame(human_blast_rows).to_csv("primer_human_blast_summary.csv", index=False)
else:
    print("[STEP 3] Пропущен (BLAST по человеку отключен)\n")

# ---------- QC НА ЧЕЛОВЕКЕ ----------

print("[STEP 4] QC неспецифических ампликонов на человеке")

if blast_human_enabled:

    for idx, row in enumerate(results):

        print(f"  [{idx+1}/{len(results)}] QC {row['Forward_primer']} + {row['Reverse_primer']}")

        fname = row['Forward_primer']
        rname = row['Reverse_primer']
        L_ref = row['Product_size']

        reason = "specific"
        best_L_human = None

        for hf in human_hits.get(fname, []):
            for hr in human_hits.get(rname, []):
                if hf['chrom'] != hr['chrom']:
                    continue

                L_human = hr['sstart'] - hf['send']
                if L_human <= 0:
                    continue

                delta = L_human - L_ref

                # запоминаем ближайший по длине человеческий продукт
                if best_L_human is None or L_human < best_L_human:
                    best_L_human = L_human

                if L_human <= L_ref:
                    reason = f"human product {L_human} <= target {L_ref}"
                elif delta <= MAX_HUMAN_DELTA:
                    reason = f"human delta {delta} <= {MAX_HUMAN_DELTA}"

        # текстовый статус (бывший human_reason)
        row['human_status'] = reason
        # компактный BLAST-результат по человеку для пары
        if best_L_human is None:
            row['Hum_BLAST'] = "not product"
        else:
            row['Hum_BLAST'] = best_L_human

    print("[STEP 4] DONE\n")
else:
    print("[STEP 4] Пропущен (BLAST по человеку отключен)\n")
    for row in results:
        row['human_status'] = "BLAST human disabled"
        row['Hum_BLAST'] = "not product"

# ---------- GENBANK ----------

print("[STEP 5] Аннотирование GenBank")

# удаляем старые продукты
record.features = [f for f in record.features if f.type != "PCR_product"]

# ----------------------------------------------------------
# ДОБАВЛЕНИЕ АННОТАЦИЙ ПРАЙМЕРОВ (ЕСЛИ ИХ НЕТ В ФАЙЛЕ)
# Формат UGENE:
#   primer 2010..2032
#     /mismatches=0
#     /note="POM; Forward primer"
#     /ugene_name="POM"
#     /ugene_group="Primers"
#   primer complement(2611..2637)
# ----------------------------------------------------------

existing_primers = set()

for f in record.features:
    if f.type == "primer" and "ugene_name" in f.qualifiers:
        existing_primers.add(f.qualifiers["ugene_name"][0])

for pname, pdata in primer_db.items():

    # если уже аннотирован — пропускаем
    if pname in existing_primers:
        continue

    for hit in pdata['binding_sites']:
        # hit: (start, end, strand, used_variant, total_mismatches, mismatches_3prime)
        start, end, strand, seq, total_mm, mm_3p = hit

        qualifiers = {
            "mismatches": "0",
            "ugene_name": pname,
            "ugene_group": "Primers",
            "note": f"{pname}; {pdata['direction']} primer"
        }

        # Forward primer
        if strand == '+':
            feature = SeqFeature(
                FeatureLocation(start, end, strand=+1),
                type="primer",
                qualifiers=qualifiers
            )

        # Reverse primer (complement)
        else:
            feature = SeqFeature(
                FeatureLocation(start, end, strand=-1),
                type="primer",
                qualifiers=qualifiers
            )

        record.features.append(feature)

# ----------------------------------------------------------
# ДОБАВЛЕНИЕ PCR ПРОДУКТОВ
# ----------------------------------------------------------

print("[STEP 5] Аннотирование GenBank")

record.features = [f for f in record.features if f.type != "PCR_product"]

for row in results:
    qualifiers = {
        "forward_primer": row["Forward_primer"],
        "reverse_primer": row["Reverse_primer"],
        "product_size": str(row["Product_size"]),
        "Tm_forward": f"{row['Tm_F']:.1f}",
        "Tm_reverse": f"{row['Tm_R']:.1f}",
        "Tm_diff": f"{row['Tm_diff']:.1f}",
        "heterodimer_TM": f"{row['Heterodimer_TM']:.1f}",
        "heterodimer_dG": f"{row['Heterodimer_dG']:.2f}",
        "specificity": row.get("human_status", "NA"),
        "specificity_note": row.get("human_reason", "")
    }

    feature = SeqFeature(
        FeatureLocation(int(row.get("Amplicon_start0", row["Amplicon_start"])), int(row.get("Amplicon_end0", row["Amplicon_end"]))),
        type="PCR_product",
        qualifiers=qualifiers
    )

    record.features.append(feature)

SeqIO.write(record, OUT_GB, "genbank")

pd.DataFrame(results).to_csv(OUT_CSV, index=False)

report.close()

print("="*60)
print("ГОТОВО")
print("="*60)
print(f"CSV:     {OUT_CSV}")
print(f"GenBank: {OUT_GB}")
print(f"Отчёт:   {OUT_REPORT}")
