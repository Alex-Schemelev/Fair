IUPAC_FROM_BASES = {
    frozenset({"A"}): "A",
    frozenset({"C"}): "C",
    frozenset({"G"}): "G",
    frozenset({"T"}): "T",
    frozenset({"A", "G"}): "R",
    frozenset({"C", "T"}): "Y",
    frozenset({"G", "C"}): "S",
    frozenset({"A", "T"}): "W",
    frozenset({"G", "T"}): "K",
    frozenset({"A", "C"}): "M",
    frozenset({"C", "G", "T"}): "B",
    frozenset({"A", "G", "T"}): "D",
    frozenset({"A", "C", "T"}): "H",
    frozenset({"A", "C", "G"}): "V",
    frozenset({"A", "C", "G", "T"}): "N",
}


def consensus_iupac_from_row(row, min_freq=0.05):
    bases = set()
    for b in ("A", "C", "G", "T"):
        try:
            if float(row.get(b, 0)) >= min_freq:
                bases.add(b)
        except Exception:
            pass
    # gap/пробел хранится в столбце '-' (иногда в таблицах консенсуса)
    gap_freq = 0.0
    try:
        gap_freq = float(row.get("-", 0) or 0)
    except Exception:
        gap_freq = 0.0

    # Пробел (gap) учитываем только если он >= min_freq.
    # Если gap значим и при этом нет ни одной базы >= min_freq — возвращаем пробел.
    if gap_freq >= min_freq and not bases:
        return " "
    if not bases:
        return "N"
    return IUPAC_FROM_BASES.get(frozenset(bases), "N")


IUPAC_COMP = {
    "A": "T", "T": "A", "C": "G", "G": "C",
    "R": "Y", "Y": "R", "S": "S", "W": "W",
    "K": "M", "M": "K", "B": "V", "V": "B",
    "D": "H", "H": "D", "N": "N", " ": " ",
}

IUPAC_SET = {
    "A": {"A"},
    "C": {"C"},
    "G": {"G"},
    "T": {"T"},
    "R": {"A", "G"},
    "Y": {"C", "T"},
    "S": {"G", "C"},
    "W": {"A", "T"},
    "K": {"G", "T"},
    "M": {"A", "C"},
    "B": {"C", "G", "T"},
    "D": {"A", "G", "T"},
    "H": {"A", "C", "T"},
    "V": {"A", "C", "G"},
    "N": {"A", "C", "G", "T"},
    " ": set(),
    "-": set(),
}


def iupac_complement(seq: str) -> str:
    return "".join(IUPAC_COMP.get(ch, "N") for ch in str(seq).upper())


def _iupac_match(a: str, b: str) -> bool:
    sa = IUPAC_SET.get(str(a).upper(), set())
    sb = IUPAC_SET.get(str(b).upper(), set())
    return bool(sa & sb)


def mismatch_stats_vs_consensus(primer_aln: str, consensus_aln: str, three_prime_window: int = 5):
    """
    primer_aln / consensus_aln: ровно те строки, которые показываем в UI (в т.ч. с '-' для инделов).
    Возвращает mismatches/matches и счётчики инделов (праймер или консенсус = '-').
    Для 3'-окна: последние three_prime_window БАЗ праймера (без '-').
    """
    p = list(str(primer_aln or ""))
    c = list(str(consensus_aln or ""))
    L = min(len(p), len(c))

    per_base_match = []
    total_matches = 0
    total_mismatches = 0

    # индекс базы праймера для каждой позиции выравнивания (None = гэп в праймере)
    primer_base_idx = []
    bidx = 0
    for i in range(L):
        if p[i] == "-":
            primer_base_idx.append(None)
        else:
            primer_base_idx.append(bidx)
            bidx += 1
    n_bases = bidx
    w = int(three_prime_window or 0)
    three_start = max(0, n_bases - w) if w > 0 else n_bases

    indels_total = 0
    indels_3p = 0

    for i in range(L):
        pb = p[i]
        cb = c[i]
        is_indel = (pb == "-") or (cb == "-")
        if is_indel:
            indels_total += 1
            bi = primer_base_idx[i]
            if bi is None:
                # инсерция в праймере: относим к 3' если соседняя база справа (или слева) в 3'-окне
                right = None
                for j in range(i + 1, L):
                    if primer_base_idx[j] is not None:
                        right = primer_base_idx[j]
                        break
                left = None
                for j in range(i - 1, -1, -1):
                    if primer_base_idx[j] is not None:
                        left = primer_base_idx[j]
                        break
                ref = right if right is not None else left
                if ref is not None and ref >= three_start:
                    indels_3p += 1
            else:
                if bi >= three_start:
                    indels_3p += 1
            continue  # инделы не считаем как match/mismatch баз

        ok = _iupac_match(pb, cb)
        per_base_match.append(ok)
        if ok:
            total_matches += 1
        else:
            total_mismatches += 1

    if len(p) != len(c):
        total_mismatches += abs(len(p) - len(c))

    if w <= 0:
        mismatches_3p = 0
        indels_3p = 0
    else:
        tail = per_base_match[-w:] if len(per_base_match) >= w else per_base_match
        mismatches_3p = sum(1 for x in tail if not x)

    return {
        "Matches": int(total_matches),
        "Mismatches_total": int(total_mismatches),
        "Mismatches_3prime": int(mismatches_3p),
        "Indels_total": int(indels_total),
        "Indels_3prime": int(indels_3p),
        "Has_indels": bool(indels_total > 0),
    }
