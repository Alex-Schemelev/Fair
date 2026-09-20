from .iupac import mismatch_stats_vs_consensus

# Defaults synchronized with web_app / In_silico_PCR.py
ANNEALING_TEMP = 60.0
MAX_3PRIME_MATCHES = 3
THREE_PRIME_WINDOW = 5


def apply_product_filters(products, conservation, homo_qc, dimers_matrix, params, filters):
    """
    Фильтрация продуктов по опциональным правилам из UI.
    Пустые числовые пороги => соответствующий фильтр не применяется.
    """
    filters = filters or {}
    params = params or {}

    def _as_float_or_none(v):
        if v is None or v == "":
            return None
        try:
            return float(v)
        except Exception:
            return None

    def _as_bool(v):
        return bool(v)

    min_s = _as_float_or_none(filters.get("min_s_index_avg", params.get("MIN_S_INDEX_AVG")))
    max_tm_diff = _as_float_or_none(filters.get("max_tm_diff", params.get("MAX_TM_DIFF")))
    filter_dimer = _as_bool(filters.get("filter_dimer_above_ta", False))
    filter_3prime = _as_bool(filters.get("filter_3prime", False))

    def _opt_int_from(src, key):
        v = src.get(key, "")
        if v is None or v == "":
            return None
        try:
            return int(v)
        except Exception:
            return None

    max_indels = _opt_int_from(filters, "max_indels_total")
    if max_indels is None:
        max_indels = _opt_int_from(params, "MAX_INDELS_TOTAL")
    max_indels_3p = _opt_int_from(filters, "max_indels_3prime")
    if max_indels_3p is None:
        max_indels_3p = _opt_int_from(params, "MAX_INDELS_3PRIME")

    ta = _as_float_or_none(params.get("ANNEALING_TEMP"))
    if ta is None:
        ta = ANNEALING_TEMP
    max_3p_matches = None
    try:
        max_3p_matches = int(params.get("MAX_3PRIME_MATCHES", MAX_3PRIME_MATCHES))
    except Exception:
        max_3p_matches = MAX_3PRIME_MATCHES

    cons_by_primer = {r.get("Primer"): r for r in (conservation or []) if r.get("Primer")}
    homo_by_primer = {r.get("Primer"): r for r in (homo_qc or []) if r.get("Primer")}
    win3 = int(params.get("THREE_PRIME_WINDOW", THREE_PRIME_WINDOW) or THREE_PRIME_WINDOW)

    # быстрый индекс гетеродимеров по паре имён
    hetero_tm = {}
    for row in dimers_matrix or []:
        if str(row.get("Type", "")).lower() != "heterodimer":
            continue
        a = row.get("Primer1")
        b = row.get("Primer2")
        if not a or not b:
            continue
        try:
            tm = float(row.get("Tm", 0) or 0)
        except Exception:
            tm = 0.0
        hetero_tm[(a, b)] = tm
        hetero_tm[(b, a)] = tm

    kept = []
    for p in products or []:
        fp = p.get("Forward_primer")
        rp = p.get("Reverse_primer")

        # 1) max ΔTm
        if max_tm_diff is not None:
            try:
                if float(p.get("Tm_diff", 0) or 0) > max_tm_diff:
                    continue
            except Exception:
                continue

        # 2) min avg S_index (оба праймера)
        if min_s is not None:
            ok_s = True
            for name in (fp, rp):
                r = cons_by_primer.get(name)
                if not r or r.get("S_index_avg") is None:
                    ok_s = False
                    break
                try:
                    if float(r.get("S_index_avg")) < min_s:
                        ok_s = False
                        break
                except Exception:
                    ok_s = False
                    break
            if not ok_s:
                continue

        # 3) димеры с Tm > Ta (гомо для F/R и гетеро для пары)
        if filter_dimer:
            bad = False
            for name in (fp, rp):
                hq = homo_by_primer.get(name) or {}
                try:
                    if float(hq.get("Homodimer_Tm", 0) or 0) > ta:
                        bad = True
                        break
                except Exception:
                    pass
            if not bad:
                # сначала колонка продукта, затем матрица
                try:
                    het = float(p.get("Heterodimer_TM", 0) or 0)
                except Exception:
                    het = 0.0
                if (fp, rp) in hetero_tm:
                    het = max(het, float(hetero_tm.get((fp, rp), 0) or 0))
                if het > ta:
                    bad = True
            if bad:
                continue

        # 4) инделы: пороги MAX_INDELS_TOTAL / MAX_INDELS_3PRIME (пусто = не фильтровать)
        if max_indels is not None or max_indels_3p is not None:
            bad_indel = False
            for name in (fp, rp):
                r = cons_by_primer.get(name)
                if not r:
                    bad_indel = True
                    break
                stats = mismatch_stats_vs_consensus(
                    r.get("Primer_seq", ""),
                    r.get("Consensus_seq", ""),
                    three_prime_window=win3,
                )
                if max_indels is not None and stats["Indels_total"] > max_indels:
                    bad_indel = True
                    break
                if max_indels_3p is not None and stats["Indels_3prime"] > max_indels_3p:
                    bad_indel = True
                    break
            if bad_indel:
                continue

        # 5) 3'-комплементарность между праймерами
        if filter_3prime and max_3p_matches is not None:
            try:
                if int(p.get("3prime_matches", 0) or 0) > max_3p_matches:
                    continue
            except Exception:
                continue

        kept.append(p)

    return kept
