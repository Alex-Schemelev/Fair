import hashlib


def stable_product_id(forward_primer, reverse_primer, start, end):
    """Стабильный ID продукта: не зависит от порядка фильтрации/индекса в списке."""
    raw = f"{forward_primer}|{reverse_primer}|{int(start)}|{int(end)}"
    digest = hashlib.sha1(raw.encode("utf-8")).hexdigest()[:12]
    return f"prod-{digest}"


def stable_primer_feature_id(name, start, end, strand="."):
    raw = f"{name}|{int(start)}|{int(end)}|{strand or '.'}"
    digest = hashlib.sha1(raw.encode("utf-8")).hexdigest()[:12]
    return f"primer-{digest}"


def rebuild_annotation_products(annotation_map, products):
    """Оставляет ref/primer фичи, продукты пересобирает из отфильтрованного списка."""
    am = annotation_map or {"length": 0, "id": "", "features": []}
    feats = [f for f in (am.get("features") or []) if f.get("category") != "product"]
    for p in products or []:
        try:
            a_start = int(p.get("Amplicon_start"))
            a_end = int(p.get("Amplicon_end"))
        except Exception:
            continue
        fp = p.get("Forward_primer", "")
        rp = p.get("Reverse_primer", "")
        name = f"{fp} + {rp}".strip(" +")
        pid = p.get("product_id") or stable_product_id(fp, rp, a_start, a_end)
        p["product_id"] = pid
        feats.append({
            "id": pid,
            "type": "PCR_product",
            "category": "product",
            "name": name,
            "start": a_start,
            "end": a_end,
            "strand": ".",
            "info": {
                "type": "PCR_product",
                "name": name,
                "start": a_start,
                "end": a_end,
                "strand": ".",
                "length": a_end - a_start + 1,
                "product_size": p.get("Product_size"),
                "Tm_forward": p.get("Tm_F"),
                "Tm_reverse": p.get("Tm_R"),
                "Tm_diff": p.get("Tm_diff"),
                "specificity": p.get("human_status"),
                "forward_primer": fp,
                "reverse_primer": rp,
                "product_id": pid,
            },
        })
    out = dict(am)
    out["features"] = feats
    return out


def extract_annotation_map(gb_path, products=None, bindings=None):
    """
    Собирает данные для интерактивной схемы аннотаций:
    - фичи из референсного/аннотированного GenBank
    - плюс (если нужно) продукты/посадки из результатов пайплайна

    Координаты в ответе: 1-based inclusive.
    """
    from Bio import SeqIO

    record = next(SeqIO.parse(gb_path, "genbank"))
    length = len(record.seq)
    features = []
    idx = 0

    def qget(qualifiers, *keys, default=""):
        for k in keys:
            if k in qualifiers and qualifiers[k]:
                v = qualifiers[k][0] if isinstance(qualifiers[k], list) else qualifiers[k]
                return str(v).replace("\\", "")
        return default

    def category_for(ftype):
        ft = (ftype or "").lower()
        if ft in ("primer",):
            return "primer"
        if ft in ("pcr_product", "amplicon"):
            return "product"
        if ft in ("source",):
            return "skip"
        return "ref"

    for feat in record.features:
        cat = category_for(feat.type)
        if cat == "skip":
            continue
        try:
            start0 = int(feat.location.start)
            end0 = int(feat.location.end)
        except Exception:
            continue
        start = start0 + 1
        end = end0
        if end < start:
            start, end = end, start

        strand = feat.location.strand
        if strand == 1:
            strand_s = "+"
        elif strand == -1:
            strand_s = "-"
        else:
            strand_s = "."

        name = qget(feat.qualifiers, "ugene_name", "gene", "product", "label", "note", default=feat.type)
        if cat == "product":
            fp = qget(feat.qualifiers, "forward_primer")
            rp = qget(feat.qualifiers, "reverse_primer")
            if fp or rp:
                name = f"{fp} + {rp}".strip(" +")

        info = {
            "type": feat.type,
            "name": name,
            "start": start,
            "end": end,
            "strand": strand_s,
            "length": end - start + 1,
        }
        for key in (
            "note", "product_size", "Tm_forward", "Tm_reverse", "Tm_diff",
            "heterodimer_TM", "heterodimer_dG", "specificity", "mismatches",
            "forward_primer", "reverse_primer", "ugene_group",
        ):
            if key in feat.qualifiers and feat.qualifiers[key]:
                val = feat.qualifiers[key][0] if isinstance(feat.qualifiers[key], list) else feat.qualifiers[key]
                info[key] = str(val).replace("\\", "")

        features.append({
            "id": f"gb-{idx}",
            "type": feat.type,
            "category": cat,
            "name": name,
            "start": start,
            "end": end,
            "strand": strand_s,
            "info": info,
        })
        idx += 1

    has_products = any(f["category"] == "product" for f in features)
    has_primers = any(f["category"] == "primer" for f in features)

    if (not has_products) and products:
        for p in products:
            try:
                a_start = int(p.get("Amplicon_start"))
                a_end = int(p.get("Amplicon_end"))
            except Exception:
                continue
            fp = p.get("Forward_primer", "")
            rp = p.get("Reverse_primer", "")
            name = f"{fp} + {rp}".strip(" +")
            pid = p.get("product_id") or stable_product_id(fp, rp, a_start, a_end)
            p["product_id"] = pid
            info = {
                "type": "PCR_product",
                "name": name,
                "start": a_start,
                "end": a_end,
                "strand": ".",
                "length": a_end - a_start + 1,
                "product_size": p.get("Product_size"),
                "Tm_forward": p.get("Tm_F"),
                "Tm_reverse": p.get("Tm_R"),
                "Tm_diff": p.get("Tm_diff"),
                "specificity": p.get("human_status"),
                "forward_primer": fp,
                "reverse_primer": rp,
                "product_id": pid,
            }
            features.append({
                "id": pid,
                "type": "PCR_product",
                "category": "product",
                "name": name,
                "start": a_start,
                "end": a_end,
                "strand": ".",
                "info": info,
            })

    if (not has_primers) and bindings:
        for b in bindings:
            try:
                start = int(b.get("Start"))
                end = int(b.get("End"))
            except Exception:
                continue
            name = b.get("Primer", "primer")
            strand_s = b.get("Strand", ".") or "."
            pid = stable_primer_feature_id(name, start, end, strand_s)
            info = {
                "type": "primer",
                "name": name,
                "start": start,
                "end": end,
                "strand": strand_s,
                "length": end - start + 1,
                "Total_mismatches": b.get("Total_mismatches"),
                "Mismatches_3prime": b.get("Mismatches_3prime"),
            }
            features.append({
                "id": pid,
                "type": "primer",
                "category": "primer",
                "name": name,
                "start": start,
                "end": end,
                "strand": strand_s,
                "info": info,
            })

    return {
        "length": length,
        "id": getattr(record, "id", "reference"),
        "features": features,
    }
