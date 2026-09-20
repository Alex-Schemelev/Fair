/* global React */
(function () {
  const e = React.createElement;

  const TYPE_COLORS = {
    PCR_product: "#0ea5e9",
    primer: "#f59e0b",
    Protein: "#8b5cf6",
    misc_feature: "#10b981",
    gene: "#3b82f6",
    CDS: "#6366f1",
    "5'UTR": "#f43f5e",
    "3'UTR": "#ec4899",
    source: "#94a3b8",
  };

  const CATEGORY_ORDER = ["ref", "primer", "product"];
  const CATEGORY_LABELS = {
    ref: "Аннотации референса",
    primer: "Праймеры",
    product: "ПЦР‑продукты",
  };

  function colorFor(type) {
    return TYPE_COLORS[type] || "#64748b";
  }

  function packLanes(items) {
    const sorted = items.slice().sort((a, b) => {
      if (a.start !== b.start) return a.start - b.start;
      return a.end - b.end;
    });
    const laneEnds = [];
    const packed = [];
    sorted.forEach((item) => {
      let lane = -1;
      for (let i = 0; i < laneEnds.length; i++) {
        if (item.start > laneEnds[i]) {
          lane = i;
          break;
        }
      }
      if (lane < 0) {
        lane = laneEnds.length;
        laneEnds.push(item.end);
      } else {
        laneEnds[lane] = item.end;
      }
      packed.push(Object.assign({}, item, { lane: lane }));
    });
    return { packed: packed, laneCount: Math.max(1, laneEnds.length) };
  }

  function genomicOverlap(a, b) {
    return a.start <= b.end && b.start <= a.end;
  }

  /**
   * Перенос блока на целевую дорожку: координаты (start/end) не меняются.
   * Перетаскиваемый блок занимает targetLane; пересекающиеся по геному блоки
   * с этой и нижележащих дорожек сдвигаются вниз, пока наложений не останется.
   */
  function moveItemToLane(items, laneById, itemId, targetLane) {
    const auto = packLanes(items);
    const lanes = {};
    auto.packed.forEach((p) => {
      lanes[p.id] = p.lane;
    });
    if (laneById) {
      Object.keys(laneById).forEach((id) => {
        if (laneById[id] !== undefined && laneById[id] !== null) {
          lanes[id] = Math.max(0, laneById[id] | 0);
        }
      });
    }
    lanes[itemId] = Math.max(0, targetLane | 0);

    let guard = 0;
    while (guard++ < 2000) {
      let changed = false;
      for (let i = 0; i < items.length; i++) {
        const a = items[i];
        if (a.id === itemId) continue;
        for (let j = 0; j < items.length; j++) {
          const b = items[j];
          if (a.id === b.id) continue;
          if (lanes[a.id] !== lanes[b.id]) continue;
          if (!genomicOverlap(a, b)) continue;
          // конфликт: сдвигаем вниз любого, кроме перетаскиваемого
          lanes[a.id] = lanes[a.id] + 1;
          changed = true;
          break;
        }
        if (changed) break;
      }
      if (!changed) break;
    }
    return lanes;
  }

  function packWithLaneMap(items, laneById) {
    if (!items || !items.length) return { packed: [], laneCount: 1 };
    if (!laneById || !Object.keys(laneById).length) {
      return packLanes(items);
    }

    const allCovered = items.every(
      (it) => laneById[it.id] !== undefined && laneById[it.id] !== null
    );
    if (allCovered) {
      let maxLane = 0;
      const packed = items.map((it) => {
        const lane = Math.max(0, laneById[it.id] | 0);
        if (lane > maxLane) maxLane = lane;
        return Object.assign({}, it, { lane: lane });
      });
      return { packed: packed, laneCount: Math.max(1, maxLane + 1) };
    }

    const auto = packLanes(items);
    const lanes = {};
    auto.packed.forEach((p) => {
      lanes[p.id] = p.lane;
    });
    Object.keys(laneById).forEach((id) => {
      if (laneById[id] !== undefined && laneById[id] !== null) {
        lanes[id] = Math.max(0, laneById[id] | 0);
      }
    });

    const pinnedIds = new Set(
      Object.keys(laneById).filter((id) => laneById[id] !== undefined && laneById[id] !== null)
    );
    let guard = 0;
    while (guard++ < 2000) {
      let changed = false;
      for (let i = 0; i < items.length; i++) {
        const a = items[i];
        for (let j = i + 1; j < items.length; j++) {
          const b = items[j];
          if (lanes[a.id] !== lanes[b.id]) continue;
          if (!genomicOverlap(a, b)) continue;
          const aPinned = pinnedIds.has(a.id);
          const bPinned = pinnedIds.has(b.id);
          let victim = null;
          if (aPinned && !bPinned) victim = b;
          else if (!aPinned && bPinned) victim = a;
          else victim = a.start >= b.start ? a : b;
          lanes[victim.id] = lanes[victim.id] + 1;
          changed = true;
        }
      }
      if (!changed) break;
    }

    let maxLane = 0;
    const packed = items.map((it) => {
      const lane = lanes[it.id] !== undefined ? lanes[it.id] : 0;
      if (lane > maxLane) maxLane = lane;
      return Object.assign({}, it, { lane: lane });
    });
    return { packed: packed, laneCount: Math.max(1, maxLane + 1) };
  }

  /**
   * Упаковка с приоритетом порядка: сначала пытаемся положить элемент в верхние дорожки,
   * не перекрывая уже лежащие. Так более приоритетные (идущие раньше в items) оказываются выше.
   */
  function packLanesByPriority(items) {
    const lanes = []; // lane -> last end
    const packed = [];
    items.forEach((item) => {
      let lane = -1;
      for (let i = 0; i < lanes.length; i++) {
        if (item.start > lanes[i]) {
          lane = i;
          break;
        }
      }
      if (lane < 0) {
        lane = lanes.length;
        lanes.push(item.end);
      } else {
        lanes[lane] = Math.max(lanes[lane], item.end);
      }
      packed.push(Object.assign({}, item, { lane: lane }));
    });
    return { packed: packed, laneCount: Math.max(1, lanes.length) };
  }

  /**
   * Жёсткая визуальная сортировка: одна дорожка = один продукт, сверху самые подходящие.
   */
  function packLanesStrictRank(items) {
    const packed = items.map((item, idx) =>
      Object.assign({}, item, { lane: idx, visualRank: idx + 1 })
    );
    return { packed: packed, laneCount: Math.max(1, items.length) };
  }

  const SCHEME_STATUS = {
    untested: {
      id: "untested",
      label: "не тестировалась",
      fill: "rgba(148,163,184,0.28)",
      stroke: "#64748b",
      chipBg: "#f1f5f9",
      chipFg: "#334155",
      chipBorder: "#94a3b8",
    },
    failed: {
      id: "failed",
      label: "неудачное тестирование",
      fill: "rgba(239,68,68,0.22)",
      stroke: "#b91c1c",
      chipBg: "#fee2e2",
      chipFg: "#991b1b",
      chipBorder: "#fca5a5",
    },
    works: {
      id: "works",
      label: "работает",
      fill: "rgba(34,197,94,0.25)",
      stroke: "#15803d",
      chipBg: "#dcfce7",
      chipFg: "#166534",
      chipBorder: "#86efac",
    },
  };

  function schemeStatusMeta(status) {
    return SCHEME_STATUS[status] || SCHEME_STATUS.untested;
  }

  function productPairInfo(feat, productsTable) {
    const info = (feat && feat.info) || {};
    let fp = info.forward_primer || "";
    let rp = info.reverse_primer || "";
    if ((!fp || !rp) && productsTable && productsTable.length) {
      const match = productsTable.find((p) => {
        const nm = String(p.Forward_primer || "") + " + " + String(p.Reverse_primer || "");
        if (nm === feat.name) return true;
        try {
          return Number(p.Amplicon_start) === feat.start && Number(p.Amplicon_end) === feat.end;
        } catch (_) {
          return false;
        }
      });
      if (match) {
        if (!fp) fp = match.Forward_primer || "";
        if (!rp) rp = match.Reverse_primer || "";
      }
    }
    if (!fp || !rp) {
      const parts = String((feat && feat.name) || "").split(/\s*\+\s*/);
      if (!fp) fp = (parts[0] || "").trim();
      if (!rp) rp = (parts[1] || "").trim();
    }
    fp = String(fp || "").trim();
    rp = String(rp || "").trim();
    const pair = fp && rp ? fp + "+" + rp : String((feat && feat.name) || "").replace(/\s*\+\s*/g, "+");
    return { fp: fp, rp: rp, pair: pair };
  }

  function formatTmRange(tmF, tmR) {
    const vals = [];
    if (tmF != null && !isNaN(Number(tmF))) vals.push(Number(tmF));
    if (tmR != null && !isNaN(Number(tmR))) vals.push(Number(tmR));
    if (!vals.length) return "—";
    const lo = Math.min.apply(null, vals);
    const hi = Math.max.apply(null, vals);
    if (Math.abs(lo - hi) < 1e-9) return lo.toFixed(1);
    return lo.toFixed(1) + "-" + hi.toFixed(1);
  }

  function formatProductSchemeLine(p) {
    const num = p.num != null ? p.num : "?";
    return (
      "#" +
      num +
      " - " +
      (p.pair || p.name || "?") +
      ": Длина " +
      p.size +
      " bp (" +
      p.start +
      "–" +
      p.end +
      ") Tm F/R: " +
      formatTmRange(p.tm_f, p.tm_r)
    );
  }

  function formatRoundTitle(rd) {
    const pairs = (rd.products || []).map((p) => p.pair || p.name).filter(Boolean);
    return "Раунд " + rd.round + " · " + (pairs.length ? pairs.join("; ") : "—");
  }

  function productTm(feat, productsTable) {
    const info = feat.info || {};
    let tmF = info.Tm_forward != null && info.Tm_forward !== "" ? Number(info.Tm_forward) : NaN;
    let tmR = info.Tm_reverse != null && info.Tm_reverse !== "" ? Number(info.Tm_reverse) : NaN;
    if ((isNaN(tmF) || isNaN(tmR)) && productsTable && productsTable.length) {
      const match = productsTable.find((p) => {
        const nm = String(p.Forward_primer || "") + " + " + String(p.Reverse_primer || "");
        if (nm === feat.name) return true;
        try {
          return Number(p.Amplicon_start) === feat.start && Number(p.Amplicon_end) === feat.end;
        } catch (_) {
          return false;
        }
      });
      if (match) {
        if (isNaN(tmF)) tmF = Number(match.Tm_F);
        if (isNaN(tmR)) tmR = Number(match.Tm_R);
      }
    }
    return { tmF: isNaN(tmF) ? null : tmF, tmR: isNaN(tmR) ? null : tmR };
  }

  function buildSchemeRounds(scheme, productFeats, productsTable, productNumById) {
    const byId = {};
    (productFeats || []).forEach((f) => {
      byId[f.id] = f;
    });
    const assignments = scheme.assignments || {};
    const roundMap = {};
    Object.keys(assignments).forEach((pid) => {
      const r = Number(assignments[pid]) || 1;
      if (!roundMap[r]) roundMap[r] = [];
      roundMap[r].push(pid);
    });
    const roundNums = Object.keys(roundMap)
      .map(Number)
      .sort((a, b) => a - b);
    return roundNums.map((r) => {
      const feats = roundMap[r].map((id) => byId[id]).filter(Boolean);
      const tms = [];
      const products = feats.map((f) => {
        const tm = productTm(f, productsTable);
        if (tm.tmF != null) tms.push(tm.tmF);
        if (tm.tmR != null) tms.push(tm.tmR);
        const pairInfo = productPairInfo(f, productsTable);
        const num =
          productNumById && productNumById[f.id] != null ? productNumById[f.id] : null;
        return {
          id: f.id,
          num: num,
          name: f.name,
          pair: pairInfo.pair,
          fp: pairInfo.fp,
          rp: pairInfo.rp,
          start: f.start,
          end: f.end,
          size: f.end - f.start + 1,
          tm_f: tm.tmF,
          tm_r: tm.tmR,
          tm_range: formatTmRange(tm.tmF, tm.tmR),
        };
      });
      // сортируем по номеру продукта, как на общей схеме
      products.sort((a, b) => {
        const na = a.num != null ? a.num : 1e9;
        const nb = b.num != null ? b.num : 1e9;
        return na - nb || a.start - b.start;
      });
      const title = formatRoundTitle({ round: r, products: products });
      return {
        round: r,
        title: title,
        product_ids: products.map((p) => p.id),
        products: products,
        tm_min: tms.length ? Math.min.apply(null, tms) : null,
        tm_max: tms.length ? Math.max.apply(null, tms) : null,
        span_start: products.length ? Math.min.apply(null, products.map((p) => p.start)) : null,
        span_end: products.length ? Math.max.apply(null, products.map((p) => p.end)) : null,
      };
    });
  }

  function primersInsideProductRanges(primerFeats, products) {
    if (!primerFeats || !primerFeats.length || !products || !products.length) return [];
    return primerFeats.filter((pr) => {
      const ps = Number(pr.start);
      const pe = Number(pr.end);
      if (isNaN(ps) || isNaN(pe)) return false;
      return products.some((p) => ps >= p.start && pe <= p.end);
    });
  }

  const SEQ_PRIMER_STATUS = {
    none: { id: "none", label: "не проводилось", fill: "#e2e8f0", stroke: "#94a3b8", text: "#334155" },
    ok: { id: "ok", label: "получилось", fill: "#86efac", stroke: "#16a34a", text: "#14532d" },
    fail: { id: "fail", label: "не получилось", fill: "#fca5a5", stroke: "#dc2626", text: "#7f1d1d" },
  };

  function nextSeqPrimerStatus(cur) {
    if (cur === "ok") return "fail";
    if (cur === "fail") return "none";
    return "ok";
  }

  function seqPrimerStatusMeta(status) {
    return SEQ_PRIMER_STATUS[status] || SEQ_PRIMER_STATUS.none;
  }

  function schemeSpan(roundsDetail) {
    const starts = [];
    const ends = [];
    (roundsDetail || []).forEach((rd) => {
      (rd.products || []).forEach((p) => {
        starts.push(p.start);
        ends.push(p.end);
      });
    });
    if (!starts.length) return null;
    return { start: Math.min.apply(null, starts), end: Math.max.apply(null, ends) };
  }

  function primerNamesFromProduct(p) {
    const info = p.info || {};
    const fp = info.forward_primer || (p.name || "").split("+")[0];
    const rp = info.reverse_primer || (p.name || "").split("+")[1];
    return {
      fp: String(fp || "").trim(),
      rp: String(rp || "").trim(),
    };
  }

  function avgConservation(primerNames, consByName) {
    if (!consByName || !primerNames.length) return null;
    const vals = [];
    primerNames.forEach((n) => {
      const r = consByName[n];
      if (r && r.S_index_avg !== null && r.S_index_avg !== undefined && !isNaN(Number(r.S_index_avg))) {
        vals.push(Number(r.S_index_avg));
      }
    });
    if (!vals.length) return null;
    return vals.reduce((a, b) => a + b, 0) / vals.length;
  }

  /** overhang для покрытия «с выходом за границы»; Infinity если не начинается/не заканчивается снаружи */
  function outsideOverhang(spanStart, spanEnd, aStart, aEnd) {
    if (!(spanStart < aStart && spanEnd > aEnd)) return Infinity;
    return (aStart - spanStart) + (spanEnd - aEnd);
  }

  /** покрывает ли набор интервалов отрезок [aStart, aEnd] */
  function unionCovers(intervals, aStart, aEnd) {
    const sorted = intervals.slice().sort((a, b) => a.start - b.start);
    let coverUntil = aStart - 1;
    for (let i = 0; i < sorted.length; i++) {
      const it = sorted[i];
      if (it.start > coverUntil + 1) return false;
      if (it.end > coverUntil) coverUntil = it.end;
      if (coverUntil >= aEnd) return true;
    }
    return coverUntil >= aEnd;
  }

  /** жадное минимальное покрытие интервала продуктами */
  function greedyCover(products, aStart, aEnd) {
    const overlapping = products
      .filter((p) => p.end >= aStart && p.start <= aEnd)
      .slice()
      .sort((a, b) => a.start - b.start);
    if (!overlapping.length) return null;

    const chosen = [];
    let coverUntil = aStart - 1;
    let i = 0;
    while (coverUntil < aEnd) {
      let best = null;
      let bestEnd = coverUntil;
      while (i < overlapping.length && overlapping[i].start <= coverUntil + 1) {
        if (overlapping[i].end > bestEnd) {
          best = overlapping[i];
          bestEnd = overlapping[i].end;
        }
        i += 1;
      }
      // также просмотрим уже пройденные, которые стартуют <= coverUntil+1
      overlapping.forEach((p) => {
        if (p.start <= coverUntil + 1 && p.end > bestEnd) {
          best = p;
          bestEnd = p.end;
        }
      });
      if (!best || bestEnd <= coverUntil) return null;
      chosen.push(best);
      coverUntil = bestEnd;
      if (chosen.length > 30) return null;
    }
    return chosen;
  }

  function rankSolutionsForAnnotation(annot, productFeats, consByName) {
    const aStart = annot.start;
    const aEnd = annot.end;
    const solutions = [];

    // 1) одиночные продукты, выходящие за границы
    productFeats.forEach((p) => {
      const oh = outsideOverhang(p.start, p.end, aStart, aEnd);
      if (!isFinite(oh)) return;
      const names = primerNamesFromProduct(p);
      const primers = [names.fp, names.rp].filter(Boolean);
      solutions.push({
        id: "single-" + p.id,
        tier: 0,
        kind: "single",
        products: [p],
        productIds: [p.id],
        primers: primers,
        spanStart: p.start,
        spanEnd: p.end,
        overhang: oh,
        avgS: avgConservation(primers, consByName),
        label: p.name || p.id,
      });
    });

    // 2) пары продуктов, вместе покрывающие участок
    const candidates = productFeats.filter((p) => p.end >= aStart && p.start <= aEnd);
    for (let i = 0; i < candidates.length; i++) {
      for (let j = i + 1; j < candidates.length; j++) {
        const pair = [candidates[i], candidates[j]];
        if (!unionCovers(pair, aStart, aEnd)) continue;
        const spanStart = Math.min(pair[0].start, pair[1].start);
        const spanEnd = Math.max(pair[0].end, pair[1].end);
        const oh = outsideOverhang(spanStart, spanEnd, aStart, aEnd);
        if (!isFinite(oh)) continue;
        // пропускаем, если один из них уже сам покрывает — это tier 0
        if (
          isFinite(outsideOverhang(pair[0].start, pair[0].end, aStart, aEnd)) ||
          isFinite(outsideOverhang(pair[1].start, pair[1].end, aStart, aEnd))
        ) {
          continue;
        }
        const primers = [];
        pair.forEach((p) => {
          const n = primerNamesFromProduct(p);
          if (n.fp) primers.push(n.fp);
          if (n.rp) primers.push(n.rp);
        });
        solutions.push({
          id: "pair-" + pair[0].id + "-" + pair[1].id,
          tier: 1,
          kind: "multi",
          products: pair,
          productIds: pair.map((p) => p.id),
          primers: Array.from(new Set(primers)),
          spanStart: spanStart,
          spanEnd: spanEnd,
          overhang: oh,
          avgS: avgConservation(Array.from(new Set(primers)), consByName),
          label: pair.map((p) => p.name).join(" | "),
        });
      }
    }

    // 3) жадное покрытие ≥3 продуктов (если одиночные/пары не закрыли логику)
    const greedy = greedyCover(productFeats, aStart, aEnd);
    if (greedy && greedy.length >= 2) {
      const anySingle = greedy.some((p) => isFinite(outsideOverhang(p.start, p.end, aStart, aEnd)));
      if (!anySingle && unionCovers(greedy, aStart, aEnd)) {
        const spanStart = Math.min.apply(null, greedy.map((p) => p.start));
        const spanEnd = Math.max.apply(null, greedy.map((p) => p.end));
        const oh = outsideOverhang(spanStart, spanEnd, aStart, aEnd);
        if (isFinite(oh)) {
          const primers = [];
          greedy.forEach((p) => {
            const n = primerNamesFromProduct(p);
            if (n.fp) primers.push(n.fp);
            if (n.rp) primers.push(n.rp);
          });
          const uniqPrimers = Array.from(new Set(primers));
          const id = "greedy-" + greedy.map((p) => p.id).join("_");
          if (!solutions.some((s) => s.id === id || (s.tier === 1 && s.productIds.join() === greedy.map((p) => p.id).join()))) {
            solutions.push({
              id: id,
              tier: 1,
              kind: "multi",
              products: greedy,
              productIds: greedy.map((p) => p.id),
              primers: uniqPrimers,
              spanStart: spanStart,
              spanEnd: spanEnd,
              overhang: oh,
              avgS: avgConservation(uniqPrimers, consByName),
              label: greedy.map((p) => p.name).join(" | "),
            });
          }
        }
      }
    }

    solutions.sort((a, b) => {
      if (a.tier !== b.tier) return a.tier - b.tier;
      if (a.overhang !== b.overhang) return a.overhang - b.overhang;
      // больше продуктов в multi — чуть хуже при равном overhang
      if (a.tier === 1 && a.products.length !== b.products.length) {
        return a.products.length - b.products.length;
      }
      const sa = a.avgS;
      const sb = b.avgS;
      if (sa !== null && sb !== null && sa !== sb) return sb - sa;
      if (sa !== null && sb === null) return -1;
      if (sa === null && sb !== null) return 1;
      return String(a.label).localeCompare(String(b.label));
    });

    return solutions;
  }

  function AnnotationMapView(props) {
    const data = props.data || { length: 0, features: [] };
    const conservation = props.conservation || [];
    const productsTable = props.products || [];
    const schemes = props.schemes || [];
    const onSchemesChange = props.onSchemesChange || (() => {});
    const onLaneMapsChange = props.onLaneMapsChange || null;
    const emptyLanes = { ref: {}, primer: {}, product: {} };
    const [showRef, setShowRef] = React.useState(true);
    const [showPrimers, setShowPrimers] = React.useState(true);
    const [showProducts, setShowProducts] = React.useState(false);
    const [showSchemes, setShowSchemes] = React.useState(true);
    const [zoom, setZoom] = React.useState(1);
    const [tooltip, setTooltip] = React.useState(null);
    const [plotWidth, setPlotWidth] = React.useState(900);
    const [selectedAnnotId, setSelectedAnnotId] = React.useState(null);
    const [activeSolutionId, setActiveSolutionId] = React.useState(null);
    const [laneMapsLocal, setLaneMapsLocal] = React.useState(emptyLanes);
    const laneMaps =
      props.laneMaps && typeof props.laneMaps === "object" ? props.laneMaps : laneMapsLocal;
    const setLaneMaps = (updater) => {
      const prev = laneMaps || emptyLanes;
      const next = typeof updater === "function" ? updater(prev) : updater;
      if (onLaneMapsChange) onLaneMapsChange(next);
      else setLaneMapsLocal(next);
    };
    const [dragPreview, setDragPreview] = React.useState(null);
    const [selectedProductIds, setSelectedProductIds] = React.useState({});
    const [schemeDialog, setSchemeDialog] = React.useState(null); // draft scheme being edited
    const dragRef = React.useRef(null);
    const wrapRef = React.useRef(null);
    const movedRef = React.useRef(false);

    const consByName = React.useMemo(() => {
      const m = {};
      (conservation || []).forEach((r) => {
        if (r && r.Primer) m[r.Primer] = r;
      });
      return m;
    }, [conservation]);

    React.useEffect(() => {
      const el = wrapRef.current;
      if (!el || typeof ResizeObserver === "undefined") {
        const onResize = () => {
          const w = (wrapRef.current && wrapRef.current.clientWidth) || window.innerWidth - 80;
          setPlotWidth(Math.max(320, w - 180));
        };
        onResize();
        window.addEventListener("resize", onResize);
        return () => window.removeEventListener("resize", onResize);
      }
      const ro = new ResizeObserver((entries) => {
        for (const entry of entries) {
          setPlotWidth(Math.max(320, entry.contentRect.width - 180));
        }
      });
      ro.observe(el);
      setPlotWidth(Math.max(320, el.clientWidth - 180));
      return () => ro.disconnect();
    }, []);

    const allProducts = React.useMemo(
      () => (data.features || []).filter((f) => f.category === "product"),
      [data]
    );

    const allPrimers = React.useMemo(
      () => (data.features || []).filter((f) => f.category === "primer"),
      [data]
    );

    // Стабильные номера продуктов на всей карте (по геномному порядку)
    const productNumById = React.useMemo(() => {
      const sorted = allProducts.slice().sort((a, b) => {
        if (a.start !== b.start) return a.start - b.start;
        if (a.end !== b.end) return a.end - b.end;
        return String(a.name || "").localeCompare(String(b.name || ""));
      });
      const map = {};
      sorted.forEach((p, i) => {
        map[p.id] = i + 1;
      });
      return map;
    }, [allProducts]);

    const selectedAnnot = React.useMemo(() => {
      if (!selectedAnnotId) return null;
      return (data.features || []).find((f) => f.id === selectedAnnotId) || null;
    }, [data, selectedAnnotId]);

    const ranked = React.useMemo(() => {
      if (!selectedAnnot) return [];
      return rankSolutionsForAnnotation(selectedAnnot, allProducts, consByName);
    }, [selectedAnnot, allProducts, consByName]);

    React.useEffect(() => {
      if (selectedAnnotId) setShowProducts(true);
      setActiveSolutionId(null);
    }, [selectedAnnotId]);

    const activeSolution =
      (activeSolutionId && ranked.find((s) => s.id === activeSolutionId)) || ranked[0] || null;

    const highlightProductIds = React.useMemo(() => {
      const set = new Set();
      if (activeSolution) activeSolution.productIds.forEach((id) => set.add(id));
      return set;
    }, [activeSolution]);

    const highlightPrimerNames = React.useMemo(() => {
      const set = new Set();
      if (activeSolution) activeSolution.primers.forEach((n) => set.add(n));
      // также все праймеры из топ-решений для мягкой подсветки
      ranked.slice(0, 8).forEach((s) => s.primers.forEach((n) => set.add(n)));
      return set;
    }, [activeSolution, ranked]);

    const productRankInfo = React.useMemo(() => {
      const bestRank = {};
      const ordered = [];
      const seen = new Set();
      ranked.forEach((s, idx) => {
        (s.products || []).forEach((p) => {
          if (bestRank[p.id] === undefined) bestRank[p.id] = idx;
          if (!seen.has(p.id)) {
            seen.add(p.id);
            ordered.push(p);
          }
        });
      });
      return { bestRank: bestRank, ordered: ordered };
    }, [ranked]);

    const features = React.useMemo(() => {
      const all = data.features || [];
      return all.filter((f) => {
        if (f.category === "ref") return showRef;
        if (f.category === "primer") return showPrimers;
        if (f.category === "product") {
          if (!showProducts) return false;
          // при выбранной аннотации показываем только ранжированные + пересекающие участок
          if (selectedAnnot) {
            if (productRankInfo.bestRank[f.id] !== undefined) return true;
            return f.end >= selectedAnnot.start && f.start <= selectedAnnot.end;
          }
          return true;
        }
        return true;
      });
    }, [data, showRef, showPrimers, showProducts, selectedAnnot, productRankInfo]);

    const enrichedSchemes = React.useMemo(() => {
      return (schemes || []).map((sch) => {
        const rounds_detail = buildSchemeRounds(sch, allProducts, productsTable, productNumById);
        const span = schemeSpan(rounds_detail);
        const statusMeta = schemeStatusMeta(sch.status);
        const round2 = rounds_detail.find((rd) => rd.round === 2);
        const seqProducts = round2
          ? round2.products
          : sch.nested
            ? []
            : ((rounds_detail[0] && rounds_detail[0].products) || []);
        const seqPrimers = primersInsideProductRanges(allPrimers, seqProducts).map((pr) => {
          const st = (sch.seq_primer_status && sch.seq_primer_status[pr.id]) || "none";
          return Object.assign({}, pr, {
            seq_status: st,
            seq_meta: seqPrimerStatusMeta(st),
          });
        });
        return Object.assign({}, sch, {
          rounds_detail: rounds_detail,
          span: span,
          status: statusMeta.id,
          color: statusMeta.fill,
          stroke: statusMeta.stroke,
          statusMeta: statusMeta,
          seq_primers: seqPrimers,
          seq_products: seqProducts,
        });
      });
    }, [schemes, allProducts, allPrimers, productsTable, productNumById]);

    const layout = React.useMemo(() => {
      const groups = [];
      const schemeBlocks = [];
      let yCursor = 36;
      // при выборе аннотации продукты поднимаем сразу под референс — но схемы всегда между референсом и остальным
      const order = selectedAnnot
        ? ["ref", "product", "primer"]
        : CATEGORY_ORDER;

      const packCategory = (cat, items) => {
        let packedResult;
        if (cat === "product" && selectedAnnot && ranked.length) {
          const rankedIds = new Set(productRankInfo.ordered.map((p) => p.id));
          const top = productRankInfo.ordered.filter((p) =>
            items.some((it) => it.id === p.id)
          );
          const byId = {};
          items.forEach((it) => {
            byId[it.id] = it;
          });
          const orderedItems = top.map((p) => byId[p.id]).filter(Boolean);
          const rest = items
            .filter((it) => !rankedIds.has(it.id))
            .slice()
            .sort((a, b) => a.start - b.start || a.end - b.end);
          const sequence = orderedItems.concat(rest).map((it) =>
            Object.assign({}, it, {
              visualRank: productNumById[it.id] != null ? productNumById[it.id] : null,
            })
          );
          const TOP_N = 50;
          if (sequence.length <= TOP_N) {
            packedResult = packLanesStrictRank(sequence);
            // packLanesStrictRank перезаписывает visualRank индексом — вернём стабильные номера
            packedResult = {
              packed: packedResult.packed.map((p) =>
                Object.assign({}, p, {
                  visualRank: productNumById[p.id] != null ? productNumById[p.id] : p.visualRank,
                })
              ),
              laneCount: packedResult.laneCount,
            };
          } else {
            const head = sequence.slice(0, TOP_N);
            const tail = sequence.slice(TOP_N);
            const headPack = packLanesStrictRank(head);
            const tailPack = packLanesByPriority(tail);
            const offset = headPack.laneCount;
            const headPacked = headPack.packed.map((p) =>
              Object.assign({}, p, {
                visualRank: productNumById[p.id] != null ? productNumById[p.id] : p.visualRank,
              })
            );
            const tailPacked = tailPack.packed.map((p) =>
              Object.assign({}, p, {
                lane: p.lane + offset,
                visualRank: productNumById[p.id] != null ? productNumById[p.id] : null,
              })
            );
            packedResult = {
              packed: headPacked.concat(tailPacked),
              laneCount: offset + tailPack.laneCount,
            };
          }
          const effectiveMap =
            dragPreview && dragPreview.category === cat
              ? dragPreview.lanes
              : laneMaps[cat];
          if (effectiveMap && Object.keys(effectiveMap).length) {
            packedResult = packWithLaneMap(sequence, effectiveMap);
            packedResult = {
              packed: packedResult.packed.map((p) =>
                Object.assign({}, p, {
                  visualRank: productNumById[p.id] != null ? productNumById[p.id] : null,
                })
              ),
              laneCount: packedResult.laneCount,
            };
          }
        } else {
          const withNums = items.map((it) =>
            Object.assign({}, it, {
              visualRank:
                cat === "product" && productNumById[it.id] != null
                  ? productNumById[it.id]
                  : null,
            })
          );
          const effectiveMap =
            dragPreview && dragPreview.category === cat
              ? dragPreview.lanes
              : laneMaps[cat];
          packedResult = packWithLaneMap(withNums, effectiveMap);
          if (cat === "product") {
            packedResult = {
              packed: packedResult.packed.map((p) =>
                Object.assign({}, p, {
                  visualRank: productNumById[p.id] != null ? productNumById[p.id] : null,
                })
              ),
              laneCount: packedResult.laneCount,
            };
          }
        }
        return packedResult;
      };

      const appendSchemes = () => {
        if (!showSchemes || !enrichedSchemes.length) return;
        const rowH = 16;
        const roundTitleH = 16;
        const gap = 3;
        const pad = 8;
        enrichedSchemes.forEach((sch) => {
          const headerH = sch.notes ? 42 : 28;
          let contentH = 0;
          (sch.rounds_detail || []).forEach((rd) => {
            contentH += roundTitleH + gap;
            contentH += ((rd.products && rd.products.length) || 0) * (rowH + gap);
          });
          if (!contentH) contentH = rowH;
          const innerH = headerH + contentH + pad;
          schemeBlocks.push({
            scheme: sch,
            y0: yCursor,
            height: innerH,
            headerH: headerH,
            rowH: rowH,
            roundTitleH: roundTitleH,
            gap: gap,
          });
          yCursor += innerH + 14;
        });
      };

      order.forEach((cat) => {
        let items = features.filter((f) => f.category === cat);
        if (!items.length) {
          if (cat === "ref") appendSchemes();
          return;
        }

        const packedResult = packCategory(cat, items);
        const blockH = cat === "product" && selectedAnnot ? 20 : 18;
        const gap = 4;
        const headerH = 22;
        const groupH = headerH + packedResult.laneCount * (blockH + gap) + 8;
        groups.push({
          category: cat,
          label:
            cat === "product" && selectedAnnot
              ? "ПЦР‑продукты (сверху — наиболее пригодные)"
              : CATEGORY_LABELS[cat] || cat,
          packed: packedResult.packed,
          laneCount: packedResult.laneCount,
          y0: yCursor,
          headerH: headerH,
          blockH: blockH,
          gap: gap,
          height: groupH,
        });
        yCursor += groupH + 10;
        if (cat === "ref") appendSchemes();
      });

      // если референса нет (скрыт), схемы всё равно показываем сверху
      if (showSchemes && enrichedSchemes.length && !groups.some((g) => g.category === "ref")) {
        // если appendSchemes ещё не вызывался (ref скрыт и не в order processing with empty)
        if (!schemeBlocks.length) {
          // reset: rebuild with schemes first is messy; call append if empty
          appendSchemes();
        }
      }

      return { groups: groups, schemeBlocks: schemeBlocks, totalH: yCursor + 20 };
    }, [
      features,
      selectedAnnot,
      ranked,
      productRankInfo,
      laneMaps,
      dragPreview,
      productNumById,
      showSchemes,
      enrichedSchemes,
    ]);

    const schemeLayout = layout;

    const seqLen = Math.max(1, Number(data.length) || 1);
    const plotW = plotWidth * zoom;
    const leftPad = 140;
    const rightPad = 24;
    const svgW = leftPad + plotW + rightPad;
    const svgH = Math.max(220, layout.totalH);

    const xOf = (pos) => leftPad + ((pos - 1) / seqLen) * plotW;
    const wOf = (start, end) => Math.max(3, ((end - start + 1) / seqLen) * plotW);

    const laneFromClientY = (group, clientY) => {
      const svg = wrapRef.current && wrapRef.current.querySelector("svg");
      if (!svg) return 0;
      const rect = svg.getBoundingClientRect();
      // учитываем viewBox масштабирование
      const vb = svg.viewBox && svg.viewBox.baseVal;
      const scaleY = vb && vb.height ? rect.height / vb.height : 1;
      const yInSvg = (clientY - rect.top) / scaleY;
      const rowH = group.blockH + group.gap;
      const local = yInSvg - group.y0 - group.headerH;
      const lane = Math.floor(local / rowH);
      return Math.max(0, lane);
    };

    const onPointerDown = (ev, feat, group) => {
      ev.preventDefault();
      ev.stopPropagation();
      movedRef.current = false;
      dragRef.current = {
        id: feat.id,
        feat: feat,
        group: group,
        category: group.category,
        startClientY: ev.clientY,
        fromLane: feat.lane,
      };
      try {
        ev.currentTarget.setPointerCapture(ev.pointerId);
      } catch (_) {}
    };

    const onPointerMove = (ev) => {
      const d = dragRef.current;
      if (!d) return;
      const dy = ev.clientY - d.startClientY;
      if (Math.abs(dy) > 4) movedRef.current = true;
      if (!movedRef.current) return;

      const targetLane = laneFromClientY(d.group, ev.clientY);
      const items = (layout.groups.find((g) => g.category === d.category) || {}).packed || [];
      const rawItems = features.filter((f) => f.category === d.category);
      const baseMap =
        laneMaps[d.category] && Object.keys(laneMaps[d.category]).length
          ? laneMaps[d.category]
          : Object.fromEntries(items.map((p) => [p.id, p.lane]));
      const nextLanes = moveItemToLane(rawItems.length ? rawItems : items, baseMap, d.id, targetLane);
      d.lastLanes = nextLanes;
      d.targetLane = targetLane;
      setDragPreview({ category: d.category, lanes: nextLanes, targetLane: targetLane });
    };

    const onPointerUp = () => {
      const d = dragRef.current;
      dragRef.current = null;
      if (!d) return;

      if (!movedRef.current) {
        setDragPreview(null);
        if (d.feat && d.feat.category === "ref") {
          setSelectedAnnotId((prev) => (prev === d.feat.id ? null : d.feat.id));
        } else if (d.feat && d.feat.category === "product") {
          setSelectedProductIds((prev) => {
            const next = Object.assign({}, prev);
            next[d.feat.id] = !next[d.feat.id];
            return next;
          });
        }
        return;
      }

      if (d.lastLanes) {
        setLaneMaps((prev) =>
          Object.assign({}, prev, {
            [d.category]: Object.assign({}, d.lastLanes),
          })
        );
      }
      setDragPreview(null);
    };

    const resetPositions = () => {
      setLaneMaps({ ref: {}, primer: {}, product: {} });
      setDragPreview(null);
    };

    const selectedCount = Object.keys(selectedProductIds).filter((k) => selectedProductIds[k]).length;

    const openCreateScheme = () => {
      const ids = Object.keys(selectedProductIds).filter((k) => selectedProductIds[k]);
      if (!ids.length) {
        alert("Отметьте галочками один или несколько продуктов (включите слой «Продукты» и кликните по ним, либо отметьте в списке ниже).");
        return;
      }
      const assignments = {};
      ids.forEach((id) => {
        assignments[id] = 1;
      });
      setSchemeDialog({
        id: null,
        name: "Схема ПЦР " + (schemes.length + 1),
        nested: false,
        status: "untested",
        notes: "",
        productIds: ids,
        assignments: assignments,
        seq_primer_status: {},
      });
    };

    const openEditScheme = (sch) => {
      const assignments = Object.assign({}, sch.assignments || {});
      const productIds = Object.keys(assignments);
      setSchemeDialog({
        id: sch.id,
        name: sch.name || "Схема ПЦР",
        nested: !!sch.nested,
        status: sch.status || "untested",
        notes: sch.notes || "",
        productIds: productIds,
        assignments: assignments,
        seq_primer_status: Object.assign({}, sch.seq_primer_status || {}),
      });
    };

    const saveSchemeDialog = () => {
      if (!schemeDialog) return;
      const name = String(schemeDialog.name || "").trim() || "Схема ПЦР";
      const assignments = {};
      (schemeDialog.productIds || []).forEach((pid) => {
        const r = schemeDialog.nested ? Number(schemeDialog.assignments[pid]) || 1 : 1;
        assignments[pid] = r;
      });
      const nextScheme = {
        id: schemeDialog.id || "scheme-" + Date.now().toString(36),
        name: name,
        nested: !!schemeDialog.nested,
        status: schemeDialog.status || "untested",
        notes: schemeDialog.notes || "",
        assignments: assignments,
        seq_primer_status: Object.assign({}, schemeDialog.seq_primer_status || {}),
      };
      let next;
      if (schemeDialog.id) {
        next = schemes.map((s) => (s.id === schemeDialog.id ? nextScheme : s));
      } else {
        next = schemes.concat([nextScheme]);
      }
      onSchemesChange(next);
      setSchemeDialog(null);
      setSelectedProductIds({});
    };

    const deleteScheme = (schemeId) => {
      if (!confirm("Удалить схему ПЦР?")) return;
      onSchemesChange(schemes.filter((s) => s.id !== schemeId));
    };

    const dialogPreviewRounds = React.useMemo(() => {
      if (!schemeDialog) return [];
      return buildSchemeRounds(
        {
          assignments: Object.fromEntries(
            (schemeDialog.productIds || []).map((pid) => [
              pid,
              schemeDialog.nested ? Number(schemeDialog.assignments[pid]) || 1 : 1,
            ])
          ),
        },
        allProducts,
        productsTable,
        productNumById
      );
    }, [schemeDialog, allProducts, productsTable, productNumById]);

    const dialogSeqPrimers = React.useMemo(() => {
      if (!schemeDialog) return [];
      const round2 = dialogPreviewRounds.find((rd) => rd.round === 2);
      const seqProducts = round2
        ? round2.products
        : schemeDialog.nested
          ? []
          : ((dialogPreviewRounds[0] && dialogPreviewRounds[0].products) || []);
      return primersInsideProductRanges(allPrimers, seqProducts);
    }, [schemeDialog, dialogPreviewRounds, allPrimers]);

    const renderSchemeBlock = (block) => {
      const sch = block.scheme;
      const meta = sch.statusMeta || schemeStatusMeta(sch.status);
      const span = sch.span;
      const x = span ? xOf(span.start) : leftPad;
      const w = span ? Math.max(wOf(span.start, span.end), 120) : plotW;
      const y = block.y0;
      const h = block.height;
      let rowY = y + block.headerH;
      const rows = [];
      (sch.rounds_detail || []).forEach((rd) => {
        rows.push(
          e(
            "text",
            {
              key: sch.id + "-rt-" + rd.round,
              x: x + 8,
              y: rowY + block.roundTitleH - 3,
              fontSize: 11,
              fill: "#0f172a",
              fontWeight: 700,
              style: { pointerEvents: "none" },
            },
            String(rd.title || formatRoundTitle(rd)).slice(0, Math.max(20, Math.floor(w / 7)))
          )
        );
        rowY += block.roundTitleH + block.gap;
        (rd.products || []).forEach((p) => {
          const line = formatProductSchemeLine(p);
          rows.push(
            e(
              "g",
              { key: sch.id + "-" + p.id + "-" + rd.round },
              e("rect", {
                x: xOf(p.start),
                y: rowY,
                width: Math.max(3, wOf(p.start, p.end)),
                height: block.rowH - 2,
                rx: 2,
                ry: 2,
                fill: meta.stroke,
                opacity: 0.35,
              }),
              e(
                "text",
                {
                  x: x + 10,
                  y: rowY + block.rowH - 4,
                  fontSize: 9,
                  fill: "#1e293b",
                  style: { pointerEvents: "none" },
                },
                line.slice(0, Math.max(24, Math.floor(w / 6)))
              )
            )
          );
          rowY += block.rowH + block.gap;
        });
      });
      return e(
        "g",
        { key: sch.id },
        e("text", { x: 8, y: y + 16, fontSize: 12, fill: "#0f172a", fontWeight: 600 }, "Схема"),
        e("rect", {
          x: x,
          y: y,
          width: Math.max(w, 40),
          height: h,
          rx: 8,
          ry: 8,
          fill: meta.fill,
          stroke: meta.stroke,
          strokeWidth: 1.5,
          strokeDasharray: "4 3",
        }),
        e(
          "text",
          {
            x: x + 8,
            y: y + 18,
            fontSize: 12,
            fill: "#0f172a",
            fontWeight: 700,
            style: { pointerEvents: "none" },
          },
          (sch.name || "Схема") +
            (sch.nested ? " · nested" : "") +
            " · " +
            meta.label
        ),
        sch.notes
          ? e(
              "text",
              {
                x: x + 8,
                y: y + 32,
                fontSize: 9,
                fill: "#475569",
                style: { pointerEvents: "none" },
              },
              String(sch.notes).slice(0, 80) + (String(sch.notes).length > 80 ? "…" : "")
            )
          : null,
        rows
      );
    };

    const renderBlock = (feat, group) => {
      const x = xOf(feat.start);
      const y = group.y0 + group.headerH + feat.lane * (group.blockH + group.gap);
      const w = wOf(feat.start, feat.end);
      const h = group.blockH;
      const fill = colorFor(feat.type);
      const isPrimer = feat.category === "primer" || feat.type === "primer";
      const isSelected = selectedAnnotId === feat.id;
      const isDragging = dragPreview && dragRef.current && dragRef.current.id === feat.id;

      let opacity = 0.9;
      let stroke = "#334155";
      let strokeWidth = 1;
      if (selectedAnnot) {
        if (feat.category === "ref") {
          opacity = isSelected ? 1 : 0.25;
          if (isSelected) {
            stroke = "#2563eb";
            strokeWidth = 2;
          }
        } else if (feat.category === "product") {
          opacity = highlightProductIds.has(feat.id) ? 0.95 : 0.12;
          if (highlightProductIds.has(feat.id)) {
            stroke = "#0369a1";
            strokeWidth = 2;
          }
        } else if (feat.category === "primer") {
          const hit = highlightPrimerNames.has(feat.name);
          opacity = hit ? 0.95 : 0.12;
          if (hit) {
            stroke = "#b45309";
            strokeWidth = 2;
          }
        }
      }

      const tipLines = [
        feat.name || feat.type,
        "type: " + feat.type,
        "coords: " + feat.start + ".." + feat.end + " (" + (feat.strand || ".") + ")",
      ];
      if (feat.category === "ref") tipLines.push("(клик — выбрать участок)");
      if (feat.category === "product" && selectedProductIds[feat.id]) {
        stroke = "#0ea5e9";
        strokeWidth = 2.5;
        opacity = Math.max(opacity, 0.95);
        tipLines.push("(выбран для схемы ПЦР)");
      }
      if (isDragging) {
        stroke = "#dc2626";
        strokeWidth = 2;
        opacity = 1;
      }
      const info = feat.info || {};
      [
        "product_size",
        "Tm_forward",
        "Tm_reverse",
        "Tm_diff",
        "specificity",
        "note",
        "ugene_group",
      ].forEach((k) => {
        if (info[k] !== undefined && info[k] !== null && String(info[k]).length) {
          tipLines.push(k + ": " + info[k]);
        }
      });

      const showTip = (clientX, clientY) => {
        setTooltip({ x: clientX + 12, y: clientY + 12, lines: tipLines });
      };

      const rankPrefix =
        feat.category === "product" && feat.visualRank
          ? "#" + feat.visualRank + " "
          : "";
      const baseLabel =
        w > 40
          ? String(feat.name || feat.type).slice(0, Math.max(4, Math.floor(w / 7)))
          : "";
      const shownLabel = (rankPrefix + baseLabel).trim() || (feat.visualRank ? "#" + feat.visualRank : "");

      const commonHandlers = {
        style: { cursor: "ns-resize" },
        onPointerDown: (ev) => onPointerDown(ev, feat, group),
        onPointerMove: onPointerMove,
        onPointerUp: onPointerUp,
        onPointerCancel: onPointerUp,
        onMouseEnter: (ev) => showTip(ev.clientX, ev.clientY),
        onMouseMove: (ev) => showTip(ev.clientX, ev.clientY),
        onMouseLeave: () => setTooltip(null),
      };

      if (isPrimer) {
        const ah = Math.min(10, Math.max(4, w * 0.35));
        let points;
        if (feat.strand === "-") {
          points = [
            [x + w, y],
            [x + ah, y],
            [x, y + h / 2],
            [x + ah, y + h],
            [x + w, y + h],
          ]
            .map((p) => p.join(","))
            .join(" ");
        } else {
          points = [
            [x, y],
            [x + Math.max(0, w - ah), y],
            [x + w, y + h / 2],
            [x + Math.max(0, w - ah), y + h],
            [x, y + h],
          ]
            .map((p) => p.join(","))
            .join(" ");
        }
        return e(
          "g",
          Object.assign({ key: feat.id }, commonHandlers),
          e("polygon", {
            points: points,
            fill: fill,
            stroke: stroke,
            strokeWidth: strokeWidth,
            opacity: opacity,
          }),
          baseLabel
            ? e(
                "text",
                {
                  x: x + w / 2,
                  y: y + h / 2 + 3,
                  textAnchor: "middle",
                  fontSize: 10,
                  fill: "#0f172a",
                  style: { pointerEvents: "none", userSelect: "none" },
                },
                baseLabel
              )
            : null
        );
      }

      return e(
        "g",
        Object.assign({ key: feat.id }, commonHandlers),
        e("rect", {
          x: x,
          y: y,
          width: w,
          height: h,
          rx: 4,
          ry: 4,
          fill: fill,
          stroke: stroke,
          strokeWidth: strokeWidth,
          opacity: opacity,
        }),
        shownLabel
          ? e(
              "text",
              {
                x: x + 4,
                y: y + h / 2 + 3,
                fontSize: 10,
                fill: "#0f172a",
                fontWeight: feat.visualRank ? 700 : 400,
                style: { pointerEvents: "none", userSelect: "none" },
              },
              shownLabel
            )
          : null
      );
    };

    const ticks = [];
    const stepGuess = Math.pow(10, Math.max(0, Math.floor(Math.log10(seqLen)) - 1));
    const step = stepGuess || 100;
    for (let p = 1; p <= seqLen; p += step) ticks.push(p);
    if (ticks[ticks.length - 1] !== seqLen) ticks.push(seqLen);

    const legendTypes = {};
    (data.features || []).forEach((f) => {
      legendTypes[f.type] = colorFor(f.type);
    });

    if (!data || !data.features || !data.features.length) {
      return e("p", { style: { color: "#64748b" } }, "Нет данных для схемы. Запустите анализ с GenBank‑референсом.");
    }

    const highlightBand =
      selectedAnnot && layout.groups.length
        ? e("rect", {
            x: xOf(selectedAnnot.start),
            y: 28,
            width: wOf(selectedAnnot.start, selectedAnnot.end),
            height: Math.max(40, svgH - 40),
            fill: "rgba(37, 99, 235, 0.14)",
            stroke: "rgba(37, 99, 235, 0.35)",
            strokeWidth: 1,
            style: { pointerEvents: "none" },
          })
        : null;

    return e(
      "div",
      { className: "annomap-wrap", ref: wrapRef },
      e(
        "div",
        { className: "annomap-controls" },
        e("label", null,
          e("input", { type: "checkbox", checked: showRef, onChange: (ev) => setShowRef(ev.target.checked) }),
          " Референс"
        ),
        e("label", null,
          e("input", { type: "checkbox", checked: showPrimers, onChange: (ev) => setShowPrimers(ev.target.checked) }),
          " Праймеры"
        ),
        e("label", null,
          e("input", { type: "checkbox", checked: showProducts, onChange: (ev) => setShowProducts(ev.target.checked) }),
          " Продукты"
        ),
        e("label", null,
          e("input", { type: "checkbox", checked: showSchemes, onChange: (ev) => setShowSchemes(ev.target.checked) }),
          " Схемы ПЦР"
        ),
        e("button", { type: "button", className: "btn-secondary", onClick: () => setZoom((z) => Math.max(0.5, +(z - 0.25).toFixed(2))) }, "−"),
        e("span", { style: { color: "#64748b", fontSize: 13 } }, "zoom " + zoom.toFixed(2)),
        e("button", { type: "button", className: "btn-secondary", onClick: () => setZoom((z) => Math.min(8, +(z + 0.25).toFixed(2))) }, "+"),
        e("button", { type: "button", className: "btn-secondary", onClick: resetPositions }, "Сбросить позиции"),
        selectedAnnot
          ? e("button", {
              type: "button",
              className: "btn-secondary",
              onClick: () => setSelectedAnnotId(null),
            }, "Снять выбор: " + (selectedAnnot.name || selectedAnnot.type))
          : null,
        e("button", {
          type: "button",
          className: "btn-primary",
          onClick: openCreateScheme,
          disabled: selectedCount === 0,
          style: selectedCount === 0 ? { opacity: 0.5, cursor: "not-allowed" } : null,
        }, "Создать схему ПЦР" + (selectedCount ? " (" + selectedCount + ")" : ""))
      ),
      allProducts.length > 0 &&
        e(
          "div",
          {
            className: "annomap-product-pick",
            style: {
              maxHeight: 140,
              overflow: "auto",
              border: "1px solid #cbd5e1",
              borderRadius: 8,
              padding: "6px 8px",
              marginBottom: 8,
              background: "#fff",
              fontSize: 12,
            },
          },
          e(
            "div",
            { style: { marginBottom: 4, color: "#475569", fontWeight: 600 } },
            "Продукты для схемы (галочки или клик по блоку на карте):"
          ),
          allProducts.map((p) =>
            e(
              "label",
              {
                key: p.id,
                style: {
                  display: "flex",
                  gap: 6,
                  alignItems: "center",
                  marginBottom: 2,
                  color: "#1e293b",
                },
              },
              e("input", {
                type: "checkbox",
                checked: !!selectedProductIds[p.id],
                onChange: (ev) => {
                  const checked = ev.target.checked;
                  setSelectedProductIds((prev) => {
                    const next = Object.assign({}, prev);
                    next[p.id] = checked;
                    return next;
                  });
                },
              }),
              "#" +
                (productNumById[p.id] != null ? productNumById[p.id] : "?") +
                "  " +
                (p.name || p.id) +
                "  [" +
                p.start +
                "–" +
                p.end +
                ", " +
                (p.end - p.start + 1) +
                " bp]"
            )
          )
        ),
      schemes.length > 0 &&
        e(
          "div",
          { className: "annomap-scheme-list" },
          schemes.map((sch) => {
            const meta = schemeStatusMeta(sch.status);
            return e(
              "span",
              {
                key: sch.id,
                className: "annomap-scheme-chip",
                style: {
                  display: "inline-flex",
                  gap: 6,
                  alignItems: "center",
                  background: meta.chipBg,
                  color: meta.chipFg,
                  borderColor: meta.chipBorder,
                },
              },
              (sch.name || "Схема") +
                (sch.nested ? " · nested" : "") +
                " · " +
                meta.label,
              e(
                "button",
                {
                  type: "button",
                  className: "btn-secondary",
                  style: { padding: "2px 8px", fontSize: 11 },
                  onClick: () => openEditScheme(sch),
                },
                "Изменить"
              ),
              e(
                "button",
                {
                  type: "button",
                  className: "btn-secondary",
                  style: { padding: "2px 8px", fontSize: 11 },
                  onClick: () => deleteScheme(sch.id),
                },
                "✕"
              )
            );
          })
        ),
      e(
        "div",
        { className: "annomap-legend" },
        Object.keys(legendTypes).map((t) =>
          e(
            "span",
            { key: t, className: "annomap-legend-item" },
            e("span", { className: "annomap-swatch", style: { background: legendTypes[t] } }),
            t
          )
        )
      ),
      e(
        "div",
        { className: "annomap-scroll" },
        e(
          "svg",
          {
            width: "100%",
            viewBox: "0 0 " + svgW + " " + svgH,
            preserveAspectRatio: "xMinYMin meet",
            className: "annomap-svg",
            style: { width: "100%", height: "auto", minHeight: Math.min(svgH, 480) + "px" },
            onClick: () => {},
          },
          highlightBand,
          e("line", {
            x1: leftPad,
            y1: 18,
            x2: leftPad + plotW,
            y2: 18,
            stroke: "#94a3b8",
            strokeWidth: 2,
          }),
          ticks.map((p) =>
            e(
              "g",
              { key: "t" + p },
              e("line", { x1: xOf(p), y1: 14, x2: xOf(p), y2: 22, stroke: "#94a3b8" }),
              e("text", { x: xOf(p), y: 12, textAnchor: "middle", fontSize: 10, fill: "#64748b" }, String(p))
            )
          ),
          layout.groups.map((group) =>
            e(
              "g",
              { key: group.category },
              e("text", { x: 8, y: group.y0 + 14, fontSize: 12, fill: "#0f172a", fontWeight: 600 }, group.label),
              e("line", {
                x1: leftPad,
                y1: group.y0 + group.headerH - 6,
                x2: leftPad + plotW,
                y2: group.y0 + group.headerH - 6,
                stroke: "#e2e8f0",
              }),
              group.packed.map((feat) => renderBlock(feat, group))
            )
          ),
          schemeLayout.schemeBlocks.map((block) => renderSchemeBlock(block))
        )
      ),
      selectedAnnot &&
        e(
          "div",
          { className: "annomap-rank" },
          e(
            "div",
            { style: { padding: "8px 10px", borderBottom: "1px solid #e2e8f0", background: "#f8fafc" } },
            e("strong", null, "Ранжирование для: "),
            selectedAnnot.name || selectedAnnot.type,
            " (",
            selectedAnnot.start,
            "–",
            selectedAnnot.end,
            ")",
            ranked.length
              ? e("span", { style: { color: "#64748b" } }, " — вариантов: " + ranked.length)
              : e("span", { style: { color: "#b91c1c" } }, " — подходящих продуктов не найдено")
          ),
          ranked.length > 0 &&
            e(
              "table",
              null,
              e(
                "thead",
                null,
                e("tr", null,
                  e("th", null, "#"),
                  e("th", null, "Тип"),
                  e("th", null, "Продукты"),
                  e("th", null, "Span"),
                  e("th", null, "Overhang"),
                  e("th", null, "Avg S"),
                  e("th", null, "Праймеры")
                )
              ),
              e(
                "tbody",
                null,
                ranked.map((s, idx) =>
                  e(
                    "tr",
                    {
                      key: s.id,
                      className: activeSolution && activeSolution.id === s.id ? "active" : "",
                      onClick: () => setActiveSolutionId(s.id),
                    },
                    e("td", null, idx + 1),
                    e("td", null, s.kind === "single" ? "один" : "набор (" + s.products.length + ")"),
                    e("td", null, s.label),
                    e("td", null, s.spanStart + "–" + s.spanEnd),
                    e("td", null, s.overhang),
                    e("td", null, s.avgS === null || s.avgS === undefined ? "—" : Number(s.avgS).toFixed(3)),
                    e("td", null, s.primers.join(", "))
                  )
                )
              )
            )
        ),
      tooltip
        ? e(
            "div",
            { className: "annomap-tooltip", style: { left: tooltip.x + "px", top: tooltip.y + "px" } },
            tooltip.lines.map((line, i) => e("div", { key: i }, line))
          )
        : null,
      schemeDialog &&
        e(
          "div",
          {
            className: "annomap-modal-backdrop",
            onClick: (ev) => {
              if (ev.target === ev.currentTarget) setSchemeDialog(null);
            },
          },
          e(
            "div",
            { className: "annomap-modal", onClick: (ev) => ev.stopPropagation() },
            e("h3", null, schemeDialog.id ? "Редактировать схему ПЦР" : "Создать схему ПЦР"),
            e(
              "div",
              { className: "field" },
              e("label", null, "Название"),
              e("input", {
                type: "text",
                value: schemeDialog.name,
                onChange: (ev) =>
                  setSchemeDialog(Object.assign({}, schemeDialog, { name: ev.target.value })),
              })
            ),
            e(
              "div",
              { className: "field" },
              e("label", null, "Статус схемы"),
              e(
                "select",
                {
                  value: schemeDialog.status || "untested",
                  onChange: (ev) =>
                    setSchemeDialog(Object.assign({}, schemeDialog, { status: ev.target.value })),
                },
                Object.keys(SCHEME_STATUS).map((key) =>
                  e("option", { key: key, value: key }, SCHEME_STATUS[key].label)
                )
              )
            ),
            e(
              "div",
              { className: "field" },
              e(
                "label",
                { style: { display: "flex", gap: 8, alignItems: "center" } },
                e("input", {
                  type: "checkbox",
                  checked: !!schemeDialog.nested,
                  onChange: (ev) => {
                    const nested = ev.target.checked;
                    const assignments = Object.assign({}, schemeDialog.assignments);
                    if (!nested) {
                      (schemeDialog.productIds || []).forEach((pid) => {
                        assignments[pid] = 1;
                      });
                    }
                    setSchemeDialog(Object.assign({}, schemeDialog, { nested: nested, assignments: assignments }));
                  },
                }),
                " Nested PCR (несколько раундов)"
              )
            ),
            e(
              "div",
              { className: "field" },
              e("label", null, "Продукты и раунды"),
              (schemeDialog.productIds || []).map((pid) => {
                const feat = allProducts.find((f) => f.id === pid);
                const pair = feat ? productPairInfo(feat, productsTable).pair : pid;
                const label = feat
                  ? pair + " [" + feat.start + "–" + feat.end + "]"
                  : pid;
                return e(
                  "div",
                  {
                    key: pid,
                    style: {
                      display: "flex",
                      gap: 8,
                      alignItems: "center",
                      marginBottom: 6,
                      fontSize: 13,
                    },
                  },
                  e("span", { style: { flex: 1 } }, label),
                  schemeDialog.nested
                    ? e(
                        "select",
                        {
                          value: String(schemeDialog.assignments[pid] || 1),
                          onChange: (ev) => {
                            const assignments = Object.assign({}, schemeDialog.assignments);
                            assignments[pid] = Number(ev.target.value) || 1;
                            setSchemeDialog(Object.assign({}, schemeDialog, { assignments: assignments }));
                          },
                          style: { width: 120 },
                        },
                        [1, 2, 3, 4, 5].map((n) =>
                          e("option", { key: n, value: String(n) }, "Раунд " + n)
                        )
                      )
                    : e("span", { style: { color: "#64748b" } }, "раунд 1")
                );
              })
            ),
            dialogPreviewRounds.map((rd) =>
              e(
                "div",
                { key: "rd-" + rd.round, className: "annomap-round-box" },
                e(
                  "div",
                  { style: { fontWeight: 700, marginBottom: 6 } },
                  rd.title || formatRoundTitle(rd)
                ),
                (rd.products || []).map((p) =>
                  e(
                    "div",
                    { key: p.id, style: { fontSize: 12, color: "#334155", marginBottom: 4 } },
                    formatProductSchemeLine(p)
                  )
                )
              )
            ),
            e(
              "div",
              { className: "field" },
              e("label", null, "Праймеры для секвенирования"),
              e(
                "p",
                { style: { fontSize: 12, color: "#64748b", margin: "0 0 8px" } },
                schemeDialog.nested
                  ? "Праймеры внутри диапазона продуктов 2-го раунда (включая концы). Клик: не проводилось → получилось → не получилось."
                  : "Для обычной ПЦР — праймеры внутри продуктов схемы. Клик: не проводилось → получилось → не получилось."
              ),
              dialogSeqPrimers.length
                ? e(
                    "div",
                    { style: { display: "flex", flexWrap: "wrap", gap: 8 } },
                    dialogSeqPrimers.map((pr) => {
                      const st =
                        (schemeDialog.seq_primer_status && schemeDialog.seq_primer_status[pr.id]) ||
                        "none";
                      const sm = seqPrimerStatusMeta(st);
                      return e(
                        "button",
                        {
                          key: pr.id,
                          type: "button",
                          title: sm.label + " (клик — сменить статус)",
                          onClick: () => {
                            const seq_primer_status = Object.assign(
                              {},
                              schemeDialog.seq_primer_status || {}
                            );
                            seq_primer_status[pr.id] = nextSeqPrimerStatus(st);
                            setSchemeDialog(
                              Object.assign({}, schemeDialog, {
                                seq_primer_status: seq_primer_status,
                              })
                            );
                          },
                          style: {
                            border: "1px solid " + sm.stroke,
                            background: sm.fill,
                            color: sm.text,
                            borderRadius: 8,
                            padding: "6px 10px",
                            cursor: "pointer",
                            fontSize: 12,
                            textAlign: "left",
                          },
                        },
                        e("div", { style: { fontWeight: 600 } }, pr.name || pr.id),
                        e(
                          "div",
                          { style: { fontSize: 11, opacity: 0.85 } },
                          pr.start + "–" + pr.end + " · " + sm.label
                        )
                      );
                    })
                  )
                : e(
                    "div",
                    { style: { fontSize: 12, color: "#94a3b8" } },
                    schemeDialog.nested
                      ? "Нет праймеров внутри продуктов 2-го раунда (назначьте продукты на раунд 2)."
                      : "Нет праймеров внутри выбранных продуктов."
                  )
            ),
            e(
              "div",
              { className: "field" },
              e("label", null, "Смесь / концентрации / примечания"),
              e("textarea", {
                value: schemeDialog.notes || "",
                placeholder: "Мастер-микс, концентрации праймеров, объёмы, программа термоциклера…",
                onChange: (ev) =>
                  setSchemeDialog(Object.assign({}, schemeDialog, { notes: ev.target.value })),
              })
            ),
            e(
              "div",
              { className: "annomap-modal-actions" },
              e(
                "button",
                { type: "button", className: "btn-secondary", onClick: () => setSchemeDialog(null) },
                "Отмена"
              ),
              e(
                "button",
                { type: "button", className: "btn-primary", onClick: saveSchemeDialog },
                "Сохранить схему"
              )
            )
          )
        ),
      e(
        "p",
        { className: "annomap-hint" },
        "Отметьте продукты галочками и нажмите «Создать схему ПЦР». Продукты остаются на общей карте и могут входить в несколько схем. Схемы сохраняются вместе с прогоном (кнопка «Сохранить прогон»)."
      )
    );
  }

  window.AnnotationMapView = AnnotationMapView;
})();
