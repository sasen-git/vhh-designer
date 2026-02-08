#!/usr/bin/env python3
"""
Correlation Results v3 (CSV + IMGT) - Complete Breakdown
=======================================================

Goal
- Read the *dedup* annotated CSV with IMGT columns (imgt_##) and CDR strings.
- Produce:
  (1) correlation_results_v3_summary.json  (your "Complete Breakdown")
  (2) correlation_rules_v3.json           (compensation-style rules JSON)

Why this exists
- Your older correlations script is NPZ-based and extracts FRs via substring search,
  which is not IMGT-position-safe.
- This version treats IMGT columns as the coordinate system and emits the same rule
  schema style as the compensation pipeline:
    {condition, position, suggested_aa, confidence, support, lift, baseline_prob, source}

Outputs
- correlation_results_v3_summary.json:
  {
    "meta": {...},
    "germline_counts": {...},
    "hallmark_distributions": {...},
    "family_cdr_dependent": {...},
    "family_conservation": {...},
    "family_cdr_length_stats": {...},
    "indel_correlations": {...},
    "family_triplet_rules": {...}
  }

- correlation_rules_v3.json:
  list[ {condition, position, suggested_aa, confidence, support, lift, baseline_prob, source, rule_type} ]

Notes on performance
- Full-pass summaries are cheap (family counts, hallmark distributions, length stats).
- Rule mining can explode combinatorially; this script uses deterministic sampling for rules by default.
  Set --sample-per-million 1000000 to use all rows (slow).

USAGE
  python analyze_cdr_framework_correlations_imgt_csv_v3.py \
    --csv data/.../vhh_full_annotated_imgt_dedup.csv \
    -o models/correlations/imgt_csv_v3 \
    --chunk-size 200000 \
    --positions fr3fr4 \
    --conditions minimal \
    --min-support 5000 \
    --min-confidence 0.70 \
    --min-lift 1.15 \
    --sample-per-million 200000

"""

import os
import csv
import json
import math
import argparse
import hashlib
from collections import Counter, defaultdict
from datetime import datetime
from typing import Dict, List, Tuple, Optional

# -----------------------------
# IMGT boundaries (VH/VHH domain)
# -----------------------------
IMGT_BOUNDARIES = {
    "FR1": (1, 26),
    "CDR1": (27, 38),
    "FR2": (39, 55),
    "CDR2": (56, 65),
    "FR3": (66, 104),
    "CDR3": (105, 117),
    "FR4": (118, 128),
}

HALLMARK_POSITIONS = [42, 49, 50, 52]

# The "vernier-ish" subset used in your compensation script (good default for faster rule mining)
DEFAULT_VERNIER_POSITIONS = [42, 49, 50, 52, 66, 67, 68, 69, 71, 76, 78, 82, 87, 89, 91, 94]

AA_POS = set("KRH")
AA_NEG = set("DE")

def norm_cell(x: Optional[str]) -> str:
    """Normalize CSV cell to a residue-like token."""
    if x is None:
        return "-"
    x = str(x).strip()
    if x == "" or x.lower() in {"na", "nan", "none"}:
        return "-"
    return x

def charge(seq: str) -> int:
    seq = seq or ""
    pos = sum(1 for c in seq if c in AA_POS)
    neg = sum(1 for c in seq if c in AA_NEG)
    return pos - neg

def bin_int(x: int, bins: List[Tuple[int, int, str]]) -> str:
    for lo, hi, name in bins:
        if lo <= x <= hi:
            return name
    return "other"

def bin_len(n: int) -> str:
    # tweak as you like
    if n <= 8: return "short"
    if n <= 14: return "mid"
    if n <= 21: return "long"
    return "xlong"

def deterministic_sample(key: str, sample_per_million: int) -> bool:
    """Deterministic sampling using md5(key) mod 1e6."""
    if sample_per_million >= 1_000_000:
        return True
    h = hashlib.md5(key.encode("utf-8")).hexdigest()[:8]
    v = int(h, 16) % 1_000_000
    return v < sample_per_million

def choose_positions(mode: str, available_imgt_positions: List[int]) -> List[int]:
    avail = set(available_imgt_positions)
    if mode == "vernier":
        return [p for p in DEFAULT_VERNIER_POSITIONS if p in avail]
    if mode == "fr3fr4":
        ps = []
        ps.extend(range(IMGT_BOUNDARIES["FR3"][0], IMGT_BOUNDARIES["FR3"][1] + 1))
        ps.extend(range(IMGT_BOUNDARIES["FR4"][0], IMGT_BOUNDARIES["FR4"][1] + 1))
        ps.extend(HALLMARK_POSITIONS)
        return [p for p in sorted(set(ps)) if p in avail]
    if mode == "all_available":
        return sorted(avail)
    raise ValueError(f"Unknown positions mode: {mode}")

def build_conditions(row: Dict[str, str], mode: str, family: str, hallmarks: str) -> List[str]:
    """
    Return a list of string conditions.
    IMPORTANT: to keep compatibility with your compensation designer, conditions are strings.
    """
    cdr1 = row.get("cdr1", "") or ""
    cdr2 = row.get("cdr2", "") or ""
    cdr3 = row.get("cdr3", "") or ""

    conds = []
    # Always include family + hallmarks as conditions (useful for summaries / archetypes)
    conds.append(f"family={family}")
    if hallmarks:
        conds.append(f"hallmarks={hallmarks}")

    # Minimal condition set (fast + high-signal)
    def safe_at(s: str, i: int) -> str:
        try:
            return s[i]
        except Exception:
            return "-"

    if mode == "minimal":
        conds.extend([
            f"cdr1[0]={safe_at(cdr1, 0)}",
            f"cdr2[0]={safe_at(cdr2, 0)}",
            f"cdr3[-1]={safe_at(cdr3, -1)}",
            f"cdr3_len_bin={bin_len(len(cdr3))}",
            f"cdr3_cys_bin={bin_int((cdr3 or '').count('C'), [(0,0,'0'), (1,1,'1'), (2,2,'2'), (3,99,'3plus')])}",
            f"cdr3_charge_bin={bin_int(charge(cdr3), [(-99,-3,'neg3plus'), (-2,-1,'neg'), (0,0,'0'), (1,2,'pos'), (3,99,'pos3plus')])}",
        ])
        return conds

    if mode == "full":
        # a few extra boundary residues and lengths/charges for CDR1/2 as well
        conds.extend([
            f"cdr1[0]={safe_at(cdr1, 0)}", f"cdr1[1]={safe_at(cdr1, 1)}", f"cdr1[-1]={safe_at(cdr1, -1)}",
            f"cdr2[0]={safe_at(cdr2, 0)}", f"cdr2[1]={safe_at(cdr2, 1)}", f"cdr2[-1]={safe_at(cdr2, -1)}",
            f"cdr3[0]={safe_at(cdr3, 0)}", f"cdr3[1]={safe_at(cdr3, 1)}", f"cdr3[-1]={safe_at(cdr3, -1)}",
            f"cdr1_len={len(cdr1)}", f"cdr2_len={len(cdr2)}", f"cdr3_len={len(cdr3)}",
            f"cdr1_charge={charge(cdr1)}", f"cdr2_charge={charge(cdr2)}", f"cdr3_charge={charge(cdr3)}",
            f"cdr3_len_bin={bin_len(len(cdr3))}",
            f"cdr3_cys_bin={bin_int((cdr3 or '').count('C'), [(0,0,'0'), (1,1,'1'), (2,2,'2'), (3,99,'3plus')])}",
            f"cdr3_charge_bin={bin_int(charge(cdr3), [(-99,-3,'neg3plus'), (-2,-1,'neg'), (0,0,'0'), (1,2,'pos'), (3,99,'pos3plus')])}",
        ])
        return conds

    raise ValueError(f"Unknown conditions mode: {mode}")

def infer_hallmarks(imgt: Dict[int, str]) -> str:
    # e.g. "FERG"
    return "".join(imgt.get(p, "-") for p in HALLMARK_POSITIONS)

def infer_family(imgt: Dict[int, str], aa_v_full: str) -> str:
    """
    Family heuristic aligned to your earlier labels:
    - VH_like: IMGT50 == L
    - F_C2 / F_C4 / Y_C2 / Y_C4 based on IMGT42 and cysteine count
    - Other_VHH otherwise

    NOTE: cysteine logic is adjustable. By default:
      C4 if total cysteines in aa_v_full >= 4 else C2
    """
    p42 = imgt.get(42, "-")
    p50 = imgt.get(50, "-")

    if p50 == "L":
        return "VH_like"

    cys_total = (aa_v_full or "").count("C")
    c_label = "C4" if cys_total >= 4 else "C2"

    if p42 in {"F", "Y"}:
        return f"{p42}_{c_label}"

    return "Other_VHH"

def update_mean_var(stats: Dict, key: Tuple, x: int):
    s = stats.setdefault(key, {"n": 0, "sum": 0.0, "sumsq": 0.0})
    s["n"] += 1
    s["sum"] += x
    s["sumsq"] += x * x

def finalize_mean_sd(s: Dict) -> Tuple[float, float]:
    n = s.get("n", 0)
    if n <= 1:
        return (float(s.get("sum", 0.0)) / max(1, n), 0.0)
    mean = s["sum"] / n
    var = max(0.0, (s["sumsq"] / n) - mean * mean)
    return (mean, math.sqrt(var))

def extract_rules_from_counts(
    counts_by_pos_cond: Dict[int, Dict[str, Counter]],
    cond_totals_by_pos: Dict[int, Counter],
    baseline_by_pos: Dict[int, Counter],
    min_support: int,
    min_confidence: float,
    min_lift: float,
    source: str,
    rule_type: str,
) -> List[dict]:
    """
    Same math as compensation:
      conf = P(aa | cond, pos)
      baseline_prob = P(aa | pos)
      lift = conf / baseline_prob
    """
    rules = []
    for pos, cond_map in counts_by_pos_cond.items():
        base_counts = baseline_by_pos.get(pos, Counter())
        base_total = sum(base_counts.values())
        if base_total < min_support:
            continue

        for cond, aa_counts in cond_map.items():
            cond_total = cond_totals_by_pos[pos].get(cond, 0)
            if cond_total < min_support:
                continue

            aa, count = max(aa_counts.items(), key=lambda kv: kv[1])
            if count < min_support:
                continue

            conf = count / cond_total
            if conf < min_confidence:
                continue

            base_prob = base_counts.get(aa, 0) / base_total if base_total else 0.0
            lift = (conf / base_prob) if base_prob > 0 else float("inf")
            if lift < min_lift:
                continue

            rules.append({
                "condition": cond,
                "position": f"IMGT{pos}",
                "suggested_aa": aa,
                "confidence": round(conf, 3),
                "support": int(count),
                "lift": round(float(lift), 2) if lift != float("inf") else float("inf"),
                "baseline_prob": round(base_prob, 3),
                "source": source,
                "rule_type": rule_type,
            })

    rules.sort(key=lambda r: (-r["support"], -r["confidence"], -(r["lift"] if r["lift"] != float("inf") else 1e9)))
    return rules

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", required=True)
    ap.add_argument("-o", "--output", required=True)
    ap.add_argument("--chunk-size", type=int, default=200000, help="Progress print frequency only (streaming csv).")
    ap.add_argument("--positions", choices=["vernier", "fr3fr4", "all_available"], default="fr3fr4")
    ap.add_argument("--conditions", choices=["minimal", "full"], default="minimal")
    ap.add_argument("--min-support", type=int, default=5000)
    ap.add_argument("--min-confidence", type=float, default=0.70)
    ap.add_argument("--min-lift", type=float, default=1.15)
    ap.add_argument("--sample-per-million", type=int, default=200000,
                    help="Deterministic sampling for rule mining (0..1,000,000). Summaries use all rows.")
    ap.add_argument("--max-rows", type=int, default=0, help="Debug: stop after N rows (0=all).")
    args = ap.parse_args()

    os.makedirs(args.output, exist_ok=True)

    # -----------------------------
    # PASS: stream CSV
    # -----------------------------
    total = 0
    used_for_rules = 0

    # Summary outputs
    source_counts = Counter()
    family_counts = Counter()

    hallmark_pos_counts = {p: Counter() for p in HALLMARK_POSITIONS}
    hallmark_pattern_counts = Counter()

    # family -> per hallmark pos -> counter
    family_hallmark_pos_counts = defaultdict(lambda: {p: Counter() for p in HALLMARK_POSITIONS})

    # CDR length stats per family
    length_stats = {}  # key=(family, cdr_name) -> {n,sum,sumsq}

    # Indel/missing stats (use IMGT grid for FR3/FR4 gaps)
    indel_stats = Counter()  # keys like ("F_C2","fr4_missing",0) etc.

    # Triplet rules counts per family
    # (family, boundary, lhs3, rhs3) -> count
    triplet_counts = Counter()

    # RULE MINING COUNTS (sampled)
    # We mine "CDR-dependent FR positions" using the same counting math as compensation,
    # but conditions include CDR features (cdr3_len_bin, cdr1[0], etc.)
    #
    # Structure:
    #   family -> pos -> cond -> Counter(aa)
    rule_counts = defaultdict(lambda: defaultdict(lambda: defaultdict(Counter)))
    rule_cond_totals = defaultdict(lambda: defaultdict(Counter))   # family -> pos -> Counter(cond -> total)
    rule_baseline = defaultdict(lambda: defaultdict(Counter))      # family -> pos -> Counter(aa -> total)

    # Conservation (sampled): family -> pos -> Counter(aa)
    conservation = defaultdict(lambda: defaultdict(Counter))

    # detect IMGT columns present
    with open(args.csv, "r", newline="") as fh:
        reader = csv.DictReader(fh)
        header = reader.fieldnames or []
        imgt_positions = []
        for h in header:
            if h.startswith("imgt_"):
                try:
                    imgt_positions.append(int(h.split("_")[1]))
                except Exception:
                    pass
        imgt_positions = sorted(set(imgt_positions))
        target_positions = choose_positions(args.positions, imgt_positions)

        needed_cols = set(["aa_v_full", "cdr1", "cdr2", "cdr3", "fr3", "fr4", "source", "family", "hallmarks"])
        present = {c for c in needed_cols if c in header}

        print("=" * 70)
        print("CORRELATION RESULTS v3 (CSV + IMGT)")
        print("=" * 70)
        print(f"Input:  {args.csv}")
        print(f"Output: {args.output}")
        print(f"IMGT columns detected: {len(imgt_positions)}")
        print(f"Target IMGT positions ({args.positions}): {len(target_positions)}")
        print(f"Sampling for rule mining: {args.sample_per_million}/1,000,000")
        print(f"Conditions: {args.conditions}")
        print(f"Thresholds: support>={args.min_support}, confidence>={args.min_confidence}, lift>={args.min_lift}")
        print(f"Columns present: {sorted(present)}")
        print()

        for row in reader:
            total += 1
            if args.max_rows and total > args.max_rows:
                break

            src = row.get("source", "")
            if src:
                source_counts[src] += 1

            aa_v_full = row.get("aa_v_full", "") or ""

            # IMGT dict for this row (only the positions we care about)
            imgt = {}
            for p in HALLMARK_POSITIONS:
                col = f"imgt_{p}"
                if col in row:
                    imgt[p] = norm_cell(row.get(col))
                else:
                    imgt[p] = "-"

            # hallmark distributions
            hallmarks = row.get("hallmarks", "") or infer_hallmarks(imgt)
            hallmark_pattern_counts[hallmarks] += 1
            for p in HALLMARK_POSITIONS:
                hallmark_pos_counts[p][imgt.get(p, "-")] += 1

            # family
            fam = row.get("family", "") or infer_family(imgt, aa_v_full)
            family_counts[fam] += 1
            for p in HALLMARK_POSITIONS:
                family_hallmark_pos_counts[fam][p][imgt.get(p, "-")] += 1

            # CDR length stats
            cdr1 = row.get("cdr1", "") or ""
            cdr2 = row.get("cdr2", "") or ""
            cdr3 = row.get("cdr3", "") or ""
            update_mean_var(length_stats, (fam, "CDR1"), len(cdr1))
            update_mean_var(length_stats, (fam, "CDR2"), len(cdr2))
            update_mean_var(length_stats, (fam, "CDR3"), len(cdr3))

            # Indel / gap counts from FR3/FR4 IMGT grid (if present)
            fr3_missing = 0
            fr4_missing = 0
            for p in range(IMGT_BOUNDARIES["FR3"][0], IMGT_BOUNDARIES["FR3"][1] + 1):
                col = f"imgt_{p}"
                if col in row and norm_cell(row.get(col)) == "-":
                    fr3_missing += 1
            for p in range(IMGT_BOUNDARIES["FR4"][0], IMGT_BOUNDARIES["FR4"][1] + 1):
                col = f"imgt_{p}"
                if col in row and norm_cell(row.get(col)) == "-":
                    fr4_missing += 1

            indel_stats[(fam, "fr3_missing", fr3_missing)] += 1
            indel_stats[(fam, "fr4_missing", fr4_missing)] += 1
            indel_stats[(fam, "cdr3_len", len(cdr3))] += 1

            # Triplet boundary rules (do on *all* rows; cheap)
            fr3 = row.get("fr3", "") or ""
            fr4 = row.get("fr4", "") or ""
            if len(cdr3) >= 3 and len(fr4) >= 3:
                lhs = cdr3[-3:]
                rhs = fr4[:3]
                triplet_counts[(fam, "cdr3_suffix3_to_fr4_prefix3", lhs, rhs)] += 1
            if len(cdr2) >= 3 and len(fr3) >= 3:
                lhs = cdr2[-3:]
                rhs = fr3[:3]
                triplet_counts[(fam, "cdr2_suffix3_to_fr3_prefix3", lhs, rhs)] += 1

            # -----------------------
            # RULE MINING (sampled)
            # -----------------------
            # Use aa_v_full as deterministic sampling key
            if deterministic_sample(aa_v_full if aa_v_full else f"row{total}", args.sample_per_million):
                used_for_rules += 1

                # build IMGT dict for target positions
                imgt_t = {}
                for p in target_positions:
                    col = f"imgt_{p}"
                    if col in row:
                        imgt_t[p] = norm_cell(row.get(col))
                    else:
                        imgt_t[p] = "-"

                # conditions for this row
                conds = build_conditions(row, args.conditions, fam, hallmarks)

                # update baseline + per-condition counts, per position
                # (this is the compensation counting pattern, but grouped by family)
                for p in target_positions:
                    aa = imgt_t.get(p, "-")
                    if aa in {"-", "X", "*"}:
                        continue
                    rule_baseline[fam][p][aa] += 1
                    for cond in conds:
                        rule_counts[fam][p][cond][aa] += 1
                        rule_cond_totals[fam][p][cond] += 1

                    # conservation counts
                    conservation[fam][p][aa] += 1

            if total % args.chunk_size == 0:
                print(f"Rows scanned: {total:,} | sampled for rules: {used_for_rules:,}", flush=True)

    # -----------------------------
    # POST: build v3 breakdown
    # -----------------------------
    # Family distribution
    total_n = sum(family_counts.values()) or 1
    germline_counts = {k: {"count": int(v), "pct": round(100*v/total_n, 2)} for k, v in family_counts.most_common()}

    # Hallmark distributions
    hallmark_distributions = {
        "positions": {},
        "patterns_top": []
    }
    for p in HALLMARK_POSITIONS:
        tot = sum(hallmark_pos_counts[p].values()) or 1
        top = hallmark_pos_counts[p].most_common(10)
        hallmark_distributions["positions"][f"imgt_{p}"] = [
            {"aa": aa, "count": int(c), "pct": round(100*c/tot, 2)} for aa, c in top
        ]
    # top patterns
    for pat, c in hallmark_pattern_counts.most_common(20):
        hallmark_distributions["patterns_top"].append({"pattern": pat, "count": int(c)})

    # Family conservation (consensus per family for available target positions)
    family_conservation = {}
    for fam, pos_map in conservation.items():
        fam_out = {"consensus": {}, "conservation_pct": {}}
        for p, aa_counts in pos_map.items():
            tot = sum(aa_counts.values()) or 1
            aa, c = max(aa_counts.items(), key=lambda kv: kv[1])
            fam_out["consensus"][f"IMGT{p}"] = aa
            fam_out["conservation_pct"][f"IMGT{p}"] = round(100*c/tot, 2)
        family_conservation[fam] = fam_out

    # CDR length stats
    family_cdr_length_stats = {}
    for (fam, cdr_name), s in length_stats.items():
        mean, sd = finalize_mean_sd(s)
        family_cdr_length_stats.setdefault(fam, {})
        family_cdr_length_stats[fam][cdr_name] = {
            "mean": round(mean, 3),
            "sd": round(sd, 3),
            "n": int(s["n"]),
        }

    # Extract compensation-style rules per family
    all_rules = []
    family_cdr_dependent = {}

    for fam in rule_counts.keys():
        rules = extract_rules_from_counts(
            counts_by_pos_cond=rule_counts[fam],
            cond_totals_by_pos=rule_cond_totals[fam],
            baseline_by_pos=rule_baseline[fam],
            min_support=args.min_support,
            min_confidence=args.min_confidence,
            min_lift=args.min_lift,
            source="correlation_imgt_csv_v3",
            rule_type="cdr_to_imgt",
        )
        all_rules.extend(rules)

        # Summarize "how many FR positions have at least one CDR-dependent rule"
        pos_with_rules = set(r["position"] for r in rules if r["condition"].startswith("cdr") or "cdr" in r["condition"])
        family_cdr_dependent[fam] = {
            "n_rules": len(rules),
            "n_positions_with_rules": len(pos_with_rules),
            "top_rules": rules[:50],
        }

    # Triplet rules: convert counts to "rules" (confidence etc.)
    # For each (fam, boundary, lhs) pick most frequent rhs and compute confidence
    family_triplet_rules = {}
    triplet_rules_out = []
    triplet_index = defaultdict(lambda: Counter())  # (fam,boundary,lhs)->Counter(rhs)
    for (fam, boundary, lhs, rhs), c in triplet_counts.items():
        triplet_index[(fam, boundary, lhs)][rhs] += c

    for (fam, boundary, lhs), rhs_counts in triplet_index.items():
        total_lhs = sum(rhs_counts.values())
        if total_lhs < args.min_support:
            continue
        rhs_top, c_top = max(rhs_counts.items(), key=lambda kv: kv[1])
        conf = c_top / total_lhs
        if conf < args.min_confidence:
            continue

        # baseline rhs probability within this family+boundary (approx)
        # Compute baseline from aggregated counts for fam+boundary
        base_counter = Counter()
        for (f2, b2, l2), rc in triplet_index.items():
            if f2 == fam and b2 == boundary:
                base_counter.update(rc)
        base_total = sum(base_counter.values()) or 1
        base_prob = base_counter.get(rhs_top, 0) / base_total
        lift = (conf / base_prob) if base_prob > 0 else float("inf")
        if lift < args.min_lift:
            continue

        rule = {
            "condition": f"family={fam} AND triplet:{boundary}:{lhs}",
            "position": f"TRIPLET:{boundary}",
            "suggested_aa": rhs_top,  # still string; here it's a 3-mer
            "confidence": round(conf, 3),
            "support": int(c_top),
            "lift": round(float(lift), 2) if lift != float("inf") else float("inf"),
            "baseline_prob": round(base_prob, 3),
            "source": "correlation_imgt_csv_v3",
            "rule_type": "triplet",
        }
        triplet_rules_out.append(rule)

    triplet_rules_out.sort(key=lambda r: (-r["confidence"], -r["support"]))
    for fam, _ in family_counts.items():
        family_triplet_rules[fam] = {
            "n_rules": sum(1 for r in triplet_rules_out if f"family={fam} " in r["condition"]),
            "top_rules": [r for r in triplet_rules_out if f"family={fam} " in r["condition"]][:50]
        }

    all_rules.extend(triplet_rules_out)

    # Indel correlations summary (distributions; you can add rule mining here later)
    indel_correlations = {"fr3_missing": {}, "fr4_missing": {}}
    for fam, _ in family_counts.items():
        fr3 = Counter()
        fr4 = Counter()
        for (f2, kind, val), c in indel_stats.items():
            if f2 != fam:
                continue
            if kind == "fr3_missing":
                fr3[val] += c
            elif kind == "fr4_missing":
                fr4[val] += c
        indel_correlations["fr3_missing"][fam] = dict(fr3.most_common(20))
        indel_correlations["fr4_missing"][fam] = dict(fr4.most_common(20))

    # -----------------------------
    # WRITE OUTPUTS
    # -----------------------------
    summary = {
        "meta": {
            "generated_at": datetime.now().isoformat(timespec="seconds"),
            "input_csv": args.csv,
            "total_sequences": int(total),
            "sampled_for_rules": int(used_for_rules),
            "sample_per_million": int(args.sample_per_million),
            "positions_mode": args.positions,
            "conditions_mode": args.conditions,
            "thresholds": {
                "min_support": args.min_support,
                "min_confidence": args.min_confidence,
                "min_lift": args.min_lift,
            },
            "sources_top": source_counts.most_common(20),
        },

        # 1) Family distribution
        "germline_counts": germline_counts,

        # 2) Hallmarks
        "hallmark_distributions": hallmark_distributions,

        # 3) CDR-dependent FR positions (rules)
        "family_cdr_dependent": family_cdr_dependent,

        # 4) Family conservation
        "family_conservation": family_conservation,

        # 5) CDR length stats
        "family_cdr_length_stats": family_cdr_length_stats,

        # 6) NEW: indel distributions (extend to rules as needed)
        "indel_correlations": indel_correlations,

        # 7) Triplet rules
        "family_triplet_rules": family_triplet_rules,
    }

    summary_path = os.path.join(args.output, "correlation_results_v3_summary.json")
    rules_path = os.path.join(args.output, "correlation_rules_v3.json")

    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)

    with open(rules_path, "w") as f:
        json.dump(all_rules, f, indent=2)

    print("\nDone.")
    print(f"Summary JSON: {summary_path}")
    print(f"Rules JSON:   {rules_path}")
    print(f"Total rules:  {len(all_rules):,}")

if __name__ == "__main__":
    main()
