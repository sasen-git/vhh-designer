#!/usr/bin/env python3
"""
VHH Compensation Analysis - IMGT Positions (v2)

What this script does
- Streams an annotated VHH CSV (fr1/cdr1/fr2/cdr2/fr3/cdr3/fr4 + family + hallmarks)
- Learns simple conditional residue rules: condition -> (IMGT position -> preferred AA)
- Outputs JSON rules compatible with the vhh_designer_v7 IMGT-rule schema:
    {condition, position: "IMGT##", suggested_aa, confidence, support, lift, baseline_prob, source}

Why v2
- Adds sanity checks to ensure the FR strings in the CSV are aligned to true IMGT boundaries
- Greatly reduces runtime/memory by:
    * restricting analysis to a configurable set of FR positions
    * using a smaller, higher-signal condition set by default
    * avoiding pandas iterrows()

USAGE (recommended starting point)
python vhh_compensation_imgt_v2.py \
  --csv data/databases/annotated/vhh_full_annotated_v3.csv \
  --output models/compensation/imgt_v1 \
  --chunk-size 200000 \
  --positions fr3fr4 \
  --conditions minimal \
  --min-support 5000 \
  --min-confidence 0.70 \
  --min-lift 1.15
"""

import os
import json
import re
import argparse
import pandas as pd
from collections import defaultdict
from datetime import datetime
import gc
from typing import Iterable, Tuple, List

# ----------------------------
# IMGT boundaries (VHH/VH domain)
# ----------------------------
IMGT_BOUNDARIES = {
    "FR1": (1, 26),      # 26
    "CDR1": (27, 38),    # 12
    "FR2": (39, 55),     # 17  (includes IMGT41 "W" at FR2 index 3)
    "CDR2": (56, 65),    # 10
    "FR3": (66, 104),    # 39
    "CDR3": (105, 117),  # 13 + insertions at 111.* in true IMGT; CSV CDR3 may be variable-length
    "FR4": (118, 128),   # 11
}

HALLMARK_POSITIONS = [42, 49, 50, 52]
VERNIER_POSITIONS = [42, 49, 50, 52, 66, 67, 68, 69, 71, 76, 78, 82, 87, 89, 91, 94]

EXPECTED_FIXED_LENGTHS = {
    "fr2": 17,
    "fr3": 39,
    "fr4": 11,
}
# CDR lengths are variable; fr1 may vary depending on extraction, and is not required for fr3fr4 mode.
# cdr3 is variable length


# ----------------------------
# Simple AA property helpers
# ----------------------------
HYDROPHOBICITY = {
    "I": 4.5, "V": 4.2, "L": 3.8, "F": 2.8, "C": 2.5, "M": 1.9, "A": 1.8,
    "G": -0.4, "T": -0.7, "S": -0.8, "W": -0.9, "Y": -1.3, "P": -1.6,
    "H": -3.2, "E": -3.5, "Q": -3.5, "D": -3.5, "N": -3.5, "K": -3.9, "R": -4.5,
}
CHARGE = {"K": 1, "R": 1, "H": 0.1, "D": -1, "E": -1}
for _aa in "ACFGILMNPQSTVWY":
    CHARGE[_aa] = 0


def get_cdr_features(cdr: str) -> dict:
    if not cdr:
        return {"length": 0, "charge": 0, "hydrophobicity": 0.0, "n_cys": 0}
    length = len(cdr)
    charge = sum(CHARGE.get(aa, 0) for aa in cdr)
    hydro = sum(HYDROPHOBICITY.get(aa, 0.0) for aa in cdr) / max(length, 1)
    return {
        "length": length,
        "charge": charge,
        "hydrophobicity": round(hydro, 2),
        "n_cys": cdr.count("C"),
    }


def categorize_cdr3(cdr3: str) -> dict:
    f = get_cdr_features(cdr3)
    L = f["length"]

    if L <= 8:
        len_cat = "short"
    elif L <= 14:
        len_cat = "medium"
    else:
        len_cat = "long"

    q = f["charge"]
    if q <= -2:
        charge_cat = "negative"
    elif q >= 2:
        charge_cat = "positive"
    else:
        charge_cat = "neutral"

    n = f["n_cys"]
    if n == 0:
        cys_cat = "no_cys"
    elif n == 1:
        cys_cat = "one_cys"
    else:
        cys_cat = "multi_cys"

    return {
        "cdr3_len": len_cat,
        "cdr3_charge": charge_cat,
        "cdr3_cys": cys_cat,
        **f,
    }


class StreamingCompensationAnalyzer:
    """
    Stores:
      counts[pos][condition][aa] = count
      cond_totals[pos][condition] = number of sequences counted for that (pos, condition)
      baseline[pos][aa] = unconditional count
    """
    def __init__(self, target_positions: Iterable[int]):
        self.target_positions = set(int(p) for p in target_positions)
        self.counts = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
        self.cond_totals = defaultdict(lambda: defaultdict(int))
        self.baseline = defaultdict(lambda: defaultdict(int))
        self.total_sequences_seen = 0
        self.total_sequences_used = 0

    def _iter_target_fr_residues(self, fr1: str, fr2: str, fr3: str, fr4: str) -> Iterable[Tuple[int, str]]:
        for i, aa in enumerate(fr1):
            pos = 1 + i
            if pos in self.target_positions:
                yield pos, aa
        for i, aa in enumerate(fr2):
            pos = 39 + i
            if pos in self.target_positions:
                yield pos, aa
        for i, aa in enumerate(fr3):
            pos = 66 + i
            if pos in self.target_positions:
                yield pos, aa
        for i, aa in enumerate(fr4):
            pos = 118 + i
            if pos in self.target_positions:
                yield pos, aa

    def _build_conditions(self, cdr1: str, cdr2: str, cdr3: str, family: str, hallmarks: str, mode: str) -> List[str]:
        cats = categorize_cdr3(cdr3)
        conds = [
            f"cdr3_len={cats['cdr3_len']}",
            f"cdr3_charge={cats['cdr3_charge']}",
            f"cdr3_cys={cats['cdr3_cys']}",
            f"family={family}",
            f"hallmarks={hallmarks}",
        ]

        if len(cdr3) >= 1:
            conds.append(f"cdr3[-1]={cdr3[-1]}")
        if len(cdr3) >= 2:
            conds.append(f"cdr3[-2]={cdr3[-2]}")
        if len(cdr3) >= 3:
            conds.append(f"cdr3[-3]={cdr3[-3]}")

        if mode == "full":
            for i, aa in enumerate(cdr1[:3]):
                conds.append(f"cdr1[{i}]={aa}")
            for i, aa in enumerate(cdr2[:2]):
                conds.append(f"cdr2[{i}]={aa}")
            for i, aa in enumerate(cdr3[:2]):
                conds.append(f"cdr3[{i}]={aa}")

        return conds

    def update(self, fr1: str, cdr1: str, fr2: str, cdr2: str, fr3: str, cdr3: str, fr4: str,
               family: str, hallmarks: str, condition_mode: str):
        self.total_sequences_seen += 1

        if not fr3 or not cdr3:
            return

        # Sanity check only the fixed framework regions we actually map to IMGT positions.
        # We do NOT require fixed CDR lengths (they are variable in IMGT).
        need_fr2 = any(39 <= p <= 55 for p in self.target_positions) or any(p in HALLMARK_POSITIONS for p in self.target_positions)
        need_fr3 = any(66 <= p <= 104 for p in self.target_positions)
        need_fr4 = any(118 <= p <= 128 for p in self.target_positions)

        if need_fr2 and (not fr2 or len(fr2) != EXPECTED_FIXED_LENGTHS['fr2']):
            return
        if need_fr3 and (not fr3 or len(fr3) != EXPECTED_FIXED_LENGTHS['fr3']):
            return
        if need_fr4 and (not fr4 or len(fr4) != EXPECTED_FIXED_LENGTHS['fr4']):
            return

        conds = self._build_conditions(cdr1, cdr2, cdr3, family, hallmarks, condition_mode)

        for pos, aa in self._iter_target_fr_residues(fr1, fr2, fr3, fr4):
            self.baseline[pos][aa] += 1
            for cond in conds:
                self.counts[pos][cond][aa] += 1
                self.cond_totals[pos][cond] += 1

        self.total_sequences_used += 1

    def extract_rules(self, min_support: int, min_confidence: float, min_lift: float, source: str) -> List[dict]:
        rules = []
        for pos, cond_dict in self.counts.items():
            base_counts = self.baseline[pos]
            base_total = sum(base_counts.values())
            if base_total < min_support:
                continue

            for cond, aa_counts in cond_dict.items():
                cond_total = self.cond_totals[pos][cond]
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
                })

        rules.sort(key=lambda r: (-r["support"], -r["confidence"], -r["lift"] if r["lift"] != float("inf") else -1e9))
        return rules

    def vernier_archetypes(self, min_support: int) -> dict:
        arch = defaultdict(dict)
        for pos in VERNIER_POSITIONS:
            if pos not in self.counts:
                continue
            for cond, aa_counts in self.counts[pos].items():
                if not (cond.startswith("family=") or cond.startswith("hallmarks=")):
                    continue
                total = sum(aa_counts.values())
                if total < min_support:
                    continue
                aa, c = max(aa_counts.items(), key=lambda kv: kv[1])
                arch[cond][f"IMGT{pos}"] = {"aa": aa, "count": int(c), "frequency": round(c / total, 3)}
        return dict(arch)


def choose_positions(mode: str) -> List[int]:
    if mode == "vernier":
        return sorted(set(VERNIER_POSITIONS))
    if mode == "fr3fr4":
        fr3 = list(range(IMGT_BOUNDARIES["FR3"][0], IMGT_BOUNDARIES["FR3"][1] + 1))
        fr4 = list(range(IMGT_BOUNDARIES["FR4"][0], IMGT_BOUNDARIES["FR4"][1] + 1))
        return sorted(set(fr3 + fr4 + HALLMARK_POSITIONS))
    if mode == "all":
        ps = []
        for region, (s, e) in IMGT_BOUNDARIES.items():
            if region.startswith("FR"):
                ps.extend(range(s, e + 1))
        return sorted(set(ps))
    raise ValueError(f"Unknown positions mode: {mode}")


def process_csv_streaming(csv_path: str, output_dir: str, chunk_size: int,
                          positions_mode: str, condition_mode: str,
                          min_support: int, min_confidence: float, min_lift: float):
    os.makedirs(output_dir, exist_ok=True)
    target_positions = choose_positions(positions_mode)
    analyzer = StreamingCompensationAnalyzer(target_positions)

    print(f"Input: {csv_path}")
    print(f"Output: {output_dir}")
    print(f"Chunk size: {chunk_size:,}")
    print(f"Positions: {positions_mode} ({len(target_positions)} IMGT positions)")
    print(f"Conditions: {condition_mode}")
    print(f"Thresholds: support>={min_support}, confidence>={min_confidence}, lift>={min_lift}")
    print()

    total_rows = 0
    for chunk_idx, chunk in enumerate(pd.read_csv(csv_path, chunksize=chunk_size, low_memory=False), start=1):
        # Normalize column names to allow FR1/FR1... uppercase headers, spaces, etc.
        def _norm_col(c: str) -> str:
            c = str(c).strip().lower()
            c = re.sub(r"[^0-9a-zA-Z_]+", "_", c)
            return c
        chunk.rename(columns={c: _norm_col(c) for c in chunk.columns}, inplace=True)

        total_rows += len(chunk)
        used_before = analyzer.total_sequences_used

        for row in chunk.itertuples(index=False):
            fr1 = str(getattr(row, "fr1", "") or "")
            cdr1 = str(getattr(row, "cdr1", "") or "")
            fr2 = str(getattr(row, "fr2", "") or "")
            cdr2 = str(getattr(row, "cdr2", "") or "")
            fr3 = str(getattr(row, "fr3", "") or "")
            cdr3 = str(getattr(row, "cdr3", "") or "")
            fr4 = str(getattr(row, "fr4", "") or "")
            family = str(getattr(row, "family", "Unknown") or "Unknown")
            hallmarks = str(getattr(row, "hallmarks", "----") or "----")

            analyzer.update(fr1, cdr1, fr2, cdr2, fr3, cdr3, fr4, family, hallmarks, condition_mode)

        used_now = analyzer.total_sequences_used
        print(f"Chunk {chunk_idx}: rows={total_rows:,} | used={used_now:,} (+{used_now-used_before:,}) | seen={analyzer.total_sequences_seen:,}")
        gc.collect()

    print("\\n" + "=" * 70)
    print("ANALYSIS COMPLETE")
    print("=" * 70)
    print(f"Seen: {analyzer.total_sequences_seen:,} | Used (passed sanity): {analyzer.total_sequences_used:,}")

    rules = analyzer.extract_rules(
        min_support=min_support,
        min_confidence=min_confidence,
        min_lift=min_lift,
        source="compensation_imgt_v2",
    )
    archetypes = analyzer.vernier_archetypes(min_support=max(10000, min_support))

    rules_path = os.path.join(output_dir, "compensation_imgt_rules.json")
    with open(rules_path, "w") as f:
        json.dump(rules, f, indent=2)
    arch_path = os.path.join(output_dir, "vernier_archetypes_imgt.json")
    with open(arch_path, "w") as f:
        json.dump(archetypes, f, indent=2)

    print(f"Saved rules: {rules_path} ({len(rules)} rules)")
    print(f"Saved archetypes: {arch_path} ({len(archetypes)} conditions)")

    print("\\nTop rules:")
    for r in rules[:20]:
        print(f"  {r['position']:7s} -> {r['suggested_aa']} | {r['condition'][:32]:32s} | conf={r['confidence']:.2f} lift={r['lift']} n={r['support']:,}")

    return rules, archetypes


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", required=True, help="Annotated CSV (vhh_full_annotated_v3.csv)")
    ap.add_argument("-o", "--output", default="models/compensation/imgt_v2", help="Output directory")
    ap.add_argument("--chunk-size", type=int, default=200000, help="Rows per chunk")
    ap.add_argument("--positions", choices=["vernier", "fr3fr4", "all"], default="fr3fr4",
                    help="Which FR positions to learn rules for")
    ap.add_argument("--conditions", choices=["minimal", "full"], default="minimal",
                    help="How many CDR conditions to include")
    ap.add_argument("--min-support", type=int, default=5000)
    ap.add_argument("--min-confidence", type=float, default=0.70)
    ap.add_argument("--min-lift", type=float, default=1.15)
    args = ap.parse_args()

    print("=" * 70)
    print("VHH COMPENSATION ANALYSIS - IMGT POSITIONS (v2)")
    print("=" * 70)
    start = datetime.now()
    process_csv_streaming(
        csv_path=args.csv,
        output_dir=args.output,
        chunk_size=args.chunk_size,
        positions_mode=args.positions,
        condition_mode=args.conditions,
        min_support=args.min_support,
        min_confidence=args.min_confidence,
        min_lift=args.min_lift,
    )
    print(f"\\nTotal time: {datetime.now() - start}")
    print("Done.")


if __name__ == "__main__":
    main()