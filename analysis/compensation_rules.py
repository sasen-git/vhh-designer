#!/usr/bin/env python3
import os, csv, json, re, argparse, hashlib
from collections import defaultdict
from datetime import datetime
from typing import Iterable, Tuple, List, Dict, Optional

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
VERNIER_POSITIONS = [42, 49, 50, 52, 66, 67, 68, 69, 71, 76, 78, 82, 87, 89, 91, 94]

HYDROPHOBICITY = {
    "I": 4.5, "V": 4.2, "L": 3.8, "F": 2.8, "C": 2.5, "M": 1.9, "A": 1.8,
    "G": -0.4, "T": -0.7, "S": -0.8, "W": -0.9, "Y": -1.3, "P": -1.6,
    "H": -3.2, "E": -3.5, "Q": -3.5, "D": -3.5, "N": -3.5, "K": -3.9, "R": -4.5,
}
CHARGE = {"K": 1, "R": 1, "H": 0.1, "D": -1, "E": -1}
for _aa in "ACFGILMNPQSTVWY":
    CHARGE[_aa] = 0

def _norm_col(c: str) -> str:
    c = str(c).strip().lower()
    c = re.sub(r"[^0-9a-zA-Z_]+", "_", c)
    return c

def norm_cell(x: Optional[str]) -> str:
    if x is None:
        return "-"
    x = str(x).strip()
    if x == "" or x.lower() in {"na", "nan", "none"}:
        return "-"
    return x

def get_cdr_features(cdr: str) -> dict:
    if not cdr:
        return {"length": 0, "charge": 0, "hydrophobicity": 0.0, "n_cys": 0}
    length = len(cdr)
    charge = sum(CHARGE.get(aa, 0) for aa in cdr)
    hydro = sum(HYDROPHOBICITY.get(aa, 0.0) for aa in cdr) / max(length, 1)
    return {"length": length, "charge": charge, "hydrophobicity": round(hydro, 2), "n_cys": cdr.count("C")}

def categorize_cdr3(cdr3: str) -> dict:
    f = get_cdr_features(cdr3)
    L = f["length"]
    if L <= 8: len_cat = "short"
    elif L <= 14: len_cat = "medium"
    else: len_cat = "long"

    q = f["charge"]
    if q <= -2: charge_cat = "negative"
    elif q >= 2: charge_cat = "positive"
    else: charge_cat = "neutral"

    n = f["n_cys"]
    if n == 0: cys_cat = "no_cys"
    elif n == 1: cys_cat = "one_cys"
    else: cys_cat = "multi_cys"

    return {"cdr3_len": len_cat, "cdr3_charge": charge_cat, "cdr3_cys": cys_cat, **f}

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

def infer_hallmarks(imgt: Dict[int, str]) -> str:
    return "".join(imgt.get(p, "-") for p in HALLMARK_POSITIONS)

def infer_family(imgt: Dict[int, str], aa_v_full: str) -> str:
    p42 = imgt.get(42, "-")
    p50 = imgt.get(50, "-")
    if p50 == "L":
        return "VH_like"
    cys_total = (aa_v_full or "").count("C")
    c_label = "C4" if cys_total >= 4 else "C2"
    if p42 in {"F", "Y"}:
        return f"{p42}_{c_label}"
    return "Other_VHH"

class StreamingCompensationAnalyzer:
    def __init__(self, target_positions: Iterable[int]):
        self.target_positions = list(sorted(set(int(p) for p in target_positions)))
        self.counts = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
        self.cond_totals = defaultdict(lambda: defaultdict(int))
        self.baseline = defaultdict(lambda: defaultdict(int))
        self.total_sequences_seen = 0
        self.total_sequences_used = 0

    def _build_conditions(self, cdr1: str, cdr2: str, cdr3: str, family: str, hallmarks: str, mode: str) -> List[str]:
        cats = categorize_cdr3(cdr3)
        conds = [
            f"cdr3_len={cats['cdr3_len']}",
            f"cdr3_charge={cats['cdr3_charge']}",
            f"cdr3_cys={cats['cdr3_cys']}",
            f"family={family}",
            f"hallmarks={hallmarks}",
        ]
        if len(cdr3) >= 1: conds.append(f"cdr3[-1]={cdr3[-1]}")
        if len(cdr3) >= 2: conds.append(f"cdr3[-2]={cdr3[-2]}")
        if len(cdr3) >= 3: conds.append(f"cdr3[-3]={cdr3[-3]}")

        if mode == "full":
            for i, aa in enumerate((cdr1 or "")[:3]): conds.append(f"cdr1[{i}]={aa}")
            for i, aa in enumerate((cdr2 or "")[:2]): conds.append(f"cdr2[{i}]={aa}")
            for i, aa in enumerate((cdr3 or "")[:2]): conds.append(f"cdr3[{i}]={aa}")
        return conds

    def update_imgt_row(self, imgt_pos_to_aa: Dict[int, str], cdr1: str, cdr2: str, cdr3: str,
                        family: str, hallmarks: str, condition_mode: str):
        self.total_sequences_seen += 1
        if not cdr3:
            return

        conds = self._build_conditions(cdr1, cdr2, cdr3, family, hallmarks, condition_mode)

        contributed = 0
        for pos in self.target_positions:
            aa = imgt_pos_to_aa.get(pos, "-")
            if not aa or aa in ("-", "X", "*"):
                continue
            self.baseline[pos][aa] += 1
            for cond in conds:
                self.counts[pos][cond][aa] += 1
                self.cond_totals[pos][cond] += 1
            contributed += 1

        if contributed > 0:
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

        rules.sort(key=lambda r: (-r["support"], -r["confidence"], -(r["lift"] if r["lift"] != float("inf") else 1e9)))
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

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", required=True)
    ap.add_argument("-o", "--output", default="models/compensation/imgt_csvcols_v6")
    ap.add_argument("--print-every", type=int, default=200000)
    ap.add_argument("--positions", choices=["vernier", "fr3fr4", "all"], default="fr3fr4")
    ap.add_argument("--conditions", choices=["minimal", "full"], default="minimal")
    ap.add_argument("--min-support", type=int, default=5000)
    ap.add_argument("--min-confidence", type=float, default=0.70)
    ap.add_argument("--min-lift", type=float, default=1.15)
    ap.add_argument("--keep-definitional", action="store_true")
    args = ap.parse_args()

    os.makedirs(args.output, exist_ok=True)

    target_positions = choose_positions(args.positions)
    analyzer = StreamingCompensationAnalyzer(target_positions)

    with open(args.csv, "r", newline="") as fh:
        r = csv.reader(fh)
        header_raw = next(r)
        header = [_norm_col(h) for h in header_raw]
        idx = {h: i for i, h in enumerate(header)}

        # IMGT column indices for target positions
        pos_idx = {}
        missing_cols = []
        for p in target_positions:
            col = f"imgt_{p}"
            if col in idx:
                pos_idx[p] = idx[col]
            else:
                missing_cols.append(col)

        if missing_cols:
            raise SystemExit(f"Missing IMGT columns in CSV header (first 20): {missing_cols[:20]}")

        aa_idx = idx.get("aa_v_full", None)
        cdr1_idx = idx.get("cdr1", None)
        cdr2_idx = idx.get("cdr2", None)
        cdr3_idx = idx.get("cdr3", None)
        fam_idx = idx.get("family", None)
        hm_idx = idx.get("hallmarks", None)

        if cdr3_idx is None:
            raise SystemExit("Missing required column: cdr3")

        total = 0
        for row in r:
            total += 1

            aa_v_full = row[aa_idx] if (aa_idx is not None and aa_idx < len(row)) else ""
            cdr1 = row[cdr1_idx] if (cdr1_idx is not None and cdr1_idx < len(row)) else ""
            cdr2 = row[cdr2_idx] if (cdr2_idx is not None and cdr2_idx < len(row)) else ""
            cdr3 = row[cdr3_idx] if (cdr3_idx is not None and cdr3_idx < len(row)) else ""

            # build hallmark+family from IMGT cols (or use provided columns if present)
            imgt_hm = {p: norm_cell(row[pos_idx[p]] if pos_idx[p] < len(row) else "-") for p in HALLMARK_POSITIONS}

            hallmarks = row[hm_idx] if (hm_idx is not None and hm_idx < len(row) and row[hm_idx]) else infer_hallmarks(imgt_hm)
            family = row[fam_idx] if (fam_idx is not None and fam_idx < len(row) and row[fam_idx]) else infer_family(imgt_hm, aa_v_full)

            imgt_map = {p: norm_cell(row[pos_idx[p]] if pos_idx[p] < len(row) else "-") for p in target_positions}

            analyzer.update_imgt_row(
                imgt_pos_to_aa=imgt_map,
                cdr1=cdr1, cdr2=cdr2, cdr3=cdr3,
                family=family, hallmarks=hallmarks,
                condition_mode=args.conditions,
            )

            if total % args.print_every == 0:
                print(f"Rows scanned: {total:,} | seen={analyzer.total_sequences_seen:,} | used={analyzer.total_sequences_used:,}", flush=True)

    print("\n" + "=" * 70)
    print("ANALYSIS COMPLETE")
    print("=" * 70)
    print(f"Seen: {analyzer.total_sequences_seen:,} | Used: {analyzer.total_sequences_used:,}")

    rules = analyzer.extract_rules(
        min_support=args.min_support,
        min_confidence=args.min_confidence,
        min_lift=args.min_lift,
        source="compensation_imgt_csvcols_v6",
    )

    if not args.keep_definitional:
        hallmark_nums = {42, 49, 50, 52}
        def _pos_num(pos_str: str) -> int:
            try: return int(pos_str.replace("IMGT", ""))
            except Exception: return -1
        before = len(rules)
        rules = [r for r in rules
                 if not ((r["condition"].startswith("hallmarks=") or r["condition"].startswith("family="))
                         and _pos_num(r["position"]) in hallmark_nums)]
        removed = before - len(rules)
        if removed:
            print(f"Filtered definitional hallmark rules: removed {removed}, kept {len(rules)}")

    archetypes = analyzer.vernier_archetypes(min_support=max(10000, args.min_support))

    rules_path = os.path.join(args.output, "compensation_imgt_rules.json")
    arch_path  = os.path.join(args.output, "vernier_archetypes_imgt.json")
    with open(rules_path, "w") as f:
        json.dump(rules, f, indent=2)
    with open(arch_path, "w") as f:
        json.dump(archetypes, f, indent=2)

    print(f"Saved rules: {rules_path} ({len(rules)} rules)")
    print(f"Saved archetypes: {arch_path} ({len(archetypes)} conditions)")
    print("\nTop rules:")
    for r in rules[:20]:
        print(f"  {r['position']:7s} -> {r['suggested_aa']} | {r['condition'][:32]:32s} | conf={r['confidence']:.2f} lift={r['lift']} n={r['support']:,}")

if __name__ == "__main__":
    main()
