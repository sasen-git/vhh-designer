#!/usr/bin/env python3
import sys, csv
from collections import Counter

FR3 = list(range(66, 105))     # 66..104 inclusive (39)
FR4 = list(range(118, 129))    # 118..128 inclusive (11)

def is_missing(x: str) -> bool:
    if x is None:
        return True
    x = x.strip()
    if not x:
        return True
    xl = x.lower()
    return xl in {"na", "nan", "none"}

def is_gap(x: str) -> bool:
    return x.strip() == "-"

def main():
    if len(sys.argv) < 2:
        print("Usage: qc_imgt_coverage_fast.py <csv> [max_rows]", file=sys.stderr)
        sys.exit(2)
    path = sys.argv[1]
    max_rows = int(sys.argv[2]) if len(sys.argv) >= 3 else 0  # 0 = all

    fr3_len = Counter()
    fr4_len = Counter()
    fr3_miss = Counter()
    fr4_miss = Counter()
    n = 0

    want_cols = [f"imgt_{i}" for i in (FR3 + FR4)]

    with open(path, "r", newline="") as f:
        r = csv.reader(f)
        header = next(r)
        col_to_idx = {h: i for i, h in enumerate(header)}

        missing_cols = [c for c in want_cols if c not in col_to_idx]
        if missing_cols:
            print("Missing expected columns (first 30):", missing_cols[:30], file=sys.stderr)
            sys.exit(1)

        fr3_idx = [col_to_idx[f"imgt_{i}"] for i in FR3]
        fr4_idx = [col_to_idx[f"imgt_{i}"] for i in FR4]

        for row in r:
            n += 1

            # FR3
            p3 = 0
            m3 = 0
            for j in fr3_idx:
                x = row[j] if j < len(row) else ""
                if is_missing(x) or is_gap(x):
                    m3 += 1
                else:
                    p3 += 1

            # FR4
            p4 = 0
            m4 = 0
            for j in fr4_idx:
                x = row[j] if j < len(row) else ""
                if is_missing(x) or is_gap(x):
                    m4 += 1
                else:
                    p4 += 1

            fr3_len[p3] += 1
            fr4_len[p4] += 1
            fr3_miss[m3] += 1
            fr4_miss[m4] += 1

            if max_rows and n >= max_rows:
                break

    exp3 = len(FR3)
    exp4 = len(FR4)

    def pct(x): return (x / n * 100.0) if n else 0.0

    print(f"Rows scanned: {n:,}")
    print(f"FR3 present=={exp3}: {fr3_len[exp3]:,} ({pct(fr3_len[exp3]):.2f}%)")
    print(f"FR4 present=={exp4}: {fr4_len[exp4]:,} ({pct(fr4_len[exp4]):.2f}%)")

    print("\nTop FR3 present-lengths:")
    for k,v in fr3_len.most_common(10):
        print(f"  {k:>2}: {v:,} ({pct(v):.2f}%)")

    print("\nTop FR4 present-lengths:")
    for k,v in fr4_len.most_common(10):
        print(f"  {k:>2}: {v:,} ({pct(v):.2f}%)")

    print("\nTop FR3 missing-counts (missing+gaps across 39 positions):")
    for k,v in fr3_miss.most_common(10):
        print(f"  miss={k:>2}: {v:,} ({pct(v):.2f}%)")

    print("\nTop FR4 missing-counts (missing+gaps across 11 positions):")
    for k,v in fr4_miss.most_common(10):
        print(f"  miss={k:>2}: {v:,} ({pct(v):.2f}%)")

if __name__ == "__main__":
    main()
