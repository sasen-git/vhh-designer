#!/usr/bin/env python3
import csv
import sys
from pathlib import Path

def dedupe_one(in_path: Path, out_path: Path):
    with in_path.open("r", newline="") as fin:
        r = csv.reader(fin)
        header = next(r)

        # keep first occurrence of each column name
        keep_idx = []
        seen = set()
        for i, h in enumerate(header):
            if h not in seen:
                seen.add(h)
                keep_idx.append(i)

        out_header = [header[i] for i in keep_idx]

        with out_path.open("w", newline="") as fout:
            w = csv.writer(fout)
            w.writerow(out_header)
            for row in r:
                # guard short rows
                out_row = [row[i] if i < len(row) else "" for i in keep_idx]
                w.writerow(out_row)

def main():
    if len(sys.argv) < 2:
        print("Usage: dedupe_csv_headers.py <csv1> [csv2 ...]  (writes *.dedup.csv)", file=sys.stderr)
        sys.exit(2)

    for p in sys.argv[1:]:
        in_path = Path(p)
        out_path = in_path.with_suffix(in_path.suffix + ".dedup.csv")
        dedupe_one(in_path, out_path)
        print(f"Wrote {out_path}")

if __name__ == "__main__":
    main()
