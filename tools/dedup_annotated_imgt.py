#!/usr/bin/env python3
import argparse, csv, glob, hashlib, os, sys
from typing import List, Dict

def find_csvs(indir: str) -> List[str]:
    files = sorted(glob.glob(os.path.join(indir, "*.csv")))
    if not files:
        raise SystemExit(f"No .csv files found in: {indir}")
    return files

def first_unique_header_indices(header: List[str]) -> List[int]:
    seen = set()
    keep = []
    for i, h in enumerate(header):
        if h not in seen:
            seen.add(h)
            keep.append(i)
    return keep

def md5_hex(s: str) -> str:
    return hashlib.md5(s.encode("utf-8")).hexdigest()

def main():
    ap = argparse.ArgumentParser(description="Deduplicate annotated IMGT CSV shards by sequence key (streaming bucket method).")
    ap.add_argument("-i", "--input-dir", required=True, help="Directory with shard CSVs (e.g., vhh_full_annotated_v4)")
    ap.add_argument("-o", "--output", required=True, help="Output deduped CSV path")
    ap.add_argument("--key", default="aa_v_full", help="Column name to dedup on (default: aa_v_full)")
    ap.add_argument("--buckets", type=int, default=256, help="Number of hash buckets (default: 256)")
    ap.add_argument("--tmpdir", default=None, help="Temp directory (default: alongside output)")
    args = ap.parse_args()

    csvs = find_csvs(args.input_dir)
    out_path = args.output
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    tmpdir = args.tmpdir or os.path.join(os.path.dirname(out_path), "_dedup_tmp")
    os.makedirs(tmpdir, exist_ok=True)

    # ---------- PASS 1: partition rows into bucket files ----------
    bucket_paths = [os.path.join(tmpdir, f"bucket_{b:03d}.csv") for b in range(args.buckets)]
    bucket_fhs = [None] * args.buckets
    bucket_writers = [None] * args.buckets
    wrote_header = [False] * args.buckets

    total_in = 0
    total_written_to_buckets = 0
    header_out = None
    key_idx_out = None

    for f in csvs:
        with open(f, "r", newline="") as fh:
            r = csv.reader(fh)
            header = next(r)
            keep_idx = first_unique_header_indices(header)
            header_kept = [header[i] for i in keep_idx]

            if header_out is None:
                header_out = header_kept
                if args.key not in header_out:
                    raise SystemExit(f"Key column '{args.key}' not found in header of {f}\nAvailable: {header_out[:40]} ...")
                key_idx_out = header_out.index(args.key)
            else:
                # minimal schema check: ensure key exists; don’t force strict identical headers
                if args.key not in header_kept:
                    raise SystemExit(f"Key column '{args.key}' not found in header of {f}")

            for row in r:
                total_in += 1
                row_kept = [row[i] if i < len(row) else "" for i in keep_idx]
                key = row_kept[key_idx_out] if key_idx_out is not None else ""
                if not key:
                    continue
                h = md5_hex(key)
                b = int(h[:2], 16) % args.buckets  # 0..255 if buckets=256

                if bucket_fhs[b] is None:
                    bucket_fhs[b] = open(bucket_paths[b], "w", newline="")
                    bucket_writers[b] = csv.writer(bucket_fhs[b])
                if not wrote_header[b]:
                    bucket_writers[b].writerow(header_out)
                    wrote_header[b] = True

                bucket_writers[b].writerow(row_kept)
                total_written_to_buckets += 1

        print(f"[pass1] {os.path.basename(f)} done | rows_seen={total_in:,} | bucketed={total_written_to_buckets:,}", file=sys.stderr)

    for fh in bucket_fhs:
        if fh:
            fh.close()

    # ---------- PASS 2: dedup within each bucket and append to final output ----------
    unique = 0
    dupes = 0

    with open(out_path, "w", newline="") as out_fh:
        w = csv.writer(out_fh)
        w.writerow(header_out)

        for b, bp in enumerate(bucket_paths):
            if not os.path.exists(bp) or os.path.getsize(bp) == 0:
                continue

            seen = set()
            with open(bp, "r", newline="") as in_fh:
                r = csv.reader(in_fh)
                header = next(r)  # bucket header
                if args.key not in header:
                    raise SystemExit(f"[bucket {b:03d}] key '{args.key}' missing from bucket header")
                key_idx = header.index(args.key)

                for row in r:
                    k = row[key_idx]
                    if k in seen:
                        dupes += 1
                        continue
                    seen.add(k)
                    w.writerow(row)
                    unique += 1

            print(f"[pass2] bucket {b:03d} done | unique_out={unique:,} | dupes={dupes:,}", file=sys.stderr)

    print("\nDONE", file=sys.stderr)
    print(f"Input rows seen (incl. blanks): {total_in:,}", file=sys.stderr)
    print(f"Rows bucketed (nonblank key):   {total_written_to_buckets:,}", file=sys.stderr)
    print(f"Unique rows written:           {unique:,}", file=sys.stderr)
    print(f"Duplicates removed:            {dupes:,}", file=sys.stderr)
    print(f"Output: {out_path}", file=sys.stderr)

if __name__ == "__main__":
    main()
