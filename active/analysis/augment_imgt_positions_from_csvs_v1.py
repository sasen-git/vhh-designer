#!/usr/bin/env python3
"""
Augment existing annotated VHH CSV shards with full IMGT position columns.

- Input: directory of per-shard CSVs (e.g., *.dedup.csv)
- Output: same shards with added imgt_1..imgt_128 columns (integer IMGT positions only)
- Uses ANARCI to compute IMGT numbering from aa_v_full (or fallback sequence column)

Notes:
- We only emit integer IMGT positions 1..128. Insertions like 27A are ignored.
- This is typically what you want for "framework-wide" and stable positional features.
"""

from __future__ import annotations

import argparse
import csv
import glob
import math
import os
import sys
import time
from typing import Dict, List, Optional, Tuple

import pandas as pd

try:
    from anarci import anarci as _anarci
    HAVE_ANARCI = True
except Exception:
    _anarci = None
    HAVE_ANARCI = False


SEQ_COL_CANDIDATES = ["aa_v_full", "v_full", "aa_v", "sequence", "seq"]


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description="Augment shard CSVs with full IMGT integer position columns (imgt_1..imgt_128)."
    )
    ap.add_argument("-i", "--input-dir", required=True, help="Directory containing shard CSVs")
    ap.add_argument("-o", "--output-dir", required=True, help="Directory to write augmented CSVs")
    ap.add_argument("--glob", default="*.dedup.csv", help="Glob pattern inside input-dir (default: *.dedup.csv)")
    ap.add_argument("--seq-col", default="", help="Sequence column to use (default: auto-detect)")
    ap.add_argument("--positions", default="1-128", help="IMGT integer positions, e.g. '1-128' or '1-26,39-55,66-104,118-128'")
    ap.add_argument("--chunk-size", type=int, default=50000, help="Rows per pandas chunk (default: 50000)")
    ap.add_argument("-b", "--batch-size", type=int, default=20000, help="Sequences per ANARCI call (default: 20000)")
    ap.add_argument("--ncpu", type=int, default=9, help="ANARCI/hmmscan CPU threads (default: 9)")
    ap.add_argument("--overwrite-imgt", action="store_true",
                    help="Overwrite existing imgt_* columns instead of only filling missing/blank values")
    ap.add_argument("--progress", choices=["compact", "none"], default="compact",
                    help="Progress display (default: compact)")
    ap.add_argument("--limit-rows", type=int, default=0, help="For testing: limit rows per file (0=all)")
    return ap.parse_args()


def parse_positions(spec: str) -> List[int]:
    """
    Parse '1-128' or '1-26,39-55,66-104,118-128' into sorted unique list.
    """
    spec = spec.strip()
    if not spec:
        raise ValueError("Empty --positions spec")
    parts = [p.strip() for p in spec.split(",") if p.strip()]
    out: List[int] = []
    for p in parts:
        if "-" in p:
            a, b = p.split("-", 1)
            a, b = int(a), int(b)
            if a > b:
                a, b = b, a
            out.extend(range(a, b + 1))
        else:
            out.append(int(p))
    out = sorted(set(out))
    # keep only sensible range
    out = [x for x in out if 1 <= x <= 128]
    return out


def pick_seq_col(cols: List[str], requested: str) -> str:
    cols_l = [c.strip().lower() for c in cols]
    colset = set(cols_l)
    if requested:
        r = requested.strip().lower()
        if r in colset:
            return r
        raise SystemExit(f"--seq-col '{requested}' not found in columns. Available include: {sorted(list(colset))[:50]} ...")
    for c in SEQ_COL_CANDIDATES:
        if c in colset:
            return c
    raise SystemExit(f"Could not auto-detect sequence column. Tried {SEQ_COL_CANDIDATES}. Add --seq-col.")


def chunks_total(nrows: Optional[int], chunk_size: int) -> Optional[int]:
    if nrows is None:
        return None
    return int(math.ceil(nrows / max(1, chunk_size)))


def format_hms(seconds: float) -> str:
    seconds = max(0.0, float(seconds))
    h = int(seconds // 3600)
    m = int((seconds % 3600) // 60)
    s = int(seconds % 60)
    if h > 0:
        return f"{h:d}:{m:02d}:{s:02d}"
    return f"{m:d}:{s:02d}"


def anarci_imgt_numbering(seqs: List[str], ncpu: int) -> List[Optional[Dict[int, str]]]:
    """
    Returns list of dicts mapping integer IMGT position -> residue (single letter).
    If ANARCI fails or sequence invalid, entry is None.
    Insertions (e.g., 27A) are ignored.
    """
    if not HAVE_ANARCI or _anarci is None:
        raise RuntimeError("ANARCI not installed but required. Install 'anarci' in this env.")

    # ANARCI expects strings; empty/too short sequences may fail — we handle by returning None later.
    inp = []
    for i, s in enumerate(seqs):
        s = "" if s is None else str(s).strip()
        inp.append((f"q{i}", s))

    try:
        numbered, _alignment_details, _hit_tables = _anarci(
            inp,
            scheme="imgt",
            assign_germline=False,
            ncpu=ncpu,
        )
    except Exception:
        return [None] * len(seqs)

    out: List[Optional[Dict[int, str]]] = [None] * len(seqs)

    for i, item in enumerate(numbered):
        # numbered[i] is typically: [ (chain_type, numbering_list) ] or [] on failure
        if not item:
            out[i] = None
            continue

        # pick first chain
        try:
            _chain_type, numbering = item[0]
        except Exception:
            out[i] = None
            continue

        # numbering is list of tuples: ((pos, ins), aa) where aa can be '-' for gaps
        d: Dict[int, str] = {}
        try:
            for (pos, ins), aa in numbering:
                # keep only integer positions; ignore insertions like 27A
                if aa is None or aa == "-" or aa == "":
                    continue
                if isinstance(pos, int) and 1 <= pos <= 128:
                    # keep first occurrence if any weird duplicates occur
                    if pos not in d:
                        d[pos] = aa
        except Exception:
            out[i] = None
            continue

        out[i] = d

    return out


def print_progress(prefix: str, done: int, total: int, t0: float, extra: str = ""):
    if total <= 0:
        return
    elapsed = time.time() - t0
    rate = done / elapsed if elapsed > 0 else 0.0
    remaining = (total - done) / rate if rate > 0 else 0.0
    msg = f"{prefix}: {done}/{total} ({done/total*100:5.1f}%) | el {format_hms(elapsed)} | eta {format_hms(remaining)}"
    if extra:
        msg += f" | {extra}"
    # single-line update
    sys.stderr.write("\r" + msg + " " * 10)
    sys.stderr.flush()
    if done >= total:
        sys.stderr.write("\n")
        sys.stderr.flush()


def ensure_dir(path: str):
    os.makedirs(path, exist_ok=True)


def main():
    args = parse_args()
    positions = parse_positions(args.positions)
    ensure_dir(args.output_dir)

    files = sorted(glob.glob(os.path.join(args.input_dir, args.glob)))
    if not files:
        raise SystemExit(f"No files matched: {os.path.join(args.input_dir, args.glob)}")

    eprint("======================================================================")
    eprint("IMGT Augment (CSV shards -> full imgt_1..imgt_128 columns)")
    eprint("======================================================================")
    eprint(f"Input dir:   {args.input_dir}")
    eprint(f"Glob:        {args.glob}")
    eprint(f"Output dir:  {args.output_dir}")
    eprint(f"Batch size:  {args.batch_size}")
    eprint(f"Chunk size:  {args.chunk_size}")
    eprint(f"NCPU:        {args.ncpu}")
    eprint(f"Positions:   {positions[0]}..{positions[-1]}  (n={len(positions)})")
    eprint(f"Overwrite:   {args.overwrite_imgt}")
    eprint(f"Progress:    {args.progress}")
    eprint("")

    for fpath in files:
        base = os.path.basename(fpath)
        outpath = os.path.join(args.output_dir, base.replace(".csv", ".imgtfull.csv"))

        eprint(f"\nProcessing: {base}")
        # Read just header to pick seq col
        header_df = pd.read_csv(fpath, nrows=1, dtype=str, low_memory=False)
        cols = [c.strip().lower() for c in header_df.columns]
        seq_col = pick_seq_col(cols, args.seq_col)

        # normalize columns in streaming reads too
        def norm_cols(df: pd.DataFrame) -> pd.DataFrame:
            df.columns = [c.strip().lower() for c in df.columns]
            return df

        existing_imgt = set([c for c in cols if c.startswith("imgt_")])
        target_imgt_cols = [f"imgt_{p}" for p in positions]
        want_imgt_set = set(target_imgt_cols)

        # If not overwriting, we will only create missing columns OR fill blanks.
        # For schema consistency, we still ensure all target cols exist in output.

        total_rows_written = 0
        total_valid = 0
        t_file0 = time.time()

        reader = pd.read_csv(
            fpath,
            dtype=str,
            low_memory=False,
            chunksize=args.chunk_size,
        )

        first_chunk = True
        chunk_idx = 0

        for chunk in reader:
            chunk_idx += 1
            chunk = norm_cols(chunk)

            if args.limit_rows and total_rows_written >= args.limit_rows:
                break

            if args.limit_rows:
                remaining = args.limit_rows - total_rows_written
                if remaining <= 0:
                    break
                if len(chunk) > remaining:
                    chunk = chunk.iloc[:remaining].copy()

            if seq_col not in chunk.columns:
                raise SystemExit(f"Sequence column '{seq_col}' missing in chunk. Columns: {list(chunk.columns)[:30]} ...")

            seqs = chunk[seq_col].fillna("").astype(str).tolist()

            # Run ANARCI in batches
            n = len(seqs)
            batch_total = int(math.ceil(n / max(1, args.batch_size)))
            imgt_maps: List[Optional[Dict[int, str]]] = [None] * n

            t0 = time.time()
            for bi in range(batch_total):
                bstart = bi * args.batch_size
                bend = min(n, (bi + 1) * args.batch_size)
                batch = seqs[bstart:bend]
                maps = anarci_imgt_numbering(batch, ncpu=args.ncpu)
                imgt_maps[bstart:bend] = maps

                if args.progress == "compact":
                    # show progress within this chunk
                    done = bend
                    extra = f"chunk {chunk_idx} rows {total_rows_written+done:,}"
                    print_progress(base, done, n, t0, extra=extra)

            # Count valid
            valid_mask = [m is not None for m in imgt_maps]
            total_valid += sum(valid_mask)

            # Build new IMGT columns DataFrame
            # (Use "" for missing / invalid)
            newcols = {}
            for p in positions:
                col = f"imgt_{p}"
                # List comp is faster than per-row df.apply
                newcols[col] = [("" if m is None else m.get(p, "")) for m in imgt_maps]
            newdf = pd.DataFrame(newcols)

            # Merge into chunk
            for col in target_imgt_cols:
                if col in chunk.columns:
                    if args.overwrite_imgt:
                        chunk[col] = newdf[col]
                    else:
                        # fill only where existing is blank/NA
                        cur = chunk[col].fillna("").astype(str)
                        chunk[col] = cur.where(cur.str.len() > 0, newdf[col])
                else:
                    chunk[col] = newdf[col]

            # Reorder columns: all non-imgt first, then imgt in numeric order
            non_imgt_cols = [c for c in chunk.columns if not c.startswith("imgt_")]
            # Make sure we don't accidentally drop extra imgt cols outside target set; append them after target if present.
            extra_imgt_cols = sorted([c for c in chunk.columns if c.startswith("imgt_") and c not in want_imgt_set],
                                     key=lambda x: int(x.split("_", 1)[1]) if x.split("_", 1)[1].isdigit() else 10**9)
            out_cols = non_imgt_cols + target_imgt_cols + extra_imgt_cols
            chunk = chunk[out_cols]

            # Write
            mode = "w" if first_chunk else "a"
            chunk.to_csv(
                outpath,
                mode=mode,
                index=False,
                header=first_chunk,
                quoting=csv.QUOTE_MINIMAL,
            )
            first_chunk = False
            total_rows_written += len(chunk)

        eprint(f"Done: wrote {total_rows_written:,} rows to {outpath}")
        eprint(f"  Valid ANARCI mappings (in written rows): {total_valid:,} ({(total_valid/total_rows_written*100 if total_rows_written else 0):.2f}%)")
        eprint(f"  Elapsed: {format_hms(time.time() - t_file0)}")


if __name__ == "__main__":
    main()
