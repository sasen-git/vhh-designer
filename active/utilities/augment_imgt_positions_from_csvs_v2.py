#!/usr/bin/env python3
import os, sys, glob, math, time, argparse, contextlib, io
from typing import List, Dict, Optional, Tuple

import pandas as pd

def eprint(*a, **k):
    print(*a, file=sys.stderr, **k)

def iter_indices(n: int, batch: int):
    for i in range(0, n, batch):
        yield i, min(n, i + batch)

def numeric_imgt_cols(cols):
    out = []
    for c in cols:
        if c.startswith("imgt_"):
            try:
                out.append((int(c.split("_", 1)[1]), c))
            except Exception:
                pass
    return [c for _, c in sorted(out)]

def get_anarci():
    try:
        from anarci import anarci
        return anarci
    except Exception as ex:
        raise RuntimeError("ANARCI not available in this environment") from ex

def anarci_imgt_numbering(
    seqs: List[str],
    ncpu: int = 1,
    quiet: bool = True,
) -> List[Optional[Dict[int, str]]]:
    """
    Returns list aligned to seqs: dict {imgt_pos(int): AA} for 1..128, insertions ignored.
    """
    anarci = get_anarci()

    # ANARCI expects list of (id, seq)
    inp = [(f"s{i}", s) for i, s in enumerate(seqs)]

    kw = dict(scheme="imgt", database="ALL")
    # ANARCI supports ncpu in many builds; safe-guard
    try:
        kw["ncpu"] = int(ncpu)
    except Exception:
        pass

    # Suppress ANARCI stdout/stderr spam
    if quiet:
        sink = io.StringIO()
        with contextlib.redirect_stdout(sink), contextlib.redirect_stderr(sink):
            numbered, alignment_details, hit_tables = anarci(inp, **kw)
    else:
        numbered, alignment_details, hit_tables = anarci(inp, **kw)

    out: List[Optional[Dict[int, str]]] = []
    for res in numbered:
        if res is None:
            out.append(None)
            continue

        # res is typically a list of chains; take the first chain record
        try:
            chain_rec = res[0]
            numbering = chain_rec[0]  # list of (((pos, ins), aa), ...)
        except Exception:
            out.append(None)
            continue

        d: Dict[int, str] = {}
        try:
            for (pos, ins), aa in numbering:
                # ignore gaps/none
                if aa is None or aa == "-" or aa == "":
                    continue
                # ignore insertion letters (e.g. 27A)
                if ins not in ("", " "):
                    continue
                if isinstance(pos, int) and 1 <= pos <= 128:
                    d[pos] = aa
        except Exception:
            out.append(None)
            continue

        out.append(d if d else None)
    return out

def out_name(infile: str) -> str:
    base = os.path.basename(infile)
    if base.endswith(".dedup.csv"):
        stem = base[:-len(".dedup.csv")]
        return f"{stem}.imgtfull.dedup.csv"
    elif base.endswith(".csv"):
        stem = base[:-len(".csv")]
        return f"{stem}.imgtfull.csv"
    else:
        return base + ".imgtfull.csv"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--indir", required=True)
    ap.add_argument("-o", "--outdir", required=True)
    ap.add_argument("--glob", default="*.csv")
    ap.add_argument("--seq-col", default="aa_v_full")
    ap.add_argument("--overwrite-imgt", action="store_true")
    ap.add_argument("--batch-size", type=int, default=20000)
    ap.add_argument("--chunk-size", type=int, default=50000)
    ap.add_argument("--ncpu", type=int, default=1)
    ap.add_argument("--progress", choices=["none", "compact"], default="compact")
    ap.add_argument("--overwrite", action="store_true")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    files = sorted(glob.glob(os.path.join(args.indir, args.glob)))
    if not files:
        raise SystemExit(f"No input files matched {args.glob} in {args.indir}")

    print("="*70)
    print("IMGT Augment (CSV shards -> full imgt_1..imgt_128 columns) [v2]")
    print("="*70)
    print(f"Input dir:   {args.indir}")
    print(f"Glob:        {args.glob}")
    print(f"Output dir:  {args.outdir}")
    print(f"Seq col:     {args.seq_col}")
    print(f"Batch size:  {args.batch_size}")
    print(f"Chunk size:  {args.chunk_size}")
    print(f"NCPU:        {args.ncpu}")
    print(f"Positions:   1..128")
    print(f"Overwrite:   {args.overwrite}")
    print(f"Overwrite IMGT cols: {args.overwrite_imgt}")
    print(f"Progress:    {args.progress}")
    print()

    imgt_cols = [f"imgt_{i}" for i in range(1, 129)]

    for infile in files:
        base = os.path.basename(infile)
        outfile = os.path.join(args.outdir, out_name(infile))

        if os.path.exists(outfile) and not args.overwrite:
            print(f"Skip (exists): {outfile}")
            continue

        print(f"Processing: {base}")
        t0 = time.time()
        total_written = 0
        total_valid = 0

        # stream read
        reader = pd.read_csv(infile, dtype=str, low_memory=False, chunksize=args.chunk_size)
        first = True
        chunk_idx = 0

        for chunk in reader:
            chunk_idx += 1
            if args.seq_col not in chunk.columns:
                raise SystemExit(f"Missing --seq-col {args.seq_col} in {infile}")

            chunk = chunk.copy()
            if args.overwrite_imgt:
                drop_cols = [c for c in chunk.columns if c.startswith("imgt_")]
                if drop_cols:
                    chunk = chunk.drop(columns=drop_cols)

            seqs = chunk[args.seq_col].fillna("").astype(str).tolist()

            # ANARCI in batches
            maps: List[Optional[Dict[int, str]]] = []
            for b0, b1 in iter_indices(len(seqs), args.batch_size):
                maps.extend(anarci_imgt_numbering(seqs[b0:b1], ncpu=args.ncpu, quiet=True))

            # build dense matrix (strings; "" for missing)
            mat = [[""] * 128 for _ in range(len(seqs))]
            valid_here = 0
            for i, d in enumerate(maps):
                if d:
                    valid_here += 1
                    for pos, aa in d.items():
                        mat[i][pos - 1] = aa

            imgt_df = pd.DataFrame(mat, columns=imgt_cols)

            out_chunk = pd.concat([chunk.reset_index(drop=True), imgt_df], axis=1)

            mode = "w" if first else "a"
            out_chunk.to_csv(outfile, index=False, mode=mode, header=first)
            first = False

            total_written += len(out_chunk)
            total_valid += valid_here

            if args.progress == "compact":
                el = time.time() - t0
                rate = total_written / max(el, 1e-9)
                print(f"{base}: wrote {total_written:,} rows | valid {total_valid:,} | {rate:,.1f} rows/s", flush=True)

        el = time.time() - t0
        print(f"Done: wrote {total_written:,} rows -> {outfile}")
        print(f"  Valid ANARCI rows: {total_valid:,}/{total_written:,} ({(100.0*total_valid/max(total_written,1)):.2f}%)")
        print(f"  Elapsed: {el/60:.1f} min\n")

if __name__ == "__main__":
    main()


