#!/usr/bin/env python3
"""
VHH Database ANARCI Annotation Pipeline (streaming, low-memory) - v3
===================================================================

Key differences vs v2:
- Writes ONE CSV per shard incrementally (no giant in-memory concat)
- Processes each NPZ shard in chunks, and ANARCI in batches inside each chunk
- Emits IMGT per-position columns for FR3/FR4 by default (optional FR2)
- Produces a summary JSON + manifest of outputs

Example:
  python annotate_all_shards_v3.py \
    -i data/databases/shards \
    -o data/databases/annotated/vhh_full_annotated_final_v2 \
    --use-anarci \
    -b 2000 \
    --chunk-size 50000
"""

from __future__ import annotations
import argparse
import os
import sys
import glob
import json
import csv
import time
import math
from collections import Counter
from datetime import datetime
from typing import Dict, List, Any, Optional, Iterable, Tuple

import numpy as np

# Optional deps
try:
    import pandas as pd
except Exception:
    pd = None

# ----------------------------- IMGT constants -----------------------------

IMGT_REGIONS = {
    "FR1": (1, 26),
    "CDR1": (27, 38),
    "FR2": (39, 55),
    "CDR2": (56, 65),
    "FR3": (66, 104),
    "CDR3": (105, 117),
    "FR4": (118, 128),
}

HALLMARK_POSITIONS = (42, 49, 50, 52)  # IMGT positions in FR2 for VHH/VH interface motifs


# ----------------------------- ANARCI wrapper -----------------------------

# We support running with or without ANARCI installed.
# If ANARCI is missing, the script can still run in "pass-through" mode
# (i.e., without IMGT re-numbering) as long as --use-anarci is not set.
try:
    from anarci import anarci as _anarci
    HAVE_ANARCI = True
    anarci = _anarci
except Exception:
    HAVE_ANARCI = False
    anarci = None  # type: ignore


# Global knob for ANARCI parallelism (set from CLI in main)
ANARCI_NCPU = 1

def get_anarci():
    try:
        from anarci import anarci
        return anarci
    except Exception as e:
        raise RuntimeError(
            "ANARCI not available. Install with `pip install anarci` (or conda) "
            "or run without --use-anarci."
        ) from e


def _looks_like_numbering_list(x: object) -> bool:
    """Heuristic: is x a list of residue-numbering tuples returned by ANARCI?"""
    if not isinstance(x, list) or not x:
        return False
    t = x[0]
    if not isinstance(t, tuple) or len(t) != 2:
        return False
    a, b = t
    # Common shapes:
    # 1) ((pos, ins), aa)
    # 2) (pos, aa)
    if isinstance(a, tuple) and len(a) == 2 and isinstance(b, str):
        return True
    if isinstance(a, (int, str)) and isinstance(b, str):
        return True
    return False


def _extract_numbering_from_domain(domain: object) -> Optional[list]:
    """Find the first numbering list inside a domain tuple/list in a version-robust way."""
    if domain is None:
        return None

    # Sometimes the domain itself is the numbering list.
    if _looks_like_numbering_list(domain):
        return domain  # type: ignore

    if isinstance(domain, (tuple, list)):
        for part in domain:
            if _looks_like_numbering_list(part):
                return part
            # Some ANARCI builds nest one more level
            if isinstance(part, list) and part and isinstance(part[0], list) and _looks_like_numbering_list(part[0]):
                return part[0]

    return None


def anarci_imgt_numbering(seqs: List[str], quiet: bool = True, ncpu: Optional[int] = None) -> List[Optional[Dict[int, str]]]:
    """Run ANARCI and return {IMGT_position:int -> AA:str} for each input sequence.

    NOTE: ANARCI's return schema differs slightly across versions; this parser is defensive.
    We only keep non-insertion IMGT positions as integers.
    """
    if not HAVE_ANARCI:
        raise RuntimeError("ANARCI not installed but --use-anarci was requested")

    inp = [(f"s{i}", s) for i, s in enumerate(seqs)]

    # ANARCI can be chatty (hmmer fallback messages). Optionally silence stderr + stdout.
    import io
    import contextlib

    try:
        if quiet:
            ctx = contextlib.ExitStack()
            ctx.enter_context(contextlib.redirect_stderr(io.StringIO()))
            ctx.enter_context(contextlib.redirect_stdout(io.StringIO()))
        else:
            ctx = contextlib.nullcontext()

        with ctx:
            # ANARCI is the bottleneck. Try to use multiple CPUs if supported, and skip
            # germline assignment (we only need IMGT numbering). Different ANARCI versions
            # expose different kwargs, so we progressively fall back.
            numbered = alignment_details = hit_tables = None
            kw_tries = [
                {"scheme": "imgt", "allow": set(["H"]), "assign_germline": False, "ncpu": ncpu},
                {"scheme": "imgt", "assign_germline": False, "ncpu": ncpu},
                {"scheme": "imgt", "allow": set(["H"]), "ncpu": ncpu},
                {"scheme": "imgt", "ncpu": ncpu},
                {"scheme": "imgt", "allow": set(["H"])},
                {"scheme": "imgt"},
            ]
            last_err = None
            for kw in kw_tries:
                try:
                    numbered, alignment_details, hit_tables = anarci(inp, **kw)
                    last_err = None
                    break
                except TypeError as e:
                    last_err = e
                    continue
            if numbered is None:
                raise last_err if last_err is not None else RuntimeError("ANARCI call failed")
    except Exception:
        # If ANARCI hard-fails for a batch, return Nones so caller can continue.
        return [None] * len(seqs)

    out: List[Optional[Dict[int, str]]] = []

    # numbered should be list-like per sequence
    if not isinstance(numbered, list):
        return [None] * len(seqs)

    for item in numbered:
        # Extremely defensive: sometimes item can be malformed
        if isinstance(item, (int, float)):
            out.append(None)
            continue
        if not item:
            out.append(None)
            continue

        # item is typically a list of domains; sometimes a single domain tuple.
        domains = item if isinstance(item, list) else [item]

        best_map: Dict[int, str] = {}
        found = False

        for domain in domains:
            numbering_list = _extract_numbering_from_domain(domain)
            if not numbering_list:
                continue

            posmap: Dict[int, str] = {}
            for entry in numbering_list:
                # entry is usually ( (pos,ins), aa ) but may vary.
                try:
                    if isinstance(entry, (tuple, list)) and len(entry) == 2:
                        key, aa = entry
                    elif isinstance(entry, (tuple, list)) and len(entry) == 3:
                        # sometimes (pos, ins, aa)
                        key = (entry[0], entry[1])
                        aa = entry[2]
                    else:
                        continue
                except Exception:
                    continue

                # key can be (pos, ins) OR just pos
                if isinstance(key, tuple) and len(key) == 2:
                    pos, ins = key
                    # Keep only non-insertion positions
                    if ins not in ("", " "):
                        continue
                else:
                    pos = key

                try:
                    pos_i = int(pos)
                except Exception:
                    continue

                if isinstance(aa, str) and aa != "-" and aa != "":
                    posmap[pos_i] = aa

            if posmap:
                # Prefer the first non-empty domain (VHHs should be single-domain)
                best_map = posmap
                found = True
                break

        out.append(best_map if found else None)

    # Ensure length matches inputs
    if len(out) != len(seqs):
        if len(out) < len(seqs):
            out.extend([None] * (len(seqs) - len(out)))
        else:
            out = out[: len(seqs)]

    return out


def safe_str(x: Any) -> str:
    if x is None:
        return ""
    if isinstance(x, (bytes, bytearray)):
        try:
            return x.decode("utf-8")
        except Exception:
            return str(x)
    return str(x)


def slice_region(imgt: Dict[int, str], start: int, end: int) -> str:
    # inclusive range
    return "".join(imgt.get(p, "") for p in range(start, end + 1))


def compute_hallmarks(imgt: Dict[int, str]) -> str:
    return "".join(imgt.get(p, "") for p in HALLMARK_POSITIONS)


def classify_family_from_hallmarks(h: str) -> str:
    """
    Very lightweight family labeling used by downstream scripts.
    You can refine later; this just preserves rough groupings.
    """
    # Common patterns seen in your dataset printouts: FERG, FERF, VGLW, YQRL, FERA, YERL, etc.
    if len(h) != 4 or any(c == "" for c in h):
        return "unknown"
    if h.startswith("Y"):
        # Yxxx are often Y_C2 / Y-related groups
        return "Y_like"
    if h[0] == "F" and h[1:3] == "ER":
        return "F_C2"
    if h.startswith("VG"):
        return "VH_like"
    return "other"


def iter_indices(n: int, chunk_size: int) -> Iterable[Tuple[int, int]]:
    for start in range(0, n, chunk_size):
        yield start, min(n, start + chunk_size)


def _fmt_hms(seconds: float) -> str:
    seconds = int(max(0, seconds))
    h = seconds // 3600
    m = (seconds % 3600) // 60
    s = seconds % 60
    if h:
        return f"{h:d}:{m:02d}:{s:02d}"
    return f"{m:02d}:{s:02d}"


def _compact_progress(label: str, i: int, total: int, start_ts: float, last: float, *, min_interval: float = 1.0) -> float:
    """Single-line progress updater (no tqdm). Returns updated last-print timestamp."""
    now = time.time()
    if (now - last) < min_interval and i < total:
        return last

    done = min(i, total)
    frac = (done / total) if total else 1.0
    elapsed = now - start_ts
    rate = (elapsed / done) if done else 0.0
    eta = rate * (total - done) if done else 0.0
    msg = f"{label}: {done}/{total} ({frac*100:4.1f}%) | el {_fmt_hms(elapsed)} | eta {_fmt_hms(eta)}"

    # Keep the line short to avoid wrapping in narrow terminals.
    if len(msg) > 110:
        msg = msg[:107] + "..."

    if sys.stderr.isatty():
        # Carriage-return update in-place.
        pad = max(0, 110 - len(msg))
        sys.stderr.write("\r" + msg + (" " * pad))
        sys.stderr.flush()
    else:
        print(msg, file=sys.stderr)
    return now


# ----------------------------- core processing -----------------------------

def process_npz_shard(
    shard_path: str,
    out_dir: str,
    use_anarci: bool,
    anarci_batch_size: int,
    anarci_ncpu: int,
    chunk_size: int,
    emit_imgt_fr3fr4: bool,
    emit_imgt_fr2: bool,
    progress: str = "compact",
    overwrite: bool = False,
) -> Dict[str, Any]:
    base = os.path.splitext(os.path.basename(shard_path))[0]
    out_csv = os.path.join(out_dir, f"{base}.csv")
    if overwrite and os.path.exists(out_csv):
        os.remove(out_csv)


    npz = np.load(shard_path, allow_pickle=True)
    keys = list(npz.files)

    # Identify sequence field
    seq_key = None
    for cand in ["aa_v_full", "seq", "sequence", "aa", "aa_full"]:
        if cand in keys:
            seq_key = cand
            break
    if seq_key is None:
        raise RuntimeError(f"No sequence column found in {shard_path}. Keys: {keys[:20]}...")

    seqs = npz[seq_key]
    n = len(seqs)

    # Progress for ANARCI (compact single-line updates; avoids tqdm spam).
    is_tty = sys.stderr.isatty()
    show_progress = (progress != "none")
    total_batches = math.ceil(n / max(1, anarci_batch_size)) if use_anarci else 0
    done_batches = 0
    prog_start = time.time()
    last_print = 0.0

    def _progress_update(force: bool = False) -> None:
        """Single-line progress update for ANARCI.

        Uses '\r' when running in a real TTY; otherwise prints occasional lines.
        """
        nonlocal last_print
        if not use_anarci or not show_progress:
            return
        now = time.time()
        if (not force) and (now - last_print < 1.0):
            return
        if total_batches <= 0:
            return
        pct = (done_batches / total_batches) * 100.0
        elapsed = now - prog_start
        per_batch = (elapsed / done_batches) if done_batches else 0.0
        eta = per_batch * max(0, (total_batches - done_batches))

        base = os.path.basename(shard_path)
        msg = f"ANARCI {base}: {done_batches}/{total_batches} ({pct:4.1f}%) | el {_fmt_hms(elapsed)} | eta {_fmt_hms(eta)}"
        # Keep the line from wrapping in typical terminals.
        if len(msg) > 110:
            msg = msg[:107] + "..."

        if is_tty:
            print("\r" + msg.ljust(110), end="", file=sys.stderr, flush=True)
        else:
            # In non-TTY contexts, don't spam.
            print(msg, file=sys.stderr, flush=True)
        last_print = now

    # Metadata keys to carry through (everything except sequences + common redundant fields)
    skip = {seq_key}
    meta_keys = [k for k in keys if k not in skip]

    # Output columns (stable order)
    core_cols = [
        "aa_v_full",
        "valid",
        "error",
        "hallmarks",
        "family",
        "fr1","cdr1","fr2","cdr2","fr3","cdr3","fr4",
        "len_v","len_fr1","len_cdr1","len_fr2","len_cdr2","len_fr3","len_cdr3","len_fr4",
        "source_shard",
    ]
    imgt_cols: List[str] = []
    if emit_imgt_fr3fr4:
        imgt_cols += [f"imgt_{p}" for p in range(66, 105)]
        imgt_cols += [f"imgt_{p}" for p in range(118, 129)]
    if emit_imgt_fr2:
        imgt_cols += [f"imgt_{p}" for p in range(39, 56)]
    # Always include hallmark columns explicitly for debugging/joins
    imgt_cols += [f"imgt_{p}" for p in HALLMARK_POSITIONS if f"imgt_{p}" not in imgt_cols]

    output_cols = core_cols + imgt_cols + meta_keys

    # Write header now so you can watch files appear immediately
    os.makedirs(out_dir, exist_ok=True)
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=output_cols, extrasaction="ignore")
        w.writeheader()

    valid_count = 0
    fam_counts = Counter()
    start_time = datetime.now()

    for cstart, cend in iter_indices(n, chunk_size):
        # Prepare chunk sequences as python strings
        chunk_seqs = [safe_str(s) for s in seqs[cstart:cend]]

        # Run ANARCI in sub-batches if requested
        imgt_dicts: List[Optional[Dict[int, str]]]
        if use_anarci:
            imgt_dicts = []
            for bstart, bend in iter_indices(len(chunk_seqs), anarci_batch_size):
                imgt_dicts.extend(anarci_imgt_numbering(chunk_seqs[bstart:bend], ncpu=anarci_ncpu))
                done_batches += 1
                if show_progress:
                    _progress_update()
        else:
            # Without ANARCI we cannot produce IMGT columns robustly
            imgt_dicts = [None] * len(chunk_seqs)

        # Gather metadata arrays for chunk
        meta_chunk = {k: npz[k][cstart:cend] for k in meta_keys}

        rows: List[Dict[str, Any]] = []
        for i, seq in enumerate(chunk_seqs):
            row: Dict[str, Any] = {"aa_v_full": seq, "source_shard": base}
            # attach metadata
            for k in meta_keys:
                v = meta_chunk[k][i]
                row[k] = safe_str(v)

            imgt = imgt_dicts[i]
            if imgt is None:
                row["valid"] = False
                row["error"] = "anarci_failed" if use_anarci else "no_anarci"
                row["hallmarks"] = ""
                row["family"] = "unknown"
                # leave regions empty
                for c in ["fr1","cdr1","fr2","cdr2","fr3","cdr3","fr4"]:
                    row[c] = ""
                for c in ["len_v","len_fr1","len_cdr1","len_fr2","len_cdr2","len_fr3","len_cdr3","len_fr4"]:
                    row[c] = ""
                rows.append(row)
                continue

            # regions
            fr1 = slice_region(imgt, *IMGT_REGIONS["FR1"])
            cdr1 = slice_region(imgt, *IMGT_REGIONS["CDR1"])
            fr2 = slice_region(imgt, *IMGT_REGIONS["FR2"])
            cdr2 = slice_region(imgt, *IMGT_REGIONS["CDR2"])
            fr3 = slice_region(imgt, *IMGT_REGIONS["FR3"])
            cdr3 = slice_region(imgt, *IMGT_REGIONS["CDR3"])
            fr4 = slice_region(imgt, *IMGT_REGIONS["FR4"])

            row["fr1"],row["cdr1"],row["fr2"],row["cdr2"],row["fr3"],row["cdr3"],row["fr4"] = fr1,cdr1,fr2,cdr2,fr3,cdr3,fr4
            row["len_v"] = str(len(seq))
            row["len_fr1"] = str(len(fr1))
            row["len_cdr1"] = str(len(cdr1))
            row["len_fr2"] = str(len(fr2))
            row["len_cdr2"] = str(len(cdr2))
            row["len_fr3"] = str(len(fr3))
            row["len_cdr3"] = str(len(cdr3))
            row["len_fr4"] = str(len(fr4))

            h = compute_hallmarks(imgt)
            row["hallmarks"] = h
            fam = classify_family_from_hallmarks(h)
            row["family"] = fam

            # per-position IMGT columns
            for p in HALLMARK_POSITIONS:
                row[f"imgt_{p}"] = imgt.get(p, "")
            if emit_imgt_fr3fr4:
                for p in range(66, 105):
                    row[f"imgt_{p}"] = imgt.get(p, "")
                for p in range(118, 129):
                    row[f"imgt_{p}"] = imgt.get(p, "")
            if emit_imgt_fr2:
                for p in range(39, 56):
                    row[f"imgt_{p}"] = imgt.get(p, "")

            row["valid"] = True
            row["error"] = ""
            rows.append(row)

            valid_count += 1
            fam_counts[fam] += 1

        # append chunk rows
        with open(out_csv, "a", newline="") as f:
            w = csv.DictWriter(f, fieldnames=output_cols, extrasaction="ignore")
            w.writerows(rows)

        if (cend % (chunk_size * 5) == 0) or (cend == n):
            elapsed = (datetime.now() - start_time).total_seconds()
            rate = (cend / elapsed) if elapsed > 0 else 0
            print(f"    {base}: wrote {cend:,}/{n:,} rows | valid so far {valid_count:,} | {rate:,.0f} seq/s")

    # Finish progress line cleanly.
    if use_anarci and show_progress and is_tty:
        _progress_update(force=True)
        print("", file=sys.stderr)

    shard_elapsed = (datetime.now() - start_time).total_seconds()
    return {
        "shard": base,
        "input": n,
        "valid": valid_count,
        "families": dict(fam_counts),
        "out_csv": out_csv,
        "elapsed_seconds": shard_elapsed,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input-dir", required=True, help="Directory containing NPZ shards")
    ap.add_argument("-o", "--output-dir", required=True, help="Output directory (CSV per shard + summary)")
    ap.add_argument("--use-anarci", action="store_true", help="Run ANARCI to derive IMGT numbering")
    ap.add_argument("-b", "--batch-size", type=int, default=2000, help="ANARCI batch size (default 2000)")
    ap.add_argument("--ncpu", type=int, default=max(1, (os.cpu_count() or 2) - 1), help="ANARCI parallel workers (if supported)")
    ap.add_argument("--overwrite", action="store_true", help="Overwrite existing per-shard CSVs instead of appending")
    ap.add_argument("--chunk-size", type=int, default=50000, help="Rows per streaming chunk (default 50000)")
    ap.add_argument(
        "--progress",
        choices=["compact", "none"],
        default="compact",
        help="Progress display for ANARCI: 'compact' (single line) or 'none' (default compact)",
    )
    ap.add_argument("--no-imgt-fr3fr4", action="store_true", help="Do NOT emit per-position IMGT columns for FR3/FR4")
    ap.add_argument("--emit-imgt-fr2", action="store_true", help="Also emit per-position IMGT columns for FR2 (39-55)")
    args = ap.parse_args()
    global ANARCI_NCPU
    ANARCI_NCPU = max(1, int(args.ncpu))

    input_dir = args.input_dir
    out_dir = args.output_dir
    use_anarci = args.use_anarci

    npz_files = sorted(glob.glob(os.path.join(input_dir, "*.npz")))
    npz_files = [f for f in npz_files if not f.endswith("Zone.Identifier")]

    print("=" * 70)
    print("VHH Database ANARCI Annotation Pipeline (streaming v3)")
    print("=" * 70)
    print(f"Input directory: {input_dir}")
    print(f"Output directory: {out_dir}")
    print(f"Use ANARCI: {use_anarci}")
    print(f"Batch size: {args.batch_size}")
    print(f"Chunk size: {args.chunk_size}")
    print(f"Emit IMGT FR3/FR4 per-position: {not args.no_imgt_fr3fr4}")
    print(f"Emit IMGT FR2 per-position: {args.emit_imgt_fr2}")
    if use_anarci:
        print(f"Progress: {args.progress}")
    print()
    print(f"Found {len(npz_files)} NPZ files:")
    for f in npz_files:
        size_mb = os.path.getsize(f) / 1024 / 1024
        print(f"  {os.path.basename(f)}: {size_mb:.1f} MB")
    print()

    os.makedirs(out_dir, exist_ok=True)

    all_valid = 0
    all_input = 0
    all_fams = Counter()
    manifest = []
    global_start = datetime.now()

    for shard_path in npz_files:
        print(f"\n  Processing: {os.path.basename(shard_path)}")
        info = process_npz_shard(
            shard_path=shard_path,
            out_dir=out_dir,
            use_anarci=use_anarci,
            anarci_batch_size=args.batch_size,
            anarci_ncpu=args.ncpu,
            chunk_size=args.chunk_size,
            emit_imgt_fr3fr4=not args.no_imgt_fr3fr4,
            emit_imgt_fr2=args.emit_imgt_fr2,
            progress=args.progress,
            overwrite=args.overwrite,
        )
        manifest.append(info)
        all_input += info["input"]
        all_valid += info["valid"]
        all_fams.update(info["families"])
        print(f"    Valid: {info['valid']:,}/{info['input']:,}  | wrote {info['out_csv']}")

    elapsed = (datetime.now() - global_start).total_seconds()
    summary = {
        "input_dir": input_dir,
        "output_dir": out_dir,
        "use_anarci": use_anarci,
        "batch_size": args.batch_size,
        "chunk_size": args.chunk_size,
        "emit_imgt_fr3fr4": (not args.no_imgt_fr3fr4),
        "emit_imgt_fr2": args.emit_imgt_fr2,
        "total_input": all_input,
        "total_valid": all_valid,
        "family_counts": dict(all_fams),
        "shards_processed": len(manifest),
        "elapsed_seconds": elapsed,
        "timestamp": datetime.now().isoformat(),
        "manifest": manifest,
    }
    summary_path = os.path.join(out_dir, "annotation_summary.json")
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"\nSummary saved: {summary_path}")
    print(f"Done. Total valid: {all_valid:,}/{all_input:,}  | elapsed {elapsed/3600:.2f} h")


if __name__ == "__main__":
    main()
