#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
NPZ FULLSCAN v3 - FIXED with ASCII decoding
- Live dual progress bars (global shards + per-shard mini-bar)
- Direct CDR extraction from numbering arrays (no ANARCI on DB entries)
- Smart FULL-region optimization
- Persistent per-shard stats lines (hits, max ID, top 25% cutoff, skipped, time)
- Per-run CSV summary (timestamped) saved to output dir

Usage example:
python -u npz_fullscan_v3_fixed.py \
  --db-root "/home/user/KA-Search/extracted/oas-paper/oasdb_20230111/Heavy/Human" \
  --query-seq "EIQLQQ...TVSA" \
  --regions full \
  --min-id-full 0.30 \
  --outdir "/home/user/KA-Search/runs_fullscan" \
  --tag "full_30_human"
"""

import os
import re
import time
import math
import glob
import json
import signal
import argparse
import datetime as dt
from typing import Dict, List, Tuple

import numpy as np
from tqdm import tqdm

# Optional: ANARCI will be imported only for query extraction
try:
    from anarci import run_anarci
    _ANARCI_AVAILABLE = True
except Exception:
    run_anarci = None
    _ANARCI_AVAILABLE = False

# ------------------------
# ASCII Decoding (from working npz_fullscan.py)
# ------------------------
AA_PAD = 124  
AA_GAP = 0

def _ints_to_aa(int_row):
    """Convert integer array to amino acid string.
    Ultra-robust: Handles integers, strings, bytes, and any weird characters.
    """
    out = []
    for v in int_row:
        # Handle None or empty
        if v is None:
            continue
            
        # Check if it's already a string/character
        if isinstance(v, (str, bytes)):
            if isinstance(v, bytes):
                try:
                    v = v.decode('utf-8', errors='ignore')
                except:
                    continue
            v = str(v).strip()
            # Only accept single valid amino acid characters
            if len(v) == 1 and v.upper() in 'ACDEFGHIKLMNPQRSTVWY':
                out.append(v.upper())
            continue
        
        # Handle numpy types that look like strings
        try:
            v_str = str(v).strip()
            if len(v_str) == 1 and v_str.upper() in 'ACDEFGHIKLMNPQRSTVWY':
                out.append(v_str.upper())
                continue
        except:
            pass
        
        # Handle integer encoding (ASCII codes)
        try:
            iv = int(v)
            if iv == AA_GAP or iv == AA_PAD:
                continue  # skip gaps/padding
            elif 65 <= iv <= 90:  # 'A'..'Z' ASCII codes
                out.append(chr(iv))
            elif 97 <= iv <= 122:  # 'a'..'z' ASCII codes (lowercase)
                out.append(chr(iv).upper())
        except (ValueError, TypeError):
            # Not a valid integer, skip silently
            continue
    return "".join(out)

def extract_cdr1_from_npz(int_row):
    """Extract CDR1 from NPZ array positions 39-45."""
    return _ints_to_aa(int_row[39:46])

def extract_cdr2_from_npz(int_row):
    """Extract CDR2 from NPZ array positions 66-90 (split by gaps)."""
    part1 = _ints_to_aa(int_row[66:75])
    part2 = _ints_to_aa(int_row[85:91])
    return part1 + part2

def extract_cdr3_from_npz(int_row):
    """Extract CDR3 from NPZ array positions 150-190."""
    cdr3_region = int_row[150:190]
    cdr3_aa = _ints_to_aa(cdr3_region)
    
    # CDR3 starts with C and goes until W
    if cdr3_aa and cdr3_aa[0] == 'C':
        w_pos = cdr3_aa.find('W')
        if w_pos > 0:
            return cdr3_aa[:w_pos]
        else:
            # No W found, return continuous sequence
            parts = []
            for char in cdr3_aa:
                if char:
                    parts.append(char)
                elif len(parts) > 0 and len(parts) < 20:
                    continue
                else:
                    break
            return ''.join(parts)
    return cdr3_aa

def extract_whole_from_npz(int_row):
    """Extract whole variable domain from NPZ array."""
    return _ints_to_aa(int_row[16:190])

def heuristic_cdrh3(full_seq: str) -> str:
    """Backup CDR3 extraction using C...W pattern."""
    s = re.sub(r"[^A-Z]", "", full_seq.upper())
    end = None
    for pat in [r'WGQG', r'W.QG', r'WG.G', r'W..G', r'WG..']:
        m = list(re.finditer(pat, s))
        if m:
            end = m[-1].start()
            break
    if end is None:
        w = s.rfind('W', max(0, len(s) - 30))
        if w == -1:
            return ""
        end = w
    c_start = s.rfind('C', max(0, end - 35), end)
    if c_start == -1:
        c_start = s.rfind('C', 0, end)
        if c_start == -1:
            return ""
    if end - c_start < 3:
        return ""
    return s[c_start:end]

# ------------------------
# Timeout helper (POSIX)
# ------------------------
class _Timeout(Exception):
    pass

def _alarm_handler(signum, frame):
    raise _Timeout()

# ---------------------------------------------
# Identity helpers with PROPER ALIGNMENT
# ---------------------------------------------
def edit_distance(s1: str, s2: str) -> int:
    """Calculate Levenshtein edit distance between two strings."""
    if len(s1) < len(s2):
        return edit_distance(s2, s1)
    
    if len(s2) == 0:
        return len(s1)
    
    previous_row = range(len(s2) + 1)
    for i, c1 in enumerate(s1):
        current_row = [i + 1]
        for j, c2 in enumerate(s2):
            insertions = previous_row[j + 1] + 1
            deletions = current_row[j] + 1
            substitutions = previous_row[j] + (c1 != c2)
            current_row.append(min(insertions, deletions, substitutions))
        previous_row = current_row
    
    return previous_row[-1]


def identity_fraction(a: str, b: str) -> float:
    """
    Calculate sequence identity using edit distance alignment.
    This properly handles insertions/deletions unlike naive 1:1 comparison.
    Returns identity as fraction 0.0-1.0
    """
    if not a or not b:
        return 0.0
    
    max_len = max(len(a), len(b))
    if max_len == 0:
        return 0.0
    
    dist = edit_distance(a, b)
    identity = 1.0 - (dist / max_len)
    
    return max(0.0, identity)

# ---------------------------------------------
# ANARCI-based IMGT CDR extraction (for query only)
# ---------------------------------------------
def extract_cdrs_from_query_sequence(seq: str, scheme: str = "imgt", timeout_s: int = 5) -> Dict[str, str]:
    """Extract IMGT CDRs from query using ANARCI, with timeout and quiet failure.
    Returns dict with keys: cdr1, cdr2, cdr3 (empty strings on failure).
    """
    # quick sanitize
    seq_clean = seq.upper().replace("-", "").replace(".", "")
    if len(seq_clean) < 70 or any(c not in "ACDEFGHIKLMNPQRSTVWY" for c in seq_clean):
        return {"cdr1": "", "cdr2": "", "cdr3": heuristic_cdrh3(seq_clean)}

    if not _ANARCI_AVAILABLE:
        return {"cdr1": "", "cdr2": "", "cdr3": heuristic_cdrh3(seq_clean)}

    # timeout protection
    signal.signal(signal.SIGALRM, _alarm_handler)
    signal.alarm(timeout_s)
    try:
        res = run_anarci([("H", seq_clean)], scheme=scheme, allowed_species=None)
    except _Timeout:
        return {"cdr1": "", "cdr2": "", "cdr3": heuristic_cdrh3(seq_clean)}
    except Exception:
        return {"cdr1": "", "cdr2": "", "cdr3": heuristic_cdrh3(seq_clean)}
    finally:
        signal.alarm(0)

    # Modern ANARCI result indexing: res[1][0][0][0] -> numbering list
    try:
        numbering = res[1][0][0][0]  # list of ((pos, ins), aa)
    except Exception:
        return {"cdr1": "", "cdr2": "", "cdr3": heuristic_cdrh3(seq_clean)}

    if not numbering or not isinstance(numbering, list):
        return {"cdr1": "", "cdr2": "", "cdr3": heuristic_cdrh3(seq_clean)}

    def pos_val(pos_tuple):
        pos = pos_tuple[0]
        ins = pos_tuple[1] if len(pos_tuple) > 1 else ""
        if isinstance(ins, str) and ins.strip():
            return float(f"{pos}.{ord(ins.upper())-64:02d}")
        return float(pos)

    def grab(start: float, end: float) -> str:
        out = []
        for entry in numbering:
            if not isinstance(entry, (tuple, list)) or len(entry) != 2:
                continue
            t, aa = entry
            if not isinstance(t, (tuple, list)) or aa in (None, "-"):
                continue
            v = pos_val(t)
            if start <= v <= end + 0.09:
                out.append(aa)
        return "".join(out)

    # IMGT heavy chain boundaries
    cdr1 = grab(27.0, 38.9)
    cdr2 = grab(56.0, 65.9)
    cdr3 = grab(105.0, 117.9)
    return {"cdr1": cdr1, "cdr2": cdr2, "cdr3": cdr3}

# ---------------------------------------------
# NPZ processing per shard (with mini progress bar)
# ---------------------------------------------
def process_npz_file(npz_path: str, cfg: Dict, query_cdrs: Dict[str, str], shard_idx: int, n_shards: int) -> Dict:
    start = time.time()
    regions = {r.strip().lower() for r in str(cfg["regions"]).split(",") if r.strip()}
    only_full = regions == {"full"}

    # Load numberings from NPZ (FIXED - was looking for 'seqs')
    arr = np.load(npz_path, allow_pickle=True)
    if "numberings" in arr:
        numberings = arr["numberings"]
    else:
        print(f"[ERROR] No 'numberings' key in {os.path.basename(npz_path)}")
        arr.close()
        return {
            "shard": shard_idx,
            "n_total": 0,
            "n_hits": 0,
            "n_skipped": 0,
            "max_id": 0.0,
            "q75": 0.0,
            "time_s": 0.0,
        }
    
    total_seqs = int(len(numberings))

    # per-shard mini bar (live), auto-hides when closed
    inner_bar = tqdm(total=total_seqs, ncols=80, position=1, leave=False)
    inner_bar.set_description(f"Shard {shard_idx:03d}/{n_shards}")

    hits_id_values: List[float] = []
    n_hits = 0
    n_skipped = 0

    # Precompute query parts
    q_full = query_cdrs.get("full", "")
    q_c1 = query_cdrs.get("cdr1", "")
    q_c2 = query_cdrs.get("cdr2", "")
    q_c3 = query_cdrs.get("cdr3", "")

    min_id_full = float(cfg.get("min_id_full", 0.0))
    min_id_c1 = float(cfg.get("min_id_cdr1", 0.0))
    min_id_c2 = float(cfg.get("min_id_cdr2", 0.0))
    min_id_c3 = float(cfg.get("min_id_cdr3", 0.0))

    # FIXED: Now using direct extraction from numbering arrays
    for i in range(total_seqs):
        try:
            numbering_row = numberings[i]
        except Exception:
            n_skipped += 1
            inner_bar.update(1)
            continue

        if only_full:
            # direct full-seq compare using ASCII decoding
            db_full = extract_whole_from_npz(numbering_row)
            
            if len(db_full) < 70:  # Skip very short sequences
                n_skipped += 1
                inner_bar.update(1)
                continue
            
            ident = identity_fraction(q_full, db_full)
            if ident >= min_id_full:
                n_hits += 1
                hits_id_values.append(ident)
        else:
            # extract DB CDRs directly from numbering array
            db_c1 = extract_cdr1_from_npz(numbering_row)
            db_c2 = extract_cdr2_from_npz(numbering_row)
            db_c3 = extract_cdr3_from_npz(numbering_row)
            
            # Fallback for CDR3 if it doesn't start with C
            if not db_c3 or not db_c3.startswith('C'):
                full = extract_whole_from_npz(numbering_row)
                db_c3 = heuristic_cdrh3(full)
            
            if not any([db_c1, db_c2, db_c3]):
                n_skipped += 1
                inner_bar.update(1)
                continue
            
            keep = True
            if "cdr1" in regions:
                if identity_fraction(q_c1, db_c1) < min_id_c1:
                    keep = False
            if keep and "cdr2" in regions:
                if identity_fraction(q_c2, db_c2) < min_id_c2:
                    keep = False
            if keep and "cdr3" in regions:
                ident3 = identity_fraction(q_c3, db_c3)
                if ident3 < min_id_c3:
                    keep = False
                else:
                    hits_id_values.append(ident3)
            if keep:
                n_hits += 1

        inner_bar.update(1)

    inner_bar.close()
    arr.close()
    elapsed = time.time() - start

    # compute stats
    if hits_id_values:
        ids = np.array(hits_id_values, dtype=float)
        max_id = float(np.max(ids))
        q75 = float(np.percentile(ids, 75))
    else:
        max_id = 0.0
        q75 = 0.0

    # Persistent stats line
    pct = (n_hits / total_seqs * 100.0) if total_seqs else 0.0
    print(f"Shard {shard_idx:03d}/{n_shards} — {total_seqs:,} sequences")
    print(
        f"✔ hits={n_hits:,} ({pct:.1f}%)\tmax ID={max_id:.2f}\ttop 25% ≥ {q75:.2f}"
        f"\tskipped={n_skipped:,}\ttime={elapsed:.1f} s"
    )
    print("─" * 79)  # separator

    return {
        "shard": shard_idx,
        "n_total": total_seqs,
        "n_hits": n_hits,
        "n_skipped": n_skipped,
        "max_id": max_id,
        "q75": q75,
        "time_s": elapsed,
    }

# ---------------------------------------------
# MAIN
# ---------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="NPZ FULLSCAN v3 (FIXED) — Direct CDR extraction with progress bars and stats")
    parser.add_argument("--db-root", required=True, help="Directory with .npz shards (e.g., .../Heavy/Human)")
    parser.add_argument("--query-seq", required=True, help="Query AA sequence (FR1..FR4)")
    parser.add_argument("--regions", default="full", help="Comma list: full or any of cdr1,cdr2,cdr3")
    parser.add_argument("--numbering-scheme", dest="scheme", default="imgt")
    parser.add_argument("--min-id-full", type=float, default=0.0)
    parser.add_argument("--min-id-cdr1", type=float, default=0.0)
    parser.add_argument("--min-id-cdr2", type=float, default=0.0)
    parser.add_argument("--min-id-cdr3", type=float, default=0.0)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--tag", default="run")
    args = parser.parse_args()

    print("=" * 80)
    print("NPZ FULLSCAN v3 (FIXED) - Direct CDR Extraction with ASCII Decoding")
    print("=" * 80)
    print(f"Started:       {dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Database:      {args.db_root}")
    print(f"Regions:       {args.regions}")
    print(f"Numbering:     {args.scheme}")
    print(f"CDR Extractor: Direct from numbering arrays (ASCII decoding)")
    print("=" * 80)

    os.makedirs(args.outdir, exist_ok=True)

    # Collect shards
    npz_files = sorted(glob.glob(os.path.join(args.db_root, "*.npz")))
    print(f"Found {len(npz_files)} shards\n")

    if not npz_files:
        print("[ERROR] No NPZ files found!")
        return

    # Extract query CDRs via ANARCI once
    print("Extracting query CDRs (imgt scheme)...\n")
    q_full = re.sub(r"[^A-Za-z]", "", args.query_seq.upper())

    query_cdrs = {"full": q_full, "cdr1": "", "cdr2": "", "cdr3": ""}
    
    # Extract query CDRs if needed
    regions_selected = {r.strip().lower() for r in str(args.regions).split(",") if r.strip()}
    if regions_selected != {"full"}:
        # Need CDRs for comparison
        qc = extract_cdrs_from_query_sequence(q_full, args.scheme, timeout_s=5)
        query_cdrs.update({k: qc.get(k, "") for k in ("cdr1", "cdr2", "cdr3")})

    print("Query sequence (cleaned):")
    print(q_full)
    print("\nQuery CDRs extracted (strict IMGT boundaries):")
    print(f"  CDR-H1: {query_cdrs.get('cdr1','')} (len={len(query_cdrs.get('cdr1',''))})")
    print(f"  CDR-H2: {query_cdrs.get('cdr2','')} (len={len(query_cdrs.get('cdr2',''))})")
    print(f"  CDR-H3: {query_cdrs.get('cdr3','')} (len={len(query_cdrs.get('cdr3',''))})")
    print(f"  Full:   {len(q_full)} aa\n")

    if regions_selected == {"full"}:
        print("[INFO] Full-sequence mode detected — using direct ASCII decoding for database entries.\n")

    cfg = {
        "scheme": args.scheme,
        "regions": args.regions,
        "min_id_full": args.min_id_full,
        "min_id_cdr1": args.min_id_cdr1,
        "min_id_cdr2": args.min_id_cdr2,
        "min_id_cdr3": args.min_id_cdr3,
    }

    # Prepare per-run CSV path (timestamp + tag)
    ts = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    base = re.sub(r"[^A-Za-z0-9_.-]", "_", args.tag or "run")
    csv_path = os.path.join(args.outdir, f"summary_stats_{base}_{ts}.csv")

    # Write CSV header
    with open(csv_path, "w", encoding="utf-8") as fh:
        fh.write("shard_id,total_seqs,hits,skipped,max_id,top25_cutoff,time_s\n")

    total_hits = 0
    outer_bar = tqdm(total=len(npz_files), desc="Processing shards", ncols=100, position=0, leave=True)

    for i, npz_path in enumerate(npz_files, start=1):
        stats = process_npz_file(npz_path, cfg, {**query_cdrs}, i, len(npz_files))
        total_hits += stats["n_hits"]
        # append to CSV
        with open(csv_path, "a", encoding="utf-8") as fh:
            fh.write(
                f"{stats['shard']},{stats['n_total']},{stats['n_hits']},{stats['n_skipped']},{stats['max_id']:.6f},{stats['q75']:.6f},{stats['time_s']:.3f}\n"
            )
        outer_bar.update(1)

    outer_bar.close()

    print(f"\n✅ Done — all {len(npz_files)} shards processed.")
    print(f"Total hits: {total_hits:,}")
    print(f"Summary saved to: {csv_path}")


if __name__ == "__main__":
    main()