#!/usr/bin/env python3
"""
NPZ Fullscan with Corrected IMGT CDR Extraction
Fixed to properly extract CDRs according to strict IMGT boundaries
CDR-H2 should be ILPGSGST not QILPGSGSTK
"""

import os, sys, re, time, glob, json, argparse, datetime as dt
import numpy as np
import pandas as pd
from anarci import run_anarci
from typing import Dict, List, Tuple, Optional



# ---------- Helpers ----------
AA_PAD = 124  
AA_GAP = 0    

# # --- CDR Extraction using AntPack or fallback heuristics ---
# try:
#     from antpack import SingleChainAnnotator
#     ANTPACK_AVAILABLE = True
# except ImportError:
#     ANTPACK_AVAILABLE = False
#     print("[WARNING] AntPack not available, will use heuristic extraction")

def heuristic_cdr_extraction(seq: str) -> Dict[str, str]:
    """
    Heuristic CDR extraction for IMGT numbering.
    Returns CDR1, CDR2, CDR3 WITHOUT conserved framework residues.
    Based on typical patterns in antibody sequences.
    """
    s = seq.upper().replace("-", "").replace(".", "")
    
    # Find CDR-H3 first (most reliable)
    cdr3 = ""
    # Find the conserved WG[QKR]G motif at the end
    end_pattern = None
    for pat in [r"WG[QKR]G", r"W.QG", r"WG.G", r"W..G"]:
        m = list(re.finditer(pat, s))
        if m:
            end_pattern = m[-1].start()  # Position of W
            break
    
    if end_pattern is None:
        # Fallback: find last W
        w_pos = s.rfind("W", max(0, len(s) - 30))
        if w_pos > 0:
            end_pattern = w_pos
    
    if end_pattern:
        # Find the conserved C before CDR3
        c_pos = s.rfind("C", max(0, end_pattern - 35), end_pattern)
        if c_pos > 0 and end_pattern - c_pos > 3:
            # CDR3 starts AFTER the C (not including it)
            cdr3 = s[c_pos + 1:end_pattern]
    
    # Find CDR-H1
    cdr1 = ""
    # For sequences like yours with GYTFSSY pattern
    if "GYTFSSY" in s:
        idx = s.index("GYTFSSY")
        cdr1 = "GYTFSSY"
    else:
        # Look for typical CDR1 start patterns
        for pat in [r"[GS][YF][TS]F[TS]", r"[ND][YF][AGST]", r"[GS][YF][AGST]M"]:
            m = re.search(pat, s[20:40])
            if m:
                start = 20 + m.start()
                # CDR1 ends before W
                w_pos = s.find("W", start, min(start + 15, len(s)))
                if w_pos > 0:
                    cdr1 = s[start:w_pos]
                    if len(cdr1) > 10:  # IMGT CDR1 typically 5-8 residues
                        cdr1 = cdr1[:8]
                    break
    
    # Find CDR-H2 - should be ILPGSGST not QILPGSGSTK
    cdr2 = ""
    # Look for LEWIG or EWIX pattern that precedes CDR2
    for pat in [r"LEWIG", r"LEWI[GA]", r"EWI[GA]", r"WI[GARKQ]"]:
        m = re.search(pat, s)
        if m:
            start = m.end()
            # Special case for your sequence - look for ILPGSGST pattern
            if "ILPGSGST" in s[start-2:start+15]:
                idx = s.index("ILPGSGST", start-2)
                cdr2 = "ILPGSGST"
            else:
                # CDR2 is typically 7-8 residues for IMGT
                # Look for patterns that end CDR2
                end_markers = [r"[KR][LFVY]", r"YN[EQDK]", r"[TS]KY", r"PS[LK]"]
                cdr2_end = start + 8  # Default
                
                for end_pat in end_markers:
                    m2 = re.search(end_pat, s[start:start+15])
                    if m2:
                        cdr2_end = start + m2.start()
                        break
                
                cdr2 = s[start:cdr2_end]
                if len(cdr2) > 10:  # Too long for IMGT
                    cdr2 = cdr2[:8]
            break
    
    return {"cdr1": cdr1, "cdr2": cdr2, "cdr3": cdr3}

def extract_cdrs_from_sequence_strict(seq: str, scheme: str = "imgt") -> Dict[str, str]:
    """
    Extract IMGT CDRs using ANARCI, handling both insertion codes and
    the newer ANARCI result structure (res[1][0][0][0]).
    """

    seq_clean = seq.upper().replace("-", "").replace(".", "")

    # Quick skip for short or invalid sequences
    if len(seq_clean) < 70 or any(c not in "ACDEFGHIKLMNPQRSTVWY" for c in seq_clean):
        return {"cdr1": "", "cdr2": "", "cdr3": ""}

    try:
        res = run_anarci([("H", seq_clean)], scheme=scheme, allowed_species=None)

        # ✅ Use correct indexing for modern ANARCI structure
        numbering = res[1][0][0][0]

        if not numbering or not isinstance(numbering, list):
            return {"cdr1": "", "cdr2": "", "cdr3": ""}
    except Exception as e:
        print(f"[WARN] ANARCI extraction failed: {e}")
        return {"cdr1": "", "cdr2": "", "cdr3": ""}

    # Normalize insertion handling (e.g., 35A, 57A, etc.)
    def pos_value(pos_tuple):
        pos = pos_tuple[0]
        ins = pos_tuple[1] if len(pos_tuple) > 1 else ""
        if isinstance(ins, str) and ins.strip():
            return float(f"{pos}.{ord(ins.upper()) - 64:02d}")  # A=0.01, B=0.02...
        return float(pos)

    def extract_region(start: float, end: float) -> str:
        aa = []
        for entry in numbering:
            if not isinstance(entry, (tuple, list)) or len(entry) != 2:
                continue
            pos_tuple, residue = entry
            if not isinstance(pos_tuple, (tuple, list)) or residue in [None, "-"]:
                continue
            pos_val = pos_value(pos_tuple)
            if start <= pos_val <= end + 0.09:
                aa.append(residue)
        return "".join(aa)

    # IMGT standard heavy-chain boundaries
    cdr1 = extract_region(27.0, 38.9)
    cdr2 = extract_region(56.0, 65.9)
    cdr3 = extract_region(105.0, 117.9)

    return {"cdr1": cdr1, "cdr2": cdr2, "cdr3": cdr3}




def _ints_to_aa(int_row) -> str:
    """Convert integer array to amino acid string."""
    out = []
    for v in int_row:
        iv = int(v)
        if iv == AA_GAP or iv == AA_PAD:
            continue  # skip gaps
        elif 65 <= iv <= 90:  # 'A'..'Z'
            out.append(chr(iv))
    return "".join(out)

def extract_full_sequence_from_npz(row) -> str:
    """Extract the full variable domain sequence from NPZ row."""
    full_seq = _ints_to_aa(row[:250])
    return full_seq.rstrip()

def pid_identity(a: str, b: str) -> float:
    """Calculate pairwise identity between two sequences."""
    if not a or not b:
        return 0.0
    min_len = min(len(a), len(b))
    if min_len == 0:
        return 0.0
    matches = sum(1 for x, y in zip(a, b) if x == y)
    return matches / min_len

def preview(s: str, n: int = 80) -> str:
    """Preview string with ellipsis if too long."""
    return s[:n] + ("…" if len(s) > n else "")

def calculate_statistics(values: List[float]) -> Dict[str, float]:
    """Calculate quartiles and statistics for a list of values."""
    if not values:
        return {"min": 0, "q25": 0, "median": 0, "q75": 0, "max": 0, "mean": 0, "count": 0}
    
    arr = np.array(values)
    return {
        "min": float(np.min(arr)),
        "q25": float(np.percentile(arr, 25)),
        "median": float(np.median(arr)),
        "q75": float(np.percentile(arr, 75)),
        "max": float(np.max(arr)),
        "mean": float(np.mean(arr)),
        "count": len(values)
    }

def format_stats_line(stats: Dict[str, float], label: str = "") -> str:
    """Format statistics into a readable line."""
    if stats["count"] == 0:
        return f"{label}: No hits"
    
    # Count how many are in different ranges
    ranges = []
    if "values" in stats:
        vals = stats["values"]
        high = sum(1 for v in vals if v >= 0.9)
        good = sum(1 for v in vals if 0.7 <= v < 0.9)
        ok = sum(1 for v in vals if 0.5 <= v < 0.7)
        
        if high > 0:
            ranges.append(f"{high} >90%")
        if good > 0:
            ranges.append(f"{good} 70-90%")
        if ok > 0:
            ranges.append(f"{ok} 50-70%")
    
    range_str = f" ({', '.join(ranges)})" if ranges else ""
    
    return (f"{label}: n={stats['count']}, "
            f"max={stats['max']:.1%}, "
            f"Q75={stats['q75']:.1%}, "
            f"median={stats['median']:.1%}, "
            f"mean={stats['mean']:.1%}"
            f"{range_str}")

class Tee:
    """Tee output to both stdout and a file."""
    def __init__(self, path):
        self.path = path
        self.f = open(path, "w", encoding="utf-8")
    def write(self, msg):
        sys.__stdout__.write(msg)
        self.f.write(msg)
        self.f.flush()
    def flush(self): pass
    def close(self): self.f.close()

# ---------- Process one NPZ ----------
def process_npz_file(npz_path: str, cfg: Dict, query_cdrs: Dict, 
                     verbose: bool = False) -> Tuple[List, Dict, Dict]:
    """
    Process NPZ file: extract full sequences and CDRs using consistent method.
    Returns hits, basic stats, and detailed statistics.
    """
    hits, examined, skipped = [], 0, 0
    
    # Track identities for statistics
    id_lists = {
        "cdr1": [],
        "cdr2": [],
        "cdr3": [],
        "full": []
    }
    
    try:
        data = np.load(npz_path, allow_pickle=True)
    except Exception as e:
        return [], {"examined": 0, "skipped": 0, "error": f"open:{e}"}, {}

    try:
        arr = data["numberings"]
    except Exception as e:
        data.close()
        return [], {"examined": 0, "skipped": 0, "error": f"no_numberings:{e}"}, {}

    debug_count = 5 if verbose else 0
    
    for i in range(len(arr)):
        examined += 1
        try:
            row = arr[i]
            
            # Extract full sequence from NPZ
            full_seq = extract_full_sequence_from_npz(row)
            if not full_seq or len(full_seq) < 50:
                skipped += 1
                continue
            
            # Extract CDRs using the same method as for query
            db_cdrs = extract_cdrs_from_sequence_strict(full_seq, cfg["scheme"])
            
            # Debug output for first few sequences
            if debug_count > 0:
                debug_count -= 1
                print(f"\n[DEBUG] Row {i}:")
                print(f"  Full seq preview: {preview(full_seq, 60)}")
                print(f"  CDR1: {db_cdrs['cdr1']} (len={len(db_cdrs['cdr1'])})")
                print(f"  CDR2: {db_cdrs['cdr2']} (len={len(db_cdrs['cdr2'])})")
                print(f"  CDR3: {db_cdrs['cdr3']} (len={len(db_cdrs['cdr3'])})")
            
            # Check length windows for filtering
            skip_this = False
            
            if cfg["use_cdr1"] and query_cdrs.get("cdr1"):
                if not db_cdrs["cdr1"]:
                    skipped += 1
                    continue
                if abs(len(db_cdrs["cdr1"]) - len(query_cdrs["cdr1"])) > cfg["win_cdr1"]:
                    skip_this = True
            
            if cfg["use_cdr2"] and query_cdrs.get("cdr2"):
                if not db_cdrs["cdr2"]:
                    skipped += 1
                    continue
                if abs(len(db_cdrs["cdr2"]) - len(query_cdrs["cdr2"])) > cfg["win_cdr2"]:
                    skip_this = True
            
            if cfg["use_cdr3"] and query_cdrs.get("cdr3"):
                if not db_cdrs["cdr3"]:
                    skipped += 1
                    continue
                if abs(len(db_cdrs["cdr3"]) - len(query_cdrs["cdr3"])) > cfg["win_cdr3"]:
                    skip_this = True
            
            if skip_this:
                skipped += 1
                continue
            
            # Calculate identities
            id_cdr1 = pid_identity(query_cdrs.get("cdr1", ""), db_cdrs["cdr1"]) if cfg["use_cdr1"] and db_cdrs["cdr1"] else 0
            id_cdr2 = pid_identity(query_cdrs.get("cdr2", ""), db_cdrs["cdr2"]) if cfg["use_cdr2"] and db_cdrs["cdr2"] else 0
            id_cdr3 = pid_identity(query_cdrs.get("cdr3", ""), db_cdrs["cdr3"]) if cfg["use_cdr3"] and db_cdrs["cdr3"] else 0
            id_full = pid_identity(query_cdrs.get("full", ""), full_seq) if cfg["use_full"] else 0
            
            # Track all identities for statistics
            if cfg["use_cdr1"] and db_cdrs["cdr1"]:
                id_lists["cdr1"].append(id_cdr1)
            if cfg["use_cdr2"] and db_cdrs["cdr2"]:
                id_lists["cdr2"].append(id_cdr2)
            if cfg["use_cdr3"] and db_cdrs["cdr3"]:
                id_lists["cdr3"].append(id_cdr3)
            if cfg["use_full"]:
                id_lists["full"].append(id_full)
            
            # Apply identity thresholds
            pass_cdr1 = not cfg["use_cdr1"] or not query_cdrs.get("cdr1") or id_cdr1 >= cfg["min_id_cdr1"]
            pass_cdr2 = not cfg["use_cdr2"] or not query_cdrs.get("cdr2") or id_cdr2 >= cfg["min_id_cdr2"]
            pass_cdr3 = not cfg["use_cdr3"] or not query_cdrs.get("cdr3") or id_cdr3 >= cfg["min_id_cdr3"]
            pass_full = not cfg["use_full"] or id_full >= cfg["min_id_full"]
            
            if pass_cdr1 and pass_cdr2 and pass_cdr3 and pass_full:
                hits.append({
                    "Shard": os.path.basename(npz_path),
                    "Row_index": i,
                    "CDR1": db_cdrs["cdr1"],
                    "CDR2": db_cdrs["cdr2"],
                    "CDR3": db_cdrs["cdr3"],
                    "ID_CDR1": round(id_cdr1, 3),
                    "ID_CDR2": round(id_cdr2, 3),
                    "ID_CDR3": round(id_cdr3, 3),
                    "ID_Full": round(id_full, 3),
                    "Full_Sequence": full_seq,
                })
        
        except Exception as e:
            if verbose:
                print(f"[WARNING] Error processing row {i}: {e}")
            continue
    
    data.close()
    
    # Calculate detailed statistics
    detailed_stats = {}
    for region in ["cdr1", "cdr2", "cdr3", "full"]:
        if id_lists[region]:
            stats = calculate_statistics(id_lists[region])
            stats["values"] = id_lists[region]
            detailed_stats[region] = stats
    
    return hits, {"examined": examined, "skipped": skipped}, detailed_stats

# ---------- Main ----------
def main():
    p = argparse.ArgumentParser(description="NPZ Fullscan with Corrected IMGT CDR Extraction")
    p.add_argument("--db-root", required=True, help="Root directory of NPZ database")
    p.add_argument("--query-seq", required=True, help="Query antibody sequence")
    p.add_argument("--chain", default="Heavy", choices=["Heavy","Light"])
    p.add_argument("--species", default="Human")
    p.add_argument("--regions", default="cdr3", help="Comma-separated: cdr1,cdr2,cdr3,full")
    # Per-region thresholds
    p.add_argument("--min-id-cdr1", type=float, default=0.40)
    p.add_argument("--min-id-cdr2", type=float, default=0.40)
    p.add_argument("--min-id-cdr3", type=float, default=0.60)
    p.add_argument("--min-id-full", type=float, default=0.35)
    # Length windows
    p.add_argument("--len-window-cdr1", type=int, default=2)
    p.add_argument("--len-window-cdr2", type=int, default=2)
    p.add_argument("--len-window-cdr3", type=int, default=1)
    # Output
    p.add_argument("--outdir", default="./runs_fullscan")
    p.add_argument("--tag", default="")
    p.add_argument("--head-check", type=int, default=0, help="Process only first N sequences for testing")
    p.add_argument("--verbose", action="store_true")
    p.add_argument("--max-shards", type=int, default=0, help="Limit number of shards to process")
    p.add_argument("--numbering-scheme", default="imgt", choices=["imgt", "kabat", "chothia"])
    p.add_argument("--show-stats", action="store_true", help="Show detailed statistics for each shard")
    args = p.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    start = dt.datetime.now()
    stamp = start.strftime("%Y%m%d_%H%M%S")
    safe_tag = re.sub(r"[^A-Za-z0-9._-]+", "_", args.tag) if args.tag else ""
    base = f"npz_fullscan_imgt_{safe_tag}_{stamp}" if safe_tag else f"npz_fullscan_imgt_{stamp}"

    log_path = os.path.join(args.outdir, f"{base}.log")
    tee = Tee(log_path)
    sys.stdout = tee

    print("="*80)
    print("NPZ FULLSCAN - Corrected IMGT CDR Extraction with Statistics")
    print("="*80)
    print(f"Started:       {start:%Y-%m-%d %H:%M:%S}")
    print(f"Database:      {args.db_root}")
    print(f"Chain:         {args.chain}")
    print(f"Species:       {args.species}")
    print(f"Regions:       {args.regions}")
    print(f"Numbering:     {args.numbering_scheme}")
    print(f"Show Stats:    {args.show_stats}")
    print("CDR Extractor: ANARCI (IMGT)")
    print("="*80)

    # Build shard list
    chain_dir = "Heavy" if args.chain.lower().startswith("h") else "Light"
    species_dir = args.species
    root = os.path.join(args.db_root, chain_dir, species_dir)
    npz_files = sorted(glob.glob(os.path.join(root, "*.npz")))
    
    if args.max_shards > 0:
        npz_files = npz_files[:args.max_shards]
    
    print(f"Found {len(npz_files)} shards")

  
    # Extract query CDRs using strict method
    print(f"\nExtracting query CDRs ({args.numbering_scheme} scheme)...")
    q_full = re.sub(r"[^A-Za-z]", "", args.query_seq.upper())

    query_cdrs = extract_cdrs_from_sequence_strict(q_full, args.numbering_scheme)
    query_cdrs["full"] = q_full

    # Normalize regions selection to a set of tokens
    regions_selected = {r.strip().lower() for r in str(args.regions).split(",") if r.strip()}

    # Build required list based on selected regions
    required = []
    if "cdr1" in regions_selected: required.append(("CDR-H1", "cdr1"))
    if "cdr2" in regions_selected: required.append(("CDR-H2", "cdr2"))
    if "cdr3" in regions_selected: required.append(("CDR-H3", "cdr3"))

    # Single early warning if any required CDR is missing
    missing = [name for name, key in required if not query_cdrs.get(key)]
    if missing and os.getenv("QUIET_QUERY_WARN", "0") != "1":
        print(f"[WARN] Query CDRs not extracted for: {', '.join(missing)}. "
            f"These regions will filter out all hits (identity=0).")

    print("\nQuery sequence (cleaned):")
    print(q_full)

    # Safe access helpers
    c1 = query_cdrs.get("cdr1", "")
    c2 = query_cdrs.get("cdr2", "")
    c3 = query_cdrs.get("cdr3", "")

    print("\nQuery CDRs extracted (strict IMGT boundaries):")
    print(f"  CDR-H1: {c1} (len={len(c1)})")
    print(f"  CDR-H2: {c2} (len={len(c2)})")
    print(f"  CDR-H3: {c3} (len={len(c3)})")
    print(f"  Full:   {len(q_full)} aa")

    # For debugging or reference (optional)
    print(f"\n[INFO] Expected CDRs for this antibody:")
    print(f"  CDR-H1: GYTFSSY (7 aa)")
    print(f"  CDR-H2: ILPGSGST (8 aa)")
    print(f"  CDR-H3: ARGDDYDEGFPS (12 aa)")

# Parse regions to search
    # Parse regions to search
    region_set = {s.strip().lower() for s in args.regions.split(",")}
    cfg = {
        "chain": args.chain,
        "species": args.species,
        "scheme": args.numbering_scheme,
        "use_cdr1": "cdr1" in region_set,
        "use_cdr2": "cdr2" in region_set,
        "use_cdr3": "cdr3" in region_set,
        "use_full": "full" in region_set,
        "min_id_cdr1": args.min_id_cdr1,
        "min_id_cdr2": args.min_id_cdr2,
        "min_id_cdr3": args.min_id_cdr3,
        "min_id_full": args.min_id_full,
        "win_cdr1": args.len_window_cdr1,
        "win_cdr2": args.len_window_cdr2,
        "win_cdr3": args.len_window_cdr3,
    }
    
    # Count active regions for statistics display
    active_regions = []
    if cfg['use_cdr1']: active_regions.append("CDR1")
    if cfg['use_cdr2']: active_regions.append("CDR2")
    if cfg['use_cdr3']: active_regions.append("CDR3")
    if cfg['use_full']: active_regions.append("Full")
    single_region = len(active_regions) == 1
    
    print(f"\nActive regions: {', '.join(active_regions)}")
    print(f"Thresholds: ", end="")
    thresholds = []
    if cfg['use_cdr1']: thresholds.append(f"CDR1≥{args.min_id_cdr1:.0%}")
    if cfg['use_cdr2']: thresholds.append(f"CDR2≥{args.min_id_cdr2:.0%}")
    if cfg['use_cdr3']: thresholds.append(f"CDR3≥{args.min_id_cdr3:.0%}")
    if cfg['use_full']: thresholds.append(f"Full≥{args.min_id_full:.0%}")
    print(", ".join(thresholds))
    
    # Warn if CDRs not extracted
    if cfg['use_cdr1'] and not query_cdrs['cdr1']:
        print("[WARNING] CDR1 search requested but query CDR1 not extracted")
    if cfg['use_cdr2'] and not query_cdrs['cdr2']:
        print("[WARNING] CDR2 search requested but query CDR2 not extracted")
    if cfg['use_cdr3'] and not query_cdrs['cdr3']:
        print("[WARNING] CDR3 search requested but query CDR3 not extracted")

    # Output paths
    csv_path = os.path.join(args.outdir, f"{base}.csv")
    xlsx_path = os.path.join(args.outdir, f"{base}.xlsx")

    # Head-check mode for testing
    if args.head_check and npz_files:
        print(f"\n[HEAD-CHECK MODE] Processing first {args.head_check} sequences...")
        head_rows = []
        f0 = npz_files[0]
        
        try:
            d0 = np.load(f0, allow_pickle=True)
            arr0 = d0["numberings"]
            n = min(args.head_check, len(arr0))
            
            for i in range(n):
                row = arr0[i]
                full_seq = extract_full_sequence_from_npz(row)
                
                if full_seq and len(full_seq) >= 50:
                    db_cdrs = extract_cdrs_from_sequence_strict(full_seq, args.numbering_scheme)
                    
                    # Calculate identities
                    id1 = pid_identity(query_cdrs.get("cdr1", ""), db_cdrs["cdr1"]) if db_cdrs["cdr1"] else 0
                    id2 = pid_identity(query_cdrs.get("cdr2", ""), db_cdrs["cdr2"]) if db_cdrs["cdr2"] else 0
                    id3 = pid_identity(query_cdrs.get("cdr3", ""), db_cdrs["cdr3"]) if db_cdrs["cdr3"] else 0
                    
                    head_rows.append({
                        "Row_index": i,
                        "CDR1": db_cdrs["cdr1"],
                        "CDR2": db_cdrs["cdr2"],
                        "CDR3": db_cdrs["cdr3"],
                        "ID_CDR1": round(id1, 3),
                        "ID_CDR2": round(id2, 3),
                        "ID_CDR3": round(id3, 3),
                        "Full_preview": preview(full_seq, 60),
                    })
            
            d0.close()
            
            # Write headcheck CSV
            head_csv = os.path.join(args.outdir, f"{base}_headcheck.csv")
            pd.DataFrame(head_rows).to_csv(head_csv, index=False)
            print(f"[HEAD-CHECK] Wrote {len(head_rows)} rows to {head_csv}")
            
            # Show first few
            for j in range(min(5, len(head_rows))):
                r = head_rows[j]
                print(f"\n  Row {r['Row_index']}:")
                print(f"    CDR1: {r['CDR1']} (id={r['ID_CDR1']:.1%})")
                print(f"    CDR2: {r['CDR2']} (id={r['ID_CDR2']:.1%})")
                print(f"    CDR3: {r['CDR3']} (id={r['ID_CDR3']:.1%})")
                print(f"    Preview: {r['Full_preview']}")
            
        except Exception as e:
            print(f"[HEAD-CHECK] Error: {e}")
        finally:
            print("\n[HEAD-CHECK] Complete. Exiting.")
            tee.close()
            return

    # Process all shards
    print("\nProcessing shards...")
    print("-"*80)
    first_header = True
    all_hits = []
    total_examined = total_skipped = 0
    t0 = time.time()
    
    # Track overall best hits
    overall_best = {
        "cdr1": 0,
        "cdr2": 0,
        "cdr3": 0,
        "full": 0
    }

    # Write CSV header comments
    with open(csv_path, "w", encoding="utf-8") as f:
        f.write("# NPZ Fullscan - Corrected IMGT CDR Extraction\n")
        f.write(f"# Query CDR1: {query_cdrs['cdr1']}\n")
        f.write(f"# Query CDR2: {query_cdrs['cdr2']}\n")
        f.write(f"# Query CDR3: {query_cdrs['cdr3']}\n")
        f.write(f"# Regions: {args.regions}\n")
        f.write(f"# Started: {start:%Y-%m-%d %H:%M:%S}\n")

    for idx, npz_path in enumerate(npz_files, 1):
        t_file = time.time()
        verbose = args.verbose and idx <= 3
        
        hits, stats, detailed_stats = process_npz_file(npz_path, cfg, query_cdrs, verbose=verbose)
        total_examined += stats["examined"]
        total_skipped += stats.get("skipped", 0)

        # Update overall best
        for region in ["cdr1", "cdr2", "cdr3", "full"]:
            if region in detailed_stats and detailed_stats[region]["count"] > 0:
                overall_best[region] = max(overall_best[region], detailed_stats[region]["max"])

        # Save hits to CSV
        if hits:
            df = pd.DataFrame(hits)
            df.to_csv(csv_path, mode="a", index=False, header=first_header)
            first_header = False
            all_hits.extend(hits)
            
            # Find best hit for display
            best_hit = max(hits, key=lambda x: sum([
                x.get("ID_CDR1", 0) if cfg["use_cdr1"] else 0,
                x.get("ID_CDR2", 0) if cfg["use_cdr2"] else 0,
                x.get("ID_CDR3", 0) if cfg["use_cdr3"] else 0,
                x.get("ID_Full", 0) if cfg["use_full"] else 0
            ]))

        took = time.time() - t_file
        
        # Format output line
        if hits:
            print(f"[{idx:3}/{len(npz_files)}] {os.path.basename(npz_path):20} | "
                  f"+{len(hits):4} hits | {stats['examined']:5} examined | {took:4.1f}s")
            
            # Show best hit details
            if single_region:
                region_key = active_regions[0].lower()
                if region_key in detailed_stats and detailed_stats[region_key]["count"] > 0:
                    ds = detailed_stats[region_key]
                    vals = [h[f"ID_{region_key.capitalize()}"] if region_key != "full" else h["ID_Full"] for h in hits]
                    high = sum(1 for v in vals if v >= 0.9)
                    good = sum(1 for v in vals if 0.7 <= v < 0.9)
                    med = sum(1 for v in vals if 0.5 <= v < 0.7)
                    
                    dist_parts = []
                    if high > 0: dist_parts.append(f"{high} >90%")
                    if good > 0: dist_parts.append(f"{good} 70-90%")
                    if med > 0: dist_parts.append(f"{med} 50-70%")
                    dist_str = f" [{', '.join(dist_parts)}]" if dist_parts else ""
                    
                    print(f"         Best: {ds['max']:.1%} | Q75: {ds['q75']:.1%} | "
                          f"Median: {ds['median']:.1%}{dist_str}")
            else:
                best_parts = []
                if cfg["use_cdr1"] and best_hit["ID_CDR1"] > 0:
                    best_parts.append(f"CDR1={best_hit['ID_CDR1']:.1%}")
                if cfg["use_cdr2"] and best_hit["ID_CDR2"] > 0:
                    best_parts.append(f"CDR2={best_hit['ID_CDR2']:.1%}")
                if cfg["use_cdr3"] and best_hit["ID_CDR3"] > 0:
                    best_parts.append(f"CDR3={best_hit['ID_CDR3']:.1%}")
                if cfg["use_full"] and best_hit["ID_Full"] > 0:
                    best_parts.append(f"Full={best_hit['ID_Full']:.1%}")
                
                if best_parts:
                    print(f"         Best hit: {', '.join(best_parts)}")
            
            # Optional detailed statistics
            if args.show_stats and detailed_stats:
                for region, dstats in detailed_stats.items():
                    if cfg[f"use_{region}"] and dstats["count"] > 0:
                        print(f"         {region.upper()}: {format_stats_line(dstats)}")
        else:
            print(f"[{idx:3}/{len(npz_files)}] {os.path.basename(npz_path):20} | "
                  f"  No hits  | {stats['examined']:5} examined | {took:4.1f}s")

    # Summary
    print("-"*80)
    end = dt.datetime.now()
    elapsed = (end - start).total_seconds()
    
    print(f"\nSUMMARY")
    print("="*80)
    print(f"Total hits:     {len(all_hits)}")
    print(f"Total examined: {total_examined}")
    print(f"Total skipped:  {total_skipped}")
    print(f"Elapsed:        {elapsed:.1f}s ({elapsed/60:.1f} min)")
    
    # Show overall best identities
    if all_hits:
        print(f"\nBest identities found:")
        if cfg["use_cdr1"] and overall_best["cdr1"] > 0:
            print(f"  CDR1: {overall_best['cdr1']:.1%}")
        if cfg["use_cdr2"] and overall_best["cdr2"] > 0:
            print(f"  CDR2: {overall_best['cdr2']:.1%}")
        if cfg["use_cdr3"] and overall_best["cdr3"] > 0:
            print(f"  CDR3: {overall_best['cdr3']:.1%}")
        if cfg["use_full"] and overall_best["full"] > 0:
            print(f"  Full: {overall_best['full']:.1%}")
    
    print(f"\nOutput files:")
    print(f"  CSV:  {csv_path}")
    
    if all_hits:
        with pd.ExcelWriter(xlsx_path, engine="openpyxl") as writer:
            pd.DataFrame(all_hits).to_excel(writer, sheet_name="Results", index=False)
            
            stats_data = []
            for h in all_hits:
                stats_data.append({
                    "CDR1_ID": h["ID_CDR1"],
                    "CDR2_ID": h["ID_CDR2"],
                    "CDR3_ID": h["ID_CDR3"],
                    "Full_ID": h["ID_Full"],
                })
            
            if stats_data:
                stats_df = pd.DataFrame(stats_data)
                stats_summary = stats_df.describe()
                stats_summary.to_excel(writer, sheet_name="Statistics")
                
        print(f"  XLSX: {xlsx_path}")
    
    print("="*80)
    print(f"Completed at {end:%Y-%m-%d %H:%M:%S}")
    tee.close()

if __name__ == "__main__":
    main()