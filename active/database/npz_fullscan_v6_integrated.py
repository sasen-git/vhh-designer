#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
NPZ FULLSCAN v6 INTEGRATED - Combined Interactive Analysis & Full Scan

This version integrates:
1. Interactive threshold analysis (sampling mode)
2. Full exhaustive scan with all v6 improvements

NEW: Asks if you want to sample first to determine optimal thresholds!

Usage:
  python3 npz_fullscan_v6_integrated.py --db-root /path/to/db --query-seq "EIQLQQ..."
  
  The script will ask if you want to:
  1. Run a quick analysis to determine optimal thresholds (recommended)
  2. Proceed directly with full scan
"""

import os
import re
import time
import glob
import signal
import argparse
import datetime as dt
import sys
import json
import platform
from typing import Dict, List, Tuple, Optional
from pathlib import Path
import random

import numpy as np
import pandas as pd
from tqdm import tqdm

# Optional: ANARCI for query extraction
try:
    from anarci import run_anarci
    _ANARCI_AVAILABLE = True
except Exception:
    run_anarci = None
    _ANARCI_AVAILABLE = False

# ------------------------
# Constants
# ------------------------
AA_PAD = 124  
AA_GAP = 0

# Amino acid similarity groups for position constraints
SIMILAR_GROUPS = {
    'H': ['H', 'K', 'R'],      # Positive charged
    'K': ['H', 'K', 'R'],      # Positive charged  
    'R': ['H', 'K', 'R'],      # Positive charged
    'D': ['D', 'E'],           # Negative charged
    'E': ['D', 'E'],           # Negative charged
    'S': ['S', 'T'],           # Polar uncharged (small)
    'T': ['S', 'T'],           # Polar uncharged (small)
    'N': ['N', 'Q'],           # Polar uncharged (amide)
    'Q': ['N', 'Q'],           # Polar uncharged (amide)
    'L': ['L', 'I', 'V'],      # Hydrophobic aliphatic
    'I': ['L', 'I', 'V'],      # Hydrophobic aliphatic
    'V': ['L', 'I', 'V'],      # Hydrophobic aliphatic
    'F': ['F', 'Y', 'W'],      # Aromatic
    'Y': ['F', 'Y', 'W'],      # Aromatic
    'W': ['F', 'Y', 'W'],      # Aromatic
    'A': ['A', 'G'],           # Small nonpolar
    'G': ['A', 'G'],           # Small nonpolar
    'C': ['C'],                # Cysteine (special)
    'M': ['M'],                # Methionine (special)
    'P': ['P'],                # Proline (special)
}

# ------------------------
# Timeout Handler
# ------------------------
class _Timeout(Exception):
    pass

def _alarm_handler(signum, frame):
    raise _Timeout()

# ------------------------
# Decoder functions
# ------------------------
def _ints_to_aa(int_row):
    """Ultra-robust decoder for NPZ integer arrays."""
    out = []
    for v in int_row:
        if v is None:
            continue
            
        if isinstance(v, (str, bytes)):
            if isinstance(v, bytes):
                try:
                    v = v.decode('utf-8', errors='ignore')
                except:
                    continue
            v = str(v).strip()
            if len(v) == 1 and v.upper() in 'ACDEFGHIKLMNPQRSTVWY':
                out.append(v.upper())
            continue
        
        try:
            v_str = str(v).strip()
            if len(v_str) == 1 and v_str.upper() in 'ACDEFGHIKLMNPQRSTVWY':
                out.append(v_str.upper())
                continue
        except:
            pass
        
        try:
            iv = int(v)
            if iv == AA_GAP or iv == AA_PAD:
                continue
            elif 65 <= iv <= 90:
                out.append(chr(iv))
            elif 97 <= iv <= 122:
                out.append(chr(iv).upper())
        except (ValueError, TypeError):
            continue
    return "".join(out)

# ------------------------
# CDR/FR Extraction
# ------------------------
def extract_cdr1_from_npz(int_row):
    return _ints_to_aa(int_row[39:46])

def extract_cdr2_from_npz(int_row):
    part1 = _ints_to_aa(int_row[66:75])
    part2 = _ints_to_aa(int_row[85:91])
    return part1 + part2

def extract_cdr3_from_npz(int_row):
    cdr3_region = int_row[150:190]
    cdr3_aa = _ints_to_aa(cdr3_region)
    if cdr3_aa and cdr3_aa[0] == 'C':
        w_pos = cdr3_aa.find('W')
        if w_pos > 0:
            return cdr3_aa[:w_pos]
    return cdr3_aa

def extract_frameworks_from_npz(int_row):
    """Extract framework regions (everything except CDRs)."""
    fr1 = _ints_to_aa(int_row[0:27])
    fr2 = _ints_to_aa(int_row[46:65])
    fr3 = _ints_to_aa(int_row[91:150])
    fr4 = _ints_to_aa(int_row[190:])
    return fr1 + fr2 + fr3 + fr4

def extract_whole_from_npz(int_row):
    return _ints_to_aa(int_row[16:190])

def heuristic_cdrh3(full_seq: str) -> str:
    """Backup CDR3 extraction using heuristics."""
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
# Query CDR extraction via ANARCI
# ------------------------
def extract_cdrs_from_query_sequence(seq: str, scheme: str = "imgt", timeout_s: int = 5) -> Dict[str, str]:
    """Extract CDRs from query using ANARCI."""
    if not _ANARCI_AVAILABLE:
        h3 = heuristic_cdrh3(seq)
        return {"cdr1": "", "cdr2": "", "cdr3": h3}
    
    seq_clean = seq.upper().replace("-", "").replace(".", "")
    
    if hasattr(signal, 'SIGALRM'):
        signal.signal(signal.SIGALRM, _alarm_handler)
        signal.alarm(timeout_s)
    
    try:
        res = run_anarci([("H", seq_clean)], scheme=scheme, allowed_species=None)
        if hasattr(signal, 'SIGALRM'):
            signal.alarm(0)
        
        if not res or not res[1] or not res[1][0] or not res[1][0][0]:
            return {"cdr1": "", "cdr2": "", "cdr3": heuristic_cdrh3(seq_clean)}
        
        numbering = res[1][0][0][0]
        
        def pos_val(pos_tuple):
            pos = pos_tuple[0]
            ins = pos_tuple[1] if len(pos_tuple) > 1 else ""
            if isinstance(ins, str) and ins.strip():
                return float(f"{pos}.{ord(ins.upper())-64:02d}")
            return float(pos)
        
        def grab(start, end):
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
        
        return {
            "cdr1": grab(27.0, 38.9),
            "cdr2": grab(56.0, 65.9),
            "cdr3": grab(105.0, 117.9)
        }
    except Exception:
        if hasattr(signal, 'SIGALRM'):
            signal.alarm(0)
        return {"cdr1": "", "cdr2": "", "cdr3": heuristic_cdrh3(seq_clean)}

# ------------------------
# Alignment & Identity
# ------------------------
def edit_dist(s1: str, s2: str) -> int:
    """Fast edit distance calculation."""
    len1, len2 = len(s1), len(s2)
    
    if len1 < len2:
        s1, s2 = s2, s1
        len1, len2 = len2, len1
    
    if len2 == 0:
        return len1
    
    prev_row = list(range(len2 + 1))
    
    for i, c1 in enumerate(s1):
        curr_row = [i + 1]
        for j, c2 in enumerate(s2):
            insertions = prev_row[j + 1] + 1
            deletions = curr_row[j] + 1
            substitutions = prev_row[j] + (c1 != c2)
            curr_row.append(min(insertions, deletions, substitutions))
        prev_row = curr_row
    
    return prev_row[-1]

def calc_id_by_region(qval: str, dval: str) -> float:
    """Calculate identity for a region."""
    if not qval or not dval:
        return 0.0
    max_len = max(len(qval), len(dval))
    if max_len == 0:
        return 0.0
    dist = edit_dist(qval, dval)
    return 1.0 - (dist / max_len)

# ------------------------
# File organization helpers
# ------------------------
def get_species_from_path(db_path: str) -> str:
    """Extract species from database path."""
    path_parts = Path(db_path).parts
    for part in reversed(path_parts):
        if part.lower() in ['human', 'camel', 'mouse', 'rabbit', 'rat', 'rhesus', 'humanised', 'humanized']:
            return part.lower()
    return 'unknown'

def create_output_structure(outdir: str, species: str, regions: str, thresholds: Dict,
                           tag: Optional[str] = None, mode: str = "fullscan") -> Tuple[str, str]:
    """Create organized output structure and generate filename."""
    date_str = dt.datetime.now().strftime("%Y%m%d")
    folder_name = f"{date_str}_{species}"
    full_folder_path = os.path.join(outdir, folder_name)
    os.makedirs(full_folder_path, exist_ok=True)
    
    timestamp = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    
    # Parse regions
    regions_list = [r.strip().lower() for r in regions.split(",")]
    if 'full' in regions_list:
        region_str = 'full'
    elif 'frameworks' in regions_list:
        region_str = 'frameworks'
    elif all(r in ['cdr1', 'cdr2', 'cdr3'] for r in regions_list):
        region_str = 'cdrs'
    else:
        region_str = '_'.join(regions_list)
    
    # Build identity string
    identity_parts = []
    for region in regions_list:
        if region == 'full' and thresholds.get('min_id_full'):
            identity_parts.append(f"full-{thresholds['min_id_full']}")
        elif region == 'cdr1' and thresholds.get('min_id_cdr1'):
            identity_parts.append(f"c1-{thresholds['min_id_cdr1']}")
        elif region == 'cdr2' and thresholds.get('min_id_cdr2'):
            identity_parts.append(f"c2-{thresholds['min_id_cdr2']}")
        elif region == 'cdr3' and thresholds.get('min_id_cdr3'):
            identity_parts.append(f"c3-{thresholds['min_id_cdr3']}")
        elif region == 'frameworks' and thresholds.get('min_id_frameworks'):
            identity_parts.append(f"fr-{thresholds['min_id_frameworks']}")
    
    identity_str = '_'.join(identity_parts) if identity_parts else 'no_threshold'
    
    # Add mode prefix
    mode_prefix = "analysis_" if mode == "analysis" else ""
    tag_str = f"_{tag}" if tag else ""
    
    # Final filename
    filename = f"{mode_prefix}{timestamp}_{species}_{region_str}_{identity_str}{tag_str}.csv"
    full_path = os.path.join(full_folder_path, filename)
    
    return full_path, full_folder_path

# ============================================================================
# INTERACTIVE ANALYSIS FUNCTIONS
# ============================================================================

def sample_and_analyze_shard(npz_path: str, query_data: dict, sample_size: int,
                             regions: List[str], use_length_filter: bool = False,
                             length_windows: Dict = None) -> pd.DataFrame:
    """Sample sequences from a shard for analysis."""
    try:
        arr = np.load(npz_path, allow_pickle=True)
        if "numberings" not in arr:
            arr.close()
            return pd.DataFrame()
        
        numberings = arr["numberings"]
        total_seqs = len(numberings)
        
        if total_seqs == 0:
            arr.close()
            return pd.DataFrame()
        
        # Random sample
        sample_indices = np.random.choice(total_seqs, min(sample_size, total_seqs), replace=False)
        
        results = []
        for idx in sample_indices:
            row = numberings[idx]
            
            # Extract regions
            db_regions = {}
            if 'full' in regions:
                db_regions['full'] = extract_whole_from_npz(row)
            if 'cdr1' in regions:
                db_regions['cdr1'] = extract_cdr1_from_npz(row)
            if 'cdr2' in regions:
                db_regions['cdr2'] = extract_cdr2_from_npz(row)
            if 'cdr3' in regions:
                db_regions['cdr3'] = extract_cdr3_from_npz(row)
            if 'frameworks' in regions:
                db_regions['frameworks'] = extract_frameworks_from_npz(row)
            
            # Apply length filters if configured
            if use_length_filter and length_windows:
                skip = False
                for cdr in ['cdr1', 'cdr2', 'cdr3']:
                    if cdr in length_windows and cdr in db_regions:
                        window = length_windows[cdr]
                        query_len = len(query_data.get(cdr, ""))
                        db_len = len(db_regions[cdr])
                        if abs(db_len - query_len) > window:
                            skip = True
                            break
                if skip:
                    continue
            
            # Calculate identities
            result = {"shard": os.path.basename(npz_path), "index": int(idx)}
            
            for region in regions:
                if region in db_regions and region in query_data:
                    result[f"id_{region}"] = calc_id_by_region(
                        query_data[region], db_regions[region]
                    )
                    result[f"{region}_seq"] = db_regions[region]
            
            # Calculate combined scores
            id_values = [v for k, v in result.items() if k.startswith("id_")]
            if id_values:
                result["avg_id"] = np.mean(id_values)
                result["min_id"] = np.min(id_values)
            
            results.append(result)
        
        arr.close()
        return pd.DataFrame(results)
        
    except Exception as e:
        print(f"Error processing shard: {e}")
        return pd.DataFrame()

def run_interactive_analysis(db_root: str, query_seq: str, regions: List[str],
                            sample_size: int = 1000, target_hits: int = 10000) -> Dict:
    """Run interactive analysis to determine optimal thresholds."""
    
    print("\n" + "="*80)
    print("INTERACTIVE THRESHOLD ANALYSIS")
    print("="*80)
    print(f"Sampling {sample_size} sequences per shard...")
    print(f"Target: Find thresholds to get ~{target_hits:,} hits")
    print("="*80 + "\n")
    
    # Extract query CDRs
    q_full = re.sub(r"[^A-Za-z]", "", query_seq.upper())
    query_cdrs = {"full": q_full}
    
    if set(regions) != {"full"}:
        qc = extract_cdrs_from_query_sequence(q_full)
        query_cdrs.update(qc)
    
    if 'frameworks' in regions:
        query_cdrs['frameworks'] = q_full  # Placeholder
    
    # Find shards
    npz_files = sorted(glob.glob(os.path.join(db_root, "*.npz")))
    
    if not npz_files:
        print("ERROR: No NPZ files found!")
        return {}
    
    print(f"Found {len(npz_files)} shards")
    print(f"Total sequences to sample: {len(npz_files) * sample_size:,}\n")
    
    # Sample from shards
    all_samples = []
    sample_bar = tqdm(npz_files[:min(10, len(npz_files))], desc="Sampling shards")
    
    for npz_file in sample_bar:
        df = sample_and_analyze_shard(npz_file, query_cdrs, sample_size, regions)
        if not df.empty:
            all_samples.append(df)
    
    if not all_samples:
        print("ERROR: No samples collected!")
        return {}
    
    combined = pd.concat(all_samples, ignore_index=True)
    
    # Analyze distributions
    print(f"\n{'='*80}")
    print("IDENTITY DISTRIBUTION ANALYSIS")
    print(f"{'='*80}\n")
    
    recommendations = {}
    
    for region in regions:
        col = f"id_{region}"
        if col not in combined.columns:
            continue
        
        values = combined[col].values
        
        print(f"\n{region.upper()} Statistics:")
        print(f"  Max:    {values.max():.3f} ({values.max()*100:.1f}%)")
        print(f"  Q95:    {np.percentile(values, 95):.3f} ({np.percentile(values, 95)*100:.1f}%)")
        print(f"  Q90:    {np.percentile(values, 90):.3f} ({np.percentile(values, 90)*100:.1f}%)")
        print(f"  Q75:    {np.percentile(values, 75):.3f} ({np.percentile(values, 75)*100:.1f}%)")
        print(f"  Median: {np.median(values):.3f} ({np.median(values)*100:.1f}%)")
        print(f"  Mean:   {values.mean():.3f} ({values.mean()*100:.1f}%)")
        
        # Estimate threshold for target hits
        total_sampled = len(combined)
        scale_factor = len(npz_files) * 5_000_000 / total_sampled  # Approximate
        needed_in_sample = target_hits / scale_factor
        
        if needed_in_sample < len(values):
            # Find threshold that would give us target_hits
            sorted_vals = np.sort(values)[::-1]
            threshold_idx = min(int(needed_in_sample), len(sorted_vals) - 1)
            recommended_threshold = sorted_vals[threshold_idx]
            
            print(f"\n  RECOMMENDED THRESHOLD: {recommended_threshold:.3f} ({recommended_threshold*100:.1f}%)")
            print(f"  This should give approximately {target_hits:,} hits")
            
            recommendations[region] = recommended_threshold
    
    # Show top matches
    print(f"\n{'='*80}")
    print("TOP 10 MATCHES IN SAMPLE")
    print(f"{'='*80}\n")
    
    if 'avg_id' in combined.columns:
        top10 = combined.nlargest(10, 'avg_id')
        for i, (_, row) in enumerate(top10.iterrows(), 1):
            print(f"#{i:2d} Shard: {row['shard'][:40]:40s} (idx {row['index']:,})")
            for region in regions:
                col = f"id_{region}"
                if col in row:
                    print(f"    {region:10s}: {row[col]*100:5.1f}%")
            print(f"    Average:    {row['avg_id']*100:5.1f}%")
            print()
    
    return recommendations

# ============================================================================
# FULL SCAN FUNCTIONS (from v6)
# ============================================================================

def write_metadata_block(file_handle, args, query_cdrs, start_time, config):
    """Write comprehensive metadata block to CSV file."""
    file_handle.write("# ============================================\n")
    file_handle.write("# NPZ FULLSCAN v6 INTEGRATED - RUN METADATA\n")
    file_handle.write("# ============================================\n")
    file_handle.write(f"# Run started: {start_time.strftime('%Y-%m-%d %H:%M:%S')}\n")
    file_handle.write(f"# Script: {os.path.basename(__file__)}\n")
    file_handle.write(f"# Command: {' '.join(sys.argv)}\n")
    file_handle.write(f"# Platform: {platform.platform()}\n")
    file_handle.write(f"# Python: {sys.version.split()[0]}\n")
    file_handle.write("#\n")
    file_handle.write("# DATABASE INFORMATION\n")
    file_handle.write(f"# Database path: {args.db_root}\n")
    file_handle.write(f"# Species: {get_species_from_path(args.db_root)}\n")
    file_handle.write("#\n")
    file_handle.write("# QUERY SEQUENCE\n")
    file_handle.write(f"# Length: {len(query_cdrs['full'])} aa\n")
    file_handle.write(f"# Sequence: {query_cdrs['full']}\n")
    file_handle.write(f"# CDR1: {query_cdrs.get('cdr1', 'N/A')} (len={len(query_cdrs.get('cdr1', ''))})\n")
    file_handle.write(f"# CDR2: {query_cdrs.get('cdr2', 'N/A')} (len={len(query_cdrs.get('cdr2', ''))})\n")
    file_handle.write(f"# CDR3: {query_cdrs.get('cdr3', 'N/A')} (len={len(query_cdrs.get('cdr3', ''))})\n")
    file_handle.write("#\n")
    file_handle.write("# SEARCH PARAMETERS\n")
    file_handle.write(f"# Regions searched: {args.regions}\n")
    file_handle.write(f"# Numbering scheme: {args.scheme}\n")
    
    # Write thresholds
    for key, value in config.items():
        if key.startswith('min_id_') and value is not None:
            region = key.replace('min_id_', '').upper()
            file_handle.write(f"# {region} min identity: {value}\n")
    
    if config.get('use_length_filter'):
        file_handle.write("#\n# LENGTH FILTERING\n")
        for cdr, window in config.get('length_windows', {}).items():
            file_handle.write(f"# {cdr.upper()} length window: Â±{window} aa\n")
    
    if config.get('interactive_analysis_run'):
        file_handle.write("#\n# INTERACTIVE ANALYSIS\n")
        file_handle.write(f"# Analysis performed: Yes\n")
        file_handle.write(f"# Recommended thresholds applied: {config.get('used_recommendations', False)}\n")
    
    if args.tag:
        file_handle.write(f"#\n# Tag: {args.tag}\n")
    
    file_handle.write("# ============================================\n")
    file_handle.write("#\n")

def process_npz_file_fullscan(npz_path: str, config: dict, query_data: dict,
                              shard_idx: int, total_shards: int) -> dict:
    """Process single NPZ file for full scan."""
    
    shard_name = os.path.basename(npz_path)
    stats = {
        "shard": shard_name,
        "n_total": 0,
        "n_hits": 0,
        "n_skipped": 0,
        "n_filtered_length": 0,
        "n_filtered_seq_length": 0,
        "n_filtered_cdr_length": 0,
        "max_id": 0.0,
        "min_id": 1.0,
        "q75": 0.0,
        "time_s": 0.0,
        "best_seq": "",
        "worst_seq": ""
    }
    
    t0 = time.time()
    
    try:
        arr = np.load(npz_path, allow_pickle=True)
        if "numberings" not in arr:
            arr.close()
            stats["n_skipped"] = -1
            stats["time_s"] = time.time() - t0
            return stats
        
        numberings = arr["numberings"]
        n_total = len(numberings)
        stats["n_total"] = n_total
        
        if n_total == 0:
            arr.close()
            stats["time_s"] = time.time() - t0
            return stats
    except Exception:
        stats["n_skipped"] = -1
        stats["time_s"] = time.time() - t0
        return stats
    
    hits = []
    regions_to_check = {r.strip().lower() for r in str(config["regions"]).split(",") if r.strip()}
    
    # Determine if we're in full mode for optimizations
    is_full_mode = regions_to_check == {"full"}
    query_full_len = len(query_data.get('full', ''))
    
    inner_bar = tqdm(total=n_total, desc=f"  Shard {shard_idx}/{total_shards}",
                     ncols=100, position=1, leave=False)
    
    for idx, row in enumerate(numberings):
        inner_bar.update(1)
        
        # Length optimizations for full mode
        if is_full_mode and config.get("enable_seq_length_prefilter", True):
            db_full_approx = extract_whole_from_npz(row)
            db_full_len = len(db_full_approx)
            
            seq_len_threshold = config.get("seq_length_threshold", 10)
            if abs(db_full_len - query_full_len) > seq_len_threshold:
                stats["n_filtered_seq_length"] += 1
                continue
            
            if config.get("use_length_filter") and config.get("enable_cdr_length_prefilter", True):
                skip = False
                temp_cdrs = {
                    'cdr1': extract_cdr1_from_npz(row),
                    'cdr2': extract_cdr2_from_npz(row),
                    'cdr3': extract_cdr3_from_npz(row)
                }
                
                for cdr in ['cdr1', 'cdr2', 'cdr3']:
                    if cdr in config.get("length_windows", {}):
                        window = config["length_windows"][cdr]
                        query_len = config["query_lengths"][cdr]
                        db_len = len(temp_cdrs[cdr])
                        if abs(db_len - query_len) > window:
                            stats["n_filtered_cdr_length"] += 1
                            skip = True
                            break
                
                if skip:
                    continue
        
        # Extract regions
        db_regions = {}
        if 'full' in regions_to_check:
            if is_full_mode and 'db_full_approx' in locals():
                db_regions['full'] = db_full_approx
            else:
                db_regions['full'] = extract_whole_from_npz(row)
        if 'cdr1' in regions_to_check:
            db_regions['cdr1'] = extract_cdr1_from_npz(row)
        if 'cdr2' in regions_to_check:
            db_regions['cdr2'] = extract_cdr2_from_npz(row)
        if 'cdr3' in regions_to_check:
            db_regions['cdr3'] = extract_cdr3_from_npz(row)
        if 'frameworks' in regions_to_check:
            db_regions['frameworks'] = extract_frameworks_from_npz(row)
        
        # Standard length filters for non-full modes
        if not is_full_mode and config.get("use_length_filter"):
            skip = False
            for cdr in ['cdr1', 'cdr2', 'cdr3']:
                if cdr in config.get("length_windows", {}):
                    window = config["length_windows"][cdr]
                    query_len = config["query_lengths"][cdr]
                    db_len = len(db_regions.get(cdr, ""))
                    if abs(db_len - query_len) > window:
                        stats["n_filtered_length"] += 1
                        skip = True
                        break
            if skip:
                continue
        
        # Calculate identities
        identities = {}
        for region in regions_to_check:
            if region in db_regions and region in query_data:
                identities[region] = calc_id_by_region(query_data[region], db_regions[region])
        
        # Check thresholds
        passes = True
        for region in regions_to_check:
            min_key = f"min_id_{region}"
            if config.get(min_key) is not None and region in identities:
                if identities[region] < config[min_key]:
                    passes = False
                    break
        
        if passes and identities:
            avg_id = np.mean(list(identities.values()))
            hits.append({
                'idx': idx,
                'identities': identities,
                'avg_id': avg_id,
                'sequence': db_regions.get('full', '')
            })
            
            if avg_id > stats["max_id"]:
                stats["max_id"] = avg_id
                stats["best_seq"] = db_regions.get('full', '')[:50] + "..."
            if avg_id < stats["min_id"]:
                stats["min_id"] = avg_id
                stats["worst_seq"] = db_regions.get('full', '')[:50] + "..."
    
    inner_bar.close()
    arr.close()
    
    stats["n_hits"] = len(hits)
    if hits:
        avg_ids = [h['avg_id'] for h in hits]
        stats["q75"] = np.percentile(avg_ids, 75)
    
    stats["time_s"] = time.time() - t0
    return stats

# ============================================================================
# MAIN INTEGRATED FUNCTION
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="NPZ Fullscan v6 INTEGRATED - Interactive Analysis + Full Scan",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
This integrated version allows you to:
1. Run a quick sampling analysis to determine optimal thresholds
2. Apply those thresholds to a full exhaustive scan

The script will ask if you want to run the analysis first (recommended).
        """
    )
    
    # Required arguments
    parser.add_argument("--db-root", required=True, help="Path to database directory")
    parser.add_argument("--query-seq", required=True, help="Query heavy chain sequence")
    
    # Search configuration
    parser.add_argument("--regions", default="cdr3",
                       help="Comma-separated regions: full,cdr1,cdr2,cdr3,frameworks")
    parser.add_argument("--scheme", default="imgt", help="Numbering scheme")
    
    # Identity thresholds (optional - can be determined by analysis)
    parser.add_argument("--min-id-full", type=float, default=None)
    parser.add_argument("--min-id-cdr1", type=float, default=None)
    parser.add_argument("--min-id-cdr2", type=float, default=None)
    parser.add_argument("--min-id-cdr3", type=float, default=None)
    parser.add_argument("--min-id-frameworks", type=float, default=None)
    
    # Length filtering
    parser.add_argument("--cdr1-len-window", type=int, default=None)
    parser.add_argument("--cdr2-len-window", type=int, default=None)
    parser.add_argument("--cdr3-len-window", type=int, default=None)
    
    # Interactive analysis parameters
    parser.add_argument("--sample-size", type=int, default=1000,
                       help="Sequences to sample per shard for analysis")
    parser.add_argument("--target-hits", type=int, default=10000,
                       help="Target number of hits for threshold recommendation")
    
    # Skip interactive prompt
    parser.add_argument("--skip-interactive", action="store_true",
                       help="Skip interactive analysis prompt and go straight to full scan")
    parser.add_argument("--analysis-only", action="store_true",
                       help="Only run analysis, don't do full scan")
    
    # Output
    parser.add_argument("--outdir", default="./results")
    parser.add_argument("--tag", default="")
    
    # Optimization controls
    parser.add_argument("--disable-seq-length-prefilter", action="store_true")
    parser.add_argument("--disable-cdr-length-prefilter", action="store_true")
    parser.add_argument("--seq-length-threshold", type=int, default=10)
    
    args = parser.parse_args()
    
    # Track start time
    start_time = dt.datetime.now()
    
    print("\n" + "="*80)
    print("NPZ FULLSCAN v6 INTEGRATED")
    print("Interactive Analysis + Full Scan")
    print("="*80)
    print(f"Started: {start_time.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Database: {args.db_root}")
    print(f"Species: {get_species_from_path(args.db_root)}")
    print("="*80 + "\n")
    
    # Parse regions
    regions_list = [r.strip().lower() for r in args.regions.split(",") if r.strip()]
    
    # Extract query sequence and CDRs
    q_full = re.sub(r"[^A-Za-z]", "", args.query_seq.upper())
    query_cdrs = {"full": q_full}
    
    if set(regions_list) != {"full"}:
        qc = extract_cdrs_from_query_sequence(q_full, args.scheme)
        query_cdrs.update(qc)
    
    if 'frameworks' in regions_list:
        query_cdrs['frameworks'] = q_full
    
    print("Query sequence analysis:")
    print(f"  Full length: {len(q_full)} aa")
    print(f"  CDR1: {query_cdrs.get('cdr1', 'N/A')} (len={len(query_cdrs.get('cdr1', ''))})")
    print(f"  CDR2: {query_cdrs.get('cdr2', 'N/A')} (len={len(query_cdrs.get('cdr2', ''))})")
    print(f"  CDR3: {query_cdrs.get('cdr3', 'N/A')} (len={len(query_cdrs.get('cdr3', ''))})")
    print()
    
    # Interactive decision
    recommendations = {}
    used_recommendations = False
    
    if not args.skip_interactive:
        print("="*80)
        print("INTERACTIVE THRESHOLD ANALYSIS")
        print("="*80)
        print("\nWould you like to run a quick analysis to determine optimal thresholds?")
        print("This will sample sequences to recommend identity thresholds before the full scan.")
        print("\nOptions:")
        print("  1. Yes - Run analysis first (recommended)")
        print("  2. No  - Proceed directly to full scan")
        print("  3. Analysis only - Just run analysis, no full scan")
        
        while True:
            choice = input("\nYour choice (1/2/3): ").strip()
            if choice in ['1', '2', '3']:
                break
            print("Please enter 1, 2, or 3")
        
        if choice == '1' or choice == '3':
            # Run interactive analysis
            recommendations = run_interactive_analysis(
                args.db_root, args.query_seq, regions_list,
                args.sample_size, args.target_hits
            )
            
            if recommendations and choice == '1':
                print("\n" + "="*80)
                print("APPLY RECOMMENDED THRESHOLDS?")
                print("="*80)
                print("\nRecommended thresholds based on analysis:")
                for region, threshold in recommendations.items():
                    print(f"  {region.upper()}: {threshold:.3f} ({threshold*100:.1f}%)")
                
                apply = input("\nApply these thresholds to full scan? (y/n): ").strip().lower()
                if apply == 'y':
                    # Apply recommendations
                    for region, threshold in recommendations.items():
                        setattr(args, f"min_id_{region}", threshold)
                    used_recommendations = True
                    print("âœ“ Recommended thresholds applied")
                else:
                    print("âœ— Using manual thresholds (if specified)")
            
            if choice == '3':  # Analysis only
                print("\n" + "="*80)
                print("ANALYSIS COMPLETE")
                print("="*80)
                print("Recommended thresholds for full scan:")
                for region, threshold in recommendations.items():
                    print(f"  --min-id-{region} {threshold:.3f}")
                print("\nAnalysis complete. Exiting without full scan.")
                return
    
    # Build config for full scan
    cfg = {
        "scheme": args.scheme,
        "regions": args.regions,
        "min_id_full": args.min_id_full,
        "min_id_cdr1": args.min_id_cdr1,
        "min_id_cdr2": args.min_id_cdr2,
        "min_id_cdr3": args.min_id_cdr3,
        "min_id_frameworks": args.min_id_frameworks,
        "use_length_filter": False,
        "length_windows": {},
        "query_lengths": {
            "cdr1": len(query_cdrs.get("cdr1", "")),
            "cdr2": len(query_cdrs.get("cdr2", "")),
            "cdr3": len(query_cdrs.get("cdr3", ""))
        },
        "enable_seq_length_prefilter": not args.disable_seq_length_prefilter,
        "enable_cdr_length_prefilter": not args.disable_cdr_length_prefilter,
        "seq_length_threshold": args.seq_length_threshold,
        "interactive_analysis_run": bool(recommendations),
        "used_recommendations": used_recommendations
    }
    
    # Add length windows
    if args.cdr1_len_window is not None:
        cfg["use_length_filter"] = True
        cfg["length_windows"]["cdr1"] = args.cdr1_len_window
    if args.cdr2_len_window is not None:
        cfg["use_length_filter"] = True
        cfg["length_windows"]["cdr2"] = args.cdr2_len_window
    if args.cdr3_len_window is not None:
        cfg["use_length_filter"] = True
        cfg["length_windows"]["cdr3"] = args.cdr3_len_window
    
    # Show full scan configuration
    print("\n" + "="*80)
    print("FULL SCAN CONFIGURATION")
    print("="*80)
    print(f"Regions: {', '.join(regions_list)}")
    
    # Show thresholds
    has_thresholds = False
    for region in regions_list:
        threshold = cfg.get(f"min_id_{region}")
        if threshold is not None:
            print(f"  {region.upper()} min identity: {threshold:.3f} ({threshold*100:.1f}%)")
            has_thresholds = True
    
    if not has_thresholds:
        print("  No identity thresholds set (will return all matches)")
    
    if cfg["use_length_filter"]:
        print("\nLength filtering:")
        for cdr, window in cfg["length_windows"].items():
            print(f"  {cdr.upper()}: Â±{window} aa")
    
    print("="*80 + "\n")
    
    # Find shards
    npz_files = sorted(glob.glob(os.path.join(args.db_root, "*.npz")))
    
    if not npz_files:
        print("ERROR: No NPZ files found!")
        return
    
    print(f"Found {len(npz_files)} shards to process\n")
    
    # Create output structure
    species = get_species_from_path(args.db_root)
    thresholds = {
        f'min_id_{r}': cfg.get(f'min_id_{r}')
        for r in regions_list
    }
    
    csv_path, output_folder = create_output_structure(
        args.outdir, species, args.regions, thresholds, args.tag
    )
    
    print(f"Output folder: {output_folder}")
    print(f"Output file: {os.path.basename(csv_path)}\n")
    
    # Ensure output directory exists
    os.makedirs(args.outdir, exist_ok=True)
    
    # Write CSV header
    with open(csv_path, "w", encoding="utf-8") as fh:
        write_metadata_block(fh, args, query_cdrs, start_time, cfg)
        fh.write("shard_id,total_seqs,hits,skipped,length_filtered,seq_length_filtered,"
                 "cdr_length_filtered,max_id,min_id,q75_id,best_seq_preview,worst_seq_preview,time_s\n")
    
    # Run full scan
    print("="*80)
    print("RUNNING FULL SCAN")
    print("="*80 + "\n")
    
    total_hits = 0
    total_filtered_length = 0
    total_filtered_seq_length = 0
    total_filtered_cdr_length = 0
    global_best_id = 0.0
    global_worst_id = 1.0
    global_best_seq = ""
    global_worst_seq = ""
    
    outer_bar = tqdm(total=len(npz_files), desc="Processing shards", ncols=100, position=0, leave=True)
    
    for i, npz_path in enumerate(npz_files, start=1):
        stats = process_npz_file_fullscan(npz_path, cfg, query_cdrs, i, len(npz_files))
        
        # Update totals
        total_hits += stats["n_hits"]
        total_filtered_length += stats["n_filtered_length"]
        total_filtered_seq_length += stats["n_filtered_seq_length"]
        total_filtered_cdr_length += stats["n_filtered_cdr_length"]
        
        # Track global best/worst
        if stats["max_id"] > global_best_id:
            global_best_id = stats["max_id"]
            global_best_seq = stats["best_seq"]
        if stats["n_hits"] > 0 and stats["min_id"] < global_worst_id:
            global_worst_id = stats["min_id"]
            global_worst_seq = stats["worst_seq"]
        
        # Append to CSV
        with open(csv_path, "a", encoding="utf-8") as fh:
            fh.write(
                f"{stats['shard']},{stats['n_total']},{stats['n_hits']},{stats['n_skipped']},"
                f"{stats['n_filtered_length']},{stats['n_filtered_seq_length']},"
                f"{stats['n_filtered_cdr_length']},{stats['max_id']:.6f},{stats['min_id']:.6f},"
                f"{stats['q75']:.6f},\"{stats['best_seq']}\",\"{stats['worst_seq']}\",{stats['time_s']:.3f}\n"
            )
        outer_bar.update(1)
    
    outer_bar.close()
    
    # Calculate elapsed time
    end_time = dt.datetime.now()
    elapsed = (end_time - start_time).total_seconds()
    
    # Write summary footer
    with open(csv_path, "a", encoding="utf-8") as fh:
        fh.write("\n# ============================================\n")
        fh.write("# RUN SUMMARY\n")
        fh.write("# ============================================\n")
        fh.write(f"# Run completed: {end_time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        fh.write(f"# Elapsed time: {elapsed:.2f} seconds ({elapsed/60:.2f} minutes)\n")
        fh.write(f"# Total shards processed: {len(npz_files)}\n")
        fh.write(f"# Total hits: {total_hits:,}\n")
        if cfg["use_length_filter"]:
            fh.write(f"# Filtered by CDR length: {total_filtered_length:,}\n")
        if total_filtered_seq_length > 0:
            fh.write(f"# Filtered by sequence length: {total_filtered_seq_length:,}\n")
        if total_filtered_cdr_length > 0:
            fh.write(f"# Filtered by CDR length pre-filter: {total_filtered_cdr_length:,}\n")
        fh.write(f"# Best identity found: {global_best_id:.6f}\n")
        fh.write(f"# Worst identity found: {global_worst_id:.6f}\n")
        fh.write("# ============================================\n")
    
    # Print summary
    print(f"\nâœ… Done â€” all {len(npz_files)} shards processed.")
    print(f"\n{'='*80}")
    print("RUN SUMMARY")
    print(f"{'='*80}")
    print(f"Completed: {end_time.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Elapsed time: {elapsed:.2f} seconds ({elapsed/60:.2f} minutes)")
    print(f"Total hits: {total_hits:,}")
    
    if used_recommendations:
        print(f"\nâœ“ Used recommended thresholds from interactive analysis")
        if total_hits > 0:
            accuracy = abs(total_hits - args.target_hits) / args.target_hits * 100
            print(f"  Target hits: {args.target_hits:,}")
            print(f"  Actual hits: {total_hits:,}")
            print(f"  Accuracy: {100-accuracy:.1f}%")
    
    print(f"\nBest identity: {global_best_id:.6f}")
    print(f"Worst identity: {global_worst_id:.6f}")
    
    print(f"\nOutput saved to:")
    print(f"  {csv_path}")
    print("="*80 + "\n")


if __name__ == "__main__":
    main()
