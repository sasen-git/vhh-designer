#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
NPZ FULLSCAN v8 - VHH Database with Target/Patent Metadata
Enhanced antibody search supporting both OAS and VHH unified databases

NEW IN V8:
- Support for VHH unified database shards (12M+ sequences)
- Displays target/antigen binding information  
- Shows patent IDs and titles for patent-derived sequences
- Source tracking (OAS_Camel, INDI_patent, SAbDab, TheraSAbDab, etc.)
- Shard-aware searching with shard_index.json support
- Annotated-first search option (search high-value sequences first)

DATABASE FORMATS SUPPORTED:
1. OAS format: NPZ with "numberings" array (IMGT-numbered integer arrays)
2. VHH format: NPZ with string arrays (aa_v_full, cdr1, cdr2, cdr3, targets, etc.)

USAGE:
1. Interactive mode:
   python3 npz_fullscan_v8_vhh.py
   
2. VHH database search:
   python3 npz_fullscan_v8_vhh.py --vhh-shards /path/to/VHH_shards --query-seq "EIQLQQ..."
   
3. Search annotated sequences only (fast target lookup):
   python3 npz_fullscan_v8_vhh.py --vhh-shards /path/to/VHH_shards --annotated-only --query-seq "..."
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
import math
from typing import Dict, List, Tuple, Optional
from pathlib import Path
import random

import numpy as np
import pandas as pd
from tqdm import tqdm

# ------------------------
# Configuration Memory
# ------------------------
LAST_QUERY_FILE = os.path.expanduser("~/.ka_search_last_query.txt")
LAST_DB_FILE = os.path.expanduser("~/.ka_search_last_db.txt")

def save_last_query_file(filepath: str):
    """Save the last used query file path."""
    try:
        with open(LAST_QUERY_FILE, 'w') as f:
            f.write(filepath)
    except:
        pass

def load_last_query_file() -> Optional[str]:
    """Load the last used query file path."""
    try:
        if os.path.exists(LAST_QUERY_FILE):
            with open(LAST_QUERY_FILE, 'r') as f:
                filepath = f.read().strip()
                if os.path.exists(filepath):
                    return filepath
    except:
        pass
    return None

def save_last_db(db_path: str):
    """Save the last used database path."""
    try:
        with open(LAST_DB_FILE, 'w') as f:
            f.write(db_path)
    except:
        pass

def load_last_db() -> Optional[str]:
    """Load the last used database path."""
    try:
        if os.path.exists(LAST_DB_FILE):
            with open(LAST_DB_FILE, 'r') as f:
                db_path = f.read().strip()
                if os.path.exists(db_path):
                    return db_path
    except:
        pass
    return None

# ------------------------
# Tag Detection and Removal
# ------------------------
# Common antibody expression tags
KNOWN_TAGS = {
    # N-terminal tags/signals
    'MSKIK': ('N', 'Signal peptide (MSKIK)'),
    'MKWVTFISLLFLFSSAYS': ('N', 'IgK signal peptide'),
    'METDTLLLWVLLLWVPGSTG': ('N', 'Human IgH signal peptide'),
    
    # HA tag (single and 3x)
    'YPYDVPDYA': ('C', 'HA tag'),
    'YPYDVPDYAYPYDVPDYAYPYDVPDYA': ('C', '3xHA tag'),
    
    # Strep tags
    'WSHPQFEK': ('C', 'Strep-tag II'),
    'SAWSHPQFEK': ('C', 'Strep-tag II with SA'),
    
    # C-tag
    'EPEA': ('C', 'C-tag'),
    
    # His tags
    'HHHHHH': ('C', '6xHis tag'),
    'HHHHHHHH': ('C', '8xHis tag'),
    
    # FLAG tag
    'DYKDDDDK': ('C', 'FLAG tag'),
    
    # Myc tag
    'EQKLISEEDL': ('C', 'Myc tag'),
    
    # V5 tag
    'GKPIPNPLLGLDST': ('C', 'V5 tag'),
}

def detect_tags(sequence: str) -> List[Tuple[str, str, str]]:
    """
    Detect known tags in a sequence.
    Returns list of (tag_sequence, position, tag_name).
    """
    detected = []
    seq_upper = sequence.upper()
    
    # Check N-terminal tags (first 30 aa)
    n_term = seq_upper[:30]
    for tag, (pos, name) in KNOWN_TAGS.items():
        if pos == 'N' and n_term.startswith(tag):
            detected.append((tag, 'N-terminal', name))
    
    # Check C-terminal tags (last 50 aa to catch 3xHA)
    c_term = seq_upper[-50:]
    for tag, (pos, name) in KNOWN_TAGS.items():
        if pos == 'C' and tag in c_term:
            # Find position to verify it's at the end
            idx = seq_upper.rfind(tag)
            if idx >= len(seq_upper) - len(tag) - 5:  # Within 5 aa of end
                detected.append((tag, 'N-terminal' if pos == 'N' else 'C-terminal', name))
    
    return detected

def remove_tags(sequence: str, tags_to_remove: List[Tuple[str, str, str]]) -> str:
    """Remove detected tags from sequence."""
    result = sequence.upper()
    
    # Sort tags by position in sequence (remove from end first)
    c_tags = [(t, p, n) for t, p, n in tags_to_remove if 'C-terminal' in p]
    n_tags = [(t, p, n) for t, p, n in tags_to_remove if 'N-terminal' in p]
    
    # Remove C-terminal tags (from end, longest first)
    c_tags.sort(key=lambda x: len(x[0]), reverse=True)
    for tag, pos, name in c_tags:
        idx = result.rfind(tag)
        if idx > len(result) // 2:  # Only if in second half
            result = result[:idx]
    
    # Remove N-terminal tags (from start)
    for tag, pos, name in n_tags:
        if result.startswith(tag):
            result = result[len(tag):]
    
    return result

def interactive_tag_removal(sequence: str) -> str:
    """Interactively detect and offer to remove tags."""
    detected = detect_tags(sequence)
    
    if not detected:
        return sequence
    
    print("\n  ⚠️  Detected expression tags:")
    for tag, pos, name in detected:
        print(f"      • {name} ({pos}): {tag}")
    
    remove = input("\n  Remove these tags? (y/n) [y]: ").strip().lower()
    
    if remove in ['', 'y', 'yes']:
        cleaned = remove_tags(sequence, detected)
        print(f"\n  ✓ Tags removed. New length: {len(cleaned)} aa (was {len(sequence)} aa)")
        return cleaned
    else:
        print("  → Keeping tags in sequence")
        return sequence

# ------------------------
# Sequence Alignment Display  
# ------------------------
def needleman_wunsch_align(seq1: str, seq2: str) -> Tuple[str, str]:
    """
    Perform Needleman-Wunsch global alignment with affine-like gap handling.
    Uses same similarity scoring as identity calculation.
    Returns (aligned_seq1, aligned_seq2) with gaps as '-'.
    """
    n, m = len(seq1), len(seq2)
    
    # Initialize scoring matrix
    # Use same costs: match=0, conservative=0.25, semi=0.5, mismatch=1, gap=1
    
    # dp[i][j] = min cost to align seq1[:i] with seq2[:j]
    dp = [[0.0] * (m + 1) for _ in range(n + 1)]
    
    # Initialize first row and column (gap penalties)
    for i in range(n + 1):
        dp[i][0] = float(i)
    for j in range(m + 1):
        dp[0][j] = float(j)
    
    # Fill the matrix
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            cost = substitution_cost(seq1[i-1], seq2[j-1])
            dp[i][j] = min(
                dp[i-1][j-1] + cost,      # match/mismatch
                dp[i-1][j] + 1.0,          # gap in seq2
                dp[i][j-1] + 1.0           # gap in seq1
            )
    
    # Traceback to get alignment
    aligned1, aligned2 = [], []
    i, j = n, m
    
    while i > 0 or j > 0:
        if i > 0 and j > 0:
            cost = substitution_cost(seq1[i-1], seq2[j-1])
            if dp[i][j] == dp[i-1][j-1] + cost:
                aligned1.append(seq1[i-1])
                aligned2.append(seq2[j-1])
                i -= 1
                j -= 1
                continue
        
        if i > 0 and dp[i][j] == dp[i-1][j] + 1.0:
            aligned1.append(seq1[i-1])
            aligned2.append('-')
            i -= 1
        elif j > 0:
            aligned1.append('-')
            aligned2.append(seq2[j-1])
            j -= 1
        else:
            break
    
    return ''.join(reversed(aligned1)), ''.join(reversed(aligned2))

def highlight_alignment(query: str, target: str) -> Tuple[str, str, str]:
    """
    Create alignment display with proper gap insertion.
    Uses Needleman-Wunsch to find optimal alignment.
    Returns (query_display, match_line, target_display).
    
    Symbols:
      | = exact match (100% credit)
      : = conservative substitution (75% credit) - D↔E, K↔R, L↔I↔V, F↔Y, S↔T, N↔Q
      . = semi-conservative (50% credit) - H↔K/R, M↔L/I/V, W↔F/Y, A↔G/S
      x = non-conservative (0% credit)
        = gap
    """
    # Get proper alignment with gaps
    aligned_q, aligned_t = needleman_wunsch_align(query, target)
    
    m_display = []
    
    for q_char, t_char in zip(aligned_q, aligned_t):
        if q_char == '-' or t_char == '-':
            m_display.append(' ')  # Gap
        elif q_char == t_char:
            m_display.append('|')  # Exact match
        elif t_char in CONSERVATIVE_GROUPS.get(q_char, set()):
            m_display.append(':')  # Conservative (75%)
        elif t_char in SEMICONSERVATIVE_GROUPS.get(q_char, set()):
            m_display.append('.')  # Semi-conservative (50%)
        else:
            m_display.append('x')  # Non-conservative
    
    return aligned_q, ''.join(m_display), aligned_t

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

# VHH Shard paths (common locations)
VHH_SHARD_PATHS = [
    os.path.expanduser("~/KA-Search/VHH_shards"),
    "/home/sasenefrem/KA-Search/VHH_shards",
    os.path.expanduser("~/VHH_shards"),
    "./VHH_shards",
]

# OAS database paths
OAS_DB_PATHS = [
    os.path.expanduser("~/KA-Search/extracted/oas-paper/oasdb_20230111/Heavy"),
    "/home/sasenefrem/KA-Search/extracted/oas-paper/oasdb_20230111/Heavy",
]

# Amino acid similarity groups for position constraints
SIMILAR_GROUPS = {
    'H': ['H', 'K', 'R'], 'K': ['H', 'K', 'R'], 'R': ['H', 'K', 'R'],
    'D': ['D', 'E'], 'E': ['D', 'E'],
    'S': ['S', 'T'], 'T': ['S', 'T'],
    'N': ['N', 'Q'], 'Q': ['N', 'Q'],
    'L': ['L', 'I', 'V'], 'I': ['L', 'I', 'V'], 'V': ['L', 'I', 'V'],
    'F': ['F', 'Y', 'W'], 'Y': ['F', 'Y', 'W'], 'W': ['F', 'Y', 'W'],
    'A': ['A', 'G'], 'G': ['A', 'G'],
    'C': ['C'], 'M': ['M'], 'P': ['P'],
}

# ------------------------
# Timeout Handler
# ------------------------
class _Timeout(Exception):
    pass

def _alarm_handler(signum, frame):
    raise _Timeout()

# ------------------------
# VHH Database Functions
# ------------------------
def load_shard_index(shard_dir: str) -> Optional[Dict]:
    """Load shard_index.json if it exists."""
    index_path = os.path.join(shard_dir, "shard_index.json")
    if os.path.exists(index_path):
        try:
            with open(index_path, 'r') as f:
                return json.load(f)
        except Exception as e:
            print(f"Warning: Could not load shard index: {e}")
    return None

def detect_npz_format(npz_path: str) -> str:
    """Detect if NPZ is OAS format or VHH format."""
    try:
        data = np.load(npz_path, allow_pickle=True)
        keys = list(data.files)
        data.close()
        
        if "numberings" in keys:
            return "oas"
        elif "aa_v_full" in keys or "cdr3" in keys:
            return "vhh"
        elif "arr_0" in keys:
            return "oas"  # Legacy OAS format
        else:
            return "unknown"
    except Exception:
        return "unknown"

def get_vhh_shards(shard_dir: str, annotated_only: bool = False) -> List[str]:
    """Get list of VHH shard files, optionally filtering to annotated-only."""
    shard_index = load_shard_index(shard_dir)
    
    if annotated_only:
        # Look for the annotated shard specifically
        annotated_path = os.path.join(shard_dir, "vhh_annotated.npz")
        if os.path.exists(annotated_path):
            return [annotated_path]
        
        # Fall back to checking shard index
        if shard_index and "shards" in shard_index:
            for shard in shard_index["shards"]:
                if "annotated" in shard.get("filename", "").lower():
                    path = os.path.join(shard_dir, shard["filename"])
                    if os.path.exists(path):
                        return [path]
    
    # Get all NPZ files
    npz_files = sorted(glob.glob(os.path.join(shard_dir, "*.npz")))
    
    # If shard index exists, use its ordering
    if shard_index and "shards" in shard_index:
        ordered_files = []
        for shard in shard_index["shards"]:
            path = os.path.join(shard_dir, shard["filename"])
            if os.path.exists(path):
                ordered_files.append(path)
        if ordered_files:
            return ordered_files
    
    return npz_files

def print_vhh_database_info(shard_dir: str):
    """Print information about VHH database."""
    print("\n" + "="*70)
    print("VHH UNIFIED DATABASE")
    print("="*70)
    
    shard_index = load_shard_index(shard_dir)
    
    if shard_index:
        print(f"\nDatabase: {shard_index.get('database_name', 'VHH Unified')}")
        print(f"Created: {shard_index.get('created', 'Unknown')}")
        print(f"Total sequences: {shard_index.get('total_sequences', 'Unknown'):,}")
        print(f"Number of shards: {len(shard_index.get('shards', []))}")
        
        print("\nShards:")
        print("-"*70)
        for shard in shard_index.get("shards", []):
            desc = shard.get("description", "")
            count = shard.get("count", 0)
            size = shard.get("size_mb", 0)
            print(f"  {shard['filename']}: {count:,} seqs ({size:.1f} MB)")
            if desc:
                print(f"    {desc}")
    else:
        # No index, just count files
        npz_files = glob.glob(os.path.join(shard_dir, "*.npz"))
        print(f"\nFound {len(npz_files)} shard files")
        for f in sorted(npz_files):
            size_mb = os.path.getsize(f) / (1024*1024)
            print(f"  {os.path.basename(f)}: {size_mb:.1f} MB")
    
    print("="*70 + "\n")

# ------------------------
# Decoder functions (OAS format)
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
# CDR Extraction (OAS format)
# ------------------------
def extract_cdr1_from_npz(int_row):
    """Extract CDR-H1 using IMGT positions (27-38)."""
    cdr1 = _ints_to_aa(int_row[26:39])
    return cdr1.replace("-", "")

def extract_cdr2_from_npz(int_row):
    """Extract CDR-H2 using IMGT positions (56-65)."""
    cdr2 = _ints_to_aa(int_row[55:66])
    return cdr2.replace("-", "")

def extract_cdr3_from_npz(int_row):
    """Extract CDR-H3 using IMGT positions (105-117)."""
    cdr3_region = int_row[104:128]
    cdr3_aa = _ints_to_aa(cdr3_region)
    
    if not cdr3_aa:
        return ""
    
    cdr3_aa = cdr3_aa.replace("-", "")
    if not cdr3_aa:
        return ""
    
    if cdr3_aa and cdr3_aa[0] == 'C':
        cdr3_aa = cdr3_aa[1:]
    
    if not cdr3_aa:
        return ""
    
    fr4_motifs = [
        'WGQG', 'WGKG', 'WGRG', 'WGAG', 'WGLG', 'WGEG',
        'WGNG', 'WGSG', 'WGTG', 'WGVG', 'WGMG', 'WGPG',
        'WGDG', 'WGHG', 'WGYG', 'WGFG', 'WGCG', 'WGIG',
        'WGWG', 'WG'
    ]
    
    fr4_start = -1
    for motif in fr4_motifs:
        pos = cdr3_aa.find(motif)
        if pos > 0:
            if fr4_start == -1 or pos < fr4_start:
                fr4_start = pos
    
    if fr4_start > 0:
        cdr3 = cdr3_aa[:fr4_start]
        if 3 <= len(cdr3) <= 25:
            return cdr3
        else:
            return ""
    else:
        return ""

def extract_frameworks_from_npz(int_row):
    """Extract framework regions."""
    fr1 = _ints_to_aa(int_row[0:26])
    fr2 = _ints_to_aa(int_row[39:55])
    fr3 = _ints_to_aa(int_row[66:104])
    fr4 = _ints_to_aa(int_row[118:128])
    frameworks = (fr1 + fr2 + fr3 + fr4).replace("-", "")
    return frameworks

def extract_whole_from_npz(int_row):
    """Extract full sequence from NPZ row."""
    full_seq = _ints_to_aa(int_row[0:128])
    return full_seq.replace("-", "")

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
# Region-by-region identity calculation for full sequences
# ------------------------
def extract_frameworks(full_seq: str, cdr1: str, cdr2: str, cdr3: str) -> Tuple[str, str, str, str]:
    """
    Extract framework regions by finding CDR positions in full sequence.
    Returns (FR1, FR2, FR3, FR4) or empty strings if CDRs not found.
    """
    if not full_seq:
        return "", "", "", ""
    
    # Check for empty CDRs (can't find position of empty string meaningfully)
    if not cdr1 or not cdr2 or not cdr3:
        return "", "", "", ""
    
    # Find CDR positions
    cdr1_start = full_seq.find(cdr1)
    cdr2_start = full_seq.find(cdr2)
    cdr3_start = full_seq.find(cdr3)
    
    # If we can't find CDRs, return empty
    if cdr1_start < 0 or cdr2_start < 0 or cdr3_start < 0:
        return "", "", "", ""
    
    cdr1_end = cdr1_start + len(cdr1)
    cdr2_end = cdr2_start + len(cdr2)
    cdr3_end = cdr3_start + len(cdr3)
    
    fr1 = full_seq[:cdr1_start]
    fr2 = full_seq[cdr1_end:cdr2_start]
    fr3 = full_seq[cdr2_end:cdr3_start]
    fr4 = full_seq[cdr3_end:]
    
    return fr1, fr2, fr3, fr4


def calc_region_by_region_identity(
    q_full: str, q_cdr1: str, q_cdr2: str, q_cdr3: str,
    db_full: str, db_cdr1: str, db_cdr2: str, db_cdr3: str
) -> Tuple[float, float]:
    """
    Calculate full sequence identity by comparing each region separately.
    This is ~5x faster than full edit distance and handles length differences properly.
    
    Returns (exact_identity, weighted_identity).
    """
    if not q_full or not db_full:
        return 0.0, 0.0
    
    # Extract frameworks
    q_fr1, q_fr2, q_fr3, q_fr4 = extract_frameworks(q_full, q_cdr1, q_cdr2, q_cdr3)
    db_fr1, db_fr2, db_fr3, db_fr4 = extract_frameworks(db_full, db_cdr1, db_cdr2, db_cdr3)
    
    # Determine which comparison mode to use
    have_q_frameworks = bool(q_fr1)
    have_db_frameworks = bool(db_fr1)
    
    if have_q_frameworks and have_db_frameworks:
        # Best case: Compare all 7 regions separately
        regions = [
            (q_fr1, db_fr1),
            (q_cdr1, db_cdr1),
            (q_fr2, db_fr2),
            (q_cdr2, db_cdr2),
            (q_fr3, db_fr3),
            (q_cdr3, db_cdr3),
            (q_fr4, db_fr4),
        ]
    elif q_cdr1 and q_cdr2 and q_cdr3 and db_cdr1 and db_cdr2 and db_cdr3:
        # Fallback: CDRs available but frameworks couldn't be extracted
        regions = [
            (q_cdr1, db_cdr1),
            (q_cdr2, db_cdr2),
            (q_cdr3, db_cdr3),
        ]
    else:
        # Last resort: Fall back to full sequence comparison
        # This is slower but handles missing CDR annotations
        max_len = max(len(q_full), len(db_full))
        weighted_dist = edit_dist_weighted(q_full, db_full)
        exact_dist = edit_dist_exact(q_full, db_full)
        return (1.0 - exact_dist / max_len, 1.0 - weighted_dist / max_len)
    
    # Calculate total edit distance across all regions
    total_weighted_dist = 0.0
    total_exact_dist = 0
    total_max_len = 0
    
    for q_region, db_region in regions:
        if q_region and db_region:
            max_len = max(len(q_region), len(db_region))
            total_max_len += max_len
            total_weighted_dist += edit_dist_weighted(q_region, db_region)
            total_exact_dist += edit_dist_exact(q_region, db_region)
        elif q_region:
            # Only query has this region - counts as deletions
            total_max_len += len(q_region)
            total_weighted_dist += len(q_region)
            total_exact_dist += len(q_region)
        elif db_region:
            # Only db has this region - counts as insertions
            total_max_len += len(db_region)
            total_weighted_dist += len(db_region)
            total_exact_dist += len(db_region)
    
    if total_max_len == 0:
        return 0.0, 0.0
    
    exact_id = 1.0 - (total_exact_dist / total_max_len)
    weighted_id = 1.0 - (total_weighted_dist / total_max_len)
    
    return exact_id, weighted_id


# ------------------------
# CDR Weighting for multi-region searches
# ------------------------
# CDR3 is most important for antigen binding specificity
CDR_WEIGHTS = {
    "cdr1": 0.2,   # 20% - contributes to binding but less variable
    "cdr2": 0.3,   # 30% - important for binding
    "cdr3": 0.5,   # 50% - most critical for specificity
    "full": 1.0,   # 100% if searching full sequence only
}

def calc_weighted_cdr_score(identities: Dict[str, float], regions: List[str]) -> float:
    """
    Calculate weighted CDR identity score.
    
    Args:
        identities: Dict of region -> identity (e.g., {"cdr1": 0.7, "cdr2": 0.8, "cdr3": 0.6})
        regions: List of regions being searched
    
    Returns:
        Weighted average identity
    """
    if not identities or not regions:
        return 0.0
    
    # Get weights for searched regions
    total_weight = 0.0
    weighted_sum = 0.0
    
    for region in regions:
        region = region.strip().lower()
        if region in identities and identities[region] > 0:
            weight = CDR_WEIGHTS.get(region, 0.33)  # Default equal weight
            weighted_sum += weight * identities[region]
            total_weight += weight
    
    if total_weight == 0:
        return 0.0
    
    return weighted_sum / total_weight


# ------------------------
# Identity calculation with similar AA scoring
# ------------------------
# Amino acid similarity groups based on biochemical properties
# Conservative: very similar properties (0.25 cost)
# Semi-conservative: somewhat similar (0.5 cost)

CONSERVATIVE_GROUPS = {
    # Strongly similar - 0.25 cost (75% credit)
    'D': {'E'},                # Acidic (negative charge)
    'E': {'D'},
    'K': {'R'},                # Basic (positive charge, similar size)
    'R': {'K'},
    'L': {'I', 'V'},           # Aliphatic hydrophobic (branched)
    'I': {'L', 'V'},
    'V': {'L', 'I'},
    'F': {'Y'},                # Aromatic (similar size)
    'Y': {'F'},
    'S': {'T'},                # Small polar (hydroxyl)
    'T': {'S'},
    'N': {'Q'},                # Amide
    'Q': {'N'},
}

SEMICONSERVATIVE_GROUPS = {
    # Somewhat similar - 0.5 cost (50% credit)
    'H': {'K', 'R'},           # Basic (different structure)
    'K': {'H'},
    'R': {'H'},
    'M': {'L', 'I', 'V'},      # Hydrophobic (sulfur-containing)
    'L': {'M'},
    'I': {'M'},
    'V': {'M'},
    'W': {'F', 'Y'},           # Aromatic (larger)
    'F': {'W'},
    'Y': {'W'},
    'A': {'G', 'S'},           # Small
    'G': {'A'},
    'S': {'A'},
}

def substitution_cost(c1: str, c2: str) -> float:
    """
    Return substitution cost between two amino acids.
    0.00 = exact match (100% credit)
    0.25 = conservative substitution (75% credit)
    0.50 = semi-conservative substitution (50% credit)
    1.00 = non-conservative (0% credit)
    """
    if c1 == c2:
        return 0.0
    if c2 in CONSERVATIVE_GROUPS.get(c1, set()):
        return 0.25
    if c2 in SEMICONSERVATIVE_GROUPS.get(c1, set()):
        return 0.5
    return 1.0

def edit_dist_weighted(s1: str, s2: str) -> float:
    """
    Calculate weighted Levenshtein distance with partial credit for similar AAs.
    Returns float distance (can be non-integer due to fractional costs).
    """
    if len(s1) < len(s2):
        return edit_dist_weighted(s2, s1)
    
    if len(s2) == 0:
        return float(len(s1))
    
    prev_row = [float(i) for i in range(len(s2) + 1)]
    
    for i, c1 in enumerate(s1):
        curr_row = [float(i + 1)]
        for j, c2 in enumerate(s2):
            insertions = prev_row[j + 1] + 1.0
            deletions = curr_row[j] + 1.0
            substitutions = prev_row[j] + substitution_cost(c1, c2)
            curr_row.append(min(insertions, deletions, substitutions))
        prev_row = curr_row
    
    return prev_row[-1]

def edit_dist_exact(s1: str, s2: str) -> int:
    """Calculate exact Levenshtein distance (no partial credit)."""
    if len(s1) < len(s2):
        return edit_dist_exact(s2, s1)
    
    if len(s2) == 0:
        return len(s1)
    
    prev_row = list(range(len(s2) + 1))
    
    for i, c1 in enumerate(s1):
        curr_row = [i + 1]
        for j, c2 in enumerate(s2):
            insertions = prev_row[j + 1] + 1
            deletions = curr_row[j] + 1
            substitutions = prev_row[j] + (c1 != c2)
            curr_row.append(min(insertions, deletions, substitutions))
        prev_row = curr_row
    
    return prev_row[-1]

def calc_identity_both(qval: str, dval: str) -> Tuple[float, float]:
    """
    Calculate both exact and weighted identity.
    Returns (exact_identity, weighted_identity).
    """
    if not qval or not dval:
        return 0.0, 0.0
    max_len = max(len(qval), len(dval))
    if max_len == 0:
        return 0.0, 0.0
    
    exact_dist = edit_dist_exact(qval, dval)
    weighted_dist = edit_dist_weighted(qval, dval)
    
    exact_id = 1.0 - (exact_dist / max_len)
    weighted_id = 1.0 - (weighted_dist / max_len)
    
    return exact_id, weighted_id

def calc_id_by_region(qval: str, dval: str) -> float:
    """Calculate weighted identity for a region (for backward compatibility)."""
    _, weighted_id = calc_identity_both(qval, dval)
    return weighted_id

def calc_exact_id_by_region(qval: str, dval: str) -> float:
    """Calculate exact identity for a region."""
    exact_id, _ = calc_identity_both(qval, dval)
    return exact_id

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
            "cdr1": grab(27, 38),
            "cdr2": grab(56, 65),
            "cdr3": grab(105, 117)
        }
    
    except _Timeout:
        if hasattr(signal, 'SIGALRM'):
            signal.alarm(0)
        return {"cdr1": "", "cdr2": "", "cdr3": heuristic_cdrh3(seq_clean)}
    except Exception:
        if hasattr(signal, 'SIGALRM'):
            signal.alarm(0)
        return {"cdr1": "", "cdr2": "", "cdr3": heuristic_cdrh3(seq_clean)}

# ------------------------
# VHH NPZ Processing
# ------------------------
def process_vhh_npz_file(npz_path: str, config: Dict, query_cdrs: Dict,
                         shard_idx: int, total_shards: int) -> Dict:
    """Process a VHH format NPZ file."""
    shard_name = os.path.basename(npz_path)
    start_time = time.time()
    
    # Print shard header
    print(f"\n  [{shard_idx}/{total_shards}] {shard_name}")
    
    try:
        data = np.load(npz_path, allow_pickle=True)
        
        # Load arrays
        ids = data.get("ids", data.get("arr_0", np.array([])))
        sources = data.get("source", np.array([]))
        aa_full = data.get("aa_v_full", np.array([]))
        cdr1_arr = data.get("cdr1", np.array([]))
        cdr2_arr = data.get("cdr2", np.array([]))
        cdr3_arr = data.get("cdr3", np.array([]))
        
        # Metadata arrays
        targets_arr = data.get("targets", np.array([]))
        patent_id_arr = data.get("patent_id", np.array([]))
        patent_title_arr = data.get("patent_title", np.array([]))
        organism_arr = data.get("organism", np.array([]))
        pdb_id_arr = data.get("pdb_id", np.array([]))
        reference_arr = data.get("reference", np.array([]))
        
        n_seqs = len(aa_full) if len(aa_full) > 0 else len(ids)
        
        if n_seqs == 0:
            data.close()
            return empty_stats(shard_name)
        
        regions = config["regions"].split(",")
        regions = [r.strip() for r in regions]
        
        # Determine what we need to extract
        need_full = "full" in regions
        need_cdr1 = "cdr1" in regions or need_full  # Also need CDRs for region-by-region identity
        need_cdr2 = "cdr2" in regions or need_full
        need_cdr3 = "cdr3" in regions or need_full
        
        # Initialize counters
        n_hits = 0
        n_skipped = 0
        n_filtered_length = 0
        n_filtered_seq_length = 0
        n_filtered_cdr_length = 0
        
        identities = []
        best_id = 0.0
        worst_id = 1.0
        best_seq = ""
        worst_seq = ""
        hits_data = []
        
        # Helper to safely get string from array
        def safe_str(arr, idx, default=""):
            if arr is None or len(arr) == 0 or idx >= len(arr):
                return default
            val = arr[idx]
            if val is None:
                return default
            if isinstance(val, bytes):
                return val.decode('utf-8', errors='ignore')
            return str(val) if val else default
        
        # Process sequences with progress bar
        print(f"      Sequences: {n_seqs:,}")
        inner_bar = tqdm(range(n_seqs), 
                        desc=f"      Progress",
                        ncols=100, position=0, leave=False,
                        mininterval=1.0)
        
        last_update = time.time()
        
        for i in inner_bar:
            # Only extract what we need for the search
            # For full sequence search, CDRs are optional (for region-by-region)
            # For CDR-only search, missing CDRs mean we skip
            
            if need_cdr3:
                db_cdr3 = safe_str(cdr3_arr, i)
                if not db_cdr3 and "cdr3" in regions:
                    # Only skip if CDR3 is explicitly being searched
                    n_skipped += 1
                    continue
            else:
                db_cdr3 = ""
            
            if need_cdr1:
                db_cdr1 = safe_str(cdr1_arr, i)
                if not db_cdr1 and "cdr1" in regions:
                    n_skipped += 1
                    continue
            else:
                db_cdr1 = ""
                
            if need_cdr2:
                db_cdr2 = safe_str(cdr2_arr, i)
                if not db_cdr2 and "cdr2" in regions:
                    n_skipped += 1
                    continue
            else:
                db_cdr2 = ""
            
            if need_full:
                full_seq = safe_str(aa_full, i)
                if not full_seq:
                    n_skipped += 1
                    continue
            else:
                full_seq = ""
            
            # CDR length pre-filter (fast rejection)
            if config.get("enable_cdr_length_prefilter"):
                skip_cdr = False
                if need_cdr3:
                    q_len = config["query_lengths"].get("cdr3", 0)
                    if q_len > 0 and abs(len(db_cdr3) - q_len) > 3:
                        skip_cdr = True
                if not skip_cdr and need_cdr1:
                    q_len = config["query_lengths"].get("cdr1", 0)
                    if q_len > 0 and abs(len(db_cdr1) - q_len) > 3:
                        skip_cdr = True
                if not skip_cdr and need_cdr2:
                    q_len = config["query_lengths"].get("cdr2", 0)
                    if q_len > 0 and abs(len(db_cdr2) - q_len) > 3:
                        skip_cdr = True
                
                if skip_cdr:
                    n_filtered_cdr_length += 1
                    continue
            
            # Calculate identities ONLY for searched regions
            region_identities = {}  # Store as dict for weighted calculation
            passes_threshold = True
            
            for region in regions:
                if region == "cdr1":
                    db_seq, q_seq = db_cdr1, query_cdrs.get("cdr1", "")
                elif region == "cdr2":
                    db_seq, q_seq = db_cdr2, query_cdrs.get("cdr2", "")
                elif region == "cdr3":
                    db_seq, q_seq = db_cdr3, query_cdrs.get("cdr3", "")
                elif region == "full":
                    db_seq, q_seq = full_seq, query_cdrs.get("full", "")
                else:
                    continue
                
                if not db_seq or not q_seq:
                    passes_threshold = False
                    break
                
                # Length filter
                if config.get("use_length_filter"):
                    window = config["length_windows"].get(region, None)
                    if window is not None:
                        if abs(len(db_seq) - len(q_seq)) > window:
                            n_filtered_length += 1
                            passes_threshold = False
                            break
                
                # Calculate identity
                if region == "full":
                    # Use region-by-region for accurate full sequence comparison
                    _, identity = calc_region_by_region_identity(
                        q_seq, query_cdrs.get("cdr1", ""), query_cdrs.get("cdr2", ""), query_cdrs.get("cdr3", ""),
                        db_seq, db_cdr1, db_cdr2, db_cdr3
                    )
                else:
                    identity = calc_id_by_region(q_seq, db_seq)
                region_identities[region] = identity
                
                # Check per-region threshold (if set)
                min_id = config.get(f"min_id_{region}")
                if min_id is not None and identity < min_id:
                    passes_threshold = False
                    break
            
            if passes_threshold and region_identities:
                # Calculate weighted score for multi-CDR searches
                if len(regions) > 1 and "full" not in regions:
                    weighted_identity = calc_weighted_cdr_score(region_identities, regions)
                else:
                    weighted_identity = np.mean(list(region_identities.values()))
                
                # Check weighted threshold (if set)
                min_weighted = config.get("min_weighted_identity")
                if min_weighted is not None and weighted_identity < min_weighted:
                    continue  # Skip this hit
                
                n_hits += 1
                identities.append(weighted_identity)
                
                # Get metadata and full sequence for output (lazy load)
                if not full_seq:
                    full_seq = safe_str(aa_full, i)
                if not db_cdr1:
                    db_cdr1 = safe_str(cdr1_arr, i)
                if not db_cdr2:
                    db_cdr2 = safe_str(cdr2_arr, i)
                if not db_cdr3:
                    db_cdr3 = safe_str(cdr3_arr, i)
                
                seq_id = safe_str(ids, i, f"seq_{i}")
                source = safe_str(sources, i)
                target = safe_str(targets_arr, i)
                patent_id = safe_str(patent_id_arr, i)
                patent_title = safe_str(patent_title_arr, i)
                organism = safe_str(organism_arr, i)
                pdb_id = safe_str(pdb_id_arr, i)
                reference = safe_str(reference_arr, i)
                
                # Calculate identities for all CDRs for output
                cdr1_id = calc_id_by_region(query_cdrs.get("cdr1", ""), db_cdr1) if db_cdr1 and query_cdrs.get("cdr1") else 0.0
                cdr2_id = calc_id_by_region(query_cdrs.get("cdr2", ""), db_cdr2) if db_cdr2 and query_cdrs.get("cdr2") else 0.0
                cdr3_id = calc_id_by_region(query_cdrs.get("cdr3", ""), db_cdr3) if db_cdr3 and query_cdrs.get("cdr3") else 0.0
                
                hit_data = {
                    "shard": shard_name,
                    "index": i,
                    "id": seq_id,
                    "source": source,
                    "weighted_identity": weighted_identity,
                    "avg_identity": weighted_identity,  # Keep for backward compat
                    "full_sequence": full_seq,
                    "full_identity": region_identities.get("full", 0.0),  # Actual full identity
                    "cdr1_seq": db_cdr1,
                    "cdr1_identity": region_identities.get("cdr1", cdr1_id),  # Use searched value if available
                    "cdr2_seq": db_cdr2,
                    "cdr2_identity": region_identities.get("cdr2", cdr2_id),
                    "cdr3_seq": db_cdr3,
                    "cdr3_identity": region_identities.get("cdr3", cdr3_id),
                    "target": target,
                    "patent_id": patent_id,
                    "patent_title": patent_title,
                    "organism": organism,
                    "pdb_id": pdb_id,
                    "reference": reference,
                }
                
                hits_data.append(hit_data)
                
                if weighted_identity > best_id:
                    best_id = weighted_identity
                    best_seq = db_cdr3 if db_cdr3 else full_seq[:50]
                
                if weighted_identity < worst_id:
                    worst_id = weighted_identity
                    worst_seq = db_cdr3 if db_cdr3 else full_seq[:50]
            
            # Update progress bar with stats
            if (i + 1) % 50000 == 0 or time.time() - last_update > 5:
                inner_bar.set_postfix({
                    "hits": n_hits,
                    "skip": n_skipped,
                    "filt": n_filtered_length + n_filtered_cdr_length
                })
                last_update = time.time()
        
        inner_bar.close()
        data.close()
        
        elapsed = time.time() - start_time
        rate = n_seqs / elapsed if elapsed > 0 else 0
        
        # Print shard summary
        print(f"      ✓ Done: {n_hits:,} hits, {n_skipped:,} skipped, {n_filtered_length + n_filtered_cdr_length:,} filtered")
        print(f"      ✓ Time: {elapsed:.1f}s ({rate:,.0f} seq/s)")
        
        q75 = np.percentile(identities, 75) if identities else 0.0
        
        return {
            "shard": shard_name,
            "n_total": n_seqs,
            "n_hits": n_hits,
            "n_skipped": n_skipped,
            "n_filtered_length": n_filtered_length,
            "n_filtered_seq_length": n_filtered_seq_length,
            "n_filtered_cdr_length": n_filtered_cdr_length,
            "max_id": best_id,
            "min_id": worst_id if n_hits > 0 else 0.0,
            "q75": q75,
            "best_seq": best_seq,
            "worst_seq": worst_seq,
            "time_s": elapsed,
            "hits": hits_data
        }
    
    except Exception as e:
        print(f"\n      ❌ Error: {e}")
        import traceback
        traceback.print_exc()
        return empty_stats(shard_name)

def empty_stats(shard_name: str) -> Dict:
    """Return empty stats dictionary."""
    return {
        "shard": shard_name,
        "n_total": 0,
        "n_hits": 0,
        "n_skipped": 0,
        "n_filtered_length": 0,
        "n_filtered_seq_length": 0,
        "n_filtered_cdr_length": 0,
        "max_id": 0.0,
        "min_id": 0.0,
        "q75": 0.0,
        "best_seq": "",
        "worst_seq": "",
        "time_s": 0.0,
        "hits": []
    }

# ------------------------
# OAS NPZ Processing (from v7)
# ------------------------
def process_oas_npz_file(npz_path: str, config: Dict, query_cdrs: Dict,
                         shard_idx: int, total_shards: int) -> Dict:
    """Process an OAS format NPZ file."""
    shard_name = os.path.basename(npz_path)
    start_time = time.time()
    
    try:
        data = np.load(npz_path, allow_pickle=True)
        
        if "numberings" not in data:
            if "arr_0" in data.files:
                arr = data["arr_0"]
            elif "sequences" in data.files:
                arr = data["sequences"]
            else:
                data.close()
                return empty_stats(shard_name)
        else:
            arr = data["numberings"]
        
        n_seqs = len(arr)
        regions = config["regions"].split(",")
        
        n_hits = 0
        n_skipped = 0
        n_filtered_length = 0
        n_filtered_seq_length = 0
        n_filtered_cdr_length = 0
        
        identities = []
        best_id = 0.0
        worst_id = 1.0
        best_seq = ""
        worst_seq = ""
        hits_data = []
        
        inner_bar = tqdm(range(n_seqs), 
                        desc=f"  Shard {shard_idx}/{total_shards}",
                        ncols=100, position=1, leave=False)
        
        for i in inner_bar:
            row = arr[i]
            
            # Sequence length pre-filter
            if config.get("enable_seq_length_prefilter") and "full" in regions:
                db_full = extract_whole_from_npz(row)
                if not db_full:
                    n_skipped += 1
                    continue
                
                len_ratio = abs(len(db_full) - len(query_cdrs["full"])) / len(query_cdrs["full"])
                if len_ratio > config.get("seq_length_threshold", 0.25):
                    n_filtered_seq_length += 1
                    continue
            
            # CDR length pre-filter
            if config.get("enable_cdr_length_prefilter"):
                skip_cdr = False
                for cdr in ["cdr1", "cdr2", "cdr3"]:
                    if cdr not in regions:
                        continue
                    
                    if cdr == "cdr1":
                        db_cdr = extract_cdr1_from_npz(row)
                    elif cdr == "cdr2":
                        db_cdr = extract_cdr2_from_npz(row)
                    else:
                        db_cdr = extract_cdr3_from_npz(row)
                    
                    q_cdr_len = config["query_lengths"].get(cdr, 0)
                    if q_cdr_len > 0:
                        len_diff = abs(len(db_cdr) - q_cdr_len)
                        if len_diff > 3:
                            skip_cdr = True
                            break
                
                if skip_cdr:
                    n_filtered_cdr_length += 1
                    continue
            
            # Extract sequences
            full_seq = extract_whole_from_npz(row)
            db_cdr1 = extract_cdr1_from_npz(row)
            db_cdr2 = extract_cdr2_from_npz(row)
            db_cdr3 = extract_cdr3_from_npz(row)
            
            all_cdrs = {'cdr1': db_cdr1, 'cdr2': db_cdr2, 'cdr3': db_cdr3}
            all_identities = {}
            
            for cdr_name, db_cdr in all_cdrs.items():
                q_cdr = query_cdrs.get(cdr_name, "")
                if db_cdr and q_cdr:
                    all_identities[cdr_name] = calc_id_by_region(q_cdr, db_cdr)
                else:
                    all_identities[cdr_name] = 0.0
            
            if full_seq and query_cdrs.get('full'):
                # Use region-by-region for accurate full sequence comparison
                _, all_identities['full'] = calc_region_by_region_identity(
                    query_cdrs['full'], query_cdrs.get("cdr1", ""), query_cdrs.get("cdr2", ""), query_cdrs.get("cdr3", ""),
                    full_seq, db_cdr1, db_cdr2, db_cdr3
                )
            else:
                all_identities['full'] = 0.0
            
            # Check thresholds
            region_identities = []
            passes_threshold = True
            region_sequences = {}
            
            for region in regions:
                region = region.strip()
                
                if region == "cdr1":
                    db_seq = db_cdr1
                    q_seq = query_cdrs.get("cdr1", "")
                elif region == "cdr2":
                    db_seq = db_cdr2
                    q_seq = query_cdrs.get("cdr2", "")
                elif region == "cdr3":
                    db_seq = db_cdr3
                    q_seq = query_cdrs.get("cdr3", "")
                elif region == "frameworks":
                    db_seq = extract_frameworks_from_npz(row)
                    q_seq = query_cdrs.get("frameworks", "")
                    region_sequences['frameworks'] = db_seq
                    if db_seq and q_seq:
                        all_identities['frameworks'] = calc_id_by_region(q_seq, db_seq)
                elif region == "full":
                    db_seq = full_seq
                    q_seq = query_cdrs.get("full", "")
                else:
                    continue
                
                if not db_seq or not q_seq:
                    passes_threshold = False
                    break
                
                if config.get("use_length_filter"):
                    window = config["length_windows"].get(region, None)
                    if window is not None:
                        if abs(len(db_seq) - len(q_seq)) > window:
                            n_filtered_length += 1
                            passes_threshold = False
                            break
                
                identity = all_identities.get(region, 0.0)
                region_identities.append(identity)
                
                min_id = config.get(f"min_id_{region}")
                if min_id is not None and identity < min_id:
                    passes_threshold = False
                    break
            
            if passes_threshold and region_identities:
                # Calculate weighted score for multi-CDR searches
                region_ids = {regions[j]: region_identities[j] for j in range(len(region_identities))}
                if len(regions) > 1 and "full" not in regions:
                    weighted_identity = calc_weighted_cdr_score(region_ids, regions)
                else:
                    weighted_identity = np.mean(region_identities)
                
                # Check weighted threshold (if set)
                min_weighted = config.get("min_weighted_identity")
                if min_weighted is not None and weighted_identity < min_weighted:
                    continue  # Skip this hit
                
                n_hits += 1
                identities.append(weighted_identity)
                
                hit_data = {
                    "shard": shard_name,
                    "index": i,
                    "id": f"oas_{i}",
                    "source": "OAS",
                    "weighted_identity": weighted_identity,
                    "avg_identity": weighted_identity,  # Keep for backward compat
                    "full_sequence": full_seq,
                    "full_identity": all_identities.get('full', 0.0),
                    "cdr1_seq": db_cdr1,
                    "cdr1_identity": all_identities['cdr1'],
                    "cdr2_seq": db_cdr2,
                    "cdr2_identity": all_identities['cdr2'],
                    "cdr3_seq": db_cdr3,
                    "cdr3_identity": all_identities['cdr3'],
                    # No metadata in OAS format
                    "target": "",
                    "patent_id": "",
                    "patent_title": "",
                    "organism": "",
                    "pdb_id": "",
                    "reference": "",
                }
                
                if 'frameworks' in region_sequences:
                    hit_data["frameworks_seq"] = region_sequences['frameworks']
                    hit_data["frameworks_identity"] = all_identities.get('frameworks', 0.0)
                
                hits_data.append(hit_data)
                
                if weighted_identity > best_id:
                    best_id = weighted_identity
                    best_seq = full_seq[:50] if "full" in regions else f"CDR3: {db_cdr3}"
                
                if weighted_identity < worst_id:
                    worst_id = weighted_identity
                    worst_seq = full_seq[:50] if "full" in regions else f"CDR3: {db_cdr3}"
            
            if (i + 1) % 1000 == 0:
                inner_bar.set_postfix({"hits": n_hits})
        
        inner_bar.close()
        data.close()
        
        q75 = np.percentile(identities, 75) if identities else 0.0
        
        return {
            "shard": shard_name,
            "n_total": n_seqs,
            "n_hits": n_hits,
            "n_skipped": n_skipped,
            "n_filtered_length": n_filtered_length,
            "n_filtered_seq_length": n_filtered_seq_length,
            "n_filtered_cdr_length": n_filtered_cdr_length,
            "max_id": best_id,
            "min_id": worst_id if n_hits > 0 else 0.0,
            "q75": q75,
            "best_seq": best_seq,
            "worst_seq": worst_seq,
            "time_s": time.time() - start_time,
            "hits": hits_data
        }
    
    except Exception as e:
        print(f"\nError processing {shard_name}: {e}")
        return empty_stats(shard_name)

# ------------------------
# Generic NPZ Processing
# ------------------------
def process_npz_file(npz_path: str, config: Dict, query_cdrs: Dict,
                     shard_idx: int, total_shards: int) -> Dict:
    """Process NPZ file, auto-detecting format."""
    fmt = detect_npz_format(npz_path)
    
    if fmt == "vhh":
        return process_vhh_npz_file(npz_path, config, query_cdrs, shard_idx, total_shards)
    else:
        return process_oas_npz_file(npz_path, config, query_cdrs, shard_idx, total_shards)

# ------------------------
# Excel Export with Metadata
# ------------------------
def get_sort_column_for_regions(regions: List[str]) -> str:
    """Determine the correct column to sort by based on search regions."""
    regions = [r.strip().lower() for r in regions]
    
    if len(regions) == 1:
        # Single region search - sort by that region's identity
        region = regions[0]
        if region == "cdr1":
            return "cdr1_identity"
        elif region == "cdr2":
            return "cdr2_identity"
        elif region == "cdr3":
            return "cdr3_identity"
        elif region == "full":
            return "full_identity"
    elif "full" in regions:
        # If full is included, sort by full
        return "full_identity"
    else:
        # Multi-CDR search - sort by average/weighted identity
        return "avg_identity"
    
    return "avg_identity"  # Default fallback


def export_to_excel(csv_path: str, results_csv_path: str, output_folder: str, 
                   query_cdrs: Dict = None, include_metadata: bool = True,
                   regions: List[str] = None) -> Optional[str]:
    """Export CSV results to formatted Excel file with metadata columns."""
    try:
        base_name = os.path.basename(csv_path).replace('.csv', '')
        excel_path = os.path.join(output_folder, base_name + '.xlsx')
        
        # Determine sort column based on search regions
        sort_col = get_sort_column_for_regions(regions) if regions else "avg_identity"
        
        # Read the summary CSV
        with open(csv_path, 'r') as f:
            lines = f.readlines()
            data_start = 0
            for i, line in enumerate(lines):
                if line.startswith('shard_id,') or line.startswith('shard,'):
                    data_start = i
                    break
        
        summary_df = pd.read_csv(csv_path, skiprows=data_start)
        
        # Read results data
        results_df = None
        if os.path.exists(results_csv_path):
            with open(results_csv_path, 'r') as f:
                lines = f.readlines()
                data_start = 0
                for i, line in enumerate(lines):
                    if not line.startswith('#'):
                        data_start = i
                        break
            results_df = pd.read_csv(results_csv_path, skiprows=data_start)
            
            # Sort by the appropriate identity column (descending - highest first)
            if sort_col in results_df.columns:
                results_df = results_df.sort_values(by=sort_col, ascending=False)
                print(f"[INFO] Results sorted by {sort_col} (descending)", file=sys.stderr)
            elif 'avg_identity' in results_df.columns:
                results_df = results_df.sort_values(by='avg_identity', ascending=False)
                print(f"[INFO] Results sorted by avg_identity (descending)", file=sys.stderr)
        
        # Create Excel writer
        with pd.ExcelWriter(excel_path, engine='openpyxl') as writer:
            # Write summary sheet
            summary_df.to_excel(writer, sheet_name='Summary', index=False)
            
            # Write results sheet
            if results_df is not None and not results_df.empty:
                # Reorder columns for better readability
                priority_cols = ['id', 'source', 'avg_identity', 'target', 'patent_id',
                                'cdr3_identity', 'cdr3_seq', 'cdr1_identity', 'cdr1_seq',
                                'cdr2_identity', 'cdr2_seq', 'full_identity', 'full_sequence']
                
                # Get columns that exist
                ordered_cols = [c for c in priority_cols if c in results_df.columns]
                other_cols = [c for c in results_df.columns if c not in ordered_cols]
                results_df = results_df[ordered_cols + other_cols]
                
                max_rows = 1048576 - 2
                if len(results_df) > max_rows:
                    for i, start_idx in enumerate(range(0, len(results_df), max_rows)):
                        end_idx = min(start_idx + max_rows, len(results_df))
                        sheet_name = f'Results_{i+1}'
                        results_df.iloc[start_idx:end_idx].to_excel(
                            writer, sheet_name=sheet_name, index=False
                        )
                else:
                    results_df.to_excel(writer, sheet_name='Results', index=False)
                
                # Create Annotated-only sheet (sequences with target info)
                if 'target' in results_df.columns:
                    annotated_df = results_df[results_df['target'].notna() & (results_df['target'] != '')]
                    if not annotated_df.empty:
                        annotated_df.to_excel(writer, sheet_name='Annotated Hits', index=False)
            
            # Add metadata sheet
            metadata = {
                'File': [os.path.basename(csv_path)],
                'Created': [dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S')],
                'Total Hits': [len(results_df) if results_df is not None else 0],
                'Annotated Hits': [len(results_df[results_df['target'].notna() & (results_df['target'] != '')]) 
                                  if results_df is not None and 'target' in results_df.columns else 0],
                'Summary Rows': [len(summary_df)],
                'Sorted By': [sort_col]
            }
            metadata_df = pd.DataFrame(metadata)
            metadata_df.to_excel(writer, sheet_name='Info', index=False)
            
            # Add Query Information sheet
            if query_cdrs:
                query_info_data = [
                    ['QUERY INFORMATION', ''],
                    ['', ''],
                    ['Full Query Sequence:', query_cdrs.get('full', 'N/A')],
                    ['Query Length:', f"{len(query_cdrs.get('full', ''))} aa" if query_cdrs.get('full') else 'N/A'],
                    ['', ''],
                    ['EXTRACTED CDRs:', ''],
                    ['CDR1 Sequence:', query_cdrs.get('cdr1', 'N/A')],
                    ['CDR1 Length:', f"{len(query_cdrs.get('cdr1', ''))} aa" if query_cdrs.get('cdr1') else 'N/A'],
                    ['CDR2 Sequence:', query_cdrs.get('cdr2', 'N/A')],
                    ['CDR2 Length:', f"{len(query_cdrs.get('cdr2', ''))} aa" if query_cdrs.get('cdr2') else 'N/A'],
                    ['CDR3 Sequence:', query_cdrs.get('cdr3', 'N/A')],
                    ['CDR3 Length:', f"{len(query_cdrs.get('cdr3', ''))} aa" if query_cdrs.get('cdr3') else 'N/A'],
                ]
                query_df = pd.DataFrame(query_info_data, columns=['Parameter', 'Value'])
                query_df.to_excel(writer, sheet_name='Query Info', index=False)
            
            # Format sheets
            for sheet_name in writer.sheets:
                worksheet = writer.sheets[sheet_name]
                for column_cells in worksheet.columns:
                    length = max(len(str(cell.value or '')) for cell in column_cells)
                    adjusted_width = min(length + 2, 50)
                    worksheet.column_dimensions[column_cells[0].column_letter].width = adjusted_width
        
        return excel_path
    
    except Exception as e:
        print(f"Warning: Could not create Excel file: {e}")
        return None

# ------------------------
# Output Functions
# ------------------------
def create_output_structure(outdir: str, db_name: str, regions: str, 
                           thresholds: Dict, tag: Optional[str] = None) -> Tuple[str, str, str]:
    """Create organized output folder and filename."""
    os.makedirs(outdir, exist_ok=True)
    
    date_str = dt.datetime.now().strftime("%Y%m%d")
    folder_name = f"{date_str}_{db_name.lower()}"
    output_folder = os.path.join(outdir, folder_name)
    os.makedirs(output_folder, exist_ok=True)
    
    csv_folder = os.path.join(output_folder, "csvs")
    os.makedirs(csv_folder, exist_ok=True)
    
    filename_parts = ["npz_scan", db_name.lower(), regions.replace(',', '_')]
    
    for region, threshold in thresholds.items():
        if threshold is not None:
            region_name = region.replace('min_id_', '')
            filename_parts.append(f"{region_name}{int(threshold*100)}")
    
    if tag:
        filename_parts.append(tag)
    
    time_str = dt.datetime.now().strftime("%H%M%S")
    filename_parts.append(time_str)
    
    filename = "_".join(filename_parts) + ".csv"
    csv_path = os.path.join(csv_folder, filename)
    
    return csv_path, output_folder, csv_folder

def write_metadata_block(fh, args, query_cdrs: Dict, start_time, config: Dict):
    """Write metadata block to CSV file."""
    fh.write("# ============================================\n")
    fh.write("# NPZ FULLSCAN v8 - VHH Database Search\n")
    fh.write("# ============================================\n")
    fh.write(f"# Run started: {start_time.strftime('%Y-%m-%d %H:%M:%S')}\n")
    fh.write(f"# Database: {config.get('db_name', 'Unknown')}\n")
    fh.write(f"# Database format: {config.get('db_format', 'Unknown')}\n")
    fh.write("#\n# QUERY INFORMATION\n")
    fh.write(f"# Query sequence: {query_cdrs.get('full', 'N/A')[:80]}...\n")
    fh.write(f"# Query CDR1: {query_cdrs.get('cdr1', 'N/A')}\n")
    fh.write(f"# Query CDR2: {query_cdrs.get('cdr2', 'N/A')}\n")
    fh.write(f"# Query CDR3: {query_cdrs.get('cdr3', 'N/A')}\n")
    fh.write("#\n# SEARCH PARAMETERS\n")
    fh.write(f"# Regions searched: {config['regions']}\n")
    
    for region in config['regions'].split(','):
        region = region.strip()
        min_id = config.get(f'min_id_{region}')
        if min_id is not None:
            fh.write(f"# {region.upper()} threshold: {min_id:.2%}\n")
    
    fh.write("# ============================================\n\n")

# ------------------------
# Histogram Visualization
# ------------------------
def print_ascii_histogram(values: List[float], bin_width: float = 1.0, label: str = "", max_bar_width: int = 50):
    """Print an ASCII histogram for numeric values."""
    if not values:
        print(f"No data for {label}")
        return
    
    clean_values = [float(v) for v in values if v is not None]
    if not clean_values:
        print(f"No valid data for {label}")
        return
    
    min_val = min(clean_values)
    max_val = max(clean_values)
    
    n_bins = int((max_val - min_val) / bin_width) + 1
    bins = [min_val + i * bin_width for i in range(n_bins + 1)]
    
    counts = [0] * n_bins
    for val in clean_values:
        bin_idx = min(int((val - min_val) / bin_width), n_bins - 1)
        counts[bin_idx] += 1
    
    max_count = max(counts) if counts else 1
    
    print(f"\n{label}")
    print("-" * 60)
    
    for i, count in enumerate(counts):
        if count == 0:
            continue
        
        bar_width = int((count / max_count) * max_bar_width)
        bar = "█" * bar_width if bar_width > 0 else ""
        
        bin_start = bins[i]
        bin_end = bins[i + 1]
        
        if bin_width == 1.0:
            bin_label = f"{int(bin_start):3d}"
        else:
            bin_label = f"{bin_start:5.1f}-{bin_end:5.1f}"
        
        print(f"  {bin_label}: {bar} ({count})")
    
    print(f"\n  Total: {len(clean_values)}")
    print(f"  Range: {min_val:.1f} - {max_val:.1f}")
    print(f"  Mean: {np.mean(clean_values):.1f}")
    print(f"  Median: {np.median(clean_values):.1f}")

# ------------------------
# Sampling and Analysis Functions
# ------------------------
def sample_vhh_shard(npz_file: str, query_cdrs: Dict, sample_size: int, 
                     regions: List[str]) -> pd.DataFrame:
    """Sample sequences from a VHH format shard and calculate identities."""
    try:
        data = np.load(npz_file, allow_pickle=True)
        
        aa_full = data.get("aa_v_full", np.array([]))
        cdr1_arr = data.get("cdr1", np.array([]))
        cdr2_arr = data.get("cdr2", np.array([]))
        cdr3_arr = data.get("cdr3", np.array([]))
        
        n_seqs = len(aa_full) if len(aa_full) > 0 else 0
        
        if n_seqs == 0:
            data.close()
            return pd.DataFrame()
        
        # Random sample
        if n_seqs <= sample_size:
            indices = list(range(n_seqs))
        else:
            indices = random.sample(range(n_seqs), sample_size)
        
        def safe_str(arr, idx):
            if arr is None or len(arr) == 0 or idx >= len(arr):
                return ""
            val = arr[idx]
            if val is None:
                return ""
            if isinstance(val, bytes):
                return val.decode('utf-8', errors='ignore')
            return str(val) if val else ""
        
        results = []
        for idx in indices:
            full_seq = safe_str(aa_full, idx)
            db_cdr1 = safe_str(cdr1_arr, idx)
            db_cdr2 = safe_str(cdr2_arr, idx)
            db_cdr3 = safe_str(cdr3_arr, idx)
            
            result = {"shard": os.path.basename(npz_file), "index": idx}
            result['len_full'] = len(full_seq) if full_seq else 0
            result['len_cdr1'] = len(db_cdr1) if db_cdr1 else 0
            result['len_cdr2'] = len(db_cdr2) if db_cdr2 else 0
            result['len_cdr3'] = len(db_cdr3) if db_cdr3 else 0
            
            # Store actual sequences for top matches display
            result['seq_full'] = full_seq
            result['seq_cdr1'] = db_cdr1
            result['seq_cdr2'] = db_cdr2
            result['seq_cdr3'] = db_cdr3
            
            # Get query CDRs - handle None values explicitly
            q_cdr1 = query_cdrs.get("cdr1") or ""
            q_cdr2 = query_cdrs.get("cdr2") or ""
            q_cdr3 = query_cdrs.get("cdr3") or ""
            q_full = query_cdrs.get("full") or ""
            
            # Calculate both exact and weighted identity for ALL CDRs
            exact_cdr1, weighted_cdr1 = calc_identity_both(q_cdr1, db_cdr1) if q_cdr1 and db_cdr1 else (0.0, 0.0)
            exact_cdr2, weighted_cdr2 = calc_identity_both(q_cdr2, db_cdr2) if q_cdr2 and db_cdr2 else (0.0, 0.0)
            exact_cdr3, weighted_cdr3 = calc_identity_both(q_cdr3, db_cdr3) if q_cdr3 and db_cdr3 else (0.0, 0.0)
            # Use region-by-region for full sequence identity
            if q_full and full_seq:
                exact_full, weighted_full = calc_region_by_region_identity(
                    q_full, q_cdr1, q_cdr2, q_cdr3,
                    full_seq, db_cdr1, db_cdr2, db_cdr3
                )
            else:
                exact_full, weighted_full = 0.0, 0.0
            
            # Store both types
            result["id_cdr1"] = weighted_cdr1
            result["id_cdr2"] = weighted_cdr2
            result["id_cdr3"] = weighted_cdr3
            result["id_full"] = weighted_full
            
            result["exact_id_cdr1"] = exact_cdr1
            result["exact_id_cdr2"] = exact_cdr2
            result["exact_id_cdr3"] = exact_cdr3
            result["exact_id_full"] = exact_full
            
            # Average identity is based on SEARCHED regions only
            searched_weighted = []
            searched_exact = []
            for region in regions:
                if region == "cdr1" and result["id_cdr1"] > 0:
                    searched_weighted.append(result["id_cdr1"])
                    searched_exact.append(result["exact_id_cdr1"])
                elif region == "cdr2" and result["id_cdr2"] > 0:
                    searched_weighted.append(result["id_cdr2"])
                    searched_exact.append(result["exact_id_cdr2"])
                elif region == "cdr3" and result["id_cdr3"] > 0:
                    searched_weighted.append(result["id_cdr3"])
                    searched_exact.append(result["exact_id_cdr3"])
                elif region == "full" and result["id_full"] > 0:
                    searched_weighted.append(result["id_full"])
                    searched_exact.append(result["exact_id_full"])
            
            result["avg_id"] = np.mean(searched_weighted) if searched_weighted else 0.0
            result["avg_exact_id"] = np.mean(searched_exact) if searched_exact else 0.0
            
            # Calculate weighted score for multi-CDR searches
            if len(regions) > 1 and "full" not in regions:
                region_ids = {r: result.get(f"id_{r}", 0) for r in regions}
                result["weighted_score"] = calc_weighted_cdr_score(region_ids, regions)
            else:
                result["weighted_score"] = result["avg_id"]
            
            results.append(result)
        
        data.close()
        return pd.DataFrame(results)
    
    except Exception as e:
        print(f"Error sampling {npz_file}: {e}")
        return pd.DataFrame()

def sample_oas_shard(npz_file: str, query_cdrs: Dict, sample_size: int, 
                     regions: List[str]) -> pd.DataFrame:
    """Sample sequences from an OAS format shard and calculate identities."""
    try:
        data = np.load(npz_file, allow_pickle=True)
        
        if "numberings" in data:
            arr = data["numberings"]
        elif "arr_0" in data.files:
            arr = data["arr_0"]
        else:
            data.close()
            return pd.DataFrame()
        
        n_seqs = len(arr)
        if n_seqs == 0:
            data.close()
            return pd.DataFrame()
        
        if n_seqs <= sample_size:
            indices = list(range(n_seqs))
        else:
            indices = random.sample(range(n_seqs), sample_size)
        
        results = []
        for idx in indices:
            row = arr[idx]
            result = {"shard": os.path.basename(npz_file), "index": idx}
            
            full_seq = extract_whole_from_npz(row)
            db_cdr1 = extract_cdr1_from_npz(row)
            db_cdr2 = extract_cdr2_from_npz(row)
            db_cdr3 = extract_cdr3_from_npz(row)
            
            result['len_full'] = len(full_seq) if full_seq else 0
            result['len_cdr1'] = len(db_cdr1) if db_cdr1 else 0
            result['len_cdr2'] = len(db_cdr2) if db_cdr2 else 0
            result['len_cdr3'] = len(db_cdr3) if db_cdr3 else 0
            
            # Store sequences for top matches display
            result['seq_full'] = full_seq
            result['seq_cdr1'] = db_cdr1
            result['seq_cdr2'] = db_cdr2
            result['seq_cdr3'] = db_cdr3
            
            # Get query CDRs - handle None values explicitly
            q_cdr1 = query_cdrs.get("cdr1") or ""
            q_cdr2 = query_cdrs.get("cdr2") or ""
            q_cdr3 = query_cdrs.get("cdr3") or ""
            q_full = query_cdrs.get("full") or ""
            
            # Calculate both exact and weighted identity for ALL CDRs
            exact_cdr1, weighted_cdr1 = calc_identity_both(q_cdr1, db_cdr1) if q_cdr1 and db_cdr1 else (0.0, 0.0)
            exact_cdr2, weighted_cdr2 = calc_identity_both(q_cdr2, db_cdr2) if q_cdr2 and db_cdr2 else (0.0, 0.0)
            exact_cdr3, weighted_cdr3 = calc_identity_both(q_cdr3, db_cdr3) if q_cdr3 and db_cdr3 else (0.0, 0.0)
            # Use region-by-region for full sequence identity
            if q_full and full_seq:
                exact_full, weighted_full = calc_region_by_region_identity(
                    q_full, q_cdr1, q_cdr2, q_cdr3,
                    full_seq, db_cdr1, db_cdr2, db_cdr3
                )
            else:
                exact_full, weighted_full = 0.0, 0.0
            
            # Store both types
            result["id_cdr1"] = weighted_cdr1
            result["id_cdr2"] = weighted_cdr2
            result["id_cdr3"] = weighted_cdr3
            result["id_full"] = weighted_full
            
            result["exact_id_cdr1"] = exact_cdr1
            result["exact_id_cdr2"] = exact_cdr2
            result["exact_id_cdr3"] = exact_cdr3
            result["exact_id_full"] = exact_full
            
            # Average identity is based on SEARCHED regions only
            searched_weighted = []
            searched_exact = []
            for region in regions:
                if region == "cdr1" and result["id_cdr1"] > 0:
                    searched_weighted.append(result["id_cdr1"])
                    searched_exact.append(result["exact_id_cdr1"])
                elif region == "cdr2" and result["id_cdr2"] > 0:
                    searched_weighted.append(result["id_cdr2"])
                    searched_exact.append(result["exact_id_cdr2"])
                elif region == "cdr3" and result["id_cdr3"] > 0:
                    searched_weighted.append(result["id_cdr3"])
                    searched_exact.append(result["exact_id_cdr3"])
                elif region == "full" and result["id_full"] > 0:
                    searched_weighted.append(result["id_full"])
                    searched_exact.append(result["exact_id_full"])
            
            result["avg_id"] = np.mean(searched_weighted) if searched_weighted else 0.0
            result["avg_exact_id"] = np.mean(searched_exact) if searched_exact else 0.0
            
            # Calculate weighted score for multi-CDR searches
            if len(regions) > 1 and "full" not in regions:
                region_ids = {r: result.get(f"id_{r}", 0) for r in regions}
                result["weighted_score"] = calc_weighted_cdr_score(region_ids, regions)
            else:
                result["weighted_score"] = result["avg_id"]
            
            results.append(result)
        
        data.close()
        return pd.DataFrame(results)
    
    except Exception as e:
        print(f"Error sampling {npz_file}: {e}")
        return pd.DataFrame()

def sample_shard(npz_file: str, query_cdrs: Dict, sample_size: int, 
                 regions: List[str]) -> pd.DataFrame:
    """Sample from shard, auto-detecting format."""
    fmt = detect_npz_format(npz_file)
    if fmt == "vhh":
        return sample_vhh_shard(npz_file, query_cdrs, sample_size, regions)
    else:
        return sample_oas_shard(npz_file, query_cdrs, sample_size, regions)

def run_interactive_analysis(npz_files: List[str], query_cdrs: Dict, regions: List[str],
                            sample_size: int = 5000, target_hits: int = 5000,
                            db_name: str = "Unknown") -> Tuple[Dict, Dict]:
    """
    Run interactive analysis to determine optimal thresholds.
    Returns (recommendations, length_windows) dictionaries.
    """
    
    print("\n" + "="*80)
    print("INTERACTIVE THRESHOLD ANALYSIS")
    print("="*80)
    print(f"Database: {db_name}")
    print(f"Sampling {sample_size:,} sequences per shard...")
    print(f"Target: Find thresholds to get ~{target_hits:,} hits")
    print("="*80 + "\n")
    
    # Estimate database size
    if "vhh" in db_name.lower():
        estimated_db_size = 12_000_000
    else:
        avg_seqs_per_shard = 500_000
        estimated_db_size = len(npz_files) * avg_seqs_per_shard
    
    print(f"Estimated database size: ~{estimated_db_size:,} sequences")
    
    # Limit shards to sample
    shards_to_sample = min(len(npz_files), 8)
    print(f"Sampling from {shards_to_sample} shards")
    print(f"Total sequences to sample: ~{shards_to_sample * sample_size:,}\n")
    
    # Sample from shards
    all_samples = []
    sample_bar = tqdm(npz_files[:shards_to_sample], desc="Sampling shards")
    
    for npz_file in sample_bar:
        df = sample_shard(npz_file, query_cdrs, sample_size, regions)
        if not df.empty:
            all_samples.append(df)
    
    if not all_samples:
        print("ERROR: No samples collected!")
        return {}, {}
    
    combined = pd.concat(all_samples, ignore_index=True)
    
    # Check if this is a multi-CDR search
    is_multi_cdr = len(regions) > 1 and "full" not in regions
    
    # Analyze distributions
    print(f"\n{'='*80}")
    print("IDENTITY DISTRIBUTION ANALYSIS")
    print(f"{'='*80}")
    print("\nWeighted identity gives partial credit for similar amino acids:")
    print("  Conservative (75%): D↔E, K↔R, L↔I↔V, F↔Y, S↔T, N↔Q")
    print("  Semi-conserv (50%): H↔K/R, M↔L/I/V, W↔F/Y, A↔G/S")
    
    if is_multi_cdr:
        print(f"\n{'─'*60}")
        print("CDR WEIGHTING (for multi-CDR searches):")
        print("  CDR1: 20% weight")
        print("  CDR2: 30% weight") 
        print("  CDR3: 50% weight (most important for binding)")
        print(f"{'─'*60}")
    
    recommendations = {}
    
    # Show per-region statistics
    for region in regions:
        col = f"id_{region}"
        exact_col = f"exact_id_{region}"
        
        if col not in combined.columns:
            continue
        
        values = combined[col].dropna().values
        exact_values = combined[exact_col].dropna().values if exact_col in combined.columns else values
        
        if len(values) == 0:
            print(f"\n{region.upper()} Statistics:")
            print(f"  No valid values extracted for {region.upper()}")
            continue
        
        weight_str = f" (weight: {CDR_WEIGHTS.get(region, 0.33)*100:.0f}%)" if is_multi_cdr else ""
        print(f"\n{region.upper()} Statistics{weight_str} (Weighted / Exact):")
        print(f"  Samples with data: {len(values):,}/{len(combined):,} ({len(values)*100/len(combined):.1f}%)")
        print(f"  Max:    {values.max()*100:5.1f}% / {exact_values.max()*100:5.1f}%")
        print(f"  Q95:    {np.percentile(values, 95)*100:5.1f}% / {np.percentile(exact_values, 95)*100:5.1f}%")
        print(f"  Q90:    {np.percentile(values, 90)*100:5.1f}% / {np.percentile(exact_values, 90)*100:5.1f}%")
        print(f"  Q75:    {np.percentile(values, 75)*100:5.1f}% / {np.percentile(exact_values, 75)*100:5.1f}%")
        print(f"  Median: {np.median(values)*100:5.1f}% / {np.median(exact_values)*100:5.1f}%")
        print(f"  Mean:   {values.mean()*100:5.1f}% / {exact_values.mean()*100:5.1f}%")
    
    # For multi-CDR searches, show WEIGHTED SCORE distribution and recommend a single threshold
    if is_multi_cdr and "weighted_score" in combined.columns:
        print(f"\n{'='*60}")
        print("WEIGHTED CDR SCORE (0.2×CDR1 + 0.3×CDR2 + 0.5×CDR3)")
        print(f"{'='*60}")
        
        weighted_vals = combined["weighted_score"].dropna().values
        
        if len(weighted_vals) > 0:
            print(f"  Samples: {len(weighted_vals):,}")
            print(f"  Max:    {weighted_vals.max()*100:5.1f}%")
            print(f"  Q95:    {np.percentile(weighted_vals, 95)*100:5.1f}%")
            print(f"  Q90:    {np.percentile(weighted_vals, 90)*100:5.1f}%")
            print(f"  Q75:    {np.percentile(weighted_vals, 75)*100:5.1f}%")
            print(f"  Median: {np.median(weighted_vals)*100:5.1f}%")
            print(f"  Mean:   {weighted_vals.mean()*100:5.1f}%")
            
            # Estimate threshold for target hits
            total_sampled = len(combined)
            scale_factor = estimated_db_size / total_sampled
            needed_in_sample = min(target_hits / scale_factor, len(weighted_vals))
            
            if needed_in_sample >= len(weighted_vals):
                recommended_weighted = np.percentile(weighted_vals, 90)
                estimated_hits = len(weighted_vals) * scale_factor * 0.1
                print(f"\n  ⚠ Warning: May not have {target_hits:,} matches at reasonable thresholds")
                print(f"\n  RECOMMENDED WEIGHTED THRESHOLD: {recommended_weighted*100:.1f}%")
                print(f"  This is the 90th percentile")
                print(f"  Estimated hits: ~{int(estimated_hits):,}")
            else:
                sorted_vals = np.sort(weighted_vals)[::-1]
                threshold_idx = int(needed_in_sample)
                recommended_weighted = sorted_vals[threshold_idx]
                percentile_rank = (1 - (threshold_idx / len(weighted_vals))) * 100
                
                print(f"\n  RECOMMENDED WEIGHTED THRESHOLD: {recommended_weighted*100:.1f}%")
                print(f"  This is the {percentile_rank:.1f}th percentile")
                print(f"  Estimated to give ~{target_hits:,} hits")
            
            recommendations["weighted"] = recommended_weighted
            
            # Ask about per-CDR minimum floors
            print(f"\n{'─'*60}")
            print("OPTIONAL: Per-CDR Minimum Floors")
            print("  These ensure no single CDR is too dissimilar")
            print(f"{'─'*60}")
            
            use_floors = input("\nSet minimum floor for each CDR? (y/n) [n]: ").strip().lower()
            if use_floors in ['y', 'yes']:
                for region in regions:
                    col = f"id_{region}"
                    if col in combined.columns:
                        vals = combined[col].dropna().values
                        suggested_floor = np.percentile(vals, 50)  # Median as default floor
                        
                        while True:
                            floor_input = input(f"  {region.upper()} minimum floor [{suggested_floor*100:.0f}%]: ").strip()
                            if not floor_input:
                                recommendations[f"floor_{region}"] = suggested_floor
                                break
                            try:
                                floor_val = float(floor_input.replace("%", ""))
                                if floor_val > 1:
                                    floor_val /= 100
                                recommendations[f"floor_{region}"] = floor_val
                                break
                            except ValueError:
                                print("    Please enter a number")
    
    else:
        # Single region search - use per-region threshold
        for region in regions:
            col = f"id_{region}"
            if col not in combined.columns:
                continue
            
            values = combined[col].dropna().values
            if len(values) == 0:
                continue
            
            # Estimate threshold for target hits
            total_sampled = len(combined)
            scale_factor = estimated_db_size / total_sampled
            needed_in_sample = min(target_hits / scale_factor, len(values))
            
            if needed_in_sample >= len(values):
                if len(values) > 10:
                    recommended_threshold = np.percentile(values, 90)
                    estimated_hits = len(values) * scale_factor * 0.1
                    print(f"\n  ⚠ Warning: May not have {target_hits:,} matches")
                    print(f"\n  RECOMMENDED {region.upper()} THRESHOLD: {recommended_threshold*100:.1f}%")
                    print(f"  Estimated hits: ~{int(estimated_hits):,}")
                else:
                    recommended_threshold = np.percentile(values, 50)
                    print(f"\n  RECOMMENDED {region.upper()} THRESHOLD: {recommended_threshold*100:.1f}%")
                
                recommendations[region] = recommended_threshold
            else:
                sorted_vals = np.sort(values)[::-1]
                threshold_idx = int(needed_in_sample)
                recommended_threshold = sorted_vals[threshold_idx]
                percentile_rank = (1 - (threshold_idx / len(values))) * 100
                
                print(f"\n  RECOMMENDED {region.upper()} THRESHOLD: {recommended_threshold*100:.1f}%")
                print(f"  This is the {percentile_rank:.1f}th percentile")
                print(f"  Estimated to give ~{target_hits:,} hits")
                
                recommendations[region] = recommended_threshold
    
    # CDR Length Analysis and Filtering
    length_windows = {}
    
    print(f"\n{'='*80}")
    print("CDR LENGTH ANALYSIS")
    print(f"{'='*80}")
    
    for cdr in ['cdr1', 'cdr2', 'cdr3']:
        if cdr not in regions:
            continue
        
        len_col = f'len_{cdr}'
        if len_col not in combined.columns:
            continue
        
        cdr_lengths = [l for l in combined[len_col].dropna().tolist() if l > 0]
        if not cdr_lengths:
            continue
        
        query_len = len(query_cdrs.get(cdr, ""))
        
        print(f"\n{cdr.upper()} Length Distribution:")
        print_ascii_histogram(cdr_lengths, bin_width=1, label="")
        
        print(f"\n  Your query {cdr.upper()} length: {query_len} aa")
        
        # Show how many sequences are within different ranges
        for window in [1, 2, 3, 5]:
            in_range = sum(1 for l in cdr_lengths if abs(l - query_len) <= window)
            pct = in_range / len(cdr_lengths) * 100
            print(f"    ±{window} aa ({query_len-window}-{query_len+window}): {in_range:,} sequences ({pct:.1f}%)")
        
        # Ask for length window
        while True:
            window_input = input(f"\n  {cdr.upper()} length window (±aa, or Enter for ±3): ").strip()
            if not window_input:
                length_windows[cdr] = 3
                break
            try:
                length_windows[cdr] = int(window_input)
                break
            except ValueError:
                print("    Please enter a number")
    
    # Show top matches with alignment
    print(f"\n{'='*80}")
    print("TOP 10 MATCHES IN SAMPLE")
    print(f"{'='*80}")
    print("Legend: | = exact (100%), : = conservative (75%), . = semi-conserv (50%), x = diff")
    
    # Determine sort column based on search type
    sort_col = 'weighted_score' if (is_multi_cdr and 'weighted_score' in combined.columns) else 'avg_id'
    
    if sort_col in combined.columns:
        top10 = combined.nlargest(10, sort_col)
        
        for i, (_, row) in enumerate(top10.iterrows(), 1):
            # Recalculate all identities at display time for consistency
            db_cdr1 = row.get('seq_cdr1', '') or ''
            db_cdr2 = row.get('seq_cdr2', '') or ''
            db_cdr3 = row.get('seq_cdr3', '') or ''
            db_full = row.get('seq_full', '') or ''
            
            q_cdr1 = query_cdrs.get('cdr1', '') or ''
            q_cdr2 = query_cdrs.get('cdr2', '') or ''
            q_cdr3 = query_cdrs.get('cdr3', '') or ''
            
            # Calculate both exact and weighted identities
            exact_cdr1, weighted_cdr1 = calc_identity_both(q_cdr1, db_cdr1) if q_cdr1 and db_cdr1 else (0.0, 0.0)
            exact_cdr2, weighted_cdr2 = calc_identity_both(q_cdr2, db_cdr2) if q_cdr2 and db_cdr2 else (0.0, 0.0)
            exact_cdr3, weighted_cdr3 = calc_identity_both(q_cdr3, db_cdr3) if q_cdr3 and db_cdr3 else (0.0, 0.0)
            
            # Calculate weighted CDR score for multi-CDR searches
            if is_multi_cdr:
                region_ids = {}
                if 'cdr1' in regions:
                    region_ids['cdr1'] = weighted_cdr1
                if 'cdr2' in regions:
                    region_ids['cdr2'] = weighted_cdr2
                if 'cdr3' in regions:
                    region_ids['cdr3'] = weighted_cdr3
                weighted_score = calc_weighted_cdr_score(region_ids, regions)
                
                # Also calc exact weighted score
                exact_region_ids = {}
                if 'cdr1' in regions:
                    exact_region_ids['cdr1'] = exact_cdr1
                if 'cdr2' in regions:
                    exact_region_ids['cdr2'] = exact_cdr2
                if 'cdr3' in regions:
                    exact_region_ids['cdr3'] = exact_cdr3
                exact_weighted_score = calc_weighted_cdr_score(exact_region_ids, regions)
            elif 'full' in regions:
                # Full sequence search - use stored identity from dataframe
                q_full = query_cdrs.get('full', '') or ''
                if q_full and db_full:
                    exact_full, weighted_full = calc_identity_both(q_full, db_full)
                else:
                    exact_full, weighted_full = row.get('exact_id_full', 0), row.get('id_full', 0)
                weighted_score = weighted_full
                exact_weighted_score = exact_full
            else:
                # Single CDR region search
                searched_weighted = []
                searched_exact = []
                for region in regions:
                    if region == 'cdr1' and weighted_cdr1 > 0:
                        searched_weighted.append(weighted_cdr1)
                        searched_exact.append(exact_cdr1)
                    elif region == 'cdr2' and weighted_cdr2 > 0:
                        searched_weighted.append(weighted_cdr2)
                        searched_exact.append(exact_cdr2)
                    elif region == 'cdr3' and weighted_cdr3 > 0:
                        searched_weighted.append(weighted_cdr3)
                        searched_exact.append(exact_cdr3)
                
                weighted_score = np.mean(searched_weighted) if searched_weighted else 0
                exact_weighted_score = np.mean(searched_exact) if searched_exact else 0
            
            print(f"\n{'─'*60}")
            if is_multi_cdr:
                print(f"#{i} - Weighted CDR Score: {weighted_score*100:.1f}% (exact: {exact_weighted_score*100:.1f}%)")
                print(f"     Scoring: 0.2×CDR1 + 0.3×CDR2 + 0.5×CDR3")
            elif 'full' in regions:
                print(f"#{i} - Full Sequence Identity: {weighted_score*100:.1f}% (exact: {exact_weighted_score*100:.1f}%)")
            else:
                print(f"#{i} - Identity: {weighted_score*100:.1f}% (exact: {exact_weighted_score*100:.1f}%)")
            print(f"   Shard: {row['shard']} (index {row['index']:,})")
            
            # Show all CDR alignments with their identities
            cdr_data = [
                ('CDR1', q_cdr1, db_cdr1, exact_cdr1, weighted_cdr1, CDR_WEIGHTS.get('cdr1', 0.2)),
                ('CDR2', q_cdr2, db_cdr2, exact_cdr2, weighted_cdr2, CDR_WEIGHTS.get('cdr2', 0.3)),
                ('CDR3', q_cdr3, db_cdr3, exact_cdr3, weighted_cdr3, CDR_WEIGHTS.get('cdr3', 0.5)),
            ]
            
            for cdr_name, q_seq, db_seq, exact_id, weighted_id, weight in cdr_data:
                if q_seq and db_seq:
                    q_disp, match_line, t_disp = highlight_alignment(q_seq, db_seq)
                    
                    # Mark if this CDR was searched, show weight for multi-CDR
                    if cdr_name.lower() in regions:
                        if is_multi_cdr:
                            searched_marker = f" ★ ({weight*100:.0f}% weight)"
                        else:
                            searched_marker = " ★"
                    else:
                        searched_marker = ""
                    
                    print(f"\n   {cdr_name} (exact: {exact_id*100:.1f}%, weighted: {weighted_id*100:.1f}%){searched_marker}:")
                    print(f"   Query:  {q_disp}")
                    print(f"           {match_line}")
                    print(f"   Match:  {t_disp}")
            
            # Show full sequence (wrapped)
            if db_full:
                print(f"\n   Full sequence ({len(db_full)} aa):")
                for j in range(0, len(db_full), 60):
                    print(f"   {db_full[j:j+60]}")
    
    print(f"\n{'='*80}\n")
    
    return recommendations, length_windows

# ------------------------
# Interactive Mode
# ------------------------
def get_user_input_interactive():
    """Get all parameters interactively from the user."""
    print("\n" + "="*80)
    print("NPZ FULLSCAN v8 - VHH DATABASE SEARCH")
    print("="*80)
    print("Welcome! Let's set up your antibody search.\n")
    
    params = {}
    
    # Check for available databases
    print("Step 1: Database Selection")
    print("-" * 40)
    
    # Find VHH shards
    vhh_path = None
    for path in VHH_SHARD_PATHS:
        if os.path.exists(path):
            npz_files = glob.glob(os.path.join(path, "*.npz"))
            if npz_files:
                vhh_path = path
                break
    
    # Find OAS databases
    oas_paths = {}
    for base_path in OAS_DB_PATHS:
        if os.path.exists(base_path):
            for species_dir in glob.glob(os.path.join(base_path, "*")):
                if os.path.isdir(species_dir):
                    species = os.path.basename(species_dir)
                    oas_paths[species] = species_dir
    
    # Check last used database
    last_db = load_last_db()
    
    # Build menu
    print("\nAvailable databases:")
    options = {}
    opt_num = 1
    
    if vhh_path:
        shard_index = load_shard_index(vhh_path)
        total_seqs = shard_index.get("total_sequences", "Unknown") if shard_index else "Unknown"
        print(f"  [{opt_num}] VHH Unified Database ({total_seqs:,} sequences)")
        print(f"      - 12M+ VHH/Nanobody sequences from multiple sources")
        print(f"      - Includes target/antigen binding info")
        print(f"      - Patent IDs and titles")
        options[str(opt_num)] = ("vhh", vhh_path)
        opt_num += 1
        
        # Add annotated-only option
        print(f"  [{opt_num}] VHH Annotated Only (~19K sequences)")
        print(f"      - High-value sequences with target/patent info")
        print(f"      - Fast search for known binders")
        options[str(opt_num)] = ("vhh_annotated", vhh_path)
        opt_num += 1
    
    for species, path in sorted(oas_paths.items()):
        npz_count = len(glob.glob(os.path.join(path, "*.npz")))
        print(f"  [{opt_num}] OAS {species} ({npz_count} shards)")
        options[str(opt_num)] = ("oas", path, species)
        opt_num += 1
    
    if not options:
        print("\n❌ No databases found!")
        print("Please specify database path with --vhh-shards or --db-root")
        sys.exit(1)
    
    # Get selection
    while True:
        choice = input(f"\nSelect database [1-{opt_num-1}]: ").strip()
        if choice in options:
            break
        print("Invalid selection, please try again.")
    
    selection = options[choice]
    
    if selection[0] == "vhh":
        params["vhh_shards"] = selection[1]
        params["db_name"] = "VHH"
        params["annotated_only"] = False
    elif selection[0] == "vhh_annotated":
        params["vhh_shards"] = selection[1]
        params["db_name"] = "VHH_Annotated"
        params["annotated_only"] = True
    else:
        params["db_root"] = selection[1]
        params["db_name"] = selection[2]
        params["annotated_only"] = False
    
    save_last_db(selection[1])
    
    # Get query sequence
    print("\n" + "-"*40)
    print("Step 2: Query Sequence")
    print("-"*40)
    
    last_query = load_last_query_file()
    
    print("\nEnter your query:")
    print("  1. Paste sequence directly")
    print("  2. Load from FASTA file")
    if last_query:
        print(f"  3. Use last file: {last_query}")
    
    while True:
        choice = input("\nChoice [1/2" + ("/3" if last_query else "") + "]: ").strip()
        
        if choice == "1":
            print("\nPaste your sequence (press Enter twice when done):")
            lines = []
            while True:
                line = input()
                if not line:
                    break
                lines.append(line)
            seq = "".join(lines).replace(" ", "").replace("\n", "")
            if seq.startswith(">"):
                seq = "".join(seq.split("\n")[1:])
            params["query_seq"] = seq.upper()
            break
        
        elif choice == "2":
            fasta_path = input("Enter FASTA file path: ").strip()
            # Handle relative paths
            if not os.path.isabs(fasta_path):
                fasta_path = os.path.abspath(fasta_path)
            
            if os.path.exists(fasta_path):
                with open(fasta_path, 'r') as f:
                    content = f.read()
                lines = content.strip().split('\n')
                seq_lines = [l for l in lines if not l.startswith('>')]
                params["query_seq"] = "".join(seq_lines).replace(" ", "").upper()
                save_last_query_file(fasta_path)
                print(f"\n  ✓ Loaded from: {fasta_path}")
                break
            else:
                print(f"  ❌ File not found: {fasta_path}")
        
        elif choice == "3" and last_query:
            with open(last_query, 'r') as f:
                content = f.read()
            lines = content.strip().split('\n')
            seq_lines = [l for l in lines if not l.startswith('>')]
            params["query_seq"] = "".join(seq_lines).replace(" ", "").upper()
            print(f"\n  ✓ Loaded from: {last_query}")
            break
        
        else:
            print("  Invalid choice")
    
    # Show full query sequence
    print(f"\n  Query sequence ({len(params['query_seq'])} aa):")
    # Show full sequence, wrapping at 60 characters
    seq = params["query_seq"]
    for i in range(0, len(seq), 60):
        print(f"  {seq[i:i+60]}")
    
    # Detect and offer to remove tags
    params["query_seq"] = interactive_tag_removal(params["query_seq"])
    
    # Extract CDRs
    print("\nExtracting CDRs from query...")
    query_cdrs = extract_cdrs_from_query_sequence(params["query_seq"])
    query_cdrs["full"] = params["query_seq"]
    
    print(f"  CDR1: {query_cdrs['cdr1']} ({len(query_cdrs['cdr1'])} aa)")
    print(f"  CDR2: {query_cdrs['cdr2']} ({len(query_cdrs['cdr2'])} aa)")
    print(f"  CDR3: {query_cdrs['cdr3']} ({len(query_cdrs['cdr3'])} aa)")
    
    params["query_cdrs"] = query_cdrs
    
    # Region selection
    print("\n" + "-"*40)
    print("Step 3: Search Regions")
    print("-"*40)
    print("\nWhich regions to search?")
    print("  1. CDR3 only (fastest, recommended for initial search)")
    print("  2. All CDRs (CDR1, CDR2, CDR3)")
    print("  3. Full sequence")
    print("  4. Custom")
    
    while True:
        choice = input("\nChoice [1-4]: ").strip()
        if choice == "1":
            params["regions"] = "cdr3"
            break
        elif choice == "2":
            params["regions"] = "cdr1,cdr2,cdr3"
            break
        elif choice == "3":
            params["regions"] = "full"
            break
        elif choice == "4":
            regions = input("Enter regions (comma-separated, e.g., cdr1,cdr3): ").strip()
            params["regions"] = regions
            break
    
    # Threshold selection
    print("\n" + "-"*40)
    print("Step 4: Identity Thresholds")
    print("-"*40)
    print("\nHow do you want to set identity thresholds?")
    print("  1. Set manually (enter specific values)")
    print("  2. Use automatic threshold detection (RECOMMENDED)")
    print("     - Samples database to find optimal thresholds")
    print("     - Shows identity distribution statistics")
    print("  3. No thresholds (return all matches)")
    
    while True:
        choice = input("\nChoice [1-3]: ").strip()
        
        if choice == "1":
            # Manual thresholds
            for region in params["regions"].split(","):
                region = region.strip()
                default = 0.7 if region == "cdr3" else 0.6
                
                while True:
                    val = input(f"  {region.upper()} minimum identity [{default:.0%}]: ").strip()
                    if not val:
                        params[f"min_id_{region}"] = default
                        break
                    try:
                        val = float(val.replace("%", "")) 
                        if val > 1:
                            val = val / 100
                        params[f"min_id_{region}"] = val
                        break
                    except ValueError:
                        print("  Invalid value, enter a number (e.g., 0.7 or 70)")
            break
        
        elif choice == "2":
            # Automatic threshold detection
            params["run_analysis"] = True
            
            # Get analysis parameters
            print("\n  Automatic threshold detection settings:")
            
            # Sample size
            default_sample = 5000
            while True:
                sample_input = input(f"  Sample size per shard [{default_sample}]: ").strip()
                if not sample_input:
                    params["sample_size"] = default_sample
                    break
                try:
                    params["sample_size"] = int(sample_input)
                    break
                except ValueError:
                    print("    Invalid number")
            
            # Target hits
            default_target = 5000
            while True:
                target_input = input(f"  Target number of hits [{default_target}]: ").strip()
                if not target_input:
                    params["target_hits"] = default_target
                    break
                try:
                    params["target_hits"] = int(target_input)
                    break
                except ValueError:
                    print("    Invalid number")
            break
        
        elif choice == "3":
            # No thresholds
            print("  No thresholds will be applied - all matches returned")
            break
        
        else:
            print("Invalid choice, please enter 1, 2, or 3")
    
    # Output directory
    print("\n" + "-"*40)
    print("Step 5: Output Location")
    print("-"*40)
    
    default_outdir = os.path.expanduser("~/ka_search_results")
    outdir = input(f"Output directory [{default_outdir}]: ").strip()
    params["outdir"] = outdir if outdir else default_outdir
    
    return params

# ------------------------
# Main Function
# ------------------------
def main():
    parser = argparse.ArgumentParser(
        description="NPZ FULLSCAN v8 - VHH Database Search with Metadata",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Interactive mode
  python3 npz_fullscan_v8_vhh.py
  
  # VHH database search
  python3 npz_fullscan_v8_vhh.py --vhh-shards ~/KA-Search/VHH_shards --query-seq "EVQL..."
  
  # Search annotated sequences only (fast target lookup)
  python3 npz_fullscan_v8_vhh.py --vhh-shards ~/KA-Search/VHH_shards --annotated-only --query-seq "..."
  
  # OAS database search (legacy)
  python3 npz_fullscan_v8_vhh.py --db-root ~/KA-Search/extracted/.../Human --query-seq "..."
        """
    )
    
    # Database options
    db_group = parser.add_mutually_exclusive_group()
    db_group.add_argument("--vhh-shards", help="Path to VHH shards directory")
    db_group.add_argument("--db-root", help="Path to OAS database directory")
    
    # VHH-specific options
    parser.add_argument("--annotated-only", action="store_true",
                       help="Search only annotated sequences (with target/patent info)")
    
    # Query options
    parser.add_argument("--query-seq", help="Query antibody sequence")
    parser.add_argument("--query-file", help="FASTA file with query sequence")
    
    # Search parameters
    parser.add_argument("--regions", default="cdr3",
                       help="Regions to search (default: cdr3)")
    parser.add_argument("--min-id-cdr1", type=float, help="Min CDR1 identity")
    parser.add_argument("--min-id-cdr2", type=float, help="Min CDR2 identity")
    parser.add_argument("--min-id-cdr3", type=float, default=0.7,
                       help="Min CDR3 identity (default: 0.7)")
    parser.add_argument("--min-id-full", type=float, help="Min full sequence identity")
    
    # Filtering options
    parser.add_argument("--length-filter", action="store_true",
                       help="Enable CDR length filtering")
    parser.add_argument("--seq-length-prefilter", action="store_true",
                       help="Enable sequence length pre-filtering")
    
    # Output options
    parser.add_argument("--outdir", default=os.path.expanduser("~/ka_search_results"),
                       help="Output directory")
    parser.add_argument("--tag", help="Tag for output filename")
    
    args = parser.parse_args()
    
    # Interactive mode if no database specified
    run_analysis = False
    sample_size = 5000
    target_hits = 5000
    
    if not args.vhh_shards and not args.db_root:
        params = get_user_input_interactive()
        
        # Update args from interactive params
        if "vhh_shards" in params:
            args.vhh_shards = params["vhh_shards"]
        if "db_root" in params:
            args.db_root = params["db_root"]
        if "query_seq" in params:
            args.query_seq = params["query_seq"]
        if "regions" in params:
            args.regions = params["regions"]
        if "outdir" in params:
            args.outdir = params["outdir"]
        if params.get("annotated_only"):
            args.annotated_only = True
        
        # Check for analysis request
        run_analysis = params.get("run_analysis", False)
        sample_size = params.get("sample_size", 5000)
        target_hits = params.get("target_hits", 5000)
        
        for region in params.get("regions", "cdr3").split(","):
            region = region.strip()
            key = f"min_id_{region}"
            if key in params:
                setattr(args, key.replace("-", "_"), params[key])
        
        query_cdrs = params.get("query_cdrs", {})
        db_name = params.get("db_name", "Unknown")
    else:
        # Command-line mode - need query sequence
        if not args.query_seq and not args.query_file:
            print("Error: --query-seq or --query-file required")
            sys.exit(1)
        
        if args.query_file:
            with open(args.query_file, 'r') as f:
                content = f.read()
            lines = content.strip().split('\n')
            seq_lines = [l for l in lines if not l.startswith('>')]
            args.query_seq = "".join(seq_lines).replace(" ", "").upper()
        
        query_cdrs = extract_cdrs_from_query_sequence(args.query_seq)
        query_cdrs["full"] = args.query_seq
        
        db_name = "VHH" if args.vhh_shards else os.path.basename(args.db_root.rstrip('/'))
    
    # Get NPZ files
    if args.vhh_shards:
        npz_files = get_vhh_shards(args.vhh_shards, args.annotated_only)
        db_format = "vhh"
        print_vhh_database_info(args.vhh_shards)
    else:
        npz_files = sorted(glob.glob(os.path.join(args.db_root, "*.npz")))
        db_format = "oas"
    
    if not npz_files:
        print("Error: No NPZ files found!")
        sys.exit(1)
    
    print(f"\nFound {len(npz_files)} shard(s) to search")
    
    # Initialize length_windows (may be updated by analysis)
    analysis_length_windows = {}
    
    # Run interactive analysis if requested
    if run_analysis:
        regions_list = args.regions.split(",")
        regions_list = [r.strip() for r in regions_list]
        
        recommendations, analysis_length_windows = run_interactive_analysis(
            npz_files, query_cdrs, regions_list,
            sample_size=sample_size, target_hits=target_hits,
            db_name=db_name
        )
        
        if recommendations:
            print(f"\n{'='*80}")
            print("THRESHOLD RECOMMENDATIONS")
            print(f"{'='*80}")
            
            # Check if this is a multi-CDR weighted search
            is_multi_cdr = len(regions_list) > 1 and "full" not in regions_list
            
            if is_multi_cdr and "weighted" in recommendations:
                print(f"\n  WEIGHTED SCORE THRESHOLD: {recommendations['weighted']*100:.1f}%")
                print(f"  (Weighted: 0.2×CDR1 + 0.3×CDR2 + 0.5×CDR3)")
                
                # Show floors if set
                floors = {k: v for k, v in recommendations.items() if k.startswith("floor_")}
                if floors:
                    print(f"\n  Per-CDR minimum floors:")
                    for key, val in floors.items():
                        region = key.replace("floor_", "").upper()
                        print(f"    {region}: ≥{val*100:.0f}%")
            else:
                for region, threshold in recommendations.items():
                    if not region.startswith("floor_") and region != "weighted":
                        print(f"  {region.upper()}: {threshold*100:.1f}%")
            
            print(f"\n{'='*80}")
            use_rec = input("\nUse recommended thresholds? (y/n) [y]: ").strip().lower()
            
            if use_rec in ['y', 'yes', '']:
                print("\n✓ Using recommended thresholds")
                
                if is_multi_cdr and "weighted" in recommendations:
                    # Set weighted threshold in config
                    args.min_weighted_identity = recommendations["weighted"]
                    
                    # Set per-CDR floors if defined
                    for key, val in recommendations.items():
                        if key.startswith("floor_"):
                            region = key.replace("floor_", "")
                            setattr(args, f"min_id_{region}", val)
                else:
                    # Single region - use direct threshold
                    for region, threshold in recommendations.items():
                        if not region.startswith("floor_") and region != "weighted":
                            setattr(args, f"min_id_{region}", threshold)
            else:
                # Let user set manually
                print("\nEnter your thresholds:")
                
                if is_multi_cdr:
                    # Ask for weighted threshold
                    rec_weighted = recommendations.get("weighted", 0.5)
                    while True:
                        val = input(f"  Weighted score minimum [{rec_weighted:.1%}]: ").strip()
                        if not val:
                            args.min_weighted_identity = rec_weighted
                            break
                        try:
                            val = float(val.replace("%", ""))
                            if val > 1:
                                val = val / 100
                            args.min_weighted_identity = val
                            break
                        except ValueError:
                            print("    Invalid value")
                    
                    # Ask for per-CDR floors
                    set_floors = input("\n  Set per-CDR minimum floors? (y/n) [n]: ").strip().lower()
                    if set_floors in ['y', 'yes']:
                        for region in regions_list:
                            rec = recommendations.get(f"floor_{region}", 0.4)
                            while True:
                                val = input(f"    {region.upper()} minimum [{rec:.0%}]: ").strip()
                                if not val:
                                    setattr(args, f"min_id_{region}", rec)
                                    break
                                try:
                                    val = float(val.replace("%", ""))
                                    if val > 1:
                                        val = val / 100
                                    setattr(args, f"min_id_{region}", val)
                                    break
                                except ValueError:
                                    print("      Invalid value")
                else:
                    for region in regions_list:
                        rec = recommendations.get(region, 0.7)
                        while True:
                            val = input(f"  {region.upper()} minimum identity [{rec:.1%}]: ").strip()
                            if not val:
                                setattr(args, f"min_id_{region}", rec)
                                break
                            try:
                                val = float(val.replace("%", ""))
                                if val > 1:
                                    val = val / 100
                                setattr(args, f"min_id_{region}", val)
                                break
                            except ValueError:
                                print("    Invalid value")
    
    # Build config - use analysis length windows if available
    default_length_windows = {"cdr1": 3, "cdr2": 3, "cdr3": 5, "full": 20}
    if analysis_length_windows:
        default_length_windows.update(analysis_length_windows)
    
    # Check if multi-CDR weighted search
    regions_list = [r.strip() for r in args.regions.split(",")]
    is_multi_cdr = len(regions_list) > 1 and "full" not in regions_list
    
    # For full sequence search, show info about region-by-region method
    if "full" in regions_list:
        print(f"\n{'─'*60}")
        print("FULL SEQUENCE SEARCH")
        print(f"{'─'*60}")
        print("Using region-by-region identity calculation:")
        print("  - Compares FR1, CDR1, FR2, CDR2, FR3, CDR3, FR4 separately")
        print("  - Handles CDR length differences properly")
        print("  - ~5x faster than full edit distance")
        print(f"  - Estimated time: ~3-4 hours for 12M sequences")
        print(f"{'─'*60}")
    
    cfg = {
        "regions": args.regions,
        "db_name": db_name,
        "db_format": db_format,
        "use_length_filter": args.length_filter or bool(analysis_length_windows),
        "enable_seq_length_prefilter": args.seq_length_prefilter,
        # Enable CDR length prefilter by default for CDR searches (huge speedup)
        "enable_cdr_length_prefilter": "cdr" in args.regions.lower(),
        "seq_length_threshold": 0.25,
        "length_windows": default_length_windows,
        "query_lengths": {
            "cdr1": len(query_cdrs.get("cdr1", "")),
            "cdr2": len(query_cdrs.get("cdr2", "")),
            "cdr3": len(query_cdrs.get("cdr3", "")),
            "full": len(query_cdrs.get("full", "")),
        },
        # Weighted threshold for multi-CDR searches
        "min_weighted_identity": getattr(args, "min_weighted_identity", None),
    }
    
    # Set thresholds
    for region in args.regions.split(","):
        region = region.strip()
        threshold = getattr(args, f"min_id_{region}", None)
        if threshold is not None:
            cfg[f"min_id_{region}"] = threshold
    
    # Create output structure
    thresholds = {f"min_id_{r.strip()}": cfg.get(f"min_id_{r.strip()}") 
                  for r in args.regions.split(",")}
    csv_path, output_folder, csv_folder = create_output_structure(
        args.outdir, db_name, args.regions, thresholds, args.tag
    )
    
    start_time = dt.datetime.now()
    
    # Create results file path
    results_path = csv_path.replace('.csv', '_results.csv')
    
    # Write CSV headers
    with open(csv_path, "w", encoding="utf-8") as fh:
        write_metadata_block(fh, args, query_cdrs, start_time, cfg)
        fh.write("shard,total_seqs,hits,skipped,length_filtered,seq_length_filtered,"
                 "cdr_length_filtered,max_id,min_id,q75_id,best_seq_preview,worst_seq_preview,time_s\n")
    
    # Write results file header with metadata columns
    with open(results_path, "w", encoding="utf-8") as fh:
        fh.write("# Detailed hit results with metadata\n")
        fh.write(f"# Query CDR1: {query_cdrs.get('cdr1', 'N/A')}\n")
        fh.write(f"# Query CDR2: {query_cdrs.get('cdr2', 'N/A')}\n")
        fh.write(f"# Query CDR3: {query_cdrs.get('cdr3', 'N/A')}\n")
        fh.write(f"# Searched regions: {args.regions}\n")
        fh.write("#\n")
        
        # Extended header with metadata columns
        header_parts = [
            "shard", "index", "id", "source", "avg_identity", "full_identity",
            "cdr1_identity", "cdr1_seq",
            "cdr2_identity", "cdr2_seq", 
            "cdr3_identity", "cdr3_seq",
            "full_sequence",
            # Metadata columns
            "target", "patent_id", "patent_title", "organism", "pdb_id", "reference"
        ]
        fh.write(",".join(header_parts) + "\n")
    
    # Run full scan
    print("="*80)
    print("RUNNING FULL SCAN")
    print("="*80)
    print(f"Searching {len(npz_files)} shards for {args.regions} matches")
    
    # Build threshold display string
    threshold_parts = []
    for r in args.regions.split(','):
        r = r.strip()
        thresh = cfg.get(f"min_id_{r}", 0)
        if thresh:
            threshold_parts.append(f"{r}≥{thresh:.0%}")
    if threshold_parts:
        print(f"Thresholds: {', '.join(threshold_parts)}")
    print("="*80)
    
    total_hits = 0
    total_annotated_hits = 0
    total_filtered_length = 0
    total_filtered_seq_length = 0
    total_filtered_cdr_length = 0
    total_seqs_processed = 0
    global_best_id = 0.0
    global_worst_id = 1.0
    global_best_seq = ""
    global_worst_seq = ""
    
    scan_start = time.time()
    
    for i, npz_path in enumerate(npz_files, start=1):
        stats = process_npz_file(npz_path, cfg, query_cdrs, i, len(npz_files))
        
        # Update totals
        total_hits += stats["n_hits"]
        total_seqs_processed += stats["n_total"]
        total_filtered_length += stats["n_filtered_length"]
        total_filtered_seq_length += stats["n_filtered_seq_length"]
        total_filtered_cdr_length += stats["n_filtered_cdr_length"]
        
        # Count annotated hits
        for hit in stats.get("hits", []):
            if hit.get("target"):
                total_annotated_hits += 1
        
        # Track global best/worst
        if stats["max_id"] > global_best_id:
            global_best_id = stats["max_id"]
            global_best_seq = stats["best_seq"]
        if stats["n_hits"] > 0 and stats["min_id"] < global_worst_id:
            global_worst_id = stats["min_id"]
            global_worst_seq = stats["worst_seq"]
        
        # Print running totals
        elapsed = time.time() - scan_start
        rate = total_seqs_processed / elapsed if elapsed > 0 else 0
        print(f"      Running total: {total_hits:,} hits from {total_seqs_processed:,} sequences ({rate:,.0f} seq/s)")
        
        # Append summary to CSV
        with open(csv_path, "a", encoding="utf-8") as fh:
            fh.write(
                f"{stats['shard']},{stats['n_total']},{stats['n_hits']},{stats['n_skipped']},"
                f"{stats['n_filtered_length']},{stats['n_filtered_seq_length']},"
                f"{stats['n_filtered_cdr_length']},{stats['max_id']:.6f},{stats['min_id']:.6f},"
                f"{stats['q75']:.6f},\"{stats['best_seq']}\",\"{stats['worst_seq']}\",{stats['time_s']:.3f}\n"
            )
            fh.flush()
        
        # Append hit details to results file
        if stats.get("hits"):
            with open(results_path, "a", encoding="utf-8") as fh:
                for hit in stats["hits"]:
                    # Escape quotes in strings
                    def escape_csv(val):
                        if val is None:
                            return ""
                        val = str(val)
                        if '"' in val or ',' in val:
                            return '"' + val.replace('"', '""') + '"'
                        return val
                    
                    row_parts = [
                        escape_csv(hit.get("shard", "")),
                        str(hit.get("index", 0)),
                        escape_csv(hit.get("id", "")),
                        escape_csv(hit.get("source", "")),
                        f"{hit.get('avg_identity', 0.0):.6f}",
                        f"{hit.get('full_identity', 0.0):.6f}",
                        f"{hit.get('cdr1_identity', 0.0):.6f}",
                        escape_csv(hit.get("cdr1_seq", "")),
                        f"{hit.get('cdr2_identity', 0.0):.6f}",
                        escape_csv(hit.get("cdr2_seq", "")),
                        f"{hit.get('cdr3_identity', 0.0):.6f}",
                        escape_csv(hit.get("cdr3_seq", "")),
                        escape_csv(hit.get("full_sequence", "")),
                        # Metadata
                        escape_csv(hit.get("target", "")),
                        escape_csv(hit.get("patent_id", "")),
                        escape_csv(hit.get("patent_title", "")),
                        escape_csv(hit.get("organism", "")),
                        escape_csv(hit.get("pdb_id", "")),
                        escape_csv(hit.get("reference", "")),
                    ]
                    fh.write(",".join(row_parts) + "\n")
                fh.flush()
    
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
        fh.write(f"# Annotated hits (with target info): {total_annotated_hits:,}\n")
        fh.write(f"# Best identity found: {global_best_id:.6f}\n")
        fh.write(f"# Worst identity found: {global_worst_id:.6f}\n")
        fh.write("# ============================================\n")
    
    # Print summary
    print(f"\n✅ Done – all {len(npz_files)} shards processed.")
    print(f"\n{'='*80}")
    print("RUN SUMMARY")
    print(f"{'='*80}")
    print(f"Completed: {end_time.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Elapsed time: {elapsed:.2f} seconds ({elapsed/60:.2f} minutes)")
    print(f"Total hits: {total_hits:,}")
    print(f"Annotated hits (with target info): {total_annotated_hits:,}")
    print(f"Best identity: {global_best_id:.2%}")
    
    # Export to Excel
    excel_path = None
    if total_hits > 0:
        print("\n📊 Creating Excel file...")
        excel_path = export_to_excel(csv_path, results_path, output_folder, query_cdrs, 
                                     regions=regions_list)
    
    # Final output summary
    print(f"\n{'='*80}")
    print("FILES SAVED")
    print(f"{'='*80}")
    
    print(f"\n📁 Output folder: {output_folder}")
    
    if excel_path and os.path.exists(excel_path):
        print(f"\n📊 Excel file:")
        print(f"   {os.path.basename(excel_path)}")
        print(f"   • Summary sheet")
        print(f"   • Results sheet ({total_hits:,} rows)")
        if total_annotated_hits > 0:
            print(f"   • Annotated Hits sheet ({total_annotated_hits:,} with target info)")
        print(f"   • Query Info sheet")
    
    print(f"\n📁 CSV files (csvs/ subfolder):")
    print(f"   1. {os.path.basename(csv_path)}")
    print(f"   2. {os.path.basename(results_path)}")
    
    print(f"\n✅ All {total_hits:,} sequences saved!")
    
    # Show sample annotated hits
    if total_annotated_hits > 0:
        print(f"\n{'='*80}")
        print("SAMPLE ANNOTATED HITS (with target info)")
        print(f"{'='*80}")
        
        # Read back some annotated hits
        try:
            results_df = pd.read_csv(results_path, comment='#')
            annotated = results_df[results_df['target'].notna() & (results_df['target'] != '')]
            
            for idx, row in annotated.head(5).iterrows():
                print(f"\n  ID: {row['id']}")
                print(f"  Source: {row['source']}")
                print(f"  Identity: {row['avg_identity']:.2%}")
                print(f"  CDR3: {row['cdr3_seq']}")
                print(f"  Target: {row['target']}")
                if row.get('patent_id'):
                    print(f"  Patent: {row['patent_id']}")
        except Exception as e:
            pass
    
    print("="*80 + "\n")


if __name__ == "__main__":
    main()
