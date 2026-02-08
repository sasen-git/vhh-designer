#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
NPZ FULLSCAN v5 - Complete Feature Set

NEW in v5:
- Sequence length pre-filtering for full scans (automatic speedup)
- CDR length pre-filtering for full scans (automatic when using full mode)
- Framework extraction and searching
- Position-specific amino acid constraints with similarity groups

Features:
- Live dual progress bars
- Direct CDR/FR extraction from numbering arrays
- Multiple filtering strategies for speed
- Position-specific matching constraints
- Persistent per-shard stats
- Per-run CSV summary

Usage examples:

# Full sequence with automatic optimizations:
python3 npz_fullscan_v5.py \
  --db-root "/path/to/Heavy/Human" \
  --query-seq "EIQLQQ..." \
  --regions full \
  --min-id-full 0.30 \
  --outdir "./results" \
  --tag "full_optimized"

# CDRs with position constraints:
python3 npz_fullscan_v5.py \
  --db-root "/path/to/Heavy/Camel" \
  --query-seq "EIQLQQ..." \
  --regions cdr1,cdr2,cdr3 \
  --min-id-cdr1 0.35 --min-id-cdr2 0.35 --min-id-cdr3 0.40 \
  --position-constraint "cdr3:4:H:similar" \
  --position-constraint "cdr3:7:C:exact" \
  --outdir "./results" \
  --tag "cdr_with_constraints"

# Framework search:
python3 npz_fullscan_v5.py \
  --db-root "/path/to/Heavy/Human" \
  --query-seq "EIQLQQ..." \
  --regions frameworks \
  --min-id-frameworks 0.35 \
  --outdir "./results" \
  --tag "framework_search"
"""

import os
import re
import time
import glob
import signal
import argparse
import datetime as dt
from typing import Dict, List, Tuple, Optional

import numpy as np
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

class _Timeout(Exception):
    pass

def _alarm_handler(signum, frame):
    raise _Timeout()


# ============================================================================
# AMINO ACID EXTRACTION
# ============================================================================

def _ints_to_aa(int_row):
    """Convert integer array to amino acid string."""
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


def extract_cdr1_from_npz(int_row):
    """Extract CDR1 from NPZ array."""
    return _ints_to_aa(int_row[39:46])


def extract_cdr2_from_npz(int_row):
    """Extract CDR2 from NPZ array."""
    part1 = _ints_to_aa(int_row[66:75])
    part2 = _ints_to_aa(int_row[85:91])
    return part1 + part2


def extract_cdr3_from_npz(int_row):
    """Extract CDR3 from NPZ array."""
    cdr3_region = int_row[150:190]
    cdr3_aa = _ints_to_aa(cdr3_region)
    
    if cdr3_aa and cdr3_aa[0] == 'C':
        w_pos = cdr3_aa.find('W')
        if w_pos > 0:
            return cdr3_aa[:w_pos]
    return cdr3_aa


def extract_whole_from_npz(int_row):
    """Extract whole variable domain from NPZ array."""
    return _ints_to_aa(int_row[16:190])


def extract_frameworks_from_npz(int_row):
    """Extract all frameworks concatenated (FR1+FR2+FR3+FR4)."""
    fr1 = _ints_to_aa(int_row[16:39])   # FR1: IMGT 1-26
    fr2 = _ints_to_aa(int_row[46:66])   # FR2: IMGT 39-55
    fr3 = _ints_to_aa(int_row[75:150])  # FR3: IMGT 66-104
    fr4 = _ints_to_aa(int_row[165:190]) # FR4: IMGT 118-128
    return fr1 + fr2 + fr3 + fr4


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


# ============================================================================
# ALIGNMENT WITH POSITION TRACKING
# ============================================================================

def needleman_wunsch_align(seq1: str, seq2: str, match=1, mismatch=-1, gap=-1) -> Tuple[str, str]:
    """
    Needleman-Wunsch global alignment returning aligned strings.
    Returns tuple of (aligned_seq1, aligned_seq2) with gaps inserted.
    """
    if not seq1 or not seq2:
        return "", ""
    
    m, n = len(seq1), len(seq2)
    
    # Initialize scoring matrix
    score = np.zeros((m + 1, n + 1))
    
    # Initialize first row and column
    for i in range(m + 1):
        score[i][0] = gap * i
    for j in range(n + 1):
        score[0][j] = gap * j
    
    # Fill scoring matrix
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match_score = score[i-1][j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch)
            delete_score = score[i-1][j] + gap
            insert_score = score[i][j-1] + gap
            score[i][j] = max(match_score, delete_score, insert_score)
    
    # Traceback
    align1, align2 = "", ""
    i, j = m, n
    
    while i > 0 or j > 0:
        if i > 0 and j > 0 and score[i][j] == score[i-1][j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch):
            align1 = seq1[i-1] + align1
            align2 = seq2[j-1] + align2
            i -= 1
            j -= 1
        elif i > 0 and score[i][j] == score[i-1][j] + gap:
            align1 = seq1[i-1] + align1
            align2 = "-" + align2
            i -= 1
        else:
            align1 = "-" + align1
            align2 = seq2[j-1] + align2
            j -= 1
    
    return align1, align2


def alignment_identity(a: str, b: str) -> float:
    """Calculate sequence identity using edit distance."""
    if not a or not b:
        return 0.0
    max_len = max(len(a), len(b))
    if max_len == 0:
        return 0.0
    
    # Use fast edit distance
    dist = edit_distance(a, b)
    identity = 1.0 - (dist / max_len)
    return max(0.0, identity)


def edit_distance(s1: str, s2: str) -> int:
    """Calculate Levenshtein edit distance."""
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


# ============================================================================
# POSITION CONSTRAINT PARSING AND CHECKING
# ============================================================================

class PositionConstraint:
    """Represents a position-specific amino acid constraint."""
    
    def __init__(self, region: str, position: int, amino_acids: List[str], use_similar: bool):
        self.region = region
        self.position = position  # 1-based position in the aligned sequence
        self.amino_acids = amino_acids
        self.use_similar = use_similar
    
    def check(self, aligned_query: str, aligned_db: str) -> bool:
        """
        Check if constraint is satisfied.
        aligned_query and aligned_db are the aligned strings (with gaps).
        Position is 1-based in the original query sequence (not counting gaps).
        """
        # Convert 1-based position to 0-based and account for gaps
        query_pos_in_alignment = self._get_alignment_position(aligned_query, self.position)
        
        if query_pos_in_alignment is None or query_pos_in_alignment >= len(aligned_db):
            return False
        
        db_aa = aligned_db[query_pos_in_alignment]
        
        # Skip if gap in database sequence
        if db_aa == '-':
            return False
        
        # Check if amino acid matches
        if self.use_similar:
            # Check similarity groups
            for allowed_aa in self.amino_acids:
                similar_group = SIMILAR_GROUPS.get(allowed_aa, [allowed_aa])
                if db_aa in similar_group:
                    return True
            return False
        else:
            # Exact match required
            return db_aa in self.amino_acids
    
    @staticmethod
    def _get_alignment_position(aligned_seq: str, query_position: int) -> Optional[int]:
        """
        Convert 1-based position in original sequence to position in aligned sequence.
        query_position is 1-based position in the original (ungapped) sequence.
        Returns 0-based index in the aligned sequence, or None if out of range.
        """
        non_gap_count = 0
        for i, char in enumerate(aligned_seq):
            if char != '-':
                non_gap_count += 1
                if non_gap_count == query_position:
                    return i
        return None
    
    def __str__(self):
        similarity = "similar" if self.use_similar else "exact"
        aas = ",".join(self.amino_acids)
        return f"{self.region}:pos{self.position}:{aas}:{similarity}"


def parse_position_constraint(constraint_str: str) -> PositionConstraint:
    """
    Parse position constraint string.
    Format: "region:position:amino_acids:mode"
    Examples:
      - "cdr3:4:H:similar"
      - "cdr3:4:H,K,R:exact"
      - "cdr1:2:S:similar"
    """
    parts = constraint_str.split(':')
    if len(parts) != 4:
        raise ValueError(f"Invalid constraint format: {constraint_str}. Expected 'region:position:amino_acids:mode'")
    
    region = parts[0].lower()
    if region not in ['cdr1', 'cdr2', 'cdr3', 'full', 'frameworks']:
        raise ValueError(f"Invalid region: {region}. Must be cdr1, cdr2, cdr3, full, or frameworks")
    
    try:
        position = int(parts[1])
        if position < 1:
            raise ValueError("Position must be >= 1")
    except ValueError:
        raise ValueError(f"Invalid position: {parts[1]}. Must be a positive integer")
    
    amino_acids = [aa.strip().upper() for aa in parts[2].split(',')]
    for aa in amino_acids:
        if aa not in 'ACDEFGHIKLMNPQRSTVWY':
            raise ValueError(f"Invalid amino acid: {aa}")
    
    mode = parts[3].lower()
    if mode not in ['similar', 'exact']:
        raise ValueError(f"Invalid mode: {mode}. Must be 'similar' or 'exact'")
    
    use_similar = (mode == 'similar')
    
    return PositionConstraint(region, position, amino_acids, use_similar)


def check_position_constraints(query_seq: str, db_seq: str, constraints: List[PositionConstraint]) -> bool:
    """
    Check if all position constraints are satisfied.
    Returns True if all constraints pass, False otherwise.
    """
    if not constraints:
        return True
    
    # Align sequences
    aligned_query, aligned_db = needleman_wunsch_align(query_seq, db_seq)
    
    # Check all constraints (AND logic)
    for constraint in constraints:
        if not constraint.check(aligned_query, aligned_db):
            return False
    
    return True



# ============================================================================
# ANARCI-BASED CDR EXTRACTION (FOR QUERY ONLY)
# ============================================================================

def extract_cdrs_from_query_sequence(seq: str, scheme: str = "imgt", timeout_s: int = 5) -> Dict[str, str]:
    """Extract IMGT CDRs from query using ANARCI."""
    seq_clean = seq.upper().replace("-", "").replace(".", "")
    if len(seq_clean) < 70 or any(c not in "ACDEFGHIKLMNPQRSTVWY" for c in seq_clean):
        return {"cdr1": "", "cdr2": "", "cdr3": heuristic_cdrh3(seq_clean)}

    if not _ANARCI_AVAILABLE:
        return {"cdr1": "", "cdr2": "", "cdr3": heuristic_cdrh3(seq_clean)}

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

    try:
        numbering = res[1][0][0][0]
    except (TypeError, IndexError, KeyError):
        return {"cdr1": "", "cdr2": "", "cdr3": heuristic_cdrh3(seq_clean)}

    def pos_val(pos_tuple):
        pos = pos_tuple[0]
        ins = pos_tuple[1] if len(pos_tuple) > 1 else ""
        if isinstance(ins, str) and ins.strip():
            return float(f"{pos}.{ord(ins.upper())-64:02d}")
        return float(pos)

    def grab_region(start, end):
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
        "cdr1": grab_region(27.0, 38.9),
        "cdr2": grab_region(56.0, 65.9),
        "cdr3": grab_region(105.0, 117.9),
    }


# ============================================================================
# PER-SHARD SCAN FUNCTION WITH ALL OPTIMIZATIONS
# ============================================================================

def process_npz_file(npz_path: str, cfg: Dict, query_cdrs: Dict, 
                      shard_idx: int, n_shards: int) -> Dict:
    """
    Process one NPZ shard with all optimizations:
    - Sequence length pre-filter for full scans
    - CDR length pre-filter for full scans
    - CDR/Framework/Full identity filtering
    - Position-specific constraints
    """
    start = time.time()
    regions = {r.strip().lower() for r in str(cfg["regions"]).split(",") if r.strip()}
    only_full = regions == {"full"}
    only_frameworks = regions == {"frameworks"}

    # Load numberings from NPZ
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
            "n_filtered_length": 0,
            "n_filtered_seq_length": 0,
            "n_filtered_cdr_length": 0,
            "n_filtered_position": 0,
            "max_id": 0.0,
            "q75": 0.0,
            "time_s": 0.0,
        }
    
    total_seqs = int(len(numberings))

    # Per-shard mini bar
    inner_bar = tqdm(total=total_seqs, ncols=80, position=1, leave=False)
    inner_bar.set_description(f"Shard {shard_idx:03d}/{n_shards}")

    hits_id_values: List[float] = []
    n_hits = 0
    n_skipped = 0
    n_filtered_length = 0
    n_filtered_seq_length = 0  # NEW: Full sequence length filter
    n_filtered_cdr_length = 0  # NEW: CDR length pre-filter
    n_filtered_position = 0    # NEW: Position constraint filter

    # Precompute query parts
    q_full = query_cdrs.get("full", "")
    q_c1 = query_cdrs.get("cdr1", "")
    q_c2 = query_cdrs.get("cdr2", "")
    q_c3 = query_cdrs.get("cdr3", "")
    q_frameworks = query_cdrs.get("frameworks", "")

    # Identity thresholds
    min_id_full = float(cfg.get("min_id_full", 0.0))
    min_id_c1 = float(cfg.get("min_id_cdr1", 0.0))
    min_id_c2 = float(cfg.get("min_id_cdr2", 0.0))
    min_id_c3 = float(cfg.get("min_id_cdr3", 0.0))
    min_id_frameworks = float(cfg.get("min_id_frameworks", 0.0))

    # Length filtering parameters
    use_length_filter = cfg.get("use_length_filter", False)
    length_windows = cfg.get("length_windows", {})
    query_lengths = cfg.get("query_lengths", {})

    # Position constraints
    position_constraints_by_region = cfg.get("position_constraints_by_region", {})

    # Optimization parameters for full scan
    seq_length_window = 10  # Skip if >10aa different (automatic optimization)
    query_full_len = len(q_full)

    # Process sequences
    for i in range(total_seqs):
        try:
            numbering_row = numberings[i]
        except Exception:
            n_skipped += 1
            inner_bar.update(1)
            continue

        # OPTIMIZATION 1: Sequence length pre-filter for full scans
        if only_full and query_full_len > 0:
            db_full = extract_whole_from_npz(numbering_row)
            db_full_len = len(db_full)
            
            # Skip if length difference > threshold
            if abs(db_full_len - query_full_len) > seq_length_window:
                n_filtered_seq_length += 1
                inner_bar.update(1)
                continue
            
            # Skip very short sequences
            if db_full_len < 70:
                n_skipped += 1
                inner_bar.update(1)
                continue

        # OPTIMIZATION 2: CDR length pre-filter for full scans
        if only_full and use_length_filter:
            # Extract CDRs for length check (cheap operation)
            db_c1_quick = extract_cdr1_from_npz(numbering_row)
            db_c2_quick = extract_cdr2_from_npz(numbering_row)
            db_c3_quick = extract_cdr3_from_npz(numbering_row)
            
            passes_cdr_length = True
            if 'cdr1' in length_windows:
                query_len = query_lengths.get('cdr1', 0)
                window = length_windows['cdr1']
                if not (query_len - window <= len(db_c1_quick) <= query_len + window):
                    passes_cdr_length = False
            
            if passes_cdr_length and 'cdr2' in length_windows:
                query_len = query_lengths.get('cdr2', 0)
                window = length_windows['cdr2']
                if not (query_len - window <= len(db_c2_quick) <= query_len + window):
                    passes_cdr_length = False
            
            if passes_cdr_length and 'cdr3' in length_windows:
                query_len = query_lengths.get('cdr3', 0)
                window = length_windows['cdr3']
                if not (query_len - window <= len(db_c3_quick) <= query_len + window):
                    passes_cdr_length = False
            
            if not passes_cdr_length:
                n_filtered_cdr_length += 1
                inner_bar.update(1)
                continue

        # Now do full processing based on mode
        if only_full:
            # Full sequence mode (already extracted above if not filtered)
            if 'db_full' not in locals():
                db_full = extract_whole_from_npz(numbering_row)
            
            # Check position constraints for full sequence
            if 'full' in position_constraints_by_region:
                if not check_position_constraints(q_full, db_full, position_constraints_by_region['full']):
                    n_filtered_position += 1
                    inner_bar.update(1)
                    continue
            
            ident = alignment_identity(q_full, db_full)
            if ident >= min_id_full:
                n_hits += 1
                hits_id_values.append(ident)

        elif only_frameworks:
            # Frameworks mode
            db_frameworks = extract_frameworks_from_npz(numbering_row)
            
            if len(db_frameworks) < 20:  # Skip very short
                n_skipped += 1
                inner_bar.update(1)
                continue
            
            # Check position constraints
            if 'frameworks' in position_constraints_by_region:
                if not check_position_constraints(q_frameworks, db_frameworks, position_constraints_by_region['frameworks']):
                    n_filtered_position += 1
                    inner_bar.update(1)
                    continue
            
            ident = alignment_identity(q_frameworks, db_frameworks)
            if ident >= min_id_frameworks:
                n_hits += 1
                hits_id_values.append(ident)

        else:
            # CDR mode (or mixed)
            db_c1 = extract_cdr1_from_npz(numbering_row)
            db_c2 = extract_cdr2_from_npz(numbering_row)
            db_c3 = extract_cdr3_from_npz(numbering_row)
            
            # CDR3 fallback
            if not db_c3 or not db_c3.startswith('C'):
                full = extract_whole_from_npz(numbering_row)
                db_c3 = heuristic_cdrh3(full)
            
            if not any([db_c1, db_c2, db_c3]):
                n_skipped += 1
                inner_bar.update(1)
                continue
            
            # Length filtering
            if use_length_filter:
                passes_filter = True
                
                if "cdr1" in regions and "cdr1" in length_windows:
                    query_len = query_lengths.get("cdr1", 0)
                    window = length_windows["cdr1"]
                    if not (query_len - window <= len(db_c1) <= query_len + window):
                        passes_filter = False
                
                if passes_filter and "cdr2" in regions and "cdr2" in length_windows:
                    query_len = query_lengths.get("cdr2", 0)
                    window = length_windows["cdr2"]
                    if not (query_len - window <= len(db_c2) <= query_len + window):
                        passes_filter = False
                
                if passes_filter and "cdr3" in regions and "cdr3" in length_windows:
                    query_len = query_lengths.get("cdr3", 0)
                    window = length_windows["cdr3"]
                    if not (query_len - window <= len(db_c3) <= query_len + window):
                        passes_filter = False
                
                if not passes_filter:
                    n_filtered_length += 1
                    inner_bar.update(1)
                    continue
            
            # Position constraints
            position_failed = False
            if 'cdr1' in position_constraints_by_region:
                if not check_position_constraints(q_c1, db_c1, position_constraints_by_region['cdr1']):
                    position_failed = True
            if not position_failed and 'cdr2' in position_constraints_by_region:
                if not check_position_constraints(q_c2, db_c2, position_constraints_by_region['cdr2']):
                    position_failed = True
            if not position_failed and 'cdr3' in position_constraints_by_region:
                if not check_position_constraints(q_c3, db_c3, position_constraints_by_region['cdr3']):
                    position_failed = True
            
            if position_failed:
                n_filtered_position += 1
                inner_bar.update(1)
                continue
            
            # Identity filtering
            keep = True
            if "cdr1" in regions:
                if alignment_identity(q_c1, db_c1) < min_id_c1:
                    keep = False
            if keep and "cdr2" in regions:
                if alignment_identity(q_c2, db_c2) < min_id_c2:
                    keep = False
            if keep and "cdr3" in regions:
                ident3 = alignment_identity(q_c3, db_c3)
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

    # Compute stats
    if hits_id_values:
        ids = np.array(hits_id_values, dtype=float)
        max_id = float(np.max(ids))
        q75 = float(np.percentile(ids, 75))
    else:
        max_id = 0.0
        q75 = 0.0

    # Persistent stats line
    pct = (n_hits / total_seqs * 100.0) if total_seqs else 0.0
    stats_str = f"Shard {shard_idx:03d}/{n_shards} — {total_seqs:,} sequences\n"
    stats_str += f"✔ hits={n_hits:,} ({pct:.1f}%)\tmax ID={max_id:.2f}\ttop 25% ≥ {q75:.2f}"
    stats_str += f"\tskipped={n_skipped:,}"
    
    if use_length_filter:
        stats_str += f"\tlen_filt={n_filtered_length:,}"
    if only_full and n_filtered_seq_length > 0:
        stats_str += f"\tseq_len_filt={n_filtered_seq_length:,}"
    if only_full and n_filtered_cdr_length > 0:
        stats_str += f"\tcdr_len_filt={n_filtered_cdr_length:,}"
    if n_filtered_position > 0:
        stats_str += f"\tpos_filt={n_filtered_position:,}"
    
    stats_str += f"\ttime={elapsed:.1f} s"
    print(stats_str)
    print("─" * 79)

    return {
        "shard": shard_idx,
        "n_total": total_seqs,
        "n_hits": n_hits,
        "n_skipped": n_skipped,
        "n_filtered_length": n_filtered_length,
        "n_filtered_seq_length": n_filtered_seq_length,
        "n_filtered_cdr_length": n_filtered_cdr_length,
        "n_filtered_position": n_filtered_position,
        "max_id": max_id,
        "q75": q75,
        "time_s": elapsed,
    }



# ============================================================================
# MAIN
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="NPZ FULLSCAN v5 — Complete feature set with optimizations",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:

# Full sequence search (automatically optimized):
python3 npz_fullscan_v5.py \\
  --db-root "/path/to/Heavy/Human" \\
  --query-seq "EIQLQQ..." \\
  --regions full \\
  --min-id-full 0.30 \\
  --outdir "./results" \\
  --tag "full_search"

# CDR search with position constraints:
python3 npz_fullscan_v5.py \\
  --db-root "/path/to/Heavy/Camel" \\
  --query-seq "EIQLQQ..." \\
  --regions cdr1,cdr2,cdr3 \\
  --min-id-cdr1 0.35 --min-id-cdr2 0.35 --min-id-cdr3 0.40 \\
  --position-constraint "cdr3:4:H:similar" \\
  --position-constraint "cdr3:7:C:exact" \\
  --outdir "./results" \\
  --tag "cdr_constrained"

# Framework search:
python3 npz_fullscan_v5.py \\
  --db-root "/path/to/Heavy/Human" \\
  --query-seq "EIQLQQ..." \\
  --regions frameworks \\
  --min-id-frameworks 0.35 \\
  --outdir "./results" \\
  --tag "framework_search"

# Position constraint syntax:
  --position-constraint "REGION:POSITION:AMINO_ACIDS:MODE"
  
  REGION: cdr1, cdr2, cdr3, full, or frameworks
  POSITION: 1-based position in the query sequence (e.g., 4)
  AMINO_ACIDS: Single AA or comma-separated list (e.g., H or H,K,R)
  MODE: "similar" or "exact"
  
  Examples:
    "cdr3:4:H:similar"     - Position 4 must be H, K, or R (positive charged)
    "cdr3:7:C:exact"       - Position 7 must be exactly C
    "cdr1:2:S,T:exact"     - Position 2 must be S or T
    "cdr2:5:F:similar"     - Position 5 must be F, Y, or W (aromatic)
        """
    )
    
    # Required arguments
    parser.add_argument("--db-root", required=True, help="Directory with .npz shards")
    parser.add_argument("--query-seq", required=True, help="Query AA sequence")
    
    # Region selection
    parser.add_argument("--regions", default="full", 
                        help="Comma list: full, frameworks, or any of cdr1,cdr2,cdr3")
    
    # Identity thresholds
    parser.add_argument("--numbering-scheme", dest="scheme", default="imgt")
    parser.add_argument("--min-id-full", type=float, default=0.0)
    parser.add_argument("--min-id-cdr1", type=float, default=0.0)
    parser.add_argument("--min-id-cdr2", type=float, default=0.0)
    parser.add_argument("--min-id-cdr3", type=float, default=0.0)
    parser.add_argument("--min-id-frameworks", type=float, default=0.0,
                        help="Identity threshold for framework regions")
    
    # Length filtering
    parser.add_argument("--cdr1-len-window", type=int, default=None,
                        help="CDR1 length filter: query_length ± window")
    parser.add_argument("--cdr2-len-window", type=int, default=None,
                        help="CDR2 length filter: query_length ± window")
    parser.add_argument("--cdr3-len-window", type=int, default=None,
                        help="CDR3 length filter: query_length ± window")
    
    # Position constraints (NEW!)
    parser.add_argument("--position-constraint", action='append', dest='position_constraints',
                        help="Position-specific AA constraint (can be used multiple times). "
                             "Format: 'region:position:amino_acids:mode' "
                             "Example: 'cdr3:4:H:similar' or 'cdr3:7:C:exact'")
    
    # Output
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--tag", default="run")
    
    args = parser.parse_args()

    print("=" * 80)
    print("NPZ FULLSCAN v5 - Complete Feature Set")
    print("=" * 80)
    print(f"Started:       {dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Database:      {args.db_root}")
    print(f"Regions:       {args.regions}")
    print(f"Numbering:     {args.scheme}")
    
    # Parse and display position constraints
    position_constraints = []
    if args.position_constraints:
        print("\nPosition Constraints:")
        for constraint_str in args.position_constraints:
            try:
                constraint = parse_position_constraint(constraint_str)
                position_constraints.append(constraint)
                print(f"  • {constraint}")
            except ValueError as e:
                print(f"  ❌ ERROR: {e}")
                return
    
    print("=" * 80)

    os.makedirs(args.outdir, exist_ok=True)

    # Collect shards
    npz_files = sorted(glob.glob(os.path.join(args.db_root, "*.npz")))
    print(f"Found {len(npz_files)} shards\n")

    if not npz_files:
        print("[ERROR] No NPZ files found!")
        return

    # Extract query CDRs via ANARCI once
    print("Extracting query CDRs...\n")
    q_full = re.sub(r"[^A-Za-z]", "", args.query_seq.upper())

    query_cdrs = {"full": q_full, "cdr1": "", "cdr2": "", "cdr3": "", "frameworks": ""}
    
    # Extract query CDRs if needed
    regions_selected = {r.strip().lower() for r in str(args.regions).split(",") if r.strip()}
    if regions_selected != {"full"}:
        qc = extract_cdrs_from_query_sequence(q_full, args.scheme, timeout_s=5)
        query_cdrs.update({k: qc.get(k, "") for k in ("cdr1", "cdr2", "cdr3")})
    
    # Extract frameworks if needed (approximate from full sequence)
    if 'frameworks' in regions_selected:
        # For query, we can estimate frameworks
        # In practice, this should use ANARCI numbering properly
        # For now, use full sequence as placeholder
        query_cdrs['frameworks'] = q_full

    print("Query sequence (cleaned):")
    print(q_full)
    print("\nQuery regions extracted:")
    print(f"  CDR-H1: {query_cdrs.get('cdr1','')} (len={len(query_cdrs.get('cdr1',''))})")
    print(f"  CDR-H2: {query_cdrs.get('cdr2','')} (len={len(query_cdrs.get('cdr2',''))})")
    print(f"  CDR-H3: {query_cdrs.get('cdr3','')} (len={len(query_cdrs.get('cdr3',''))})")
    print(f"  Full:   {len(q_full)} aa")
    if 'frameworks' in regions_selected:
        print(f"  Frameworks: {len(query_cdrs['frameworks'])} aa (combined)")
    print()

    # Show position constraints with actual query sequences
    if position_constraints:
        print("Position Constraint Details:")
        for constraint in position_constraints:
            region_seq = query_cdrs.get(constraint.region, "")
            if region_seq:
                if constraint.position <= len(region_seq):
                    query_aa = region_seq[constraint.position - 1]
                    mode_str = "similar" if constraint.use_similar else "exact"
                    allowed_str = ",".join(constraint.amino_acids)
                    
                    print(f"  • {constraint.region.upper()} position {constraint.position} (query has '{query_aa}')")
                    print(f"    Must be: {allowed_str} ({mode_str})")
                    
                    if constraint.use_similar:
                        # Show what "similar" means
                        all_allowed = set()
                        for aa in constraint.amino_acids:
                            all_allowed.update(SIMILAR_GROUPS.get(aa, [aa]))
                        print(f"    Accepts: {','.join(sorted(all_allowed))}")
                else:
                    print(f"  ⚠️  WARNING: {constraint.region.upper()} position {constraint.position} is out of range (sequence length: {len(region_seq)})")
        print()

    # Organize constraints by region for faster lookup
    position_constraints_by_region = {}
    for constraint in position_constraints:
        if constraint.region not in position_constraints_by_region:
            position_constraints_by_region[constraint.region] = []
        position_constraints_by_region[constraint.region].append(constraint)

    # Build config
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
        "position_constraints_by_region": position_constraints_by_region
    }
    
    # Check if any length windows specified
    if args.cdr1_len_window is not None:
        cfg["use_length_filter"] = True
        cfg["length_windows"]["cdr1"] = args.cdr1_len_window
    if args.cdr2_len_window is not None:
        cfg["use_length_filter"] = True
        cfg["length_windows"]["cdr2"] = args.cdr2_len_window
    if args.cdr3_len_window is not None:
        cfg["use_length_filter"] = True
        cfg["length_windows"]["cdr3"] = args.cdr3_len_window
    
    if cfg["use_length_filter"]:
        print("=" * 80)
        print("LENGTH FILTERING ENABLED")
        print("=" * 80)
        for cdr, window in cfg["length_windows"].items():
            query_len = cfg["query_lengths"][cdr]
            print(f"  {cdr.upper()}: Query length = {query_len}, Window = ±{window} aa")
            print(f"           Accepting lengths: {query_len - window} to {query_len + window}")
        print("=" * 80)
        print()

    # Show optimization info
    if regions_selected == {"full"}:
        print("=" * 80)
        print("FULL SEQUENCE MODE - Automatic Optimizations Enabled")
        print("=" * 80)
        print("  • Sequence length pre-filter: Skip if >10aa different")
        if cfg["use_length_filter"]:
            print("  • CDR length pre-filter: Check CDR lengths before full alignment")
        print("=" * 80)
        print()

    # Prepare per-run CSV path
    ts = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    base = re.sub(r"[^A-Za-z0-9_.-]", "_", args.tag or "run")
    csv_path = os.path.join(args.outdir, f"summary_stats_{base}_{ts}.csv")

    # Write CSV header
    with open(csv_path, "w", encoding="utf-8") as fh:
        fh.write("shard_id,total_seqs,hits,skipped,length_filtered,seq_length_filtered,"
                 "cdr_length_filtered,position_filtered,max_id,top25_cutoff,time_s\n")

    total_hits = 0
    total_filtered_length = 0
    total_filtered_seq_length = 0
    total_filtered_cdr_length = 0
    total_filtered_position = 0
    
    outer_bar = tqdm(total=len(npz_files), desc="Processing shards", ncols=100, position=0, leave=True)

    for i, npz_path in enumerate(npz_files, start=1):
        stats = process_npz_file(npz_path, cfg, {**query_cdrs}, i, len(npz_files))
        total_hits += stats["n_hits"]
        total_filtered_length += stats["n_filtered_length"]
        total_filtered_seq_length += stats["n_filtered_seq_length"]
        total_filtered_cdr_length += stats["n_filtered_cdr_length"]
        total_filtered_position += stats["n_filtered_position"]
        
        # Append to CSV
        with open(csv_path, "a", encoding="utf-8") as fh:
            fh.write(
                f"{stats['shard']},{stats['n_total']},{stats['n_hits']},{stats['n_skipped']},"
                f"{stats['n_filtered_length']},{stats['n_filtered_seq_length']},"
                f"{stats['n_filtered_cdr_length']},{stats['n_filtered_position']},"
                f"{stats['max_id']:.6f},{stats['q75']:.6f},{stats['time_s']:.3f}\n"
            )
        outer_bar.update(1)

    outer_bar.close()

    print(f"\n✅ Done — all {len(npz_files)} shards processed.")
    print(f"Total hits: {total_hits:,}")
    
    if cfg["use_length_filter"]:
        print(f"Total filtered by CDR length: {total_filtered_length:,}")
    if total_filtered_seq_length > 0:
        print(f"Total filtered by sequence length: {total_filtered_seq_length:,}")
    if total_filtered_cdr_length > 0:
        print(f"Total filtered by CDR length pre-filter: {total_filtered_cdr_length:,}")
    if total_filtered_position > 0:
        print(f"Total filtered by position constraints: {total_filtered_position:,}")
    
    print(f"Summary saved to: {csv_path}")


if __name__ == "__main__":
    main()
