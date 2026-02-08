#!/usr/bin/env python3
"""
VHH Database ANARCI Annotation Pipeline
========================================

Processes all VHH database shards with ANARCI to create a unified, fully-annotated database.

Features:
- Runs ANARCI for IMGT numbering (if needed)
- Extracts CDR1, CDR2, CDR3 and FR1, FR2, FR3, FR4
- Classifies VHH family (F_C2, Y_C2, VH_like, etc.)
- Preserves metadata (targets, source, patent info, etc.)
- Outputs both CSV and NPZ formats

Estimated time:
- If shards have CDRs already: ~30-60 min for 12M sequences (no ANARCI needed)
- If ANARCI needed: ~8-12 hours for 12M sequences (batch processing)

Usage:
    python annotate_all_shards.py \
        --input-dir data/databases/shards/ \
        --output-prefix data/databases/annotated/vhh_annotated_full \
        --batch-size 2000

Author: Claude (Anthropic)
Date: January 2026
"""

import os
import sys
import glob
import argparse
import warnings
import json
import gc
from typing import List, Dict, Tuple, Optional
from pathlib import Path
from datetime import datetime
from collections import Counter

import numpy as np
import pandas as pd
from tqdm import tqdm

# Try to import ANARCI
try:
    from anarci import run_anarci
    ANARCI_AVAILABLE = True
except ImportError:
    ANARCI_AVAILABLE = False
    print("Warning: ANARCI not available. Install with: pip install anarci")


# =============================================================================
# CONSTANTS
# =============================================================================

# IMGT region boundaries (inclusive)
IMGT_REGIONS = {
    'FR1': (1, 26),
    'CDR1': (27, 38),
    'FR2': (39, 55),      # Note: Actual FR2 starts at 41 for W anchor
    'CDR2': (56, 65),
    'FR3': (66, 104),
    'CDR3': (105, 117),
    'FR4': (118, 128),
}

# Hallmark positions for family classification (IMGT numbering)
# IMPORTANT: Do NOT derive these by indexing into the extracted FR2 string.
# We always compute hallmarks directly from the IMGT-numbered dict (imgt_dict).
HALLMARK_POSITIONS = {
    42: 'pos42',  # F or Y in classical VHH
    49: 'pos49',  # E or Q in classical VHH, G in humanized
    50: 'pos50',  # R in classical VHH, L in human VH
    52: 'pos52',  # Variable
}


# =============================================================================
# ANARCI FUNCTIONS
# =============================================================================

def pos_val(pos_tuple: Tuple) -> float:
    """Convert IMGT position with insertion code to a float value."""
    pos = pos_tuple[0]
    ins = pos_tuple[1] if len(pos_tuple) > 1 else " "
    
    if isinstance(ins, str) and ins.strip():
        ins_val = ord(ins.upper()) - 64
        return float(f"{pos}.{ins_val:02d}")
    return float(pos)


def extract_region(numbering: List, start: float, end: float) -> str:
    """Extract amino acids from a region defined by IMGT positions."""
    aa_list = []
    for entry in numbering:
        if not isinstance(entry, (tuple, list)) or len(entry) != 2:
            continue
        pos_tuple, residue = entry
        if residue in (None, '-', '.'):
            continue
        if not isinstance(pos_tuple, (tuple, list)):
            continue
        v = pos_val(pos_tuple)
        if start <= v <= end + 0.99:
            aa_list.append(residue)
    return ''.join(aa_list)


def extract_full_v_domain(numbering: List) -> str:
    """Extract full V-domain sequence from ANARCI numbering."""
    aa_list = []
    for entry in numbering:
        if not isinstance(entry, (tuple, list)) or len(entry) != 2:
            continue
        pos_tuple, residue = entry
        if residue in (None, '-', '.'):
            continue
        aa_list.append(residue)
    return ''.join(aa_list)


def numbering_to_imgt_dict(numbering: List) -> Dict[int, str]:
    """Convert ANARCI numbering to dict of IMGT position -> amino acid."""
    imgt_dict = {}
    for entry in numbering:
        if not isinstance(entry, (tuple, list)) or len(entry) != 2:
            continue
        pos_tuple, residue = entry
        if residue in (None, '-', '.'):
            continue
        if not isinstance(pos_tuple, (tuple, list)):
            continue
        
        pos = pos_tuple[0]
        ins = pos_tuple[1] if len(pos_tuple) > 1 else " "
        
        # Store with insertion code handling
        if isinstance(ins, str) and ins.strip():
            # Has insertion code - store as float key or string
            key = f"{pos}{ins}"
        else:
            key = int(pos)
        
        imgt_dict[key] = residue
    
    return imgt_dict


def run_anarci_batch(sequences: List[Tuple[str, str]], scheme: str = "imgt") -> List[Optional[List]]:
    """Run ANARCI on a batch of sequences."""
    if not ANARCI_AVAILABLE:
        return [None] * len(sequences)
    
    if not sequences:
        return []
    
    try:
        result = run_anarci(sequences, scheme=scheme, allowed_species=None)
        numberings_list = result[1]
        
        parsed_results = []
        for i, entry in enumerate(numberings_list):
            try:
                if entry is None or not entry:
                    parsed_results.append(None)
                    continue
                
                domain_data = entry[0]
                if domain_data is None:
                    parsed_results.append(None)
                    continue
                
                numbering = domain_data[0]
                if numbering and isinstance(numbering, list):
                    parsed_results.append(numbering)
                else:
                    parsed_results.append(None)
                    
            except (IndexError, TypeError):
                parsed_results.append(None)
        
        return parsed_results
        
    except Exception as e:
        warnings.warn(f"ANARCI batch failed: {e}")
        return [None] * len(sequences)


# =============================================================================
# REGION EXTRACTION (WITHOUT ANARCI)
# =============================================================================

def extract_regions_from_sequence(aa_v_full: str, cdr1: str, cdr2: str, cdr3: str) -> Dict:
    """
    Extract FR regions and build IMGT position dict from sequence + CDRs.
    Used when ANARCI is not available or CDRs are pre-extracted.
    
    FR2 anchors W at IMGT 41 (corrected mapping).
    """
    result = {
        'fr1': '', 'fr2': '', 'fr3': '', 'fr4': '',
        'cdr1': cdr1, 'cdr2': cdr2, 'cdr3': cdr3,
        'imgt_positions': {},
        'valid': False,
        'error': None,
    }
    
    if not aa_v_full or not cdr1 or not cdr2 or not cdr3:
        result['error'] = 'missing_data'
        return result
    
    try:
        # Find CDR positions in full sequence
        cdr1_start = aa_v_full.find(cdr1)
        if cdr1_start == -1:
            result['error'] = 'cdr1_not_found'
            return result
        cdr1_end = cdr1_start + len(cdr1)
        
        cdr2_start = aa_v_full.find(cdr2, cdr1_end)
        if cdr2_start == -1:
            result['error'] = 'cdr2_not_found'
            return result
        cdr2_end = cdr2_start + len(cdr2)
        
        cdr3_start = aa_v_full.find(cdr3, cdr2_end)
        if cdr3_start == -1:
            result['error'] = 'cdr3_not_found'
            return result
        cdr3_end = cdr3_start + len(cdr3)
        
        # Extract framework regions
        fr1 = aa_v_full[:cdr1_start]
        fr2 = aa_v_full[cdr1_end:cdr2_start]
        fr3 = aa_v_full[cdr2_end:cdr3_start]
        fr4 = aa_v_full[cdr3_end:]
        
        # Validate lengths (typical VHH ranges)
        if not (15 <= len(fr1) <= 30):
            result['error'] = f'fr1_length_{len(fr1)}'
            return result
        if not (10 <= len(fr2) <= 20):
            result['error'] = f'fr2_length_{len(fr2)}'
            return result
        if not (30 <= len(fr3) <= 50):
            result['error'] = f'fr3_length_{len(fr3)}'
            return result
        if not (5 <= len(fr4) <= 15):
            result['error'] = f'fr4_length_{len(fr4)}'
            return result
        
        result['fr1'] = fr1
        result['fr2'] = fr2
        result['fr3'] = fr3
        result['fr4'] = fr4
        result['valid'] = True
        
        # Build IMGT position dict
        # Standard IMGT numbering for VHH:
        # FR1: 1-26, CDR1: 27-38, FR2: 41-54 (W anchor at 41), CDR2: 56-65, FR3: 66-104, CDR3: 105-117, FR4: 118-128
        positions = {}
        
        # FR1 (IMGT 1-26)
        for i, aa in enumerate(fr1):
            imgt_pos = i + 1
            if imgt_pos <= 26:
                positions[imgt_pos] = aa
        
        # CDR1 (IMGT 27-38) - variable length
        cdr1_imgt_start = 27
        for i, aa in enumerate(cdr1):
            positions[cdr1_imgt_start + i] = aa
        
        # FR2 (IMGT 41-54) - W anchor at 41
        fr2_imgt_start = 41
        for i, aa in enumerate(fr2):
            imgt_pos = fr2_imgt_start + i
            if imgt_pos <= 54:
                positions[imgt_pos] = aa
        
        # CDR2 (IMGT 56-65)
        cdr2_imgt_start = 56
        for i, aa in enumerate(cdr2):
            positions[cdr2_imgt_start + i] = aa
        
        # FR3 (IMGT 66-104)
        fr3_imgt_start = 66
        for i, aa in enumerate(fr3):
            imgt_pos = fr3_imgt_start + i
            if imgt_pos <= 104:
                positions[imgt_pos] = aa
        
        # CDR3 (IMGT 105-117) - highly variable
        cdr3_imgt_start = 105
        for i, aa in enumerate(cdr3):
            positions[cdr3_imgt_start + i] = aa
        
        # FR4 (IMGT 118-128)
        fr4_imgt_start = 118
        for i, aa in enumerate(fr4):
            imgt_pos = fr4_imgt_start + i
            if imgt_pos <= 128:
                positions[imgt_pos] = aa
        
        result['imgt_positions'] = positions
        
    except Exception as e:
        result['error'] = str(e)
    
    return result


# =============================================================================
# FAMILY CLASSIFICATION
# =============================================================================

def classify_family(imgt_dict: Dict, full_seq: str) -> str:
    """
    Classify VHH family based on hallmark positions.
    
    Hallmarks (IMGT positions):
    - pos42: F (llama) or Y (alpaca) in classical VHH
    - pos49: E in classical VHH, G in humanized
    - pos50: R in classical VHH, L in human VH
    - pos52: varies
    """
    pos42 = imgt_dict.get(42, '-')
    pos49 = imgt_dict.get(49, '-')
    pos50 = imgt_dict.get(50, '-')
    pos52 = imgt_dict.get(52, '-')
    
    # Count cysteines
    n_cys = full_seq.count('C') if full_seq else 2
    
    # VH-like: pos50 = L (human-like)
    if pos50 == 'L':
        return 'VH_like'
    
    # Classical VHH criteria
    is_classical = (
        pos42 in ['F', 'Y'] and
        pos49 in ['E', 'Q'] and
        pos50 == 'R'
    )
    
    if is_classical:
        if pos42 == 'Y' and n_cys == 2:
            return 'Y_C2'
        elif pos42 == 'F' and n_cys == 2:
            return 'F_C2'
        elif pos42 == 'F' and n_cys == 4:
            return 'F_C4'
        elif pos42 == 'Y' and n_cys == 4:
            return 'Y_C4'
        else:
            return 'Classical_other'
    
    return 'Non_classical'


def get_hallmark_string(imgt_dict: Dict) -> str:
    """Get hallmark positions as string like 'FERG' or 'YERL'."""
    pos42 = imgt_dict.get(42, '-')
    pos49 = imgt_dict.get(49, '-')
    pos50 = imgt_dict.get(50, '-')
    pos52 = imgt_dict.get(52, '-')
    return f"{pos42}{pos49}{pos50}{pos52}"


# =============================================================================
# NPZ LOADING
# =============================================================================

def load_npz_shard(path: str) -> Tuple[pd.DataFrame, Dict]:
    """
    Load a single NPZ shard into a DataFrame.
    
    Handles various formats:
    - Standard: aa_v_full, cdr1, cdr2, cdr3, ids
    - With metadata: targets, source, patent_id, etc.
    """
    data = np.load(path, allow_pickle=True)
    
    df_dict = {}
    metadata = {}
    
    # Load all arrays
    for key in data.files:
        arr = data[key]
        
        # Handle metadata separately
        if key == 'metadata':
            try:
                metadata = arr.item() if hasattr(arr, 'item') else dict(arr)
            except:
                pass
            continue
        
        # Convert to list for DataFrame
        if hasattr(arr, 'tolist'):
            df_dict[key] = arr.tolist()
        else:
            df_dict[key] = list(arr)
    
    # Create DataFrame
    if df_dict:
        # Make sure all arrays have same length
        max_len = max(len(v) for v in df_dict.values())
        for k in df_dict:
            if len(df_dict[k]) < max_len:
                df_dict[k] = df_dict[k] + [None] * (max_len - len(df_dict[k]))
        
        df = pd.DataFrame(df_dict)
    else:
        df = pd.DataFrame()
    
    return df, metadata


# =============================================================================
# MAIN PROCESSING
# =============================================================================

def process_shard(
    shard_path: str,
    use_anarci: bool = False,
    batch_size: int = 2000,
    *,
    emit_imgt_fr3fr4: bool = True,
    emit_imgt_fr2: bool = False,
) -> pd.DataFrame:
    """
    Process a single NPZ shard.
    
    If shard has CDRs, derive FRs and IMGT positions without ANARCI.
    If shard only has sequences, run ANARCI.
    """
    shard_name = os.path.basename(shard_path)
    print(f"\n  Processing: {shard_name}")
    
    # Load shard
    df, metadata = load_npz_shard(shard_path)
    n_seqs = len(df)
    print(f"    Loaded {n_seqs:,} sequences")
    
    if n_seqs == 0:
        return pd.DataFrame()
    
    # Check what columns we have
    has_aa_v_full = 'aa_v_full' in df.columns
    has_cdrs = 'cdr1' in df.columns and 'cdr2' in df.columns and 'cdr3' in df.columns
    has_sequences = 'sequences' in df.columns
    
    print(f"    Has aa_v_full: {has_aa_v_full}, Has CDRs: {has_cdrs}")
    
    # Determine processing mode
    want_imgt_positions = emit_imgt_fr3fr4 or emit_imgt_fr2
    if has_aa_v_full and has_cdrs and not use_anarci and want_imgt_positions:
        # If we want per-position IMGT columns, we MUST run ANARCI to get true IMGT numbering.
        print("    Mode: CDRs present but per-position IMGT requested -> forcing ANARCI for correct numbering")
        use_anarci = True

    if has_aa_v_full and has_cdrs and not use_anarci:
        # Fast path: derive FRs from existing CDRs (no IMGT per-position columns)
        print(f"    Mode: Deriving FRs from CDRs (no ANARCI)")
        df = process_with_existing_cdrs(df)
    elif has_aa_v_full:
        # Run ANARCI to extract CDRs + true IMGT positions
        print(f"    Mode: Running ANARCI")
        df = process_with_anarci(df, batch_size, emit_imgt_fr3fr4=emit_imgt_fr3fr4, emit_imgt_fr2=emit_imgt_fr2)
    elif has_sequences:
        # sequences column instead of aa_v_full
        df['aa_v_full'] = df['sequences']
        print(f"    Mode: Running ANARCI on 'sequences' column")
        df = process_with_anarci(df, batch_size, emit_imgt_fr3fr4=emit_imgt_fr3fr4, emit_imgt_fr2=emit_imgt_fr2)
    else:
        print(f"    ERROR: No sequence column found!")
        return pd.DataFrame()
    
    # Add source shard info
    df['source_shard'] = shard_name
    
    # Keep track of original columns for metadata
    return df


def process_with_existing_cdrs(df: pd.DataFrame) -> pd.DataFrame:
    """Process sequences that already have CDRs extracted."""
    
    results = []
    
    for idx, row in tqdm(df.iterrows(), total=len(df), desc="    Extracting"):
        aa_v_full = str(row.get('aa_v_full', '')) if row.get('aa_v_full') else ''
        cdr1 = str(row.get('cdr1', '')) if row.get('cdr1') else ''
        cdr2 = str(row.get('cdr2', '')) if row.get('cdr2') else ''
        cdr3 = str(row.get('cdr3', '')) if row.get('cdr3') else ''
        
        # Extract regions
        regions = extract_regions_from_sequence(aa_v_full, cdr1, cdr2, cdr3)
        
        result = {
            'aa_v_full': aa_v_full,
            'cdr1': cdr1,
            'cdr2': cdr2,
            'cdr3': cdr3,
            'fr1': regions['fr1'],
            'fr2': regions['fr2'],
            'fr3': regions['fr3'],
            'fr4': regions['fr4'],
            'len_v': len(aa_v_full),
            'len_cdr1': len(cdr1),
            'len_cdr2': len(cdr2),
            'len_cdr3': len(cdr3),
            'len_fr1': len(regions['fr1']),
            'len_fr2': len(regions['fr2']),
            'len_fr3': len(regions['fr3']),
            'len_fr4': len(regions['fr4']),
            'valid': regions['valid'],
            'error': regions['error'],
        }
        
        # Classify family
        if regions['valid']:
            result['family'] = classify_family(regions['imgt_positions'], aa_v_full)
            result['hallmarks'] = get_hallmark_string(regions['imgt_positions'])
            
            # Store key IMGT positions
            for pos in [42, 49, 50, 52]:
                result[f'imgt_{pos}'] = regions['imgt_positions'].get(pos, '-')
        else:
            result['family'] = 'Unknown'
            result['hallmarks'] = '----'
            for pos in [42, 49, 50, 52]:
                result[f'imgt_{pos}'] = '-'
        
        # Preserve original metadata columns
        for col in df.columns:
            if col not in result and col not in ['aa_v_full', 'cdr1', 'cdr2', 'cdr3']:
                result[col] = row.get(col)
        
        results.append(result)
    
    return pd.DataFrame(results)


def process_with_anarci(df: pd.DataFrame, batch_size: int = 2000, *, emit_imgt_fr3fr4: bool = True, emit_imgt_fr2: bool = False) -> pd.DataFrame:
    """Process sequences using ANARCI for numbering."""
    
    if not ANARCI_AVAILABLE:
        print("    ERROR: ANARCI not available!")
        return pd.DataFrame()
    
    sequences = df['aa_v_full'].tolist()
    n_seqs = len(sequences)
    
    # Initialize results
    results = []
    
    # Process in batches
    for start_idx in tqdm(range(0, n_seqs, batch_size), desc="    ANARCI"):
        end_idx = min(start_idx + batch_size, n_seqs)
        batch_seqs = sequences[start_idx:end_idx]
        batch_df = df.iloc[start_idx:end_idx]
        
        # Prepare ANARCI input
        anarci_input = [(f"seq_{i}", seq) for i, seq in enumerate(batch_seqs) if seq]
        
        # Handle empty sequences
        if not anarci_input:
            for idx in range(start_idx, end_idx):
                results.append({'valid': False, 'error': 'empty_sequence'})
            continue
        
        # Run ANARCI
        numberings = run_anarci_batch(anarci_input, scheme="imgt")
        
        # Process results
        for i, numbering in enumerate(numberings):
            row = batch_df.iloc[i]
            aa_v_full = batch_seqs[i] if batch_seqs[i] else ''
            
            result = {
                'aa_v_full': aa_v_full,
                'len_v': len(aa_v_full),
            }
            
            if numbering is None:
                result['valid'] = False
                result['error'] = 'anarci_failed'
                result['cdr1'] = ''
                result['cdr2'] = ''
                result['cdr3'] = ''
                result['fr1'] = ''
                result['fr2'] = ''
                result['fr3'] = ''
                result['fr4'] = ''
                result['family'] = 'Unknown'
                result['hallmarks'] = '----'
            else:
                try:
                    # Extract regions
                    cdr1 = extract_region(numbering, 27.0, 38.99)
                    cdr2 = extract_region(numbering, 56.0, 65.99)
                    cdr3 = extract_region(numbering, 105.0, 117.99)
                    fr1 = extract_region(numbering, 1.0, 26.99)
                    fr2 = extract_region(numbering, 39.0, 55.99)
                    fr3 = extract_region(numbering, 66.0, 104.99)
                    fr4 = extract_region(numbering, 118.0, 128.99)
                    
                    result['cdr1'] = cdr1
                    result['cdr2'] = cdr2
                    result['cdr3'] = cdr3
                    result['fr1'] = fr1
                    result['fr2'] = fr2
                    result['fr3'] = fr3
                    result['fr4'] = fr4
                    result['len_cdr1'] = len(cdr1)
                    result['len_cdr2'] = len(cdr2)
                    result['len_cdr3'] = len(cdr3)
                    result['len_fr1'] = len(fr1)
                    result['len_fr2'] = len(fr2)
                    result['len_fr3'] = len(fr3)
                    result['len_fr4'] = len(fr4)
                    
                    # Get IMGT dict and classify
                    imgt_dict = numbering_to_imgt_dict(numbering)
                    result['family'] = classify_family(imgt_dict, aa_v_full)
                    result['hallmarks'] = get_hallmark_string(imgt_dict)
                    
                    for pos in [42, 49, 50, 52]:
                        result[f'imgt_{pos}'] = imgt_dict.get(pos, '-')


                    # Optional: emit per-position IMGT columns for max-coverage rule learning.
                    # These are "portable" positions (stable across sequences), unlike FR-substring indices.
                    if emit_imgt_fr2:
                        for pos in range(39, 56):   # FR2 span
                            result[f'imgt_{pos}'] = imgt_dict.get(pos, '-')
                    if emit_imgt_fr3fr4:
                        for pos in range(66, 105):  # FR3 span
                            result[f'imgt_{pos}'] = imgt_dict.get(pos, '-')
                        for pos in range(118, 129): # FR4 span
                            result[f'imgt_{pos}'] = imgt_dict.get(pos, '-')
                    
                    result['valid'] = True
                    result['error'] = None
                    
                except Exception as e:
                    result['valid'] = False
                    result['error'] = str(e)
            
            # Preserve original metadata
            for col in df.columns:
                if col not in result and col != 'aa_v_full':
                    result[col] = row.get(col)
            
            results.append(result)
    
    return pd.DataFrame(results)


def process_all_shards(
    input_dir: str,
    output_prefix: str,
    use_anarci: bool = False,
    batch_size: int = 2000,
    checkpoint_interval: int = 1000000,
    *,
    emit_imgt_fr3fr4: bool = True,
    emit_imgt_fr2: bool = False,
) -> None:
    """Process all NPZ shards in a directory."""
    
    print("=" * 70)
    print("VHH Database ANARCI Annotation Pipeline")
    print("=" * 70)
    print(f"Input directory: {input_dir}")
    print(f"Output prefix: {output_prefix}")
    print(f"Use ANARCI: {use_anarci}")
    print(f"Batch size: {batch_size}")
    print()
    
    # Find all NPZ files
    npz_files = sorted(glob.glob(os.path.join(input_dir, "*.npz")))
    npz_files = [f for f in npz_files if not f.endswith('Zone.Identifier')]
    
    print(f"Found {len(npz_files)} NPZ files:")
    for f in npz_files:
        size_mb = os.path.getsize(f) / 1024 / 1024
        print(f"  {os.path.basename(f)}: {size_mb:.1f} MB")
    print()
    
    # Create output directory
    output_dir = os.path.dirname(output_prefix)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    
    # Process each shard
    all_dfs = []
    total_processed = 0
    total_valid = 0
    family_counts = Counter()
    
    start_time = datetime.now()
    
    for shard_path in npz_files:
        df = process_shard(shard_path, use_anarci=use_anarci, batch_size=batch_size, emit_imgt_fr3fr4=emit_imgt_fr3fr4, emit_imgt_fr2=emit_imgt_fr2)
        
        if len(df) > 0:
            valid_mask = df['valid'] == True
            n_valid = valid_mask.sum()
            total_processed += len(df)
            total_valid += n_valid
            
            # Count families
            for fam in df.loc[valid_mask, 'family'].value_counts().to_dict().items():
                family_counts[fam[0]] += fam[1]
            
            print(f"    Valid: {n_valid:,} / {len(df):,} ({n_valid/len(df)*100:.1f}%)")
            
            all_dfs.append(df)
            
            # Free memory
            gc.collect()
    
    elapsed = (datetime.now() - start_time).total_seconds()
    
    print(f"\n{'='*70}")
    print("MERGING AND SAVING")
    print(f"{'='*70}")
    
    # Merge all DataFrames
    print(f"\nMerging {len(all_dfs)} shards...")
    df_all = pd.concat(all_dfs, ignore_index=True)
    
    print(f"Total sequences: {len(df_all):,}")
    print(f"Valid sequences: {total_valid:,} ({total_valid/len(df_all)*100:.1f}%)")
    
    # Generate unique IDs
    df_all = df_all.reset_index(drop=True)
    df_all['id'] = [f"vhh_{i:08d}" for i in range(len(df_all))]
    
    # Filter to valid only for output
    df_valid = df_all[df_all['valid'] == True].copy()
    
    print(f"\nFamily distribution:")
    for fam, count in sorted(family_counts.items(), key=lambda x: -x[1]):
        pct = 100 * count / total_valid if total_valid > 0 else 0
        print(f"  {fam}: {count:,} ({pct:.1f}%)")
    
    # Define output columns order
    core_cols = [
        'id', 'aa_v_full',
        'fr1', 'cdr1', 'fr2', 'cdr2', 'fr3', 'cdr3', 'fr4',
        'family', 'hallmarks',
        'imgt_42', 'imgt_49', 'imgt_50', 'imgt_52',

        # Optional: per-position IMGT columns (only present if emitted during ANARCI processing)
        # FR2 (39-55) is only emitted if emit_imgt_fr2=True
        *[f'imgt_{p}' for p in range(39, 56)],
        # FR3 (66-104) and FR4 (118-128) are emitted by default (emit_imgt_fr3fr4=True)
        *[f'imgt_{p}' for p in range(66, 105)],
        *[f'imgt_{p}' for p in range(118, 129)],

        'len_v', 'len_fr1', 'len_cdr1', 'len_fr2', 'len_cdr2', 'len_fr3', 'len_cdr3', 'len_fr4',
        'source_shard',
    ]
    
    # Add any metadata columns
    metadata_cols = [c for c in df_valid.columns if c not in core_cols and c not in ['valid', 'error']]
    output_cols = core_cols + metadata_cols
    
    # Filter to existing columns
    output_cols = [c for c in output_cols if c in df_valid.columns]
    
    # Save CSV
    csv_path = f"{output_prefix}.csv"
    print(f"\nSaving CSV: {csv_path}")
    df_valid[output_cols].to_csv(csv_path, index=False)
    csv_size = os.path.getsize(csv_path) / 1024 / 1024
    print(f"  Size: {csv_size:.1f} MB")
    
    # Save NPZ
    npz_path = f"{output_prefix}.npz"
    print(f"\nSaving NPZ: {npz_path}")
    
    npz_dict = {}
    for col in output_cols:
        arr = df_valid[col].values
        # Convert object arrays to proper types
        if arr.dtype == object:
            npz_dict[col] = np.array([str(x) if x is not None else '' for x in arr])
        else:
            npz_dict[col] = arr
    
    np.savez_compressed(npz_path, **npz_dict)
    npz_size = os.path.getsize(npz_path) / 1024 / 1024
    print(f"  Size: {npz_size:.1f} MB")
    
    # Save summary
    summary = {
        'total_input': int(len(df_all)),
        'total_valid': int(total_valid),
        'family_counts': {k: int(v) for k, v in family_counts.items()},
        'shards_processed': len(npz_files),
        'elapsed_seconds': float(elapsed),
        'output_csv': csv_path,
        'output_npz': npz_path,
        'timestamp': datetime.now().isoformat(),
    }
    
    summary_path = f"{output_prefix}_summary.json"
    with open(summary_path, 'w') as f:
        json.dump(summary, f, indent=2)
    print(f"\nSummary saved: {summary_path}")
    
    print(f"\n{'='*70}")
    print("COMPLETE")
    print(f"{'='*70}")
    print(f"Processed {total_processed:,} sequences in {elapsed/60:.1f} minutes")
    print(f"Valid output: {total_valid:,} sequences")
    print(f"Output files:")
    print(f"  {csv_path}")
    print(f"  {npz_path}")
    print(f"  {summary_path}")


# =============================================================================
# MAIN
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Annotate VHH database shards with ANARCI/IMGT numbering'
    )
    parser.add_argument('--input-dir', '-i', required=True,
                        help='Directory containing NPZ shards')
    parser.add_argument('--output-prefix', '-o', required=True,
                        help='Output file prefix (will create .csv and .npz)')
    parser.add_argument('--use-anarci', '-a', action='store_true',
                        help='Force ANARCI even if CDRs exist (slower but more accurate)')
    parser.add_argument('--batch-size', '-b', type=int, default=2000,
                        help='ANARCI batch size (default: 2000)')
    parser.add_argument('--checkpoint-interval', '-c', type=int, default=1000000,
                        help='Checkpoint interval (default: 1000000)')
    parser.add_argument('--no-imgt-fr3fr4', action='store_true',
                        help='Disable per-position IMGT columns for FR3/FR4 (saves space). Default: emit FR3/FR4.')
    parser.add_argument('--emit-imgt-fr2', action='store_true',
                        help='Also emit per-position IMGT columns for FR2 (IMGT39-55). Default: off (saves space).')

    
    args = parser.parse_args()
    
    process_all_shards(
        input_dir=args.input_dir,
        output_prefix=args.output_prefix,
        use_anarci=args.use_anarci,
        batch_size=args.batch_size,
        checkpoint_interval=args.checkpoint_interval,
        emit_imgt_fr3fr4=(not args.no_imgt_fr3fr4),
        emit_imgt_fr2=args.emit_imgt_fr2,
    )


if __name__ == '__main__':
    main()
