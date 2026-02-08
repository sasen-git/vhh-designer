#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Process VHH Sequences from GenBank/PLAbDab and Merge with OAS Camel Data

This script:
1. Loads VHH sequences from CSV (GenBank/PLAbDab format)
2. Extracts CDRs (uses provided CDRs or re-runs ANARCI)
3. Preserves antigen/target information
4. Optionally merges with OAS camel data
5. Creates unified CSV and NPZ outputs

Key feature: Preserves antigen info for display in search results!

Usage:
    # Process VHH CSV only
    python process_vhh_with_antigens.py \
        --vhh-csv /path/to/vhh_sequences.csv \
        --output-prefix vhh_annotated

    # Merge with OAS camel data
    python process_vhh_with_antigens.py \
        --vhh-csv /path/to/vhh_sequences.csv \
        --oas-csv /path/to/camel_vhh_clean_anarci_renumbered.csv \
        --output-prefix VHH_db_unified

Author: Claude (Anthropic)
Date: 2025
"""

import os
import sys
import ast
import argparse
import re
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from tqdm import tqdm

# Try to import ANARCI for re-numbering if needed
try:
    from anarci import run_anarci
    ANARCI_AVAILABLE = True
except ImportError:
    ANARCI_AVAILABLE = False


# =============================================================================
# CDR EXTRACTION FUNCTIONS
# =============================================================================

def parse_cdr_dict(cdr_string: str) -> Dict[str, str]:
    """Parse CDR dictionary from string format."""
    if pd.isna(cdr_string) or not cdr_string:
        return {}
    
    try:
        # Try ast.literal_eval first (handles Python dict strings)
        return ast.literal_eval(cdr_string)
    except (ValueError, SyntaxError):
        pass
    
    # Try manual parsing for format like "{'CDRH1': 'XXX', 'CDRH2': 'YYY'}"
    try:
        # Extract key-value pairs
        pattern = r"'(\w+)':\s*'([^']+)'"
        matches = re.findall(pattern, cdr_string)
        return {k: v for k, v in matches}
    except Exception:
        return {}


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


def run_anarci_single(sequence: str) -> Optional[Tuple[str, str, str, str]]:
    """Run ANARCI on a single sequence and extract CDRs."""
    if not ANARCI_AVAILABLE:
        return None
    
    try:
        result = run_anarci([("seq", sequence)], scheme="imgt")
        numberings = result[1]
        
        if not numberings or not numberings[0]:
            return None
        
        numbering = numberings[0][0][0]
        if not numbering:
            return None
        
        cdr1 = extract_region(numbering, 27.0, 38.99)
        cdr2 = extract_region(numbering, 56.0, 65.99)
        cdr3 = extract_region(numbering, 105.0, 117.99)
        full_v = extract_full_v_domain(numbering)
        
        return (full_v, cdr1, cdr2, cdr3)
    except Exception:
        return None


# =============================================================================
# MAIN PROCESSING
# =============================================================================

def load_vhh_csv(file_path: str) -> pd.DataFrame:
    """Load VHH CSV file."""
    print(f"Loading VHH CSV: {file_path}")
    
    # Try different separators
    for sep in [',', '\t']:
        try:
            df = pd.read_csv(file_path, sep=sep, low_memory=False)
            if len(df.columns) > 1:
                break
        except Exception:
            continue
    
    print(f"  Loaded {len(df):,} rows, {len(df.columns)} columns")
    print(f"  Columns: {list(df.columns)}")
    
    return df


def process_vhh_sequences(df: pd.DataFrame, use_anarci: bool = False) -> pd.DataFrame:
    """Process VHH sequences and extract CDRs."""
    print(f"\nProcessing {len(df):,} VHH sequences...")
    
    # Find sequence column
    seq_col = None
    for col in ['sequence', 'Sequence', 'aa_seq', 'VHH_sequence']:
        if col in df.columns:
            seq_col = col
            break
    
    if seq_col is None:
        raise ValueError(f"Could not find sequence column. Available: {list(df.columns)}")
    
    print(f"  Sequence column: '{seq_col}'")
    
    # Find CDR column (if pre-computed)
    cdr_col = None
    for col in ['cdr_sequences', 'CDR_sequences', 'cdrs']:
        if col in df.columns:
            cdr_col = col
            break
    
    # Find target/antigen column
    target_cols = []
    for col in ['targets_mentioned', 'target', 'antigen', 'Target', 'Antigen', 
                'antigen_name', 'target_name', 'linked_targets']:
        if col in df.columns:
            target_cols.append(col)
    
    if target_cols:
        print(f"  Target/antigen columns found: {target_cols}")
    
    # Process each sequence
    results = []
    
    for idx, row in tqdm(df.iterrows(), total=len(df), desc="  Processing"):
        seq = row[seq_col]
        
        if pd.isna(seq) or len(str(seq)) < 50:
            continue
        
        seq = str(seq).upper().replace(' ', '').replace('-', '')
        
        # Get CDRs - either from pre-computed or via ANARCI
        cdr1 = cdr2 = cdr3 = full_v = ""
        
        if cdr_col and not use_anarci:
            # Use pre-computed CDRs
            cdr_dict = parse_cdr_dict(row.get(cdr_col, ""))
            cdr1 = cdr_dict.get('CDRH1', cdr_dict.get('cdr1', ''))
            cdr2 = cdr_dict.get('CDRH2', cdr_dict.get('cdr2', ''))
            cdr3 = cdr_dict.get('CDRH3', cdr_dict.get('cdr3', ''))
            full_v = seq  # Use full sequence as V-domain approximation
        
        if use_anarci or not cdr3:
            # Run ANARCI for proper numbering
            anarci_result = run_anarci_single(seq)
            if anarci_result:
                full_v, cdr1, cdr2, cdr3 = anarci_result
            elif not cdr3:
                # Skip if no CDRs available
                continue
        
        # Build result record
        record = {
            'source': row.get('source', 'GenBank'),
            'seq_id': row.get('ID', row.get('model', f"vhh_{idx}")),
            'aa_v_full': full_v if full_v else seq,
            'cdr1': cdr1,
            'cdr2': cdr2,
            'cdr3': cdr3,
            'len_v': len(full_v) if full_v else len(seq),
            'len_cdr1': len(cdr1),
            'len_cdr2': len(cdr2),
            'len_cdr3': len(cdr3),
            'aa_raw': seq,
        }
        
        # Add target/antigen info
        for tcol in target_cols:
            val = row.get(tcol, '')
            if pd.notna(val) and val:
                record['targets'] = str(val)
                break
        
        # Add other useful metadata
        for meta_col in ['definition', 'organism', 'reference_title', 'reference_authors', 
                         'update_date', 'type', 'model']:
            if meta_col in df.columns:
                val = row.get(meta_col, '')
                if pd.notna(val):
                    record[meta_col] = str(val)
        
        results.append(record)
    
    df_out = pd.DataFrame(results)
    print(f"  Processed: {len(df_out):,} sequences")
    
    return df_out


def load_oas_csv(file_path: str) -> pd.DataFrame:
    """Load OAS camel processed CSV."""
    print(f"\nLoading OAS camel CSV: {file_path}")
    df = pd.read_csv(file_path)
    print(f"  Loaded {len(df):,} rows")
    
    # Standardize column names to match VHH format
    df['source'] = 'OAS_Camel'
    
    # Add empty targets column if not present
    if 'targets' not in df.columns:
        df['targets'] = ''
    
    return df


def merge_databases(vhh_df: pd.DataFrame, oas_df: pd.DataFrame) -> pd.DataFrame:
    """Merge VHH and OAS databases."""
    print(f"\nMerging databases...")
    print(f"  VHH sequences: {len(vhh_df):,}")
    print(f"  OAS sequences: {len(oas_df):,}")
    
    # Standardize columns
    common_cols = ['source', 'seq_id', 'aa_v_full', 'cdr1', 'cdr2', 'cdr3',
                   'len_v', 'len_cdr1', 'len_cdr2', 'len_cdr3', 'targets']
    
    # Ensure all columns exist
    for col in common_cols:
        if col not in vhh_df.columns:
            vhh_df[col] = ''
        if col not in oas_df.columns:
            oas_df[col] = ''
    
    # Select and order columns
    vhh_subset = vhh_df[common_cols].copy()
    oas_subset = oas_df[common_cols].copy()
    
    # Concatenate
    merged = pd.concat([vhh_subset, oas_subset], ignore_index=True)
    
    # Generate unique IDs
    merged['unified_id'] = [f"vhh_{i:07d}" for i in range(len(merged))]
    
    print(f"  Merged total: {len(merged):,} sequences")
    
    # Source breakdown
    print(f"\n  Source breakdown:")
    for src, count in merged['source'].value_counts().items():
        print(f"    {src}: {count:,}")
    
    return merged


def save_outputs(df: pd.DataFrame, output_prefix: str):
    """Save output files."""
    print(f"\n{'='*60}")
    print("SAVING OUTPUTS")
    print(f"{'='*60}")
    
    # Save full CSV
    csv_path = f"{output_prefix}.csv"
    df.to_csv(csv_path, index=False)
    print(f"\n  CSV saved: {csv_path}")
    print(f"    {len(df):,} sequences, {len(df.columns)} columns")
    
    # Save NPZ with antigen info
    npz_path = f"{output_prefix}.npz"
    
    # Prepare arrays
    arrays = {
        'ids': np.array(df['unified_id'].values if 'unified_id' in df.columns 
                        else df['seq_id'].values, dtype=object),
        'source': np.array(df['source'].values, dtype=object),
        'aa_v_full': np.array(df['aa_v_full'].values, dtype=object),
        'cdr1': np.array(df['cdr1'].values, dtype=object),
        'cdr2': np.array(df['cdr2'].values, dtype=object),
        'cdr3': np.array(df['cdr3'].values, dtype=object),
        'len_v': np.array(df['len_v'].values, dtype=np.int32),
        'len_cdr3': np.array(df['len_cdr3'].values, dtype=np.int32),
    }
    
    # Add targets/antigen info if present
    if 'targets' in df.columns:
        arrays['targets'] = np.array(df['targets'].fillna('').values, dtype=object)
    
    np.savez_compressed(npz_path, **arrays)
    print(f"  NPZ saved: {npz_path}")
    print(f"    Arrays: {list(arrays.keys())}")
    
    return csv_path, npz_path


def print_summary(df: pd.DataFrame):
    """Print summary statistics."""
    print(f"\n{'='*60}")
    print("SUMMARY STATISTICS")
    print(f"{'='*60}")
    
    print(f"\nTotal sequences: {len(df):,}")
    
    # Source breakdown
    if 'source' in df.columns:
        print(f"\nBy source:")
        for src, count in df['source'].value_counts().items():
            print(f"  {src}: {count:,}")
    
    # Sequences with target info
    if 'targets' in df.columns:
        has_targets = (df['targets'].fillna('') != '').sum()
        print(f"\nWith antigen/target info: {has_targets:,} ({has_targets/len(df)*100:.1f}%)")
    
    # Length statistics
    if 'len_v' in df.columns:
        v_lens = df['len_v'].values
        print(f"\nV-domain length:")
        print(f"  Min:    {np.min(v_lens)}")
        print(f"  Max:    {np.max(v_lens)}")
        print(f"  Mean:   {np.mean(v_lens):.1f}")
        print(f"  Median: {np.median(v_lens):.1f}")
    
    if 'len_cdr3' in df.columns:
        cdr3_lens = df['len_cdr3'].values
        print(f"\nCDR3 length:")
        print(f"  Min:    {np.min(cdr3_lens)}")
        print(f"  Max:    {np.max(cdr3_lens)}")
        print(f"  Mean:   {np.mean(cdr3_lens):.1f}")
    
    # Example with targets
    if 'targets' in df.columns:
        with_targets = df[df['targets'].fillna('') != '']
        if len(with_targets) > 0:
            print(f"\nExample sequences with target info:")
            for i, (_, row) in enumerate(with_targets.head(3).iterrows()):
                print(f"\n  [{i+1}] {row.get('seq_id', 'N/A')}")
                print(f"      CDR3: {row['cdr3']}")
                targets = str(row['targets'])[:80]
                print(f"      Targets: {targets}{'...' if len(str(row['targets'])) > 80 else ''}")


def main():
    parser = argparse.ArgumentParser(
        description="Process VHH sequences with antigen info and merge with OAS data",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Process VHH CSV only
    python process_vhh_with_antigens.py \\
        --vhh-csv /path/to/vhh_sequences.csv \\
        --output-prefix vhh_annotated

    # Merge with OAS camel data
    python process_vhh_with_antigens.py \\
        --vhh-csv /path/to/vhh_sequences.csv \\
        --oas-csv /path/to/camel_vhh_clean_anarci_renumbered.csv \\
        --output-prefix VHH_db_unified
        """
    )
    
    parser.add_argument(
        "--vhh-csv", required=True,
        help="Path to VHH sequences CSV (GenBank/PLAbDab format)"
    )
    parser.add_argument(
        "--oas-csv", default=None,
        help="Path to OAS camel processed CSV (optional, for merging)"
    )
    parser.add_argument(
        "--output-prefix", required=True,
        help="Output file prefix"
    )
    parser.add_argument(
        "--use-anarci", action="store_true",
        help="Re-run ANARCI even if CDRs are pre-computed"
    )
    
    args = parser.parse_args()
    
    # Header
    print("=" * 60)
    print("VHH SEQUENCE PROCESSOR WITH ANTIGEN INFO")
    print("=" * 60)
    print(f"VHH CSV:     {args.vhh_csv}")
    print(f"OAS CSV:     {args.oas_csv or 'None (VHH only)'}")
    print(f"Output:      {args.output_prefix}.*")
    print(f"Use ANARCI:  {args.use_anarci}")
    print(f"ANARCI:      {'Available' if ANARCI_AVAILABLE else 'Not available'}")
    print("=" * 60)
    
    # Load and process VHH CSV
    vhh_df = load_vhh_csv(args.vhh_csv)
    vhh_processed = process_vhh_sequences(vhh_df, use_anarci=args.use_anarci)
    
    # Optionally merge with OAS data
    if args.oas_csv:
        oas_df = load_oas_csv(args.oas_csv)
        merged_df = merge_databases(vhh_processed, oas_df)
        output_df = merged_df
    else:
        # Just use VHH data
        vhh_processed['unified_id'] = vhh_processed['seq_id']
        output_df = vhh_processed
    
    # Save outputs
    csv_path, npz_path = save_outputs(output_df, args.output_prefix)
    
    # Print summary
    print_summary(output_df)
    
    print(f"\n{'='*60}")
    print("COMPLETE!")
    print(f"{'='*60}")
    print(f"\nOutput files:")
    print(f"  CSV: {csv_path}")
    print(f"  NPZ: {npz_path}")
    print(f"\nThe NPZ includes 'targets' array for antigen info in search results!")
    print("=" * 60)


if __name__ == "__main__":
    main()
