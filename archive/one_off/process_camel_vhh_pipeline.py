#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
OAS Camel VHH Complete Pipeline

This script processes OAS CSV/TSV files containing camel heavy-chain nucleotide 
sequences and creates a clean VHH database with proper ANARCI numbering.

Pipeline steps:
1. Load all CSV/TSV files from a directory (or single file)
2. Extract nucleotide sequences from first column or 'sequence' column
3. Translate nucleotides → amino acids
4. Run ANARCI for IMGT numbering
5. Extract CDR1, CDR2, CDR3 and full V-domain
6. Save clean CSV and NPZ outputs

Usage:
    # Process all CSV files in a directory
    python process_camel_vhh_pipeline.py \
        --input-dir /path/to/csv/files/ \
        --output-prefix camel_vhh_clean \
        --batch-size 2000

    # Process a single file
    python process_camel_vhh_pipeline.py \
        --input-file /path/to/file.csv \
        --output-prefix camel_vhh_clean

Author: Claude (Anthropic)
Date: 2025
"""

import os
import sys
import glob
import argparse
import warnings
from typing import List, Dict, Tuple, Optional
from pathlib import Path

import numpy as np
import pandas as pd
from tqdm import tqdm

# Try to import Biopython for translation
try:
    from Bio.Seq import Seq
    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False
    warnings.warn("Biopython not available. Using fallback codon table.")

# Try to import ANARCI
try:
    from anarci import run_anarci
    ANARCI_AVAILABLE = True
except ImportError:
    ANARCI_AVAILABLE = False


# =============================================================================
# CONSTANTS
# =============================================================================

# Standard genetic code (fallback if Biopython not available)
CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}


# =============================================================================
# TRANSLATION FUNCTIONS
# =============================================================================

def clean_nt_sequence(nt_seq: str) -> str:
    """Clean nucleotide sequence."""
    if not nt_seq or not isinstance(nt_seq, str):
        return ""
    
    # Uppercase and remove whitespace, gaps
    nt_seq = nt_seq.upper().replace(' ', '').replace('\n', '').replace('-', '').replace('.', '')
    
    # Keep only valid nucleotide characters
    nt_seq = ''.join(c for c in nt_seq if c in 'ATCGN')
    
    return nt_seq


def translate_nt_to_aa(nt_seq: str) -> str:
    """Translate nucleotide sequence to amino acid."""
    nt_seq = clean_nt_sequence(nt_seq)
    
    if len(nt_seq) < 3:
        return ""
    
    # Trim to multiple of 3
    trim = len(nt_seq) % 3
    if trim:
        nt_seq = nt_seq[:-trim]
    
    if BIOPYTHON_AVAILABLE:
        try:
            # Replace N with A for translation (or could use random)
            nt_seq_clean = nt_seq.replace('N', 'A')
            aa = str(Seq(nt_seq_clean).translate(to_stop=False))
            return aa
        except Exception:
            pass
    
    # Fallback translation
    aa_list = []
    for i in range(0, len(nt_seq), 3):
        codon = nt_seq[i:i+3]
        if 'N' in codon:
            aa_list.append('X')  # Unknown
        elif codon in CODON_TABLE:
            aa_list.append(CODON_TABLE[codon])
        else:
            aa_list.append('X')
    
    return ''.join(aa_list)


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
# FILE LOADING
# =============================================================================

def detect_separator(file_path: str) -> str:
    """Detect CSV/TSV separator."""
    with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
        first_line = f.readline()
    
    if '\t' in first_line:
        return '\t'
    return ','


def find_nt_column(df: pd.DataFrame) -> str:
    """Find the nucleotide sequence column."""
    
    # Check for 'sequence' column (preferred)
    if 'sequence' in df.columns:
        return 'sequence'
    
    # Case-insensitive check for 'sequence'
    for col in df.columns:
        if col.lower() == 'sequence':
            return col
    
    # Check first column - if it looks like nucleotides, use it
    first_col = df.columns[0]
    sample = df[first_col].dropna().iloc[0] if len(df) > 0 else ""
    
    if isinstance(sample, str):
        sample_upper = sample.upper().replace('-', '').replace('.', '')[:100]
        nt_chars = sum(1 for c in sample_upper if c in 'ATCGN')
        if len(sample_upper) > 0 and nt_chars / len(sample_upper) > 0.8:
            print(f"  Using first column '{first_col}' as nucleotide sequence")
            return first_col
    
    # Try other common names
    for name in ['nt_sequence', 'nucleotide', 'nt_seq', 'sequence_nt']:
        if name in df.columns:
            return name
        for col in df.columns:
            if col.lower() == name.lower():
                return col
    
    raise ValueError(f"Could not find nucleotide sequence column. Columns: {list(df.columns)[:10]}")


def load_single_file(file_path: str) -> pd.DataFrame:
    """Load a single CSV/TSV file, handling OAS format with JSON metadata on first line."""
    sep = detect_separator(file_path)
    
    # First, try loading with skiprows=1 (OAS format has JSON metadata on line 1)
    # This is safer for OAS files
    try:
        df = pd.read_csv(file_path, sep=sep, low_memory=False, skiprows=1)
        
        # Check if we got the expected 'sequence' column
        if 'sequence' in df.columns:
            print(f"    (OAS format - skipped metadata line)")
            return df
    except Exception:
        pass
    
    # Fallback: try loading without skipping rows
    try:
        df = pd.read_csv(file_path, sep=sep, low_memory=False)
        
        # Check if first column looks like JSON metadata
        if len(df.columns) > 0:
            first_col = str(df.columns[0])
            if '{' in first_col or 'Run' in first_col:
                # First row was JSON, need to skip it
                df = pd.read_csv(file_path, sep=sep, low_memory=False, skiprows=1)
                print(f"    (OAS format detected, reloading with skiprows=1)")
        
        return df
    except Exception as e:
        # Try with different encoding
        try:
            df = pd.read_csv(file_path, sep=sep, low_memory=False, skiprows=1, encoding='latin-1')
            return df
        except Exception as e2:
            raise ValueError(f"Failed to load {file_path}: {e2}")


def load_all_files(input_dir: str = None, input_file: str = None, 
                   file_pattern: str = "*.csv") -> pd.DataFrame:
    """Load all CSV/TSV files from directory or single file."""
    
    if input_file:
        print(f"Loading single file: {input_file}")
        df = load_single_file(input_file)
        df['_source_file'] = os.path.basename(input_file)
        return df
    
    if input_dir:
        # Find all matching files - handle both regular files and nested directories
        all_files = []
        
        # Method 1: Direct CSV/TSV files in directory
        csv_files = glob.glob(os.path.join(input_dir, "*.csv"))
        tsv_files = glob.glob(os.path.join(input_dir, "*.tsv"))
        
        for f in csv_files + tsv_files:
            if os.path.isfile(f):
                all_files.append(f)
        
        # Method 2: Nested structure (OAS download format) - directories named *.csv containing files
        for item in os.listdir(input_dir):
            item_path = os.path.join(input_dir, item)
            if os.path.isdir(item_path) and item.endswith('.csv'):
                # Look for CSV file inside
                nested_csv = os.path.join(item_path, item)
                if os.path.isfile(nested_csv):
                    all_files.append(nested_csv)
                else:
                    # Try to find any CSV file in the directory
                    for nested_file in os.listdir(item_path):
                        if nested_file.endswith('.csv') and not nested_file.endswith('Zone.Identifier'):
                            all_files.append(os.path.join(item_path, nested_file))
        
        # Method 3: Recursive search for all CSV files
        if not all_files:
            print("  No direct CSV files found, searching recursively...")
            for root, dirs, files in os.walk(input_dir):
                for f in files:
                    if f.endswith('.csv') and not f.endswith('Zone.Identifier'):
                        all_files.append(os.path.join(root, f))
        
        # Remove duplicates while preserving order
        all_files = list(dict.fromkeys(all_files))
        
        if not all_files:
            raise FileNotFoundError(f"No CSV/TSV files found in {input_dir}")
        
        print(f"Found {len(all_files)} files to process:")
        for f in all_files:
            print(f"  - {os.path.relpath(f, input_dir)}")
        
        # Load and combine all files
        dfs = []
        for file_path in tqdm(all_files, desc="Loading files"):
            try:
                df = load_single_file(file_path)
                df['_source_file'] = os.path.basename(file_path)
                dfs.append(df)
                print(f"    Loaded {len(df):,} rows from {os.path.basename(file_path)}")
            except Exception as e:
                print(f"    Warning: Failed to load {file_path}: {e}")
        
        if not dfs:
            raise ValueError("No files could be loaded successfully")
        
        # Combine all dataframes
        combined = pd.concat(dfs, ignore_index=True)
        print(f"\nCombined: {len(combined):,} total rows from {len(dfs)} files")
        
        return combined
    
    raise ValueError("Must specify either --input-dir or --input-file")


# =============================================================================
# MAIN PIPELINE
# =============================================================================

def process_sequences(df: pd.DataFrame, nt_col: str, batch_size: int, 
                      min_nt_len: int = 150, min_aa_len: int = 70) -> pd.DataFrame:
    """Process all sequences: translate, ANARCI, extract CDRs."""
    
    print(f"\n{'='*60}")
    print("PROCESSING SEQUENCES")
    print(f"{'='*60}")
    
    # Step 1: Filter and translate
    print("\n[Step 1] Filtering and translating nucleotide sequences...")
    
    # Remove null sequences
    df = df[df[nt_col].notna()].copy()
    print(f"  Non-null sequences: {len(df):,}")
    
    # Clean and get lengths
    df['nt_clean'] = df[nt_col].apply(clean_nt_sequence)
    df['nt_len'] = df['nt_clean'].str.len()
    
    # Filter by minimum nucleotide length
    df = df[df['nt_len'] >= min_nt_len].copy()
    print(f"  After min NT length ({min_nt_len}): {len(df):,}")
    
    # Translate
    print("  Translating to amino acids...")
    tqdm.pandas(desc="  Translating")
    df['aa_raw'] = df['nt_clean'].progress_apply(translate_nt_to_aa)
    
    # Filter out sequences with internal stop codons
    def has_internal_stop(aa: str) -> bool:
        if not aa:
            return True
        if '*' in aa[:-1]:  # Stop codon not at end
            return True
        return False
    
    df = df[~df['aa_raw'].apply(has_internal_stop)].copy()
    print(f"  After removing internal stops: {len(df):,}")
    
    # Remove trailing stop
    df['aa_raw'] = df['aa_raw'].str.rstrip('*')
    
    # Filter by minimum AA length
    df['aa_len'] = df['aa_raw'].str.len()
    df = df[df['aa_len'] >= min_aa_len].copy()
    print(f"  After min AA length ({min_aa_len}): {len(df):,}")
    
    if len(df) == 0:
        raise ValueError("No sequences passed filtering!")
    
    # Step 2: Run ANARCI
    print(f"\n[Step 2] Running ANARCI numbering (batch size: {batch_size})...")
    
    if not ANARCI_AVAILABLE:
        print("ERROR: ANARCI is not installed!")
        print("Install with: pip install anarci")
        sys.exit(1)
    
    sequences = df['aa_raw'].tolist()
    n_seqs = len(sequences)
    
    # Initialize result columns
    results = {
        'anarci_status': ['pending'] * n_seqs,
        'aa_v_full': [''] * n_seqs,
        'cdr1': [''] * n_seqs,
        'cdr2': [''] * n_seqs,
        'cdr3': [''] * n_seqs,
        'len_v': [0] * n_seqs,
        'len_cdr1': [0] * n_seqs,
        'len_cdr2': [0] * n_seqs,
        'len_cdr3': [0] * n_seqs,
    }
    
    n_success = 0
    n_fail = 0
    
    for start_idx in tqdm(range(0, n_seqs, batch_size), desc="  ANARCI"):
        end_idx = min(start_idx + batch_size, n_seqs)
        batch_seqs = sequences[start_idx:end_idx]
        
        # Prepare ANARCI input
        anarci_input = [(f"seq_{i}", seq) for i, seq in enumerate(batch_seqs)]
        
        # Run ANARCI
        numberings = run_anarci_batch(anarci_input, scheme="imgt")
        
        # Process results
        for i, numbering in enumerate(numberings):
            global_idx = start_idx + i
            
            if numbering is None or len(numbering) == 0:
                results['anarci_status'][global_idx] = 'fail'
                n_fail += 1
                continue
            
            try:
                # Extract regions (IMGT boundaries for heavy chain)
                cdr1 = extract_region(numbering, 27.0, 38.99)
                cdr2 = extract_region(numbering, 56.0, 65.99)
                cdr3 = extract_region(numbering, 105.0, 117.99)
                full_v = extract_full_v_domain(numbering)
                
                # Check minimum V-domain length
                if len(full_v) < min_aa_len:
                    results['anarci_status'][global_idx] = 'short'
                    n_fail += 1
                    continue
                
                # Store results
                results['anarci_status'][global_idx] = 'ok'
                results['aa_v_full'][global_idx] = full_v
                results['cdr1'][global_idx] = cdr1
                results['cdr2'][global_idx] = cdr2
                results['cdr3'][global_idx] = cdr3
                results['len_v'][global_idx] = len(full_v)
                results['len_cdr1'][global_idx] = len(cdr1)
                results['len_cdr2'][global_idx] = len(cdr2)
                results['len_cdr3'][global_idx] = len(cdr3)
                n_success += 1
                
            except Exception as e:
                results['anarci_status'][global_idx] = 'error'
                n_fail += 1
    
    # Add results to dataframe
    for col, values in results.items():
        df[col] = values
    
    print(f"\n  ANARCI results:")
    print(f"    Success: {n_success:,} ({n_success/n_seqs*100:.1f}%)")
    print(f"    Failed:  {n_fail:,} ({n_fail/n_seqs*100:.1f}%)")
    
    return df


def save_outputs(df: pd.DataFrame, output_prefix: str, nt_col: str):
    """Save output files."""
    
    print(f"\n{'='*60}")
    print("SAVING OUTPUTS")
    print(f"{'='*60}")
    
    # Filter to successful entries
    df_ok = df[df['anarci_status'] == 'ok'].copy()
    print(f"\nSuccessfully processed sequences: {len(df_ok):,}")
    
    if len(df_ok) == 0:
        print("WARNING: No sequences were successfully processed!")
        return
    
    # Generate sequence IDs
    df_ok = df_ok.reset_index(drop=True)
    df_ok['seq_id'] = [f"camel_vhh_{i:07d}" for i in range(len(df_ok))]
    
    # Core output columns (always included)
    output_cols = ['seq_id', 'aa_v_full', 'cdr1', 'cdr2', 'cdr3', 
                   'len_v', 'len_cdr1', 'len_cdr2', 'len_cdr3',
                   'aa_raw', 'nt_len', 'anarci_status']
    
    # Add original nucleotide sequence
    if nt_col in df_ok.columns:
        df_ok['nt_sequence'] = df_ok[nt_col]
        output_cols.insert(1, 'nt_sequence')
    
    # Add useful OAS metadata columns if present
    oas_metadata_cols = [
        'v_call', 'd_call', 'j_call',  # Gene calls
        'cdr3_aa',  # Original OAS CDR3 for comparison
        'productive',
        'sequence_alignment_aa',  # Original OAS AA alignment for comparison
        '_source_file'  # Which file this came from
    ]
    
    for col in oas_metadata_cols:
        if col in df_ok.columns:
            output_cols.append(col)
    
    # Filter to columns that exist
    output_cols = [c for c in output_cols if c in df_ok.columns]
    
    df_out = df_ok[output_cols].copy()
    
    # Save CSV
    csv_path = f"{output_prefix}_anarci_renumbered.csv"
    df_out.to_csv(csv_path, index=False)
    print(f"\n  CSV saved: {csv_path}")
    print(f"    {len(df_out):,} sequences, {len(df_out.columns)} columns")
    
    # Save NPZ for downstream search
    npz_path = f"{output_prefix}_db.npz"
    
    np.savez_compressed(
        npz_path,
        ids=np.array(df_out['seq_id'].values, dtype=object),
        aa_v_full=np.array(df_out['aa_v_full'].values, dtype=object),
        cdr1=np.array(df_out['cdr1'].values, dtype=object),
        cdr2=np.array(df_out['cdr2'].values, dtype=object),
        cdr3=np.array(df_out['cdr3'].values, dtype=object),
        len_v=np.array(df_out['len_v'].values, dtype=np.int32),
        len_cdr3=np.array(df_out['len_cdr3'].values, dtype=np.int32),
    )
    print(f"  NPZ saved: {npz_path}")
    
    return df_out


def print_summary(df: pd.DataFrame):
    """Print summary statistics."""
    
    df_ok = df[df['anarci_status'] == 'ok']
    
    print(f"\n{'='*60}")
    print("SUMMARY STATISTICS")
    print(f"{'='*60}")
    
    print(f"\nTotal input sequences: {len(df):,}")
    print(f"Successfully processed: {len(df_ok):,} ({len(df_ok)/len(df)*100:.1f}%)")
    
    if len(df_ok) == 0:
        return
    
    # V-domain lengths
    v_lens = df_ok['len_v'].values
    print(f"\nV-domain length:")
    print(f"  Min:    {np.min(v_lens)}")
    print(f"  Max:    {np.max(v_lens)}")
    print(f"  Mean:   {np.mean(v_lens):.1f}")
    print(f"  Median: {np.median(v_lens):.1f}")
    
    # CDR3 lengths
    cdr3_lens = df_ok['len_cdr3'].values
    print(f"\nCDR3 length:")
    print(f"  Min:    {np.min(cdr3_lens)}")
    print(f"  Max:    {np.max(cdr3_lens)}")
    print(f"  Mean:   {np.mean(cdr3_lens):.1f}")
    print(f"  Median: {np.median(cdr3_lens):.1f}")
    
    # CDR1 lengths
    cdr1_lens = df_ok['len_cdr1'].values
    print(f"\nCDR1 length:")
    print(f"  Min:    {np.min(cdr1_lens)}")
    print(f"  Max:    {np.max(cdr1_lens)}")
    print(f"  Mean:   {np.mean(cdr1_lens):.1f}")
    
    # CDR2 lengths
    cdr2_lens = df_ok['len_cdr2'].values
    print(f"\nCDR2 length:")
    print(f"  Min:    {np.min(cdr2_lens)}")
    print(f"  Max:    {np.max(cdr2_lens)}")
    print(f"  Mean:   {np.mean(cdr2_lens):.1f}")
    
    # Quality check
    print(f"\n{'='*60}")
    print("QUALITY CHECK")
    print(f"{'='*60}")
    
    # V-domain in expected range (100-140 aa for VHH)
    v_good = ((v_lens >= 100) & (v_lens <= 140)).sum()
    print(f"\nV-domain 100-140 aa: {v_good:,} ({v_good/len(df_ok)*100:.1f}%)")
    
    # CDR3 in expected range (5-30 aa)
    cdr3_good = ((cdr3_lens >= 5) & (cdr3_lens <= 30)).sum()
    print(f"CDR3 5-30 aa: {cdr3_good:,} ({cdr3_good/len(df_ok)*100:.1f}%)")
    
    # Show examples
    print(f"\nExample sequences (first 3):")
    for i, (_, row) in enumerate(df_ok.head(3).iterrows()):
        print(f"\n  [{i+1}] V-domain ({row['len_v']} aa):")
        print(f"      {row['aa_v_full'][:60]}...")
        print(f"      CDR1 ({row['len_cdr1']} aa): {row['cdr1']}")
        print(f"      CDR2 ({row['len_cdr2']} aa): {row['cdr2']}")
        print(f"      CDR3 ({row['len_cdr3']} aa): {row['cdr3']}")


def main():
    parser = argparse.ArgumentParser(
        description="Process OAS camel VHH sequences: translate, ANARCI, extract CDRs",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
    # Process all CSV files in a directory
    python process_camel_vhh_pipeline.py \\
        --input-dir /home/sasenefrem/KA-Search/extracted/oas-paper/ \\
        --output-prefix camel_vhh_clean \\
        --batch-size 2000

    # Process a single file
    python process_camel_vhh_pipeline.py \\
        --input-file /path/to/file.csv \\
        --output-prefix camel_vhh_clean
        """
    )
    
    # Input options (mutually exclusive)
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        "--input-dir",
        help="Directory containing CSV/TSV files to process"
    )
    input_group.add_argument(
        "--input-file",
        help="Single CSV/TSV file to process"
    )
    
    parser.add_argument(
        "--output-prefix", required=True,
        help="Output file prefix (will create _anarci_renumbered.csv and _db.npz)"
    )
    parser.add_argument(
        "--batch-size", type=int, default=2000,
        help="ANARCI batch size (default: 2000)"
    )
    parser.add_argument(
        "--min-nt-len", type=int, default=150,
        help="Minimum nucleotide sequence length (default: 150)"
    )
    parser.add_argument(
        "--min-aa-len", type=int, default=70,
        help="Minimum amino acid sequence length (default: 70)"
    )
    parser.add_argument(
        "--max-sequences", type=int, default=None,
        help="Maximum sequences to process (for testing)"
    )
    
    args = parser.parse_args()
    
    # Header
    print("=" * 60)
    print("OAS CAMEL VHH PROCESSING PIPELINE")
    print("=" * 60)
    print(f"Input:       {args.input_dir or args.input_file}")
    print(f"Output:      {args.output_prefix}_*")
    print(f"Batch size:  {args.batch_size}")
    print(f"Min NT len:  {args.min_nt_len}")
    print(f"Min AA len:  {args.min_aa_len}")
    print(f"Biopython:   {'Available' if BIOPYTHON_AVAILABLE else 'Fallback'}")
    print(f"ANARCI:      {'Available' if ANARCI_AVAILABLE else 'NOT AVAILABLE!'}")
    print("=" * 60)
    
    if not ANARCI_AVAILABLE:
        print("\nERROR: ANARCI is required but not installed.")
        print("Install with: pip install anarci")
        sys.exit(1)
    
    # Create output directory if needed
    output_dir = os.path.dirname(args.output_prefix)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Step 1: Load files
    print(f"\n{'='*60}")
    print("LOADING DATA")
    print(f"{'='*60}")
    
    df = load_all_files(
        input_dir=args.input_dir,
        input_file=args.input_file
    )
    
    # Find nucleotide column
    nt_col = find_nt_column(df)
    print(f"\nNucleotide column: '{nt_col}'")
    
    # Limit for testing
    if args.max_sequences:
        df = df.head(args.max_sequences)
        print(f"Limited to first {args.max_sequences} sequences (testing)")
    
    # Step 2: Process sequences
    df = process_sequences(
        df, 
        nt_col=nt_col,
        batch_size=args.batch_size,
        min_nt_len=args.min_nt_len,
        min_aa_len=args.min_aa_len
    )
    
    # Step 3: Save outputs
    df_out = save_outputs(df, args.output_prefix, nt_col)
    
    # Step 4: Print summary
    print_summary(df)
    
    print(f"\n{'='*60}")
    print("PIPELINE COMPLETE!")
    print(f"{'='*60}")
    print(f"\nOutput files:")
    print(f"  CSV: {args.output_prefix}_anarci_renumbered.csv")
    print(f"  NPZ: {args.output_prefix}_db.npz")
    print(f"\nThe NPZ file can be used for similarity searching.")
    print("=" * 60)


if __name__ == "__main__":
    main()