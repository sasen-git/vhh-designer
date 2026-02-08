#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Build Clean Camelid VHH Database from OAS Nucleotide Data

This script:
1. Loads camel heavy-chain nucleotide sequences from an OAS TSV/CSV export
2. Translates them to amino acids (properly, from raw nucleotides)
3. Renumbers them with ANARCI (scheme="imgt")
4. Extracts full V-domain and CDR1/2/3
5. Saves a clean output (CSV and optionally NPZ) for downstream similarity search

This bypasses the truncated sequence_alignment_aa column in OAS and gives you
complete V-domains with proper IMGT numbering.

Usage:
    python build_camel_vhh_db.py \
        --input-tsv camel_heavy_oas.tsv \
        --output-prefix camel_vhh_anarci \
        --batch-size 2000 \
        --min-nt-len 150

Author: Claude (Anthropic)
Date: 2025
"""

import os
import sys
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
    warnings.warn("ANARCI not available. Cannot perform IMGT numbering.")


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

# Valid amino acids
VALID_AA = set('ACDEFGHIKLMNPQRSTVWY')

# Possible nucleotide column names in OAS (priority order)
NT_COLUMN_PRIORITY = [
    'sequence',
    'sequence_alignment', 
    'nt_sequence',
    'sequence_nt',
    'nucleotide_sequence',
]

# IMGT CDR boundaries (heavy chain)
IMGT_CDR_BOUNDARIES = {
    'cdr1': (27, 38),
    'cdr2': (56, 65),
    'cdr3': (105, 117),
}


# =============================================================================
# TRANSLATION FUNCTIONS
# =============================================================================

def translate_nt_to_aa_biopython(nt_seq: str) -> str:
    """Translate nucleotide sequence to amino acid using Biopython."""
    # Clean the sequence
    nt_seq = nt_seq.upper().replace(' ', '').replace('\n', '').replace('-', '').replace('.', '')
    
    # Remove non-ATCG characters
    nt_seq = ''.join(c for c in nt_seq if c in 'ATCGN')
    
    # Trim to multiple of 3
    trim = len(nt_seq) % 3
    if trim:
        nt_seq = nt_seq[:-trim]
    
    if len(nt_seq) < 3:
        return ""
    
    # Translate using Biopython
    try:
        aa = str(Seq(nt_seq).translate(to_stop=False))
        return aa
    except Exception as e:
        warnings.warn(f"Biopython translation failed: {e}")
        return ""


def translate_nt_to_aa_fallback(nt_seq: str) -> str:
    """Translate nucleotide sequence to amino acid using fallback codon table."""
    # Clean the sequence
    nt_seq = nt_seq.upper().replace(' ', '').replace('\n', '').replace('-', '').replace('.', '')
    
    # Remove non-ATCG characters (replace N with random choice or skip)
    nt_seq = ''.join(c for c in nt_seq if c in 'ATCG')
    
    # Trim to multiple of 3
    trim = len(nt_seq) % 3
    if trim:
        nt_seq = nt_seq[:-trim]
    
    if len(nt_seq) < 3:
        return ""
    
    # Translate codon by codon
    aa_list = []
    for i in range(0, len(nt_seq), 3):
        codon = nt_seq[i:i+3]
        if codon in CODON_TABLE:
            aa_list.append(CODON_TABLE[codon])
        else:
            aa_list.append('X')  # Unknown codon
    
    return ''.join(aa_list)


def nt_to_aa(nt_seq: str) -> str:
    """Translate nucleotide to amino acid sequence."""
    if not nt_seq or not isinstance(nt_seq, str):
        return ""
    
    if BIOPYTHON_AVAILABLE:
        return translate_nt_to_aa_biopython(nt_seq)
    else:
        return translate_nt_to_aa_fallback(nt_seq)


# =============================================================================
# ANARCI FUNCTIONS
# =============================================================================

def pos_val(pos_tuple: Tuple) -> float:
    """Convert IMGT position with insertion code to a float value."""
    pos = pos_tuple[0]
    ins = pos_tuple[1] if len(pos_tuple) > 1 else " "
    
    if isinstance(ins, str) and ins.strip():
        # A→.01, B→.02, etc
        ins_val = ord(ins.upper()) - 64  # A=1, B=2, etc.
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
        if start <= v <= end + 0.99:  # Include insertions up to end.99
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
    """
    Run ANARCI on a batch of sequences.
    
    Args:
        sequences: List of (name, sequence) tuples
        scheme: Numbering scheme (default: "imgt")
    
    Returns:
        List of numbering results (None if failed)
    """
    if not ANARCI_AVAILABLE:
        return [None] * len(sequences)
    
    if not sequences:
        return []
    
    try:
        # Run ANARCI
        result = run_anarci(sequences, scheme=scheme, allowed_species=None)
        
        # Extract numberings
        # Result structure: (numberings_list, alignment_details, hit_tables)
        numberings_list = result[1]
        
        parsed_results = []
        for i, entry in enumerate(numberings_list):
            try:
                if entry is None:
                    parsed_results.append(None)
                    continue
                
                # entry should be a list of (numbering, start, end) for each domain
                if not entry or not isinstance(entry, list):
                    parsed_results.append(None)
                    continue
                
                # Get first domain numbering
                domain_data = entry[0]
                if domain_data is None or not isinstance(domain_data, (list, tuple)):
                    parsed_results.append(None)
                    continue
                
                # domain_data = (numbering_list, start_pos, end_pos)
                numbering = domain_data[0]
                if numbering and isinstance(numbering, list):
                    parsed_results.append(numbering)
                else:
                    parsed_results.append(None)
                    
            except (IndexError, TypeError) as e:
                parsed_results.append(None)
        
        return parsed_results
        
    except Exception as e:
        warnings.warn(f"ANARCI batch failed: {e}")
        return [None] * len(sequences)


# =============================================================================
# DATA LOADING AND FILTERING
# =============================================================================

def find_nucleotide_column(df: pd.DataFrame) -> Optional[str]:
    """Find the nucleotide sequence column in the dataframe."""
    for col_name in NT_COLUMN_PRIORITY:
        if col_name in df.columns:
            return col_name
    
    # Try case-insensitive match
    lower_cols = {c.lower(): c for c in df.columns}
    for col_name in NT_COLUMN_PRIORITY:
        if col_name.lower() in lower_cols:
            return lower_cols[col_name.lower()]
    
    return None


def load_oas_data(input_path: str) -> pd.DataFrame:
    """Load OAS TSV/CSV file."""
    print(f"Loading data from: {input_path}")
    
    # Detect separator
    with open(input_path, 'r') as f:
        first_line = f.readline()
    
    if '\t' in first_line:
        sep = '\t'
    else:
        sep = ','
    
    # Load with pandas
    df = pd.read_csv(input_path, sep=sep, low_memory=False)
    print(f"  Loaded {len(df):,} rows, {len(df.columns)} columns")
    
    return df


def filter_sequences(df: pd.DataFrame, nt_col: str, min_nt_len: int) -> pd.DataFrame:
    """Apply quality filters to sequences."""
    original_count = len(df)
    
    # Filter 1: Non-null nucleotide sequence
    df = df[df[nt_col].notna()].copy()
    print(f"  After removing null nt_seq: {len(df):,} ({len(df)/original_count*100:.1f}%)")
    
    # Filter 2: Minimum nucleotide length
    df['nt_len'] = df[nt_col].str.len()
    df = df[df['nt_len'] >= min_nt_len].copy()
    print(f"  After min_nt_len >= {min_nt_len}: {len(df):,} ({len(df)/original_count*100:.1f}%)")
    
    # Filter 3: Productive (if column exists)
    if 'productive' in df.columns:
        before = len(df)
        df['productive_str'] = df['productive'].astype(str).str.upper()
        df = df[df['productive_str'].isin(['TRUE', 'T', '1', 'YES'])].copy()
        print(f"  After productive filter: {len(df):,} ({len(df)/before*100:.1f}% of previous)")
    
    # Filter 4: Heavy chain locus (if column exists)
    if 'locus' in df.columns:
        before = len(df)
        df = df[df['locus'].str.upper() == 'IGH'].copy()
        print(f"  After locus=IGH filter: {len(df):,} ({len(df)/before*100:.1f}% of previous)")
    
    return df


def translate_and_filter_aa(df: pd.DataFrame, nt_col: str) -> pd.DataFrame:
    """Translate nucleotides to amino acids and filter."""
    print("\nTranslating nucleotide sequences to amino acids...")
    
    # Translate
    tqdm.pandas(desc="Translating")
    df['aa_raw'] = df[nt_col].progress_apply(nt_to_aa)
    
    original_count = len(df)
    
    # Filter: Remove sequences with internal stop codons
    def has_internal_stop(aa: str) -> bool:
        if not aa:
            return True
        # Check for stop codon not at the end
        if '*' in aa[:-1]:
            return True
        return False
    
    df = df[~df['aa_raw'].apply(has_internal_stop)].copy()
    print(f"  After removing internal stops: {len(df):,} ({len(df)/original_count*100:.1f}%)")
    
    # Filter: Minimum amino acid length
    df['aa_len'] = df['aa_raw'].str.len()
    df = df[df['aa_len'] >= 70].copy()
    print(f"  After min aa_len >= 70: {len(df):,} ({len(df)/original_count*100:.1f}%)")
    
    # Remove trailing stop codon if present
    df['aa_raw'] = df['aa_raw'].str.rstrip('*')
    
    return df


# =============================================================================
# MAIN PROCESSING
# =============================================================================

def process_with_anarci(df: pd.DataFrame, batch_size: int) -> pd.DataFrame:
    """Run ANARCI numbering on all sequences and extract regions."""
    print(f"\nRunning ANARCI numbering in batches of {batch_size}...")
    
    if not ANARCI_AVAILABLE:
        print("ERROR: ANARCI is not available. Cannot proceed.")
        sys.exit(1)
    
    # Prepare input
    sequences = list(df['aa_raw'].values)
    n_seqs = len(sequences)
    
    # Initialize result columns
    results = {
        'anarci_status': ['pending'] * n_seqs,
        'aa_v_full': [''] * n_seqs,
        'cdr1': [''] * n_seqs,
        'cdr2': [''] * n_seqs,
        'cdr3': [''] * n_seqs,
        'len_v': [0] * n_seqs,
        'len_cdr3': [0] * n_seqs,
    }
    
    # Process in batches
    n_success = 0
    n_fail = 0
    
    for start_idx in tqdm(range(0, n_seqs, batch_size), desc="ANARCI batches"):
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
                # Extract regions
                cdr1 = extract_region(numbering, 27.0, 38.99)
                cdr2 = extract_region(numbering, 56.0, 65.99)
                cdr3 = extract_region(numbering, 105.0, 117.99)
                full_v = extract_full_v_domain(numbering)
                
                # Check minimum lengths
                if len(full_v) < 70:
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
                results['len_cdr3'][global_idx] = len(cdr3)
                n_success += 1
                
            except Exception as e:
                results['anarci_status'][global_idx] = 'error'
                n_fail += 1
    
    # Add results to dataframe
    for col, values in results.items():
        df[col] = values
    
    print(f"\nANARCI completed:")
    print(f"  Success: {n_success:,} ({n_success/n_seqs*100:.1f}%)")
    print(f"  Failed:  {n_fail:,} ({n_fail/n_seqs*100:.1f}%)")
    
    return df


def save_outputs(df: pd.DataFrame, output_prefix: str):
    """Save output files."""
    print(f"\nSaving outputs with prefix: {output_prefix}")
    
    # Filter to successful entries only
    df_ok = df[df['anarci_status'] == 'ok'].copy()
    print(f"  Saving {len(df_ok):,} successfully processed sequences")
    
    # Determine ID column
    if 'sequence_id' in df.columns:
        id_col = 'sequence_id'
    elif 'sequence_id_original' in df.columns:
        id_col = 'sequence_id_original'
    else:
        df_ok['sequence_id'] = range(len(df_ok))
        id_col = 'sequence_id'
    
    # Select columns for output
    output_cols = [id_col]
    
    # Add optional columns if they exist
    for col in ['nt_seq', 'aa_raw', 'aa_v_full', 'cdr1', 'cdr2', 'cdr3', 
                'len_v', 'len_cdr3', 'anarci_status', 'v_call', 'productive']:
        if col in df_ok.columns:
            output_cols.append(col)
    
    # Rename nt column if needed
    if 'nt_seq' not in df_ok.columns:
        for col in NT_COLUMN_PRIORITY:
            if col in df_ok.columns:
                df_ok['nt_seq'] = df_ok[col]
                if 'nt_seq' not in output_cols:
                    output_cols.insert(1, 'nt_seq')
                break
    
    df_out = df_ok[output_cols].copy()
    
    # Save CSV
    csv_path = f"{output_prefix}_anarci_renumbered.csv"
    df_out.to_csv(csv_path, index=False)
    print(f"  CSV saved: {csv_path}")
    
    # Save NPZ for downstream search
    npz_path = f"{output_prefix}_db.npz"
    
    np.savez_compressed(
        npz_path,
        ids=np.array(df_out[id_col].values, dtype=object),
        aa_v_full=np.array(df_out['aa_v_full'].values, dtype=object),
        cdr1=np.array(df_out['cdr1'].values, dtype=object),
        cdr2=np.array(df_out['cdr2'].values, dtype=object),
        cdr3=np.array(df_out['cdr3'].values, dtype=object),
    )
    print(f"  NPZ saved: {npz_path}")
    
    return df_out


def print_sanity_checks(df: pd.DataFrame):
    """Print sanity check statistics."""
    df_ok = df[df['anarci_status'] == 'ok']
    
    print("\n" + "="*70)
    print("SANITY CHECK STATISTICS")
    print("="*70)
    
    print(f"\nTotal sequences processed: {len(df):,}")
    print(f"ANARCI success: {len(df_ok):,} ({len(df_ok)/len(df)*100:.1f}%)")
    
    if len(df_ok) == 0:
        print("\nNo successful sequences - cannot compute statistics.")
        return
    
    # V-domain length stats
    v_lens = df_ok['len_v'].values
    print(f"\nV-domain length statistics:")
    print(f"  Min:    {np.min(v_lens)}")
    print(f"  Max:    {np.max(v_lens)}")
    print(f"  Mean:   {np.mean(v_lens):.1f}")
    print(f"  Median: {np.median(v_lens):.1f}")
    print(f"  Q25:    {np.percentile(v_lens, 25):.1f}")
    print(f"  Q75:    {np.percentile(v_lens, 75):.1f}")
    
    # CDR3 length stats
    cdr3_lens = df_ok['len_cdr3'].values
    print(f"\nCDR3 length statistics:")
    print(f"  Min:    {np.min(cdr3_lens)}")
    print(f"  Max:    {np.max(cdr3_lens)}")
    print(f"  Mean:   {np.mean(cdr3_lens):.1f}")
    print(f"  Median: {np.median(cdr3_lens):.1f}")
    print(f"  Q25:    {np.percentile(cdr3_lens, 25):.1f}")
    print(f"  Q75:    {np.percentile(cdr3_lens, 75):.1f}")
    
    # Expected ranges check
    print(f"\nExpected ranges check:")
    
    # VHH V-domains should be ~110-130 aa
    v_in_range = ((v_lens >= 100) & (v_lens <= 140)).sum()
    print(f"  V-domain 100-140 aa: {v_in_range:,} ({v_in_range/len(df_ok)*100:.1f}%)")
    
    # VHH CDR3 typically 8-24 aa
    cdr3_in_range = ((cdr3_lens >= 5) & (cdr3_lens <= 30)).sum()
    print(f"  CDR3 5-30 aa: {cdr3_in_range:,} ({cdr3_in_range/len(df_ok)*100:.1f}%)")
    
    # Show some examples
    print(f"\nExample sequences (first 3):")
    for i, (_, row) in enumerate(df_ok.head(3).iterrows()):
        print(f"\n  [{i+1}]")
        print(f"    V-domain ({row['len_v']} aa): {row['aa_v_full'][:50]}...")
        print(f"    CDR1: {row['cdr1']}")
        print(f"    CDR2: {row['cdr2']}")
        print(f"    CDR3 ({row['len_cdr3']} aa): {row['cdr3']}")
    
    print("\n" + "="*70)


# =============================================================================
# MAIN
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Build clean camelid VHH database from OAS nucleotide data",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
    python build_camel_vhh_db.py \\
        --input-tsv camel_heavy_oas.tsv \\
        --output-prefix camel_vhh_anarci \\
        --batch-size 2000 \\
        --min-nt-len 150

This script translates nucleotide sequences to amino acids and renumbers
them with ANARCI (IMGT scheme) to get complete V-domains and CDR regions.
        """
    )
    
    parser.add_argument(
        "--input-tsv", required=True,
        help="Path to OAS camel heavy-chain export (TSV/CSV)"
    )
    parser.add_argument(
        "--output-prefix", required=True,
        help="Base path for output files (no extension)"
    )
    parser.add_argument(
        "--batch-size", type=int, default=2000,
        help="Number of sequences per ANARCI batch (default: 2000)"
    )
    parser.add_argument(
        "--min-nt-len", type=int, default=150,
        help="Minimum nucleotide length to keep (default: 150)"
    )
    parser.add_argument(
        "--max-sequences", type=int, default=None,
        help="Maximum sequences to process (for testing)"
    )
    
    args = parser.parse_args()
    
    # Print header
    print("="*70)
    print("BUILD CAMELID VHH DATABASE FROM OAS NUCLEOTIDE DATA")
    print("="*70)
    print(f"Input:        {args.input_tsv}")
    print(f"Output:       {args.output_prefix}_*")
    print(f"Batch size:   {args.batch_size}")
    print(f"Min nt len:   {args.min_nt_len}")
    print(f"Biopython:    {'Available' if BIOPYTHON_AVAILABLE else 'Not available (using fallback)'}")
    print(f"ANARCI:       {'Available' if ANARCI_AVAILABLE else 'NOT AVAILABLE - REQUIRED!'}")
    print("="*70)
    
    if not ANARCI_AVAILABLE:
        print("\nERROR: ANARCI is required but not installed.")
        print("Install with: pip install anarci")
        sys.exit(1)
    
    # Check input file exists
    if not os.path.exists(args.input_tsv):
        print(f"\nERROR: Input file not found: {args.input_tsv}")
        sys.exit(1)
    
    # Create output directory if needed
    output_dir = os.path.dirname(args.output_prefix)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created output directory: {output_dir}")
    
    # Step 1: Load data
    print("\n" + "="*70)
    print("STEP 1: Loading OAS data")
    print("="*70)
    df = load_oas_data(args.input_tsv)
    
    # Find nucleotide column
    nt_col = find_nucleotide_column(df)
    if nt_col is None:
        print(f"\nERROR: Could not find nucleotide sequence column.")
        print(f"Expected one of: {NT_COLUMN_PRIORITY}")
        print(f"Found columns: {list(df.columns)}")
        sys.exit(1)
    print(f"  Using nucleotide column: '{nt_col}'")
    
    # Limit for testing
    if args.max_sequences:
        df = df.head(args.max_sequences)
        print(f"  Limited to first {args.max_sequences} sequences (testing mode)")
    
    # Step 2: Filter sequences
    print("\n" + "="*70)
    print("STEP 2: Filtering sequences")
    print("="*70)
    df = filter_sequences(df, nt_col, args.min_nt_len)
    
    if len(df) == 0:
        print("\nERROR: No sequences passed filtering!")
        sys.exit(1)
    
    # Step 3: Translate to amino acids
    print("\n" + "="*70)
    print("STEP 3: Translating nucleotides to amino acids")
    print("="*70)
    df = translate_and_filter_aa(df, nt_col)
    
    if len(df) == 0:
        print("\nERROR: No sequences passed translation filtering!")
        sys.exit(1)
    
    # Step 4: Run ANARCI
    print("\n" + "="*70)
    print("STEP 4: Running ANARCI numbering")
    print("="*70)
    df = process_with_anarci(df, args.batch_size)
    
    # Step 5: Save outputs
    print("\n" + "="*70)
    print("STEP 5: Saving outputs")
    print("="*70)
    df_out = save_outputs(df, args.output_prefix)
    
    # Step 6: Sanity checks
    print_sanity_checks(df)
    
    print("\n" + "="*70)
    print("COMPLETE!")
    print("="*70)
    print(f"Output files:")
    print(f"  CSV: {args.output_prefix}_anarci_renumbered.csv")
    print(f"  NPZ: {args.output_prefix}_db.npz")
    print("="*70)


if __name__ == "__main__":
    main()
