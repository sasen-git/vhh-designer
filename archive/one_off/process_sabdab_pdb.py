#!/usr/bin/env python3
"""
Extract VHH sequences from SAbDab IMGT-numbered PDB files and merge with unified database.

This script:
1. Reads sabdab_nano_summary_all.tsv for metadata
2. Extracts sequences from IMGT-numbered PDB files
3. Extracts CDR regions using IMGT positions
4. Optionally merges with existing unified VHH database
"""

import os
import sys
import argparse
import pandas as pd
import numpy as np
from pathlib import Path
from collections import defaultdict

# 3-letter to 1-letter amino acid mapping
AA_MAP = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
    'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
    'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
    'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y',
    # Non-standard (map to X)
    'MSE': 'M',  # Selenomethionine
    'SEC': 'U',  # Selenocysteine
    'PYL': 'O',  # Pyrrolysine
}

# IMGT CDR definitions (positions are inclusive)
IMGT_CDR1_START, IMGT_CDR1_END = 27, 38
IMGT_CDR2_START, IMGT_CDR2_END = 56, 65
IMGT_CDR3_START, IMGT_CDR3_END = 105, 117
IMGT_V_END = 128  # End of V-domain


def extract_sequence_from_pdb(pdb_path: str, chain_id: str) -> tuple:
    """
    Extract amino acid sequence from a PDB file for a specific chain.
    Uses CA atoms to get one residue per position.
    
    Returns: (full_sequence, position_dict) where position_dict maps IMGT position -> AA
    """
    position_to_aa = {}
    
    try:
        with open(pdb_path, 'r') as f:
            for line in f:
                if not line.startswith('ATOM'):
                    continue
                
                # PDB format columns:
                # 1-6: Record name (ATOM)
                # 13-16: Atom name
                # 17-20: Residue name (right-justified, often 17-19 used)
                # 22: Chain ID
                # 23-26: Residue sequence number
                # 27: Insertion code
                
                atom_name = line[12:16].strip()
                if atom_name != 'CA':  # Only use CA atoms (one per residue)
                    continue
                
                chain = line[21]
                if chain != chain_id:
                    continue
                
                res_name = line[17:20].strip()
                try:
                    res_num = int(line[22:26].strip())
                except ValueError:
                    continue
                
                # Skip if already seen this position (handles alt conformations)
                if res_num in position_to_aa:
                    continue
                
                # Convert to 1-letter code
                aa = AA_MAP.get(res_name, 'X')
                position_to_aa[res_num] = aa
    
    except FileNotFoundError:
        return None, None
    except Exception as e:
        print(f"Error reading {pdb_path}: {e}", file=sys.stderr)
        return None, None
    
    if not position_to_aa:
        return None, None
    
    # Build full sequence string from sorted positions
    sorted_positions = sorted(position_to_aa.keys())
    full_sequence = ''.join(position_to_aa[pos] for pos in sorted_positions)
    
    return full_sequence, position_to_aa


def extract_cdr_from_positions(position_dict: dict, start: int, end: int) -> str:
    """Extract CDR sequence from position dictionary given IMGT start/end."""
    if position_dict is None:
        return ''
    
    cdr_residues = []
    for pos in range(start, end + 1):
        if pos in position_dict:
            cdr_residues.append(position_dict[pos])
    
    return ''.join(cdr_residues)


def extract_v_domain(position_dict: dict) -> str:
    """Extract V-domain (positions 1-128) as contiguous sequence."""
    if position_dict is None:
        return ''
    
    v_residues = []
    for pos in sorted(position_dict.keys()):
        if pos <= IMGT_V_END:
            v_residues.append(position_dict[pos])
    
    return ''.join(v_residues)


def process_sabdab_data(tsv_path: str, pdb_dir: str, use_imgt: bool = True) -> pd.DataFrame:
    """
    Process SAbDab TSV and extract sequences from PDB files.
    
    Args:
        tsv_path: Path to sabdab_nano_summary_all.tsv
        pdb_dir: Directory containing PDB files (or imgt/ subdirectory)
        use_imgt: If True, use imgt/ subdirectory for IMGT-numbered structures
    
    Returns:
        DataFrame with sequences and metadata
    """
    # Read TSV
    print(f"Reading TSV: {tsv_path}")
    df = pd.read_csv(tsv_path, sep='\t')
    print(f"Loaded {len(df)} entries from TSV")
    
    # Determine PDB directory
    if use_imgt:
        actual_pdb_dir = os.path.join(pdb_dir, 'imgt')
        if not os.path.exists(actual_pdb_dir):
            print(f"IMGT directory not found, using main directory: {pdb_dir}")
            actual_pdb_dir = pdb_dir
    else:
        actual_pdb_dir = pdb_dir
    
    print(f"Using PDB directory: {actual_pdb_dir}")
    
    # Process each entry
    results = []
    seen_sequences = set()  # Track unique sequences
    
    for idx, row in df.iterrows():
        if idx % 500 == 0:
            print(f"Processing entry {idx}/{len(df)}...")
        
        pdb_id = str(row['pdb']).lower()
        chain_id = row['Hchain']
        
        if pd.isna(chain_id) or chain_id == 'NA':
            continue
        
        # Find PDB file
        pdb_path = os.path.join(actual_pdb_dir, f"{pdb_id}.pdb")
        if not os.path.exists(pdb_path):
            continue
        
        # Extract sequence
        full_seq, pos_dict = extract_sequence_from_pdb(pdb_path, chain_id)
        
        if full_seq is None or len(full_seq) < 50:  # Skip very short sequences
            continue
        
        # Extract V-domain and CDRs
        v_domain = extract_v_domain(pos_dict)
        cdr1 = extract_cdr_from_positions(pos_dict, IMGT_CDR1_START, IMGT_CDR1_END)
        cdr2 = extract_cdr_from_positions(pos_dict, IMGT_CDR2_START, IMGT_CDR2_END)
        cdr3 = extract_cdr_from_positions(pos_dict, IMGT_CDR3_START, IMGT_CDR3_END)
        
        # Skip if V-domain too short
        if len(v_domain) < 70:
            continue
        
        # Create unique ID
        seq_id = f"SAbDab_{pdb_id}_{chain_id}"
        
        # Skip duplicate sequences (keep first occurrence)
        if v_domain in seen_sequences:
            continue
        seen_sequences.add(v_domain)
        
        # Build target string from antigen info
        targets = []
        if pd.notna(row.get('antigen_name')) and row['antigen_name'] != 'NA':
            targets.append(str(row['antigen_name']))
        
        # Collect metadata
        result = {
            'id': seq_id,
            'source': 'SAbDab_PDB',
            'pdb_id': pdb_id,
            'chain': chain_id,
            'aa_v_full': v_domain,
            'cdr1': cdr1,
            'cdr2': cdr2,
            'cdr3': cdr3,
            'len_v': len(v_domain),
            'len_cdr3': len(cdr3),
            'targets': '; '.join(targets) if targets else '',
            'antigen_type': row.get('antigen_type', ''),
            'species': row.get('heavy_species', ''),
            'antigen_species': row.get('antigen_species', ''),
            'resolution': row.get('resolution', ''),
            'method': row.get('method', ''),
            'heavy_subclass': row.get('heavy_subclass', ''),
            'pmid': row.get('pmid', ''),
            'authors': row.get('authors', ''),
            'compound': row.get('compound', ''),
        }
        results.append(result)
    
    result_df = pd.DataFrame(results)
    print(f"\nExtracted {len(result_df)} unique VHH sequences from PDB files")
    
    return result_df


def deduplicate_database(df: pd.DataFrame) -> pd.DataFrame:
    """
    Remove duplicate sequences, keeping the entry with the most metadata.
    Priority: entries with target info > entries with longer IDs (more info) > first occurrence
    """
    print(f"\n=== Scanning for Duplicate Sequences ===")
    initial_count = len(df)
    
    # Find duplicates based on aa_v_full
    dup_mask = df.duplicated(subset=['aa_v_full'], keep=False)
    num_in_dup_groups = dup_mask.sum()
    
    if num_in_dup_groups == 0:
        print("No duplicate sequences found!")
        return df
    
    # Get duplicate groups
    dup_sequences = df[dup_mask]['aa_v_full'].unique()
    print(f"Found {len(dup_sequences)} sequences with duplicates ({num_in_dup_groups} total entries)")
    
    # Show sample duplicates
    print("\nSample duplicate groups:")
    for seq in dup_sequences[:3]:
        dups = df[df['aa_v_full'] == seq][['id', 'source', 'targets', 'len_v']].head(5)
        print(f"\n  Sequence (first 50 aa): {seq[:50]}...")
        for _, row in dups.iterrows():
            target_preview = str(row['targets'])[:40] + '...' if len(str(row['targets'])) > 40 else row['targets']
            print(f"    - {row['id']} [{row['source']}] targets: {target_preview}")
    
    # Strategy: Keep the entry with best metadata
    # Score: has_target (10 pts) + source priority + id length
    source_priority = {
        'TheraSAbDab': 5,    # Therapeutic - highest value
        'SAbDab_PDB': 4,     # Structural data
        'SAbDab': 4,
        'GenBank': 3,        # Annotated
        'OAS_Camel': 1,      # Repertoire (least metadata)
    }
    
    def score_entry(row):
        score = 0
        # Has target info
        if pd.notna(row['targets']) and str(row['targets']).strip():
            score += 10
        # Source priority
        score += source_priority.get(row['source'], 0)
        return score
    
    df['_score'] = df.apply(score_entry, axis=1)
    
    # Sort by score descending, then keep first occurrence of each sequence
    df_sorted = df.sort_values('_score', ascending=False)
    df_dedup = df_sorted.drop_duplicates(subset=['aa_v_full'], keep='first')
    df_dedup = df_dedup.drop(columns=['_score']).reset_index(drop=True)
    
    removed = initial_count - len(df_dedup)
    print(f"\nRemoved {removed} duplicate entries")
    print(f"Final count: {len(df_dedup)} unique sequences")
    
    # Show what was kept for sample duplicates
    print("\nKept entries for sample duplicates:")
    for seq in dup_sequences[:3]:
        kept = df_dedup[df_dedup['aa_v_full'] == seq][['id', 'source', 'targets']].iloc[0]
        target_preview = str(kept['targets'])[:50] + '...' if len(str(kept['targets'])) > 50 else kept['targets']
        print(f"  {kept['id']} [{kept['source']}] - {target_preview}")
    
    return df_dedup


def merge_with_existing(sabdab_df: pd.DataFrame, existing_csv: str, deduplicate: bool = True) -> pd.DataFrame:
    """Merge SAbDab data with existing unified database."""
    print(f"\nLoading existing database: {existing_csv}")
    existing_df = pd.read_csv(existing_csv)
    print(f"Existing database has {len(existing_df)} sequences")
    
    # Select common columns for merge
    common_cols = ['id', 'source', 'aa_v_full', 'cdr1', 'cdr2', 'cdr3', 'len_v', 'len_cdr3', 'targets']
    
    # Ensure columns exist
    for col in common_cols:
        if col not in existing_df.columns:
            existing_df[col] = ''
        if col not in sabdab_df.columns:
            sabdab_df[col] = ''
    
    # Combine all entries first
    merged_df = pd.concat([
        existing_df[common_cols],
        sabdab_df[common_cols]
    ], ignore_index=True)
    
    print(f"Combined database has {len(merged_df)} total entries (before deduplication)")
    
    # Deduplicate if requested
    if deduplicate:
        merged_df = deduplicate_database(merged_df)
    
    return merged_df


def save_outputs(df: pd.DataFrame, output_prefix: str):
    """Save CSV and NPZ outputs."""
    # Save CSV
    csv_path = f"{output_prefix}.csv"
    df.to_csv(csv_path, index=False)
    print(f"Saved CSV: {csv_path}")
    
    # Save NPZ with arrays for search
    npz_path = f"{output_prefix}.npz"
    
    # Convert to arrays
    ids = np.array(df['id'].values, dtype=object)
    sources = np.array(df['source'].values, dtype=object)
    aa_v_full = np.array(df['aa_v_full'].values, dtype=object)
    cdr1 = np.array(df['cdr1'].values, dtype=object)
    cdr2 = np.array(df['cdr2'].values, dtype=object)
    cdr3 = np.array(df['cdr3'].values, dtype=object)
    len_v = np.array(df['len_v'].values, dtype=np.int32)
    len_cdr3 = np.array(df['len_cdr3'].values, dtype=np.int32)
    targets = np.array(df['targets'].fillna('').values, dtype=object)
    
    np.savez_compressed(
        npz_path,
        ids=ids,
        source=sources,
        aa_v_full=aa_v_full,
        cdr1=cdr1,
        cdr2=cdr2,
        cdr3=cdr3,
        len_v=len_v,
        len_cdr3=len_cdr3,
        targets=targets
    )
    print(f"Saved NPZ: {npz_path}")
    
    # Print statistics
    print("\n=== Database Statistics ===")
    print(f"Total sequences: {len(df)}")
    print(f"\nSource breakdown:")
    print(df['source'].value_counts().to_string())
    print(f"\nWith target/antigen info: {(df['targets'] != '').sum()}")
    print(f"\nV-domain length: {df['len_v'].min()}-{df['len_v'].max()} aa (mean {df['len_v'].mean():.1f})")
    print(f"CDR3 length: {df['len_cdr3'].min()}-{df['len_cdr3'].max()} aa (mean {df['len_cdr3'].mean():.1f})")


def main():
    parser = argparse.ArgumentParser(
        description='Extract VHH sequences from SAbDab PDB files and merge with database'
    )
    parser.add_argument('--tsv', required=True,
                        help='Path to sabdab_nano_summary_all.tsv')
    parser.add_argument('--pdb-dir', required=True,
                        help='Directory containing PDB files (with imgt/ subdirectory)')
    parser.add_argument('--existing-csv', default=None,
                        help='Existing unified database CSV to merge with')
    parser.add_argument('--output-prefix', required=True,
                        help='Output prefix for CSV and NPZ files')
    parser.add_argument('--no-imgt', action='store_true',
                        help='Use main PDB directory instead of imgt/ subdirectory')
    parser.add_argument('--no-dedup', action='store_true',
                        help='Skip deduplication (keep all entries)')
    
    args = parser.parse_args()
    
    # Process SAbDab data
    sabdab_df = process_sabdab_data(
        args.tsv,
        args.pdb_dir,
        use_imgt=not args.no_imgt
    )
    
    if len(sabdab_df) == 0:
        print("No sequences extracted from PDB files!")
        return 1
    
    # Print sample entries
    print("\n=== Sample Extracted Sequences ===")
    for _, row in sabdab_df.head(3).iterrows():
        print(f"\n{row['id']}:")
        print(f"  V-domain: {row['aa_v_full'][:50]}... ({row['len_v']} aa)")
        print(f"  CDR1: {row['cdr1']}")
        print(f"  CDR2: {row['cdr2']}")
        print(f"  CDR3: {row['cdr3']} ({row['len_cdr3']} aa)")
        print(f"  Target: {row['targets'][:80]}..." if len(str(row['targets'])) > 80 else f"  Target: {row['targets']}")
    
    # Merge with existing if provided
    if args.existing_csv:
        final_df = merge_with_existing(sabdab_df, args.existing_csv, deduplicate=not args.no_dedup)
    else:
        # Even without existing, deduplicate the SAbDab data itself
        if not args.no_dedup:
            final_df = deduplicate_database(sabdab_df)
        else:
            final_df = sabdab_df
    
    # Save outputs
    save_outputs(final_df, args.output_prefix)
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
