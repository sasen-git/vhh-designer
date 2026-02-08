#!/usr/bin/env python3
"""
Process INDI database and merge with existing VHH unified database.

Creates a comprehensive VHH database with:
- Detailed source tracking (OAS_Camel, INDI_patent, INDI_structure, SAbDab_PDB, etc.)
- Target/binding partner information where available
- Patent numbers for patent-derived sequences
- Deduplication keeping entries with most metadata

INDI files:
- patent_sequence.tsv + patent_meta.tsv (targets from biomolecules)
- structure_sequence.tsv + structure_meta.tsv (PDB info)
- abgenbank_sequence.tsv (no metadata file)
- manual_sequence.tsv + manual_meta.tsv (curated entries)
- ngs_sequence.tsv (SKIP - NGS repertoire data, no functional annotation)
"""

import os
import sys
import argparse
import pandas as pd
import numpy as np
from pathlib import Path


def load_indi_sequences(seq_file: str, meta_file: str = None, source_name: str = 'INDI') -> pd.DataFrame:
    """
    Load INDI sequence file and optionally join with metadata.
    
    Args:
        seq_file: Path to sequence TSV (sequence, CDR1, CDR2, CDR3, IDs)
        meta_file: Path to metadata TSV (optional)
        source_name: Source label (e.g., 'INDI_patent')
    
    Returns:
        DataFrame with sequences and metadata
    """
    print(f"\nLoading {seq_file}...")
    
    # Load sequences
    seq_df = pd.read_csv(seq_file, sep='\t')
    print(f"  Loaded {len(seq_df)} sequences")
    
    # Rename columns for consistency
    seq_df = seq_df.rename(columns={
        'sequence': 'aa_v_full',
        'CDR1': 'cdr1',
        'CDR2': 'cdr2', 
        'CDR3': 'cdr3',
        'Comma-separated db-specific IDS': 'db_ids'
    })
    
    # Add source
    seq_df['source'] = source_name
    
    # Calculate lengths
    seq_df['len_v'] = seq_df['aa_v_full'].str.len()
    seq_df['len_cdr3'] = seq_df['cdr3'].fillna('').str.len()
    
    # Create ID from first db_id
    seq_df['id'] = seq_df['db_ids'].astype(str).str.split(',').str[0]
    seq_df['id'] = source_name + '_' + seq_df['id'].astype(str)
    
    # Initialize target columns
    seq_df['targets'] = ''
    seq_df['patent_id'] = ''
    seq_df['patent_title'] = ''
    seq_df['organism'] = ''
    seq_df['pdb_id'] = ''
    seq_df['reference'] = ''
    
    # Join with metadata if available
    if meta_file and os.path.exists(meta_file):
        print(f"  Joining with metadata: {meta_file}")
        meta_df = pd.read_csv(meta_file, sep='\t')
        meta_df['id'] = meta_df['id'].astype(str)
        
        # Create lookup dict for metadata
        meta_dict = {}
        for _, row in meta_df.iterrows():
            meta_dict[str(row['id'])] = row.to_dict()
        
        # Apply metadata to sequences
        def apply_metadata(row):
            # Get all IDs for this sequence
            ids = str(row['db_ids']).split(',')
            
            targets = []
            patent_ids = []
            patent_titles = []
            organisms = []
            pdb_ids = []
            references = []
            
            for id_str in ids:
                id_str = id_str.strip()
                if id_str in meta_dict:
                    meta = meta_dict[id_str]
                    
                    # Extract targets based on source type
                    if 'comma-separated identified biomolecules' in meta:
                        biomol = meta['comma-separated identified biomolecules']
                        if pd.notna(biomol) and biomol:
                            targets.append(str(biomol))
                    
                    if 'patent name' in meta:
                        pat = meta['patent name']
                        if pd.notna(pat) and pat:
                            patent_ids.append(str(pat))
                    
                    if 'patent title' in meta:
                        title = meta['patent title']
                        if pd.notna(title) and title:
                            patent_titles.append(str(title).strip('"'))
                    
                    if 'organism' in meta:
                        org = meta['organism']
                        if pd.notna(org) and org:
                            organisms.append(str(org))
                    
                    if 'chain title' in meta:
                        chain = meta['chain title']
                        if pd.notna(chain) and chain:
                            targets.append(str(chain))
                    
                    if 'pdb title' in meta:
                        pdb_title = meta['pdb title']
                        if pd.notna(pdb_title) and pdb_title:
                            references.append(str(pdb_title))
                    
                    if 'title' in meta:
                        title = meta['title']
                        if pd.notna(title) and title:
                            references.append(str(title))
                    
                    if 'url' in meta:
                        url = meta['url']
                        if pd.notna(url) and url:
                            # Extract PDB ID from URL if present
                            if 'rcsb.org' in str(url):
                                pdb = str(url).split('/')[-1].upper()
                                pdb_ids.append(pdb)
            
            # Combine and deduplicate
            row['targets'] = '; '.join(list(dict.fromkeys(targets)))[:500]  # Limit length
            row['patent_id'] = '; '.join(list(dict.fromkeys(patent_ids)))
            row['patent_title'] = '; '.join(list(dict.fromkeys(patent_titles)))[:300]
            row['organism'] = '; '.join(list(dict.fromkeys(organisms)))
            row['pdb_id'] = '; '.join(list(dict.fromkeys(pdb_ids)))
            row['reference'] = '; '.join(list(dict.fromkeys(references)))[:300]
            
            return row
        
        seq_df = seq_df.apply(apply_metadata, axis=1)
        
        # Count how many got metadata
        has_targets = (seq_df['targets'] != '').sum()
        print(f"  Sequences with target info: {has_targets}")
    
    return seq_df


def process_all_indi(indi_dir: str, skip_ngs: bool = False) -> pd.DataFrame:
    """Process all INDI source files."""
    all_dfs = []
    
    # Patent sequences (highest value - has targets!)
    patent_seq = os.path.join(indi_dir, 'patent_sequence.tsv')
    patent_meta = os.path.join(indi_dir, 'patent_meta.tsv')
    if os.path.exists(patent_seq):
        df = load_indi_sequences(patent_seq, patent_meta, 'INDI_patent')
        all_dfs.append(df)
    
    # Structure sequences
    struct_seq = os.path.join(indi_dir, 'structure_sequence.tsv')
    struct_meta = os.path.join(indi_dir, 'structure_meta.tsv')
    if os.path.exists(struct_seq):
        df = load_indi_sequences(struct_seq, struct_meta, 'INDI_structure')
        all_dfs.append(df)
    
    # Manual sequences
    manual_seq = os.path.join(indi_dir, 'manual_sequence.tsv')
    manual_meta = os.path.join(indi_dir, 'manual_meta.tsv')
    if os.path.exists(manual_seq):
        df = load_indi_sequences(manual_seq, manual_meta, 'INDI_manual')
        all_dfs.append(df)
    
    # GenBank sequences (no metadata file)
    genbank_seq = os.path.join(indi_dir, 'abgenbank_sequence.tsv')
    if os.path.exists(genbank_seq):
        df = load_indi_sequences(genbank_seq, None, 'INDI_genbank')
        all_dfs.append(df)
    
    # NGS sequences (11M+ sequences - camelid repertoires from various species)
    if not skip_ngs:
        ngs_seq = os.path.join(indi_dir, 'ngs_sequence.tsv')
        if os.path.exists(ngs_seq):
            print(f"\nLoading NGS data (this may take a few minutes)...")
            df = load_indi_sequences(ngs_seq, None, 'INDI_NGS')
            all_dfs.append(df)
    else:
        print("\nSkipping NGS data (use without --skip-ngs to include 11M+ sequences)")
    
    if not all_dfs:
        print("No INDI files found!")
        return pd.DataFrame()
    
    combined = pd.concat(all_dfs, ignore_index=True)
    print(f"\nTotal INDI sequences: {len(combined)}")
    print(f"Source breakdown:")
    print(combined['source'].value_counts().to_string())
    
    return combined


def load_existing_database(csv_path: str) -> pd.DataFrame:
    """Load existing unified database."""
    print(f"\nLoading existing database: {csv_path}")
    df = pd.read_csv(csv_path, low_memory=False)
    print(f"  Loaded {len(df)} sequences")
    
    # Ensure required columns exist
    for col in ['id', 'source', 'aa_v_full', 'cdr1', 'cdr2', 'cdr3', 
                'len_v', 'len_cdr3', 'targets', 'patent_id', 'patent_title',
                'organism', 'pdb_id', 'reference']:
        if col not in df.columns:
            df[col] = ''
    
    return df


def deduplicate_database(df: pd.DataFrame) -> pd.DataFrame:
    """
    Remove duplicate sequences, keeping the entry with the most metadata.
    Priority: entries with target info > patent info > structural info > repertoire
    """
    print(f"\n=== Deduplicating Database ===")
    initial_count = len(df)
    
    # Score each entry based on metadata richness
    source_priority = {
        'TheraSAbDab': 10,       # Therapeutic - highest value
        'INDI_patent': 9,        # Patent with targets
        'SAbDab_PDB': 8,         # Structural with antigen
        'INDI_structure': 7,     # Structural
        'SAbDab': 7,
        'INDI_manual': 6,        # Curated
        'GenBank': 5,            # Annotated
        'INDI_genbank': 4,
        'OAS_Camel': 1,          # Repertoire (least metadata)
    }
    
    def score_entry(row):
        score = 0
        # Has target info (most important)
        if pd.notna(row.get('targets')) and str(row.get('targets', '')).strip():
            score += 20
        # Has patent info
        if pd.notna(row.get('patent_id')) and str(row.get('patent_id', '')).strip():
            score += 10
        # Has organism info
        if pd.notna(row.get('organism')) and str(row.get('organism', '')).strip():
            score += 5
        # Has PDB info
        if pd.notna(row.get('pdb_id')) and str(row.get('pdb_id', '')).strip():
            score += 5
        # Source priority
        score += source_priority.get(row.get('source', ''), 0)
        return score
    
    df['_score'] = df.apply(score_entry, axis=1)
    
    # Find duplicates
    dup_mask = df.duplicated(subset=['aa_v_full'], keep=False)
    num_in_dup_groups = dup_mask.sum()
    
    if num_in_dup_groups > 0:
        dup_sequences = df[dup_mask]['aa_v_full'].nunique()
        print(f"Found {dup_sequences} sequences with duplicates ({num_in_dup_groups} total entries)")
        
        # Show sample duplicates
        print("\nSample duplicate groups:")
        sample_seqs = df[dup_mask]['aa_v_full'].unique()[:3]
        for seq in sample_seqs:
            dups = df[df['aa_v_full'] == seq][['id', 'source', 'targets', '_score']].head(4)
            print(f"\n  Sequence: {seq[:50]}...")
            for _, row in dups.iterrows():
                target_preview = str(row['targets'])[:40] if row['targets'] else 'N/A'
                print(f"    [{row['source']}] score={row['_score']}: {target_preview}")
    
    # Sort by score descending, keep first
    df_sorted = df.sort_values('_score', ascending=False)
    df_dedup = df_sorted.drop_duplicates(subset=['aa_v_full'], keep='first')
    df_dedup = df_dedup.drop(columns=['_score']).reset_index(drop=True)
    
    removed = initial_count - len(df_dedup)
    print(f"\nRemoved {removed} duplicate entries")
    print(f"Final count: {len(df_dedup)} unique sequences")
    
    return df_dedup


def save_final_database(df: pd.DataFrame, output_prefix: str):
    """Save the final database as CSV and NPZ."""
    
    # Select and order columns
    output_cols = [
        'id', 'source', 'aa_v_full', 'cdr1', 'cdr2', 'cdr3',
        'len_v', 'len_cdr3', 'targets', 'patent_id', 'patent_title',
        'organism', 'pdb_id', 'reference'
    ]
    
    # Ensure all columns exist
    for col in output_cols:
        if col not in df.columns:
            df[col] = ''
    
    df_out = df[output_cols].copy()
    
    # Fill NaN with empty string
    df_out = df_out.fillna('')
    
    # Save CSV
    csv_path = f"{output_prefix}.csv"
    df_out.to_csv(csv_path, index=False)
    print(f"\nSaved CSV: {csv_path}")
    
    # Save NPZ
    npz_path = f"{output_prefix}.npz"
    np.savez_compressed(
        npz_path,
        ids=np.array(df_out['id'].values, dtype=object),
        source=np.array(df_out['source'].values, dtype=object),
        aa_v_full=np.array(df_out['aa_v_full'].values, dtype=object),
        cdr1=np.array(df_out['cdr1'].values, dtype=object),
        cdr2=np.array(df_out['cdr2'].values, dtype=object),
        cdr3=np.array(df_out['cdr3'].values, dtype=object),
        len_v=np.array(df_out['len_v'].values, dtype=np.int32),
        len_cdr3=np.array(df_out['len_cdr3'].values, dtype=np.int32),
        targets=np.array(df_out['targets'].values, dtype=object),
        patent_id=np.array(df_out['patent_id'].values, dtype=object),
        patent_title=np.array(df_out['patent_title'].values, dtype=object),
        organism=np.array(df_out['organism'].values, dtype=object),
        pdb_id=np.array(df_out['pdb_id'].values, dtype=object),
        reference=np.array(df_out['reference'].values, dtype=object),
    )
    print(f"Saved NPZ: {npz_path}")
    
    # Print comprehensive statistics
    print("\n" + "="*60)
    print("FINAL DATABASE STATISTICS")
    print("="*60)
    
    print(f"\nTotal unique sequences: {len(df_out)}")
    
    print(f"\n--- Source Breakdown ---")
    print(df_out['source'].value_counts().to_string())
    
    print(f"\n--- Metadata Coverage ---")
    has_targets = (df_out['targets'] != '').sum()
    has_patent = (df_out['patent_id'] != '').sum()
    has_organism = (df_out['organism'] != '').sum()
    has_pdb = (df_out['pdb_id'] != '').sum()
    
    print(f"With target/binding info: {has_targets} ({100*has_targets/len(df_out):.2f}%)")
    print(f"With patent ID: {has_patent} ({100*has_patent/len(df_out):.2f}%)")
    print(f"With organism info: {has_organism} ({100*has_organism/len(df_out):.2f}%)")
    print(f"With PDB ID: {has_pdb} ({100*has_pdb/len(df_out):.2f}%)")
    
    print(f"\n--- Sequence Statistics ---")
    print(f"V-domain length: {df_out['len_v'].min()}-{df_out['len_v'].max()} aa (mean {df_out['len_v'].mean():.1f})")
    print(f"CDR3 length: {df_out['len_cdr3'].min()}-{df_out['len_cdr3'].max()} aa (mean {df_out['len_cdr3'].mean():.1f})")
    
    # Sample entries with rich metadata
    print(f"\n--- Sample Entries with Target Info ---")
    rich_entries = df_out[df_out['targets'] != ''].head(5)
    for _, row in rich_entries.iterrows():
        print(f"\n  {row['id']} [{row['source']}]")
        print(f"    CDR3: {row['cdr3']}")
        print(f"    Target: {str(row['targets'])[:80]}{'...' if len(str(row['targets'])) > 80 else ''}")
        if row['patent_id']:
            print(f"    Patent: {row['patent_id']}")


def main():
    parser = argparse.ArgumentParser(
        description='Process INDI database and merge with existing VHH database'
    )
    parser.add_argument('--indi-dir', required=True,
                        help='Directory containing INDI TSV files')
    parser.add_argument('--existing-csv', default=None,
                        help='Existing unified database CSV to merge with')
    parser.add_argument('--output-prefix', required=True,
                        help='Output prefix for final CSV and NPZ files')
    parser.add_argument('--no-dedup', action='store_true',
                        help='Skip deduplication')
    parser.add_argument('--skip-ngs', action='store_true',
                        help='Skip NGS data (11M sequences) for faster processing')
    
    args = parser.parse_args()
    
    # Process INDI data
    indi_df = process_all_indi(args.indi_dir, skip_ngs=args.skip_ngs)
    
    if len(indi_df) == 0:
        print("No INDI sequences found!")
        return 1
    
    # Show sample INDI entries
    print("\n=== Sample INDI Sequences ===")
    for source in indi_df['source'].unique():
        sample = indi_df[indi_df['source'] == source].head(1)
        for _, row in sample.iterrows():
            print(f"\n{row['id']} [{row['source']}]:")
            print(f"  Sequence: {row['aa_v_full'][:50]}... ({row['len_v']} aa)")
            print(f"  CDR3: {row['cdr3']} ({row['len_cdr3']} aa)")
            if row['targets']:
                print(f"  Targets: {str(row['targets'])[:60]}...")
            if row['patent_id']:
                print(f"  Patent: {row['patent_id']}")
    
    # Merge with existing if provided
    if args.existing_csv:
        existing_df = load_existing_database(args.existing_csv)
        
        # Standardize columns in existing
        col_mapping = {
            'id': 'id', 'source': 'source', 'aa_v_full': 'aa_v_full',
            'cdr1': 'cdr1', 'cdr2': 'cdr2', 'cdr3': 'cdr3',
            'len_v': 'len_v', 'len_cdr3': 'len_cdr3', 'targets': 'targets'
        }
        
        combined_df = pd.concat([existing_df, indi_df], ignore_index=True)
        print(f"\nCombined database: {len(combined_df)} entries")
    else:
        combined_df = indi_df
    
    # Deduplicate
    if not args.no_dedup:
        final_df = deduplicate_database(combined_df)
    else:
        final_df = combined_df
    
    # Save final database
    save_final_database(final_df, args.output_prefix)
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
