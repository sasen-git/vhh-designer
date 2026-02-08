#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
OAS/KA-Search Data Diagnostic Script

This script examines your data directory and determines:
1. What format the data is in (NPZ vs TSV/CSV)
2. Whether it contains raw nucleotide sequences
3. What the sequence lengths look like
4. Whether the data is truncated or complete

Run this on your system:
    python diagnose_oas_data.py --data-dir /path/to/your/data
"""

import os
import sys
import glob
import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def check_npz_files(data_dir: str) -> dict:
    """Analyze NPZ files in directory."""
    npz_files = glob.glob(os.path.join(data_dir, "*.npz"))
    
    if not npz_files:
        return {"found": False, "count": 0}
    
    results = {
        "found": True,
        "count": len(npz_files),
        "files": [],
        "sample_analysis": None
    }
    
    # Analyze first NPZ file
    sample_file = npz_files[0]
    try:
        data = np.load(sample_file, allow_pickle=True)
        arrays = list(data.keys())
        
        file_info = {
            "path": sample_file,
            "arrays": arrays,
            "shapes": {k: data[k].shape for k in arrays},
            "dtypes": {k: str(data[k].dtype) for k in arrays}
        }
        
        # Check for numberings array (KA-Search format)
        if "numberings" in arrays:
            numberings = data["numberings"]
            file_info["is_kasearch_format"] = True
            file_info["num_sequences"] = len(numberings)
            
            # Check sequence lengths
            if len(numberings) > 0:
                # Sample some sequences
                sample_indices = np.random.choice(
                    len(numberings), 
                    min(100, len(numberings)), 
                    replace=False
                )
                
                lengths = []
                for idx in sample_indices:
                    row = numberings[idx]
                    # Count non-zero, non-pad positions
                    non_empty = sum(1 for v in row if v not in (0, 124, None))
                    lengths.append(non_empty)
                
                file_info["sample_seq_lengths"] = {
                    "min": min(lengths),
                    "max": max(lengths),
                    "mean": sum(lengths) / len(lengths),
                    "sample_size": len(lengths)
                }
        else:
            file_info["is_kasearch_format"] = False
        
        data.close()
        results["sample_analysis"] = file_info
        
    except Exception as e:
        results["error"] = str(e)
    
    return results


def check_tabular_files(data_dir: str) -> dict:
    """Analyze TSV/CSV files in directory."""
    tsv_files = glob.glob(os.path.join(data_dir, "*.tsv"))
    csv_files = glob.glob(os.path.join(data_dir, "*.csv"))
    
    all_files = tsv_files + csv_files
    
    if not all_files:
        return {"found": False, "count": 0}
    
    results = {
        "found": True,
        "tsv_count": len(tsv_files),
        "csv_count": len(csv_files),
        "files": all_files[:10],  # First 10 files
        "sample_analysis": None
    }
    
    # Analyze first file
    sample_file = all_files[0]
    try:
        # Detect separator
        sep = '\t' if sample_file.endswith('.tsv') else ','
        
        # Read just the header and first few rows
        df_sample = pd.read_csv(sample_file, sep=sep, nrows=100, low_memory=False)
        
        file_info = {
            "path": sample_file,
            "columns": list(df_sample.columns),
            "num_rows_sample": len(df_sample),
        }
        
        # Check for key columns
        has_sequence = 'sequence' in df_sample.columns
        has_sequence_alignment = 'sequence_alignment' in df_sample.columns
        has_sequence_aa = 'sequence_alignment_aa' in df_sample.columns
        
        file_info["has_sequence_column"] = has_sequence
        file_info["has_sequence_alignment"] = has_sequence_alignment
        file_info["has_sequence_aa"] = has_sequence_aa
        
        # Analyze sequence column if present
        if has_sequence:
            seqs = df_sample['sequence'].dropna()
            if len(seqs) > 0:
                lengths = seqs.str.len()
                
                # Check if it's nucleotide (ATCGN) or amino acid
                sample_seq = str(seqs.iloc[0])[:100].upper()
                is_nucleotide = all(c in 'ATCGN-.' for c in sample_seq.replace(' ', ''))
                
                file_info["sequence_info"] = {
                    "type": "nucleotide" if is_nucleotide else "amino_acid",
                    "length_min": int(lengths.min()),
                    "length_max": int(lengths.max()),
                    "length_mean": float(lengths.mean()),
                    "sample_preview": sample_seq[:50] + "..."
                }
        
        # Analyze sequence_alignment_aa if present
        if has_sequence_aa:
            seqs_aa = df_sample['sequence_alignment_aa'].dropna()
            if len(seqs_aa) > 0:
                lengths_aa = seqs_aa.str.len()
                file_info["sequence_aa_info"] = {
                    "length_min": int(lengths_aa.min()),
                    "length_max": int(lengths_aa.max()),
                    "length_mean": float(lengths_aa.mean()),
                }
        
        results["sample_analysis"] = file_info
        
    except Exception as e:
        results["error"] = str(e)
    
    return results


def diagnose_directory(data_dir: str):
    """Run full diagnostics on a data directory."""
    
    print("=" * 70)
    print("OAS/KA-SEARCH DATA DIAGNOSTIC")
    print("=" * 70)
    print(f"\nDirectory: {data_dir}")
    print()
    
    if not os.path.exists(data_dir):
        print(f"ERROR: Directory does not exist!")
        return
    
    # List all files
    all_files = os.listdir(data_dir)
    print(f"Total items in directory: {len(all_files)}")
    
    # Count by extension
    extensions = {}
    for f in all_files:
        ext = os.path.splitext(f)[1].lower()
        extensions[ext] = extensions.get(ext, 0) + 1
    
    print("\nFiles by extension:")
    for ext, count in sorted(extensions.items()):
        print(f"  {ext or '(no extension)'}: {count}")
    
    # Check NPZ files
    print("\n" + "-" * 70)
    print("NPZ FILES (KA-Search pre-processed format)")
    print("-" * 70)
    
    npz_results = check_npz_files(data_dir)
    
    if npz_results["found"]:
        print(f"\nFound {npz_results['count']} NPZ files")
        
        if npz_results.get("sample_analysis"):
            info = npz_results["sample_analysis"]
            print(f"\nSample file: {os.path.basename(info['path'])}")
            print(f"  Arrays: {info['arrays']}")
            
            if info.get("is_kasearch_format"):
                print(f"  ✅ This is KA-Search format")
                print(f"  Number of sequences: {info['num_sequences']:,}")
                
                if info.get("sample_seq_lengths"):
                    sl = info["sample_seq_lengths"]
                    print(f"\n  Sequence lengths (sampled {sl['sample_size']} seqs):")
                    print(f"    Min:  {sl['min']}")
                    print(f"    Max:  {sl['max']}")
                    print(f"    Mean: {sl['mean']:.1f}")
                    
                    if sl['max'] < 100:
                        print(f"\n  ⚠️  WARNING: Sequences appear TRUNCATED!")
                        print(f"     Expected VHH length: 110-130 aa")
                        print(f"     Observed max: {sl['max']} positions")
                    else:
                        print(f"\n  ✅ Sequence lengths look reasonable")
            else:
                print(f"  ❓ Unknown NPZ format (not standard KA-Search)")
    else:
        print("\nNo NPZ files found")
    
    # Check tabular files
    print("\n" + "-" * 70)
    print("TABULAR FILES (TSV/CSV - raw OAS format)")
    print("-" * 70)
    
    tab_results = check_tabular_files(data_dir)
    
    if tab_results["found"]:
        print(f"\nFound {tab_results['tsv_count']} TSV + {tab_results['csv_count']} CSV files")
        
        if tab_results.get("sample_analysis"):
            info = tab_results["sample_analysis"]
            print(f"\nSample file: {os.path.basename(info['path'])}")
            print(f"  Columns ({len(info['columns'])} total):")
            
            # Show key columns
            key_cols = ['sequence', 'sequence_alignment', 'sequence_alignment_aa', 
                       'v_call', 'productive', 'locus', 'cdr3']
            for col in key_cols:
                status = "✅" if col in info['columns'] else "❌"
                print(f"    {status} {col}")
            
            if info.get("has_sequence_column"):
                print(f"\n  'sequence' column analysis:")
                seq_info = info.get("sequence_info", {})
                print(f"    Type: {seq_info.get('type', 'unknown')}")
                print(f"    Length range: {seq_info.get('length_min')} - {seq_info.get('length_max')}")
                print(f"    Mean length: {seq_info.get('length_mean', 0):.1f}")
                print(f"    Preview: {seq_info.get('sample_preview', 'N/A')}")
                
                if seq_info.get('type') == 'nucleotide':
                    print(f"\n  ✅ This file has RAW NUCLEOTIDE sequences!")
                    print(f"     This is what you need for rebuilding the database.")
            
            if info.get("sequence_aa_info"):
                aa_info = info["sequence_aa_info"]
                print(f"\n  'sequence_alignment_aa' column:")
                print(f"    Length range: {aa_info['length_min']} - {aa_info['length_max']}")
                print(f"    Mean length: {aa_info['length_mean']:.1f}")
                
                if aa_info['length_max'] < 100:
                    print(f"    ⚠️  WARNING: AA sequences appear truncated!")
    else:
        print("\nNo TSV/CSV files found")
    
    # Recommendation
    print("\n" + "=" * 70)
    print("RECOMMENDATION")
    print("=" * 70)
    
    if tab_results.get("found") and tab_results.get("sample_analysis", {}).get("has_sequence_column"):
        seq_info = tab_results["sample_analysis"].get("sequence_info", {})
        if seq_info.get("type") == "nucleotide":
            print("""
✅ GOOD NEWS: You have raw nucleotide sequences!

You can rebuild the database using:

    python build_camel_vhh_db.py \\
        --input-tsv "{}" \\
        --output-prefix /home/sasenefrem/KA-Search/camel_vhh_clean \\
        --batch-size 2000
""".format(tab_results["sample_analysis"]["path"]))
        else:
            print("""
⚠️ The 'sequence' column appears to contain amino acids, not nucleotides.
   You may need to find a different data source with raw nucleotide sequences.
""")
    elif npz_results.get("found"):
        print("""
⚠️ This directory only contains NPZ files (pre-processed KA-Search format).
   These files likely have the truncated sequences.

   You need to find the ORIGINAL OAS data with nucleotide sequences.
   
   Options:
   1. Check if raw OAS TSV/CSV files exist elsewhere on your system
   2. Download fresh from OAS: https://opig.stats.ox.ac.uk/webapps/oas/
   
   Look for files with a 'sequence' column containing nucleotide strings (ATCGN...).
""")
    else:
        print("""
❓ No recognizable data files found in this directory.
   Please check the path and try again.
""")


def main():
    parser = argparse.ArgumentParser(
        description="Diagnose OAS/KA-Search data directory"
    )
    parser.add_argument(
        "--data-dir", 
        default="/home/sasenefrem/KA-Search/extracted/oas-paper/oasdb_20230111/Heavy/Camel",
        help="Path to data directory"
    )
    
    args = parser.parse_args()
    diagnose_directory(args.data_dir)


if __name__ == "__main__":
    main()
