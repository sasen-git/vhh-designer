#!/usr/bin/env python3
"""
Shard the VHH database for efficient KA-Search queries.

Creates multiple NPZ shards from a single large database file,
plus an index file for coordinating searches across shards.
"""

import os
import sys
import json
import argparse
import numpy as np
import pandas as pd
from pathlib import Path


def load_database(npz_path: str) -> dict:
    """Load NPZ database into memory."""
    print(f"Loading database: {npz_path}")
    data = np.load(npz_path, allow_pickle=True)
    
    arrays = {}
    for key in data.files:
        arrays[key] = data[key]
        print(f"  {key}: {len(arrays[key])} entries, dtype={arrays[key].dtype}")
    
    return arrays


def create_shards(arrays: dict, output_dir: str, shard_size: int = 2000000, 
                  prefix: str = "shard") -> dict:
    """
    Split arrays into shards and save as separate NPZ files.
    
    Args:
        arrays: Dict of numpy arrays from the database
        output_dir: Directory to save shards
        shard_size: Number of sequences per shard
        prefix: Filename prefix for shards
    
    Returns:
        Index dict with shard information
    """
    os.makedirs(output_dir, exist_ok=True)
    
    total_seqs = len(arrays['ids'])
    num_shards = (total_seqs + shard_size - 1) // shard_size
    
    print(f"\nCreating {num_shards} shards of up to {shard_size:,} sequences each")
    print(f"Total sequences: {total_seqs:,}")
    print(f"Output directory: {output_dir}")
    
    index = {
        'total_sequences': int(total_seqs),
        'shard_size': shard_size,
        'num_shards': num_shards,
        'shards': []
    }
    
    for i in range(num_shards):
        start_idx = i * shard_size
        end_idx = min((i + 1) * shard_size, total_seqs)
        
        shard_name = f"{prefix}_{i:03d}.npz"
        shard_path = os.path.join(output_dir, shard_name)
        
        # Extract slice for each array
        shard_arrays = {}
        for key, arr in arrays.items():
            shard_arrays[key] = arr[start_idx:end_idx]
        
        # Save shard
        np.savez_compressed(shard_path, **shard_arrays)
        
        # Get file size
        file_size = os.path.getsize(shard_path)
        
        # Record shard info
        shard_info = {
            'filename': shard_name,
            'start_idx': int(start_idx),
            'end_idx': int(end_idx),
            'count': int(end_idx - start_idx),
            'size_mb': round(file_size / (1024 * 1024), 2)
        }
        index['shards'].append(shard_info)
        
        print(f"  Created {shard_name}: {end_idx - start_idx:,} sequences ({shard_info['size_mb']} MB)")
    
    return index


def create_source_based_shards(arrays: dict, output_dir: str, prefix: str = "shard") -> dict:
    """
    Create shards based on source type for more targeted searching.
    Splits into: annotated (with targets), oas_camel, indi_ngs
    """
    os.makedirs(output_dir, exist_ok=True)
    
    sources = arrays['source']
    total_seqs = len(sources)
    
    print(f"\nCreating source-based shards")
    print(f"Total sequences: {total_seqs:,}")
    
    # Define source groups
    annotated_sources = {'INDI_patent', 'INDI_structure', 'INDI_manual', 
                         'SAbDab_PDB', 'SAbDab', 'GenBank', 'TheraSAbDab', 'INDI_genbank'}
    
    # Create masks
    annotated_mask = np.isin(sources, list(annotated_sources))
    oas_mask = sources == 'OAS_Camel'
    indi_ngs_mask = sources == 'INDI_NGS'
    
    index = {
        'total_sequences': int(total_seqs),
        'shard_type': 'source_based',
        'shards': []
    }
    
    # Helper to save a shard
    def save_shard(mask, name, description):
        indices = np.where(mask)[0]
        if len(indices) == 0:
            return
        
        shard_arrays = {}
        for key, arr in arrays.items():
            shard_arrays[key] = arr[indices]
        
        shard_path = os.path.join(output_dir, f"{name}.npz")
        np.savez_compressed(shard_path, **shard_arrays)
        
        file_size = os.path.getsize(shard_path)
        
        shard_info = {
            'filename': f"{name}.npz",
            'description': description,
            'count': int(len(indices)),
            'size_mb': round(file_size / (1024 * 1024), 2),
            'sources': list(np.unique(shard_arrays['source']))
        }
        index['shards'].append(shard_info)
        
        print(f"  Created {name}.npz: {len(indices):,} sequences ({shard_info['size_mb']} MB)")
        print(f"    Sources: {shard_info['sources']}")
    
    # Save annotated shard (patents, structures, curated - with target info)
    save_shard(annotated_mask, f"{prefix}_annotated", 
               "Sequences with annotations (patents, structures, curated)")
    
    # Save OAS Camel shard
    save_shard(oas_mask, f"{prefix}_oas_camel",
               "OAS Camel repertoire sequences")
    
    # For INDI_NGS, split into sub-shards of 2M each
    indi_ngs_indices = np.where(indi_ngs_mask)[0]
    if len(indi_ngs_indices) > 0:
        shard_size = 2000000
        num_sub_shards = (len(indi_ngs_indices) + shard_size - 1) // shard_size
        
        for i in range(num_sub_shards):
            start = i * shard_size
            end = min((i + 1) * shard_size, len(indi_ngs_indices))
            sub_indices = indi_ngs_indices[start:end]
            
            shard_arrays = {}
            for key, arr in arrays.items():
                shard_arrays[key] = arr[sub_indices]
            
            shard_name = f"{prefix}_indi_ngs_{i:03d}"
            shard_path = os.path.join(output_dir, f"{shard_name}.npz")
            np.savez_compressed(shard_path, **shard_arrays)
            
            file_size = os.path.getsize(shard_path)
            
            shard_info = {
                'filename': f"{shard_name}.npz",
                'description': f"INDI NGS repertoire (part {i+1}/{num_sub_shards})",
                'count': int(len(sub_indices)),
                'size_mb': round(file_size / (1024 * 1024), 2),
                'sources': ['INDI_NGS']
            }
            index['shards'].append(shard_info)
            
            print(f"  Created {shard_name}.npz: {len(sub_indices):,} sequences ({shard_info['size_mb']} MB)")
    
    return index


def save_index(index: dict, output_dir: str):
    """Save the shard index as JSON."""
    index_path = os.path.join(output_dir, "shard_index.json")
    with open(index_path, 'w') as f:
        json.dump(index, f, indent=2)
    print(f"\nSaved index: {index_path}")


def main():
    parser = argparse.ArgumentParser(
        description='Shard VHH database for KA-Search'
    )
    parser.add_argument('--input', required=True,
                        help='Input NPZ file (VHH_db_final.npz)')
    parser.add_argument('--output-dir', required=True,
                        help='Output directory for shards')
    parser.add_argument('--shard-size', type=int, default=2000000,
                        help='Sequences per shard (default: 2M)')
    parser.add_argument('--prefix', default='vhh',
                        help='Shard filename prefix (default: vhh)')
    parser.add_argument('--by-source', action='store_true',
                        help='Create source-based shards instead of size-based')
    
    args = parser.parse_args()
    
    # Load database
    arrays = load_database(args.input)
    
    # Create shards
    if args.by_source:
        index = create_source_based_shards(arrays, args.output_dir, args.prefix)
    else:
        index = create_shards(arrays, args.output_dir, args.shard_size, args.prefix)
    
    # Save index
    save_index(index, args.output_dir)
    
    # Summary
    print("\n" + "="*60)
    print("SHARDING COMPLETE")
    print("="*60)
    print(f"Total shards: {len(index['shards'])}")
    print(f"Total sequences: {index['total_sequences']:,}")
    total_size = sum(s['size_mb'] for s in index['shards'])
    print(f"Total size: {total_size:.1f} MB")
    
    print(f"\nShards created in: {args.output_dir}/")
    for shard in index['shards']:
        print(f"  {shard['filename']}: {shard['count']:,} seqs ({shard['size_mb']} MB)")
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
