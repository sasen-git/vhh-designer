#!/usr/bin/env python3
"""
VHH CDR-Framework Correlation Analyzer
======================================
Processes NPZ shards to find correlations between:
- Single amino acid positions (CDR pos X ↔ FW pos Y)
- Triplet motifs (CDR triplet ↔ FW motif)
- Lengths (CDR length ↔ FW length)

Usage:
    python analyze_cdr_framework_correlations.py /path/to/shards/*.npz
    
    # Or specify output directory:
    python analyze_cdr_framework_correlations.py --output ./results /path/to/shards/*.npz

Output:
    - correlation_results.pkl: Raw counts for all correlations
    - correlation_summary.txt: Human-readable summary
    - (Optional) Visualizations if --visualize flag is used
"""

import numpy as np
from collections import defaultdict
import pickle
import sys
import os
import time
import argparse
from glob import glob

# ============================================================
# EXTRACTION FUNCTION
# ============================================================

def extract_frameworks(full_seq, cdr1, cdr2, cdr3):
    """Extract framework regions from full sequence given CDRs"""
    if not full_seq or not cdr1 or not cdr2 or not cdr3:
        return None
    
    full_seq, cdr1, cdr2, cdr3 = str(full_seq), str(cdr1), str(cdr2), str(cdr3)
    
    cdr1_start = full_seq.find(cdr1)
    cdr2_start = full_seq.find(cdr2)
    cdr3_start = full_seq.find(cdr3)
    
    if cdr1_start == -1 or cdr2_start == -1 or cdr3_start == -1:
        return None
    if not (cdr1_start < cdr2_start < cdr3_start):
        return None
    
    fr1 = full_seq[:cdr1_start]
    fr2 = full_seq[cdr1_start + len(cdr1):cdr2_start]
    fr3 = full_seq[cdr2_start + len(cdr2):cdr3_start]
    fr4 = full_seq[cdr3_start + len(cdr3):]
    
    # Filter extremes
    if len(fr1) < 10 or len(fr2) < 10 or len(fr3) < 20 or len(fr4) < 5:
        return None
    if len(fr1) > 40 or len(fr2) > 25 or len(fr3) > 50 or len(fr4) > 20:
        return None
    
    return {
        'FR1': fr1, 'FR2': fr2, 'FR3': fr3, 'FR4': fr4,
        'cdr1': cdr1, 'cdr2': cdr2, 'cdr3': cdr3
    }

# ============================================================
# SHARD PROCESSING
# ============================================================

def process_shard(npz_path, verbose=True):
    """Process a single NPZ shard and return count dictionaries"""
    
    if verbose:
        print(f"  Processing: {os.path.basename(npz_path)}", end=" ", flush=True)
    start_time = time.time()
    
    try:
        data = np.load(npz_path, allow_pickle=True)
        full_seqs = data['aa_v_full']
        cdr1s = data['cdr1']
        cdr2s = data['cdr2']
        cdr3s = data['cdr3']
    except Exception as e:
        print(f"ERROR: {e}")
        return None
    
    n_total = len(full_seqs)
    
    # Initialize counts
    # Single AA: (cdr_name, cdr_pos, cdr_aa, fw_name, fw_pos, fw_aa) -> count
    single_aa_counts = defaultdict(int)
    
    # Triplets: (junction_type, cdr_triplet, fw_motif) -> count
    triplet_counts = defaultdict(int)
    
    # Lengths: (cdr_name, cdr_len, fw_name, fw_len) -> count
    length_counts = defaultdict(int)
    
    # Framework sequences for consensus
    fw_sequences = {'FR1': [], 'FR2': [], 'FR3': [], 'FR4': []}
    
    n_valid = 0
    for i in range(n_total):
        fw = extract_frameworks(full_seqs[i], cdr1s[i], cdr2s[i], cdr3s[i])
        if fw is None:
            continue
        
        n_valid += 1
        
        # Store framework sequences (sample for memory efficiency)
        if n_valid <= 10000:
            for fw_name in ['FR1', 'FR2', 'FR3', 'FR4']:
                fw_sequences[fw_name].append(fw[fw_name])
        
        # === SINGLE AA CORRELATIONS ===
        # Analyze boundary positions (most important for grafting)
        cdr_positions = {
            'cdr1': [0, 1, 2, -3, -2, -1],
            'cdr2': [0, 1, 2, -3, -2, -1],
            'cdr3': [0, 1, 2, -3, -2, -1]
        }
        
        fw_positions = {
            'FR1': list(range(-10, 0)),  # Last 10
            'FR2': list(range(10)) + list(range(-10, 0)),  # First 10 + last 10
            'FR3': list(range(10)) + list(range(-10, 0)),
            'FR4': list(range(10))  # First 10
        }
        
        for cdr_name in ['cdr1', 'cdr2', 'cdr3']:
            cdr_seq = fw[cdr_name]
            for cdr_pos in cdr_positions[cdr_name]:
                try:
                    cdr_aa = cdr_seq[cdr_pos]
                except IndexError:
                    continue
                
                for fw_name in ['FR1', 'FR2', 'FR3', 'FR4']:
                    fw_seq = fw[fw_name]
                    for fw_pos in fw_positions[fw_name]:
                        try:
                            fw_aa = fw_seq[fw_pos]
                        except IndexError:
                            continue
                        
                        key = (cdr_name, cdr_pos, cdr_aa, fw_name, fw_pos, fw_aa)
                        single_aa_counts[key] += 1
        
        # === TRIPLET ASSOCIATIONS ===
        junctions = [
            ('CDR1_start', 'FR1_end', fw['cdr1'][:3] if len(fw['cdr1']) >= 3 else '', fw['FR1'][-5:] if len(fw['FR1']) >= 5 else ''),
            ('CDR1_end', 'FR2_start', fw['cdr1'][-3:] if len(fw['cdr1']) >= 3 else '', fw['FR2'][:5] if len(fw['FR2']) >= 5 else ''),
            ('CDR2_start', 'FR2_end', fw['cdr2'][:3] if len(fw['cdr2']) >= 3 else '', fw['FR2'][-5:] if len(fw['FR2']) >= 5 else ''),
            ('CDR2_end', 'FR3_start', fw['cdr2'][-3:] if len(fw['cdr2']) >= 3 else '', fw['FR3'][:5] if len(fw['FR3']) >= 5 else ''),
            ('CDR3_start', 'FR3_end', fw['cdr3'][:3] if len(fw['cdr3']) >= 3 else '', fw['FR3'][-5:] if len(fw['FR3']) >= 5 else ''),
            ('CDR3_end', 'FR4_start', fw['cdr3'][-3:] if len(fw['cdr3']) >= 3 else '', fw['FR4'][:5] if len(fw['FR4']) >= 5 else ''),
        ]
        
        for junction_from, junction_to, cdr_trip, fw_motif in junctions:
            if len(cdr_trip) == 3 and len(fw_motif) == 5:
                triplet_counts[(junction_from, junction_to, cdr_trip, fw_motif)] += 1
        
        # === LENGTH CORRELATIONS ===
        for cdr_name in ['cdr1', 'cdr2', 'cdr3']:
            cdr_len = len(fw[cdr_name])
            for fw_name in ['FR1', 'FR2', 'FR3', 'FR4']:
                fw_len = len(fw[fw_name])
                length_counts[(cdr_name, cdr_len, fw_name, fw_len)] += 1
    
    elapsed = time.time() - start_time
    if verbose:
        print(f"→ {n_valid:,}/{n_total:,} valid ({elapsed:.1f}s)")
    
    return {
        'single_aa': dict(single_aa_counts),
        'triplets': dict(triplet_counts),
        'lengths': dict(length_counts),
        'fw_sequences': {k: v[:1000] for k, v in fw_sequences.items()},  # Keep sample
        'n_valid': n_valid,
        'n_total': n_total
    }

# ============================================================
# MERGE RESULTS
# ============================================================

def merge_results(results_list):
    """Merge results from multiple shards"""
    merged = {
        'single_aa': defaultdict(int),
        'triplets': defaultdict(int),
        'lengths': defaultdict(int),
        'fw_sequences': {'FR1': [], 'FR2': [], 'FR3': [], 'FR4': []},
        'n_valid': 0,
        'n_total': 0
    }
    
    for results in results_list:
        if results is None:
            continue
        
        merged['n_valid'] += results['n_valid']
        merged['n_total'] += results['n_total']
        
        for key, count in results['single_aa'].items():
            merged['single_aa'][key] += count
        
        for key, count in results['triplets'].items():
            merged['triplets'][key] += count
        
        for key, count in results['lengths'].items():
            merged['lengths'][key] += count
        
        # Keep sample of FW sequences
        for fw_name in ['FR1', 'FR2', 'FR3', 'FR4']:
            if len(merged['fw_sequences'][fw_name]) < 5000:
                merged['fw_sequences'][fw_name].extend(results['fw_sequences'].get(fw_name, []))
    
    return {
        'single_aa': dict(merged['single_aa']),
        'triplets': dict(merged['triplets']),
        'lengths': dict(merged['lengths']),
        'fw_sequences': merged['fw_sequences'],
        'n_valid': merged['n_valid'],
        'n_total': merged['n_total']
    }

# ============================================================
# COMPUTE STATISTICS
# ============================================================

def compute_statistics(merged_results, min_n=50):
    """Compute confidence percentages and organize results"""
    
    stats = {
        'single_aa_associations': [],
        'triplet_associations': [],
        'length_correlations': [],
        'summary': {}
    }
    
    # === SINGLE AA STATISTICS ===
    # Group by (cdr_name, cdr_pos, cdr_aa, fw_name, fw_pos)
    single_aa = merged_results['single_aa']
    
    cdr_fw_pairs = defaultdict(lambda: defaultdict(int))
    for key, count in single_aa.items():
        cdr_name, cdr_pos, cdr_aa, fw_name, fw_pos, fw_aa = key
        pair_key = (cdr_name, cdr_pos, cdr_aa, fw_name, fw_pos)
        cdr_fw_pairs[pair_key][fw_aa] += count
    
    for pair_key, fw_aa_counts in cdr_fw_pairs.items():
        cdr_name, cdr_pos, cdr_aa, fw_name, fw_pos = pair_key
        total = sum(fw_aa_counts.values())
        
        if total < min_n:
            continue
        
        top_fw_aa = max(fw_aa_counts, key=fw_aa_counts.get)
        confidence = 100 * fw_aa_counts[top_fw_aa] / total
        
        stats['single_aa_associations'].append({
            'cdr': cdr_name,
            'cdr_pos': cdr_pos,
            'cdr_aa': cdr_aa,
            'fw': fw_name,
            'fw_pos': fw_pos,
            'predicted_fw_aa': top_fw_aa,
            'confidence': confidence,
            'n': total,
            'distribution': dict(fw_aa_counts)
        })
    
    # === TRIPLET STATISTICS ===
    triplets = merged_results['triplets']
    
    junction_pairs = defaultdict(lambda: defaultdict(int))
    for key, count in triplets.items():
        junction_from, junction_to, cdr_trip, fw_motif = key
        pair_key = (junction_from, junction_to, cdr_trip)
        junction_pairs[pair_key][fw_motif] += count
    
    for pair_key, fw_motif_counts in junction_pairs.items():
        junction_from, junction_to, cdr_trip = pair_key
        total = sum(fw_motif_counts.values())
        
        if total < min_n:
            continue
        
        top_fw_motif = max(fw_motif_counts, key=fw_motif_counts.get)
        confidence = 100 * fw_motif_counts[top_fw_motif] / total
        
        stats['triplet_associations'].append({
            'junction_from': junction_from,
            'junction_to': junction_to,
            'cdr_triplet': cdr_trip,
            'predicted_fw_motif': top_fw_motif,
            'confidence': confidence,
            'n': total,
            'distribution': dict(fw_motif_counts)
        })
    
    # === LENGTH STATISTICS ===
    lengths = merged_results['lengths']
    
    # Group by (cdr_name, cdr_len, fw_name)
    length_pairs = defaultdict(lambda: defaultdict(int))
    for key, count in lengths.items():
        cdr_name, cdr_len, fw_name, fw_len = key
        pair_key = (cdr_name, cdr_len, fw_name)
        length_pairs[pair_key][fw_len] += count
    
    for pair_key, fw_len_counts in length_pairs.items():
        cdr_name, cdr_len, fw_name = pair_key
        total = sum(fw_len_counts.values())
        
        if total < min_n:
            continue
        
        top_fw_len = max(fw_len_counts, key=fw_len_counts.get)
        confidence = 100 * fw_len_counts[top_fw_len] / total
        
        stats['length_correlations'].append({
            'cdr': cdr_name,
            'cdr_len': cdr_len,
            'fw': fw_name,
            'predicted_fw_len': top_fw_len,
            'confidence': confidence,
            'n': total,
            'distribution': dict(fw_len_counts)
        })
    
    # Summary
    stats['summary'] = {
        'n_sequences': merged_results['n_valid'],
        'n_single_aa_associations': len(stats['single_aa_associations']),
        'n_triplet_associations': len(stats['triplet_associations']),
        'n_length_correlations': len(stats['length_correlations'])
    }
    
    return stats

# ============================================================
# WRITE SUMMARY
# ============================================================

def write_summary(stats, output_path):
    """Write human-readable summary"""
    
    with open(output_path, 'w') as f:
        f.write("="*70 + "\n")
        f.write("CDR-FRAMEWORK CORRELATION ANALYSIS RESULTS\n")
        f.write("="*70 + "\n\n")
        
        f.write(f"Total sequences analyzed: {stats['summary']['n_sequences']:,}\n")
        f.write(f"Single AA associations: {stats['summary']['n_single_aa_associations']:,}\n")
        f.write(f"Triplet associations: {stats['summary']['n_triplet_associations']:,}\n")
        f.write(f"Length correlations: {stats['summary']['n_length_correlations']:,}\n\n")
        
        # Top single AA associations
        f.write("-"*70 + "\n")
        f.write("TOP 50 SINGLE AA ASSOCIATIONS (by confidence)\n")
        f.write("-"*70 + "\n\n")
        
        sorted_aa = sorted(stats['single_aa_associations'], 
                          key=lambda x: (x['confidence'], x['n']), reverse=True)
        
        for assoc in sorted_aa[:50]:
            f.write(f"{assoc['cdr']}[{assoc['cdr_pos']}]={assoc['cdr_aa']} → "
                   f"{assoc['fw']}[{assoc['fw_pos']}]={assoc['predicted_fw_aa']}: "
                   f"{assoc['confidence']:.1f}% (n={assoc['n']:,})\n")
        
        # Top triplet associations
        f.write("\n" + "-"*70 + "\n")
        f.write("TOP 50 TRIPLET ASSOCIATIONS (by confidence)\n")
        f.write("-"*70 + "\n\n")
        
        sorted_trip = sorted(stats['triplet_associations'],
                            key=lambda x: (x['confidence'], x['n']), reverse=True)
        
        for assoc in sorted_trip[:50]:
            f.write(f"{assoc['junction_from']} {assoc['cdr_triplet']} → "
                   f"{assoc['junction_to']} {assoc['predicted_fw_motif']}: "
                   f"{assoc['confidence']:.1f}% (n={assoc['n']:,})\n")
        
        # Length correlations by CDR-FW pair
        f.write("\n" + "-"*70 + "\n")
        f.write("LENGTH CORRELATIONS\n")
        f.write("-"*70 + "\n\n")
        
        for cdr in ['cdr1', 'cdr2', 'cdr3']:
            for fw in ['FR1', 'FR2', 'FR3', 'FR4']:
                relevant = [x for x in stats['length_correlations'] 
                           if x['cdr'] == cdr and x['fw'] == fw]
                if relevant:
                    avg_conf = np.mean([x['confidence'] for x in relevant])
                    f.write(f"{cdr.upper()} vs {fw}: avg confidence {avg_conf:.1f}%\n")
    
    print(f"  Summary written to: {output_path}")

# ============================================================
# MAIN
# ============================================================

def main():
    parser = argparse.ArgumentParser(description='Analyze CDR-Framework correlations')
    parser.add_argument('npz_files', nargs='+', help='NPZ files to process (supports wildcards)')
    parser.add_argument('--output', '-o', default='.', help='Output directory')
    parser.add_argument('--min-n', type=int, default=50, help='Minimum sample size for associations')
    parser.add_argument('--visualize', action='store_true', help='Generate visualizations')
    
    args = parser.parse_args()
    
    # Expand wildcards
    npz_files = []
    for pattern in args.npz_files:
        if '*' in pattern:
            npz_files.extend(glob(pattern))
        else:
            npz_files.append(pattern)
    
    npz_files = [f for f in npz_files if f.endswith('.npz') and os.path.exists(f)]
    
    if not npz_files:
        print("ERROR: No valid NPZ files found")
        sys.exit(1)
    
    print("="*70)
    print("CDR-FRAMEWORK CORRELATION ANALYZER")
    print("="*70)
    print(f"\nFound {len(npz_files)} NPZ files to process\n")
    
    # Process all shards
    all_results = []
    for npz_file in npz_files:
        result = process_shard(npz_file)
        if result:
            all_results.append(result)
    
    if not all_results:
        print("ERROR: No valid results")
        sys.exit(1)
    
    # Merge
    print(f"\nMerging {len(all_results)} shard results...")
    merged = merge_results(all_results)
    print(f"  Total: {merged['n_valid']:,} valid / {merged['n_total']:,} total sequences")
    
    # Compute statistics
    print("\nComputing statistics...")
    stats = compute_statistics(merged, min_n=args.min_n)
    print(f"  Single AA associations: {stats['summary']['n_single_aa_associations']:,}")
    print(f"  Triplet associations: {stats['summary']['n_triplet_associations']:,}")
    print(f"  Length correlations: {stats['summary']['n_length_correlations']:,}")
    
    # Save outputs
    os.makedirs(args.output, exist_ok=True)
    
    # Raw results (for passing to Claude or further analysis)
    pkl_path = os.path.join(args.output, 'correlation_results.pkl')
    with open(pkl_path, 'wb') as f:
        pickle.dump({
            'raw_counts': merged,
            'statistics': stats
        }, f)
    print(f"\n✓ Raw results saved to: {pkl_path}")
    
    # Human-readable summary
    summary_path = os.path.join(args.output, 'correlation_summary.txt')
    write_summary(stats, summary_path)
    
    # Visualizations
    if args.visualize:
        try:
            from visualize_correlations import create_visualizations
            create_visualizations(stats, args.output)
        except ImportError:
            print("\n  Note: Run visualize_correlations.py separately for visualizations")
            print("  Or pass correlation_results.pkl to Claude for visualization")
    
    print(f"\n{'='*70}")
    print("COMPLETE!")
    print(f"{'='*70}")
    print(f"\nTo visualize results:")
    print(f"  1. Run: python visualize_correlations.py {pkl_path}")
    print(f"  2. Or upload {pkl_path} to Claude for visualization")

if __name__ == '__main__':
    main()
