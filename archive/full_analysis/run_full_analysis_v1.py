#!/usr/bin/env python3
"""
VHH CDR-Framework Complete Analysis Pipeline
=============================================
Master script that runs all analyses sequentially.

Usage:
    python run_full_analysis.py /path/to/shards/*.npz --output ./results
    
    # With custom framework to score:
    python run_full_analysis.py /path/to/shards/*.npz --output ./results \
        --fr1 "EVQLVESGGGLVQPGGSLRLSCAAS" \
        --fr2 "WFRQAPGKGREFVA" \
        --fr3 "YYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYC" \
        --fr4 "WGQGTQVTVSS"
    
    # With specific CDRs to score against:
    python run_full_analysis.py /path/to/shards/*.npz --output ./results \
        --cdr1 "GFTFSSYA" --cdr2 "ISYDGSNK" --cdr3 "ARDLLVRY"

Output structure:
    ./results/
    ├── 1_correlations/
    │   ├── correlation_results.pkl
    │   ├── correlation_summary.txt
    │   └── figures/
    ├── 2_advanced_correlations/
    │   ├── correlation_results_advanced.pkl
    │   ├── correlation_summary_advanced.txt
    │   └── germline_vs_universal.txt
    ├── 3_pfr_cdr_models/
    │   ├── pfr_cdr_models.pkl
    │   ├── model_summary.txt
    │   └── figures/
    ├── 4_framework_scoring/
    │   └── framework_report.txt
    └── FINAL_REPORT.md
"""

import numpy as np
from collections import defaultdict, Counter
import pickle
import sys
import os
import time
import argparse
from glob import glob
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

# ============================================================
# CONFIGURATION
# ============================================================

DEFAULT_FRAMEWORK = {
    'FR1': 'EVQLVESGGGLVQPGGSLRLSCAAS',
    'FR2': 'WFRQAPGKGREFVA',
    'FR3': 'YYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYC',
    'FR4': 'WGQGTQVTVSS'
}

# Amino acid properties
HYDROPHOBICITY = {
    'I': 4.5, 'V': 4.2, 'L': 3.8, 'F': 2.8, 'C': 2.5,
    'M': 1.9, 'A': 1.8, 'G': -0.4, 'T': -0.7, 'S': -0.8,
    'W': -0.9, 'Y': -1.3, 'P': -1.6, 'H': -3.2, 'E': -3.5,
    'Q': -3.5, 'D': -3.5, 'N': -3.5, 'K': -3.9, 'R': -4.5
}

CHARGE = {
    'K': 1, 'R': 1, 'H': 0.1, 'D': -1, 'E': -1,
    'A': 0, 'C': 0, 'F': 0, 'G': 0, 'I': 0, 'L': 0, 'M': 0,
    'N': 0, 'P': 0, 'Q': 0, 'S': 0, 'T': 0, 'V': 0, 'W': 0, 'Y': 0
}

AROMATIC = {'F': 1, 'W': 1, 'Y': 1, 'H': 0.5}
POLAR = {'S': 1, 'T': 1, 'N': 1, 'Q': 1, 'Y': 0.5, 'H': 0.5, 'C': 0.5}

AA_LIST = list('ACDEFGHIKLMNPQRSTVWY')
AA_TO_IDX = {aa: i for i, aa in enumerate(AA_LIST)}

VERNIER_POSITIONS = {
    'FR1': [-3, -2, -1],
    'FR2': [0, 1, 2, -3, -2, -1],
    'FR3': [0, 1, 2, 3, 4, -5, -4, -3, -2, -1],
    'FR4': [0, 1, 2, 3]
}

HALLMARK_POSITIONS = {
    'FR1': [0, 5, 10, 15, 20],
    'FR2': [0, 5, 10, 14],
    'FR3': [0, 5, 10, 15, 20, 25, 30, 35],
    'FR4': [0, 5]
}

# ============================================================
# UTILITY FUNCTIONS
# ============================================================

def print_header(title):
    print("\n" + "="*70)
    print(f"  {title}")
    print("="*70 + "\n")

def print_step(step_num, title):
    print(f"\n{'─'*70}")
    print(f"  STEP {step_num}: {title}")
    print(f"{'─'*70}\n")

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
    
    if len(fr1) < 10 or len(fr2) < 10 or len(fr3) < 20 or len(fr4) < 5:
        return None
    if len(fr1) > 40 or len(fr2) > 25 or len(fr3) > 50 or len(fr4) > 20:
        return None
    
    return {
        'FR1': fr1, 'FR2': fr2, 'FR3': fr3, 'FR4': fr4,
        'cdr1': cdr1, 'cdr2': cdr2, 'cdr3': cdr3
    }

def get_framework_signature(fw_dict):
    sig = []
    for fw_name in ['FR1', 'FR2', 'FR3', 'FR4']:
        fw_seq = fw_dict[fw_name]
        for pos in HALLMARK_POSITIONS[fw_name]:
            if pos < len(fw_seq):
                sig.append(fw_seq[pos])
            else:
                sig.append('X')
    return ''.join(sig)

def compute_sequence_properties(seq):
    if not seq or len(seq) == 0:
        return {'length': 0, 'charge': 0, 'hydrophobicity': 0, 
                'aromaticity': 0, 'polar_fraction': 0}
    seq = str(seq).upper()
    n = len(seq)
    return {
        'length': n,
        'charge': sum(CHARGE.get(aa, 0) for aa in seq),
        'hydrophobicity': sum(HYDROPHOBICITY.get(aa, 0) for aa in seq) / n,
        'aromaticity': sum(AROMATIC.get(aa, 0) for aa in seq) / n,
        'polar_fraction': sum(POLAR.get(aa, 0) for aa in seq) / n
    }

def extract_cdr_features(cdr1, cdr2, cdr3):
    features = []
    feature_names = []
    
    features.extend([len(cdr1), len(cdr2), len(cdr3)])
    feature_names.extend(['cdr1_len', 'cdr2_len', 'cdr3_len'])
    
    for cdr_name, cdr_seq in [('cdr1', cdr1), ('cdr2', cdr2), ('cdr3', cdr3)]:
        props = compute_sequence_properties(cdr_seq)
        features.extend([props['charge'], props['hydrophobicity'], 
                        props['aromaticity'], props['polar_fraction']])
        feature_names.extend([f'{cdr_name}_charge', f'{cdr_name}_hydro',
                             f'{cdr_name}_aromatic', f'{cdr_name}_polar'])
    
    for cdr_name, cdr_seq in [('cdr1', cdr1), ('cdr2', cdr2), ('cdr3', cdr3)]:
        for i in range(2):
            features.append(AA_TO_IDX.get(cdr_seq[i], 10) if i < len(cdr_seq) else 10)
            feature_names.append(f'{cdr_name}_pos{i}')
        for i in [-2, -1]:
            features.append(AA_TO_IDX.get(cdr_seq[i], 10) if abs(i) <= len(cdr_seq) else 10)
            feature_names.append(f'{cdr_name}_pos{i}')
    
    return np.array(features, dtype=np.float32), feature_names

# ============================================================
# STEP 1: BASIC CORRELATION ANALYSIS
# ============================================================

def run_basic_correlations(npz_files, output_dir):
    """Run basic CDR-FW correlation analysis"""
    
    print_step(1, "BASIC CORRELATION ANALYSIS")
    os.makedirs(output_dir, exist_ok=True)
    
    # Initialize counts
    single_aa_counts = defaultdict(int)
    triplet_counts = defaultdict(int)
    length_counts = defaultdict(int)
    
    n_valid_total = 0
    n_total_total = 0
    
    cdr_positions = {'cdr1': [0, 1, 2, -3, -2, -1], 'cdr2': [0, 1, 2, -3, -2, -1], 'cdr3': [0, 1, 2, -3, -2, -1]}
    fw_positions = {'FR1': list(range(-10, 0)), 'FR2': list(range(10)) + list(range(-10, 0)),
                   'FR3': list(range(10)) + list(range(-10, 0)), 'FR4': list(range(10))}
    
    for npz_path in npz_files:
        print(f"  Processing: {os.path.basename(npz_path)}", end=" ", flush=True)
        start = time.time()
        
        try:
            data = np.load(npz_path, allow_pickle=True)
            full_seqs, cdr1s, cdr2s, cdr3s = data['aa_v_full'], data['cdr1'], data['cdr2'], data['cdr3']
        except Exception as e:
            print(f"ERROR: {e}")
            continue
        
        n_valid = 0
        for i in range(len(full_seqs)):
            fw = extract_frameworks(full_seqs[i], cdr1s[i], cdr2s[i], cdr3s[i])
            if fw is None:
                continue
            n_valid += 1
            
            # Single AA
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
                                single_aa_counts[(cdr_name, cdr_pos, cdr_aa, fw_name, fw_pos, fw_aa)] += 1
                            except IndexError:
                                continue
            
            # Triplets
            junctions = [
                ('CDR1_start', 'FR1_end', fw['cdr1'][:3], fw['FR1'][-5:]),
                ('CDR1_end', 'FR2_start', fw['cdr1'][-3:], fw['FR2'][:5]),
                ('CDR2_start', 'FR2_end', fw['cdr2'][:3], fw['FR2'][-5:]),
                ('CDR2_end', 'FR3_start', fw['cdr2'][-3:], fw['FR3'][:5]),
                ('CDR3_start', 'FR3_end', fw['cdr3'][:3], fw['FR3'][-5:]),
                ('CDR3_end', 'FR4_start', fw['cdr3'][-3:], fw['FR4'][:5]),
            ]
            for jf, jt, ct, fm in junctions:
                if len(ct) == 3 and len(fm) == 5:
                    triplet_counts[(jf, jt, ct, fm)] += 1
            
            # Lengths
            for cdr_name in ['cdr1', 'cdr2', 'cdr3']:
                for fw_name in ['FR1', 'FR2', 'FR3', 'FR4']:
                    length_counts[(cdr_name, len(fw[cdr_name]), fw_name, len(fw[fw_name]))] += 1
        
        n_valid_total += n_valid
        n_total_total += len(full_seqs)
        print(f"→ {n_valid:,} valid ({time.time()-start:.1f}s)")
    
    # Compute statistics
    print("\n  Computing statistics...")
    stats = {'single_aa': [], 'triplets': [], 'lengths': [], 'summary': {}}
    
    # Single AA
    aa_pairs = defaultdict(lambda: defaultdict(int))
    for key, count in single_aa_counts.items():
        pair_key = key[:5]
        aa_pairs[pair_key][key[5]] += count
    
    for pair_key, dist in aa_pairs.items():
        total = sum(dist.values())
        if total < 50:
            continue
        top_aa = max(dist, key=dist.get)
        stats['single_aa'].append({
            'cdr': pair_key[0], 'cdr_pos': pair_key[1], 'cdr_aa': pair_key[2],
            'fw': pair_key[3], 'fw_pos': pair_key[4], 'predicted_fw_aa': top_aa,
            'confidence': 100 * dist[top_aa] / total, 'n': total
        })
    
    # Triplets
    trip_pairs = defaultdict(lambda: defaultdict(int))
    for key, count in triplet_counts.items():
        pair_key = key[:3]
        trip_pairs[pair_key][key[3]] += count
    
    for pair_key, dist in trip_pairs.items():
        total = sum(dist.values())
        if total < 50:
            continue
        top_motif = max(dist, key=dist.get)
        stats['triplets'].append({
            'junction_from': pair_key[0], 'junction_to': pair_key[1],
            'cdr_triplet': pair_key[2], 'predicted_fw_motif': top_motif,
            'confidence': 100 * dist[top_motif] / total, 'n': total
        })
    
    stats['summary'] = {
        'n_sequences': n_valid_total,
        'n_single_aa': len(stats['single_aa']),
        'n_triplets': len(stats['triplets'])
    }
    
    # Save
    with open(os.path.join(output_dir, 'correlation_results.pkl'), 'wb') as f:
        pickle.dump({'statistics': stats, 'raw_counts': {
            'single_aa': dict(single_aa_counts),
            'triplets': dict(triplet_counts),
            'lengths': dict(length_counts)
        }}, f)
    
    # Write summary
    with open(os.path.join(output_dir, 'correlation_summary.txt'), 'w') as f:
        f.write(f"Basic Correlation Analysis\n{'='*50}\n\n")
        f.write(f"Sequences analyzed: {n_valid_total:,}\n")
        f.write(f"Single AA associations: {len(stats['single_aa']):,}\n")
        f.write(f"Triplet associations: {len(stats['triplets']):,}\n")
    
    print(f"\n  ✓ Results saved to {output_dir}")
    print(f"    Sequences: {n_valid_total:,}")
    print(f"    Single AA rules: {len(stats['single_aa']):,}")
    print(f"    Triplet rules: {len(stats['triplets']):,}")
    
    return stats

# ============================================================
# STEP 2: ADVANCED CORRELATION ANALYSIS (GERMLINE CONTROL)
# ============================================================

def run_advanced_correlations(npz_files, output_dir):
    """Run advanced analysis with germline control"""
    
    print_step(2, "ADVANCED CORRELATION ANALYSIS (GERMLINE CONTROL)")
    os.makedirs(output_dir, exist_ok=True)
    
    global_single_aa = defaultdict(int)
    cluster_single_aa = defaultdict(lambda: defaultdict(int))
    fw_signatures = Counter()
    mi_data = defaultdict(lambda: defaultdict(int))
    
    n_valid_total = 0
    
    cdr_positions = {'cdr1': [0, 1, 2, -3, -2, -1], 'cdr2': [0, 1, 2, -3, -2, -1], 'cdr3': [0, 1, 2, -3, -2, -1]}
    fw_positions = {'FR1': list(range(-10, 0)), 'FR2': list(range(10)) + list(range(-10, 0)),
                   'FR3': list(range(10)) + list(range(-10, 0)), 'FR4': list(range(10))}
    
    for npz_path in npz_files:
        print(f"  Processing: {os.path.basename(npz_path)}", end=" ", flush=True)
        start = time.time()
        
        try:
            data = np.load(npz_path, allow_pickle=True)
            full_seqs, cdr1s, cdr2s, cdr3s = data['aa_v_full'], data['cdr1'], data['cdr2'], data['cdr3']
        except:
            print("ERROR")
            continue
        
        n_valid = 0
        for i in range(len(full_seqs)):
            fw = extract_frameworks(full_seqs[i], cdr1s[i], cdr2s[i], cdr3s[i])
            if fw is None:
                continue
            n_valid += 1
            
            fw_sig = get_framework_signature(fw)
            fw_signatures[fw_sig] += 1
            
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
                                key = (cdr_name, cdr_pos, cdr_aa, fw_name, fw_pos, fw_aa)
                                global_single_aa[key] += 1
                                cluster_single_aa[fw_sig][key] += 1
                                mi_data[(cdr_name, cdr_pos, fw_name, fw_pos)][(cdr_aa, fw_aa)] += 1
                            except IndexError:
                                continue
        
        n_valid_total += n_valid
        n_clusters = len(set(fw_signatures.keys()))
        print(f"→ {n_valid:,} valid, {n_clusters} clusters ({time.time()-start:.1f}s)")
    
    # Classify universal vs germline-specific
    print("\n  Classifying universal vs germline-specific rules...")
    
    top_clusters = [c[0] for c in sorted(fw_signatures.items(), key=lambda x: x[1], reverse=True)[:50] if c[1] >= 100]
    
    universal_rules = []
    germline_rules = []
    
    aa_pairs = defaultdict(lambda: defaultdict(int))
    for key, count in global_single_aa.items():
        pair_key = key[:5]
        aa_pairs[pair_key][key[5]] += count
    
    for pair_key, dist in aa_pairs.items():
        total = sum(dist.values())
        if total < 50:
            continue
        
        top_aa = max(dist, key=dist.get)
        global_conf = 100 * dist[top_aa] / total
        
        # Check cluster consistency
        cluster_agrees = 0
        cluster_tested = 0
        
        for cluster_id in top_clusters:
            cluster_dist = defaultdict(int)
            for key, count in cluster_single_aa[cluster_id].items():
                if key[:5] == pair_key:
                    cluster_dist[key[5]] += count
            
            if sum(cluster_dist.values()) < 20:
                continue
            
            cluster_tested += 1
            if max(cluster_dist, key=cluster_dist.get) == top_aa:
                cluster_agrees += 1
        
        if cluster_tested < 3:
            continue
        
        consistency = cluster_agrees / cluster_tested
        is_vernier = pair_key[4] in VERNIER_POSITIONS.get(pair_key[3], [])
        
        rule = {
            'cdr': pair_key[0], 'cdr_pos': pair_key[1], 'cdr_aa': pair_key[2],
            'fw': pair_key[3], 'fw_pos': pair_key[4], 'predicted_fw_aa': top_aa,
            'confidence': global_conf, 'n': total, 'consistency': consistency,
            'is_vernier': is_vernier
        }
        
        if consistency >= 0.8:
            universal_rules.append(rule)
        else:
            germline_rules.append(rule)
    
    # Compute MI
    print("  Computing mutual information...")
    mi_results = []
    for pos_pair, aa_counts in mi_data.items():
        total = sum(aa_counts.values())
        if total < 1000:
            continue
        
        cdr_marg = defaultdict(int)
        fw_marg = defaultdict(int)
        for (ca, fa), c in aa_counts.items():
            cdr_marg[ca] += c
            fw_marg[fa] += c
        
        mi = 0.0
        for (ca, fa), c in aa_counts.items():
            p_joint = c / total
            p_cdr = cdr_marg[ca] / total
            p_fw = fw_marg[fa] / total
            if p_joint > 0 and p_cdr > 0 and p_fw > 0:
                mi += p_joint * np.log2(p_joint / (p_cdr * p_fw))
        
        h_cdr = -sum((c/total) * np.log2(c/total) for c in cdr_marg.values() if c > 0)
        h_fw = -sum((c/total) * np.log2(c/total) for c in fw_marg.values() if c > 0)
        nmi = mi / max(h_cdr, h_fw) if max(h_cdr, h_fw) > 0 else 0
        
        mi_results.append({
            'cdr': pos_pair[0], 'cdr_pos': pos_pair[1],
            'fw': pos_pair[2], 'fw_pos': pos_pair[3],
            'nmi': nmi, 'n': total,
            'is_vernier': pos_pair[3] in VERNIER_POSITIONS.get(pos_pair[2], [])
        })
    
    mi_results = sorted(mi_results, key=lambda x: x['nmi'], reverse=True)
    
    # Save
    stats = {
        'universal_rules': sorted(universal_rules, key=lambda x: x['confidence'], reverse=True),
        'germline_rules': sorted(germline_rules, key=lambda x: x['confidence'], reverse=True),
        'mutual_information': mi_results,
        'summary': {
            'n_sequences': n_valid_total,
            'n_clusters': len(fw_signatures),
            'n_universal': len(universal_rules),
            'n_germline': len(germline_rules)
        }
    }
    
    with open(os.path.join(output_dir, 'correlation_results_advanced.pkl'), 'wb') as f:
        pickle.dump(stats, f)
    
    # Write summaries
    with open(os.path.join(output_dir, 'correlation_summary_advanced.txt'), 'w') as f:
        f.write(f"Advanced Correlation Analysis\n{'='*50}\n\n")
        f.write(f"Sequences: {n_valid_total:,}\n")
        f.write(f"FW clusters: {len(fw_signatures):,}\n")
        f.write(f"Universal rules: {len(universal_rules):,}\n")
        f.write(f"Germline-specific rules: {len(germline_rules):,}\n\n")
        f.write("TOP 20 UNIVERSAL RULES\n" + "-"*50 + "\n")
        for r in universal_rules[:20]:
            v = "[VERNIER]" if r['is_vernier'] else ""
            f.write(f"{r['cdr']}[{r['cdr_pos']}]={r['cdr_aa']} → {r['fw']}[{r['fw_pos']}]={r['predicted_fw_aa']}: "
                   f"{r['confidence']:.1f}% {v}\n")
    
    with open(os.path.join(output_dir, 'germline_vs_universal.txt'), 'w') as f:
        f.write("UNIVERSAL vs GERMLINE-SPECIFIC RULES\n")
        f.write("="*50 + "\n\n")
        f.write(f"UNIVERSAL (consistent across germlines): {len(universal_rules)}\n")
        f.write(f"GERMLINE-SPECIFIC (vary by cluster): {len(germline_rules)}\n")
    
    print(f"\n  ✓ Results saved to {output_dir}")
    print(f"    Universal rules: {len(universal_rules):,}")
    print(f"    Germline-specific rules: {len(germline_rules):,}")
    
    return stats

# ============================================================
# STEP 3: BUILD p(FR|CDR) MODELS
# ============================================================

def build_pfr_cdr_models(npz_files, output_dir, max_samples=500000):
    """Build conditional probability models"""
    
    print_step(3, "BUILD p(FR|CDR) MODELS")
    os.makedirs(output_dir, exist_ok=True)
    
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler
    
    target_positions = {
        'FR1': list(range(-10, 0)),
        'FR2': list(range(15)),
        'FR3': list(range(10)) + list(range(-10, 0)),
        'FR4': list(range(11))
    }
    
    position_data = defaultdict(lambda: {'X': [], 'y': []})
    n_collected = 0
    feature_names = None
    
    print(f"  Collecting training data (max {max_samples:,} samples)...")
    
    for npz_path in npz_files:
        if n_collected >= max_samples:
            break
        print(f"  Loading: {os.path.basename(npz_path)}", end=" ", flush=True)
        
        try:
            data = np.load(npz_path, allow_pickle=True)
            full_seqs, cdr1s, cdr2s, cdr3s = data['aa_v_full'], data['cdr1'], data['cdr2'], data['cdr3']
        except:
            print("ERROR")
            continue
        
        n_valid = 0
        for i in range(len(full_seqs)):
            if n_collected >= max_samples:
                break
            
            fw = extract_frameworks(full_seqs[i], cdr1s[i], cdr2s[i], cdr3s[i])
            if fw is None:
                continue
            
            features, feature_names = extract_cdr_features(fw['cdr1'], fw['cdr2'], fw['cdr3'])
            
            for fw_name, positions in target_positions.items():
                fw_seq = fw[fw_name]
                for pos in positions:
                    try:
                        aa = fw_seq[pos]
                        if aa not in AA_TO_IDX:
                            continue
                        position_data[(fw_name, pos)]['X'].append(features)
                        position_data[(fw_name, pos)]['y'].append(AA_TO_IDX[aa])
                    except IndexError:
                        continue
            
            n_valid += 1
            n_collected += 1
        
        print(f"→ {n_valid:,}")
    
    print(f"\n  Training {len(position_data)} position models...")
    
    models = {}
    scalers = {}
    metrics = {}
    
    for idx, ((fw_name, fw_pos), data) in enumerate(position_data.items()):
        X = np.array(data['X'])
        y = np.array(data['y'])
        
        if len(X) < 100:
            continue
        
        unique = np.unique(y)
        if len(unique) < 2:
            models[(fw_name, fw_pos)] = {'type': 'constant', 'aa': AA_LIST[unique[0]]}
            metrics[(fw_name, fw_pos)] = {'accuracy': 1.0, 'n_samples': len(X), 'top_class': AA_LIST[unique[0]], 'top_class_freq': 1.0}
            continue
        
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X)
        scalers[(fw_name, fw_pos)] = scaler
        
        model = LogisticRegression(max_iter=500, multi_class='multinomial', solver='lbfgs', class_weight='balanced', n_jobs=-1)
        
        try:
            model.fit(X_scaled, y)
            acc = model.score(X_scaled, y)
            counts = Counter(y)
            top_class = max(counts, key=counts.get)
            top_freq = counts[top_class] / len(y)
            
            models[(fw_name, fw_pos)] = {'type': 'logistic', 'model': model, 'classes': model.classes_}
            metrics[(fw_name, fw_pos)] = {
                'accuracy': acc, 'n_samples': len(X), 'n_classes': len(unique),
                'top_class': AA_LIST[top_class], 'top_class_freq': top_freq,
                'improvement': acc - top_freq
            }
        except:
            continue
        
        if (idx + 1) % 20 == 0:
            print(f"    Trained {idx+1}/{len(position_data)}...")
    
    # Compute importances
    importances = {}
    for (fw_name, fw_pos), m in models.items():
        if m['type'] == 'logistic':
            coefs = np.abs(m['model'].coef_).mean(axis=0)
            importances[(fw_name, fw_pos)] = [(fn, coefs[i]) for i, fn in enumerate(feature_names)]
    
    # Summary
    summary = {
        'total_positions': len(models),
        'fw_accuracy': {fw: np.mean([m['accuracy'] for (f, p), m in metrics.items() if f == fw]) 
                       for fw in ['FR1', 'FR2', 'FR3', 'FR4']},
        'most_cdr_dependent': sorted([(k, v.get('improvement', 0)) for k, v in metrics.items()], 
                                     key=lambda x: x[1], reverse=True)[:20]
    }
    
    # Save
    with open(os.path.join(output_dir, 'pfr_cdr_models.pkl'), 'wb') as f:
        pickle.dump({
            'models': models, 'scalers': scalers, 'metrics': metrics,
            'feature_names': feature_names, 'importances': importances, 'summary': summary
        }, f)
    
    with open(os.path.join(output_dir, 'model_summary.txt'), 'w') as f:
        f.write(f"p(FR|CDR) Model Summary\n{'='*50}\n\n")
        f.write(f"Training samples: {n_collected:,}\n")
        f.write(f"Positions modeled: {len(models)}\n\n")
        f.write("Avg accuracy by framework:\n")
        for fw, acc in summary['fw_accuracy'].items():
            f.write(f"  {fw}: {acc:.1%}\n")
    
    print(f"\n  ✓ Models saved to {output_dir}")
    print(f"    Positions modeled: {len(models)}")
    
    return {'models': models, 'scalers': scalers, 'metrics': metrics, 'feature_names': feature_names}

# ============================================================
# STEP 4: SCORE FRAMEWORK
# ============================================================

def score_framework(model_data, framework, cdr1, cdr2, cdr3, output_dir):
    """Score a framework against CDRs"""
    
    print_step(4, "FRAMEWORK SCORING")
    os.makedirs(output_dir, exist_ok=True)
    
    models = model_data['models']
    scalers = model_data['scalers']
    
    print(f"  Framework: {framework['FR1'][:15]}.../{framework['FR2']}/.../{framework['FR4']}")
    print(f"  CDRs: {cdr1} / {cdr2} / {cdr3}")
    
    features = extract_cdr_features(cdr1, cdr2, cdr3)[0].reshape(1, -1)
    
    results = []
    
    for (fw_name, fw_pos), model_info in models.items():
        try:
            fw_aa = framework[fw_name][fw_pos]
        except:
            continue
        
        if model_info['type'] == 'constant':
            expected = model_info['aa']
            prob = 1.0 if fw_aa == expected else 0.0
        else:
            model = model_info['model']
            scaler = scalers.get((fw_name, fw_pos))
            if scaler is None:
                continue
            
            X_scaled = scaler.transform(features)
            probs = model.predict_proba(X_scaled)[0]
            classes = model_info['classes']
            
            prob_dist = {AA_LIST[c]: probs[i] for i, c in enumerate(classes)}
            prob = prob_dist.get(fw_aa, 0.0)
            expected = AA_LIST[classes[np.argmax(probs)]]
        
        results.append({
            'fw': fw_name, 'pos': fw_pos, 'current': fw_aa,
            'expected': expected, 'prob': prob, 'match': fw_aa == expected
        })
    
    # Sort by probability
    results = sorted(results, key=lambda x: x['prob'])
    
    # Write report
    report_path = os.path.join(output_dir, 'framework_report.txt')
    with open(report_path, 'w') as f:
        f.write("FRAMEWORK SCORING REPORT\n")
        f.write("="*60 + "\n\n")
        f.write(f"CDR1: {cdr1}\n")
        f.write(f"CDR2: {cdr2}\n")
        f.write(f"CDR3: {cdr3}\n\n")
        
        n_match = sum(1 for r in results if r['match'])
        f.write(f"Matching: {n_match}/{len(results)} ({100*n_match/len(results):.1f}%)\n\n")
        
        f.write("PROBLEM POSITIONS (prob < 10%)\n")
        f.write("-"*60 + "\n")
        f.write(f"{'Position':<12} {'Current':<8} {'Expected':<8} {'Prob':<8}\n")
        f.write("-"*60 + "\n")
        
        problems = [r for r in results if r['prob'] < 0.1]
        for r in problems[:20]:
            f.write(f"{r['fw']}[{r['pos']}]".ljust(12))
            f.write(f"{r['current']:<8} {r['expected']:<8} {r['prob']:.1%}\n")
        
        f.write("\n\nRECOMMENDED MUTATIONS\n")
        f.write("-"*60 + "\n")
        for i, r in enumerate(problems[:10], 1):
            f.write(f"  {i}. {r['fw']}[{r['pos']}]: {r['current']} → {r['expected']}\n")
    
    print(f"\n  ✓ Report saved to {report_path}")
    print(f"    Problem positions: {len(problems)}")
    
    # Print summary
    print("\n  TOP PROBLEM POSITIONS:")
    for r in problems[:5]:
        print(f"    {r['fw']}[{r['pos']}]: {r['current']} → {r['expected']} (prob: {r['prob']:.1%})")
    
    return results

# ============================================================
# STEP 5: GENERATE VISUALIZATIONS
# ============================================================

def generate_visualizations(corr_path, model_path, output_dir):
    """Generate all visualizations"""
    
    print_step(5, "GENERATING VISUALIZATIONS")
    os.makedirs(output_dir, exist_ok=True)
    
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import seaborn as sns
        plt.style.use('seaborn-v0_8-whitegrid')
    except Exception as e:
        print(f"  ⚠️  Matplotlib not available: {e}")
        return
    
    # Load model data
    if os.path.exists(model_path):
        with open(model_path, 'rb') as f:
            model_data = pickle.load(f)
        
        metrics = model_data['metrics']
        
        # Conservation plot
        positions = sorted(metrics.items(), key=lambda x: x[1].get('top_class_freq', 0), reverse=True)
        
        fig, ax = plt.subplots(figsize=(16, 6))
        labels = [f"{p[0][0]}[{p[0][1]}]" for p in positions]
        values = [p[1].get('top_class_freq', 0) for p in positions]
        colors = ['#27ae60' if v > 0.95 else '#f1c40f' if v > 0.8 else '#e74c3c' for v in values]
        
        ax.bar(range(len(labels)), values, color=colors)
        ax.axhline(0.95, color='green', linestyle='--', alpha=0.5)
        ax.axhline(0.80, color='orange', linestyle='--', alpha=0.5)
        ax.set_xticks(range(0, len(labels), 3))
        ax.set_xticklabels([labels[i] for i in range(0, len(labels), 3)], rotation=90, fontsize=7)
        ax.set_ylabel('Conservation')
        ax.set_title('Framework Position Conservation')
        
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'position_conservation.png'), dpi=150)
        plt.close()
        print("    Saved: position_conservation.png")
        
        # Feature importance
        if model_data.get('importances'):
            feature_names = model_data['feature_names']
            avg_imp = defaultdict(float)
            count = 0
            for pos, imps in model_data['importances'].items():
                for fn, val in imps:
                    avg_imp[fn] += val
                count += 1
            avg_imp = {k: v/count for k, v in avg_imp.items()}
            
            fig, ax = plt.subplots(figsize=(10, 8))
            sorted_imp = sorted(avg_imp.items(), key=lambda x: x[1], reverse=True)[:15]
            ax.barh([x[0] for x in sorted_imp], [x[1] for x in sorted_imp])
            ax.set_xlabel('Average Importance')
            ax.set_title('Top Predictive CDR Features')
            ax.invert_yaxis()
            
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, 'feature_importance.png'), dpi=150)
            plt.close()
            print("    Saved: feature_importance.png")
    
    print(f"\n  ✓ Visualizations saved to {output_dir}")

# ============================================================
# FINAL REPORT
# ============================================================

def write_final_report(output_dir, results):
    """Write comprehensive final report"""
    
    report_path = os.path.join(output_dir, 'FINAL_REPORT.md')
    
    with open(report_path, 'w') as f:
        f.write("# VHH CDR-Framework Analysis - Final Report\n\n")
        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        f.write("## Summary\n\n")
        
        if 'basic' in results:
            s = results['basic'].get('summary', {})
            f.write(f"- **Sequences analyzed**: {s.get('n_sequences', 'N/A'):,}\n")
            f.write(f"- **Single AA associations**: {s.get('n_single_aa', 'N/A'):,}\n")
            f.write(f"- **Triplet associations**: {s.get('n_triplets', 'N/A'):,}\n\n")
        
        if 'advanced' in results:
            s = results['advanced'].get('summary', {})
            f.write(f"- **FW clusters (germline proxies)**: {s.get('n_clusters', 'N/A'):,}\n")
            f.write(f"- **Universal rules**: {s.get('n_universal', 'N/A'):,}\n")
            f.write(f"- **Germline-specific rules**: {s.get('n_germline', 'N/A'):,}\n\n")
        
        f.write("## Output Files\n\n")
        f.write("```\n")
        f.write(f"{output_dir}/\n")
        f.write("├── 1_correlations/\n")
        f.write("│   ├── correlation_results.pkl\n")
        f.write("│   └── correlation_summary.txt\n")
        f.write("├── 2_advanced_correlations/\n")
        f.write("│   ├── correlation_results_advanced.pkl\n")
        f.write("│   └── germline_vs_universal.txt\n")
        f.write("├── 3_pfr_cdr_models/\n")
        f.write("│   ├── pfr_cdr_models.pkl\n")
        f.write("│   └── model_summary.txt\n")
        f.write("├── 4_framework_scoring/\n")
        f.write("│   └── framework_report.txt\n")
        f.write("├── 5_figures/\n")
        f.write("│   ├── position_conservation.png\n")
        f.write("│   └── feature_importance.png\n")
        f.write("└── FINAL_REPORT.md\n")
        f.write("```\n\n")
        
        f.write("## Next Steps\n\n")
        f.write("1. Review `4_framework_scoring/framework_report.txt` for problem positions\n")
        f.write("2. Upload `3_pfr_cdr_models/pfr_cdr_models.pkl` to Claude for custom analysis\n")
        f.write("3. Design focused library at problem positions\n")
        f.write("4. Validate experimentally\n")
    
    print(f"\n  ✓ Final report: {report_path}")

# ============================================================
# MAIN
# ============================================================

def main():
    parser = argparse.ArgumentParser(description='Complete VHH CDR-Framework Analysis Pipeline')
    parser.add_argument('npz_files', nargs='+', help='NPZ files to process')
    parser.add_argument('--output', '-o', default='./cdr_fw_analysis', help='Output directory')
    parser.add_argument('--max-samples', type=int, default=500000, help='Max samples for model training')
    parser.add_argument('--fr1', default=DEFAULT_FRAMEWORK['FR1'], help='FR1 sequence')
    parser.add_argument('--fr2', default=DEFAULT_FRAMEWORK['FR2'], help='FR2 sequence')
    parser.add_argument('--fr3', default=DEFAULT_FRAMEWORK['FR3'], help='FR3 sequence')
    parser.add_argument('--fr4', default=DEFAULT_FRAMEWORK['FR4'], help='FR4 sequence')
    parser.add_argument('--cdr1', default='GFTFSSYA', help='CDR1 for scoring')
    parser.add_argument('--cdr2', default='ISYDGSNK', help='CDR2 for scoring')
    parser.add_argument('--cdr3', default='ARDLYYYYGMDV', help='CDR3 for scoring')
    parser.add_argument('--skip-basic', action='store_true', help='Skip basic correlation analysis')
    parser.add_argument('--skip-advanced', action='store_true', help='Skip advanced correlation analysis')
    
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
    
    framework = {'FR1': args.fr1, 'FR2': args.fr2, 'FR3': args.fr3, 'FR4': args.fr4}
    
    print_header("VHH CDR-FRAMEWORK COMPLETE ANALYSIS PIPELINE")
    print(f"  NPZ files: {len(npz_files)}")
    print(f"  Output: {args.output}")
    print(f"  Max samples: {args.max_samples:,}")
    
    start_time = time.time()
    os.makedirs(args.output, exist_ok=True)
    
    results = {}
    
    # Step 1: Basic correlations
    if not args.skip_basic:
        results['basic'] = run_basic_correlations(npz_files, os.path.join(args.output, '1_correlations'))
    
    # Step 2: Advanced correlations
    if not args.skip_advanced:
        results['advanced'] = run_advanced_correlations(npz_files, os.path.join(args.output, '2_advanced_correlations'))
    
    # Step 3: Build models
    model_data = build_pfr_cdr_models(npz_files, os.path.join(args.output, '3_pfr_cdr_models'), args.max_samples)
    
    # Step 4: Score framework
    score_framework(model_data, framework, args.cdr1, args.cdr2, args.cdr3, 
                   os.path.join(args.output, '4_framework_scoring'))
    
    # Step 5: Visualizations
    generate_visualizations(
        os.path.join(args.output, '1_correlations', 'correlation_results.pkl'),
        os.path.join(args.output, '3_pfr_cdr_models', 'pfr_cdr_models.pkl'),
        os.path.join(args.output, '5_figures')
    )
    
    # Final report
    write_final_report(args.output, results)
    
    elapsed = time.time() - start_time
    
    print_header("ANALYSIS COMPLETE")
    print(f"  Total time: {elapsed/60:.1f} minutes")
    print(f"  Results: {args.output}")
    print(f"\n  Upload pkl files to Claude for further analysis!")

if __name__ == '__main__':
    main()
