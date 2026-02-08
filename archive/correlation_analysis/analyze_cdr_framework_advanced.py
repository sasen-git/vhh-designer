#!/usr/bin/env python3
"""
Advanced VHH CDR-Framework Correlation Analyzer
================================================
Accounts for germline preferences by:
1. Clustering sequences by framework similarity (proxy for V-gene family)
2. Computing associations WITHIN each cluster
3. Identifying UNIVERSAL associations (consistent across clusters)
4. Flagging GERMLINE-SPECIFIC associations (only in certain clusters)

Additional analyses:
- Vernier zone / contact position flagging
- Mutual information for position pairs
- Conditional probability p(FR | CDR) estimates

Usage:
    python analyze_cdr_framework_advanced.py /path/to/shards/*.npz --output ./results

Output:
    - correlation_results_advanced.pkl: Full analysis with germline control
    - correlation_summary_advanced.txt: Human-readable summary
    - germline_vs_universal.txt: Which rules are germline-specific vs universal
"""

import numpy as np
from collections import defaultdict, Counter
import pickle
import sys
import os
import time
import argparse
from glob import glob

# ============================================================
# VERNIER ZONE POSITIONS (FR residues that contact CDRs)
# Based on canonical VHH structures - IMGT numbering approximated
# ============================================================

VERNIER_POSITIONS = {
    'FR1': [-3, -2, -1],  # Last 3 positions contact CDR1
    'FR2': [0, 1, 2, -3, -2, -1],  # Both ends contact CDR1/CDR2
    'FR3': [0, 1, 2, 3, 4, -5, -4, -3, -2, -1],  # Extended contact zone
    'FR4': [0, 1, 2, 3]  # First positions contact CDR3
}

# Key "hallmark" positions that differ between VHH germlines
# These help cluster by germline family
HALLMARK_POSITIONS = {
    'FR1': [0, 5, 10, 15, 20],
    'FR2': [0, 5, 10, 14],
    'FR3': [0, 5, 10, 15, 20, 25, 30, 35],
    'FR4': [0, 5]
}

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
    
    if len(fr1) < 10 or len(fr2) < 10 or len(fr3) < 20 or len(fr4) < 5:
        return None
    if len(fr1) > 40 or len(fr2) > 25 or len(fr3) > 50 or len(fr4) > 20:
        return None
    
    return {
        'FR1': fr1, 'FR2': fr2, 'FR3': fr3, 'FR4': fr4,
        'cdr1': cdr1, 'cdr2': cdr2, 'cdr3': cdr3
    }

def get_framework_signature(fw_dict):
    """Get a simplified signature for germline clustering"""
    sig = []
    for fw_name in ['FR1', 'FR2', 'FR3', 'FR4']:
        fw_seq = fw_dict[fw_name]
        for pos in HALLMARK_POSITIONS[fw_name]:
            if pos < len(fw_seq):
                sig.append(fw_seq[pos])
            else:
                sig.append('X')
    return ''.join(sig)

# ============================================================
# SHARD PROCESSING WITH GERMLINE CLUSTERING
# ============================================================

def process_shard_advanced(npz_path, verbose=True):
    """Process shard with framework clustering for germline control"""
    
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
    
    # Initialize counts - now keyed by (cluster_id, association)
    # Global counts (across all clusters)
    global_single_aa = defaultdict(int)
    global_triplets = defaultdict(int)
    global_lengths = defaultdict(int)
    
    # Per-cluster counts
    cluster_single_aa = defaultdict(lambda: defaultdict(int))
    cluster_triplets = defaultdict(lambda: defaultdict(int))
    
    # Framework signature counts (for clustering)
    fw_signatures = Counter()
    
    # Mutual information data (position pairs)
    mi_data = defaultdict(lambda: defaultdict(int))
    
    n_valid = 0
    
    for i in range(n_total):
        fw = extract_frameworks(full_seqs[i], cdr1s[i], cdr2s[i], cdr3s[i])
        if fw is None:
            continue
        
        n_valid += 1
        
        # Get framework signature for germline clustering
        fw_sig = get_framework_signature(fw)
        fw_signatures[fw_sig] += 1
        
        # === SINGLE AA CORRELATIONS ===
        cdr_positions = {
            'cdr1': [0, 1, 2, -3, -2, -1],
            'cdr2': [0, 1, 2, -3, -2, -1],
            'cdr3': [0, 1, 2, -3, -2, -1]
        }
        
        fw_positions = {
            'FR1': list(range(-10, 0)),
            'FR2': list(range(10)) + list(range(-10, 0)),
            'FR3': list(range(10)) + list(range(-10, 0)),
            'FR4': list(range(10))
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
                        global_single_aa[key] += 1
                        cluster_single_aa[fw_sig][key] += 1
                        
                        # For MI calculation
                        mi_key = (cdr_name, cdr_pos, fw_name, fw_pos)
                        mi_data[mi_key][(cdr_aa, fw_aa)] += 1
        
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
                key = (junction_from, junction_to, cdr_trip, fw_motif)
                global_triplets[key] += 1
                cluster_triplets[fw_sig][key] += 1
        
        # === LENGTH CORRELATIONS ===
        for cdr_name in ['cdr1', 'cdr2', 'cdr3']:
            cdr_len = len(fw[cdr_name])
            for fw_name in ['FR1', 'FR2', 'FR3', 'FR4']:
                fw_len = len(fw[fw_name])
                global_lengths[(cdr_name, cdr_len, fw_name, fw_len)] += 1
    
    elapsed = time.time() - start_time
    if verbose:
        n_clusters = len(fw_signatures)
        print(f"→ {n_valid:,}/{n_total:,} valid, {n_clusters} FW clusters ({elapsed:.1f}s)")
    
    return {
        'global_single_aa': dict(global_single_aa),
        'global_triplets': dict(global_triplets),
        'global_lengths': dict(global_lengths),
        'cluster_single_aa': {k: dict(v) for k, v in cluster_single_aa.items()},
        'cluster_triplets': {k: dict(v) for k, v in cluster_triplets.items()},
        'fw_signatures': dict(fw_signatures),
        'mi_data': {k: dict(v) for k, v in mi_data.items()},
        'n_valid': n_valid,
        'n_total': n_total
    }

# ============================================================
# MERGE RESULTS
# ============================================================

def merge_results_advanced(results_list):
    """Merge results from multiple shards"""
    merged = {
        'global_single_aa': defaultdict(int),
        'global_triplets': defaultdict(int),
        'global_lengths': defaultdict(int),
        'cluster_single_aa': defaultdict(lambda: defaultdict(int)),
        'cluster_triplets': defaultdict(lambda: defaultdict(int)),
        'fw_signatures': Counter(),
        'mi_data': defaultdict(lambda: defaultdict(int)),
        'n_valid': 0,
        'n_total': 0
    }
    
    for results in results_list:
        if results is None:
            continue
        
        merged['n_valid'] += results['n_valid']
        merged['n_total'] += results['n_total']
        
        for key, count in results['global_single_aa'].items():
            merged['global_single_aa'][key] += count
        
        for key, count in results['global_triplets'].items():
            merged['global_triplets'][key] += count
        
        for key, count in results['global_lengths'].items():
            merged['global_lengths'][key] += count
        
        for sig, counts in results['cluster_single_aa'].items():
            for key, count in counts.items():
                merged['cluster_single_aa'][sig][key] += count
        
        for sig, counts in results['cluster_triplets'].items():
            for key, count in counts.items():
                merged['cluster_triplets'][sig][key] += count
        
        merged['fw_signatures'].update(results['fw_signatures'])
        
        for pos_pair, aa_counts in results['mi_data'].items():
            for aa_pair, count in aa_counts.items():
                merged['mi_data'][pos_pair][aa_pair] += count
    
    return {
        'global_single_aa': dict(merged['global_single_aa']),
        'global_triplets': dict(merged['global_triplets']),
        'global_lengths': dict(merged['global_lengths']),
        'cluster_single_aa': {k: dict(v) for k, v in merged['cluster_single_aa'].items()},
        'cluster_triplets': {k: dict(v) for k, v in merged['cluster_triplets'].items()},
        'fw_signatures': dict(merged['fw_signatures']),
        'mi_data': {k: dict(v) for k, v in merged['mi_data'].items()},
        'n_valid': merged['n_valid'],
        'n_total': merged['n_total']
    }

# ============================================================
# COMPUTE MUTUAL INFORMATION
# ============================================================

def compute_mutual_information(mi_data, min_n=1000):
    """Compute mutual information between CDR and FW positions"""
    mi_results = []
    
    for pos_pair, aa_counts in mi_data.items():
        cdr_name, cdr_pos, fw_name, fw_pos = pos_pair
        
        total = sum(aa_counts.values())
        if total < min_n:
            continue
        
        # Marginals
        cdr_marginal = defaultdict(int)
        fw_marginal = defaultdict(int)
        for (cdr_aa, fw_aa), count in aa_counts.items():
            cdr_marginal[cdr_aa] += count
            fw_marginal[fw_aa] += count
        
        # Compute MI
        mi = 0.0
        for (cdr_aa, fw_aa), count in aa_counts.items():
            p_joint = count / total
            p_cdr = cdr_marginal[cdr_aa] / total
            p_fw = fw_marginal[fw_aa] / total
            
            if p_joint > 0 and p_cdr > 0 and p_fw > 0:
                mi += p_joint * np.log2(p_joint / (p_cdr * p_fw))
        
        # Normalized MI (0-1 scale)
        h_cdr = -sum((c/total) * np.log2(c/total) for c in cdr_marginal.values() if c > 0)
        h_fw = -sum((c/total) * np.log2(c/total) for c in fw_marginal.values() if c > 0)
        
        nmi = mi / max(h_cdr, h_fw) if max(h_cdr, h_fw) > 0 else 0
        
        # Check if position is in Vernier zone
        is_vernier = fw_pos in VERNIER_POSITIONS.get(fw_name, [])
        
        mi_results.append({
            'cdr': cdr_name,
            'cdr_pos': cdr_pos,
            'fw': fw_name,
            'fw_pos': fw_pos,
            'mi': mi,
            'nmi': nmi,
            'n': total,
            'is_vernier': is_vernier
        })
    
    return sorted(mi_results, key=lambda x: x['nmi'], reverse=True)

# ============================================================
# IDENTIFY UNIVERSAL VS GERMLINE-SPECIFIC ASSOCIATIONS
# ============================================================

def classify_associations(merged, min_n=50, min_clusters=3):
    """
    Classify associations as:
    - UNIVERSAL: consistent across multiple germline clusters
    - GERMLINE-SPECIFIC: only present in certain clusters
    """
    
    # Get top clusters (by sequence count)
    top_clusters = sorted(merged['fw_signatures'].items(), 
                         key=lambda x: x[1], reverse=True)[:50]
    top_cluster_ids = set(c[0] for c in top_clusters if c[1] >= 100)
    
    print(f"  Analyzing {len(top_cluster_ids)} major FW clusters...")
    
    # For each association, check consistency across clusters
    results = {
        'universal': [],
        'germline_specific': [],
        'cluster_info': {}
    }
    
    # Analyze single AA associations
    global_aa = merged['global_single_aa']
    cluster_aa = merged['cluster_single_aa']
    
    # Group by (cdr, cdr_pos, cdr_aa, fw, fw_pos) -> fw_aa distribution
    aa_associations = defaultdict(lambda: defaultdict(int))
    for key, count in global_aa.items():
        cdr, cdr_pos, cdr_aa, fw, fw_pos, fw_aa = key
        assoc_key = (cdr, cdr_pos, cdr_aa, fw, fw_pos)
        aa_associations[assoc_key][fw_aa] += count
    
    for assoc_key, fw_aa_dist in aa_associations.items():
        total = sum(fw_aa_dist.values())
        if total < min_n:
            continue
        
        top_fw_aa = max(fw_aa_dist, key=fw_aa_dist.get)
        global_conf = 100 * fw_aa_dist[top_fw_aa] / total
        
        # Check consistency across clusters
        cluster_agreements = 0
        cluster_total = 0
        cluster_details = []
        
        for cluster_id in top_cluster_ids:
            cluster_counts = cluster_aa.get(cluster_id, {})
            
            # Get distribution for this association in this cluster
            cluster_dist = defaultdict(int)
            for key, count in cluster_counts.items():
                if key[:5] == assoc_key:
                    cluster_dist[key[5]] += count
            
            cluster_n = sum(cluster_dist.values())
            if cluster_n < 20:
                continue
            
            cluster_total += 1
            cluster_top = max(cluster_dist, key=cluster_dist.get) if cluster_dist else None
            
            if cluster_top == top_fw_aa:
                cluster_agreements += 1
                cluster_conf = 100 * cluster_dist[cluster_top] / cluster_n
                cluster_details.append({
                    'cluster': cluster_id[:20],
                    'agrees': True,
                    'conf': cluster_conf,
                    'n': cluster_n
                })
            else:
                cluster_details.append({
                    'cluster': cluster_id[:20],
                    'agrees': False,
                    'top_aa': cluster_top,
                    'n': cluster_n
                })
        
        if cluster_total < min_clusters:
            continue
        
        consistency = cluster_agreements / cluster_total if cluster_total > 0 else 0
        
        # Check if in Vernier zone
        fw_name, fw_pos = assoc_key[3], assoc_key[4]
        is_vernier = fw_pos in VERNIER_POSITIONS.get(fw_name, [])
        
        assoc_info = {
            'cdr': assoc_key[0],
            'cdr_pos': assoc_key[1],
            'cdr_aa': assoc_key[2],
            'fw': fw_name,
            'fw_pos': fw_pos,
            'predicted_fw_aa': top_fw_aa,
            'global_confidence': global_conf,
            'n': total,
            'consistency': consistency,
            'n_clusters_tested': cluster_total,
            'n_clusters_agree': cluster_agreements,
            'is_vernier': is_vernier,
            'cluster_details': cluster_details
        }
        
        if consistency >= 0.8:  # 80%+ clusters agree
            results['universal'].append(assoc_info)
        else:
            results['germline_specific'].append(assoc_info)
    
    # Sort by confidence
    results['universal'] = sorted(results['universal'], 
                                  key=lambda x: x['global_confidence'], reverse=True)
    results['germline_specific'] = sorted(results['germline_specific'],
                                          key=lambda x: x['global_confidence'], reverse=True)
    
    # Store cluster info
    results['cluster_info'] = {
        'n_clusters': len(top_cluster_ids),
        'top_clusters': [(c[0][:30], c[1]) for c in top_clusters[:20]]
    }
    
    return results

# ============================================================
# COMPUTE FULL STATISTICS
# ============================================================

def compute_statistics_advanced(merged, min_n=50):
    """Compute comprehensive statistics"""
    
    stats = {
        'single_aa_associations': [],
        'triplet_associations': [],
        'length_correlations': [],
        'mutual_information': [],
        'universal_vs_germline': {},
        'summary': {}
    }
    
    # Standard single AA stats (same as before)
    single_aa = merged['global_single_aa']
    
    cdr_fw_pairs = defaultdict(lambda: defaultdict(int))
    for key, count in single_aa.items():
        cdr, cdr_pos, cdr_aa, fw, fw_pos, fw_aa = key
        pair_key = (cdr, cdr_pos, cdr_aa, fw, fw_pos)
        cdr_fw_pairs[pair_key][fw_aa] += count
    
    for pair_key, fw_aa_counts in cdr_fw_pairs.items():
        cdr, cdr_pos, cdr_aa, fw, fw_pos = pair_key
        total = sum(fw_aa_counts.values())
        
        if total < min_n:
            continue
        
        top_fw_aa = max(fw_aa_counts, key=fw_aa_counts.get)
        confidence = 100 * fw_aa_counts[top_fw_aa] / total
        
        is_vernier = fw_pos in VERNIER_POSITIONS.get(fw, [])
        
        stats['single_aa_associations'].append({
            'cdr': cdr,
            'cdr_pos': cdr_pos,
            'cdr_aa': cdr_aa,
            'fw': fw,
            'fw_pos': fw_pos,
            'predicted_fw_aa': top_fw_aa,
            'confidence': confidence,
            'n': total,
            'is_vernier': is_vernier,
            'distribution': dict(fw_aa_counts)
        })
    
    # Triplet stats
    triplets = merged['global_triplets']
    
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
    
    # Mutual information
    print("  Computing mutual information...")
    stats['mutual_information'] = compute_mutual_information(merged['mi_data'])
    
    # Universal vs germline-specific
    print("  Classifying universal vs germline-specific associations...")
    stats['universal_vs_germline'] = classify_associations(merged)
    
    # Summary
    n_universal = len(stats['universal_vs_germline']['universal'])
    n_germline = len(stats['universal_vs_germline']['germline_specific'])
    
    stats['summary'] = {
        'n_sequences': merged['n_valid'],
        'n_fw_clusters': len(merged['fw_signatures']),
        'n_single_aa_associations': len(stats['single_aa_associations']),
        'n_triplet_associations': len(stats['triplet_associations']),
        'n_universal_rules': n_universal,
        'n_germline_specific_rules': n_germline,
        'n_vernier_associations': len([a for a in stats['single_aa_associations'] if a['is_vernier']])
    }
    
    return stats

# ============================================================
# WRITE SUMMARY
# ============================================================

def write_summary_advanced(stats, output_dir):
    """Write comprehensive summary"""
    
    # Main summary
    summary_path = os.path.join(output_dir, 'correlation_summary_advanced.txt')
    with open(summary_path, 'w') as f:
        f.write("="*70 + "\n")
        f.write("ADVANCED CDR-FRAMEWORK CORRELATION ANALYSIS\n")
        f.write("="*70 + "\n\n")
        
        s = stats['summary']
        f.write(f"Total sequences: {s['n_sequences']:,}\n")
        f.write(f"Framework clusters (germline proxies): {s['n_fw_clusters']:,}\n")
        f.write(f"Single AA associations: {s['n_single_aa_associations']:,}\n")
        f.write(f"  - In Vernier zone: {s['n_vernier_associations']:,}\n")
        f.write(f"Triplet associations: {s['n_triplet_associations']:,}\n")
        f.write(f"\nUNIVERSAL rules (consistent across germlines): {s['n_universal_rules']:,}\n")
        f.write(f"GERMLINE-SPECIFIC rules: {s['n_germline_specific_rules']:,}\n")
        
        # Top universal rules
        f.write("\n" + "="*70 + "\n")
        f.write("TOP 30 UNIVERSAL RULES (True CDR-FW Compatibility)\n")
        f.write("These are consistent across germline families\n")
        f.write("="*70 + "\n\n")
        
        for assoc in stats['universal_vs_germline']['universal'][:30]:
            vernier = "[VERNIER]" if assoc['is_vernier'] else ""
            f.write(f"{assoc['cdr']}[{assoc['cdr_pos']}]={assoc['cdr_aa']} → "
                   f"{assoc['fw']}[{assoc['fw_pos']}]={assoc['predicted_fw_aa']}: "
                   f"{assoc['global_confidence']:.1f}% "
                   f"(n={assoc['n']:,}, {assoc['n_clusters_agree']}/{assoc['n_clusters_tested']} clusters) "
                   f"{vernier}\n")
        
        # Top germline-specific rules
        f.write("\n" + "="*70 + "\n")
        f.write("TOP 30 GERMLINE-SPECIFIC RULES (May be confounded)\n")
        f.write("These vary by germline family - use with caution\n")
        f.write("="*70 + "\n\n")
        
        for assoc in stats['universal_vs_germline']['germline_specific'][:30]:
            vernier = "[VERNIER]" if assoc['is_vernier'] else ""
            f.write(f"{assoc['cdr']}[{assoc['cdr_pos']}]={assoc['cdr_aa']} → "
                   f"{assoc['fw']}[{assoc['fw_pos']}]={assoc['predicted_fw_aa']}: "
                   f"{assoc['global_confidence']:.1f}% "
                   f"(n={assoc['n']:,}, {assoc['n_clusters_agree']}/{assoc['n_clusters_tested']} clusters) "
                   f"{vernier}\n")
        
        # Top MI pairs
        f.write("\n" + "="*70 + "\n")
        f.write("TOP 30 POSITION PAIRS BY MUTUAL INFORMATION\n")
        f.write("="*70 + "\n\n")
        
        for mi in stats['mutual_information'][:30]:
            vernier = "[VERNIER]" if mi['is_vernier'] else ""
            f.write(f"{mi['cdr']}[{mi['cdr_pos']}] ↔ {mi['fw']}[{mi['fw_pos']}]: "
                   f"NMI={mi['nmi']:.3f} (n={mi['n']:,}) {vernier}\n")
    
    print(f"  Summary: {summary_path}")
    
    # Germline vs Universal detailed report
    report_path = os.path.join(output_dir, 'germline_vs_universal.txt')
    with open(report_path, 'w') as f:
        f.write("GERMLINE VS UNIVERSAL ASSOCIATION REPORT\n")
        f.write("="*70 + "\n\n")
        
        f.write("INTERPRETATION GUIDE:\n")
        f.write("-"*70 + "\n")
        f.write("UNIVERSAL rules: Same CDR→FW pattern in 80%+ of germline clusters\n")
        f.write("  → These reflect TRUE structural/functional requirements\n")
        f.write("  → SAFE to use for framework design\n\n")
        f.write("GERMLINE-SPECIFIC rules: Pattern varies by germline cluster\n")
        f.write("  → May be evolutionary accident, not functional requirement\n")
        f.write("  → Use with CAUTION - may not generalize\n\n")
        f.write("VERNIER positions: Framework positions that contact CDRs\n")
        f.write("  → Rules at these positions are more likely to be real\n")
        f.write("-"*70 + "\n\n")
        
        f.write(f"UNIVERSAL: {len(stats['universal_vs_germline']['universal']):,} rules\n")
        f.write(f"GERMLINE-SPECIFIC: {len(stats['universal_vs_germline']['germline_specific']):,} rules\n")
    
    print(f"  Report: {report_path}")

# ============================================================
# MAIN
# ============================================================

def main():
    parser = argparse.ArgumentParser(description='Advanced CDR-Framework correlation analysis')
    parser.add_argument('npz_files', nargs='+', help='NPZ files to process')
    parser.add_argument('--output', '-o', default='.', help='Output directory')
    parser.add_argument('--min-n', type=int, default=50, help='Minimum sample size')
    
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
    print("ADVANCED CDR-FRAMEWORK CORRELATION ANALYZER")
    print("(With Germline Control & Mutual Information)")
    print("="*70)
    print(f"\nFound {len(npz_files)} NPZ files\n")
    
    # Process shards
    all_results = []
    for npz_file in npz_files:
        result = process_shard_advanced(npz_file)
        if result:
            all_results.append(result)
    
    if not all_results:
        print("ERROR: No valid results")
        sys.exit(1)
    
    # Merge
    print(f"\nMerging {len(all_results)} shards...")
    merged = merge_results_advanced(all_results)
    print(f"  Total: {merged['n_valid']:,} sequences")
    print(f"  Framework clusters: {len(merged['fw_signatures']):,}")
    
    # Compute statistics
    print("\nComputing statistics...")
    stats = compute_statistics_advanced(merged, min_n=args.min_n)
    
    # Save
    os.makedirs(args.output, exist_ok=True)
    
    pkl_path = os.path.join(args.output, 'correlation_results_advanced.pkl')
    with open(pkl_path, 'wb') as f:
        pickle.dump({
            'raw_counts': merged,
            'statistics': stats
        }, f)
    print(f"\n✓ Results: {pkl_path}")
    
    write_summary_advanced(stats, args.output)
    
    print(f"\n{'='*70}")
    print("SUMMARY")
    print(f"{'='*70}")
    print(f"Sequences: {stats['summary']['n_sequences']:,}")
    print(f"FW clusters (germline proxies): {stats['summary']['n_fw_clusters']:,}")
    print(f"\nUNIVERSAL rules (true compatibility): {stats['summary']['n_universal_rules']:,}")
    print(f"GERMLINE-SPECIFIC rules (use caution): {stats['summary']['n_germline_specific_rules']:,}")
    print(f"\nVernier zone associations: {stats['summary']['n_vernier_associations']:,}")
    print(f"\nUpload {pkl_path} to Claude for visualization!")

if __name__ == '__main__':
    main()
