#!/usr/bin/env python3
"""
LOW MEMORY VERSION v3 - Germline-Aware CDR-Framework Analysis

Features:
- Germline classification using FR2 hallmark positions (42/49/50/52)
- Cysteine counting for C2 vs C4 classification
- Species/germline family grouping (Alpaca Y_C2/F_C2/F_C4, Llama VHH1-5, Camel, VH-like)
- TRUE CDR-dependence analysis WITHIN germline families
- Separates germline effects from CDR-framework compatibility

Usage:
    python run_full_analysis_lowmem_v3.py /path/to/shards/*.npz --output ./results_v3
"""

import argparse
import numpy as np
import pickle
import gc
import sys
import time
from pathlib import Path
from collections import defaultdict
from datetime import datetime

# ============================================================
# CONFIGURATION
# ============================================================

AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY")
AA_TO_IDX = {aa: i for i, aa in enumerate(AMINO_ACIDS)}

# FR2 Hallmark positions (0-indexed in our FR2)
# IMGT 42 → FR2[3], IMGT 49 → FR2[10], IMGT 50 → FR2[11], IMGT 52 → FR2[13]
HALLMARK_POSITIONS = {
    'pos42': 3,   # Y/F - key VHH hallmark
    'pos49': 10,  # E/Q
    'pos50': 11,  # R/C
    'pos52': 13,  # L/G/F/W - distinguishes VHH from VH
}

# Germline family definitions based on FR2 hallmarks
GERMLINE_FAMILIES = {
    # Alpaca families
    'Alpaca_Y_C2': {
        'description': 'IGHV3S53-like (Y-E-R-G/L)',
        'fr2_42': ['Y'],
        'fr2_52': ['G', 'L', 'F'],  # NOT W
        'cys_pattern': 'C2',
    },
    'Alpaca_F_C2': {
        'description': 'IGHV3-3-like (F-E-R-L)',
        'fr2_42': ['F'],
        'fr2_52': ['L', 'G', 'F'],
        'cys_pattern': 'C2',
    },
    'Alpaca_F_C4': {
        'description': 'IGHV3S65-like (F-E-R-L, extra Cys)',
        'fr2_42': ['F'],
        'fr2_52': ['L', 'G', 'F'],
        'cys_pattern': 'C4',
    },
    # Llama families
    'Llama_VHH1': {
        'description': 'IGHV3S1 (F-E-R-G)',
        'fr2_42': ['F'],
        'fr2_52': ['G'],
        'cys_pattern': 'C2',
    },
    'Llama_VHH3': {
        'description': 'IGHV3S3 (Y-E-R-W)',
        'fr2_42': ['Y'],
        'fr2_52': ['W'],
        'cys_pattern': 'C2',
    },
    'Llama_VHH4': {
        'description': 'IGHV3S4 (Y-E-R-L, extra Cys)',
        'fr2_42': ['Y'],
        'fr2_52': ['L'],
        'cys_pattern': 'C4',
    },
    # VH-like (non-classical)
    'VH_like': {
        'description': 'VH-derived (V-G-L-W)',
        'fr2_42': ['V'],
        'fr2_52': ['W'],
        'cys_pattern': 'any',
    },
}

# AA similarity groups
AA_GROUPS = {
    'hydrophobic': set('AILMVFWP'),
    'aromatic': set('FWY'),
    'polar': set('STNQ'),
    'positive': set('KRH'),
    'negative': set('DE'),
    'small': set('GAS'),
    'aliphatic': set('AILV'),
}

CDR_POSITIONS = {
    'cdr1': [0, 1, 2, -3, -2, -1],
    'cdr2': [0, 1, 2, -3, -2, -1],
    'cdr3': [0, 1, 2, -3, -2, -1],
}

FR_POSITIONS = {
    'FR1': list(range(25)),
    'FR2': list(range(14)),
    'FR3': list(range(32)),
    'FR4': list(range(11)),
}

ALL_CDR_FW_PAIRS = [
    (cdr, fw) 
    for cdr in ['cdr1', 'cdr2', 'cdr3'] 
    for fw in ['FR1', 'FR2', 'FR3', 'FR4']
]

# ============================================================
# PROGRESS BAR
# ============================================================

class ProgressBar:
    def __init__(self, total, width=40, prefix=''):
        self.total = total
        self.width = width
        self.prefix = prefix
        self.current = 0
        self.start_time = time.time()
        
    def update(self, n=1):
        self.current += n
        self._display()
    
    def _display(self):
        if self.total == 0:
            return
        pct = self.current / self.total
        filled = int(self.width * pct)
        bar = '█' * filled + '░' * (self.width - filled)
        elapsed = time.time() - self.start_time
        if self.current > 0:
            eta = (elapsed / self.current) * (self.total - self.current)
            eta_str = f"ETA {eta:.0f}s"
        else:
            eta_str = "ETA --"
        sys.stdout.write(f'\r    {self.prefix}|{bar}| {self.current:,}/{self.total:,} ({pct*100:.1f}%) {eta_str}   ')
        sys.stdout.flush()
    
    def finish(self):
        elapsed = time.time() - self.start_time
        sys.stdout.write(f'\r    {self.prefix}|{"█" * self.width}| {self.total:,}/{self.total:,} (100%) Done in {elapsed:.1f}s\n')
        sys.stdout.flush()

# ============================================================
# GERMLINE CLASSIFICATION
# ============================================================

def count_cysteines(seq_data):
    """Count total cysteines across all regions."""
    total_cys = 0
    for region in ['FR1', 'FR2', 'FR3', 'FR4', 'cdr1', 'cdr2', 'cdr3']:
        seq = seq_data.get(region, '')
        if seq:
            total_cys += seq.count('C')
    return total_cys

def classify_germline(seq_data):
    """
    Classify sequence into germline family based on FR2 hallmarks and Cys count.
    
    Returns: (family_name, hallmark_dict)
    """
    fr2 = seq_data.get('FR2', '')
    if len(fr2) < 14:
        return ('Unknown', {})
    
    # Extract hallmark residues
    hallmarks = {
        'pos42': fr2[HALLMARK_POSITIONS['pos42']] if len(fr2) > HALLMARK_POSITIONS['pos42'] else '?',
        'pos49': fr2[HALLMARK_POSITIONS['pos49']] if len(fr2) > HALLMARK_POSITIONS['pos49'] else '?',
        'pos50': fr2[HALLMARK_POSITIONS['pos50']] if len(fr2) > HALLMARK_POSITIONS['pos50'] else '?',
        'pos52': fr2[HALLMARK_POSITIONS['pos52']] if len(fr2) > HALLMARK_POSITIONS['pos52'] else '?',
    }
    
    # Count cysteines
    n_cys = count_cysteines(seq_data)
    cys_pattern = 'C4' if n_cys >= 4 else 'C2'
    hallmarks['n_cys'] = n_cys
    hallmarks['cys_pattern'] = cys_pattern
    
    # Motif-based classification
    pos42 = hallmarks['pos42']
    pos52 = hallmarks['pos52']
    
    # Check for VH-like first (V at 42, W at 52)
    if pos42 == 'V' and pos52 == 'W':
        return ('VH_like', hallmarks)
    
    # Check for VHH with W at 52 (like Llama VHH3)
    if pos42 in ['Y', 'F'] and pos52 == 'W':
        return ('VHH_W52', hallmarks)
    
    # Classical VHH patterns
    if pos42 == 'Y':
        if cys_pattern == 'C4':
            return ('Y_C4', hallmarks)
        else:
            return ('Y_C2', hallmarks)
    elif pos42 == 'F':
        if cys_pattern == 'C4':
            return ('F_C4', hallmarks)
        else:
            return ('F_C2', hallmarks)
    
    # Other/Unknown
    return ('Other_VHH', hallmarks)

def get_motif_string(hallmarks):
    """Create a motif string like 'F-E-R-L' from hallmarks."""
    return f"{hallmarks.get('pos42', '?')}-{hallmarks.get('pos49', '?')}-{hallmarks.get('pos50', '?')}-{hallmarks.get('pos52', '?')}"

# ============================================================
# HELPER FUNCTIONS
# ============================================================

def get_position_value(seq, pos):
    if seq is None or len(seq) == 0:
        return None
    try:
        if pos < 0:
            if len(seq) >= abs(pos):
                return seq[pos]
        else:
            if len(seq) > pos:
                return seq[pos]
    except:
        pass
    return None

def extract_frameworks(full_seq, cdr1, cdr2, cdr3):
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
    
    return {'FR1': fr1, 'FR2': fr2, 'FR3': fr3, 'FR4': fr4,
            'cdr1': cdr1, 'cdr2': cdr2, 'cdr3': cdr3}

def process_shard_streaming(npz_path, accumulator, show_progress=True):
    try:
        data = np.load(npz_path, allow_pickle=True)
        
        if 'aa_v_full' in data and 'cdr1' in data and 'cdr2' in data and 'cdr3' in data:
            full_seqs = data['aa_v_full']
            cdr1s = data['cdr1']
            cdr2s = data['cdr2']
            cdr3s = data['cdr3']
            
            n_total = len(full_seqs)
            n_valid = 0
            
            if show_progress:
                pbar = ProgressBar(n_total)
            
            for i in range(n_total):
                fw = extract_frameworks(full_seqs[i], cdr1s[i], cdr2s[i], cdr3s[i])
                if fw:
                    accumulator.add_sequence(fw)
                    n_valid += 1
                
                if show_progress and i % 10000 == 0:
                    pbar.update(10000 if i > 0 else 0)
            
            if show_progress:
                pbar.current = n_total
                pbar.finish()
            
            data.close()
            return n_valid, n_total
        
        data.close()
        return 0, 0
        
    except Exception as e:
        print(f"\n  Warning: Error processing {npz_path}: {e}")
        return 0, 0

# ============================================================
# GERMLINE-AWARE ACCUMULATOR
# ============================================================

class GermlineAwareAccumulator:
    """Accumulates statistics separately for each germline family."""
    
    def __init__(self):
        # Global stats
        self.n_sequences = 0
        self.germline_counts = defaultdict(int)
        self.motif_counts = defaultdict(int)
        
        # Per-germline accumulators
        # Structure: family -> cdr -> cdr_pos -> fw -> fw_pos -> cdr_aa -> fw_aa -> count
        self.family_single_counts = defaultdict(
            lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(
                lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(int)))))))
        
        # Per-germline triplet counts
        # Structure: family -> cdr -> fw -> cdr_triplet -> fw_triplet -> count
        self.family_triplet_counts = defaultdict(
            lambda: defaultdict(lambda: defaultdict(
                lambda: defaultdict(lambda: defaultdict(int)))))
        
        # Global (all families) for comparison
        self.global_single_counts = defaultdict(lambda: defaultdict(lambda: defaultdict(
            lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(int))))))
        
        # FW marginals per family
        self.family_fw_marginals = defaultdict(
            lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(int))))
        
        # Hallmark distribution tracking
        self.hallmark_distributions = {
            'pos42': defaultdict(int),
            'pos49': defaultdict(int),
            'pos50': defaultdict(int),
            'pos52': defaultdict(int),
        }
        
        # CDR length distributions per family
        self.family_cdr_lengths = defaultdict(lambda: defaultdict(list))
    
    def add_sequence(self, seq_data):
        self.n_sequences += 1
        
        # Classify germline
        family, hallmarks = classify_germline(seq_data)
        self.germline_counts[family] += 1
        
        # Track motif
        motif = get_motif_string(hallmarks)
        self.motif_counts[motif] += 1
        
        # Track hallmark distributions
        for pos in ['pos42', 'pos49', 'pos50', 'pos52']:
            aa = hallmarks.get(pos, '?')
            self.hallmark_distributions[pos][aa] += 1
        
        # Track CDR lengths per family
        for cdr in ['cdr1', 'cdr2', 'cdr3']:
            seq = seq_data.get(cdr, '')
            if seq:
                self.family_cdr_lengths[family][cdr].append(len(seq))
        
        # FW marginals per family
        for fw in ['FR1', 'FR2', 'FR3', 'FR4']:
            fw_seq = seq_data.get(fw, '')
            if not fw_seq:
                continue
            for fw_pos in FR_POSITIONS[fw]:
                fw_aa = get_position_value(fw_seq, fw_pos)
                if fw_aa in AA_TO_IDX:
                    self.family_fw_marginals[family][fw][fw_pos][fw_aa] += 1
        
        # Process CDR-FW correlations
        for cdr, fw in ALL_CDR_FW_PAIRS:
            cdr_seq = seq_data.get(cdr, '')
            fw_seq = seq_data.get(fw, '')
            
            if not cdr_seq or not fw_seq:
                continue
            
            # Single position correlations
            for cdr_pos in CDR_POSITIONS[cdr]:
                cdr_aa = get_position_value(cdr_seq, cdr_pos)
                if cdr_aa not in AA_TO_IDX:
                    continue
                
                for fw_pos in FR_POSITIONS[fw]:
                    fw_aa = get_position_value(fw_seq, fw_pos)
                    if fw_aa in AA_TO_IDX:
                        # Per-family
                        self.family_single_counts[family][cdr][cdr_pos][fw][fw_pos][cdr_aa][fw_aa] += 1
                        # Global
                        self.global_single_counts[cdr][cdr_pos][fw][fw_pos][cdr_aa][fw_aa] += 1
            
            # Triplet correlations
            if len(cdr_seq) >= 3 and len(fw_seq) >= 3:
                cdr_triplet = cdr_seq[-3:]
                fw_triplet = fw_seq[:3]
                
                if (all(aa in AA_TO_IDX for aa in cdr_triplet) and 
                    all(aa in AA_TO_IDX for aa in fw_triplet)):
                    self.family_triplet_counts[family][cdr][fw][cdr_triplet][fw_triplet] += 1
    
    def compute_results(self, min_count=50, min_confidence=75.0):
        results = {
            'n_sequences': self.n_sequences,
            'germline_counts': dict(self.germline_counts),
            'motif_counts': dict(self.motif_counts),
            'hallmark_distributions': {k: dict(v) for k, v in self.hallmark_distributions.items()},
            
            # Per-family analysis
            'family_conservation': {},
            'family_cdr_dependent': {},
            'family_triplet_rules': {},
            'family_cdr_length_stats': {},
            
            # Global analysis
            'global_rules': [],
            'true_cdr_dependent': [],  # Rules that hold WITHIN families
            
            # Cross-family comparison
            'germline_dependent_positions': [],  # Positions that differ BETWEEN families
        }
        
        # === Germline family summaries ===
        for family in self.germline_counts:
            # CDR length stats
            length_stats = {}
            for cdr in ['cdr1', 'cdr2', 'cdr3']:
                lengths = self.family_cdr_lengths[family][cdr]
                if lengths:
                    length_stats[cdr] = {
                        'mean': np.mean(lengths),
                        'std': np.std(lengths),
                        'min': min(lengths),
                        'max': max(lengths),
                        'n': len(lengths),
                    }
            results['family_cdr_length_stats'][family] = length_stats
            
            # FW conservation per family
            family_cons = {}
            for fw in ['FR1', 'FR2', 'FR3', 'FR4']:
                family_cons[fw] = {}
                for fw_pos in self.family_fw_marginals[family][fw]:
                    counts = self.family_fw_marginals[family][fw][fw_pos]
                    total = sum(counts.values())
                    if total < min_count:
                        continue
                    best_aa = max(counts, key=counts.get)
                    best_pct = 100.0 * counts[best_aa] / total
                    family_cons[fw][fw_pos] = {
                        'best_aa': best_aa,
                        'best_pct': best_pct,
                        'total': total,
                        'distribution': dict(counts),
                    }
            results['family_conservation'][family] = family_cons
        
        # === Find TRUE CDR-dependent positions (within each family) ===
        for family in self.family_single_counts:
            family_cdr_dep = {}
            
            for cdr in self.family_single_counts[family]:
                for cdr_pos in self.family_single_counts[family][cdr]:
                    for fw in self.family_single_counts[family][cdr][cdr_pos]:
                        for fw_pos in self.family_single_counts[family][cdr][cdr_pos][fw]:
                            joint = self.family_single_counts[family][cdr][cdr_pos][fw][fw_pos]
                            
                            # Collect predictions for this FW position
                            predictions = {}
                            for cdr_aa in joint:
                                total = sum(joint[cdr_aa].values())
                                if total < min_count:
                                    continue
                                best_fw_aa = max(joint[cdr_aa], key=joint[cdr_aa].get)
                                best_count = joint[cdr_aa][best_fw_aa]
                                confidence = 100.0 * best_count / total
                                
                                if confidence >= min_confidence:
                                    key = f"{cdr}[{cdr_pos}]={cdr_aa}"
                                    predictions[key] = {
                                        'fw_aa': best_fw_aa,
                                        'confidence': confidence,
                                        'n': total,
                                    }
                            
                            # Check if different CDRs predict different FW AAs
                            if predictions:
                                predicted_aas = set(p['fw_aa'] for p in predictions.values())
                                if len(predicted_aas) > 1:
                                    # TRUE CDR-dependence within this family!
                                    pos_key = f"{fw}[{fw_pos}]"
                                    if pos_key not in family_cdr_dep:
                                        family_cdr_dep[pos_key] = {
                                            'variants': predicted_aas,
                                            'rules': {}
                                        }
                                    family_cdr_dep[pos_key]['rules'].update(predictions)
            
            results['family_cdr_dependent'][family] = family_cdr_dep
        
        # === Find GERMLINE-dependent positions (differ between families) ===
        # For each FW position, check if different families have different consensus AAs
        for fw in ['FR1', 'FR2', 'FR3', 'FR4']:
            for fw_pos in FR_POSITIONS[fw]:
                family_consensus = {}
                for family in results['family_conservation']:
                    cons = results['family_conservation'][family].get(fw, {}).get(fw_pos)
                    if cons and cons['best_pct'] >= 75:
                        family_consensus[family] = cons['best_aa']
                
                if len(family_consensus) >= 2:
                    unique_aas = set(family_consensus.values())
                    if len(unique_aas) > 1:
                        results['germline_dependent_positions'].append({
                            'position': f"{fw}[{fw_pos}]",
                            'family_consensus': family_consensus,
                            'variants': unique_aas,
                        })
        
        # === Family triplet rules ===
        for family in self.family_triplet_counts:
            family_rules = []
            for cdr in self.family_triplet_counts[family]:
                for fw in self.family_triplet_counts[family][cdr]:
                    for cdr_triplet in self.family_triplet_counts[family][cdr][fw]:
                        counts = self.family_triplet_counts[family][cdr][fw][cdr_triplet]
                        total = sum(counts.values())
                        if total < min_count:
                            continue
                        
                        best_fw_triplet = max(counts, key=counts.get)
                        confidence = 100.0 * counts[best_fw_triplet] / total
                        
                        if confidence >= min_confidence:
                            family_rules.append({
                                'cdr': cdr,
                                'cdr_triplet': cdr_triplet,
                                'fw': fw,
                                'fw_triplet': best_fw_triplet,
                                'confidence': confidence,
                                'n': total,
                            })
            
            family_rules.sort(key=lambda x: x['confidence'], reverse=True)
            results['family_triplet_rules'][family] = family_rules
        
        return results

# ============================================================
# MAIN
# ============================================================

def format_time(seconds):
    if seconds < 60:
        return f"{seconds:.1f}s"
    elif seconds < 3600:
        return f"{seconds/60:.1f}m"
    else:
        return f"{seconds/3600:.1f}h"

def main():
    parser = argparse.ArgumentParser(description='Germline-Aware CDR-Framework Analysis v3')
    parser.add_argument('npz_files', nargs='+', help='NPZ shard files')
    parser.add_argument('--output', '-o', required=True, help='Output directory')
    parser.add_argument('--max-shards', type=int, default=None, help='Max shards to process')
    parser.add_argument('--min-count', type=int, default=50, help='Min count for rules')
    parser.add_argument('--min-confidence', type=float, default=75.0, help='Min confidence %%')
    
    args = parser.parse_args()
    
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    npz_files = sorted(args.npz_files)
    if args.max_shards:
        npz_files = npz_files[:args.max_shards]
    
    print()
    print(f"{'='*70}")
    print(f"  GERMLINE-AWARE CDR-Framework Analysis v3")
    print(f"{'='*70}")
    print(f"  Shards: {len(npz_files)}")
    print(f"  Min confidence: {args.min_confidence}%")
    print(f"  Output: {output_dir}")
    print(f"{'='*70}")
    print()
    print("  Features:")
    print("    • Germline classification (FR2 hallmarks: 42/49/50/52)")
    print("    • Cysteine counting (C2 vs C4)")
    print("    • Family grouping (Y_C2, F_C2, F_C4, VH-like, etc.)")
    print("    • Per-family CDR-FW correlation analysis")
    print("    • TRUE CDR-dependence (within family)")
    print("    • Germline-dependence (between families)")
    print()
    
    overall_start = time.time()
    shard_times = []
    
    accumulator = GermlineAwareAccumulator()
    
    print("STEP 1: Processing shards")
    print("-" * 70)
    
    total_valid = 0
    total_seqs = 0
    
    for i, npz_file in enumerate(npz_files):
        shard_start = time.time()
        shard_name = Path(npz_file).name
        
        print(f"\n  [{i+1}/{len(npz_files)}] {shard_name}")
        
        n_valid, n_total = process_shard_streaming(npz_file, accumulator, show_progress=True)
        
        shard_time = time.time() - shard_start
        shard_times.append(shard_time)
        
        total_valid += n_valid
        total_seqs += n_total
        
        avg_time = sum(shard_times) / len(shard_times)
        remaining = len(npz_files) - (i + 1)
        eta = avg_time * remaining
        
        print(f"    ✓ {n_valid:,}/{n_total:,} valid | Time: {format_time(shard_time)} | ETA: {format_time(eta)}")
        
        gc.collect()
    
    step1_time = time.time() - overall_start
    
    print(f"\n  Total: {total_valid:,} valid sequences from {total_seqs:,}")
    print(f"  Step 1 time: {format_time(step1_time)}")
    
    # Print germline distribution
    print()
    print("-" * 70)
    print("GERMLINE DISTRIBUTION")
    print("-" * 70)
    for family, count in sorted(accumulator.germline_counts.items(), key=lambda x: -x[1]):
        pct = 100.0 * count / total_valid
        print(f"  {family:20s}: {count:>10,} ({pct:5.1f}%)")
    
    # Print hallmark distribution
    print()
    print("-" * 70)
    print("FR2 HALLMARK DISTRIBUTIONS")
    print("-" * 70)
    for pos in ['pos42', 'pos49', 'pos50', 'pos52']:
        dist = accumulator.hallmark_distributions[pos]
        total = sum(dist.values())
        print(f"\n  Position {pos.replace('pos', '')} (IMGT):")
        for aa, count in sorted(dist.items(), key=lambda x: -x[1])[:5]:
            pct = 100.0 * count / total
            print(f"    {aa}: {count:>10,} ({pct:5.1f}%)")
    
    # Compute results
    print()
    print("-" * 70)
    print("STEP 2: Computing correlations")
    print("-" * 70)
    
    step2_start = time.time()
    results = accumulator.compute_results(
        min_count=args.min_count,
        min_confidence=args.min_confidence
    )
    step2_time = time.time() - step2_start
    
    # Summary
    n_germline_dep = len(results['germline_dependent_positions'])
    total_true_cdr_dep = sum(len(v) for v in results['family_cdr_dependent'].values())
    
    print(f"\n  ✓ Germline-dependent positions: {n_germline_dep}")
    print(f"  ✓ True CDR-dependent positions (across all families): {total_true_cdr_dep}")
    print(f"  Step 2 time: {format_time(step2_time)}")
    
    # Print germline-dependent positions
    print()
    print("-" * 70)
    print("GERMLINE-DEPENDENT POSITIONS (differ between families)")
    print("-" * 70)
    for pos_info in results['germline_dependent_positions'][:15]:
        print(f"\n  {pos_info['position']}: {pos_info['variants']}")
        for fam, aa in sorted(pos_info['family_consensus'].items()):
            print(f"    {fam}: {aa}")
    
    # Print true CDR-dependent per family
    print()
    print("-" * 70)
    print("TRUE CDR-DEPENDENT POSITIONS (within families)")
    print("-" * 70)
    for family in sorted(results['family_cdr_dependent'].keys()):
        cdr_dep = results['family_cdr_dependent'][family]
        if cdr_dep:
            print(f"\n  {family}:")
            for pos, info in list(cdr_dep.items())[:5]:
                print(f"    {pos}: {info['variants']}")
    
    # Save results
    print()
    print("-" * 70)
    print("STEP 3: Saving results")
    print("-" * 70)
    
    corr_dir = output_dir / "1_correlations"
    corr_dir.mkdir(exist_ok=True)
    
    with open(corr_dir / "correlation_results_v3.pkl", 'wb') as f:
        pickle.dump(results, f)
    
    # Write comprehensive summary
    with open(corr_dir / "correlation_summary_v3.txt", 'w') as f:
        f.write(f"Germline-Aware CDR-Framework Analysis v3\n")
        f.write(f"=" * 70 + "\n\n")
        f.write(f"Total sequences: {total_valid:,}\n")
        f.write(f"Min confidence: {args.min_confidence}%\n\n")
        
        f.write(f"GERMLINE DISTRIBUTION\n")
        f.write(f"-" * 50 + "\n")
        for family, count in sorted(results['germline_counts'].items(), key=lambda x: -x[1]):
            pct = 100.0 * count / total_valid
            f.write(f"  {family:20s}: {count:>10,} ({pct:5.1f}%)\n")
        
        f.write(f"\nFR2 MOTIF DISTRIBUTION (top 20)\n")
        f.write(f"-" * 50 + "\n")
        for motif, count in sorted(results['motif_counts'].items(), key=lambda x: -x[1])[:20]:
            pct = 100.0 * count / total_valid
            f.write(f"  {motif}: {count:>10,} ({pct:5.1f}%)\n")
        
        f.write(f"\nGERMLINE-DEPENDENT POSITIONS\n")
        f.write(f"-" * 50 + "\n")
        f.write(f"(These positions differ BETWEEN germline families - not CDR-dependent)\n\n")
        for pos_info in results['germline_dependent_positions']:
            f.write(f"  {pos_info['position']}: {pos_info['variants']}\n")
            for fam, aa in sorted(pos_info['family_consensus'].items()):
                f.write(f"    {fam}: {aa}\n")
            f.write("\n")
        
        f.write(f"\nTRUE CDR-DEPENDENT POSITIONS (per family)\n")
        f.write(f"-" * 50 + "\n")
        f.write(f"(These positions differ based on CDR WITHIN the same germline family)\n\n")
        for family in sorted(results['family_cdr_dependent'].keys()):
            cdr_dep = results['family_cdr_dependent'][family]
            if cdr_dep:
                f.write(f"  {family}:\n")
                for pos, info in cdr_dep.items():
                    f.write(f"    {pos}: {info['variants']}\n")
                    for rule, details in list(info['rules'].items())[:3]:
                        f.write(f"      {rule} → {details['fw_aa']} ({details['confidence']:.1f}%)\n")
                f.write("\n")
        
        f.write(f"\nCDR LENGTH STATISTICS BY FAMILY\n")
        f.write(f"-" * 50 + "\n")
        for family in sorted(results['family_cdr_length_stats'].keys()):
            stats = results['family_cdr_length_stats'][family]
            f.write(f"\n  {family}:\n")
            for cdr in ['cdr1', 'cdr2', 'cdr3']:
                if cdr in stats:
                    s = stats[cdr]
                    f.write(f"    {cdr}: mean={s['mean']:.1f} ± {s['std']:.1f}, "
                           f"range={s['min']}-{s['max']}, n={s['n']:,}\n")
    
    print(f"  ✓ Saved: {corr_dir}/correlation_results_v3.pkl")
    print(f"  ✓ Saved: {corr_dir}/correlation_summary_v3.txt")
    
    # Final report
    total_time = time.time() - overall_start
    
    print()
    print(f"{'='*70}")
    print(f"  COMPLETE")
    print(f"{'='*70}")
    print(f"  Total sequences:  {total_valid:,}")
    print(f"  Germline families: {len(results['germline_counts'])}")
    print(f"  Germline-dependent positions: {n_germline_dep}")
    print(f"  True CDR-dependent positions: {total_true_cdr_dep}")
    print(f"  Total time: {format_time(total_time)}")
    print(f"  Output: {output_dir}")
    print(f"{'='*70}")
    print()

if __name__ == '__main__':
    main()
