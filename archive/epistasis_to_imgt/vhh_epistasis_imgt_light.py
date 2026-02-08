#!/usr/bin/env python3
"""
VHH IMGT Epistasis Analysis - LIGHTWEIGHT VERSION
==================================================

Skips heavy computations to run faster and use less memory:
- SKIPS: Random Forest models (analysis_3) - saves ~10GB RAM, hours of compute
- SKIPS: Mutual Information matrices (analysis_5) - saves ~5GB RAM
- KEEPS: Compensation stats (analysis_1) - lightweight streaming
- KEEPS: Vernier clusters (analysis_2) - essential for designer
- KEEPS: Higher-order rules (analysis_4) - essential for designer

Expected runtime: ~2-4 hours for 12M sequences (vs 2-5 days for full)
Expected memory: ~4GB peak (vs 20GB+ for full)

Usage:
    python vhh_epistasis_imgt_light.py \
        --npz data/processed/full_12m/full_12m_imgt.npz \
        --output models/epistasis/imgt_light_v1 \
        --checkpoint-interval 1000000
"""

import numpy as np
import pickle
import os
import sys
import re
import gc
from collections import defaultdict, Counter
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional
import time
import argparse

try:
    from tqdm import tqdm
    HAS_TQDM = True
except ImportError:
    HAS_TQDM = False
    def tqdm(iterable, desc=None, total=None, **kwargs):
        if desc:
            print(f"{desc}...")
        return iterable

# ============================================================
# AMINO ACID PROPERTIES
# ============================================================

AA_PROPERTIES = {
    'A': {'hydrophobic': 1, 'volume': 88.6, 'charge': 0, 'aromatic': 0},
    'R': {'hydrophobic': 0, 'volume': 173.4, 'charge': 1, 'aromatic': 0},
    'N': {'hydrophobic': 0, 'volume': 114.1, 'charge': 0, 'aromatic': 0},
    'D': {'hydrophobic': 0, 'volume': 111.1, 'charge': -1, 'aromatic': 0},
    'C': {'hydrophobic': 1, 'volume': 108.5, 'charge': 0, 'aromatic': 0},
    'Q': {'hydrophobic': 0, 'volume': 143.8, 'charge': 0, 'aromatic': 0},
    'E': {'hydrophobic': 0, 'volume': 138.4, 'charge': -1, 'aromatic': 0},
    'G': {'hydrophobic': 0, 'volume': 60.1, 'charge': 0, 'aromatic': 0},
    'H': {'hydrophobic': 0, 'volume': 153.2, 'charge': 0.5, 'aromatic': 1},
    'I': {'hydrophobic': 1, 'volume': 166.7, 'charge': 0, 'aromatic': 0},
    'L': {'hydrophobic': 1, 'volume': 166.7, 'charge': 0, 'aromatic': 0},
    'K': {'hydrophobic': 0, 'volume': 168.6, 'charge': 1, 'aromatic': 0},
    'M': {'hydrophobic': 1, 'volume': 162.9, 'charge': 0, 'aromatic': 0},
    'F': {'hydrophobic': 1, 'volume': 189.9, 'charge': 0, 'aromatic': 1},
    'P': {'hydrophobic': 0, 'volume': 112.7, 'charge': 0, 'aromatic': 0},
    'S': {'hydrophobic': 0, 'volume': 89.0, 'charge': 0, 'aromatic': 0},
    'T': {'hydrophobic': 0, 'volume': 116.1, 'charge': 0, 'aromatic': 0},
    'W': {'hydrophobic': 1, 'volume': 227.8, 'charge': 0, 'aromatic': 1},
    'Y': {'hydrophobic': 0, 'volume': 193.6, 'charge': 0, 'aromatic': 1},
    'V': {'hydrophobic': 1, 'volume': 140.0, 'charge': 0, 'aromatic': 0},
}

# IMGT positions for VHH hallmarks
HALLMARK_POSITIONS = [42, 49, 50, 52]

# Vernier positions (IMGT numbering)
VERNIER_IMGT = {
    'FR2': [39, 45, 47, 49],  # Positions around CDR1/CDR2
    'FR3': [66, 67, 69, 71, 73, 76, 78, 80, 82, 84, 85, 87, 89, 91, 93, 94],
}

def get_cdr_features(cdr_seq: str) -> dict:
    """Extract features from a CDR sequence."""
    if not cdr_seq:
        return {}
    
    length = len(cdr_seq)
    n_cys = cdr_seq.count('C')
    charge = sum(AA_PROPERTIES.get(aa, {}).get('charge', 0) for aa in cdr_seq)
    hydrophobic = sum(1 for aa in cdr_seq if AA_PROPERTIES.get(aa, {}).get('hydrophobic', 0))
    aromatic = sum(1 for aa in cdr_seq if AA_PROPERTIES.get(aa, {}).get('aromatic', 0))
    
    return {
        'length': length,
        'n_cys': n_cys,
        'charge': charge,
        'frac_hydrophobic': hydrophobic / length if length > 0 else 0,
        'frac_aromatic': aromatic / length if length > 0 else 0,
    }

def categorize_cdr3(cdr3: str) -> dict:
    """Categorize CDR3 for rule conditions."""
    length = len(cdr3)
    charge = sum(AA_PROPERTIES.get(aa, {}).get('charge', 0) for aa in cdr3)
    n_cys = cdr3.count('C')
    
    # Length category
    if length < 10:
        len_cat = 'short'
    elif length < 14:
        len_cat = 'medium_short'
    elif length < 18:
        len_cat = 'medium'
    elif length < 22:
        len_cat = 'medium_long'
    else:
        len_cat = 'long'
    
    # Charge category
    if charge < -2:
        chg_cat = 'very_neg'
    elif charge < -0.5:
        chg_cat = 'negative'
    elif charge < 1:
        chg_cat = 'neutral'
    elif charge < 3:
        chg_cat = 'positive'
    else:
        chg_cat = 'very_pos'
    
    # Cysteine
    cys_cat = 'with_cys' if n_cys > 0 else 'no_cys'
    
    return {
        'cdr3_len': len_cat,
        'cdr3_charge': chg_cat,
        'cdr3_cys': cys_cat,
    }

def classify_family(imgt_dict: dict, full_seq: str) -> str:
    """Classify VHH family based on IMGT positions."""
    pos42 = imgt_dict.get(42, '-')
    pos49 = imgt_dict.get(49, '-')
    pos50 = imgt_dict.get(50, '-')
    
    n_cys = full_seq.count('C') if full_seq else 2
    
    # VH-like
    if pos50 == 'L':
        return 'VH_like'
    
    # Classical VHH
    if pos42 in ['F', 'Y'] and pos49 in ['E', 'Q'] and pos50 == 'R':
        if pos42 == 'Y' and n_cys == 2:
            return 'Y_C2'
        elif pos42 == 'F' and n_cys == 2:
            return 'F_C2'
        elif pos42 == 'F' and n_cys == 4:
            return 'F_C4'
        elif pos42 == 'Y' and n_cys == 4:
            return 'Y_C4'
        else:
            return 'Classical_other'
    
    return 'Non_classical'


# ============================================================
# ANALYSIS 1: COMPENSATION STATISTICS (Streaming)
# ============================================================

class StreamingCompensationAnalyzer:
    """
    Memory-efficient compensation analysis using Welford's online algorithm.
    Computes mean/variance of CDR features for each FR residue.
    """
    
    def __init__(self):
        self.stats = defaultdict(lambda: defaultdict(lambda: {
            'count': 0, 'mean': 0.0, 'M2': 0.0
        }))
        self.total = 0
    
    def update(self, imgt_positions: dict, cdr1: str, cdr2: str, cdr3: str):
        """Update stats with one sequence."""
        self.total += 1
        
        # Get CDR features
        cdr1_features = get_cdr_features(cdr1)
        cdr2_features = get_cdr_features(cdr2)
        cdr3_features = get_cdr_features(cdr3)
        
        # For each FR position, update stats
        for imgt_pos, aa in imgt_positions.items():
            if not isinstance(imgt_pos, int):
                continue
            if imgt_pos < 1 or imgt_pos > 128:
                continue
            
            # Determine region
            if 1 <= imgt_pos <= 26:
                region = 'FR1'
            elif 39 <= imgt_pos <= 55:
                region = 'FR2'
            elif 66 <= imgt_pos <= 104:
                region = 'FR3'
            elif 118 <= imgt_pos <= 128:
                region = 'FR4'
            else:
                continue
            
            key = f"IMGT{imgt_pos}"
            
            # Update running stats for each CDR feature
            for feature_name, value in [
                ('cdr1_len', cdr1_features.get('length', 0)),
                ('cdr1_charge', cdr1_features.get('charge', 0)),
                ('cdr2_len', cdr2_features.get('length', 0)),
                ('cdr3_len', cdr3_features.get('length', 0)),
                ('cdr3_charge', cdr3_features.get('charge', 0)),
            ]:
                self._update_welford(self.stats[key][aa], feature_name, value)
    
    def _update_welford(self, stat_dict, feature, value):
        """Welford's online algorithm for mean/variance."""
        if feature not in stat_dict:
            stat_dict[feature] = {'count': 0, 'mean': 0.0, 'M2': 0.0}
        
        s = stat_dict[feature]
        s['count'] += 1
        delta = value - s['mean']
        s['mean'] += delta / s['count']
        delta2 = value - s['mean']
        s['M2'] += delta * delta2
    
    def get_results(self) -> dict:
        """Finalize and return results."""
        results = {}
        for pos_key, aa_dict in self.stats.items():
            results[pos_key] = {}
            for aa, features in aa_dict.items():
                results[pos_key][aa] = {}
                for feat, s in features.items():
                    if isinstance(s, dict) and 'count' in s:
                        variance = s['M2'] / s['count'] if s['count'] > 1 else 0
                        results[pos_key][aa][feat] = {
                            'mean': s['mean'],
                            'std': variance ** 0.5,
                            'count': s['count'],
                        }
        return results


# ============================================================
# ANALYSIS 2: VERNIER CLUSTERING
# ============================================================

class VernierClusterAnalyzer:
    """
    Cluster sequences by Vernier zone residues.
    Uses IMGT positions for consistent numbering.
    """
    
    def __init__(self):
        # Key vernier positions (IMGT)
        self.vernier_positions = [42, 49, 50, 52] + [66, 71, 76, 82, 87, 89, 91, 94]
        self.clusters = defaultdict(lambda: {
            'count': 0,
            'families': Counter(),
            'cdr3_lengths': [],
            'cdr3_charges': [],
        })
        self.max_cdr3_samples = 10000  # Limit memory for CDR3 stats
    
    def update(self, imgt_positions: dict, family: str, cdr3: str):
        """Add sequence to cluster."""
        # Build vernier signature
        sig_parts = []
        for pos in self.vernier_positions:
            aa = imgt_positions.get(pos, '-')
            sig_parts.append(f"{pos}{aa}")
        
        sig = "_".join(sig_parts)
        
        cluster = self.clusters[sig]
        cluster['count'] += 1
        cluster['families'][family] += 1
        
        # Sample CDR3 stats (don't store all)
        if len(cluster['cdr3_lengths']) < self.max_cdr3_samples:
            cluster['cdr3_lengths'].append(len(cdr3))
            charge = sum(AA_PROPERTIES.get(aa, {}).get('charge', 0) for aa in cdr3)
            cluster['cdr3_charges'].append(charge)
    
    def get_results(self) -> dict:
        """Finalize clusters."""
        results = {}
        
        for sig, data in self.clusters.items():
            if data['count'] < 100:  # Skip rare clusters
                continue
            
            # Parse signature back to pattern
            pattern = {}
            for part in sig.split('_'):
                if len(part) >= 2:
                    pos = int(part[:-1])
                    aa = part[-1]
                    pattern[pos] = aa
            
            # Compute CDR3 stats
            cdr3_stats = {}
            if data['cdr3_lengths']:
                cdr3_stats = {
                    'length': {
                        'mean': np.mean(data['cdr3_lengths']),
                        'std': np.std(data['cdr3_lengths']),
                    },
                    'charge': {
                        'mean': np.mean(data['cdr3_charges']),
                        'std': np.std(data['cdr3_charges']),
                    }
                }
            
            # Determine dominant family
            families = dict(data['families'])
            dominant = max(families, key=families.get) if families else 'Unknown'
            
            results[sig] = {
                'pattern': pattern,
                'count': data['count'],
                'families': families,
                'dominant_family': dominant,
                'cdr3_stats': cdr3_stats,
            }
        
        return results


# ============================================================
# ANALYSIS 4: HIGHER-ORDER RULES (Streaming)
# ============================================================

class StreamingRuleFinder:
    """
    Find conditional rules: "If CDR3 is short AND FR2[pos]=X, then FR3[pos]=Y"
    Uses streaming counts to minimize memory.
    """
    
    def __init__(self, min_support: int = 1000, min_confidence: float = 0.7):
        self.min_support = min_support
        self.min_confidence = min_confidence
        
        # Count tables: condition -> consequence -> count
        self.counts = defaultdict(lambda: defaultdict(int))
        self.condition_totals = defaultdict(int)
        self.total = 0
    
    def update(self, imgt_positions: dict, cdr3: str, family: str):
        """Update counts with one sequence."""
        self.total += 1
        
        cdr3_cats = categorize_cdr3(cdr3)
        
        # Build conditions from CDR3 categories + some FR positions
        conditions = []
        
        # CDR3 conditions
        for cat_name, cat_val in cdr3_cats.items():
            conditions.append(f"{cat_name}={cat_val}")
        
        # FR2 conditions (key positions)
        for pos in [42, 49, 50]:
            aa = imgt_positions.get(pos)
            if aa:
                conditions.append(f"IMGT{pos}={aa}")
        
        # For each FR3 position, record consequence given conditions
        for target_pos in range(66, 105):
            target_aa = imgt_positions.get(target_pos)
            if not target_aa:
                continue
            
            consequence = f"IMGT{target_pos}={target_aa}"
            
            # Single conditions
            for cond in conditions:
                self.counts[cond][consequence] += 1
                self.condition_totals[cond] += 1
            
            # Pairs of conditions (limited to avoid explosion)
            for i in range(min(len(conditions), 3)):
                for j in range(i + 1, min(len(conditions), 4)):
                    cond = f"{conditions[i]} AND {conditions[j]}"
                    self.counts[cond][consequence] += 1
                    self.condition_totals[cond] += 1
    
    def get_rules(self) -> List[dict]:
        """Extract high-confidence rules."""
        rules = []
        
        for condition, consequences in self.counts.items():
            total = self.condition_totals[condition]
            if total < self.min_support:
                continue
            
            for consequence, count in consequences.items():
                confidence = count / total
                if confidence >= self.min_confidence:
                    rules.append({
                        'condition': condition,
                        'result': consequence,
                        'support': count,
                        'confidence': confidence,
                        'condition_total': total,
                    })
        
        # Sort by confidence
        rules.sort(key=lambda x: -x['confidence'])
        return rules


# ============================================================
# MAIN PROCESSOR
# ============================================================

class LightEpistasisProcessor:
    """
    Lightweight epistasis processor - skips RF and MI.
    """
    
    def __init__(self, output_dir: str, checkpoint_interval: int = 500000):
        self.output_dir = output_dir
        self.checkpoint_interval = checkpoint_interval
        os.makedirs(output_dir, exist_ok=True)
        
        # Initialize analyzers
        self.compensation = StreamingCompensationAnalyzer()
        self.vernier = VernierClusterAnalyzer()
        self.rules = StreamingRuleFinder()
        
        self.processed = 0
        self.family_counts = Counter()
        self.start_time = time.time()
    
    def process_sequence(self, imgt_dict: dict, cdr1: str, cdr2: str, cdr3: str, 
                         full_seq: str):
        """Process one sequence through all analyzers."""
        self.processed += 1
        
        # Classify family
        family = classify_family(imgt_dict, full_seq)
        self.family_counts[family] += 1
        
        # Update analyzers
        self.compensation.update(imgt_dict, cdr1, cdr2, cdr3)
        self.vernier.update(imgt_dict, family, cdr3)
        self.rules.update(imgt_dict, cdr3, family)
        
        # Checkpoint
        if self.processed % self.checkpoint_interval == 0:
            self._save_checkpoint()
    
    def _save_checkpoint(self):
        """Save intermediate results."""
        elapsed = time.time() - self.start_time
        rate = self.processed / elapsed if elapsed > 0 else 0
        print(f"  Checkpoint at {self.processed:,} sequences ({rate:.0f}/sec)")
        
        checkpoint = {
            'processed': self.processed,
            'family_counts': dict(self.family_counts),
            'elapsed_seconds': elapsed,
        }
        
        path = os.path.join(self.output_dir, f'checkpoint_{self.processed}.pkl')
        with open(path, 'wb') as f:
            pickle.dump(checkpoint, f)
    
    def finalize(self) -> dict:
        """Finalize and return all results."""
        elapsed = time.time() - self.start_time
        print(f"\nFinalizing {self.processed:,} sequences...")
        
        results = {
            'analysis_1_compensation': self.compensation.get_results(),
            'analysis_2_vernier_clusters': self.vernier.get_results(),
            'analysis_3_conditional_models': {},  # SKIPPED
            'analysis_4_higher_order_rules': self.rules.get_rules(),
            'analysis_5_mutual_information': {},  # SKIPPED
            'family_counts': dict(self.family_counts),
            'total_sequences': self.processed,
            'elapsed_seconds': elapsed,
            'parameters': {
                'mode': 'imgt_light',
                'skipped': ['RF_models', 'MI_matrix'],
            }
        }
        
        return results


def extract_regions_from_sequence(aa_v_full: str, cdr1: str, cdr2: str, cdr3: str) -> dict:
    """
    Extract FR1, FR2, FR3, FR4 and build position dict from full sequence + CDRs.
    
    Returns dict with:
        - fr1, fr2, fr3, fr4: framework sequences
        - cdr1, cdr2, cdr3: CDR sequences  
        - positions: dict mapping IMGT position -> amino acid
    """
    result = {
        'fr1': '', 'fr2': '', 'fr3': '', 'fr4': '',
        'cdr1': cdr1, 'cdr2': cdr2, 'cdr3': cdr3,
        'positions': {},
        'valid': False,
    }
    
    if not aa_v_full or not cdr1 or not cdr2 or not cdr3:
        return result
    
    # Find CDR positions in full sequence
    try:
        cdr1_start = aa_v_full.find(cdr1)
        if cdr1_start == -1:
            return result
        cdr1_end = cdr1_start + len(cdr1)
        
        cdr2_start = aa_v_full.find(cdr2, cdr1_end)
        if cdr2_start == -1:
            return result
        cdr2_end = cdr2_start + len(cdr2)
        
        cdr3_start = aa_v_full.find(cdr3, cdr2_end)
        if cdr3_start == -1:
            return result
        cdr3_end = cdr3_start + len(cdr3)
        
        # Extract framework regions
        fr1 = aa_v_full[:cdr1_start]
        fr2 = aa_v_full[cdr1_end:cdr2_start]
        fr3 = aa_v_full[cdr2_end:cdr3_start]
        fr4 = aa_v_full[cdr3_end:]
        
        # Validate lengths (typical VHH ranges)
        if not (15 <= len(fr1) <= 30):
            return result
        if not (10 <= len(fr2) <= 20):
            return result
        if not (30 <= len(fr3) <= 45):
            return result
        if not (8 <= len(fr4) <= 15):
            return result
        
        result['fr1'] = fr1
        result['fr2'] = fr2
        result['fr3'] = fr3
        result['fr4'] = fr4
        result['valid'] = True
        
        # Build IMGT position dict
        # Standard IMGT numbering for VHH:
        # FR1: 1-26, CDR1: 27-38, FR2: 39-55, CDR2: 56-65, FR3: 66-104, CDR3: 105-117, FR4: 118-128
        positions = {}
        
        # FR1 (IMGT 1-26)
        for i, aa in enumerate(fr1):
            imgt_pos = i + 1  # 1-based
            if imgt_pos <= 26:
                positions[imgt_pos] = aa
        
        # CDR1 (IMGT 27-38) - variable length, center-aligned
        cdr1_imgt_start = 27
        for i, aa in enumerate(cdr1):
            positions[cdr1_imgt_start + i] = aa
        
        # FR2 (IMGT 39-55) - key hallmark region
        fr2_imgt_start = 39
        for i, aa in enumerate(fr2):
            imgt_pos = fr2_imgt_start + i
            if imgt_pos <= 55:
                positions[imgt_pos] = aa
        
        # CDR2 (IMGT 56-65)
        cdr2_imgt_start = 56
        for i, aa in enumerate(cdr2):
            positions[cdr2_imgt_start + i] = aa
        
        # FR3 (IMGT 66-104)
        fr3_imgt_start = 66
        for i, aa in enumerate(fr3):
            imgt_pos = fr3_imgt_start + i
            if imgt_pos <= 104:
                positions[imgt_pos] = aa
        
        # CDR3 (IMGT 105-117) - highly variable
        cdr3_imgt_start = 105
        for i, aa in enumerate(cdr3):
            positions[cdr3_imgt_start + i] = aa
        
        # FR4 (IMGT 118-128)
        fr4_imgt_start = 118
        for i, aa in enumerate(fr4):
            imgt_pos = fr4_imgt_start + i
            if imgt_pos <= 128:
                positions[imgt_pos] = aa
        
        result['positions'] = positions
        
    except Exception:
        return result
    
    return result


def load_npz_shard(path: str) -> Tuple[List[dict], dict]:
    """
    Load a single NPZ shard and return list of sequence records.
    
    Each record is a dict with: aa_v_full, cdr1, cdr2, cdr3, id
    """
    data = np.load(path, allow_pickle=True)
    
    records = []
    
    # Check what columns we have
    has_aa_v_full = 'aa_v_full' in data.files
    has_cdrs = 'cdr1' in data.files and 'cdr2' in data.files and 'cdr3' in data.files
    has_sequences = 'sequences' in data.files
    
    if has_aa_v_full and has_cdrs:
        # Standard format: aa_v_full + cdr1/cdr2/cdr3
        aa_seqs = data['aa_v_full']
        cdr1s = data['cdr1']
        cdr2s = data['cdr2']
        cdr3s = data['cdr3']
        ids = data['ids'] if 'ids' in data.files else [f'seq_{i}' for i in range(len(aa_seqs))]
        
        for i in range(len(aa_seqs)):
            records.append({
                'aa_v_full': str(aa_seqs[i]) if aa_seqs[i] is not None else '',
                'cdr1': str(cdr1s[i]) if cdr1s[i] is not None else '',
                'cdr2': str(cdr2s[i]) if cdr2s[i] is not None else '',
                'cdr3': str(cdr3s[i]) if cdr3s[i] is not None else '',
                'id': str(ids[i]) if i < len(ids) else f'seq_{i}',
            })
    
    elif has_sequences:
        # Simple format: just sequences array
        seqs = data['sequences']
        for i, seq in enumerate(seqs):
            if isinstance(seq, dict):
                records.append(seq)
            else:
                records.append({'aa_v_full': str(seq), 'cdr1': '', 'cdr2': '', 'cdr3': ''})
    
    else:
        # Try first array
        key = data.files[0]
        arr = data[key]
        for i, item in enumerate(arr):
            if isinstance(item, dict):
                records.append(item)
            elif isinstance(item, str):
                records.append({'aa_v_full': item, 'cdr1': '', 'cdr2': '', 'cdr3': ''})
    
    metadata = {}
    if 'metadata' in data.files:
        metadata = data['metadata'].item() if hasattr(data['metadata'], 'item') else {}
    
    return records, metadata


def main():
    parser = argparse.ArgumentParser(description='Lightweight IMGT Epistasis Analysis')
    parser.add_argument('--npz', '-n', required=True, help='NPZ file or glob pattern')
    parser.add_argument('--output', '-o', required=True, help='Output directory')
    parser.add_argument('--max-sequences', '-m', type=int, default=None, help='Max sequences to process')
    parser.add_argument('--checkpoint-interval', '-c', type=int, default=500000)
    args = parser.parse_args()
    
    print("=" * 70)
    print("VHH IMGT Epistasis Analysis - LIGHTWEIGHT VERSION")
    print("=" * 70)
    print(f"Input: {args.npz}")
    print(f"Output: {args.output}")
    print(f"Checkpoint every: {args.checkpoint_interval:,}")
    print("\nSKIPPING: RF models (analysis_3), MI matrix (analysis_5)")
    print("=" * 70)
    
    # Find input files
    import glob
    if '*' in args.npz:
        npz_files = sorted(glob.glob(args.npz))
    else:
        npz_files = [args.npz]
    
    if not npz_files:
        print(f"ERROR: No files found matching {args.npz}")
        sys.exit(1)
    
    print(f"\nFound {len(npz_files)} NPZ file(s)")
    
    # Initialize processor
    processor = LightEpistasisProcessor(args.output, args.checkpoint_interval)
    
    total_processed = 0
    max_seqs = args.max_sequences or float('inf')
    
    for npz_path in npz_files:
        if total_processed >= max_seqs:
            break
        
        print(f"\nProcessing: {npz_path}")
        
        try:
            records, metadata = load_npz_shard(npz_path)
            print(f"  Loaded {len(records):,} sequences")
            
            valid_count = 0
            for record in tqdm(records, desc="  Processing"):
                if total_processed >= max_seqs:
                    break
                
                # Extract regions and build IMGT positions
                aa_v_full = record.get('aa_v_full', '')
                cdr1 = record.get('cdr1', '')
                cdr2 = record.get('cdr2', '')
                cdr3 = record.get('cdr3', '')
                
                regions = extract_regions_from_sequence(aa_v_full, cdr1, cdr2, cdr3)
                
                if not regions['valid']:
                    continue
                
                processor.process_sequence(
                    regions['positions'], 
                    regions['cdr1'], 
                    regions['cdr2'], 
                    regions['cdr3'], 
                    aa_v_full
                )
                total_processed += 1
                valid_count += 1
            
            print(f"  Valid sequences: {valid_count:,}")
            
            # Free memory after each shard
            del records
            gc.collect()
            
        except Exception as e:
            print(f"  ERROR: {e}")
            import traceback
            traceback.print_exc()
            continue
    
    # Finalize
    print("\n" + "=" * 70)
    results = processor.finalize()
    
    # Save results
    output_path = os.path.join(args.output, 'epistasis_imgt_light.pkl')
    print(f"\nSaving to {output_path}...")
    with open(output_path, 'wb') as f:
        pickle.dump(results, f)
    
    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"Total sequences: {results['total_sequences']:,}")
    print(f"Elapsed time: {results['elapsed_seconds']/60:.1f} minutes")
    print(f"\nFamily distribution:")
    for family, count in sorted(results['family_counts'].items(), key=lambda x: -x[1]):
        pct = 100 * count / results['total_sequences']
        print(f"  {family}: {count:,} ({pct:.1f}%)")
    
    print(f"\nAnalysis 1 (compensation): {len(results['analysis_1_compensation'])} positions")
    print(f"Analysis 2 (vernier): {len(results['analysis_2_vernier_clusters'])} clusters")
    print(f"Analysis 3 (RF models): SKIPPED")
    print(f"Analysis 4 (rules): {len(results['analysis_4_higher_order_rules'])} rules")
    print(f"Analysis 5 (MI): SKIPPED")
    
    print(f"\nOutput saved to: {output_path}")
    print("Done!")


if __name__ == '__main__':
    main()