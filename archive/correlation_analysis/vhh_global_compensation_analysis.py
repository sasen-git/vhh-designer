#!/usr/bin/env python3
"""
VHH Global Compensation Analysis
================================

Following ChatGPT's recommendations for finding "counterbalance" patterns:
1. Feature-level compensation scans (FW position vs CDR global features)
2. Vernier zone clustering + CDR profiling per cluster
3. Mutual information between positions

Input: NPZ files with ANARCI-numbered VHH sequences
Output: Compensation rules, Vernier clusters, MI matrices
"""

import numpy as np
import pickle
import os
import sys
from collections import defaultdict, Counter
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional
import json

# ============================================================
# AMINO ACID PROPERTIES
# ============================================================

# Kyte-Doolittle hydrophobicity scale
HYDROPHOBICITY = {
    'I': 4.5, 'V': 4.2, 'L': 3.8, 'F': 2.8, 'C': 2.5, 'M': 1.9, 'A': 1.8,
    'G': -0.4, 'T': -0.7, 'S': -0.8, 'W': -0.9, 'Y': -1.3, 'P': -1.6,
    'H': -3.2, 'E': -3.5, 'Q': -3.5, 'D': -3.5, 'N': -3.5, 'K': -3.9, 'R': -4.5,
    '-': 0, 'X': 0, '': 0
}

# Charge at pH 7
CHARGE = {
    'K': 1, 'R': 1, 'H': 0.1,  # Basic (H is ~10% charged at pH 7)
    'D': -1, 'E': -1,          # Acidic
    'A': 0, 'C': 0, 'F': 0, 'G': 0, 'I': 0, 'L': 0, 'M': 0,
    'N': 0, 'P': 0, 'Q': 0, 'S': 0, 'T': 0, 'V': 0, 'W': 0, 'Y': 0,
    '-': 0, 'X': 0, '': 0
}

# Aromatic residues
AROMATICS = set('FWY')

# Categories for residues
def residue_category(aa):
    if aa in 'AILMFVW':
        return 'hydrophobic'
    elif aa in 'STNQ':
        return 'polar'
    elif aa in 'DE':
        return 'acidic'
    elif aa in 'KRH':
        return 'basic'
    elif aa in 'GP':
        return 'special'
    elif aa in 'CY':
        return 'reactive'
    else:
        return 'other'

# ============================================================
# CDR FEATURE EXTRACTION
# ============================================================

@dataclass
class CDRFeatures:
    """Global features for a CDR."""
    length: int
    net_charge: float
    hydrophobicity: float
    n_aromatics: int
    n_cys: int
    n_basic: int
    n_acidic: int
    n_glycine: int
    n_proline: int
    
def compute_cdr_features(cdr_seq: str) -> CDRFeatures:
    """Compute global features for a CDR sequence."""
    seq = cdr_seq.replace('-', '').replace('X', '')
    
    if not seq:
        return CDRFeatures(0, 0, 0, 0, 0, 0, 0, 0, 0)
    
    return CDRFeatures(
        length=len(seq),
        net_charge=sum(CHARGE.get(aa, 0) for aa in seq),
        hydrophobicity=sum(HYDROPHOBICITY.get(aa, 0) for aa in seq) / len(seq),
        n_aromatics=sum(1 for aa in seq if aa in AROMATICS),
        n_cys=seq.count('C'),
        n_basic=sum(1 for aa in seq if aa in 'KRH'),
        n_acidic=sum(1 for aa in seq if aa in 'DE'),
        n_glycine=seq.count('G'),
        n_proline=seq.count('P'),
    )

# ============================================================
# IMGT POSITION DEFINITIONS
# ============================================================

# Key Vernier zone positions (IMGT numbering)
# These are FW positions that contact or influence CDR conformation
VERNIER_POSITIONS = {
    'FR1': [2, 4, 25, 26],              # Near CDR1
    'FR2': [36, 39, 40, 41, 42, 44, 46, 47, 48, 49, 50, 51, 52],  # VHH hallmarks + CDR1/2 contacts
    'FR3': [66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 82, 82, 83, 84, 85, 87, 89, 91, 93, 94],  # CDR2/3 contacts
    'FR4': [103, 104, 105, 106, 107, 108],  # Near CDR3
}

# Flatten for easy iteration
ALL_VERNIER = []
for region, positions in VERNIER_POSITIONS.items():
    for pos in positions:
        ALL_VERNIER.append((region, pos))

# CDR positions (IMGT)
CDR_RANGES = {
    'CDR1': list(range(27, 39)),   # 27-38
    'CDR2': list(range(56, 66)),   # 56-65
    'CDR3': list(range(105, 118)), # 105-117 (variable)
}

# ============================================================
# DATA LOADING (from NPZ)
# ============================================================

def load_sequences_from_npz(npz_path: str, max_seqs: int = None) -> Tuple[List[Dict], Dict]:
    """
    Load sequences from NPZ file.
    Returns list of sequence dicts with IMGT positions.
    """
    data = np.load(npz_path, allow_pickle=True)
    
    sequences = []
    imgt_positions = data.get('imgt_positions', None)
    if imgt_positions is not None:
        imgt_positions = imgt_positions.tolist() if hasattr(imgt_positions, 'tolist') else list(imgt_positions)
    
    seq_array = data.get('sequences', data.get('aligned_sequences', None))
    if seq_array is None:
        return [], {}
    
    for i, seq in enumerate(seq_array):
        if max_seqs and i >= max_seqs:
            break
        
        if isinstance(seq, bytes):
            seq = seq.decode('utf-8')
        elif isinstance(seq, np.ndarray):
            seq = ''.join(seq.astype(str))
        
        sequences.append(str(seq))
    
    return sequences, {'imgt_positions': imgt_positions}

# ============================================================
# ANALYSIS 1: FW POSITION vs CDR GLOBAL FEATURES
# ============================================================

class CompensationScanner:
    """
    For each FW position, analyze how CDR global features vary by residue.
    Finds patterns like: "When FR71 is hydrophobic, CDR3 tends to be more polar"
    """
    
    def __init__(self):
        self.position_buckets = defaultdict(lambda: defaultdict(list))
        # position_buckets[fw_position][residue] = list of CDRFeatures
        
        self.family_buckets = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
        # family_buckets[family][fw_position][residue] = list of CDRFeatures
    
    def add_sequence(self, seq_dict: Dict[int, str], family: str):
        """Add one sequence to the analysis."""
        # Extract CDR features
        cdr1 = ''.join(seq_dict.get(pos, '-') for pos in CDR_RANGES['CDR1'])
        cdr2 = ''.join(seq_dict.get(pos, '-') for pos in CDR_RANGES['CDR2'])
        cdr3 = ''.join(seq_dict.get(pos, '-') for pos in CDR_RANGES['CDR3'])
        
        cdr1_feat = compute_cdr_features(cdr1)
        cdr2_feat = compute_cdr_features(cdr2)
        cdr3_feat = compute_cdr_features(cdr3)
        
        combined_features = {
            'cdr1': cdr1_feat,
            'cdr2': cdr2_feat,
            'cdr3': cdr3_feat,
            'total_charge': cdr1_feat.net_charge + cdr2_feat.net_charge + cdr3_feat.net_charge,
            'total_hydro': (cdr1_feat.hydrophobicity + cdr2_feat.hydrophobicity + cdr3_feat.hydrophobicity) / 3,
            'total_cys': cdr1_feat.n_cys + cdr2_feat.n_cys + cdr3_feat.n_cys,
        }
        
        # Bucket by each Vernier position
        for region, pos in ALL_VERNIER:
            residue = seq_dict.get(pos, '-')
            if residue not in '-X':
                self.position_buckets[pos][residue].append(combined_features)
                self.family_buckets[family][pos][residue].append(combined_features)
    
    def compute_compensation_effects(self, min_count: int = 1000) -> Dict:
        """
        For each position, compute how CDR features vary by residue.
        Returns effect sizes and statistical tests.
        """
        results = {}
        
        for pos, residue_buckets in self.position_buckets.items():
            if len(residue_buckets) < 2:
                continue
            
            # Get feature distributions for each residue
            residue_stats = {}
            for residue, features_list in residue_buckets.items():
                if len(features_list) < min_count:
                    continue
                
                # Compute mean and std for each feature
                cdr3_charges = [f['cdr3'].net_charge for f in features_list]
                cdr3_hydro = [f['cdr3'].hydrophobicity for f in features_list]
                cdr3_lengths = [f['cdr3'].length for f in features_list]
                total_charges = [f['total_charge'] for f in features_list]
                total_cys = [f['total_cys'] for f in features_list]
                
                residue_stats[residue] = {
                    'n': len(features_list),
                    'cdr3_charge_mean': np.mean(cdr3_charges),
                    'cdr3_charge_std': np.std(cdr3_charges),
                    'cdr3_hydro_mean': np.mean(cdr3_hydro),
                    'cdr3_hydro_std': np.std(cdr3_hydro),
                    'cdr3_length_mean': np.mean(cdr3_lengths),
                    'cdr3_length_std': np.std(cdr3_lengths),
                    'total_charge_mean': np.mean(total_charges),
                    'total_cys_mean': np.mean(total_cys),
                }
            
            if len(residue_stats) >= 2:
                results[pos] = residue_stats
        
        return results
    
    def find_compensation_rules(self, results: Dict, effect_threshold: float = 0.5) -> List[Dict]:
        """
        Find specific compensation rules.
        Returns rules like: "At FR71, hydrophobic residues (F/I/V) require CDR3 charge > -2"
        """
        rules = []
        
        for pos, residue_stats in results.items():
            # Group residues by category
            categories = defaultdict(list)
            for residue, stats in residue_stats.items():
                cat = residue_category(residue)
                categories[cat].append((residue, stats))
            
            # Compare categories
            for cat1, residues1 in categories.items():
                for cat2, residues2 in categories.items():
                    if cat1 >= cat2:
                        continue
                    
                    # Compare CDR3 charge
                    charges1 = [s['cdr3_charge_mean'] for r, s in residues1]
                    charges2 = [s['cdr3_charge_mean'] for r, s in residues2]
                    
                    if charges1 and charges2:
                        diff = np.mean(charges1) - np.mean(charges2)
                        if abs(diff) > effect_threshold:
                            rules.append({
                                'position': pos,
                                'category1': cat1,
                                'residues1': [r for r, s in residues1],
                                'category2': cat2,
                                'residues2': [r for r, s in residues2],
                                'feature': 'cdr3_charge',
                                'effect': diff,
                                'interpretation': f"At position {pos}, {cat1} residues ({','.join([r for r,s in residues1])}) associate with CDR3 charge {diff:+.2f} vs {cat2}"
                            })
        
        return sorted(rules, key=lambda x: -abs(x['effect']))

# ============================================================
# ANALYSIS 2: VERNIER ZONE CLUSTERING
# ============================================================

class VernierClusterer:
    """
    Cluster sequences by their Vernier zone pattern.
    Profile CDRs for each cluster to find "Vernier archetypes".
    """
    
    def __init__(self, vernier_positions: List[int] = None):
        self.vernier_positions = vernier_positions or [42, 49, 50, 52, 71, 73, 78]
        self.sequences = []
        self.families = []
        self.vernier_patterns = []
        self.cdr_features = []
    
    def add_sequence(self, seq_dict: Dict[int, str], family: str):
        """Add a sequence for clustering."""
        # Extract Vernier pattern
        pattern = tuple(seq_dict.get(pos, '-') for pos in self.vernier_positions)
        
        # Extract CDR features
        cdr1 = ''.join(seq_dict.get(pos, '-') for pos in CDR_RANGES['CDR1'])
        cdr2 = ''.join(seq_dict.get(pos, '-') for pos in CDR_RANGES['CDR2'])
        cdr3 = ''.join(seq_dict.get(pos, '-') for pos in CDR_RANGES['CDR3'])
        
        features = {
            'cdr1': compute_cdr_features(cdr1),
            'cdr2': compute_cdr_features(cdr2),
            'cdr3': compute_cdr_features(cdr3),
        }
        
        self.vernier_patterns.append(pattern)
        self.families.append(family)
        self.cdr_features.append(features)
    
    def cluster_by_pattern(self, min_cluster_size: int = 1000) -> Dict:
        """
        Group sequences by exact Vernier pattern.
        Returns pattern -> CDR feature statistics.
        """
        pattern_groups = defaultdict(list)
        
        for pattern, features, family in zip(self.vernier_patterns, self.cdr_features, self.families):
            pattern_groups[pattern].append((features, family))
        
        # Compute stats for each pattern
        cluster_stats = {}
        for pattern, members in pattern_groups.items():
            if len(members) < min_cluster_size:
                continue
            
            features_list = [m[0] for m in members]
            families = [m[1] for m in members]
            
            # CDR3 statistics
            cdr3_lengths = [f['cdr3'].length for f in features_list]
            cdr3_charges = [f['cdr3'].net_charge for f in features_list]
            cdr3_hydro = [f['cdr3'].hydrophobicity for f in features_list]
            cdr3_cys = [f['cdr3'].n_cys for f in features_list]
            
            cluster_stats[pattern] = {
                'n': len(members),
                'family_dist': dict(Counter(families)),
                'vernier_residues': dict(zip(self.vernier_positions, pattern)),
                'cdr3': {
                    'length_mean': np.mean(cdr3_lengths),
                    'length_std': np.std(cdr3_lengths),
                    'charge_mean': np.mean(cdr3_charges),
                    'charge_std': np.std(cdr3_charges),
                    'hydro_mean': np.mean(cdr3_hydro),
                    'hydro_std': np.std(cdr3_hydro),
                    'cys_mean': np.mean(cdr3_cys),
                },
            }
        
        return cluster_stats
    
    def find_archetype_constraints(self, cluster_stats: Dict) -> List[Dict]:
        """
        For each Vernier archetype, define the CDR "envelope".
        Returns constraints like: "Vernier pattern FERGL only works with CDR3 length 14-20"
        """
        constraints = []
        
        for pattern, stats in cluster_stats.items():
            # Define envelope as mean ± 2*std
            cdr3 = stats['cdr3']
            
            constraint = {
                'vernier_pattern': stats['vernier_residues'],
                'n_sequences': stats['n'],
                'families': stats['family_dist'],
                'cdr3_length_range': (
                    max(5, cdr3['length_mean'] - 2 * cdr3['length_std']),
                    cdr3['length_mean'] + 2 * cdr3['length_std']
                ),
                'cdr3_charge_range': (
                    cdr3['charge_mean'] - 2 * cdr3['charge_std'],
                    cdr3['charge_mean'] + 2 * cdr3['charge_std']
                ),
                'cdr3_hydro_range': (
                    cdr3['hydro_mean'] - 2 * cdr3['hydro_std'],
                    cdr3['hydro_mean'] + 2 * cdr3['hydro_std']
                ),
            }
            constraints.append(constraint)
        
        return sorted(constraints, key=lambda x: -x['n_sequences'])

# ============================================================
# ANALYSIS 3: MUTUAL INFORMATION
# ============================================================

def compute_mutual_information(pos1_residues: List[str], pos2_residues: List[str]) -> float:
    """
    Compute mutual information between two position columns.
    Higher MI = stronger epistatic coupling.
    """
    # Joint distribution
    joint = Counter(zip(pos1_residues, pos2_residues))
    n = len(pos1_residues)
    
    # Marginal distributions
    p1 = Counter(pos1_residues)
    p2 = Counter(pos2_residues)
    
    mi = 0.0
    for (a, b), count in joint.items():
        if count == 0:
            continue
        p_ab = count / n
        p_a = p1[a] / n
        p_b = p2[b] / n
        if p_a > 0 and p_b > 0:
            mi += p_ab * np.log2(p_ab / (p_a * p_b))
    
    return mi

class MutualInformationAnalyzer:
    """
    Compute MI between all position pairs to find epistatic couplings.
    Focus on FR-CDR pairs.
    """
    
    def __init__(self, positions: List[int]):
        self.positions = positions
        self.columns = {pos: [] for pos in positions}
    
    def add_sequence(self, seq_dict: Dict[int, str]):
        """Add a sequence."""
        for pos in self.positions:
            self.columns[pos].append(seq_dict.get(pos, '-'))
    
    def compute_mi_matrix(self) -> np.ndarray:
        """Compute pairwise MI for all positions."""
        n_pos = len(self.positions)
        mi_matrix = np.zeros((n_pos, n_pos))
        
        for i, pos1 in enumerate(self.positions):
            for j, pos2 in enumerate(self.positions):
                if i < j:
                    mi = compute_mutual_information(self.columns[pos1], self.columns[pos2])
                    mi_matrix[i, j] = mi
                    mi_matrix[j, i] = mi
        
        return mi_matrix
    
    def get_top_couplings(self, mi_matrix: np.ndarray, n_top: int = 50) -> List[Dict]:
        """Get top MI couplings."""
        couplings = []
        n_pos = len(self.positions)
        
        for i in range(n_pos):
            for j in range(i + 1, n_pos):
                couplings.append({
                    'pos1': self.positions[i],
                    'pos2': self.positions[j],
                    'mi': mi_matrix[i, j],
                })
        
        return sorted(couplings, key=lambda x: -x['mi'])[:n_top]

# ============================================================
# MAIN PROCESSING PIPELINE
# ============================================================

def classify_family(seq_dict: Dict[int, str]) -> str:
    """Classify into germline family based on FR2 hallmarks."""
    pos42 = seq_dict.get(42, '-')
    pos49 = seq_dict.get(49, '-')
    pos50 = seq_dict.get(50, '-')
    pos52 = seq_dict.get(52, '-')
    
    if pos50 == 'L':
        return 'VH_like'
    elif pos52 == 'W':
        return 'VHH_W52'
    elif pos42 == 'Y':
        cdr3 = ''.join(seq_dict.get(p, '') for p in CDR_RANGES['CDR3'])
        n_cys = cdr3.count('C')
        return 'Y_C4' if n_cys >= 2 else 'Y_C2'
    elif pos42 == 'F':
        cdr3 = ''.join(seq_dict.get(p, '') for p in CDR_RANGES['CDR3'])
        n_cys = cdr3.count('C')
        return 'F_C4' if n_cys >= 2 else 'F_C2'
    else:
        return 'Other'

def process_npz_file(npz_path: str, scanner: CompensationScanner, 
                     clusterer: VernierClusterer, mi_analyzer: MutualInformationAnalyzer,
                     imgt_positions: List[int] = None) -> int:
    """Process one NPZ file, adding sequences to all analyzers."""
    
    try:
        data = np.load(npz_path, allow_pickle=True)
    except Exception as e:
        print(f"  Error loading {npz_path}: {e}")
        return 0
    
    # Get IMGT positions
    if imgt_positions is None:
        imgt_positions = data.get('imgt_positions', None)
        if imgt_positions is not None:
            imgt_positions = list(imgt_positions)
    
    if imgt_positions is None:
        print(f"  No IMGT positions in {npz_path}")
        return 0
    
    # Get sequences
    seq_array = data.get('sequences', data.get('aligned_sequences', None))
    if seq_array is None:
        return 0
    
    count = 0
    for seq in seq_array:
        if isinstance(seq, bytes):
            seq = seq.decode('utf-8')
        elif isinstance(seq, np.ndarray):
            seq = ''.join(seq.astype(str))
        
        # Build position dict
        seq_dict = {pos: seq[i] for i, pos in enumerate(imgt_positions) if i < len(seq)}
        
        # Classify family
        family = classify_family(seq_dict)
        
        # Add to analyzers
        scanner.add_sequence(seq_dict, family)
        clusterer.add_sequence(seq_dict, family)
        mi_analyzer.add_sequence(seq_dict)
        
        count += 1
    
    return count

def main():
    import argparse
    
    parser = argparse.ArgumentParser(description='VHH Global Compensation Analysis')
    parser.add_argument('--npz-dir', '-d', required=True, help='Directory with NPZ files')
    parser.add_argument('--output', '-o', required=True, help='Output pickle file')
    parser.add_argument('--max-files', type=int, default=None, help='Max NPZ files to process')
    parser.add_argument('--sample-rate', type=float, default=0.1, help='Sample rate (0.1 = 10%)')
    
    args = parser.parse_args()
    
    # Initialize analyzers
    print("Initializing analyzers...")
    
    # Key positions for analysis
    key_vernier = [42, 49, 50, 52, 71, 73, 78, 82, 84, 94]
    key_cdr = list(range(27, 39)) + list(range(56, 66)) + list(range(105, 118))
    all_positions = sorted(set(key_vernier + key_cdr))
    
    scanner = CompensationScanner()
    clusterer = VernierClusterer(vernier_positions=key_vernier)
    mi_analyzer = MutualInformationAnalyzer(positions=all_positions)
    
    # Find NPZ files
    npz_files = sorted([f for f in os.listdir(args.npz_dir) if f.endswith('.npz')])
    if args.max_files:
        npz_files = npz_files[:args.max_files]
    
    print(f"Found {len(npz_files)} NPZ files")
    
    # Process files
    total_seqs = 0
    for i, npz_file in enumerate(npz_files):
        npz_path = os.path.join(args.npz_dir, npz_file)
        
        # Sample
        if args.sample_rate < 1.0 and np.random.random() > args.sample_rate:
            continue
        
        count = process_npz_file(npz_path, scanner, clusterer, mi_analyzer)
        total_seqs += count
        
        if (i + 1) % 100 == 0:
            print(f"  Processed {i+1}/{len(npz_files)} files, {total_seqs:,} sequences")
    
    print(f"\nTotal sequences: {total_seqs:,}")
    
    # Compute results
    print("\nComputing compensation effects...")
    compensation_results = scanner.compute_compensation_effects(min_count=100)
    compensation_rules = scanner.find_compensation_rules(compensation_results)
    
    print(f"  Found {len(compensation_rules)} compensation rules")
    
    print("\nClustering Vernier patterns...")
    cluster_stats = clusterer.cluster_by_pattern(min_cluster_size=100)
    archetype_constraints = clusterer.find_archetype_constraints(cluster_stats)
    
    print(f"  Found {len(cluster_stats)} Vernier archetypes")
    
    print("\nComputing mutual information...")
    mi_matrix = mi_analyzer.compute_mi_matrix()
    top_couplings = mi_analyzer.get_top_couplings(mi_matrix, n_top=100)
    
    print(f"  Top coupling: positions {top_couplings[0]['pos1']}-{top_couplings[0]['pos2']} (MI={top_couplings[0]['mi']:.3f})")
    
    # Save results
    results = {
        'total_sequences': total_seqs,
        'compensation': {
            'position_stats': compensation_results,
            'rules': compensation_rules,
        },
        'vernier_clusters': {
            'stats': {str(k): v for k, v in cluster_stats.items()},  # Tuple keys -> str
            'constraints': archetype_constraints,
        },
        'mutual_information': {
            'positions': mi_analyzer.positions,
            'matrix': mi_matrix.tolist(),
            'top_couplings': top_couplings,
        },
    }
    
    with open(args.output, 'wb') as f:
        pickle.dump(results, f)
    
    print(f"\nSaved results to {args.output}")
    
    # Print summary
    print("\n" + "="*80)
    print("TOP COMPENSATION RULES")
    print("="*80)
    for rule in compensation_rules[:10]:
        print(f"  {rule['interpretation']}")
    
    print("\n" + "="*80)
    print("TOP VERNIER ARCHETYPES")
    print("="*80)
    for constraint in archetype_constraints[:5]:
        v = constraint['vernier_pattern']
        print(f"  Pattern: {v}")
        print(f"    N = {constraint['n_sequences']:,}")
        print(f"    CDR3 length: {constraint['cdr3_length_range'][0]:.0f}-{constraint['cdr3_length_range'][1]:.0f}")
        print(f"    CDR3 charge: {constraint['cdr3_charge_range'][0]:.1f} to {constraint['cdr3_charge_range'][1]:.1f}")
    
    print("\n" + "="*80)
    print("TOP EPISTATIC COUPLINGS (MI)")
    print("="*80)
    for c in top_couplings[:15]:
        pos1, pos2 = c['pos1'], c['pos2']
        region1 = 'CDR' if pos1 in key_cdr else 'FR'
        region2 = 'CDR' if pos2 in key_cdr else 'FR'
        print(f"  {region1}-{pos1} <-> {region2}-{pos2}: MI = {c['mi']:.4f}")

if __name__ == '__main__':
    main()
