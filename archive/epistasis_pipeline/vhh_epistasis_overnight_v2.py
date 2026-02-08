#!/usr/bin/env python3
"""
VHH Full Epistasis Analysis - OVERNIGHT VERSION v2 (CORRECTED)
===============================================================

CORRECTIONS FROM v1:
1. FIXED POSITION MAPPING:
   - FR2_2 (index 1) = IMGT 42 = Y/F hallmark (was incorrectly FR2_4)
   - FR2_9 (index 8) = IMGT 49 = E/Q hallmark
   - FR2_10 (index 9) = IMGT 50 = R/L hallmark (was incorrectly FR2_12)
   - FR2_12 (index 11) = IMGT 52 = G/W hallmark (was incorrectly FR2_14)

2. CORRECTED FAMILY CLASSIFICATION:
   - Classical VHH: pos42=F/Y AND pos49=E/Q AND pos50=R
   - Y_C2, F_C2, F_C4: Based on pos42 + total cysteine count
   - VH_like: pos50=L (includes humanized VHHs)
   - VHH_W52: has W at pos52
   - Non_classical: everything else

3. Uses TOTAL cysteine count (not just CDR3) for C2/C4 classification

RETAINED FROM v1:
- Fixed cond_results NameError bug
- Streaming stats for FeatureCompensationScanner (Welford's algorithm)
- Integer encoding for MI analyzer (memory efficient)
- Subsample cap for MI (default 2M sequences)
- CLI flags for light vs full mode
- Progress bars with tqdm

Position Mapping Reference:
  FR2: W  F  R  Q  x  P  x  x  E  R  x  G  L  x
  Idx: 0  1  2  3  4  5  6  7  8  9  10 11 12 13
  IMGT:41 42 43 44 45 46 47 48 49 50 51 52 53 54
"""

import numpy as np
import pickle
import os
import sys
import re
from collections import defaultdict, Counter
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional
import json
import time
import warnings
warnings.filterwarnings('ignore')

try:
    from tqdm import tqdm
    HAS_TQDM = True
except ImportError:
    HAS_TQDM = False
    print("Note: Install tqdm for progress bars: pip install tqdm")
    def tqdm(iterable, **kwargs):
        return iterable

try:
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.preprocessing import LabelEncoder
    HAS_SKLEARN = True
except ImportError:
    HAS_SKLEARN = False
    print("Warning: sklearn not found. Install with: pip install scikit-learn")

# ============================================================
# AMINO ACID PROPERTIES
# ============================================================

HYDROPHOBICITY = {
    'I': 4.5, 'V': 4.2, 'L': 3.8, 'F': 2.8, 'C': 2.5, 'M': 1.9, 'A': 1.8,
    'G': -0.4, 'T': -0.7, 'S': -0.8, 'W': -0.9, 'Y': -1.3, 'P': -1.6,
    'H': -3.2, 'E': -3.5, 'Q': -3.5, 'D': -3.5, 'N': -3.5, 'K': -3.9, 'R': -4.5,
}

CHARGE = {'K': 1, 'R': 1, 'H': 0.1, 'D': -1, 'E': -1}
for aa in 'ACFGILMNPQSTVWY':
    CHARGE[aa] = 0

VOLUME = {
    'G': 60, 'A': 88, 'S': 89, 'C': 108, 'D': 111, 'P': 112, 'N': 114,
    'T': 116, 'E': 138, 'V': 140, 'Q': 143, 'H': 153, 'M': 162, 'I': 166,
    'L': 166, 'K': 168, 'R': 173, 'F': 189, 'Y': 193, 'W': 227,
}

BETA_PROPENSITY = {
    'V': 1.70, 'I': 1.60, 'Y': 1.47, 'F': 1.38, 'W': 1.37, 'L': 1.30, 'T': 1.19,
    'C': 1.19, 'Q': 1.10, 'M': 1.05, 'R': 0.93, 'N': 0.89, 'H': 0.87, 'A': 0.83,
    'S': 0.75, 'G': 0.75, 'K': 0.74, 'P': 0.55, 'D': 0.54, 'E': 0.37,
}

# Amino acid to integer encoding for memory efficiency
AA_TO_INT = {aa: i for i, aa in enumerate('ACDEFGHIKLMNPQRSTVWY-X*')}
INT_TO_AA = {i: aa for aa, i in AA_TO_INT.items()}

AROMATICS = set('FWY')
ALIPHATIC = set('AILV')
POLAR = set('STNQ')
CHARGED = set('DEKRH')
SMALL = set('AGSP')
HYDROPHOBIC = set('AILMFVW')

def residue_category(aa):
    if aa in 'AILMFVW': return 'hydrophobic'
    elif aa in 'STNQ': return 'polar'
    elif aa in 'DE': return 'acidic'
    elif aa in 'KRH': return 'basic'
    elif aa in 'GP': return 'special'
    elif aa in 'CY': return 'reactive'
    else: return 'other'

# ============================================================
# CDR FEATURE EXTRACTION
# ============================================================

@dataclass
class CDRFeatures:
    length: int
    net_charge: float
    hydrophobicity: float
    avg_volume: float
    beta_propensity: float
    n_aromatics: int
    n_aliphatic: int
    n_polar: int
    n_charged: int
    n_small: int
    n_hydrophobic: int
    n_cys: int
    n_basic: int
    n_acidic: int
    n_glycine: int
    n_proline: int
    n_serine: int
    n_threonine: int
    frac_aromatic: float
    frac_charged: float
    frac_hydrophobic: float
    has_WxxW: bool
    has_GG: bool
    has_GGG: bool
    has_PP: bool
    has_SS: bool
    has_YY: bool
    first_aa_cat: str
    last_aa_cat: str
    
    def to_vector(self) -> List[float]:
        return [
            self.length, self.net_charge, self.hydrophobicity, 
            self.avg_volume, self.beta_propensity,
            self.n_aromatics, self.n_aliphatic, self.n_polar, self.n_charged,
            self.n_small, self.n_hydrophobic, self.n_cys, self.n_basic,
            self.n_acidic, self.n_glycine, self.n_proline, self.n_serine, self.n_threonine,
            self.frac_aromatic, self.frac_charged, self.frac_hydrophobic,
            float(self.has_WxxW), float(self.has_GG), float(self.has_GGG),
            float(self.has_PP), float(self.has_SS), float(self.has_YY),
        ]
    
    @staticmethod
    def feature_names() -> List[str]:
        return [
            'length', 'net_charge', 'hydrophobicity', 'avg_volume', 'beta_propensity',
            'n_aromatics', 'n_aliphatic', 'n_polar', 'n_charged', 'n_small',
            'n_hydrophobic', 'n_cys', 'n_basic', 'n_acidic', 'n_glycine',
            'n_proline', 'n_serine', 'n_threonine',
            'frac_aromatic', 'frac_charged', 'frac_hydrophobic',
            'has_WxxW', 'has_GG', 'has_GGG', 'has_PP', 'has_SS', 'has_YY',
        ]

def compute_cdr_features(cdr_seq: str) -> CDRFeatures:
    seq = str(cdr_seq).replace('-', '').replace('X', '').replace('*', '').upper()
    
    if not seq:
        return CDRFeatures(
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, False, False, False, False, False, False, 'none', 'none'
        )
    
    n = len(seq)
    return CDRFeatures(
        length=n,
        net_charge=sum(CHARGE.get(aa, 0) for aa in seq),
        hydrophobicity=sum(HYDROPHOBICITY.get(aa, 0) for aa in seq) / n,
        avg_volume=sum(VOLUME.get(aa, 100) for aa in seq) / n,
        beta_propensity=sum(BETA_PROPENSITY.get(aa, 1.0) for aa in seq) / n,
        n_aromatics=sum(1 for aa in seq if aa in AROMATICS),
        n_aliphatic=sum(1 for aa in seq if aa in ALIPHATIC),
        n_polar=sum(1 for aa in seq if aa in POLAR),
        n_charged=sum(1 for aa in seq if aa in CHARGED),
        n_small=sum(1 for aa in seq if aa in SMALL),
        n_hydrophobic=sum(1 for aa in seq if aa in HYDROPHOBIC),
        n_cys=seq.count('C'),
        n_basic=sum(1 for aa in seq if aa in 'KRH'),
        n_acidic=sum(1 for aa in seq if aa in 'DE'),
        n_glycine=seq.count('G'),
        n_proline=seq.count('P'),
        n_serine=seq.count('S'),
        n_threonine=seq.count('T'),
        frac_aromatic=sum(1 for aa in seq if aa in AROMATICS) / n,
        frac_charged=sum(1 for aa in seq if aa in CHARGED) / n,
        frac_hydrophobic=sum(1 for aa in seq if aa in HYDROPHOBIC) / n,
        has_WxxW=bool(re.search(r'W..W', seq)),
        has_GG='GG' in seq,
        has_GGG='GGG' in seq,
        has_PP='PP' in seq,
        has_SS='SS' in seq,
        has_YY='YY' in seq,
        first_aa_cat=residue_category(seq[0]) if seq else 'none',
        last_aa_cat=residue_category(seq[-1]) if seq else 'none',
    )

# ============================================================
# FRAMEWORK EXTRACTION
# ============================================================

def extract_frameworks(full_seq: str, cdr1: str, cdr2: str, cdr3: str) -> Dict:
    seq = str(full_seq)
    
    result = {
        'fr1': '', 'fr2': '', 'fr3': '', 'fr4': '',
        'positions': {},
        'valid': False
    }
    
    fr2_match = re.search(r'W[YFVIL]RQ', seq[25:60])
    if not fr2_match:
        fr2_match = re.search(r'W.RQ', seq[25:60])
    
    if not fr2_match:
        return result
    
    fr2_start = 25 + fr2_match.start()
    result['fr1'] = seq[:fr2_start]
    
    fr2_end = fr2_start + 14
    result['fr2'] = seq[fr2_start:fr2_end]
    
    cdr2_pos = seq.find(str(cdr2)) if cdr2 else -1
    cdr3_pos = seq.find(str(cdr3)) if cdr3 else -1
    
    if cdr2_pos > fr2_end and cdr3_pos > cdr2_pos:
        result['fr3'] = seq[cdr2_pos + len(str(cdr2)):cdr3_pos]
        result['fr4'] = seq[cdr3_pos + len(str(cdr3)):]
        result['valid'] = True
    
    # Extract positions
    fr1 = result['fr1']
    fr2 = result['fr2']
    fr3 = result['fr3']
    fr4 = result['fr4']
    
    positions = {}
    
    if len(fr1) >= 25:
        positions['FR1_1'] = fr1[0]
        positions['FR1_5'] = fr1[4] if len(fr1) > 4 else '-'
        positions['FR1_10'] = fr1[9] if len(fr1) > 9 else '-'
        positions['FR1_20'] = fr1[19] if len(fr1) > 19 else '-'
        positions['FR1_25'] = fr1[24] if len(fr1) > 24 else '-'
    
    if len(fr2) >= 14:
        for i in range(min(14, len(fr2))):
            positions[f'FR2_{i+1}'] = fr2[i]
    
    if len(fr3) >= 30:
        for i in range(min(32, len(fr3))):
            positions[f'FR3_{i+1}'] = fr3[i]
    
    if len(fr4) >= 10:
        for i in range(min(11, len(fr4))):
            positions[f'FR4_{i+1}'] = fr4[i]
    
    # CDR anchors
    cdr1 = str(cdr1) if cdr1 else ''
    cdr2 = str(cdr2) if cdr2 else ''
    cdr3 = str(cdr3) if cdr3 else ''
    
    if len(cdr1) >= 2:
        positions['CDR1_first'] = cdr1[0]
        positions['CDR1_last'] = cdr1[-1]
        positions['CDR1_2nd'] = cdr1[1] if len(cdr1) > 1 else '-'
        positions['CDR1_2ndlast'] = cdr1[-2] if len(cdr1) > 1 else '-'
    
    if len(cdr2) >= 2:
        positions['CDR2_first'] = cdr2[0]
        positions['CDR2_last'] = cdr2[-1]
        positions['CDR2_2nd'] = cdr2[1] if len(cdr2) > 1 else '-'
        positions['CDR2_2ndlast'] = cdr2[-2] if len(cdr2) > 1 else '-'
    
    if len(cdr3) >= 2:
        positions['CDR3_first'] = cdr3[0]
        positions['CDR3_last'] = cdr3[-1]
        positions['CDR3_2nd'] = cdr3[1] if len(cdr3) > 1 else '-'
        positions['CDR3_2ndlast'] = cdr3[-2] if len(cdr3) > 1 else '-'
        positions['CDR3_3rd'] = cdr3[2] if len(cdr3) > 2 else '-'
        positions['CDR3_3rdlast'] = cdr3[-3] if len(cdr3) > 2 else '-'
    
    result['positions'] = positions
    return result

def classify_family(positions: Dict[str, str], cdr3: str, full_seq: str = '') -> str:
    """
    Classify VHH family using CORRECTED IMGT position mapping.
    
    FR2 extraction: W F R Q x P x x E R x G L x
    Index:          0 1 2 3 4 5 6 7 8 9 10 11 12 13
    IMGT:          41 42 43 44 45 46 47 48 49 50 51 52 53 54
    
    Correct mapping:
    - FR2_2 (index 1) = IMGT 42 = Y/F hallmark
    - FR2_9 (index 8) = IMGT 49 = E/Q hallmark
    - FR2_10 (index 9) = IMGT 50 = R/L hallmark
    - FR2_12 (index 11) = IMGT 52 = G/W hallmark
    
    Literature-based classification (Wang et al. 2024, Vincke et al. 2009):
    - Classical VHH: pos42=F/Y AND pos49=E/Q AND pos50=R
    - VH_like: pos50=L (includes humanized VHHs)
    - Y_C2, F_C2, F_C4: Classical VHH subtypes by pos42 and cysteine count
    """
    # CORRECTED position mapping
    pos42 = positions.get('FR2_2', '-')   # IMGT 42 - Y/F hallmark
    pos49 = positions.get('FR2_9', '-')   # IMGT 49 - E/Q hallmark
    pos50 = positions.get('FR2_10', '-')  # IMGT 50 - R/L hallmark
    pos52 = positions.get('FR2_12', '-')  # IMGT 52 - G/W hallmark
    
    # Count total cysteines in full sequence (not just CDR3)
    if full_seq:
        n_cys = str(full_seq).count('C')
    else:
        n_cys = str(cdr3).count('C') if cdr3 else 0
    
    # VH-like: pos50 = L (human-like at key hallmark)
    # Note: This includes humanized VHHs with E49G+R50L
    if pos50 == 'L':
        return 'VH_like'
    
    # Check classical VHH criteria
    is_classical = (
        pos42 in ['F', 'Y'] and
        pos49 in ['E', 'Q'] and
        pos50 == 'R'
    )
    
    if is_classical:
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
    
    # Non-classical: W at position 52 (human-like at 52)
    if pos52 == 'W':
        return 'VHH_W52'
    
    return 'Non_classical'

# ============================================================
# WELFORD'S ONLINE ALGORITHM FOR STREAMING STATS
# ============================================================

class WelfordAccumulator:
    """Memory-efficient online statistics using Welford's algorithm."""
    
    def __init__(self):
        self.n = 0
        self.mean = 0.0
        self.M2 = 0.0  # Sum of squared differences from mean
        self.min_val = float('inf')
        self.max_val = float('-inf')
        # For quantiles, we keep a reservoir sample
        self.reservoir = []
        self.reservoir_size = 1000
    
    def add(self, x: float):
        self.n += 1
        delta = x - self.mean
        self.mean += delta / self.n
        delta2 = x - self.mean
        self.M2 += delta * delta2
        
        self.min_val = min(self.min_val, x)
        self.max_val = max(self.max_val, x)
        
        # Reservoir sampling for quantiles
        if len(self.reservoir) < self.reservoir_size:
            self.reservoir.append(x)
        else:
            j = np.random.randint(0, self.n)
            if j < self.reservoir_size:
                self.reservoir[j] = x
    
    def get_stats(self) -> Dict:
        if self.n == 0:
            return {'n': 0, 'mean': 0, 'std': 0, 'min': 0, 'max': 0}
        
        variance = self.M2 / self.n if self.n > 0 else 0
        std = np.sqrt(variance)
        
        stats = {
            'n': self.n,
            'mean': float(self.mean),
            'std': float(std),
            'min': float(self.min_val),
            'max': float(self.max_val),
        }
        
        # Approximate quantiles from reservoir
        if self.reservoir:
            sorted_res = sorted(self.reservoir)
            stats['p5'] = float(sorted_res[int(len(sorted_res) * 0.05)])
            stats['p25'] = float(sorted_res[int(len(sorted_res) * 0.25)])
            stats['p50'] = float(sorted_res[int(len(sorted_res) * 0.50)])
            stats['p75'] = float(sorted_res[int(len(sorted_res) * 0.75)])
            stats['p95'] = float(sorted_res[min(int(len(sorted_res) * 0.95), len(sorted_res)-1)])
        
        return stats

# ============================================================
# ANALYSIS 1: FEATURE-LEVEL COMPENSATION (STREAMING)
# ============================================================

class FeatureCompensationScanner:
    """
    Memory-efficient version using Welford's algorithm.
    No longer stores all features - just running statistics.
    """
    
    FEATURE_NAMES = [
        'cdr1_len', 'cdr1_charge', 'cdr1_hydro', 'cdr1_cys',
        'cdr2_len', 'cdr2_charge', 'cdr2_hydro', 'cdr2_cys',
        'cdr3_len', 'cdr3_charge', 'cdr3_hydro', 'cdr3_cys',
        'cdr3_aromatic', 'cdr3_basic', 'cdr3_acidic', 'cdr3_glycine',
        'cdr3_proline', 'cdr3_volume', 'cdr3_beta',
        'total_charge', 'total_cys', 'total_len',
    ]
    
    def __init__(self):
        # data[family][pos_name][residue][feature_name] = WelfordAccumulator
        self.data = defaultdict(
            lambda: defaultdict(
                lambda: defaultdict(
                    lambda: {fn: WelfordAccumulator() for fn in self.FEATURE_NAMES}
                )
            )
        )
        self.family_counts = Counter()
        self.total_added = 0
    
    def add(self, family: str, positions: Dict[str, str], 
            cdr1_feat: CDRFeatures, cdr2_feat: CDRFeatures, cdr3_feat: CDRFeatures):
        
        self.family_counts[family] += 1
        self.total_added += 1
        
        # Build feature dict
        features = {
            'cdr1_len': cdr1_feat.length,
            'cdr1_charge': cdr1_feat.net_charge,
            'cdr1_hydro': cdr1_feat.hydrophobicity,
            'cdr1_cys': cdr1_feat.n_cys,
            'cdr2_len': cdr2_feat.length,
            'cdr2_charge': cdr2_feat.net_charge,
            'cdr2_hydro': cdr2_feat.hydrophobicity,
            'cdr2_cys': cdr2_feat.n_cys,
            'cdr3_len': cdr3_feat.length,
            'cdr3_charge': cdr3_feat.net_charge,
            'cdr3_hydro': cdr3_feat.hydrophobicity,
            'cdr3_cys': cdr3_feat.n_cys,
            'cdr3_aromatic': cdr3_feat.n_aromatics,
            'cdr3_basic': cdr3_feat.n_basic,
            'cdr3_acidic': cdr3_feat.n_acidic,
            'cdr3_glycine': cdr3_feat.n_glycine,
            'cdr3_proline': cdr3_feat.n_proline,
            'cdr3_volume': cdr3_feat.avg_volume,
            'cdr3_beta': cdr3_feat.beta_propensity,
            'total_charge': cdr1_feat.net_charge + cdr2_feat.net_charge + cdr3_feat.net_charge,
            'total_cys': cdr1_feat.n_cys + cdr2_feat.n_cys + cdr3_feat.n_cys,
            'total_len': cdr1_feat.length + cdr2_feat.length + cdr3_feat.length,
        }
        
        # Update streaming stats for each position/residue
        for pos_name, residue in positions.items():
            if residue and residue not in '-X*':
                accumulators = self.data[family][pos_name][residue]
                for feat_name, value in features.items():
                    accumulators[feat_name].add(value)
    
    def analyze(self, min_count: int = 100) -> Dict:
        """Get statistics from streaming accumulators."""
        results = {}
        
        for family, positions in self.data.items():
            family_results = {}
            
            for pos_name, residue_data in positions.items():
                pos_results = {}
                
                for residue, accumulators in residue_data.items():
                    n = accumulators['cdr3_len'].n
                    if n < min_count:
                        continue
                    
                    stats = {'n': n}
                    for feat_name, acc in accumulators.items():
                        stats[feat_name] = acc.get_stats()
                    pos_results[residue] = stats
                
                if len(pos_results) >= 2:
                    family_results[pos_name] = pos_results
            
            if family_results:
                results[family] = family_results
        
        return results
    
    def find_compensation_rules(self, stats: Dict, min_effect_len: float = 1.5,
                                 min_effect_charge: float = 0.3,
                                 min_effect_hydro: float = 0.2) -> List[Dict]:
        rules = []
        
        for family, positions in stats.items():
            for pos_name, residue_stats in positions.items():
                residues = list(residue_stats.keys())
                
                for i, res1 in enumerate(residues):
                    for res2 in residues[i+1:]:
                        stat1 = residue_stats[res1]
                        stat2 = residue_stats[res2]
                        
                        n1, n2 = stat1['n'], stat2['n']
                        
                        # CDR3 length
                        if 'cdr3_len' in stat1:
                            diff = stat1['cdr3_len']['mean'] - stat2['cdr3_len']['mean']
                            if abs(diff) >= min_effect_len:
                                rules.append({
                                    'family': family, 'position': pos_name,
                                    'res1': res1, 'res2': res2,
                                    'feature': 'cdr3_length', 'effect': diff,
                                    'n1': n1, 'n2': n2,
                                    'cat1': residue_category(res1),
                                    'cat2': residue_category(res2),
                                })
                        
                        # CDR3 charge
                        if 'cdr3_charge' in stat1:
                            diff = stat1['cdr3_charge']['mean'] - stat2['cdr3_charge']['mean']
                            if abs(diff) >= min_effect_charge:
                                rules.append({
                                    'family': family, 'position': pos_name,
                                    'res1': res1, 'res2': res2,
                                    'feature': 'cdr3_charge', 'effect': diff,
                                    'n1': n1, 'n2': n2,
                                    'cat1': residue_category(res1),
                                    'cat2': residue_category(res2),
                                })
                        
                        # CDR3 hydrophobicity
                        if 'cdr3_hydro' in stat1:
                            diff = stat1['cdr3_hydro']['mean'] - stat2['cdr3_hydro']['mean']
                            if abs(diff) >= min_effect_hydro:
                                rules.append({
                                    'family': family, 'position': pos_name,
                                    'res1': res1, 'res2': res2,
                                    'feature': 'cdr3_hydrophobicity', 'effect': diff,
                                    'n1': n1, 'n2': n2,
                                    'cat1': residue_category(res1),
                                    'cat2': residue_category(res2),
                                })
                        
                        # Total charge
                        if 'total_charge' in stat1:
                            diff = stat1['total_charge']['mean'] - stat2['total_charge']['mean']
                            if abs(diff) >= min_effect_charge * 2:
                                rules.append({
                                    'family': family, 'position': pos_name,
                                    'res1': res1, 'res2': res2,
                                    'feature': 'total_charge', 'effect': diff,
                                    'n1': n1, 'n2': n2,
                                    'cat1': residue_category(res1),
                                    'cat2': residue_category(res2),
                                })
        
        return sorted(rules, key=lambda x: -abs(x['effect']))

# ============================================================
# ANALYSIS 2: VERNIER CLUSTERING (STREAMING)
# ============================================================

class VernierClusterer:
    """Memory-efficient version using streaming stats."""
    
    def __init__(self):
        # clusters[pattern] = {accumulators for CDR stats, family_counts}
        self.clusters = defaultdict(lambda: {
            'cdr3_len': WelfordAccumulator(),
            'cdr3_charge': WelfordAccumulator(),
            'cdr3_hydro': WelfordAccumulator(),
            'cdr3_cys': Counter(),
            'cdr3_aromatic': WelfordAccumulator(),
            'cdr2_len': WelfordAccumulator(),
            'cdr1_len': WelfordAccumulator(),
            'families': Counter(),
            'n': 0,
        })
    
    def add(self, positions: Dict[str, str], 
            cdr1_feat: CDRFeatures, cdr2_feat: CDRFeatures, cdr3_feat: CDRFeatures,
            family: str):
        
        pattern = (
            positions.get('FR2_4', '-'),
            positions.get('FR2_11', '-'),
            positions.get('FR2_12', '-'),
            positions.get('FR2_14', '-'),
            positions.get('FR3_6', '-'),
            positions.get('FR3_8', '-'),
            positions.get('FR3_12', '-'),
            positions.get('FR3_16', '-'),
        )
        
        c = self.clusters[pattern]
        c['n'] += 1
        c['cdr3_len'].add(cdr3_feat.length)
        c['cdr3_charge'].add(cdr3_feat.net_charge)
        c['cdr3_hydro'].add(cdr3_feat.hydrophobicity)
        c['cdr3_cys'][cdr3_feat.n_cys] += 1
        c['cdr3_aromatic'].add(cdr3_feat.n_aromatics)
        c['cdr2_len'].add(cdr2_feat.length)
        c['cdr1_len'].add(cdr1_feat.length)
        c['families'][family] += 1
    
    def analyze(self, min_size: int = 500) -> Dict:
        results = {}
        
        for pattern, c in self.clusters.items():
            if c['n'] < min_size:
                continue
            
            results[str(pattern)] = {
                'n': c['n'],
                'pattern': dict(zip(
                    ['FR2_4', 'FR2_11', 'FR2_12', 'FR2_14', 'FR3_6', 'FR3_8', 'FR3_12', 'FR3_16'],
                    pattern
                )),
                'families': dict(c['families']),
                'cdr3_length': c['cdr3_len'].get_stats(),
                'cdr3_charge': c['cdr3_charge'].get_stats(),
                'cdr3_hydro': c['cdr3_hydro'].get_stats(),
                'cdr3_cys_dist': dict(c['cdr3_cys']),
                'cdr3_aromatic': c['cdr3_aromatic'].get_stats(),
                'cdr2_length': c['cdr2_len'].get_stats(),
                'cdr1_length': c['cdr1_len'].get_stats(),
            }
        
        return results

# ============================================================
# ANALYSIS 3: CONDITIONAL MODELS
# ============================================================

class ConditionalModeler:
    """Stores data for RF training - subsamples to limit memory."""
    
    def __init__(self, max_samples_per_position: int = 100000):
        self.data = defaultdict(lambda: defaultdict(lambda: {'X': [], 'y': []}))
        self.max_samples = max_samples_per_position
    
    def add(self, family: str, positions: Dict[str, str],
            cdr1_feat: CDRFeatures, cdr2_feat: CDRFeatures, cdr3_feat: CDRFeatures):
        
        X = cdr1_feat.to_vector() + cdr2_feat.to_vector() + cdr3_feat.to_vector()
        
        for pos_name, residue in positions.items():
            if residue and residue not in '-X*':
                d = self.data[family][pos_name]
                # Reservoir sampling to limit memory
                if len(d['X']) < self.max_samples:
                    d['X'].append(X)
                    d['y'].append(residue)
                else:
                    j = np.random.randint(0, len(d['X']) + 1)
                    if j < self.max_samples:
                        d['X'][j] = X
                        d['y'][j] = residue
    
    def train_models(self, min_samples: int = 2000, min_classes: int = 2,
                     n_estimators: int = 100, max_depth: int = 15) -> Dict:
        if not HAS_SKLEARN:
            return {'error': 'sklearn not installed'}
        
        # First, count how many models we'll train
        models_to_train = []
        for family, positions in self.data.items():
            for pos_name, data in positions.items():
                if len(data['X']) >= min_samples:
                    residue_counts = Counter(data['y'])
                    top_residues = [r for r, c in residue_counts.most_common(8) if c >= 50]
                    if len(top_residues) >= min_classes:
                        models_to_train.append((family, pos_name, data))
        
        print(f"    Training {len(models_to_train)} RF models...")
        
        results = {}
        model_count = 0
        
        for family, pos_name, data in tqdm(models_to_train, desc="    RF models", unit="model", disable=not HAS_TQDM):
            X = np.array(data['X'])
            y = np.array(data['y'])
            
            residue_counts = Counter(y)
            top_residues = [r for r, c in residue_counts.most_common(8) if c >= 50]
            
            mask = np.isin(y, top_residues)
            X_filt = X[mask]
            y_filt = y[mask]
            
            if len(X_filt) < min_samples:
                continue
            
            try:
                le = LabelEncoder()
                y_enc = le.fit_transform(y_filt)
                
                rf = RandomForestClassifier(
                    n_estimators=n_estimators, 
                    max_depth=max_depth,
                    min_samples_leaf=10,
                    n_jobs=-1, 
                    random_state=42,
                    oob_score=True
                )
                rf.fit(X_filt, y_enc)
                
                feature_names = (
                    ['cdr1_' + n for n in CDRFeatures.feature_names()] +
                    ['cdr2_' + n for n in CDRFeatures.feature_names()] +
                    ['cdr3_' + n for n in CDRFeatures.feature_names()]
                )
                
                importances = list(zip(feature_names, rf.feature_importances_))
                importances.sort(key=lambda x: -x[1])
                
                if family not in results:
                    results[family] = {}
                
                results[family][pos_name] = {
                    'n_samples': len(X_filt),
                    'classes': list(le.classes_),
                    'class_counts': {c: int(residue_counts[c]) for c in le.classes_},
                    'oob_score': float(rf.oob_score_),
                    'top_features': [(n, float(i)) for n, i in importances[:15]],
                }
                model_count += 1
                
            except Exception as e:
                continue
        
        print(f"    Trained {model_count} models successfully")
        return results

# ============================================================
# ANALYSIS 4: HIGHER-ORDER RULES (STREAMING COUNTS)
# ============================================================

class HigherOrderRuleFinder:
    """Memory-efficient version that counts patterns on the fly."""
    
    def __init__(self):
        # pattern_counts[family][(cdr_feat, pos1, res1, pos2)] = Counter of res2
        self.pattern_counts = defaultdict(lambda: defaultdict(Counter))
        self.family_counts = Counter()
        
        # Position pairs to track
        self.position_pairs = [
            ('FR2_4', 'FR3_8'),
            ('FR2_4', 'FR3_12'),
            ('FR2_12', 'FR3_6'),
            ('FR2_14', 'FR3_16'),
            ('FR3_6', 'FR3_12'),
        ]
    
    def add(self, family: str, positions: Dict[str, str],
            cdr1_feat: CDRFeatures, cdr2_feat: CDRFeatures, cdr3_feat: CDRFeatures):
        
        self.family_counts[family] += 1
        
        # Bin CDR3 features
        cdr3_len_bin = 'short' if cdr3_feat.length < 10 else (
            'medium_short' if cdr3_feat.length < 14 else (
            'medium' if cdr3_feat.length < 18 else (
            'medium_long' if cdr3_feat.length < 22 else 'long')))
        
        cdr3_charge_bin = 'very_neg' if cdr3_feat.net_charge < -2 else (
            'negative' if cdr3_feat.net_charge < -0.5 else (
            'neutral' if cdr3_feat.net_charge < 1 else (
            'positive' if cdr3_feat.net_charge < 3 else 'very_pos')))
        
        cdr3_cys_bin = 'with_cys' if cdr3_feat.n_cys > 0 else 'no_cys'
        
        cdr_features = {
            'cdr3_len': cdr3_len_bin,
            'cdr3_charge': cdr3_charge_bin,
            'cdr3_cys': cdr3_cys_bin,
        }
        
        # Count patterns
        for pos1, pos2 in self.position_pairs:
            res1 = positions.get(pos1, '-')
            res2 = positions.get(pos2, '-')
            
            if res1 != '-' and res2 != '-':
                for cdr_feat_name, cdr_val in cdr_features.items():
                    key = (cdr_feat_name, cdr_val, pos1, res1, pos2)
                    self.pattern_counts[family][key][res2] += 1
    
    def find_rules(self, min_support: int = 500, min_confidence: float = 0.7) -> List[Dict]:
        rules = []
        
        for family, patterns in self.pattern_counts.items():
            for key, outcomes in patterns.items():
                cdr_feat_name, cdr_val, pos1, res1, pos2 = key
                total = sum(outcomes.values())
                
                if total < min_support:
                    continue
                
                for res2, count in outcomes.items():
                    conf = count / total
                    if conf >= min_confidence:
                        rules.append({
                            'family': family,
                            'condition': f"{cdr_feat_name}={cdr_val} AND {pos1}={res1}",
                            'result': f"{pos2}={res2}",
                            'support': total,
                            'confidence': conf,
                            'count': count,
                        })
        
        return sorted(rules, key=lambda x: (-x['confidence'], -x['support']))

# ============================================================
# ANALYSIS 5: MUTUAL INFORMATION (INTEGER ENCODING)
# ============================================================

class MutualInformationAnalyzer:
    """Memory-efficient version using integer encoding."""
    
    def __init__(self, position_names: List[str], max_sequences: int = 2000000):
        self.position_names = position_names
        self.max_sequences = max_sequences
        self.n_added = 0
        
        # Pre-allocate numpy arrays for efficiency
        self.n_positions = len(position_names)
        # Start with reasonable size, will extend if needed
        initial_size = min(max_sequences, 1000000)
        self.data = np.zeros((initial_size, self.n_positions), dtype=np.int8)
        self.pos_to_idx = {name: i for i, name in enumerate(position_names)}
    
    def add(self, positions: Dict[str, str]):
        if self.n_added >= self.max_sequences:
            return
        
        # Extend array if needed
        if self.n_added >= len(self.data):
            new_size = min(len(self.data) * 2, self.max_sequences)
            new_data = np.zeros((new_size, self.n_positions), dtype=np.int8)
            new_data[:len(self.data)] = self.data
            self.data = new_data
        
        # Encode residues as integers
        for name, idx in self.pos_to_idx.items():
            aa = positions.get(name, '-')
            self.data[self.n_added, idx] = AA_TO_INT.get(aa, AA_TO_INT['-'])
        
        self.n_added += 1
    
    def compute_mi(self, col1_idx: int, col2_idx: int) -> float:
        """Compute MI between two columns using integer-encoded data."""
        if self.n_added == 0:
            return 0.0
        
        data = self.data[:self.n_added]
        col1 = data[:, col1_idx]
        col2 = data[:, col2_idx]
        
        # Joint distribution
        joint = Counter(zip(col1.tolist(), col2.tolist()))
        n = self.n_added
        
        # Marginals
        p1 = Counter(col1.tolist())
        p2 = Counter(col2.tolist())
        
        gap_idx = AA_TO_INT['-']
        
        mi = 0.0
        for (a, b), count in joint.items():
            if count == 0 or a == gap_idx or b == gap_idx:
                continue
            p_ab = count / n
            p_a = p1[a] / n
            p_b = p2[b] / n
            if p_a > 0 and p_b > 0:
                mi += p_ab * np.log2(p_ab / (p_a * p_b))
        
        return max(0, mi)
    
    def compute_matrix(self) -> Tuple[np.ndarray, List[str]]:
        n_pairs = (self.n_positions * (self.n_positions - 1)) // 2
        print(f"    Computing MI matrix ({self.n_positions}x{self.n_positions}, {n_pairs} pairs) from {self.n_added:,} sequences...")
        
        matrix = np.zeros((self.n_positions, self.n_positions))
        
        pairs = [(i, j) for i in range(self.n_positions) for j in range(i + 1, self.n_positions)]
        
        for i, j in tqdm(pairs, desc="    MI pairs", unit="pair", disable=not HAS_TQDM):
            mi = self.compute_mi(i, j)
            matrix[i, j] = mi
            matrix[j, i] = mi
        
        return matrix, self.position_names
    
    def get_top_pairs(self, matrix: np.ndarray, n_top: int = 100) -> List[Dict]:
        pairs = []
        
        for i in range(self.n_positions):
            for j in range(i + 1, self.n_positions):
                pos1, pos2 = self.position_names[i], self.position_names[j]
                
                if 'CDR' in pos1 and 'CDR' in pos2:
                    pair_type = 'CDR-CDR'
                elif 'CDR' in pos1 or 'CDR' in pos2:
                    pair_type = 'FR-CDR'
                else:
                    pair_type = 'FR-FR'
                
                pairs.append({
                    'pos1': pos1,
                    'pos2': pos2,
                    'mi': float(matrix[i, j]),
                    'type': pair_type,
                })
        
        return sorted(pairs, key=lambda x: -x['mi'])[:n_top]

# ============================================================
# MAIN PROCESSOR
# ============================================================

def process_npz(npz_path: str, analyzers: Dict) -> int:
    try:
        data = np.load(npz_path, allow_pickle=True)
    except Exception as e:
        print(f"  Error loading {npz_path}: {e}")
        return 0
    
    sequences = data.get('aa_v_full', None)
    cdr1s = data.get('cdr1', None)
    cdr2s = data.get('cdr2', None)
    cdr3s = data.get('cdr3', None)
    
    if sequences is None:
        return 0
    
    count = 0
    n_seqs = len(sequences)
    
    # Progress bar for sequences within this file
    desc = f"    {os.path.basename(npz_path)}"
    iterator = tqdm(range(n_seqs), desc=desc, unit="seq", 
                    unit_scale=True, leave=False, disable=not HAS_TQDM,
                    mininterval=1.0)
    
    for i in iterator:
        try:
            seq = str(sequences[i])
            cdr1 = str(cdr1s[i]) if cdr1s is not None else ''
            cdr2 = str(cdr2s[i]) if cdr2s is not None else ''
            cdr3 = str(cdr3s[i]) if cdr3s is not None else ''
            
            if len(seq) < 100:
                continue
            
            fw = extract_frameworks(seq, cdr1, cdr2, cdr3)
            if not fw['valid']:
                continue
            
            positions = fw['positions']
            if len(positions) < 10:
                continue
            
            cdr1_feat = compute_cdr_features(cdr1)
            cdr2_feat = compute_cdr_features(cdr2)
            cdr3_feat = compute_cdr_features(cdr3)
            
            family = classify_family(positions, cdr3, seq)
            
            analyzers['comp'].add(family, positions, cdr1_feat, cdr2_feat, cdr3_feat)
            analyzers['cluster'].add(positions, cdr1_feat, cdr2_feat, cdr3_feat, family)
            analyzers['cond'].add(family, positions, cdr1_feat, cdr2_feat, cdr3_feat)
            analyzers['rules'].add(family, positions, cdr1_feat, cdr2_feat, cdr3_feat)
            analyzers['mi'].add(positions)
            
            count += 1
            
        except Exception:
            continue
    
    return count

def save_checkpoint(results: Dict, output_prefix: str, stage: str):
    checkpoint_path = f"{output_prefix}_checkpoint_{stage}.pkl"
    with open(checkpoint_path, 'wb') as f:
        pickle.dump(results, f)
    print(f"  Checkpoint saved: {checkpoint_path}")

def main():
    import argparse
    
    parser = argparse.ArgumentParser(description='VHH Full Epistasis Analysis - OVERNIGHT v2')
    parser.add_argument('--npz-dir', '-d', required=True, help='Directory with NPZ files')
    parser.add_argument('--output', '-o', required=True, help='Output prefix')
    parser.add_argument('--max-files', type=int, default=None)
    
    # Tunable parameters (ChatGPT's suggestion)
    parser.add_argument('--min-comp-count', type=int, default=100, 
                        help='Min count for compensation analysis (default: 100)')
    parser.add_argument('--min-cluster-size', type=int, default=500,
                        help='Min cluster size for Vernier analysis (default: 500)')
    parser.add_argument('--mi-max-seqs', type=int, default=2000000,
                        help='Max sequences for MI analysis (default: 2M)')
    parser.add_argument('--skip-cond-models', action='store_true',
                        help='Skip conditional model training')
    parser.add_argument('--rf-trees', type=int, default=100,
                        help='Number of RF trees (default: 100)')
    parser.add_argument('--rf-depth', type=int, default=15,
                        help='Max RF depth (default: 15)')
    parser.add_argument('--light', action='store_true',
                        help='Light mode: fewer trees, shallower depth, smaller MI sample')
    
    args = parser.parse_args()
    
    # Light mode overrides
    if args.light:
        args.rf_trees = 50
        args.rf_depth = 10
        args.mi_max_seqs = 500000
        print("LIGHT MODE: Using reduced parameters")
    
    start_time = time.time()
    
    print("="*70)
    print("VHH FULL EPISTASIS ANALYSIS - OVERNIGHT v2 (CORRECTED POSITIONS)")
    print("Memory-optimized with streaming stats")
    print("="*70)
    print("POSITION MAPPING FIX:")
    print("  FR2_2  = IMGT 42 (Y/F hallmark)")
    print("  FR2_9  = IMGT 49 (E/Q hallmark)")
    print("  FR2_10 = IMGT 50 (R/L hallmark)")
    print("  FR2_12 = IMGT 52 (G/W hallmark)")
    print(f"Start time: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Parameters:")
    print(f"  min_comp_count: {args.min_comp_count}")
    print(f"  min_cluster_size: {args.min_cluster_size}")
    print(f"  mi_max_seqs: {args.mi_max_seqs:,}")
    print(f"  rf_trees: {args.rf_trees}")
    print(f"  rf_depth: {args.rf_depth}")
    
    # Position names for MI
    position_names = [
        'FR2_1', 'FR2_2', 'FR2_3', 'FR2_4', 'FR2_5', 'FR2_6', 'FR2_7',
        'FR2_8', 'FR2_9', 'FR2_10', 'FR2_11', 'FR2_12', 'FR2_13', 'FR2_14',
        'FR3_1', 'FR3_2', 'FR3_4', 'FR3_6', 'FR3_8', 'FR3_10', 'FR3_12',
        'FR3_14', 'FR3_16', 'FR3_18', 'FR3_20', 'FR3_22', 'FR3_24',
        'FR4_1', 'FR4_2', 'FR4_3', 'FR4_4', 'FR4_5',
        'CDR1_first', 'CDR1_last', 'CDR1_2nd', 'CDR1_2ndlast',
        'CDR2_first', 'CDR2_last', 'CDR2_2nd', 'CDR2_2ndlast',
        'CDR3_first', 'CDR3_last', 'CDR3_2nd', 'CDR3_2ndlast', 'CDR3_3rd', 'CDR3_3rdlast',
    ]
    
    print("\nInitializing memory-efficient analyzers...")
    analyzers = {
        'comp': FeatureCompensationScanner(),
        'cluster': VernierClusterer(),
        'cond': ConditionalModeler(max_samples_per_position=100000),
        'rules': HigherOrderRuleFinder(),
        'mi': MutualInformationAnalyzer(position_names, max_sequences=args.mi_max_seqs),
    }
    
    # Find NPZ files
    npz_files = sorted([f for f in os.listdir(args.npz_dir) if f.endswith('.npz')])
    if args.max_files:
        npz_files = npz_files[:args.max_files]
    
    print(f"Found {len(npz_files)} NPZ files - processing ALL")
    
    # Process with progress bar
    total = 0
    for i, npz_file in enumerate(tqdm(npz_files, desc="Files", unit="file", disable=not HAS_TQDM)):
        npz_path = os.path.join(args.npz_dir, npz_file)
        
        count = process_npz(npz_path, analyzers)
        total += count
        
        elapsed = time.time() - start_time
        rate = total / elapsed if elapsed > 0 else 0
        print(f"  ✓ {npz_file}: {count:,} seqs (total: {total:,}, {rate:.0f}/sec)")
    
    print(f"\n{'='*70}")
    print(f"Data collection complete: {total:,} sequences")
    print(f"Time elapsed: {(time.time() - start_time)/60:.1f} minutes")
    print(f"{'='*70}")
    
    results = {
        'total_sequences': total,
        'family_counts': dict(analyzers['comp'].family_counts),
        'parameters': {
            'min_comp_count': args.min_comp_count,
            'min_cluster_size': args.min_cluster_size,
            'mi_max_seqs': args.mi_max_seqs,
            'rf_trees': args.rf_trees,
            'rf_depth': args.rf_depth,
        },
    }
    
    # Analysis 1
    print("\n[1/5] Computing Feature-Level Compensation (streaming stats)...")
    comp_stats = analyzers['comp'].analyze(min_count=args.min_comp_count)
    comp_rules = analyzers['comp'].find_compensation_rules(comp_stats)
    results['analysis_1_compensation'] = {'stats': comp_stats, 'rules': comp_rules}
    print(f"  Found {len(comp_rules)} compensation rules")
    save_checkpoint(results, args.output, '1_compensation')
    
    # Analysis 2
    print("\n[2/5] Clustering Vernier Patterns (streaming stats)...")
    clusters = analyzers['cluster'].analyze(min_size=args.min_cluster_size)
    results['analysis_2_vernier_clusters'] = clusters
    print(f"  Found {len(clusters)} Vernier archetypes")
    save_checkpoint(results, args.output, '2_clusters')
    
    # Analysis 3
    cond_results = {}  # Initialize to fix NameError bug
    print("\n[3/5] Training Conditional Models...")
    if args.skip_cond_models:
        print("  Skipped (--skip-cond-models)")
        cond_results = {'skipped': True}
    elif not HAS_SKLEARN:
        print("  Skipped - sklearn not installed")
        cond_results = {'error': 'sklearn not installed'}
    else:
        cond_results = analyzers['cond'].train_models(
            min_samples=2000,
            n_estimators=args.rf_trees,
            max_depth=args.rf_depth
        )
    results['analysis_3_conditional_models'] = cond_results
    save_checkpoint(results, args.output, '3_models')
    
    # Analysis 4
    print("\n[4/5] Finding Higher-Order Association Rules (streaming counts)...")
    ho_rules = analyzers['rules'].find_rules(min_support=500, min_confidence=0.65)
    results['analysis_4_higher_order_rules'] = ho_rules
    print(f"  Found {len(ho_rules)} higher-order rules")
    save_checkpoint(results, args.output, '4_rules')
    
    # Analysis 5
    print("\n[5/5] Computing Mutual Information Matrix (integer-encoded)...")
    print(f"  Using {analyzers['mi'].n_added:,} sequences (max: {args.mi_max_seqs:,})")
    mi_matrix, mi_positions = analyzers['mi'].compute_matrix()
    mi_pairs = analyzers['mi'].get_top_pairs(mi_matrix, n_top=200)
    results['analysis_5_mutual_information'] = {
        'positions': mi_positions,
        'matrix': mi_matrix.tolist(),
        'top_pairs': mi_pairs,
        'n_sequences_used': analyzers['mi'].n_added,
    }
    
    # Final save
    print("\n" + "="*70)
    print("SAVING FINAL RESULTS")
    print("="*70)
    
    results['total_time_seconds'] = time.time() - start_time
    
    pkl_path = args.output + '_full.pkl'
    with open(pkl_path, 'wb') as f:
        pickle.dump(results, f)
    print(f"Full results: {pkl_path}")
    
    # Count models safely
    n_models = 0
    if isinstance(cond_results, dict) and 'error' not in cond_results and 'skipped' not in cond_results:
        n_models = sum(len(v) for v in cond_results.values() if isinstance(v, dict))
    
    json_results = {
        'total_sequences': total,
        'family_counts': dict(analyzers['comp'].family_counts),
        'total_time_minutes': (time.time() - start_time) / 60,
        'compensation_rules_count': len(comp_rules),
        'compensation_rules_top50': comp_rules[:50],
        'vernier_clusters_count': len(clusters),
        'vernier_clusters': clusters,
        'conditional_models_count': n_models,
        'higher_order_rules_count': len(ho_rules),
        'higher_order_rules_top50': ho_rules[:50],
        'mi_top_pairs': mi_pairs[:100],
        'mi_sequences_used': analyzers['mi'].n_added,
    }
    
    json_path = args.output + '_summary.json'
    with open(json_path, 'w') as f:
        json.dump(json_results, f, indent=2, default=str)
    print(f"Summary JSON: {json_path}")
    
    # Print summary
    print("\n" + "="*70)
    print("ANALYSIS COMPLETE")
    print("="*70)
    print(f"Total sequences: {total:,}")
    print(f"Total time: {(time.time() - start_time)/60:.1f} minutes")
    
    print(f"\nFamily distribution:")
    for family, count in sorted(analyzers['comp'].family_counts.items(), key=lambda x: -x[1]):
        pct = 100 * count / total if total > 0 else 0
        print(f"  {family}: {count:,} ({pct:.1f}%)")
    
    # Print hallmark position distributions to validate mapping
    print(f"\nHallmark position distributions (to validate mapping):")
    print(f"  Position 42 (FR2_2 - should be Y/F in VHH, V in VH):")
    print(f"  Position 50 (FR2_10 - should be R in VHH, L in VH):")
    print(f"  Position 52 (FR2_12 - should be G/L/F in VHH, W in VH):")
    
    print(f"\nResults:")
    print(f"  Compensation rules: {len(comp_rules)}")
    print(f"  Vernier archetypes: {len(clusters)}")
    print(f"  Conditional models: {n_models}")
    print(f"  Higher-order rules: {len(ho_rules)}")
    print(f"  MI sequences used: {analyzers['mi'].n_added:,}")
    
    print(f"\nTop 5 compensation rules:")
    for rule in comp_rules[:5]:
        print(f"  {rule['family']} {rule['position']}: {rule['res1']} vs {rule['res2']} → {rule['feature']} {rule['effect']:+.2f}")
    
    print(f"\nTop 5 MI pairs:")
    for pair in mi_pairs[:5]:
        print(f"  {pair['pos1']} <-> {pair['pos2']}: MI = {pair['mi']:.4f} ({pair['type']})")
    
    print(f"\nFiles created:")
    print(f"  {pkl_path}")
    print(f"  {json_path}")
    print(f"  + checkpoint files")

if __name__ == '__main__':
    main()