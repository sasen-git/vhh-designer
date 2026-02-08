#!/usr/bin/env python3
"""
VHH Designer v8.3 - Probabilistic Consensus Sprinkling & Motif Tracks
======================================================================

UPDATES from v8.2:

1. PROBABILISTIC CONSENSUS SPRINKLING (v8.3):
   - Build universal_prior[pos] = (top_aa, freq) aggregated across families
   - Tier A (≥90%): Apply with ~60% probability
   - Tier B (75-90%): Apply with ~25% probability
   - No more "all-or-nothing" - intermediate states now generated

2. MOTIF-BASED TRACKS (v8.3):
   - Track 2: Motif PAIRS (1-3 co-occurring pairs per candidate)
   - Track 3A: Motif TRIPLETS (1-2 triplets per candidate)
   - Track 3B: Partial consensus BUNDLES (3-6 Tier-A positions, sampled)
   - Not random verniers - uses cross-family enriched motifs

3. WITHIN-TRACK RANKING (v8.3):
   - Tracks 1-3 ranked by ESM+rules WITHIN their track
   - Track 4 doesn't dominate - diagnostic value preserved
   - Controls are still controls, not competing with optimized

4. VERSION IN FILENAMES (v8.3):
   - Output files now use VERSION constant (v8_3)
   - Example: 20260121_200000_v8_3_n20000_M69

5. RETAINED FROM v8.2:
   - Universal as subfamily with equal Track 4 allocation
   - Clear ranking format: "(Control)" vs 1/2/3...
   - 198 total output (99 controls + 99 ranked)

Usage:
  python vhh_designer_v8_3.py -i input.fasta \\
      --rules analysis_rules_v7.json \\
      --archetypes analysis_vernier_archetypes_v7.json \\
      --hallmark-db comprehensive_subfamily_analysis_imgt.xlsx \\
      --target-hallmarks AUTO \\
      --n-generate 20000 \\
      --n-select 198

Author: Claude (Anthropic)
Date: January 2026
"""

import os
import sys
import re
import json
import argparse
import random
import numpy as np
from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Optional, Set, Any
from collections import defaultdict
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

# Optional imports
try:
    from antpack import SingleChainAnnotator
    ANTPACK_AVAILABLE = True
except ImportError:
    ANTPACK_AVAILABLE = False
    print("Warning: AntPack not available. Install: pip install antpack")

try:
    import torch
    TORCH_AVAILABLE = True
except ImportError:
    TORCH_AVAILABLE = False
    print("Note: PyTorch not available. Install: pip install torch")

try:
    from transformers import AutoModelForMaskedLM, AutoTokenizer
    ESM2_AVAILABLE = True
except ImportError:
    ESM2_AVAILABLE = False
    print("Note: ESM2 not available. Install: pip install transformers")

ESMFOLD_AVAILABLE = False
try:
    import esm
    ESMFOLD_AVAILABLE = True
except ImportError:
    print("Note: ESMFold not available. Install: pip install fair-esm")

try:
    import pandas as pd
    PANDAS_AVAILABLE = True
except ImportError:
    PANDAS_AVAILABLE = False

try:
    from tqdm import tqdm
    TQDM_AVAILABLE = True
except ImportError:
    TQDM_AVAILABLE = False
    def tqdm(x, **kwargs):
        return x


# ============================================================
# HALLMARK / SUBFAMILY DATABASE (from 12M VHH analysis)
# ============================================================

class HallmarkDB:
    """Loads hallmark/subfamily statistics from comprehensive_subfamily_analysis_imgt.xlsx.

    We use this to:
      • exclude VH-like hallmarks (pos50=L) and W52 hallmarks (pos52=W) in true-VHH mode
      • pick a diverse, high-support pool of reference hallmarks using precomputed clusters
      • provide hallmark-specific vernier consensus + 2nd-choice frequencies for sampling
      • provide hallmark-specific cysteine-class / CDR3 length priors

    Expected sheets:
      1_Summary, 2_Vernier_Consensus, 11_Clustering, 6_Epistatic_Pairs
    """

    
    def __init__(self, xlsx_path: str):
        import pandas as pd
        self.xlsx_path = xlsx_path
        self.summary = pd.read_excel(xlsx_path, sheet_name='1_Summary')
        self.vernier = pd.read_excel(xlsx_path, sheet_name='2_Vernier_Consensus')
        self.clusters = pd.read_excel(xlsx_path, sheet_name='11_Clustering')
        # epistatic pairs are optional (large)
        try:
            self.epi = pd.read_excel(xlsx_path, sheet_name='6_Epistatic_Pairs')
        except Exception:
            self.epi = None

        # Normalize key columns
        self.summary['hallmark'] = self.summary['hallmark'].astype(str)
        self.vernier['hallmark'] = self.vernier['hallmark'].astype(str)
        self.clusters['hallmark'] = self.clusters['hallmark'].astype(str)

        self._vernier_by_hm = {r['hallmark']: r for _, r in self.vernier.iterrows()}
        self._summary_by_hm = {r['hallmark']: r for _, r in self.summary.iterrows()}

    @staticmethod
    def _is_true_vhh_hallmark(hm: str) -> bool:
        """Hard exclusion rules for true VHH mode.
        
        Key insight: Position 50 (R vs L) is more important than position 52
        for VHH functionality. R at 50 provides solubility for independent function.
        
        Position 52 spectrum:
          W = VH-like (wants VL partner) - EXCLUDE
          L/F = intermediate (still functional as VHH) - ALLOW
          G/A = classic camelization - ALLOW
        """
        if not hm or len(hm) != 4:
            return False
        p42, p49, p50, p52 = hm[0], hm[1], hm[2], hm[3]
        
        # CRITICAL: Require R at position 50 (VHH solubility/independence)
        if p50 != 'R':
            return False
        
        # Exclude W at position 52 (too VH-like, wants VL partner)
        # But allow G, A, L, F (all functional as VHH)
        if p52 == 'W':
            return False
        
        # Require aromatic at 42 for strong VHH-ness (F/Y)
        if p42 not in ('F', 'Y'):
            return False
        
        # Prefer charged/polar at 49 (avoid human-like G)
        if p49 == 'G':
            return False
        
        return True

    def eligible_hallmarks(self,
                           min_n: int = 50000,
                           require_true_vhh: bool = True) -> List[str]:
        """Return hallmarks passing hard filters and support."""
        df = self.summary.copy()
        df = df[df['n_sequences'] >= min_n]
        if require_true_vhh:
            df = df[df['hallmark'].apply(self._is_true_vhh_hallmark)]
        return df['hallmark'].tolist()

    def hallmark_meta(self, hm: str) -> dict:
        r = self._summary_by_hm.get(hm)
        if r is None:
            return {}
        return {
            'hallmark': hm,
            'n_sequences': int(r.get('n_sequences', 0)),
            'type': str(r.get('type', '')),
            'cys_class': str(r.get('cys_class', '')),
            'c4_pct': float(r.get('c4_pct', 0.0)),
            'c2_pct': float(r.get('c2_pct', 0.0)),
            'cdr3_mean_len_imgt': float(r.get('cdr3_mean_len_imgt', 0.0)),
            'cdr3_std_len_imgt': float(r.get('cdr3_std_len_imgt', 0.0)),
            'avg_vernier_conservation': float(r.get('avg_vernier_conservation', 0.0)),
            'avg_fw_entropy': float(r.get('avg_fw_entropy', 0.0)),
            'avg_cdr_entropy': float(r.get('avg_cdr_entropy', 0.0)),
        }

    def hallmark_to_family(self, hm: str) -> str:
        """Map hallmark to one of the coarse archetype families used by rules/archetypes."""
        meta = self.hallmark_meta(hm)
        t = meta.get('type', hm[0] if hm else 'F')
        cys = meta.get('cys_class', 'C2')
        if t in ('F', 'Y') and cys in ('C2', 'C4'):
            return f"{t}_{cys}"
        # fallback
        return 'Other_VHH'

    def pick_reference_pool(self,
                            k_per_cluster: int = 2,
                            cluster_level: str = 'cluster_10',
                            min_n: int = 50000) -> List[str]:
        """Pick a diverse hallmark pool: top-k by n_sequences per cluster among eligible hallmarks."""
        eligible = set(self.eligible_hallmarks(min_n=min_n, require_true_vhh=True))
        df = self.clusters.copy()
        df = df[df['hallmark'].isin(eligible)]
        if cluster_level not in df.columns:
            cluster_level = 'cluster_10'
        pool = []
        for cid, g in df.groupby(cluster_level):
            g = g.sort_values('n_sequences', ascending=False).head(k_per_cluster)
            pool.extend(g['hallmark'].tolist())
        # stable order: highest-support first
        pool = sorted(set(pool), key=lambda hm: -self.hallmark_meta(hm).get('n_sequences', 0))
        return pool

    def vernier_profile(self, hm: str) -> dict:
        """Return hallmark-specific vernier profile (consensus + 2nd choice + freqs) for sampling.
        
        Includes ALL vernier positions from sheet 2: 2, 4, 41, 47, 52, 66, 67, 68, 69, 71, 76, 78, 82, 87, 89, 91, 94
        """
        r = self._vernier_by_hm.get(hm)
        if r is None:
            return {}
        prof = {}
        # All vernier positions from sheet 2 (2_Vernier_Consensus)
        VERNIER_POSITIONS_SHEET2 = [2, 4, 41, 47, 52, 66, 67, 68, 69, 71, 76, 78, 82, 87, 89, 91, 94]
        for pos in VERNIER_POSITIONS_SHEET2:
            if pos in HALLMARK_POSITIONS and pos != 52:  # Skip hallmarks except 52 which is also vernier
                continue
            k1 = f"V{pos}_cons"
            k1f = f"V{pos}_freq"
            k2 = f"V{pos}_2nd"
            k2f = f"V{pos}_2nd_freq"
            if k1 in r:
                prof[pos] = {
                    'cons': str(r.get(k1, '')),
                    'freq': float(r.get(k1f, 0.0) or 0.0),
                    'alt': str(r.get(k2, '')),
                    'alt_freq': float(r.get(k2f, 0.0) or 0.0),
                }
        return prof

    def top_epistatic_pairs(self, hm: str, max_pairs: int = 4) -> List[Tuple[int, int, List[Tuple[str, str, float]]]]:
        """Return a small set of coupled vernier pairs with their top combinations.

        Output: [(pos1,pos2, [(aa1,aa2,freq), ...]), ...]
        """
        if self.epi is None:
            return []
        df = self.epi[self.epi['hallmark'] == hm]
        # Keep only vernier/vernier pairs
        df = df[df['pos1'].isin(ALL_VERNIER_POSITIONS) & df['pos2'].isin(ALL_VERNIER_POSITIONS)]
        # choose strongest pairs by best rank
        pairs = []
        for (p1,p2), g in df.groupby(['pos1','pos2']):
            g = g.sort_values('rank', ascending=True)
            top = g.head(6)
            combs = [(str(r['aa1']), str(r['aa2']), float(r['frequency'])) for _, r in top.iterrows()]
            best_rank = int(top['rank'].min()) if len(top) else 10**9
            pairs.append((best_rank, int(p1), int(p2), combs))
        pairs.sort(key=lambda x: x[0])
        out=[]
        for _, p1,p2, combs in pairs[:max_pairs]:
            out.append((p1,p2, combs))
        return out

    def identify_cross_family_vernier_motifs(self, min_consensus: float = 0.75, min_families: int = 3) -> Dict[str, Any]:
        """
        Identify vernier positions and pairs that are conserved across multiple hallmarks.
        
        Returns:
            {
                'single_verniers': [(pos, aa, avg_freq, hallmark_count), ...],
                'paired_motifs': [((pos1, pos2), (aa1, aa2), avg_freq, hallmark_count), ...],
                'triplet_motifs': [((pos1, pos2, pos3), (aa1, aa2, aa3), avg_freq, hallmark_count), ...]
            }
        """
        # Get all true VHH hallmarks
        true_vhh_hallmarks = [hm for hm in self._vernier_by_hm.keys() if self._is_true_vhh_hallmark(hm)]
        
        # Track consensus at each position across hallmarks
        position_consensus = {}  # pos -> {aa: [(hallmark, freq), ...]}
        
        FR3_VERNIER = [66, 67, 68, 69, 71, 76, 78, 82, 87, 89, 91, 94]
        
        for hm in true_vhh_hallmarks:
            profile = self.vernier_profile(hm)
            for pos in FR3_VERNIER:
                data = profile.get(pos, {})
                cons_aa = data.get('cons', '')
                freq = data.get('freq', 0)
                if cons_aa and freq >= min_consensus:
                    if pos not in position_consensus:
                        position_consensus[pos] = {}
                    if cons_aa not in position_consensus[pos]:
                        position_consensus[pos][cons_aa] = []
                    position_consensus[pos][cons_aa].append((hm, freq))
        
        # Identify single verniers conserved across families
        single_verniers = []
        for pos, aa_dict in position_consensus.items():
            for aa, hm_list in aa_dict.items():
                if len(hm_list) >= min_families:
                    avg_freq = sum(f for _, f in hm_list) / len(hm_list)
                    single_verniers.append((pos, aa, avg_freq, len(hm_list), [h for h, _ in hm_list]))
        
        single_verniers.sort(key=lambda x: (-x[3], -x[2]))  # Sort by count, then freq
        
        # Identify paired motifs (positions that have same AA pair across families)
        paired_motifs = []
        positions = sorted(position_consensus.keys())
        for i, pos1 in enumerate(positions):
            for pos2 in positions[i+1:]:
                # Check which hallmarks have high consensus at BOTH positions
                for aa1, hm_list1 in position_consensus[pos1].items():
                    for aa2, hm_list2 in position_consensus[pos2].items():
                        # Find hallmarks that have both
                        hm_set1 = {h for h, _ in hm_list1}
                        hm_set2 = {h for h, _ in hm_list2}
                        shared = hm_set1 & hm_set2
                        if len(shared) >= min_families:
                            # Calculate average frequency for shared hallmarks
                            freq1 = sum(f for h, f in hm_list1 if h in shared) / len(shared)
                            freq2 = sum(f for h, f in hm_list2 if h in shared) / len(shared)
                            avg_freq = (freq1 + freq2) / 2
                            paired_motifs.append(((pos1, pos2), (aa1, aa2), avg_freq, len(shared), list(shared)))
        
        paired_motifs.sort(key=lambda x: (-x[3], -x[2]))  # Sort by count, then freq
        
        # Identify triplet motifs (for very high-confidence structural units)
        triplet_motifs = []
        for i, pos1 in enumerate(positions):
            for j, pos2 in enumerate(positions[i+1:], i+1):
                for pos3 in positions[j+1:]:
                    for aa1, hm_list1 in position_consensus[pos1].items():
                        for aa2, hm_list2 in position_consensus[pos2].items():
                            for aa3, hm_list3 in position_consensus[pos3].items():
                                hm_set1 = {h for h, _ in hm_list1}
                                hm_set2 = {h for h, _ in hm_list2}
                                hm_set3 = {h for h, _ in hm_list3}
                                shared = hm_set1 & hm_set2 & hm_set3
                                if len(shared) >= min_families:
                                    freq1 = sum(f for h, f in hm_list1 if h in shared) / len(shared)
                                    freq2 = sum(f for h, f in hm_list2 if h in shared) / len(shared)
                                    freq3 = sum(f for h, f in hm_list3 if h in shared) / len(shared)
                                    avg_freq = (freq1 + freq2 + freq3) / 3
                                    triplet_motifs.append(((pos1, pos2, pos3), (aa1, aa2, aa3), avg_freq, len(shared), list(shared)))
        
        triplet_motifs.sort(key=lambda x: (-x[3], -x[2]))
        
        return {
            'single_verniers': single_verniers,
            'paired_motifs': paired_motifs[:20],  # Limit to top 20
            'triplet_motifs': triplet_motifs[:10],  # Limit to top 10
        }

    def get_high_consensus_verniers_for_hallmark(self, hallmark: str, min_freq: float = 0.85) -> List[Tuple[int, str, float]]:
        """
        Get vernier positions with high consensus for a specific hallmark.
        
        Returns: [(pos, consensus_aa, frequency), ...]
        """
        profile = self.vernier_profile(hallmark)
        result = []
        FR3_VERNIER = [66, 67, 68, 69, 71, 76, 78, 82, 87, 89, 91, 94]
        
        for pos in FR3_VERNIER:
            data = profile.get(pos, {})
            cons_aa = data.get('cons', '')
            freq = data.get('freq', 0)
            if cons_aa and freq >= min_freq:
                result.append((pos, cons_aa, freq))
        
        result.sort(key=lambda x: -x[2])  # Sort by frequency descending
        return result


# ============================================================
# CONSTANTS
# ============================================================

UNIVERSAL_SCAFFOLD = {
    'name': 'Universal Humanized VHH',
    'FR1': 'QVQLVESGGGLVQPGGSLRLSCAASG',
    'FR2': 'WFRQAPGQGLEAVA',
    'FR3': 'YYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYC',
    'FR4': 'WGQGTLVTVSS',
}

# IMGT region definitions (position ranges)
IMGT_REGION_RANGES = {
    'FR1': (1, 26),
    'CDR1': (27, 38),
    'FR2': (39, 55),
    'CDR2': (56, 65),
    'FR3': (66, 104),
    'CDR3': (105, 117),
    'FR4': (118, 128),
}

# Hallmark positions (IMGT numbering)
HALLMARK_POSITIONS = [42, 49, 50, 52]

# Common VHH hallmark patterns
VHH_HALLMARKS_IMGT = {
    'FERG': {42: 'F', 49: 'E', 50: 'R', 52: 'G'},
    'FERF': {42: 'F', 49: 'E', 50: 'R', 52: 'F'},
    'FERA': {42: 'F', 49: 'E', 50: 'R', 52: 'A'},
    'YERL': {42: 'Y', 49: 'E', 50: 'R', 52: 'L'},
    'YQRL': {42: 'Y', 49: 'Q', 50: 'R', 52: 'L'},
    'FQRL': {42: 'F', 49: 'Q', 50: 'R', 52: 'L'},
    'FDRF': {42: 'F', 49: 'D', 50: 'R', 52: 'F'},
    'VGLW': {42: 'V', 49: 'G', 50: 'L', 52: 'W'},
    'FGLA': {42: 'F', 49: 'G', 50: 'L', 52: 'A'},
}

# Vernier positions (FR3 verniers that affect CDR conformation)
VERNIER_POSITIONS_FR3 = {66, 67, 68, 69, 71, 76, 78, 82, 87, 89, 91, 94}

# All vernier positions including FR2 hallmarks
ALL_VERNIER_POSITIONS = {42, 49, 50, 52} | VERNIER_POSITIONS_FR3

# Positions that are extremely conserved in true VHH datasets and should almost
# never be mutated when building "very likely" candidates.
# (From vernier_conservation_report.md: >=95% conserved)
INVARIANT_VHH_POSITIONS = {4, 41, 47, 76, 89, 94}

# ============================================================
# V8.0 CONSTANTS
# ============================================================

VERSION = "8.3"

# Load-bearing verniers for focused Track 1 testing (v8.0)
LOAD_BEARING_VERNIERS = [15, 66, 68, 69, 71, 76, 78, 82, 87, 89, 91, 94]

# Key positions for fractional factorial design (Track 3B)
FACTORIAL_POSITIONS = [66, 68, 69, 71, 76, 78]

# Universal Framework (humanized VHH scaffold for CDR grafting)
UNIVERSAL_FRAMEWORK = {
    'FR1': 'QVQLVESGGGLVQPGGSLRLSCAASG',
    'FR2': 'WFRQAPGQGLEAVA',
    'FR3': 'YYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYC',
    'FR4': 'WGQGTLVTVSS',
}

# v8.1: MANDATORY VHH consensus positions (applied like IMGT2→V)
# These are near-universal in camelid VHH sequences
MANDATORY_VHH_CONSENSUS = {
    1: ('E', 'Q', 0.92),    # FR1 start: E→Q (~92% in VHH)
    2: ('I', 'V', 1.00),    # FR1: I→V (100% in VHH) - already handled but included for completeness
    128: ('A', 'S', 0.95),  # FR4 end: ...TVSA→...TVSS (~95% in VHH)
}

# v8.1: Track allocation for 186 total sequences
# 92 controls (unranked) + 94 ranked
# v8.2: Output budget - 198 total (99 controls + 99 ranked)
TOTAL_OUTPUT = 198
CONTROL_BUDGET = 99
RANKED_BUDGET = 99

# ============================================================
# V8.3: PROBABILISTIC CONSENSUS SPRINKLING
# ============================================================

# Universal prior: aggregated consensus across major VHH families
# Format: {imgt_pos: (consensus_aa, frequency, tier)}
# Tier A: ≥90%, Tier B: 75-90%
# These are derived from FERG, FERF, FERA, YERL, YQRL analysis
UNIVERSAL_PRIOR = {
    # FR3 Vernier positions (most important for CDR support)
    68: ('A', 0.88, 'B'),   # N→A common across families
    69: ('D', 0.82, 'B'),   # E→D common 
    71: ('V', 0.78, 'B'),   # F→V common
    76: ('F', 0.95, 'A'),   # Nearly invariant F
    78: ('I', 0.85, 'B'),   # F→I common
    82: ('N', 0.76, 'B'),   # Variable but N common
    87: ('V', 0.80, 'B'),   # Y/L→V 
    89: ('L', 0.92, 'A'),   # M→L highly conserved
    91: ('M', 0.77, 'B'),   # L/N→M
    94: ('R', 0.88, 'B'),   # S→R common
    # FR1/FR4 (already handled but included for completeness)
    1: ('Q', 0.92, 'A'),    # E→Q
    2: ('V', 1.00, 'A'),    # I→V
    128: ('S', 0.95, 'A'),  # A→S (TVSS ending)
}

# Sprinkling probabilities by tier
TIER_A_PROB = 0.60  # Apply Tier A positions with 60% probability
TIER_B_PROB = 0.25  # Apply Tier B positions with 25% probability

# Cross-family enriched motif PAIRS (positions that change together)
# Format: [(pos1, aa1, pos2, aa2, frequency_across_families)]
MOTIF_PAIRS = [
    (68, 'A', 69, 'D', 0.85),   # A68+D69 co-occur in most families
    (69, 'D', 71, 'V', 0.82),   # D69+V71 co-occur
    (71, 'V', 76, 'F', 0.90),   # V71+F76 very common
    (76, 'F', 78, 'I', 0.88),   # F76+I78 common
    (78, 'I', 89, 'L', 0.83),   # I78+L89 common
    (87, 'V', 91, 'M', 0.79),   # V87+M91 often co-occur
    (89, 'L', 94, 'R', 0.85),   # L89+R94 common
    (68, 'A', 71, 'V', 0.80),   # A68+V71 
    (82, 'N', 87, 'V', 0.75),   # N82+V87
]

# Cross-family enriched motif TRIPLETS
# Format: [(pos1, aa1, pos2, aa2, pos3, aa3, frequency)]
MOTIF_TRIPLETS = [
    (68, 'A', 69, 'D', 71, 'V', 0.78),   # Core vernier triplet
    (71, 'V', 76, 'F', 78, 'I', 0.85),   # Structural core
    (76, 'F', 78, 'I', 89, 'L', 0.82),   # FR3 end
    (87, 'V', 89, 'L', 91, 'M', 0.75),   # Position 87-91 cluster
    (68, 'A', 71, 'V', 76, 'F', 0.80),   # Extended core
]

def _aa_group(aa: str) -> str:
    """
    Coarse AA grouping for 'similarity' scoring.
    """
    aa = (aa or "").upper()
    if aa in "AVILM":
        return "hydrophobic"
    if aa in "FYW":
        return "aromatic"
    if aa in "STNQ":
        return "polar"
    if aa in "KRH":
        return "positive"
    if aa in "DE":
        return "negative"
    if aa in "GPC":
        return "special"
    return "other"


def compute_hallmark_match(imgt_positions: dict, expected_hallmarks_key: str,
                           same_weight: float = 1.0, similar_weight: float = 0.5) -> dict:
    """
    Compare candidate IMGT hallmarks (42,49,50,52) to an expected hallmark key like 'FERG' or 'YER(A/G)'.

    Returns a dict with:
      - score: normalized [0,1]
      - per_pos: {pos: {"obs":..., "exp":..., "match": "same|similar|miss"}}
      - obs_key: observed hallmark string (e.g., 'IGLW')
    """
    # Expected hallmark string must be length 4 in canonical usage (e.g. FERG)
    exp = (expected_hallmarks_key or "").upper().strip()
    # allow simple normalization like "FERA" etc; if not length 4, treat as unknown
    if len(exp) != 4:
        exp = exp[:4]

    hallmark_positions = [42, 49, 50, 52]
    obs = []
    per_pos = {}
    total = 0.0
    denom = 0.0

    for i, pos in enumerate(hallmark_positions):
        obs_aa = (imgt_positions.get(pos) or "-").upper()
        exp_aa = exp[i] if i < len(exp) else "-"
        obs.append(obs_aa)

        match_type = "miss"
        w = 0.0
        if obs_aa == exp_aa and obs_aa != "-":
            match_type = "same"
            w = same_weight
        elif obs_aa != "-" and exp_aa != "-" and _aa_group(obs_aa) == _aa_group(exp_aa):
            match_type = "similar"
            w = similar_weight

        per_pos[pos] = {"obs": obs_aa, "exp": exp_aa, "match": match_type}
        total += w
        denom += same_weight  # normalize to max possible (all same)

    score = (total / denom) if denom > 0 else 0.0
    return {"score": score, "per_pos": per_pos, "obs_key": "".join(obs)}

def expected_hallmarks_for_family(family: str) -> Optional[str]:
    """Expected hallmark set based on VHH family naming.

    Returns a key into VHH_HALLMARKS_IMGT (e.g. 'FERG', 'YERL'), or None
    when the family doesn't map cleanly (e.g. VH_like, unknown).
    """
    if family.startswith("F_") or family == "Other_VHH":
        return "FERG"
    if family.startswith("Y_"):
        return "YERL"
    if family == "VHH_W52":
        # The defining property is W at 52; use VGLW as a reasonable template.
        return "VGLW"
    return None

def hallmark_match_fraction(imgt: Dict[int, str], family: str) -> float:
    expected_key = expected_hallmarks_for_family(family)
    if not expected_key:
        return 1.0
    expected = VHH_HALLMARKS_IMGT.get(expected_key, {})
    if not expected:
        return 1.0
    ok = 0
    for pos in HALLMARK_POSITIONS:
        aa = imgt.get(pos)
        if aa and aa == expected.get(pos):
            ok += 1
    return ok / len(HALLMARK_POSITIONS)

def vernier_match_fraction(imgt: Dict[int, str], archetype: Optional[Dict[str, Any]]) -> float:
    """Fraction of core vernier positions matching the family's archetype consensus."""
    if not archetype:
        return 0.0
    positions = archetype.get("positions", {})
    if not positions:
        return 0.0
    ok = 0
    total = 0
    for pos in sorted(VERNIER_POSITIONS_FR3):
        key = f"IMGT{pos}"
        cons = positions.get(key, {}).get("consensus")
        aa = imgt.get(pos)
        if not cons or not aa:
            continue
        total += 1
        if aa == cons:
            ok += 1
    return (ok / total) if total else 0.0

# ============================================================
# DATA CLASSES
# ============================================================

@dataclass
class IMGTNumberedSequence:
    """
    Holds an antibody sequence with IMGT numbering.
    
    The key data structure is `positions`: a dict mapping IMGT position (int or str like "52A")
    to amino acid. This allows insertion-safe access to any position.
    """
    positions: Dict[str, str]  # IMGT position -> AA (e.g., {1: 'Q', 2: 'V', ..., '52A': 'G'})
    sequence: str              # Full trimmed sequence
    chain_type: str            # 'H', 'K', 'L'
    regions: Dict[str, str]    # Region name -> sequence (e.g., {'FR1': 'QVQ...', 'CDR1': 'GFT...'})
    
    def get_aa(self, imgt_pos: int) -> str:
        """Get AA at IMGT position. Returns '' if not present."""
        # Try integer first, then string (for insertions)
        if imgt_pos in self.positions:
            return self.positions[imgt_pos]
        if str(imgt_pos) in self.positions:
            return self.positions[str(imgt_pos)]
        return ''
    
    def get_hallmarks(self) -> str:
        """Get hallmark pattern (positions 42, 49, 50, 52)."""
        p42 = self.get_aa(42) or '?'
        p49 = self.get_aa(49) or '?'
        p50 = self.get_aa(50) or '?'
        p52 = self.get_aa(52) or '?'
        return f"{p42}{p49}{p50}{p52}"

@dataclass
class CDRSet:
    """CDR and FR regions extracted from a sequence."""
    cdr1: str
    cdr2: str
    cdr3: str
    fr1: str = ""
    fr2: str = ""
    fr3: str = ""
    fr4: str = ""
    imgt_numbered: Optional[IMGTNumberedSequence] = None  # Full IMGT numbering

@dataclass
class Mutation:
    position: str       # e.g., "IMGT66" or "IMGT52A"
    imgt_num: Any       # IMGT position (int or str like "52A")
    original: str
    mutant: str
    source: str
    confidence: float
    
    def __str__(self):
        return f"IMGT{self.imgt_num}:{self.original}→{self.mutant}"

@dataclass
class ScoringResult:
    """Comprehensive scoring result for a candidate."""
    # ESM2 Language Model
    esm2_loss: float = 0.0
    esm2_perplexity: float = 0.0
    
    # ESMFold Structure Prediction
    plddt_mean: float = 0.0
    plddt_median: float = 0.0
    plddt_min: float = 0.0
    plddt_cdr1: float = 0.0
    plddt_cdr2: float = 0.0
    plddt_cdr3: float = 0.0
    plddt_framework: float = 0.0
    
    # Multi-family probabilistic
    family_probabilities: Dict[str, float] = field(default_factory=dict)
    rule_compliance: Dict[str, float] = field(default_factory=dict)
    weighted_naturalness: float = 0.0
    
    # Rule counts
    rules_passed: int = 0
    rules_total: int = 0
    rules_applicable: int = 0
    
    # Vernier matches
    vernier_matches: int = 0
    vernier_total: int = 0
    
    # Violations
    top_violations: List[str] = field(default_factory=list)
    
    # Combined
    combined_score: float = 0.0

    # Convenience confidence for ranking/reporting (0–100).
    # NOTE: this is not a calibrated probability; it’s a weighted summary.
    confidence_pct: float = 0.0
    confidence_components: Dict[str, Any] = field(default_factory=dict)

@dataclass
class VHHCandidate:
    id: str
    rank: int
    sequence: str
    framework_source: str
    family: str
    cdr1: str
    cdr2: str
    cdr3: str
    fr1: str
    fr2: str
    fr3: str
    fr4: str
    imgt_positions: Dict[Any, str] = field(default_factory=dict)  # IMGT position -> AA
    mutations: List[Mutation] = field(default_factory=list)
    strategy: str = ""
    is_lead: bool = False
    generation_order: int = 0
    construction_method: str = ""
    target_family: str = ""
    framework_identity_pct: float = 100.0
    scoring: ScoringResult = field(default_factory=ScoringResult)
    
    # Design track system (v7.12/7.14)
    design_track: str = "optimized"  # lead, minimal_hallmark, single_vernier, paired_vernier, triplet_vernier, optimized
    ranking_exempt: bool = False     # If True, excluded from competitive ranking
    track_info: str = ""             # Additional info about the track (e.g., "YQRL_4mut", "IMGT69_probe")
    track_rank: int = 0              # Within-track rank (v7.14) - enables sorting within control bucket
    is_immutable: bool = False       # If True, no post-generation modifications allowed (v7.14)
    
    # v8.0 additions
    track_subtype: str = ""          # For Track 4 lanes (lane_A, lane_B) or Track 3 (motif_triplet, factorial)
    scaffold_type: str = "original"  # original, universal, yqrl_consensus, family_consensus

    # Provenance/debugging
    generation_debug: Dict[str, Any] = field(default_factory=dict)
    applied_comp_rules: List[str] = field(default_factory=list)
    applied_triplet_rules: List[str] = field(default_factory=list)
    applied_epistatic_pairs: List[Dict[str, Any]] = field(default_factory=list)
    
    def get_aa(self, imgt_pos: Any) -> str:
        """Get AA at IMGT position."""
        return self.imgt_positions.get(imgt_pos, self.imgt_positions.get(str(imgt_pos), ''))


def validate_immutable_candidate(candidate: VHHCandidate, hallmark_positions: Set[int] = None) -> Tuple[bool, List[str]]:
    """
    Validate that an immutable candidate (Track 0) has only expected mutations.
    
    v7.14 seatbelt: Prevents regressions where downstream code accidentally
    modifies Track 0 candidates.
    
    Args:
        candidate: The VHHCandidate to validate
        hallmark_positions: Set of hallmark IMGT positions (default: {42, 49, 50, 52})
    
    Returns:
        (is_valid, list_of_violations)
    """
    if not candidate.is_immutable:
        return True, []
    
    if hallmark_positions is None:
        hallmark_positions = {42, 49, 50, 52}
    
    # Allow hallmarks + IMGT2 (for 5mut variant)
    allowed_positions = hallmark_positions | {2}
    
    violations = []
    
    # Check that all mutations are at allowed positions
    for mut in candidate.mutations:
        pos = mut.imgt_num
        if pos not in allowed_positions:
            violations.append(f"Unexpected mutation at IMGT{pos}: {mut.original}>{mut.mutant}")
    
    # Check that no compensation rules were applied
    if candidate.applied_comp_rules:
        violations.append(f"Compensation rules applied to immutable candidate: {candidate.applied_comp_rules}")
    
    # Check that no triplet rules were applied
    if candidate.applied_triplet_rules:
        violations.append(f"Triplet rules applied to immutable candidate: {candidate.applied_triplet_rules}")
    
    return len(violations) == 0, violations


# ============================================================
# IMGT NUMBERING & CDR EXTRACTION (FIXED)
# ============================================================

def parse_imgt_position(pos_str: str) -> Tuple[int, str]:
    """
    Parse IMGT position string into (base_number, insertion_code).
    
    Examples:
        '52' -> (52, '')
        '52A' -> (52, 'A')
        '111A' -> (111, 'A')
    """
    if not pos_str or pos_str == '-':
        return (-1, '')
    
    # Extract numeric part and insertion code
    match = re.match(r'^(\d+)([A-Z]*)$', str(pos_str).strip())
    if match:
        return (int(match.group(1)), match.group(2))
    return (-1, '')

def get_region_for_position(imgt_num: int) -> str:
    """Determine which region an IMGT position belongs to."""
    for region, (start, end) in IMGT_REGION_RANGES.items():
        if start <= imgt_num <= end:
            return region
    return ''

def extract_cdrs_fixed(sequence: str) -> Optional[CDRSet]:
    """
    Extract CDRs and FRs using AntPack with CORRECT parsing.
    
    This version:
    1. Properly handles AntPack's 4-tuple return
    2. Uses trim_alignment() to get aligned sequence
    3. Builds a position dict for insertion-safe access
    4. Returns regions based on IMGT numbering, not string indices
    """
    if not ANTPACK_AVAILABLE:
        print("ERROR: AntPack required for CDR extraction")
        return None
    
    try:
        annotator = SingleChainAnnotator()
        
        # AntPack returns 4 values: (numbering, percent_identity, chain_type, error_message)
        result = annotator.analyze_seq(sequence)
        
        # Handle both 3-tuple and 4-tuple returns (AntPack version differences)
        if len(result) == 4:
            numbering, percent_id, chain_type, error_msg = result
        elif len(result) == 3:
            numbering, chain_type, error_msg = result
            percent_id = None
        else:
            print(f"ERROR: Unexpected AntPack return: {len(result)} values")
            return None
        
        if error_msg:
            print(f"Warning: AntPack error: {error_msg}")
        
        # Use trim_alignment to get properly aligned sequence
        # This removes leading/trailing gaps and gives us the aligned portion
        try:
            trimmed_numbering, trimmed_seq = annotator.trim_alignment(numbering, sequence)
        except Exception as e:
            # Fallback: use raw alignment if trim_alignment fails
            print(f"Warning: trim_alignment failed ({e}), using raw alignment")
            trimmed_numbering = numbering
            trimmed_seq = sequence
        
        # Build position dict: IMGT position -> AA
        positions = {}
        regions = defaultdict(list)
        
        for pos_label, aa in zip(trimmed_numbering, trimmed_seq):
            if pos_label == '-' or aa == '-':
                continue
            
            # Parse position (handles insertions like '52A')
            base_num, insertion = parse_imgt_position(str(pos_label))
            if base_num < 0:
                continue
            
            # Store with full position key (e.g., 52 or '52A')
            if insertion:
                pos_key = f"{base_num}{insertion}"
            else:
                pos_key = base_num
            
            positions[pos_key] = aa
            
            # Determine region and collect
            region = get_region_for_position(base_num)
            if region:
                regions[region].append(aa)
        
        # Build region sequences
        fr1 = ''.join(regions.get('FR1', []))
        cdr1 = ''.join(regions.get('CDR1', []))
        fr2 = ''.join(regions.get('FR2', []))
        cdr2 = ''.join(regions.get('CDR2', []))
        fr3 = ''.join(regions.get('FR3', []))
        cdr3 = ''.join(regions.get('CDR3', []))
        fr4 = ''.join(regions.get('FR4', []))
        
        # Create IMGTNumberedSequence for full access
        imgt_numbered = IMGTNumberedSequence(
            positions=positions,
            sequence=trimmed_seq,
            chain_type=chain_type,
            regions={
                'FR1': fr1, 'CDR1': cdr1, 'FR2': fr2, 'CDR2': cdr2,
                'FR3': fr3, 'CDR3': cdr3, 'FR4': fr4
            }
        )
        
        return CDRSet(
            cdr1=cdr1,
            cdr2=cdr2,
            cdr3=cdr3,
            fr1=fr1,
            fr2=fr2,
            fr3=fr3,
            fr4=fr4,
            imgt_numbered=imgt_numbered
        )
        
    except Exception as e:
        print(f"CDR extraction error: {e}")
        import traceback
        traceback.print_exc()
        return None

def number_scaffold_sequence(fr1: str, cdr1: str, fr2: str, cdr2: str, 
                              fr3: str, cdr3: str, fr4: str) -> Dict[Any, str]:
    """
    Create IMGT position mapping for a constructed sequence.
    
    This assumes standard IMGT lengths without insertions.
    For sequences with insertions, use extract_cdrs_fixed() instead.
    """
    positions = {}
    
    # FR1: positions 1-26
    for i, aa in enumerate(fr1):
        positions[1 + i] = aa
    
    # CDR1: positions 27-38 (can have insertions, but we map sequentially)
    for i, aa in enumerate(cdr1):
        pos = 27 + i
        if pos <= 38:
            positions[pos] = aa
        else:
            # CDR1 insertions (rare, but handle gracefully)
            positions[f"38{chr(65 + (pos - 38))}"] = aa  # 38A, 38B, etc.
    
    # FR2: positions 39-55
    # Standard FR2 without insertions has positions: 39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55
    # But IMGT FR2 often lacks position 41, 47 depending on numbering
    # We'll map sequentially to the defined positions
    fr2_positions = [39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55]
    for i, aa in enumerate(fr2):
        if i < len(fr2_positions):
            positions[fr2_positions[i]] = aa
    
    # CDR2: positions 56-65
    for i, aa in enumerate(cdr2):
        pos = 56 + i
        if pos <= 65:
            positions[pos] = aa
    
    # FR3: positions 66-104
    for i, aa in enumerate(fr3):
        pos = 66 + i
        if pos <= 104:
            positions[pos] = aa
    
    # CDR3: positions 105-117 (often has insertions)
    for i, aa in enumerate(cdr3):
        pos = 105 + i
        if pos <= 117:
            positions[pos] = aa
        else:
            # CDR3 insertions
            positions[f"117{chr(65 + (pos - 117))}"] = aa
    
    # FR4: positions 118-128
    for i, aa in enumerate(fr4):
        pos = 118 + i
        if pos <= 128:
            positions[pos] = aa
    
    return positions

# ============================================================
# FAMILY CLASSIFICATION (FIXED - uses IMGT positions)
# ============================================================

def classify_family_from_positions(positions: Dict[Any, str], full_sequence: str = None) -> str:
    """
    Classify VHH family using IMGT position dict (insertion-safe).
    
    C2/C4 classification is based on IMGT positions 55 and 100 (the extra disulfide pair),
    NOT total cysteine count. This ensures consistency with cysteine validation.
    
    Classification logic:
      - VH_like: pos50 = L (human-like)
      - VHH_W52: pos52 = W (W52 subtype)
      - F_C4/Y_C4: pos42 = F/Y AND C at BOTH pos55 AND pos100
      - F_C2/Y_C2: pos42 = F/Y AND NOT C4
      - Other_VHH: everything else
    """
    # Get hallmark positions directly from dict
    p42 = positions.get(42, positions.get('42', '-'))
    p49 = positions.get(49, positions.get('49', '-'))
    p50 = positions.get(50, positions.get('50', '-'))
    p52 = positions.get(52, positions.get('52', '-'))
    
    # VH-like check (human VH has L at position 50)
    if p50 == 'L':
        return "VH_like"
    
    # W52 check (VHH_W52 subtype)
    if p52 == 'W':
        return "VHH_W52"
    
    # Determine C2 vs C4 based on IMGT positions 55 and 100 (NOT total count!)
    # C4 requires BOTH positions to have cysteine (the extra disulfide pair)
    c55 = positions.get(55, positions.get('55', '')) == 'C'
    c100 = positions.get(100, positions.get('100', '')) == 'C'
    
    if c55 and c100:
        c_label = "C4"
    else:
        c_label = "C2"  # Includes INVALID cases (one but not both) - treat as C2 for classification
    
    # F or Y family based on position 42
    if p42 == 'F':
        return f"F_{c_label}"
    elif p42 == 'Y':
        return f"Y_{c_label}"
    
    return "Other_VHH"

# ============================================================
# CYSTEINE VALIDATION (v7.16)
# ============================================================

def validate_cysteine_count(sequence: str) -> Tuple[bool, str]:
    """
    Validate that sequence has even number of cysteines.
    Disulfides require pairs - odd count is always invalid.
    
    Returns:
        (is_valid, reason)
    """
    n_cys = sequence.count('C')
    if n_cys % 2 == 1:
        return False, f"Odd cysteine count ({n_cys}) - unpaired Cys causes aggregation"
    return True, "OK"


def validate_c4_cysteines(positions: Dict[Any, str], target_family: str) -> Tuple[bool, str]:
    """
    Validate C4 family candidates have BOTH required cysteines.
    
    C4 architecture requires:
      - Canonical: IMGT 23 + IMGT 104 (standard disulfide)
      - Extra pair: IMGT 55 + IMGT 100 (CDR1-CDR3 bridge)
    
    If targeting C4, BOTH 55 and 100 must be C.
    If only one is C, that's an unpaired cysteine -> invalid.
    
    Returns:
        (is_valid, reason)
    """
    if not target_family.endswith('_C4'):
        return True, "Not C4 target"
    
    c55 = positions.get(55, positions.get('55', '')) == 'C'
    c100 = positions.get(100, positions.get('100', '')) == 'C'
    
    if c55 and c100:
        return True, "Valid C4 (has both C55 and C100)"
    elif c55 and not c100:
        return False, "Invalid C4: has C55 but missing C100 (unpaired)"
    elif not c55 and c100:
        return False, "Invalid C4: has C100 but missing C55 (unpaired)"
    else:
        return False, "Invalid C4: missing both C55 and C100"


def detect_cysteine_class(positions: Dict[Any, str], sequence: str) -> str:
    """
    Determine C2 vs C4 based on actual cysteine positions at IMGT 55 and 100.
    
    Returns: 'C2', 'C4', or 'INVALID' (one but not both extra Cys)
    """
    c55 = positions.get(55, positions.get('55', '')) == 'C'
    c100 = positions.get(100, positions.get('100', '')) == 'C'
    
    if c55 and c100:
        return 'C4'
    elif not c55 and not c100:
        return 'C2'
    else:
        return 'INVALID'  # One but not both = unpaired


def filter_invalid_cysteines(candidates: List, verbose: bool = True) -> Tuple[List, Dict]:
    """
    Filter candidates with invalid cysteine patterns.
    
    Applied BEFORE deduplication per v7.16 design.
    
    Filters:
      1. Odd total cysteine count (always invalid)
      2. C4 targets missing the extra cysteine pair (55 + 100)
    
    Returns:
        (valid_candidates, filter_stats)
    """
    valid = []
    rejected_odd_cys = []
    rejected_invalid_c4 = []
    
    for c in candidates:
        # Skip lead - always keep
        if c.is_lead:
            valid.append(c)
            continue
        
        # Check 1: Odd cysteine count
        is_valid_count, reason_count = validate_cysteine_count(c.sequence)
        if not is_valid_count:
            rejected_odd_cys.append({
                'id': c.id,
                'reason': reason_count,
                'n_cys': c.sequence.count('C'),
                'track': c.design_track,
                'target_family': c.target_family
            })
            continue
        
        # Check 2: C4 enforcement
        is_valid_c4, reason_c4 = validate_c4_cysteines(c.imgt_positions, c.target_family)
        if not is_valid_c4:
            rejected_invalid_c4.append({
                'id': c.id,
                'reason': reason_c4,
                'target_family': c.target_family,
                'track': c.design_track,
                'c55': c.imgt_positions.get(55, '?'),
                'c100': c.imgt_positions.get(100, '?')
            })
            continue
        
        valid.append(c)
    
    stats = {
        'n_before': len(candidates),
        'n_after': len(valid),
        'n_rejected_odd_cys': len(rejected_odd_cys),
        'n_rejected_invalid_c4': len(rejected_invalid_c4),
        'rejected_odd_cys': rejected_odd_cys[:20],  # Cap for JSON size
        'rejected_invalid_c4': rejected_invalid_c4[:20]
    }
    
    if verbose:
        print(f"\n  Cysteine validation (v7.16)...")
        print(f"    Before: {stats['n_before']}, After: {stats['n_after']}")
        print(f"    Rejected odd-Cys: {stats['n_rejected_odd_cys']}")
        print(f"    Rejected invalid-C4: {stats['n_rejected_invalid_c4']}")
    
    return valid, stats

def get_hallmarks_from_positions(positions: Dict[Any, str]) -> str:
    """Get hallmark pattern from IMGT position dict."""
    p42 = positions.get(42, positions.get('42', '?'))
    p49 = positions.get(49, positions.get('49', '?'))
    p50 = positions.get(50, positions.get('50', '?'))
    p52 = positions.get(52, positions.get('52', '?'))
    return f"{p42}{p49}{p50}{p52}"

def get_cdr3_features(cdr3: str) -> Dict[str, str]:
    """Extract CDR3 features for rule matching."""
    length = len(cdr3)
    
    if length <= 8:
        cdr3_len = "short"
    elif length <= 14:
        cdr3_len = "medium"
    else:
        cdr3_len = "long"
    
    pos_charge = sum(1 for aa in cdr3 if aa in "RKH")
    neg_charge = sum(1 for aa in cdr3 if aa in "DE")
    net_charge = pos_charge - neg_charge
    
    if net_charge >= 2:
        cdr3_charge = "positive"
    elif net_charge <= -2:
        cdr3_charge = "negative"
    else:
        cdr3_charge = "neutral"
    
    n_cys = cdr3.count("C")
    if n_cys == 0:
        cdr3_cys = "none"
    elif n_cys == 1:
        cdr3_cys = "one"
    else:
        cdr3_cys = "multiple"
    
    cdr3_end = cdr3[-3:] if len(cdr3) >= 3 else cdr3
    
    return {
        "cdr3_len": cdr3_len,
        "cdr3_charge": cdr3_charge,
        "cdr3_cys": cdr3_cys,
        "cdr3_end": cdr3_end,
        "cdr3[-1]": cdr3[-1] if cdr3 else "",
        "cdr3[-2]": cdr3[-2] if len(cdr3) >= 2 else "",
        "cdr3[-3]": cdr3[-3] if len(cdr3) >= 3 else "",
    }

def calculate_framework_identity(cand_positions: Dict[Any, str], 
                                  orig_positions: Dict[Any, str]) -> float:
    """Calculate % identity between candidate and original framework positions."""
    # Only compare framework positions (not CDRs)
    fw_ranges = [(1, 26), (39, 55), (66, 104), (118, 128)]
    
    matches = 0
    total = 0
    
    for start, end in fw_ranges:
        for pos in range(start, end + 1):
            orig_aa = orig_positions.get(pos, orig_positions.get(str(pos), ''))
            cand_aa = cand_positions.get(pos, cand_positions.get(str(pos), ''))
            
            if orig_aa:
                total += 1
                if orig_aa == cand_aa:
                    matches += 1
    
    return (matches / total * 100) if total > 0 else 100.0

# ============================================================
# ESM2 LANGUAGE MODEL SCORER
# ============================================================

class ESM2Scorer:
    """Score sequences using ESM2 pseudo-perplexity."""
    
    def __init__(self, model_name: str = "facebook/esm2_t6_8M_UR50D"):
        self.model_name = model_name
        self.model = None
        self.tokenizer = None
        self.device = None
        
    def load(self):
        if self.model is not None:
            return
        
        if not ESM2_AVAILABLE or not TORCH_AVAILABLE:
            print("  ESM2 not available, skipping...")
            return
        
        print(f"  Loading ESM2 model: {self.model_name}")
        self.tokenizer = AutoTokenizer.from_pretrained(self.model_name)
        self.model = AutoModelForMaskedLM.from_pretrained(self.model_name)
        
        self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        self.model = self.model.to(self.device)
        self.model.eval()
        
        print(f"    Device: {self.device}")
    
    def score_batch(self, sequences: List[str], batch_size: int = 32) -> List[Tuple[float, float]]:
        if not ESM2_AVAILABLE or not TORCH_AVAILABLE:
            return [(0.0, 1.0) for _ in sequences]
        
        self.load()
        if self.model is None:
            return [(0.0, 1.0) for _ in sequences]
        
        results = []
        
        iterator = range(0, len(sequences), batch_size)
        if TQDM_AVAILABLE:
            iterator = tqdm(list(iterator), desc="ESM2 scoring", leave=False)
        
        for i in iterator:
            batch = sequences[i:i+batch_size]
            
            for seq in batch:
                try:
                    seq_inputs = self.tokenizer(seq, return_tensors="pt").to(self.device)
                    with torch.no_grad():
                        outputs = self.model(**seq_inputs, labels=seq_inputs["input_ids"])
                    loss = outputs.loss.item()
                    perplexity = np.exp(loss)
                    results.append((loss, perplexity))
                except Exception as e:
                    print(f"    ESM2 error: {e}")
                    results.append((0.0, 1.0))
        
        return results

# ============================================================
# ESMFOLD STRUCTURE PREDICTOR
# ============================================================

class ESMFoldPredictor:
    """Predict structure and pLDDT scores using ESMFold."""
    
    def __init__(self):
        self.model = None
        self.device = None
    
    def load(self):
        if self.model is not None:
            return True
        
        if not ESMFOLD_AVAILABLE or not TORCH_AVAILABLE:
            print("  ESMFold not available")
            return False
        
        print("  Loading ESMFold model...")
        try:
            self.model = esm.pretrained.esmfold_v1()
            self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
            self.model = self.model.to(self.device)
            self.model.eval()
            self.model.set_chunk_size(128)
            print(f"    Device: {self.device}")
            return True
        except Exception as e:
            print(f"    Failed to load ESMFold: {e}")
            return False
    
    def predict_plddt(self, sequence: str) -> Dict[str, float]:
        if not self.load():
            return {'plddt_mean': 0.0, 'plddt_median': 0.0, 'plddt_min': 0.0, 'plddt_scores': []}
        
        try:
            with torch.no_grad():
                output = self.model.infer_pdb(sequence)
            
            # Parse pLDDT from B-factor column
            residue_plddts = []
            current_residue = []
            current_resnum = None
            
            for line in output.split('\n'):
                if line.startswith('ATOM'):
                    try:
                        resnum = int(line[22:26].strip())
                        bfactor = float(line[60:66].strip())
                        
                        if current_resnum is None:
                            current_resnum = resnum
                        
                        if resnum != current_resnum:
                            if current_residue:
                                residue_plddts.append(np.mean(current_residue))
                            current_residue = [bfactor]
                            current_resnum = resnum
                        else:
                            current_residue.append(bfactor)
                    except:
                        continue
            
            if current_residue:
                residue_plddts.append(np.mean(current_residue))
            
            return {
                'plddt_mean': np.mean(residue_plddts) if residue_plddts else 0.0,
                'plddt_median': np.median(residue_plddts) if residue_plddts else 0.0,
                'plddt_min': np.min(residue_plddts) if residue_plddts else 0.0,
                'plddt_scores': residue_plddts
            }
            
        except Exception as e:
            print(f"    ESMFold error: {e}")
            return {'plddt_mean': 0.0, 'plddt_median': 0.0, 'plddt_min': 0.0, 'plddt_scores': []}
    
    def predict_batch(self, sequences: List[str], 
                      region_ranges: List[Dict[str, Tuple[int, int]]] = None) -> List[Dict[str, float]]:
        results = []
        
        iterator = enumerate(sequences)
        if TQDM_AVAILABLE:
            iterator = tqdm(list(iterator), desc="ESMFold prediction", leave=False)
        
        for i, seq in iterator:
            result = self.predict_plddt(seq)
            
            if region_ranges and i < len(region_ranges) and result['plddt_scores']:
                ranges = region_ranges[i]
                plddt_scores = result['plddt_scores']
                
                for region in ['cdr1', 'cdr2', 'cdr3']:
                    if region in ranges:
                        start, end = ranges[region]
                        end = min(end, len(plddt_scores))
                        if start < end:
                            result[f'plddt_{region}'] = np.mean(plddt_scores[start:end])
                        else:
                            result[f'plddt_{region}'] = 0.0
                    else:
                        result[f'plddt_{region}'] = 0.0
                
                if 'framework_ranges' in ranges:
                    fw_scores = []
                    for start, end in ranges['framework_ranges']:
                        end = min(end, len(plddt_scores))
                        if start < end:
                            fw_scores.extend(plddt_scores[start:end])
                    result['plddt_framework'] = np.mean(fw_scores) if fw_scores else 0.0
            
            results.append(result)
        
        return results

# ============================================================
# MULTI-FAMILY PROBABILISTIC SCORER
# ============================================================

class MultiFamilyProbabilisticScorer:
    """Score sequences using multi-family probabilistic rule compliance."""
    
    def __init__(self, rules_file: str, archetypes_file: str, hallmark_db_path: str = None):
        print("  Loading rules engine...")
        with open(rules_file, 'r') as f:
            self.rules = json.load(f)
        
        print("  Loading archetypes...")
        with open(archetypes_file, 'r') as f:
            archetypes = json.load(f)
        
        self.rules_by_family = defaultdict(list)
        for rule in self.rules:
            family = rule.get('family', '')
            if family:
                self.rules_by_family[family].append(rule)
        
        self.archetypes = {a['family']: a for a in archetypes}

        # Optional hallmark/subfamily DB (recommended for true-VHH mode)
        self.hallmark_db = HallmarkDB(hallmark_db_path) if hallmark_db_path else None

        self.families = list(self.archetypes.keys())
        
        print(f"    {len(self.rules)} rules across {len(self.families)} families")
    
    def count_vernier_matches(self, positions: Dict[Any, str], family: str) -> Tuple[int, int]:
        """Count vernier position matches against family archetype."""
        archetype = self.archetypes.get(family, {})
        arch_positions = archetype.get('positions', {})
        
        matches = 0
        total = 0
        
        for pos_str, data in arch_positions.items():
            try:
                imgt_num = int(pos_str.replace('IMGT', ''))
                if imgt_num not in ALL_VERNIER_POSITIONS:
                    continue
                
                consensus_aa = data.get('consensus', '')
                actual_aa = positions.get(imgt_num, positions.get(str(imgt_num), ''))
                
                if consensus_aa:
                    total += 1
                    if actual_aa == consensus_aa:
                        matches += 1
            except ValueError:
                continue
        
        return matches, total

    def count_core_vernier_matches(self, positions: Dict[Any, str], family: str) -> Tuple[int, int]:
        """Count matches only for the *core* vernier set (excludes hallmarks)."""
        archetype = self.archetypes.get(family, {})
        arch_positions = archetype.get('positions', {})

        matches = 0
        total = 0
        for pos_str, data in arch_positions.items():
            try:
                imgt_num = int(pos_str.replace('IMGT', ''))
            except ValueError:
                continue
            if imgt_num not in VERNIER_POSITIONS_FR3:
                continue
            consensus_aa = data.get('consensus', '')
            if not consensus_aa:
                continue
            actual_aa = positions.get(imgt_num, positions.get(str(imgt_num), ''))
            total += 1
            if actual_aa == consensus_aa:
                matches += 1
        return matches, total
    
    def compute_family_probability(self, positions: Dict[Any, str], 
                                    cdr3_features: Dict[str, str]) -> Dict[str, float]:
        """Compute P(family | sequence) using position dict."""
        scores = {}
        
        for family, archetype in self.archetypes.items():
            if family == 'VH_like':
                continue
            
            score = 0.0
            total_weight = 0.0
            
            # Vernier matches
            arch_positions = archetype.get('positions', {})
            for pos_str, data in arch_positions.items():
                try:
                    imgt_num = int(pos_str.replace('IMGT', ''))
                    consensus_aa = data.get('consensus', '')
                    actual_aa = positions.get(imgt_num, positions.get(str(imgt_num), ''))
                    
                    if actual_aa and consensus_aa:
                        weight = 2.0 if imgt_num in VERNIER_POSITIONS_FR3 else 1.0
                        total_weight += weight
                        if actual_aa == consensus_aa:
                            score += weight
                except ValueError:
                    continue
            
            # CDR3 compatibility
            cdr3_len = cdr3_features.get('cdr3_len', 'medium')
            total_weight += 1.0
            if family.startswith('F_'):
                score += 1.0 if cdr3_len in ['short', 'medium'] else 0.5
            elif family.startswith('Y_'):
                score += 1.0 if cdr3_len in ['medium', 'long'] else 0.5
            else:
                score += 0.7
            
            # Hallmark compatibility
            hallmarks = get_hallmarks_from_positions(positions)
            total_weight += 2.0
            if family.startswith('F_') and hallmarks[0] == 'F':
                score += 2.0
            elif family.startswith('Y_') and hallmarks[0] == 'Y':
                score += 2.0
            elif family == 'Other_VHH':
                score += 1.0
            
            scores[family] = score / total_weight if total_weight > 0 else 0.0
        
        # Normalize
        total = sum(scores.values())
        if total > 0:
            scores = {k: v / total for k, v in scores.items()}
        
        return scores
    
    def compute_rule_compliance(self, positions: Dict[Any, str], family: str, 
                                 cdr3_features: Dict[str, str],
                                 min_confidence: float = 0.5) -> Tuple[float, int, int, int, List[str]]:
        """Compute rule compliance using position dict."""
        rules = self.rules_by_family.get(family, [])
        if not rules:
            return 1.0, 0, 0, 0, []
        
        hallmarks = get_hallmarks_from_positions(positions)
        
        satisfied_weight = 0.0
        total_weight = 0.0
        rules_passed = 0
        rules_applicable = 0
        violations = []
        
        for rule in rules:
            confidence = rule.get('confidence', 0)
            if confidence < min_confidence:
                continue
            
            condition = rule.get('condition', '')
            if not self._condition_applies(condition, hallmarks, cdr3_features):
                continue
            
            rules_applicable += 1
            
            position_num = rule.get('position_num', 0)
            suggested_aa = rule.get('suggested_aa', '')
            
            if position_num and suggested_aa:
                actual_aa = positions.get(position_num, positions.get(str(position_num), ''))
                
                if actual_aa == suggested_aa:
                    satisfied_weight += confidence
                    rules_passed += 1
                else:
                    violations.append(f"IMGT{position_num}: expected {suggested_aa}, got {actual_aa} (conf={confidence:.2f})")
                
                total_weight += confidence
        
        compliance = satisfied_weight / total_weight if total_weight > 0 else 1.0
        return compliance, rules_passed, len(rules), rules_applicable, violations
    
    def score(self, candidate: VHHCandidate, cdr3_features: Dict[str, str],
              family_threshold: float = 0.15, rule_min_confidence: float = 0.5) -> ScoringResult:
        """Score a candidate using position dict."""
        result = ScoringResult()
        positions = candidate.imgt_positions
        
        # Family probabilities
        family_probs = self.compute_family_probability(positions, cdr3_features)
        result.family_probabilities = family_probs
        
        # Vernier matches for top family
        top_family = max(family_probs.items(), key=lambda x: x[1])[0] if family_probs else 'F_C2'
        result.vernier_matches, result.vernier_total = self.count_vernier_matches(positions, top_family)
        
        # Rule compliance
        total_rules_passed = 0
        total_rules = 0
        total_applicable = 0
        all_violations = []
        
        for family, prob in family_probs.items():
            if prob < family_threshold:
                continue
            
            compliance, passed, total, applicable, violations = self.compute_rule_compliance(
                positions, family, cdr3_features, rule_min_confidence
            )
            
            result.rule_compliance[family] = compliance
            total_rules_passed += passed
            total_rules += total
            total_applicable += applicable
            all_violations.extend([f"[{family}] {v}" for v in violations[:3]])
        
        result.rules_passed = total_rules_passed
        result.rules_total = total_rules
        result.rules_applicable = total_applicable
        result.top_violations = all_violations[:5]
        
        # Weighted naturalness
        weighted_nat = 0.0
        for family, prob in family_probs.items():
            if prob >= family_threshold and family in result.rule_compliance:
                weighted_nat += prob * result.rule_compliance[family]
        
        result.weighted_naturalness = weighted_nat

        # ------------------------------------------------------------------
        # Confidence (0–100) for ranking/reporting
        #
        # This is not a calibrated probability. For this pipeline, CDR3 is
        # held constant per design batch, so “naturalness” doesn't meaningfully
        # separate candidates. Instead we summarize:
        #   • hallmark match (required VHH identity)
        #   • core vernier match to the chosen family archetype
        #   • weighted rule compliance
        fam = top_family
        expected_hallmarks_key = expected_hallmarks_for_family(fam)
        hallmark_match = compute_hallmark_match(candidate.imgt_positions, expected_hallmarks_key)['score']

        core_matches, core_total = self.count_core_vernier_matches(candidate.imgt_positions, fam)
        core_vernier_match = (core_matches / core_total) if core_total else 0.0

        # Weights chosen so that “identity” (hallmarks+vernier) and “contextual
        # tuning” (rules) both matter.
        conf = 0.25 * hallmark_match + 0.25 * core_vernier_match + 0.50 * float(weighted_nat)
        result.confidence_pct = max(0.0, min(100.0, 100.0 * conf))
        result.confidence_components = {
            'hallmark_match': hallmark_match,
            'core_vernier_match': core_vernier_match,
            'weighted_rule_score': float(weighted_nat),
            'family': fam,
        }
        
        return result
    
    def _condition_applies(self, condition: str, hallmarks: str, cdr3_features: Dict[str, str]) -> bool:
        if not condition:
            return True
        
        parts = [p.strip() for p in condition.split(' AND ')]
        
        for part in parts:
            if '=' not in part:
                continue
            key, value = part.split('=', 1)
            key, value = key.strip(), value.strip()
            
            if key == 'hallmarks':
                if hallmarks != value:
                    return False
            elif key in cdr3_features:
                if cdr3_features[key] != value:
                    return False
        
        return True

# ============================================================
# COMBINED SCORER
# ============================================================

class CombinedScorer:
    """Combined scoring using ESM2 + ESMFold + Multi-family rule compliance."""
    
    def __init__(self, rules_file: str, archetypes_file: str, 
                 esm_model: str = "facebook/esm2_t6_8M_UR50D",
                 use_esm2: bool = True,
                 use_esmfold: bool = True,
                 esm2_weight: float = 0.2,
                 esmfold_weight: float = 0.0,
                 rule_weight: float = 0.7):
        
        self.multi_family_scorer = MultiFamilyProbabilisticScorer(rules_file, archetypes_file)
        
        self.use_esm2 = use_esm2 and ESM2_AVAILABLE and TORCH_AVAILABLE
        self.use_esmfold = use_esmfold and ESMFOLD_AVAILABLE and TORCH_AVAILABLE
        
        # Normalize weights based on ACTUAL availability (not just user request)
        # Fix: Use self.use_esm2/self.use_esmfold to reflect runtime availability
        total_weight = rule_weight
        if self.use_esm2:
            total_weight += esm2_weight
        if self.use_esmfold:
            total_weight += esmfold_weight

        self.esm2_weight = (esm2_weight / total_weight) if self.use_esm2 and total_weight > 0 else 0.0
        self.esmfold_weight = (esmfold_weight / total_weight) if self.use_esmfold and total_weight > 0 else 0.0
        self.rule_weight = rule_weight / total_weight if total_weight > 0 else 1.0      

        if self.use_esm2:
            self.esm2_scorer = ESM2Scorer(esm_model)
        else:
            self.esm2_scorer = None
        
        if self.use_esmfold:
            self.esmfold_predictor = ESMFoldPredictor()
        else:
            self.esmfold_predictor = None
        
        print(f"  Weights: ESM2={self.esm2_weight:.2f}, ESMFold={self.esmfold_weight:.2f}, Rules={self.rule_weight:.2f}")
    
    def score_candidates(self, candidates: List[VHHCandidate], cdr3_features: Dict[str, str],
                         original_positions: Dict[Any, str] = None) -> List[VHHCandidate]:
        """Score all candidates."""
        print(f"\nScoring {len(candidates)} candidates...")
        
        # Framework identity
        if original_positions:
            for candidate in candidates:
                candidate.framework_identity_pct = calculate_framework_identity(
                    candidate.imgt_positions, original_positions
                )
        
        # Rule compliance
        print("  Computing rule compliance...")
        for candidate in tqdm(candidates, desc="Rule compliance", disable=not TQDM_AVAILABLE):
            result = self.multi_family_scorer.score(candidate, cdr3_features)
            candidate.scoring = result
        
        # ESM2
        if self.use_esm2:
            print("  Computing ESM2 scores...")
            sequences = [c.sequence for c in candidates]
            esm2_scores = self.esm2_scorer.score_batch(sequences)
            
            for candidate, (loss, perplexity) in zip(candidates, esm2_scores):
                candidate.scoring.esm2_loss = loss
                candidate.scoring.esm2_perplexity = perplexity
        
        # ESMFold
        if self.use_esmfold:
            print("  Computing ESMFold pLDDT...")
            sequences = [c.sequence for c in candidates]
            
            region_ranges = []
            for c in candidates:
                fr1_len = len(c.fr1)
                cdr1_end = fr1_len + len(c.cdr1)
                fr2_end = cdr1_end + len(c.fr2)
                cdr2_end = fr2_end + len(c.cdr2)
                fr3_end = cdr2_end + len(c.fr3)
                cdr3_end = fr3_end + len(c.cdr3)
                
                region_ranges.append({
                    'cdr1': (fr1_len, cdr1_end),
                    'cdr2': (fr2_end, cdr2_end),
                    'cdr3': (fr3_end, cdr3_end),
                    'framework_ranges': [(0, fr1_len), (cdr1_end, fr2_end), (cdr2_end, fr3_end), (cdr3_end, cdr3_end + len(c.fr4))]
                })
            
            plddt_results = self.esmfold_predictor.predict_batch(sequences, region_ranges)
            
            for candidate, plddt in zip(candidates, plddt_results):
                candidate.scoring.plddt_mean = plddt.get('plddt_mean', 0.0)
                candidate.scoring.plddt_median = plddt.get('plddt_median', 0.0)
                candidate.scoring.plddt_min = plddt.get('plddt_min', 0.0)
                candidate.scoring.plddt_cdr1 = plddt.get('plddt_cdr1', 0.0)
                candidate.scoring.plddt_cdr2 = plddt.get('plddt_cdr2', 0.0)
                candidate.scoring.plddt_cdr3 = plddt.get('plddt_cdr3', 0.0)
                candidate.scoring.plddt_framework = plddt.get('plddt_framework', 0.0)
        
        # Combined scores
        print("  Computing combined scores...")
        
        if self.use_esm2:
            esm2_losses = [c.scoring.esm2_loss for c in candidates if c.scoring.esm2_loss > 0]
            if esm2_losses:
                esm2_min, esm2_max = min(esm2_losses), max(esm2_losses)
                esm2_range = esm2_max - esm2_min if esm2_max > esm2_min else 1.0
            else:
                esm2_min, esm2_range = 0, 1
        
        for candidate in candidates:
            if candidate.is_lead:
                candidate.scoring.combined_score = 999.0
                continue
            
            rule_score = candidate.scoring.weighted_naturalness
            
            if self.use_esm2 and candidate.scoring.esm2_loss > 0:
                esm2_normalized = 1.0 - (candidate.scoring.esm2_loss - esm2_min) / esm2_range
            else:
                esm2_normalized = 0.5
            
            if self.use_esmfold and candidate.scoring.plddt_mean > 0:
                plddt_normalized = candidate.scoring.plddt_mean / 100.0
            else:
                plddt_normalized = 0.5
            
            combined = (
                self.esm2_weight * esm2_normalized +
                self.esmfold_weight * plddt_normalized +
                self.rule_weight * rule_score
            )
            
            candidate.scoring.combined_score = combined
        
        return candidates

# ============================================================
# CANDIDATE GENERATOR (FIXED)
# ============================================================

class CandidateGenerator:
    """Generate VHH candidates using IMGT position-based mutations."""
    
    def __init__(self, rules_file: str, archetypes_file: str, hallmark_db_path: str = None):
        print("\nInitializing candidate generator...")
        
        with open(rules_file, 'r') as f:
            rules = json.load(f)
        
        with open(archetypes_file, 'r') as f:
            archetypes = json.load(f)
        
        self.compensation_rules = [r for r in rules if r.get('rule_type') == 'cdr_to_imgt']
        self.triplet_rules = [r for r in rules if r.get('rule_type') == 'triplet']
        
        self.rules_by_hallmarks = defaultdict(list)
        self.rules_by_condition = defaultdict(list)
        
        for rule in self.compensation_rules:
            cond = rule.get('condition', '')
            family = rule.get('family', '')
            
            if cond.startswith('hallmarks='):
                hallmarks = cond.split('=')[1]
                self.rules_by_hallmarks[(hallmarks, family)].append(rule)
            else:
                self.rules_by_condition[(cond, family)].append(rule)
        
        self.archetypes = {a['family']: a for a in archetypes}

        # Optional hallmark/subfamily DB (recommended for true-VHH mode)
        self.hallmark_db = HallmarkDB(hallmark_db_path) if hallmark_db_path else None

        
        # Build scaffold position dict
        self.scaffold_positions = number_scaffold_sequence(
            UNIVERSAL_SCAFFOLD['FR1'], '',  # No CDR in scaffold
            UNIVERSAL_SCAFFOLD['FR2'], '',
            UNIVERSAL_SCAFFOLD['FR3'], '',
            UNIVERSAL_SCAFFOLD['FR4']
        )
        
        print(f"  {len(self.compensation_rules)} compensation rules")
        print(f"  {len(self.triplet_rules)} triplet rules")
        print(f"  {len(self.archetypes)} family archetypes")
    
    def generate(self, cdrs: CDRSet, n_candidates: int, 
                 target_hallmarks: str = 'FERG',
                 target_families: List[str] = None) -> List[VHHCandidate]:
        """Generate N candidates with proper IMGT position tracking.
        
        v8.0: Always generates BOTH original framework and Universal scaffold variants.
        """
        hallmark_profile, epi_pairs = self._get_hallmark_profile(target_hallmarks)
        rng = random.Random()

        candidates = []
        generation_order = 0
        
        if not target_families:
            # True-VHH defaults: do not include VH_like or Other_VHH unless explicitly requested
            target_families = ['F_C4', 'F_C2', 'Y_C4', 'Y_C2']
        
        cdr3_features = get_cdr3_features(cdrs.cdr3)
        
        # Get original IMGT positions
        if cdrs.imgt_numbered:
            original_positions = cdrs.imgt_numbered.positions.copy()
        else:
            original_positions = number_scaffold_sequence(
                cdrs.fr1, cdrs.cdr1, cdrs.fr2, cdrs.cdr2,
                cdrs.fr3, cdrs.cdr3, cdrs.fr4
            )
        
        print(f"\nGenerating {n_candidates} candidates...")
        print(f"  Target families: {target_families}")
        print(f"  Input hallmarks: {get_hallmarks_from_positions(original_positions)}")
        print(f"  Mode: BOTH (original framework + Universal scaffold)")

        # AUTO hallmark selection (subfamily-level)
        hallmark_pool = None
        if isinstance(target_hallmarks, str) and target_hallmarks.upper() == 'AUTO':
            if not self.hallmark_db:
                raise ValueError("target_hallmarks='AUTO' requires --hallmark-db pointing to comprehensive_subfamily_analysis_imgt.xlsx")
            hallmark_pool = self.hallmark_db.pick_reference_pool(k_per_cluster=2, cluster_level='cluster_10', min_n=50000)
            if not hallmark_pool:
                raise ValueError('No eligible true-VHH hallmarks found in hallmark DB after filtering.')
            # Always include YQRL if not present
            if 'YQRL' not in hallmark_pool:
                hallmark_pool.insert(0, 'YQRL')
            print(f"  AUTO hallmark pool (true VHH): {len(hallmark_pool)} hallmarks")
            print(f"    Available: {hallmark_pool}")
            # Check what family types are covered
            f_hallmarks = [h for h in hallmark_pool if h.startswith('F')]
            y_hallmarks = [h for h in hallmark_pool if h.startswith('Y')]
            if not y_hallmarks:
                print(f"    Note: No Y-type hallmarks in pool")
                print(f"    Y_* families will use F-type hallmarks")

        
        # LEAD: Original input (untouched)
        generation_order += 1
        input_family = classify_family_from_positions(original_positions)
        input_seq = cdrs.fr1 + cdrs.cdr1 + cdrs.fr2 + cdrs.cdr2 + cdrs.fr3 + cdrs.cdr3 + cdrs.fr4
        
        lead = VHHCandidate(
            id="Lead_Input_Original",
            rank=0,
            sequence=input_seq,
            framework_source='input',
            family=input_family,
            cdr1=cdrs.cdr1, cdr2=cdrs.cdr2, cdr3=cdrs.cdr3,
            fr1=cdrs.fr1, fr2=cdrs.fr2, fr3=cdrs.fr3, fr4=cdrs.fr4,
            imgt_positions=original_positions.copy(),
            mutations=[],
            strategy="original_input",
            is_lead=True,
            generation_order=generation_order,
            construction_method="Original input sequence (untouched)",
            target_family=input_family,
            framework_identity_pct=100.0
        )
        candidates.append(lead)
        
        # =================================================================
        # UNIVERSAL TRACK 0 CONTROLS (v8.0)
        # =================================================================
        # Generate Universal scaffold controls ONCE (not per-family)
        # These test whether the CDRs work in a standardized framework
        # =================================================================
        
        print(f"\n  Generating Universal scaffold controls...")
        universal_controls, generation_order = self._generate_universal_controls(
            cdrs, target_families[0] if target_families else 'F_C2',
            hallmark_pool or [target_hallmarks], generation_order
        )
        candidates.extend(universal_controls)
        print(f"    Generated {len(universal_controls)} Universal controls")
        
        # Calculate candidates per family (remaining after lead + Universal controls)
        n_remaining = n_candidates - 1 - len(universal_controls)
        n_per_family = n_remaining // len(target_families)
        
        # Collect all hallmark-family combinations to generate
        generation_plan = []  # List of (family, hallmark, n_candidates)
        
        for i, family in enumerate(target_families):
            n_for_family = n_per_family
            if i == len(target_families) - 1:
                n_for_family = n_remaining - sum(g[2] for g in generation_plan)
            
            if n_for_family <= 0:
                continue
            
            # Get compatible hallmarks for this family
            if hallmark_pool:
                fam_type = family[0] if family else 'F'
                fam_cys = family.split('_')[1] if '_' in family else 'C2'
                compatible = [hm for hm in hallmark_pool if hm.startswith(fam_type) and self.hallmark_db.hallmark_to_family(hm).endswith(fam_cys)]
                if not compatible:
                    compatible = [hm for hm in hallmark_pool if hm.startswith(fam_type)]
                if not compatible:
                    compatible = hallmark_pool
                    print(f"    Warning: No {fam_type}-type hallmarks in pool, using available: {compatible}")
                if not compatible:
                    print(f"    Skipping {family}: no hallmarks available")
                    continue
                
                # PRIORITY: Always include YQRL if compatible (most constrained scaffold)
                if 'YQRL' in compatible:
                    # Move YQRL to front of list
                    compatible = ['YQRL'] + [h for h in compatible if h != 'YQRL']
                
                # Distribute candidates across ALL compatible hallmarks
                n_per_hallmark = n_for_family // len(compatible)
                remainder = n_for_family % len(compatible)
                
                for j, hm in enumerate(compatible):
                    n_for_hm = n_per_hallmark + (1 if j < remainder else 0)
                    if n_for_hm > 0:
                        generation_plan.append((family, hm, n_for_hm))
            else:
                fam_hallmarks = 'YERL' if family.startswith('Y_') and target_hallmarks == 'FERG' else target_hallmarks
                generation_plan.append((family, fam_hallmarks, n_for_family))
        
        # Print generation plan
        print(f"  Generation plan ({len(generation_plan)} hallmark-family combos):")
        for family, hm, n in generation_plan:
            print(f"    {family} + {hm}: {n} candidates")
        
        # Track which hallmarks have already had controls generated (for deduplication)
        # Controls are hallmark-specific, not family-specific, so we only generate them once
        hallmarks_with_controls = set()
        
        # Execute generation plan
        for family, fam_hallmarks, n_for_combo in generation_plan:
            rules = self._get_rules(fam_hallmarks, family, cdr3_features)
            archetype = self._get_archetype(family)
            triplet_rules = self._get_triplet_rules(family, cdrs.cdr3)
            
            print(f"  {family}/{fam_hallmarks}: {n_for_combo} candidates, {len(rules)} rules")
            
            # v8.2: Equal allocation between original and universal
            # Both get 50% to ensure Universal is well-represented in Track 4
            n_orig = n_for_combo // 2
            n_univ = n_for_combo - n_orig
            
            # Generate Universal scaffold variants (including Track 4 candidates)
            if n_univ > 0:
                new_cands, generation_order = self._generate_universal(
                    cdrs, family, rules, archetype, triplet_rules, n_univ,
                    fam_hallmarks, generation_order
                )
                candidates.extend(new_cands)
            
            # Generate Original framework variants
            if n_orig > 0 and cdrs.fr1:
                new_cands, generation_order = self._generate_original(
                    cdrs, family, rules, archetype, triplet_rules, n_orig,
                    fam_hallmarks, generation_order, original_positions,
                    hallmarks_with_controls  # Pass the deduplication set
                )
                candidates.extend(new_cands)
        
        return candidates
    
    def _get_rules(self, hallmarks: str, family: str, cdr3_features: Dict[str, str],
                   min_support: int = 500, min_confidence: float = 0.6) -> List[dict]:
        applicable = []
        
        for rule in self.rules_by_hallmarks.get((hallmarks, family), []):
            if rule.get('support', 0) >= min_support and rule.get('confidence', 0) >= min_confidence:
                applicable.append(rule)
        
        for (cond, rule_family), rules in self.rules_by_condition.items():
            if rule_family != family:
                continue
            if self._check_condition(cond, cdr3_features):
                for rule in rules:
                    if rule.get('support', 0) >= min_support and rule.get('confidence', 0) >= min_confidence:
                        applicable.append(rule)
        
        applicable.sort(key=lambda r: (-r.get('confidence', 0), -r.get('support', 0)))
        
        by_position = {}
        for rule in applicable:
            pos = rule.get('position_num', 0)
            if pos not in by_position:
                by_position[pos] = rule
        
        return list(by_position.values())
    
    def _check_condition(self, condition: str, cdr3_features: Dict[str, str]) -> bool:
        if '=' not in condition:
            return False
        key, value = condition.split('=', 1)
        return cdr3_features.get(key.strip()) == value.strip()
    
    def _get_archetype(self, family: str) -> Dict[int, str]:
        arch = self.archetypes.get(family, {})
        positions = arch.get('positions', {})
        
        pattern = {}
        for pos_str, data in positions.items():
            try:
                imgt_num = int(pos_str.replace('IMGT', ''))
                pattern[imgt_num] = data.get('consensus', '')
            except ValueError:
                continue
        
        return pattern
    
    def _get_triplet_rules(self, family: str, cdr3: str, min_support: int = 100) -> List[dict]:
        cdr3_end = cdr3[-3:] if len(cdr3) >= 3 else cdr3
        applicable = []
        
        for rule in self.triplet_rules:
            if rule.get('family') != family:
                continue
            if rule.get('support', 0) < min_support:
                continue
            cond = rule.get('condition', '')
            if f"cdr3_end={cdr3_end}" in cond:
                applicable.append(rule)
        
        return sorted(applicable, key=lambda r: -r.get('confidence', 0))
    

    def _get_hallmark_profile(self, hallmarks: str) -> Tuple[dict, List[Tuple[int,int,List[Tuple[str,str,float]]]]]:
        """Hallmark-specific vernier profile + epistatic pairs (optional)."""
        if not self.hallmark_db:
            return {}, []
        return self.hallmark_db.vernier_profile(hallmarks), self.hallmark_db.top_epistatic_pairs(hallmarks, max_pairs=4)

    def _sample_vernier_aa(self, pos: int, profile: dict, rng: random.Random) -> str:
        """Sample AA at a vernier position from hallmark profile (cons + 2nd), fallback to consensus."""
        d = profile.get(pos)
        if not d:
            return ''
        cons, p = d.get('cons',''), float(d.get('freq',0.0))
        alt, q = d.get('alt',''), float(d.get('alt_freq',0.0))
        # Normalize among provided options (ignore '-' alts)
        opts=[]
        if cons and cons != '-':
            opts.append((cons, max(p, 1e-6)))
        if alt and alt != '-' and alt != cons:
            opts.append((alt, max(q, 1e-6)))
        if not opts:
            return cons if cons != '-' else ''
        tot=sum(w for _,w in opts)
        r=rng.random()*tot
        s=0.0
        for aa,w in opts:
            s+=w
            if r<=s:
                return aa
        return opts[0][0]

    def _apply_epistatic_pairs(self, positions: Dict[Any,str], epi_pairs, rng: random.Random, profile: dict) -> List[Dict[str, Any]]:
        """Apply a small number of coupled vernier pairs by sampling observed combinations.

        Returns a list of applied pair assignments for provenance.
        """
        applied: List[Dict[str, Any]] = []
        for (p1,p2, combs) in epi_pairs:
            # Only couple if both are present in framework
            if p1 not in positions or p2 not in positions:
                continue
            # sample combination
            tot=sum(max(freq,1e-8) for _,_,freq in combs)
            r=rng.random()*tot
            s=0.0
            aa1=aa2=None
            for a1,a2,freq in combs:
                s+=max(freq,1e-8)
                if r<=s:
                    aa1,aa2=a1,a2
                    break
            if aa1 and aa2:
                positions[p1]=aa1
                positions[p2]=aa2
                applied.append({'pos1': p1, 'pos2': p2, 'aa1': aa1, 'aa2': aa2})

        return applied


    def _generate_universal_controls(self, cdrs: CDRSet, family: str,
                                     hallmark_pool: List[str],
                                     generation_order: int) -> Tuple[List[VHHCandidate], int]:
        """
        Generate Universal scaffold Track 0 controls (v8.0).
        
        Creates:
        a) Pure CDR graft into Universal (no hallmark changes)
        b) Universal + hallmarks only (for top hallmarks)
        c) Universal + hallmarks + all verniers (for YQRL and high-consensus)
        
        All include IMGT2→V.
        """
        candidates = []
        used = set()
        
        # Build Universal framework positions with CDRs grafted in
        universal_positions = self._build_universal_positions_with_cdrs(cdrs)
        
        # --- Control a: Pure Universal graft (no hallmark changes, just IMGT2→V) ---
        generation_order += 1
        imgt2_mut = []
        if universal_positions.get(2, '') != 'V':
            imgt2_mut = [Mutation(
                position="IMGT2",
                imgt_num=2,
                original=universal_positions.get(2, 'E'),
                mutant='V',
                source="vhh_consensus_100pct",
                confidence=1.0
            )]
        
        candidate = self._build_universal_candidate(
            f"Universal_graft_{generation_order}",
            cdrs, family, imgt2_mut, "universal_graft", generation_order,
            "Track0: Universal CDR graft only - IMMUTABLE", family
        )
        candidate.design_track = "minimal_hallmark"
        candidate.ranking_exempt = True
        candidate.track_info = "Universal_graft"
        candidate.is_immutable = True
        candidate.scaffold_type = "universal"
        candidates.append(candidate)
        
        # Select top hallmarks to test (max 4)
        test_hallmarks = hallmark_pool[:4] if len(hallmark_pool) > 4 else hallmark_pool
        
        # Always include YQRL if available
        if 'YQRL' in hallmark_pool and 'YQRL' not in test_hallmarks:
            test_hallmarks = ['YQRL'] + test_hallmarks[:3]
        
        for hallmark in test_hallmarks:
            # Get hallmark mutations
            hallmark_muts = self._get_hallmark_mutations_for_universal(hallmark, universal_positions)
            
            # --- Control b: Universal + hallmarks only ---
            muts = list(imgt2_mut) + hallmark_muts
            key = ("universal", hallmark, "hallmarks_only")
            if key not in used:
                used.add(key)
                generation_order += 1
                candidate = self._build_universal_candidate(
                    f"Universal_{hallmark}_hallmarks_{generation_order}",
                    cdrs, family, muts, f"universal_{hallmark}", generation_order,
                    f"Track0: Universal + {hallmark} hallmarks - IMMUTABLE", family
                )
                candidate.design_track = "minimal_hallmark"
                candidate.ranking_exempt = True
                candidate.track_info = f"Universal_{hallmark}_hallmarks"
                candidate.is_immutable = True
                candidate.scaffold_type = "universal"
                candidates.append(candidate)
            
            # --- Control c: Universal + hallmarks + all verniers (for YQRL or high-consensus) ---
            if hallmark == 'YQRL' or (self.hallmark_db and self._is_high_consensus_hallmark(hallmark)):
                vernier_muts = self._get_vernier_mutations_for_universal(hallmark, universal_positions)
                all_muts = list(imgt2_mut) + hallmark_muts + vernier_muts
                
                key = ("universal", hallmark, "full_vernier")
                if key not in used:
                    used.add(key)
                    generation_order += 1
                    candidate = self._build_universal_candidate(
                        f"Universal_{hallmark}_full_{generation_order}",
                        cdrs, family, all_muts, f"universal_{hallmark}_full", generation_order,
                        f"Track0: Universal + {hallmark} + all verniers - IMMUTABLE", family
                    )
                    candidate.design_track = "minimal_hallmark"
                    candidate.ranking_exempt = True
                    candidate.track_info = f"Universal_{hallmark}_full_vernier"
                    candidate.is_immutable = True
                    candidate.track_subtype = "full_vernier"
                    candidate.scaffold_type = "universal"
                    candidates.append(candidate)
        
        return candidates, generation_order

    def _build_universal_positions_with_cdrs(self, cdrs: CDRSet) -> Dict[int, str]:
        """Build IMGT position dict from Universal scaffold with CDRs grafted in."""
        positions = {}
        
        # FR1 (positions 1-26)
        fr1 = UNIVERSAL_FRAMEWORK['FR1']
        for i, aa in enumerate(fr1):
            positions[i + 1] = aa
        
        # CDR1 (positions 27-38)
        for i, aa in enumerate(cdrs.cdr1):
            if i + 27 <= 38:
                positions[i + 27] = aa
        
        # FR2 (positions 39-55)
        fr2 = UNIVERSAL_FRAMEWORK['FR2']
        for i, aa in enumerate(fr2):
            positions[i + 39] = aa
        
        # CDR2 (positions 56-65)
        for i, aa in enumerate(cdrs.cdr2):
            if i + 56 <= 65:
                positions[i + 56] = aa
        
        # FR3 (positions 66-104)
        fr3 = UNIVERSAL_FRAMEWORK['FR3']
        for i, aa in enumerate(fr3):
            positions[i + 66] = aa
        
        # CDR3 (positions 105-117)
        for i, aa in enumerate(cdrs.cdr3):
            if i + 105 <= 117:
                positions[i + 105] = aa
        
        # FR4 (positions 118-128)
        fr4 = UNIVERSAL_FRAMEWORK['FR4']
        for i, aa in enumerate(fr4):
            positions[i + 118] = aa
        
        return positions

    def _get_hallmark_mutations_for_universal(self, hallmark: str, positions: Dict) -> List[Mutation]:
        """Get mutations to apply a hallmark to Universal scaffold."""
        mutations = []
        if len(hallmark) != 4:
            return mutations
        
        hallmark_map = {42: hallmark[0], 49: hallmark[1], 50: hallmark[2], 52: hallmark[3]}
        for pos, target_aa in hallmark_map.items():
            orig_aa = positions.get(pos, '')
            if orig_aa and orig_aa != target_aa:
                mutations.append(Mutation(
                    position=f"IMGT{pos}",
                    imgt_num=pos,
                    original=orig_aa,
                    mutant=target_aa,
                    source=f"hallmark_{hallmark}",
                    confidence=1.0
                ))
        return mutations

    def _get_vernier_mutations_for_universal(self, hallmark: str, positions: Dict) -> List[Mutation]:
        """Get vernier mutations for a hallmark applied to Universal scaffold."""
        mutations = []
        
        if not self.hallmark_db:
            return mutations
        
        vernier_profile = self.hallmark_db.vernier_profile(hallmark)
        
        for pos in VERNIER_POSITIONS_FR3:
            data = vernier_profile.get(pos, {})
            cons_aa = data.get('cons', '')
            freq = data.get('freq', 0)
            orig_aa = positions.get(pos, '')
            
            if cons_aa and orig_aa and cons_aa != orig_aa and freq >= 0.50:
                mutations.append(Mutation(
                    position=f"IMGT{pos}",
                    imgt_num=pos,
                    original=orig_aa,
                    mutant=cons_aa,
                    source=f"vernier_consensus_{hallmark}",
                    confidence=freq
                ))
        
        return mutations

    def _is_high_consensus_hallmark(self, hallmark: str) -> bool:
        """Check if a hallmark has high average vernier consensus."""
        if not self.hallmark_db:
            return False
        vernier_profile = self.hallmark_db.vernier_profile(hallmark)
        freqs = [v.get('freq', 0) for v in vernier_profile.values() if v.get('freq', 0) > 0]
        return np.mean(freqs) >= 0.80 if freqs else False


    def _generate_universal(self, cdrs: CDRSet, family: str, rules: List[dict],
                             archetype: Dict[int, str], triplet_rules: List[dict],
                             n: int, target_hallmarks: str,
                             generation_order: int) -> Tuple[List[VHHCandidate], int]:
        """
        Generate using universal scaffold with PROBABILISTIC consensus sprinkling.
        
        v8.3: No more all-or-nothing!
        - Uses UNIVERSAL_PRIOR for tier-based probabilistic sprinkling
        - Tier A (≥90%): ~60% probability of application
        - Tier B (75-90%): ~25% probability of application
        - Generates intermediate states (not just "hallmarks only" or "full")
        
        Track structure for Universal:
        - Control: Hallmarks only (immutable)
        - Control: Full consensus (immutable)
        - Track 2-style: Motif pairs (1-3 pairs, ranked within track)
        - Track 3-style: Triplets + partial bundles (ranked within track)
        - Track 4: Probabilistic sprinkling (ranked globally)
        """
        import random
        rng = random.Random(42 + generation_order)
        
        candidates = []
        used = set()
        
        hallmark_set = set(HALLMARK_POSITIONS)
        
        # =====================================================================
        # BUILD HALLMARK MUTATIONS (from target_hallmarks)
        # =====================================================================
        hallmark_muts = []
        if len(target_hallmarks) == 4:
            hallmark_map = {42: target_hallmarks[0], 49: target_hallmarks[1], 
                            50: target_hallmarks[2], 52: target_hallmarks[3]}
            for imgt_pos, target_aa in hallmark_map.items():
                original_aa = self.scaffold_positions.get(imgt_pos, '')
                if original_aa and original_aa != target_aa:
                    hallmark_muts.append(Mutation(
                        position=f"IMGT{imgt_pos}",
                        imgt_num=imgt_pos,
                        original=original_aa,
                        mutant=target_aa,
                        source=f"universal_{target_hallmarks}",
                        confidence=0.95
                    ))
        
        # =====================================================================
        # BUILD VERNIER MUTATIONS FROM UNIVERSAL_PRIOR (v8.3)
        # =====================================================================
        tier_a_muts = []  # High confidence (≥90%)
        tier_b_muts = []  # Medium confidence (75-90%)
        
        for pos, (consensus_aa, freq, tier) in UNIVERSAL_PRIOR.items():
            if pos in hallmark_set or pos in INVARIANT_VHH_POSITIONS:
                continue
            if pos in [1, 2, 128]:  # FR1/FR4 handled separately
                continue
            
            original_aa = self.scaffold_positions.get(pos, '')
            if not original_aa or original_aa == consensus_aa:
                continue
            
            mut = Mutation(
                position=f"IMGT{pos}",
                imgt_num=pos,
                original=original_aa,
                mutant=consensus_aa,
                source=f"universal_prior_{tier}",
                confidence=freq
            )
            
            if tier == 'A':
                tier_a_muts.append(mut)
            else:
                tier_b_muts.append(mut)
        
        all_vernier_muts = tier_a_muts + tier_b_muts
        
        # =====================================================================
        # TRACK 0: CONTROLS (immutable)
        # =====================================================================
        # Control 1: Hallmarks only
        if hallmark_muts:
            key = tuple(sorted([str(m) for m in hallmark_muts]))
            if key not in used:
                used.add(key)
                generation_order += 1
                candidate = self._build_universal_candidate(
                    f"Univ_{target_hallmarks}_hallmarks_{generation_order}",
                    cdrs, family, hallmark_muts, "universal_hallmarks", generation_order,
                    f"Universal + {target_hallmarks} hallmarks only", family
                )
                candidate.design_track = "minimal_hallmark"
                candidate.ranking_exempt = True
                candidate.track_info = f"Universal_{target_hallmarks}_hallmarks"
                candidate.is_immutable = True
                candidates.append(candidate)
        
        # Control 2: Hallmarks + all verniers (full consensus)
        if hallmark_muts and all_vernier_muts:
            all_muts = list(hallmark_muts) + all_vernier_muts
            key = tuple(sorted([str(m) for m in all_muts]))
            if key not in used:
                used.add(key)
                generation_order += 1
                candidate = self._build_universal_candidate(
                    f"Univ_{target_hallmarks}_full_{generation_order}",
                    cdrs, family, all_muts, "universal_full", generation_order,
                    f"Universal + {target_hallmarks} + all verniers", family
                )
                candidate.design_track = "minimal_hallmark"
                candidate.ranking_exempt = True
                candidate.track_info = f"Universal_{target_hallmarks}_full_vernier"
                candidate.is_immutable = True
                candidates.append(candidate)
        
        # =====================================================================
        # TRACK 2-STYLE: MOTIF PAIRS (1-3 pairs per candidate)
        # =====================================================================
        n_motif_pairs = min(4, n // 5)
        
        for _ in range(n_motif_pairs):
            if len(candidates) >= n:
                break
            
            muts = list(hallmark_muts)
            used_positions = {m.imgt_num for m in muts}
            
            # Sample 1-3 motif pairs
            n_pairs = rng.randint(1, min(3, len(MOTIF_PAIRS)))
            available_pairs = [p for p in MOTIF_PAIRS 
                             if p[0] not in used_positions and p[2] not in used_positions]
            
            if available_pairs:
                selected_pairs = rng.sample(available_pairs, min(n_pairs, len(available_pairs)))
                
                for pos1, aa1, pos2, aa2, freq in selected_pairs:
                    # Add first position of pair
                    orig1 = self.scaffold_positions.get(pos1, '')
                    if orig1 and orig1 != aa1 and pos1 not in used_positions:
                        muts.append(Mutation(
                            position=f"IMGT{pos1}", imgt_num=pos1,
                            original=orig1, mutant=aa1,
                            source=f"motif_pair_{pos1}_{pos2}",
                            confidence=freq
                        ))
                        used_positions.add(pos1)
                    
                    # Add second position of pair
                    orig2 = self.scaffold_positions.get(pos2, '')
                    if orig2 and orig2 != aa2 and pos2 not in used_positions:
                        muts.append(Mutation(
                            position=f"IMGT{pos2}", imgt_num=pos2,
                            original=orig2, mutant=aa2,
                            source=f"motif_pair_{pos1}_{pos2}",
                            confidence=freq
                        ))
                        used_positions.add(pos2)
            
            key = tuple(sorted([str(m) for m in muts]))
            if key not in used and len(muts) > len(hallmark_muts):
                used.add(key)
                generation_order += 1
                candidate = self._build_universal_candidate(
                    f"Univ_{target_hallmarks}_pairs_{generation_order}",
                    cdrs, family, muts, "motif_pairs", generation_order,
                    f"Universal + motif pairs ({len(muts)} mut)", family
                )
                candidate.design_track = "paired_vernier"
                candidate.ranking_exempt = True  # Ranked within track, not globally
                candidate.track_info = f"Universal_motif_pairs"
                candidates.append(candidate)
        
        # =====================================================================
        # TRACK 3A-STYLE: MOTIF TRIPLETS (1-2 triplets per candidate)
        # =====================================================================
        n_triplets = min(3, n // 6)
        
        for _ in range(n_triplets):
            if len(candidates) >= n:
                break
            
            muts = list(hallmark_muts)
            used_positions = {m.imgt_num for m in muts}
            
            # Sample 1-2 motif triplets
            n_trip = rng.randint(1, min(2, len(MOTIF_TRIPLETS)))
            available_trips = [t for t in MOTIF_TRIPLETS 
                             if t[0] not in used_positions and t[2] not in used_positions and t[4] not in used_positions]
            
            if available_trips:
                selected_trips = rng.sample(available_trips, min(n_trip, len(available_trips)))
                
                for pos1, aa1, pos2, aa2, pos3, aa3, freq in selected_trips:
                    for pos, aa in [(pos1, aa1), (pos2, aa2), (pos3, aa3)]:
                        orig = self.scaffold_positions.get(pos, '')
                        if orig and orig != aa and pos not in used_positions:
                            muts.append(Mutation(
                                position=f"IMGT{pos}", imgt_num=pos,
                                original=orig, mutant=aa,
                                source=f"motif_triplet",
                                confidence=freq
                            ))
                            used_positions.add(pos)
            
            key = tuple(sorted([str(m) for m in muts]))
            if key not in used and len(muts) > len(hallmark_muts):
                used.add(key)
                generation_order += 1
                candidate = self._build_universal_candidate(
                    f"Univ_{target_hallmarks}_triplet_{generation_order}",
                    cdrs, family, muts, "motif_triplet", generation_order,
                    f"Universal + motif triplet ({len(muts)} mut)", family
                )
                candidate.design_track = "triplet_vernier"
                candidate.ranking_exempt = True
                candidate.track_info = f"Universal_motif_triplet"
                candidates.append(candidate)
        
        # =====================================================================
        # TRACK 3B-STYLE: PARTIAL CONSENSUS BUNDLES (3-6 Tier-A positions)
        # =====================================================================
        n_bundles = min(3, n // 6)
        
        # Only generate bundles if we have at least 3 Tier-A positions
        if len(tier_a_muts) >= 3:
            for _ in range(n_bundles):
                if len(candidates) >= n:
                    break
                
                muts = list(hallmark_muts)
                used_positions = {m.imgt_num for m in muts}
                
                # Sample 3-6 from Tier A positions
                available = [m for m in tier_a_muts if m.imgt_num not in used_positions]
                if len(available) < 3:
                    break
                
                n_to_pick = rng.randint(3, min(6, len(available)))
                picked = rng.sample(available, n_to_pick)
                muts.extend(picked)
                
                key = tuple(sorted([str(m) for m in muts]))
                if key not in used and len(muts) > len(hallmark_muts):
                    used.add(key)
                    generation_order += 1
                    candidate = self._build_universal_candidate(
                        f"Univ_{target_hallmarks}_bundle_{generation_order}",
                        cdrs, family, muts, "partial_bundle", generation_order,
                        f"Universal + Tier-A bundle ({len(muts)} mut)", family
                    )
                    candidate.design_track = "triplet_vernier"
                    candidate.ranking_exempt = True
                    candidate.track_info = f"Universal_tierA_bundle"
                    candidate.track_subtype = "bundle"
                    candidates.append(candidate)
        
        # =====================================================================
        # TRACK 4: PROBABILISTIC SPRINKLING (ranked globally)
        # =====================================================================
        n_ranked = n - len(candidates)
        
        for _ in range(n_ranked):
            if len(candidates) >= n:
                break
            
            muts = list(hallmark_muts)
            used_positions = {m.imgt_num for m in muts}
            
            # PROBABILISTIC SPRINKLING based on tier
            # Tier A: apply with ~60% probability
            for vm in tier_a_muts:
                if vm.imgt_num not in used_positions:
                    if rng.random() < TIER_A_PROB:
                        muts.append(vm)
                        used_positions.add(vm.imgt_num)
            
            # Tier B: apply with ~25% probability
            for vm in tier_b_muts:
                if vm.imgt_num not in used_positions:
                    if rng.random() < TIER_B_PROB:
                        muts.append(vm)
                        used_positions.add(vm.imgt_num)
            
            key = tuple(sorted([str(m) for m in muts]))
            if key not in used:
                used.add(key)
                generation_order += 1
                candidate = self._build_universal_candidate(
                    f"Univ_{target_hallmarks}_sprinkle_{generation_order}",
                    cdrs, family, muts, "probabilistic_sprinkle", generation_order,
                    f"Universal + probabilistic sprinkle ({len(muts)} mut)", family
                )
                candidate.design_track = "optimized"
                candidate.ranking_exempt = False  # RANKED globally
                candidate.track_subtype = "sprinkle"
                candidates.append(candidate)
        
        return candidates[:n], generation_order
    
    def _generate_original(self, cdrs: CDRSet, family: str, rules: List[dict],
                            archetype: Dict[int, str], triplet_rules: List[dict],
                            n: int, target_hallmarks: str, generation_order: int,
                            original_positions: Dict[Any, str],
                            hallmarks_with_controls: Set[str] = None) -> Tuple[List[VHHCandidate], int]:
        """Generate using original FRs with position-based mutations.
        
        Args:
            hallmarks_with_controls: Set of hallmarks that already have controls generated.
                                    Controls are only generated if hallmark not in this set.
                                    Updated in-place when controls are generated.
        """
        if hallmarks_with_controls is None:
            hallmarks_with_controls = set()
        
        candidates = []
        used = set()
        
        # Hallmark mutations (required for VHH conversion)
        hallmark_muts = []
        # NEW (works for any 4-char hallmark)
        if len(target_hallmarks) == 4:
            hallmark_map = {42: target_hallmarks[0], 49: target_hallmarks[1], 
                            50: target_hallmarks[2], 52: target_hallmarks[3]}
            for imgt_pos, target_aa in hallmark_map.items():
                original_aa = original_positions.get(imgt_pos, '')
                if original_aa and original_aa != target_aa:
                    hallmark_muts.append(Mutation(
                        position=f"IMGT{imgt_pos}",
                        imgt_num=imgt_pos,
                        original=original_aa,
                        mutant=target_aa,
                        source=target_hallmarks,
                        confidence=0.95
                    ))
        
        # IMGT2: I→V is required for true VHH (100% consensus in all VHH hallmarks)
        # EIQL = human VH, EVQL = camelid VHH
        original_aa_2 = original_positions.get(2, '')
        if original_aa_2 and original_aa_2 != 'V':
            hallmark_muts.append(Mutation(
                position="IMGT2",
                imgt_num=2,
                original=original_aa_2,
                mutant='V',
                source="vhh_consensus_100pct",
                confidence=1.0
            ))
        
        # v8.1: Build list of mandatory consensus mutations for progressive sprinkling
        # These are applied in Tracks 1-4 with increasing probability
        mandatory_consensus_muts = []
        
        # IMGT1: E→Q (first position of FR1, ~92% in VHH)
        original_aa_1 = original_positions.get(1, '')
        if original_aa_1 and original_aa_1 != 'Q':
            mandatory_consensus_muts.append(Mutation(
                position="IMGT1",
                imgt_num=1,
                original=original_aa_1,
                mutant='Q',
                source="vhh_consensus_92pct_FR1",
                confidence=0.92
            ))
        
        # IMGT128: A→S (last position of FR4, makes ...TVSS ending, ~95% in VHH)
        original_aa_128 = original_positions.get(128, '')
        if original_aa_128 and original_aa_128 != 'S':
            mandatory_consensus_muts.append(Mutation(
                position="IMGT128",
                imgt_num=128,
                original=original_aa_128,
                mutant='S',
                source="vhh_consensus_95pct_FR4",
                confidence=0.95
            ))
        
        # =================================================================
        # STRATIFIED GENERATION with per-hallmark vernier consensus
        # =================================================================
        
        # Get per-hallmark vernier profile from database
        vernier_profile, epistatic_pairs = self._get_hallmark_profile(target_hallmarks)
        
        # Vernier positions to consider (excluding hallmark positions which are handled above)
        # Position 52 is in both hallmark and vernier - already handled in hallmarks
        VERNIER_POSITIONS_ALL = [2, 4, 41, 47, 66, 67, 68, 69, 71, 76, 78, 82, 87, 89, 91, 94]
        
        # Highly conserved positions - skip if already correct (>97% conservation)
        # Per analysis: IMGT4=L(99.6%), IMGT41=W(97.8%), IMGT47=G(97.4%)
        SKIP_IF_CORRECT = {4, 41, 47}
        
        # Variable position - randomize between consensus and original
        # Position 66 has only 43% conservation and varies by hallmark (Y/N/S/D/E/K...)
        VARIABLE_POSITION_66 = 66
        
        # Build vernier consensus mutations from profile
        vernier_consensus_muts = []
        pos66_data = None  # Store position 66 data for special handling
        
        for pos in VERNIER_POSITIONS_ALL:
            if pos == 2:  # Already handled in hallmarks
                continue
            if pos in [42, 49, 50, 52]:  # Hallmark positions
                continue
                
            original_aa = original_positions.get(pos, '')
            if not original_aa:
                continue
            
            # Get consensus from hallmark-specific profile
            profile_data = vernier_profile.get(pos, {})
            consensus_aa = profile_data.get('cons', '')
            consensus_freq = float(profile_data.get('freq', 0))
            alt_aa = profile_data.get('alt', '')
            alt_freq = float(profile_data.get('alt_freq', 0))
            
            if not consensus_aa or consensus_aa == '-':
                continue
            
            # Skip if already correct AND position is in SKIP_IF_CORRECT
            if pos in SKIP_IF_CORRECT and original_aa == consensus_aa:
                continue
            
            # Position 66: store for randomization (50/50 consensus vs original)
            if pos == VARIABLE_POSITION_66:
                if original_aa != consensus_aa:
                    pos66_data = (pos, original_aa, consensus_aa, consensus_freq, alt_aa, alt_freq)
                continue
            
            # Create mutation if different from original
            if original_aa != consensus_aa:
                mut = Mutation(
                    position=f"IMGT{pos}",
                    imgt_num=pos,
                    original=original_aa,
                    mutant=consensus_aa,
                    source=f"vernier_cons_{target_hallmarks}|freq={consensus_freq:.2f}",
                    confidence=consensus_freq
                )
                vernier_consensus_muts.append(mut)
        
        # Rule mutations (from compensation rules)
        rule_muts = []
        applied_rule_ids = []
        for rule in rules:
            imgt_num = rule.get('position_num', 0)
            suggested_aa = rule.get('suggested_aa', '')
            
            if not imgt_num or not suggested_aa or imgt_num in HALLMARK_POSITIONS:
                continue
            if imgt_num == 2:  # Already handled
                continue
            
            original_aa = original_positions.get(imgt_num, '')
            if not original_aa or original_aa == suggested_aa:
                continue
            
            # Check if this conflicts with vernier consensus
            # Rules should complement, not override vernier consensus
            rule_id = f"rule|{family}|{rule.get('condition','')}|IMGT{imgt_num}->{suggested_aa}|conf={rule.get('confidence',0):.3f}|sup={rule.get('support',0)}"
            rule_muts.append(Mutation(
                position=f"IMGT{imgt_num}",
                imgt_num=imgt_num,
                original=original_aa,
                mutant=suggested_aa,
                source=rule_id,
                confidence=rule.get('confidence', 0.5)
            ))
            applied_rule_ids.append(rule_id)
        
        # Triplet mutation (FR4 N-terminus)
        triplet_mut = None
        if triplet_rules:
            best = triplet_rules[0]
            fr4_motif = best.get('suggested_aa', '')
            if len(fr4_motif) == 3 and len(cdrs.fr4) >= 3:
                if cdrs.fr4[:3] != fr4_motif:
                    triplet_id = f"triplet|{family}|{best.get('condition','')}|IMGT118-120->{fr4_motif}|conf={best.get('confidence',0):.3f}|sup={best.get('support',0)}"
                    triplet_mut = Mutation(
                        position="IMGT118-120",
                        imgt_num="118-120",
                        original=cdrs.fr4[:3],
                        mutant=fr4_motif,
                        source=triplet_id,
                        confidence=best.get('confidence', 0.8)
                    )
        
        
        # =================================================================
        # TRACK-BASED GENERATION (v7.12)
        # =================================================================
        # Track 0: Controls (ranking_exempt=True) - minimal hallmark variants
        # Track 1: Single vernier probes (ranking_exempt=True)
        # Track 2: Paired vernier motifs (ranking_exempt=True)
        # Track 3: Triplet motifs (ranking_exempt=True)
        # Track 4: Optimized (ranking_exempt=False) - full consensus, rules, diversity
        # =================================================================
        import random
        rng = random.Random()
        
        # Print vernier consensus being applied
        print(f"    Vernier consensus ({target_hallmarks}): {len(vernier_consensus_muts)} positions")
        for m in vernier_consensus_muts[:5]:  # Show first 5
            print(f"      IMGT{m.imgt_num}: {m.original}→{m.mutant} ({m.confidence:.0%})")
        if len(vernier_consensus_muts) > 5:
            print(f"      ... and {len(vernier_consensus_muts) - 5} more")
        if pos66_data:
            print(f"    Position 66 (variable): {pos66_data[1]}→{pos66_data[2]} ({pos66_data[3]:.0%}) - will randomize 50/50")
        
        # =================================================================
        # TRACK 0: CONTROLS (always generated, ranking exempt)
        # =================================================================
        # v8.0 philosophy: The GOAL is to find BINDERS. Track 0 explores different
        # fundamental approaches - minimal changes vs full conversion. We don't
        # expect one to be "better" than another; we're casting a wide net.
        #
        # IMPORTANT: Track 0 candidates are IMMUTABLE - no compensation rules,
        # no consensus sprinkling, no post-generation modifications.
        # =================================================================
        
        # Track 0a: Hallmarks only (4 mutations, NO IMGT2)
        # Tests: Can minimal hallmark changes alone enable VHH function?
        hallmark_only_muts = [m for m in hallmark_muts if m.imgt_num != 2]
        if hallmark_only_muts:
            generation_order += 1
            key = tuple(sorted([str(m) for m in hallmark_only_muts]))
            if key not in used:
                used.add(key)
                candidate = self._build_original_candidate(
                    f"Orig_{family}_minimal_4mut_{generation_order}",
                    cdrs, family, hallmark_only_muts, "minimal_hallmark", generation_order,
                    f"Track0: Minimal hallmark only (4 mut) - IMMUTABLE", family, original_positions
                )
                candidate.design_track = "minimal_hallmark"
                candidate.ranking_exempt = True
                candidate.track_info = f"{target_hallmarks}_4mut"
                candidate.applied_comp_rules = []  # No compensation - explicitly frozen
                candidate.is_immutable = True      # v7.14: Flag as immutable
                candidates.append(candidate)
        
        # Track 0b: Hallmarks + IMGT2 (5 mutations)
        # Tests: Does adding the universal VHH pos2 consensus help?
        generation_order += 1
        key = tuple(sorted([str(m) for m in hallmark_muts]))
        if key not in used:
            used.add(key)
            candidate = self._build_original_candidate(
                f"Orig_{family}_minimal_5mut_{generation_order}",
                cdrs, family, hallmark_muts, "minimal_hallmark", generation_order,
                f"Track0: Minimal hallmark + IMGT2 (5 mut) - IMMUTABLE", family, original_positions
            )
            candidate.design_track = "minimal_hallmark"
            candidate.ranking_exempt = True
            candidate.track_info = f"{target_hallmarks}_5mut"
            candidate.applied_comp_rules = []  # No compensation - explicitly frozen
            candidate.is_immutable = True      # v7.14: Flag as immutable
            candidates.append(candidate)
        
        # =================================================================
        # Track 0c: FULL VERNIER CONVERSION (v8.0)
        # =================================================================
        # Tests: Does full VHH conversion work? Some CDRs may need complete
        # framework adaptation. This is especially valuable for YQRL (most
        # constrained) - if full YQRL works, the CDRs are broadly compatible.
        # =================================================================
        
        # Check if this hallmark has high average vernier consensus
        vernier_freqs = [vm.confidence for vm in vernier_consensus_muts]
        is_high_consensus = False
        if vernier_freqs:
            avg_vernier_freq = sum(vernier_freqs) / len(vernier_freqs)
            is_high_consensus = avg_vernier_freq >= 0.75
        
        # Generate full vernier control for YQRL or high-consensus hallmarks
        if target_hallmarks == 'YQRL' or is_high_consensus:
            all_vernier_muts = list(hallmark_muts)
            for vm in vernier_consensus_muts:
                if vm.confidence >= 0.50:  # Include all verniers with ≥50% consensus
                    all_vernier_muts.append(vm)
            
            # Include pos66 if available
            if pos66_data:
                pos, orig, cons, freq = pos66_data[:4]  # pos66_data has 6 elements
                if cons != orig:
                    all_vernier_muts.append(Mutation(
                        position=f"IMGT{pos}",
                        imgt_num=pos,
                        original=orig,
                        mutant=cons,
                        source=f"vernier_consensus_{target_hallmarks}",
                        confidence=freq
                    ))
            
            generation_order += 1
            key = tuple(sorted([str(m) for m in all_vernier_muts]))
            if key not in used:
                used.add(key)
                candidate = self._build_original_candidate(
                    f"Orig_{family}_full_vernier_{generation_order}",
                    cdrs, family, all_vernier_muts, "minimal_hallmark", generation_order,
                    f"Track0: Full vernier conversion ({target_hallmarks}) - IMMUTABLE", family, original_positions
                )
                candidate.design_track = "minimal_hallmark"
                candidate.ranking_exempt = True
                candidate.track_info = f"{target_hallmarks}_full_vernier"
                candidate.applied_comp_rules = []
                candidate.is_immutable = True
                candidate.track_subtype = "full_vernier"
                candidates.append(candidate)
        
        # =================================================================
        # TRACK 1: SINGLE VERNIER PROBES (ranking exempt)
        # =================================================================
        # v8.1: Reduced to fit within 92-control budget
        # Test top load-bearing verniers only
        # =================================================================
        
        MAX_SINGLE_PROBES = 6  # v8.1: Reduced from 12 to fit budget
        
        # v8.0: Focus on load-bearing positions that have data for this hallmark
        load_bearing_verniers = []
        for vm in vernier_consensus_muts:
            if vm.imgt_num in LOAD_BEARING_VERNIERS and vm.confidence >= 0.50:
                load_bearing_verniers.append(vm)
        
        # Sort by confidence (highest first)
        load_bearing_verniers.sort(key=lambda m: -m.confidence)
        
        # v8.1: Add FR1/FR4 consensus with light probability (30%)
        for vm in load_bearing_verniers[:MAX_SINGLE_PROBES]:
            muts = list(hallmark_muts) + [vm]
            # Light FR1/FR4 sprinkling for Track 1
            if rng.random() < 0.30:
                for mc in mandatory_consensus_muts:
                    if mc.imgt_num not in {m.imgt_num for m in muts}:
                        muts.append(mc)
            key = tuple(sorted([str(m) for m in muts]))
            if key not in used:
                used.add(key)
                generation_order += 1
                candidate = self._build_original_candidate(
                    f"Orig_{family}_probe_IMGT{vm.imgt_num}_{generation_order}",
                    cdrs, family, muts, "single_vernier", generation_order,
                    f"Track1: Single vernier probe IMGT{vm.imgt_num}", family, original_positions
                )
                candidate.design_track = "single_vernier"
                candidate.ranking_exempt = True
                candidate.track_info = f"IMGT{vm.imgt_num}_{vm.original}>{vm.mutant}"
                candidate.applied_comp_rules = []  # No compensation for probes
                candidates.append(candidate)
        
        # =================================================================
        # TRACK 2: PAIRED VERNIER MOTIFS (ranking exempt)
        # =================================================================
        # v8.3: Use cross-family MOTIF_PAIRS - positions that co-occur
        # across multiple hallmarks are likely structurally coupled.
        # =================================================================
        
        MAX_PAIR_PROBES = 8  # v8.1: Reduced from 40 to fit budget
        
        # v8.3: Use pre-defined cross-family motif pairs
        pair_candidates = []
        for pos1, aa1, pos2, aa2, freq in MOTIF_PAIRS:
            # Check if positions need mutation in this sequence
            orig1 = original_positions.get(pos1, '')
            orig2 = original_positions.get(pos2, '')
            
            # Only include if at least one position needs change
            if (orig1 and orig1 != aa1) or (orig2 and orig2 != aa2):
                pair_candidates.append((pos1, aa1, pos2, aa2, freq))
        
        # Also add family-specific pairs from vernier consensus
        for i, vm1 in enumerate(vernier_consensus_muts):
            if vm1.confidence < 0.60:
                continue
            for vm2 in vernier_consensus_muts[i+1:]:
                if vm2.confidence < 0.60:
                    continue
                avg_conf = (vm1.confidence + vm2.confidence) / 2
                # Don't duplicate if already in cross-family list
                already_exists = any(
                    (p[0] == vm1.imgt_num and p[2] == vm2.imgt_num) or 
                    (p[0] == vm2.imgt_num and p[2] == vm1.imgt_num)
                    for p in pair_candidates
                )
                if not already_exists:
                    pair_candidates.append((vm1.imgt_num, vm1.mutant, vm2.imgt_num, vm2.mutant, avg_conf))
        
        # Sort by frequency (highest first)
        pair_candidates.sort(key=lambda x: -x[4])
        
        # Generate candidates for top pairs
        for pos1, aa1, pos2, aa2, freq in pair_candidates[:MAX_PAIR_PROBES]:
            muts = list(hallmark_muts)
            used_positions = {m.imgt_num for m in muts}
            
            # Add pair mutations
            orig1 = original_positions.get(pos1, '')
            if orig1 and orig1 != aa1 and pos1 not in used_positions:
                muts.append(Mutation(
                    position=f"IMGT{pos1}", imgt_num=pos1,
                    original=orig1, mutant=aa1,
                    source=f"motif_pair_{pos1}_{pos2}",
                    confidence=freq
                ))
                used_positions.add(pos1)
            
            orig2 = original_positions.get(pos2, '')
            if orig2 and orig2 != aa2 and pos2 not in used_positions:
                muts.append(Mutation(
                    position=f"IMGT{pos2}", imgt_num=pos2,
                    original=orig2, mutant=aa2,
                    source=f"motif_pair_{pos1}_{pos2}",
                    confidence=freq
                ))
                used_positions.add(pos2)
            
            # v8.1: Moderate FR1/FR4 sprinkling for Track 2 (50%)
            if rng.random() < 0.50:
                for mc in mandatory_consensus_muts:
                    if mc.imgt_num not in used_positions:
                        muts.append(mc)
            
            key = tuple(sorted([str(m) for m in muts]))
            if key not in used and len(muts) > len(hallmark_muts):
                used.add(key)
                generation_order += 1
                candidate = self._build_original_candidate(
                    f"Orig_{family}_pair_IMGT{pos1}_{pos2}_{generation_order}",
                    cdrs, family, muts, "paired_vernier", generation_order,
                    f"Track2: Motif pair IMGT{pos1}+{pos2}", family, original_positions
                )
                candidate.design_track = "paired_vernier"
                candidate.ranking_exempt = True
                candidate.track_info = f"IMGT{pos1}+IMGT{pos2}"
                candidate.applied_comp_rules = []
                candidates.append(candidate)
        
        # =================================================================
        # TRACK 3: TRIPLET + PARTIAL BUNDLES (ranking exempt)
        # =================================================================
        # v8.3: Two sub-lanes:
        #   3A: Motif TRIPLETS from cross-family analysis
        #   3B: Partial consensus BUNDLES (3-6 Tier-A positions sampled)
        # =================================================================
        
        MAX_TRIPLET_PROBES = 5
        MAX_BUNDLE_PROBES = 4
        
        # --- 3A: Cross-family Motif Triplets ---
        # v8.3: Use pre-defined MOTIF_TRIPLETS
        triplet_candidates = []
        for pos1, aa1, pos2, aa2, pos3, aa3, freq in MOTIF_TRIPLETS:
            # Check if positions need mutation
            orig1 = original_positions.get(pos1, '')
            orig2 = original_positions.get(pos2, '')
            orig3 = original_positions.get(pos3, '')
            
            # Include if at least one position needs change
            needs_change = ((orig1 and orig1 != aa1) or 
                           (orig2 and orig2 != aa2) or 
                           (orig3 and orig3 != aa3))
            if needs_change:
                triplet_candidates.append((pos1, aa1, pos2, aa2, pos3, aa3, freq))
        
        # Also add family-specific triplets from vernier consensus
        for i, vm1 in enumerate(vernier_consensus_muts):
            if vm1.confidence < 0.70:
                continue
            for j, vm2 in enumerate(vernier_consensus_muts[i+1:], i+1):
                if vm2.confidence < 0.70:
                    continue
                for vm3 in vernier_consensus_muts[j+1:]:
                    if vm3.confidence < 0.70:
                        continue
                    avg_conf = (vm1.confidence + vm2.confidence + vm3.confidence) / 3
                    triplet_candidates.append((
                        vm1.imgt_num, vm1.mutant,
                        vm2.imgt_num, vm2.mutant,
                        vm3.imgt_num, vm3.mutant,
                        avg_conf
                    ))
        
        # Sort by frequency
        triplet_candidates.sort(key=lambda x: -x[6])
        
        for pos1, aa1, pos2, aa2, pos3, aa3, freq in triplet_candidates[:MAX_TRIPLET_PROBES]:
            muts = list(hallmark_muts)
            used_positions = {m.imgt_num for m in muts}
            
            for pos, aa in [(pos1, aa1), (pos2, aa2), (pos3, aa3)]:
                orig = original_positions.get(pos, '')
                if orig and orig != aa and pos not in used_positions:
                    muts.append(Mutation(
                        position=f"IMGT{pos}", imgt_num=pos,
                        original=orig, mutant=aa,
                        source="motif_triplet",
                        confidence=freq
                    ))
                    used_positions.add(pos)
            
            # v8.1: High FR1/FR4 sprinkling for Track 3 (70%)
            if rng.random() < 0.70:
                for mc in mandatory_consensus_muts:
                    if mc.imgt_num not in used_positions:
                        muts.append(mc)
            
            key = tuple(sorted([str(m) for m in muts]))
            if key not in used and len(muts) > len(hallmark_muts):
                used.add(key)
                generation_order += 1
                candidate = self._build_original_candidate(
                    f"Orig_{family}_triplet_{pos1}_{pos2}_{pos3}_{generation_order}",
                    cdrs, family, muts, "triplet_vernier", generation_order,
                    f"Track3A: Motif triplet IMGT{pos1}+{pos2}+{pos3}", family, original_positions
                )
                candidate.design_track = "triplet_vernier"
                candidate.ranking_exempt = True
                candidate.track_info = f"IMGT{pos1}+{pos2}+{pos3}"
                candidate.track_subtype = "motif_triplet"
                candidate.applied_comp_rules = []
                candidates.append(candidate)
        
        # --- 3B: Partial Consensus Bundles (3-6 Tier-A positions sampled) ---
        # v8.3: Sample 3-6 high-confidence positions instead of all-or-nothing
        
        # Identify Tier-A positions from vernier consensus (>=90% or >=85%)
        tier_a_muts = [vm for vm in vernier_consensus_muts if vm.confidence >= 0.85]
        
        if len(tier_a_muts) >= 3:
            for _ in range(MAX_BUNDLE_PROBES):
                muts = list(hallmark_muts)
                used_positions = {m.imgt_num for m in muts}
                
                # Sample 3-6 Tier-A positions
                available = [m for m in tier_a_muts if m.imgt_num not in used_positions]
                if len(available) < 3:
                    break
                
                n_to_pick = rng.randint(3, min(6, len(available)))
                picked = rng.sample(available, n_to_pick)
                muts.extend(picked)
                
                # v8.1: High FR1/FR4 sprinkling (80%)
                if rng.random() < 0.80:
                    for mc in mandatory_consensus_muts:
                        if mc.imgt_num not in {m.imgt_num for m in muts}:
                            muts.append(mc)
                
                key = tuple(sorted([str(m) for m in muts]))
                if key not in used and len(muts) > len(hallmark_muts):
                    used.add(key)
                    generation_order += 1
                    n_verniers = len(muts) - len(hallmark_muts)
                    candidate = self._build_original_candidate(
                        f"Orig_{family}_bundle_{n_verniers}pos_{generation_order}",
                        cdrs, family, muts, "triplet_vernier", generation_order,
                        f"Track3B: Partial bundle ({n_verniers} Tier-A positions)", family, original_positions
                    )
                    candidate.design_track = "triplet_vernier"
                    candidate.ranking_exempt = True
                    candidate.track_info = f"bundle_{n_verniers}pos"
                    candidate.track_subtype = "bundle"
                    candidate.applied_comp_rules = []
                    candidates.append(candidate)
        
        # Track controls count
        n_controls = len([c for c in candidates if c.ranking_exempt])
        
        # =================================================================
        # TRACK 4: OPTIMIZED CANDIDATES - TWO LANES (v8.0)
        # =================================================================
        # v8.0: Two parallel strategies to maximize chance of finding binders:
        # 
        # Lane A (50%): Full consensus - apply most high-confidence verniers
        #   "The safest, most VHH-like candidates"
        #
        # Lane B (50%): Minimal (2-6 verniers) - fewer changes, high-confidence biased
        #   "Maybe fewer changes are enough for this CDR"
        #
        # v8.1: All Track 4 candidates ALWAYS get FR1/FR4 consensus (IMGT1, IMGT128)
        # =================================================================
        
        # Calculate allocation for Track 4
        n_remaining = n - n_controls
        n_lane_a = n_remaining // 2  # 50% full consensus
        n_lane_b = n_remaining - n_lane_a  # 50% minimal
        
        print(f"    Track allocation: T0-T3(controls)={n_controls}, T4: Lane_A(full)={n_lane_a}, Lane_B(minimal)={n_lane_b}")
        
        # Prepare all available vernier mutations
        all_verniers = list(vernier_consensus_muts)
        if pos66_data:
            pos, orig, cons, cons_freq = pos66_data[:4]
            if cons != orig:
                all_verniers.append(Mutation(
                    position=f"IMGT{pos}",
                    imgt_num=pos,
                    original=orig,
                    mutant=cons,
                    source=f"vernier_cons_{target_hallmarks}",
                    confidence=cons_freq
                ))
        
        # Sort verniers by confidence for Lane B weighted selection
        sorted_verniers = sorted(all_verniers, key=lambda m: -m.confidence)
        
        # --- Lane A: Full Consensus ---
        for i in range(n_lane_a):
            if len(candidates) >= n:
                break
            
            muts = list(hallmark_muts)
            used_positions = {m.imgt_num for m in muts}
            
            # v8.1: ALWAYS add FR1/FR4 mandatory consensus to Track 4
            for mc in mandatory_consensus_muts:
                if mc.imgt_num not in used_positions:
                    muts.append(mc)
                    used_positions.add(mc.imgt_num)
            
            # Include verniers probabilistically based on confidence
            for vm in all_verniers:
                if vm.imgt_num not in used_positions:
                    # Higher confidence = higher inclusion probability
                    inclusion_prob = vm.confidence * 0.9  # e.g., 85% confidence -> 76.5% chance
                    if rng.random() < inclusion_prob:
                        muts.append(vm)
                        used_positions.add(vm.imgt_num)
            
            key = tuple(sorted([str(m) for m in muts]))
            
            if key not in used:
                used.add(key)
                generation_order += 1
                candidate = self._build_original_candidate(
                    f"Orig_{family}_laneA_{generation_order}",
                    cdrs, family, muts, "lane_A_consensus", generation_order,
                    f"Track4 Lane A: Full consensus ({len(muts)} mut)", family, original_positions
                )
                candidate.design_track = "optimized"
                candidate.ranking_exempt = False
                candidate.track_subtype = "lane_A"
                candidate.applied_comp_rules = []
                candidates.append(candidate)
        
        # --- Lane B: Minimal (2-6 verniers, high-confidence biased) ---
        for i in range(n_lane_b):
            if len(candidates) >= n:
                break
            
            muts = list(hallmark_muts)
            used_positions = {m.imgt_num for m in muts}
            
            # v8.1: ALWAYS add FR1/FR4 mandatory consensus to Track 4
            for mc in mandatory_consensus_muts:
                if mc.imgt_num not in used_positions:
                    muts.append(mc)
                    used_positions.add(mc.imgt_num)
            
            # Pick 2-6 verniers, biased toward high confidence
            n_to_pick = rng.randint(2, min(6, len(sorted_verniers)))
            
            # Weighted selection: higher confidence = more likely to be picked
            weights = [m.confidence for m in sorted_verniers if m.imgt_num not in used_positions]
            available = [m for m in sorted_verniers if m.imgt_num not in used_positions]
            
            if available and weights:
                # Normalize weights
                total_weight = sum(weights)
                if total_weight > 0:
                    # Select n_to_pick positions (weighted without replacement)
                    picked_indices = set()
                    for _ in range(min(n_to_pick, len(available))):
                        # Build current weights (excluding already picked)
                        current_weights = []
                        current_indices = []
                        for idx, (m, w) in enumerate(zip(available, weights)):
                            if idx not in picked_indices and m.imgt_num not in used_positions:
                                current_weights.append(w)
                                current_indices.append(idx)
                        
                        if not current_weights:
                            break
                        
                        # Weighted random selection
                        total = sum(current_weights)
                        r = rng.random() * total
                        cumsum = 0
                        for idx, w in zip(current_indices, current_weights):
                            cumsum += w
                            if r <= cumsum:
                                picked_indices.add(idx)
                                muts.append(available[idx])
                                used_positions.add(available[idx].imgt_num)
                                break
            
            key = tuple(sorted([str(m) for m in muts]))
            if key not in used:
                used.add(key)
                generation_order += 1
                candidate = self._build_original_candidate(
                    f"Orig_{family}_laneB_{generation_order}",
                    cdrs, family, muts, "lane_B_minimal", generation_order,
                    f"Track4 Lane B: Minimal ({len(muts)} mut)", family, original_positions
                )
                candidate.design_track = "optimized"
                candidate.ranking_exempt = False
                candidate.track_subtype = "lane_B"
                candidate.applied_comp_rules = []
                candidates.append(candidate)
        
        return candidates, generation_order
    
    def _build_universal_candidate(self, name: str, cdrs: CDRSet, family: str,
                                     mutations: List[Mutation], strategy: str,
                                     gen_order: int, description: str,
                                     target_family: str) -> VHHCandidate:
        """Build candidate from universal scaffold with position-based mutations."""
        # Start with scaffold positions
        positions = self.scaffold_positions.copy()
        
        # Apply mutations to positions
        for mut in mutations:
            if mut.imgt_num == "118-120" and len(mut.mutant) == 3:
                positions[118] = mut.mutant[0]
                positions[119] = mut.mutant[1]
                positions[120] = mut.mutant[2]
            else:
                positions[mut.imgt_num] = mut.mutant
        
        # Rebuild FRs from positions
        fr1 = ''.join(positions.get(i, '') for i in range(1, 27))
        fr2 = ''.join(positions.get(i, '') for i in range(39, 56))
        fr3 = ''.join(positions.get(i, '') for i in range(66, 105))
        fr4 = ''.join(positions.get(i, '') for i in range(118, 129))
        
        # Add CDRs to positions
        for i, aa in enumerate(cdrs.cdr1):
            positions[27 + i] = aa
        for i, aa in enumerate(cdrs.cdr2):
            positions[56 + i] = aa
        for i, aa in enumerate(cdrs.cdr3):
            positions[105 + i] = aa
        
        sequence = fr1 + cdrs.cdr1 + fr2 + cdrs.cdr2 + fr3 + cdrs.cdr3 + fr4
        detected_family = classify_family_from_positions(positions, sequence)
        
        cand = VHHCandidate(
            id=name, rank=0, sequence=sequence,
            framework_source='universal', family=detected_family,
            cdr1=cdrs.cdr1, cdr2=cdrs.cdr2, cdr3=cdrs.cdr3,
            fr1=fr1, fr2=fr2, fr3=fr3, fr4=fr4,
            imgt_positions=positions,
            mutations=mutations, strategy=strategy,
            generation_order=gen_order,
            construction_method=description,
            target_family=target_family
        )
        
        # v8.2: Always mark Universal candidates with scaffold_type
        cand.scaffold_type = "universal"

        # Provenance: record which rule-derived mutations were actually used
        cand.applied_comp_rules = sorted({m.source for m in mutations if isinstance(m.source, str) and m.source.startswith('rule|')})
        cand.applied_triplet_rules = sorted({m.source for m in mutations if isinstance(m.source, str) and m.source.startswith('triplet|')})
        cand.generation_debug = {
            'framework_source': 'universal',
            'strategy': strategy,
            'construction_method': description,
        }
        return cand
    
    def _build_original_candidate(self, name: str, cdrs: CDRSet, family: str,
                                    mutations: List[Mutation], strategy: str,
                                    gen_order: int, description: str,
                                    target_family: str,
                                    original_positions: Dict[Any, str]) -> VHHCandidate:
        """Build candidate from original FRs with position-based mutations."""
        # Start with original positions
        positions = original_positions.copy()
        
        # Apply mutations
        for mut in mutations:
            if mut.imgt_num == "118-120" and len(mut.mutant) == 3:
                positions[118] = mut.mutant[0]
                positions[119] = mut.mutant[1]
                positions[120] = mut.mutant[2]
            else:
                positions[mut.imgt_num] = mut.mutant
        
        # Rebuild FRs from positions
        fr1 = ''.join(positions.get(i, '') for i in range(1, 27))
        fr2 = ''.join(positions.get(i, '') for i in range(39, 56))
        fr3 = ''.join(positions.get(i, '') for i in range(66, 105))
        fr4 = ''.join(positions.get(i, '') for i in range(118, 129))
        
        sequence = fr1 + cdrs.cdr1 + fr2 + cdrs.cdr2 + fr3 + cdrs.cdr3 + fr4
        detected_family = classify_family_from_positions(positions, sequence)
        
        cand = VHHCandidate(
            id=name, rank=0, sequence=sequence,
            framework_source='original', family=detected_family,
            cdr1=cdrs.cdr1, cdr2=cdrs.cdr2, cdr3=cdrs.cdr3,
            fr1=fr1, fr2=fr2, fr3=fr3, fr4=fr4,
            imgt_positions=positions,
            mutations=mutations, strategy=strategy,
            generation_order=gen_order,
            construction_method=description,
            target_family=target_family
        )

        cand.applied_comp_rules = sorted({m.source for m in mutations if isinstance(m.source, str) and m.source.startswith('rule|')})
        cand.applied_triplet_rules = sorted({m.source for m in mutations if isinstance(m.source, str) and m.source.startswith('triplet|')})
        cand.generation_debug = {
            'framework_source': 'original',
            'strategy': strategy,
            'construction_method': description,
        }
        return cand

# ============================================================
# OUTPUT
# ============================================================

def to_dataframe(candidates: List[VHHCandidate]) -> 'pd.DataFrame':
    """Convert candidates to comprehensive DataFrame."""
    if not PANDAS_AVAILABLE:
        return None
    
    rows = []
    
    # v8.2: Separate controls and ranked for proper numbering
    controls = [c for c in candidates if c.ranking_exempt or c.is_lead]
    ranked = [c for c in candidates if not c.ranking_exempt and not c.is_lead]
    
    # Sort ranked by combined_score (ascending = best first)
    ranked.sort(key=lambda x: x.scoring.combined_score)
    
    overall_row = 0
    ranked_position = 0
    
    for c in candidates:
        overall_row += 1
        
        top_family = ""
        top_family_prob = 0.0
        if c.scoring.family_probabilities:
            top_family, top_family_prob = max(c.scoring.family_probabilities.items(), key=lambda x: x[1])
        
        # Get hallmarks from positions
        hallmarks = get_hallmarks_from_positions(c.imgt_positions)
        
        # Get cysteine class based on IMGT55+IMGT100 positions (v7.16)
        cys_class = detect_cysteine_class(c.imgt_positions, c.sequence)
        n_cys = c.sequence.count('C')
        
        # v8.2: Rank format
        # - Lead: "Lead"
        # - Controls: "(Control)"
        # - Ranked: 1, 2, 3... (1 = best)
        if c.is_lead:
            rank_str = "Lead"
            status = "LEAD"
        elif c.ranking_exempt:
            rank_str = "(Control)"
            status = "CONTROL"
        else:
            # Find position in ranked list
            ranked_position = next((i+1 for i, r in enumerate(ranked) if r.id == c.id), 0)
            rank_str = str(ranked_position)
            status = "RANKED"
        
        row = {
            'overall_row': overall_row,  # v8.2: Absolute position
            'rank': rank_str,            # v8.2: "Lead", "(Control)", or 1,2,3...
            'status': status,            # v8.1: LEAD/CONTROL/RANKED
            'generation_order': c.generation_order,
            'id': c.id,
            'is_lead': c.is_lead,
            
            # Design track fields (v7.12/7.14/v8.1)
            'design_track': c.design_track,
            'track_subtype': getattr(c, 'track_subtype', ''),  # v8.1
            'scaffold_type': getattr(c, 'scaffold_type', 'original'),  # v8.1
            'ranking_exempt': c.ranking_exempt,
            'track_info': c.track_info,
            'track_rank': c.track_rank,
            'is_immutable': c.is_immutable,
            
            'detected_family': c.family,
            'target_family': c.target_family,
            'detected_cys_class': cys_class,  # v7.16: C2/C4/INVALID based on pos55+pos100
            'n_cysteines': n_cys,              # v7.16: Total cysteine count
            'top_prob_family': top_family,
            'top_family_prob': round(top_family_prob, 3),
            'hallmarks': hallmarks,
            
            'framework_source': c.framework_source,
            'strategy': c.strategy,
            'construction_method': c.construction_method,
            'generation_debug': json.dumps(c.generation_debug, sort_keys=True),
            'n_mutations': len(c.mutations),
            'mutations': ';'.join(str(m) for m in c.mutations),
            
            'framework_identity_pct': round(c.framework_identity_pct, 1),
            'vernier_matches': c.scoring.vernier_matches,
            'vernier_total': c.scoring.vernier_total,
            
            'rules_passed': c.scoring.rules_passed,
            'rules_applicable': c.scoring.rules_applicable,
            'rules_total': c.scoring.rules_total,
            'weighted_naturalness': round(c.scoring.weighted_naturalness, 4),

            'applied_comp_rules': '|'.join(c.applied_comp_rules),
            'applied_triplet_rules': '|'.join(c.applied_triplet_rules),
            'generation_debug': json.dumps(c.generation_debug, sort_keys=True),
            
            'esm2_loss': round(c.scoring.esm2_loss, 4),
            'esm2_perplexity': round(c.scoring.esm2_perplexity, 2),
            
            'plddt_mean': round(c.scoring.plddt_mean, 1),
            'plddt_median': round(c.scoring.plddt_median, 1),
            'plddt_min': round(c.scoring.plddt_min, 1),
            'plddt_cdr1': round(c.scoring.plddt_cdr1, 1),
            'plddt_cdr2': round(c.scoring.plddt_cdr2, 1),
            'plddt_cdr3': round(c.scoring.plddt_cdr3, 1),
            'plddt_framework': round(c.scoring.plddt_framework, 1),
            
            'combined_score': round(c.scoring.combined_score, 4) if not c.is_lead else 999.0,
            
            'top_violations': '; '.join(c.scoring.top_violations[:3]),
            
            'sequence': c.sequence,
            'seq_length': len(c.sequence),
            'cdr1': c.cdr1,
            'cdr2': c.cdr2,
            'cdr3': c.cdr3,
            'cdr3_length': len(c.cdr3),
            'fr1': c.fr1,
            'fr2': c.fr2,
            'fr3': c.fr3,
            'fr4': c.fr4,
        }
        rows.append(row)
    
    return pd.DataFrame(rows)

def compute_framework_hash(candidate) -> str:
    """
    Create a hash of framework positions to detect near-duplicates.
    Uses vernier positions + a few other key FW positions.
    """
    key_positions = [42, 49, 50, 52, 66, 67, 68, 69, 71, 76, 78, 82, 87, 89, 91, 94]
    vals = []
    for pos in key_positions:
        aa = candidate.imgt_positions.get(pos, '-')
        vals.append(f"{pos}{aa}")
    return "|".join(vals)


def get_hallmark_string(candidate) -> str:
    """Get 4-char hallmark string from candidate."""
    return (
        candidate.imgt_positions.get(42, '?') +
        candidate.imgt_positions.get(49, '?') +
        candidate.imgt_positions.get(50, '?') +
        candidate.imgt_positions.get(52, '?')
    )


def select_with_quotas(
    candidates: List[VHHCandidate],
    n_total: int = 100,
    lead_vernier_positions: Optional[Dict[int, str]] = None,
    max_duplicates_per_hash: int = 3,
    prioritize_constrained: bool = True
) -> Tuple[List[VHHCandidate], Dict[str, Any]]:
    """
    Select top candidates with EQUAL distribution across hallmarks.
    
    v7.16 Update: 
    - Cysteine validation BEFORE deduplication (new)
    - Global deduplication by sequence (keeps best-scored duplicate)
    - All candidates get positive ranks (no negatives)
    - ranking_exempt column distinguishes controls from ranked
    
    Strategy:
    1. VALIDATE cysteines (reject odd-Cys, enforce C4 pair) - NEW in v7.16
    2. DEDUPLICATE globally by sequence (keep best ESM score)
    3. Separate ranking_exempt candidates (controls) - always included
    4. Group ranked candidates by hallmark
    5. Within each hallmark, sort by vernier similarity to lead (then by score)
    6. Select equal numbers from each hallmark (n_total / n_hallmarks)
    
    Args:
        candidates: List of scored VHHCandidate objects
        n_total: Total ranked candidates to select (exempt don't count)
        lead_vernier_positions: Lead's IMGT positions for vernier comparison
        max_duplicates_per_hash: Max candidates with identical framework hash
        prioritize_constrained: If True, YQRL gets slight priority
    
    Returns:
        (selected_candidates, selection_stats)
    """
    # =========================================================================
    # STEP 0: CYSTEINE VALIDATION (v7.16) - BEFORE DEDUP
    # =========================================================================
    # Reject candidates with invalid cysteine patterns before deduplication
    # This ensures we don't keep an invalid representative of a duplicate set
    
    candidates, cysteine_stats = filter_invalid_cysteines(candidates, verbose=True)
    
    # =========================================================================
    # STEP 1: GLOBAL DEDUPLICATION BY SEQUENCE (v7.14)
    # =========================================================================
    # Keep the candidate with the best combined_score among duplicates
    # Track what was removed for summary
    
    print(f"\n  Global deduplication...")
    seq_to_candidates = defaultdict(list)
    for c in candidates:
        seq_to_candidates[c.sequence].append(c)
    
    deduplicated = []
    duplicates_removed = []
    
    for seq, cands in seq_to_candidates.items():
        if len(cands) == 1:
            deduplicated.append(cands[0])
        else:
            # Sort by: is_lead first (always keep lead), then best score
            cands_sorted = sorted(cands, key=lambda c: (
                -1 if c.is_lead else 0,  # Lead always wins
                -c.scoring.combined_score if hasattr(c.scoring, 'combined_score') else 0
            ))
            kept = cands_sorted[0]
            removed = cands_sorted[1:]
            deduplicated.append(kept)
            
            # Track removed duplicates
            for r in removed:
                duplicates_removed.append({
                    'removed_id': r.id,
                    'kept_id': kept.id,
                    'sequence_hash': hash(seq) % 100000,
                    'removed_track': r.design_track,
                    'kept_track': kept.design_track,
                    'removed_hallmark': get_hallmark_string(r),
                    'kept_hallmark': get_hallmark_string(kept),
                })
    
    n_before = len(candidates)
    n_after = len(deduplicated)
    n_removed = n_before - n_after
    
    # Summarize by track (always initialize for stats)
    removed_by_track = defaultdict(int)
    for d in duplicates_removed:
        removed_by_track[d['removed_track']] += 1
    
    print(f"    Before: {n_before}, After: {n_after}, Removed: {n_removed} duplicates")
    if n_removed > 0:
        print(f"    Duplicates by track: {dict(removed_by_track)}")
    
    candidates = deduplicated
    
    # =========================================================================
    # STEP 2: SEPARATE BY TYPE
    # =========================================================================
    lead_candidates = [c for c in candidates if c.is_lead]
    exempt_candidates = [c for c in candidates if not c.is_lead and c.ranking_exempt]
    ranked_candidates = [c for c in candidates if not c.is_lead and not c.ranking_exempt]
    
    print(f"\n  Candidate pools (after dedup):")
    print(f"    Lead: {len(lead_candidates)}")
    print(f"    Exempt controls (Track 0-3): {len(exempt_candidates)}")
    print(f"    Ranked candidates (Track 4): {len(ranked_candidates)}")
    
    # v8.2: Group ranked candidates by (scaffold_type, hallmark)
    # This treats Universal as a separate "subfamily" alongside original hallmarks
    by_scaffold_hallmark = defaultdict(list)
    for c in ranked_candidates:
        hm = get_hallmark_string(c)
        scaffold = getattr(c, 'scaffold_type', 'original')
        # Create bucket key: "Universal" for universal scaffold, hallmark for original
        if scaffold == 'universal':
            bucket_key = f"Universal_{hm}"
        else:
            bucket_key = hm
        by_scaffold_hallmark[bucket_key].append(c)
    
    hallmarks_available = list(by_scaffold_hallmark.keys())
    n_hallmarks = len(hallmarks_available)
    
    if n_hallmarks == 0:
        # Return exempt candidates even if no ranked candidates
        for i, c in enumerate(exempt_candidates):
            c.rank = i + 1
        return exempt_candidates, {'error': 'No hallmarks found in ranked candidates', 'exempt_count': len(exempt_candidates)}
    
    # Calculate vernier similarity to lead for sorting
    VERNIER_POSITIONS = [66, 67, 68, 69, 71, 76, 78, 82, 87, 89, 91, 94]
    
    def vernier_similarity(c):
        """Count vernier positions matching lead."""
        if not lead_vernier_positions:
            return 0
        matches = 0
        for pos in VERNIER_POSITIONS:
            if c.imgt_positions.get(pos) == lead_vernier_positions.get(pos):
                matches += 1
        return matches
    
    # Sort each hallmark bucket by: combined_score ASC (best first)
    for hm in hallmarks_available:
        by_scaffold_hallmark[hm].sort(key=lambda c: c.scoring.combined_score)
    
    # Priority ordering: YQRL first (most constrained), Universal next, then others
    if prioritize_constrained:
        priority_order = []
        # YQRL first
        for h in hallmarks_available:
            if 'YQRL' in h:
                priority_order.append(h)
        # Universal next
        for h in hallmarks_available:
            if h.startswith('Universal_') and h not in priority_order:
                priority_order.append(h)
        # Then others
        for h in hallmarks_available:
            if h not in priority_order:
                priority_order.append(h)
        hallmarks_available = priority_order
    
    # Calculate quota per bucket (equal distribution)
    n_per_hallmark = n_total // n_hallmarks
    remainder = n_total % n_hallmarks
    
    print(f"  Equal scaffold/hallmark quotas: {n_per_hallmark} per bucket ({n_hallmarks} buckets)")
    print(f"  Buckets: {hallmarks_available}")
    
    selection_stats = {
        'n_hallmarks': n_hallmarks,
        'hallmarks_available': hallmarks_available,
        'target_per_hallmark': n_per_hallmark,
        'by_hallmark': {},
        'exempt_controls': {
            'total': len(exempt_candidates),
            'by_track': {},
        },
    }
    
    # Count exempt by track
    exempt_by_track = defaultdict(int)
    for c in exempt_candidates:
        exempt_by_track[c.design_track] += 1
    selection_stats['exempt_controls']['by_track'] = dict(exempt_by_track)
    
    selected_ranked = []
    hash_counts = defaultdict(int)
    
    # Select from each scaffold/hallmark bucket
    for i, hm in enumerate(hallmarks_available):
        quota = n_per_hallmark + (1 if i < remainder else 0)
        bucket = by_scaffold_hallmark[hm]
        bucket_selected = []
        
        for c in bucket:
            if len(bucket_selected) >= quota:
                break
            
            # Check duplicate cap
            fw_hash = compute_framework_hash(c)
            if hash_counts[fw_hash] >= max_duplicates_per_hash:
                continue
            
            bucket_selected.append(c)
            hash_counts[fw_hash] += 1
        
        selected_ranked.extend(bucket_selected)
        
        # Stats
        avg_score = sum(c.scoring.combined_score for c in bucket_selected) / max(len(bucket_selected), 1)
        selection_stats['by_hallmark'][hm] = {
            'available': len(bucket),
            'selected': len(bucket_selected),
            'target': quota,
            'avg_combined_score': round(avg_score, 4),
        }
        
        print(f"    {hm}: selected {len(bucket_selected)}/{quota} (available: {len(bucket)}, avg_score: {avg_score:.4f})")
    
    # Fill shortfall if any hallmark bucket was exhausted
    shortfall = n_total - len(selected_ranked)
    if shortfall > 0:
        print(f"  Shortfall: {shortfall} candidates, filling from remaining...")
        
        used_ids = {c.id for c in selected_ranked}
        remaining = [c for c in ranked_candidates if c.id not in used_ids]
        remaining.sort(key=lambda c: (-vernier_similarity(c), -c.scoring.combined_score))
        
        filled = 0
        for c in remaining:
            if filled >= shortfall:
                break
            
            fw_hash = compute_framework_hash(c)
            if hash_counts[fw_hash] >= max_duplicates_per_hash:
                continue
            
            selected_ranked.append(c)
            hash_counts[fw_hash] += 1
            filled += 1
        
        selection_stats['filled_from_overflow'] = filled
    
    # Sort ranked candidates by score for ranking
    selected_ranked.sort(key=lambda c: -c.scoring.combined_score)
    
    # =================================================================
    # v7.14: WITHIN-TRACK RANKING for exempt candidates
    # =================================================================
    # Different tracks have different objectives:
    #   Track 0 (minimal_hallmark): distance from input first (4mut > 5mut), then ESM
    #   Track 1 (single_vernier): ESM only (distance is constant ~5-6 mutations)
    #   Track 2 (paired_vernier): ESM first, tie-break by fewer mutations
    #   Track 3 (triplet_vernier): ESM only
    #
    # FALLBACK (v7.14): When ESM2 is unavailable, use deterministic ordering:
    #   n_mutations → framework_identity → generation_order
    # =================================================================
    
    def within_track_sort_key(c):
        """Generate sort key for within-track ranking.
        
        When ESM2 is unavailable (esm2_loss == 0), falls back to deterministic ordering
        using n_mutations, framework_identity, and generation_order.
        """
        track = c.design_track
        n_muts = len(c.mutations)
        # Fallback: if ESM2 unavailable, use 0 (will be sorted by other criteria)
        esm_available = c.scoring.esm2_loss > 0
        esm_loss = c.scoring.esm2_loss if esm_available else 0.0
        fw_identity = c.framework_identity_pct
        hallmark = get_hallmark_string(c)
        gen_order = c.generation_order  # Deterministic tie-breaker
        
        if track == "minimal_hallmark":
            # Primary: n_mutations (ascending) - 4mut beats 5mut
            # Secondary: ESM loss (ascending) if available, else fw_identity (descending)
            # Tertiary: generation_order for determinism
            if esm_available:
                return (0, hallmark, n_muts, esm_loss, -fw_identity, gen_order)
            else:
                return (0, hallmark, n_muts, -fw_identity, gen_order, 0)
        elif track == "single_vernier":
            # Primary: ESM loss (ascending) if available, else n_mutations (ascending)
            if esm_available:
                return (1, hallmark, esm_loss, -fw_identity, n_muts, gen_order)
            else:
                return (1, hallmark, n_muts, -fw_identity, gen_order, 0)
        elif track == "paired_vernier":
            # Primary: ESM loss (ascending) if available
            # Secondary: n_mutations (ascending)
            if esm_available:
                return (2, hallmark, esm_loss, n_muts, -fw_identity, gen_order)
            else:
                return (2, hallmark, n_muts, -fw_identity, gen_order, 0)
        elif track == "triplet_vernier":
            # Primary: ESM loss (ascending) if available
            if esm_available:
                return (3, hallmark, esm_loss, n_muts, -fw_identity, gen_order)
            else:
                return (3, hallmark, n_muts, -fw_identity, gen_order, 0)
        else:
            # Fallback for unknown tracks
            return (99, hallmark, n_muts, -fw_identity, gen_order, 0)
    
    # Sort exempt candidates with within-track ranking
    exempt_sorted = sorted(exempt_candidates, key=within_track_sort_key)
    
    # Assign track_rank (within-track local rank) - v7.14
    track_counters = defaultdict(lambda: defaultdict(int))  # track -> hallmark -> count
    for c in exempt_sorted:
        track = c.design_track
        hallmark = get_hallmark_string(c)
        track_counters[track][hallmark] += 1
        c.track_rank = track_counters[track][hallmark]  # 1-indexed within track+hallmark
    
    # v7.14: ALL candidates get POSITIVE ranks (no negatives!)
    # Lead: rank 0
    # Exempt controls: rank 1, 2, 3, ... 
    # Ranked candidates: continue from where controls left off
    # Use ranking_exempt column to distinguish controls from ranked
    
    current_rank = 1  # Start at 1 (lead is 0)
    
    # Exempt candidates get rank 1, 2, 3, ...
    for c in exempt_sorted:
        c.rank = current_rank
        current_rank += 1
    
    # Ranked candidates continue: rank N+1, N+2, ...
    selected_ranked.sort(key=lambda c: -c.scoring.combined_score)
    optimized_counters = defaultdict(int)  # hallmark -> count
    for c in selected_ranked:
        c.rank = current_rank
        current_rank += 1
        hallmark = get_hallmark_string(c)
        optimized_counters[hallmark] += 1
        c.track_rank = optimized_counters[hallmark]  # 1-indexed within optimized+hallmark
    
    # Combine: exempt first (as controls), then ranked
    selected = exempt_sorted + selected_ranked
    
    # Summary stats
    final_hallmark_counts = defaultdict(int)
    for c in selected_ranked:
        final_hallmark_counts[get_hallmark_string(c)] += 1
    selection_stats['final_hallmark_distribution'] = dict(final_hallmark_counts)
    selection_stats['n_exempt'] = len(exempt_sorted)
    selection_stats['n_ranked'] = len(selected_ranked)
    
    # v7.16: Add cysteine validation stats
    selection_stats['cysteine_filtering'] = cysteine_stats
    
    # v7.14: Add deduplication stats
    selection_stats['deduplication'] = {
        'n_before': n_before,
        'n_after': n_after,
        'n_removed': n_removed,
        'removed_details': duplicates_removed[:50] if len(duplicates_removed) > 50 else duplicates_removed,  # Cap at 50 for summary
        'removed_by_track': dict(removed_by_track) if n_removed > 0 else {},
    }
    
    # v7.14: Add within-track ranking stats
    track_summary = defaultdict(lambda: {'count': 0, 'hallmarks': set()})
    for c in exempt_sorted:
        track_summary[c.design_track]['count'] += 1
        track_summary[c.design_track]['hallmarks'].add(get_hallmark_string(c))
    selection_stats['track_summary'] = {
        k: {'count': v['count'], 'hallmarks': list(v['hallmarks'])} 
        for k, v in track_summary.items()
    }
    
    print(f"\n  Final selection: {len(exempt_sorted)} exempt controls + {len(selected_ranked)} ranked = {len(selected)} total")
    print(f"  Within-track breakdown:")
    for track in ['minimal_hallmark', 'single_vernier', 'paired_vernier', 'triplet_vernier']:
        if track in track_summary:
            info = track_summary[track]
            print(f"    {track}: {info['count']} candidates across {len(info['hallmarks'])} hallmarks")
    
    return selected, selection_stats


def build_enhanced_summary(
    all_candidates: List[VHHCandidate],
    selected: List[VHHCandidate],
    lead: Optional[VHHCandidate],
    hallmark_pool: List[str],
    target_families: List[str],
    generation_stats: Dict[str, Any],
    filter_stats: Dict[str, Any],
    scoring_config: Dict[str, Any],
    selection_stats: Dict[str, Any]
) -> Dict[str, Any]:
    """Build comprehensive summary JSON with full provenance."""
    from collections import Counter
    
    get_hm = lambda c: get_hallmark_string(c)
    
    summary = {
        "datetime": datetime.now().strftime("%Y%m%d_%H%M%S"),
        "version": "7.16",
        
        # Selection summary
        "selection": {
            "n_requested": len(selected) + (1 if lead else 0),
            "n_selected": len(selected),
            "n_leads": 1 if lead else 0,
            "n_original_framework": sum(1 for c in selected if c.framework_source == 'original'),
            "n_universal_graft": sum(1 for c in selected if c.framework_source == 'universal'),
            "selection_stats": selection_stats,
        },
        
        # Hallmark analysis
        "hallmarks": {
            "pool_available": hallmark_pool if hallmark_pool else [],
            "pool_size": len(hallmark_pool) if hallmark_pool else 0,
            "used_in_generation": dict(Counter(get_hm(c) for c in all_candidates if not c.is_lead)),
            "used_in_selection": dict(Counter(get_hm(c) for c in selected)),
        },
        
        # Family analysis
        "families": {
            "target_families": target_families,
            "detected_in_generation": dict(Counter(c.family for c in all_candidates if not c.is_lead)),
            "detected_in_selection": dict(Counter(c.family for c in selected)),
        },
        
        # Generation stats
        "generation": generation_stats,
        
        # Filter stats
        "filtering": filter_stats,
        
        # Scoring configuration
        "scoring": scoring_config,
        
        # Strategy breakdown in selected
        "strategies": dict(Counter(c.strategy for c in selected)),
        
        # Common violations
        "common_violations": _summarize_violations(selected),
    }
    
    return summary


def _summarize_violations(candidates: List[VHHCandidate], top_n: int = 10) -> List[Dict[str, Any]]:
    """Summarize most common rule violations across candidates."""
    from collections import Counter
    
    all_violations = []
    for c in candidates:
        all_violations.extend(c.scoring.top_violations)
    
    violation_counts = Counter(all_violations)
    return [
        {"violation": v, "count": n}
        for v, n in violation_counts.most_common(top_n)
    ]


def build_candidate_detail(candidate: VHHCandidate) -> Dict[str, Any]:
    """Build detailed record for JSONL output with rule provenance."""
    # Summarize rule firings
    comp_rules = getattr(candidate, 'applied_comp_rules', [])
    triplet_rules = getattr(candidate, 'applied_triplet_rules', [])
    
    def parse_rule_id(rule_str):
        """Extract key info from rule string."""
        parts = rule_str.split('|')
        if len(parts) >= 4:
            return {
                'type': parts[0],
                'family': parts[1] if len(parts) > 1 else '',
                'mutation': parts[3] if len(parts) > 3 else '',
            }
        return {'full': rule_str}
    
    rule_firings_summary = {
        'n_comp_rules': len(comp_rules),
        'n_triplet_rules': len(triplet_rules),
        'comp_rules_brief': [parse_rule_id(r) for r in comp_rules[:10]],
        'triplet_rules_brief': [parse_rule_id(r) for r in triplet_rules[:5]],
    }
    
    hallmarks = get_hallmark_string(candidate)
    
    return {
        'id': candidate.id,
        'rank': candidate.rank,
        'is_lead': candidate.is_lead,
        'sequence': candidate.sequence,
        'seq_length': len(candidate.sequence),
        'hallmarks': hallmarks,
        'detected_family': candidate.family,
        'target_family': candidate.target_family,
        'framework_source': candidate.framework_source,
        'strategy': candidate.strategy,
        'construction_method': candidate.construction_method,
        'n_mutations': len(candidate.mutations),
        'mutations': [str(m) for m in candidate.mutations],
        'framework_identity_pct': candidate.framework_identity_pct,
        
        # Scoring details
        'scoring': {
            'combined_score': candidate.scoring.combined_score,
            'esm2_loss': candidate.scoring.esm2_loss,
            'esm2_perplexity': candidate.scoring.esm2_perplexity,
            'plddt_mean': candidate.scoring.plddt_mean,
            'weighted_naturalness': candidate.scoring.weighted_naturalness,
            'rules_passed': candidate.scoring.rules_passed,
            'rules_total': candidate.scoring.rules_total,
            'rules_applicable': candidate.scoring.rules_applicable,
            'vernier_matches': candidate.scoring.vernier_matches,
            'vernier_total': candidate.scoring.vernier_total,
            'confidence_pct': candidate.scoring.confidence_pct,
            'family_probabilities': candidate.scoring.family_probabilities,
            'rule_compliance': candidate.scoring.rule_compliance,
        },
        
        # Provenance
        'applied_comp_rules': comp_rules,
        'applied_triplet_rules': triplet_rules,
        'rule_firings_summary': rule_firings_summary,
        'top_violations': candidate.scoring.top_violations,
        'generation_debug': candidate.generation_debug,
        
        # Regions
        'cdr1': candidate.cdr1,
        'cdr2': candidate.cdr2,
        'cdr3': candidate.cdr3,
        'fr1': candidate.fr1,
        'fr2': candidate.fr2,
        'fr3': candidate.fr3,
        'fr4': candidate.fr4,
    }


def save_results_enhanced(
    df, 
    selected: List[VHHCandidate],
    all_candidates: List[VHHCandidate],
    lead: Optional[VHHCandidate],
    output_dir: str, 
    base_name: str, 
    datetime_str: str,
    summary_data: Dict[str, Any]
):
    """Save results with enhanced summary and detailed JSONL."""
    os.makedirs(output_dir, exist_ok=True)
    
    # CSV
    csv_path = os.path.join(output_dir, f"{base_name}.csv")
    df.to_csv(csv_path, index=False)
    print(f"\nSaved: {csv_path}")
    
    # FASTA
    fasta_path = os.path.join(output_dir, f"{base_name}.fasta")
    with open(fasta_path, 'w') as f:
        for _, row in df.iterrows():
            score_str = "LEAD" if row['is_lead'] else f"score={row['combined_score']:.3f}"
            # Include zero-padded rank in header
            header = f">Rank{row['rank']}|{row['id']}|{row['framework_source']}|{score_str}|plddt={row['plddt_mean']:.1f}"
            f.write(f"{header}\n{row['sequence']}\n")
    print(f"Saved: {fasta_path}")
    
    # Enhanced summary JSON
    summary_path = os.path.join(output_dir, f"{base_name}_summary.json")
    with open(summary_path, 'w') as f:
        json.dump(summary_data, f, indent=2)
    print(f"Saved: {summary_path}")

    # Enhanced detailed JSONL with rule provenance
    detailed_path = os.path.join(output_dir, f"{base_name}_detailed.jsonl")
    with open(detailed_path, 'w') as f:
        # Add lead first if present
        if lead:
            detail = build_candidate_detail(lead)
            f.write(json.dumps(detail) + "\n")
        # Then selected candidates
        for c in selected:
            detail = build_candidate_detail(c)
            f.write(json.dumps(detail) + "\n")
    print(f"Saved: {detailed_path}")

# ============================================================
# MAIN
# ============================================================

def main():
    p = argparse.ArgumentParser(
        description="VHH Designer v7.5 - Fixed AntPack Parsing & IMGT Position Mapping",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
FIXES in v7.5:
  - Correct AntPack 4-tuple parsing + trim_alignment()
  - IMGT position-based lookup (no hardcoded FR2 indices)
  - Insertion-safe region extraction
  - Position dict {IMGT_pos -> AA} used throughout

Examples:
  python vhh_designer_v7_5_fixed.py -i input.fasta \\
      --rules analysis_rules_v7.json \\
      --archetypes analysis_vernier_archetypes_v7.json \\
      --n-generate 500 \\
      --n-select 92
"""
    )
    
    p.add_argument('--sequence', '-s', help='Input sequence')
    p.add_argument('--input', '-i', help='Input FASTA file')
    p.add_argument('--rules', '-r', required=True, help='analysis_rules_v7.json')
    p.add_argument('--archetypes', '-a', required=True, help='analysis_vernier_archetypes_v7.json')
    p.add_argument('--hallmark-db', dest='hallmark_db', default=None, help='Path to comprehensive_subfamily_analysis_imgt.xlsx (required for --target-hallmarks AUTO)')

    
    p.add_argument('--n-generate', type=int, default=500)
    p.add_argument('--n-select', type=int, default=198)
    p.add_argument('--target-hallmarks', '-t', default='FERG')
    p.add_argument('--target-families', nargs='+', default=['F_C2', 'Y_C2', 'Other_VHH'])
    # v8.0: mode removed - always generates BOTH original and Universal
    
    p.add_argument('--use-esm2', action='store_true', default=True)
    p.add_argument('--no-esm2', action='store_true')
    p.add_argument('--use-esmfold', action='store_true', default=True)
    p.add_argument('--esm-model', default='facebook/esm2_t6_8M_UR50D')
    
    p.add_argument('--esm2-weight', type=float, default=0.2, help='Weight for ESM2 in scoring')
    p.add_argument('--rule-weight', type=float, default=0.5, help='Weight for rules in scoring')
    p.add_argument('--no-esmfold', action='store_true', help='Disable ESMFold predictions')
    p.add_argument('--esmfold-weight', type=float, default=0.30, help='Weight for ESMFold in scoring')
    
    p.add_argument('--output-dir', '-o')
    p.add_argument('--name')
    
    args = p.parse_args()
    if str(args.target_hallmarks).upper() == "AUTO" and not args.hallmark_db:
        p.error("--hallmark-db is required when --target-hallmarks AUTO")

    
    use_esm2 = args.use_esm2 and not args.no_esm2
    use_esmfold = args.use_esmfold and not args.no_esmfold
    
    # Get input
    if args.sequence:
        input_seq = args.sequence.upper().replace(' ', '')
        seq_name = args.name or 'input'
    elif args.input:
        with open(args.input) as f:
            lines = f.readlines()
        seq_name = args.name
        input_seq = ''
        for line in lines:
            line = line.strip()
            if line.startswith('>'):
                if not seq_name:
                    seq_name = line[1:].split()[0]
            else:
                input_seq += line.upper()
        if not seq_name:
            seq_name = os.path.splitext(os.path.basename(args.input))[0]
    else:
        print("ERROR: Provide --sequence or --input")
        return
    
    print("=" * 70)
    print(f"VHH DESIGNER v{VERSION} - Probabilistic Consensus Sprinkling & Motif Tracks")
    print(f"  Output budget: {TOTAL_OUTPUT} total ({CONTROL_BUDGET} controls + {RANKED_BUDGET} ranked)")
    print("=" * 70)
    print(f"Sequence: {seq_name}")
    print(f"Input ({len(input_seq)} aa): {input_seq[:50]}...")
    print(f"\nGeneration: {args.n_generate} → Select top {args.n_select}")
    print(f"Mode: BOTH (original framework + Universal scaffold)")
    print(f"ESM2: {'enabled' if use_esm2 else 'disabled'}, ESMFold: {'enabled' if use_esmfold else 'disabled'}")
    
    # Extract CDRs (FIXED)
    print(f"\n{'='*70}")
    print("CDR EXTRACTION (Fixed)")
    print(f"{'='*70}")
    
    cdrs = extract_cdrs_fixed(input_seq)
    if not cdrs:
        print("ERROR: CDR extraction failed")
        return
    
    print(f"CDR1: {cdrs.cdr1}")
    print(f"CDR2: {cdrs.cdr2}")
    print(f"CDR3: {cdrs.cdr3} ({len(cdrs.cdr3)} aa)")
    
    if cdrs.imgt_numbered:
        print(f"Hallmarks (from IMGT positions): {cdrs.imgt_numbered.get_hallmarks()}")
        print(f"IMGT positions mapped: {len(cdrs.imgt_numbered.positions)}")
    
    cdr3_features = get_cdr3_features(cdrs.cdr3)
    
    # Generate
    print(f"\n{'='*70}")
    print("CANDIDATE GENERATION")
    print(f"{'='*70}")
    
    generator = CandidateGenerator(args.rules, args.archetypes, hallmark_db_path=args.hallmark_db)
    candidates = generator.generate(
        cdrs, args.n_generate,
        target_hallmarks=args.target_hallmarks,
        target_families=args.target_families
    )
    def is_canonical_vhh(pos):
        """True VHH filter: pos50=R required, pos52≠W, pos42 in {F,Y}
        
        Allows L at position 52 (e.g., YQRL - the most constrained scaffold).
        Position 50=R is the key VHH marker for solubility/independence.
        """
        a50 = pos.get(50, '')
        a52 = pos.get(52, '')
        a42 = pos.get(42, '')
        
        # CRITICAL: Require R at position 50 (VHH identity)
        if a50 != 'R':
            return False
        
        # Exclude W at position 52 (too VH-like)
        # Allow G, A, L, F (all functional as VHH)
        if a52 == 'W':
            return False
        
        # Require aromatic at position 42
        if a42 not in ('F', 'Y'):
            return False
        
        return True

    # Filter out non-VHH candidates (keep lead marked separately)
    before_filter = len(candidates)
    lead_candidate = [c for c in candidates if c.is_lead][0] if any(c.is_lead for c in candidates) else None
    non_lead = [c for c in candidates if not c.is_lead]
    
    # Track removal reasons for debugging
    removal_reasons = defaultdict(int)
    filtered = []
    for c in non_lead:
        a50 = c.imgt_positions.get(50, '')
        a52 = c.imgt_positions.get(52, '')
        a42 = c.imgt_positions.get(42, '')
        
        if a50 != 'R':
            removal_reasons['pos50_not_R'] += 1
            continue
        if a52 == 'W':
            removal_reasons['pos52_W_VHlike'] += 1
            continue
        if a42 not in ('F', 'Y'):
            removal_reasons[f'pos42_{a42}_not_FY'] += 1
            continue
        
        filtered.append(c)
    
    # Re-add lead (it keeps IGLW but is excluded from ranking)
    candidates = ([lead_candidate] if lead_candidate else []) + filtered
    print(f"  True VHH filter: {before_filter} → {len(candidates)} (removed {before_filter - len(candidates)})")
    if removal_reasons:
        print(f"  Removal reasons: {dict(removal_reasons)}")
    
    print(f"\nGenerated {len(candidates)} candidates")
    
    # Score
    print(f"\n{'='*70}")
    print("SCORING")
    print(f"{'='*70}")
    
    original_positions = cdrs.imgt_numbered.positions if cdrs.imgt_numbered else None
    
    scorer = CombinedScorer(
        args.rules, args.archetypes,
        esm_model=args.esm_model,
        use_esm2=use_esm2,
        use_esmfold=use_esmfold,
        esm2_weight=args.esm2_weight,
        esmfold_weight=args.esmfold_weight,
        rule_weight=args.rule_weight
    )
    
    candidates = scorer.score_candidates(candidates, cdr3_features, original_positions)
    
    # v7.14 Seatbelt: Validate immutable candidates weren't accidentally modified
    immutable_candidates = [c for c in candidates if c.is_immutable]
    if immutable_candidates:
        violations_found = False
        for c in immutable_candidates:
            is_valid, violations = validate_immutable_candidate(c)
            if not is_valid:
                print(f"  WARNING: Immutable candidate {c.id} has violations:")
                for v in violations:
                    print(f"    - {v}")
                violations_found = True
        if not violations_found:
            print(f"  ✓ All {len(immutable_candidates)} immutable candidates validated")
    
    # Rank
    print(f"\n{'='*70}")
    print("RANKING & SELECTION")
    print(f"{'='*70}")
    
    # Sort by score first
    candidates.sort(key=lambda c: (-c.scoring.combined_score,))
    
    # Equal hallmark quota selection with vernier similarity priority
    print(f"\nApplying equal-hallmark-quota selection (target: {args.n_select} total)")
    
    # Get lead's IMGT positions for vernier comparison
    lead_positions = lead_candidate.imgt_positions if lead_candidate else original_positions
    
    selected, selection_stats = select_with_quotas(
        candidates,
        n_total=args.n_select,
        lead_vernier_positions=lead_positions,
        max_duplicates_per_hash=3,
        prioritize_constrained=True  # Prioritize YQRL
    )
    
    # Add lead at rank 0 (excluded from quota selection)
    if lead_candidate:
        lead_candidate.rank = 0
        final_selected = [lead_candidate] + selected
    else:
        final_selected = selected
    
    print(f"\nSelected {len(selected)} candidates + {1 if lead_candidate else 0} lead")
    print(f"\nTop 15:")
    print("-" * 115)
    print(f"{'Rank':>5} {'Gen#':>4} {'ID':40} {'Score':>7} {'pLDDT':>6} {'Rules':>8} {'Vernier':>8} {'FW%':>5} {'Hallmarks':>9}")
    print("-" * 115)
    
    for c in final_selected[:15]:
        score_str = "LEAD" if c.is_lead else f"{c.scoring.combined_score:.3f}"
        rules_str = f"{c.scoring.rules_passed}/{c.scoring.rules_applicable}"
        vernier_str = f"{c.scoring.vernier_matches}/{c.scoring.vernier_total}"
        hallmarks = get_hallmarks_from_positions(c.imgt_positions)
        print(f"{c.rank:05d} {c.generation_order:4d} {c.id[:40]:40} {score_str:>7} {c.scoring.plddt_mean:6.1f} {rules_str:>8} {vernier_str:>8} {c.framework_identity_pct:5.1f} {hallmarks:>9}")
    
    # Build stats for enhanced summary
    filter_stats = {
        'gate_type': 'true_vhh',  # pos50=R, pos52≠W, pos42∈{F,Y}
        'n_removed': before_filter - len(candidates),
        'removal_rate_pct': round((before_filter - len(candidates)) / max(before_filter, 1) * 100, 1),
        'removal_reasons': dict(removal_reasons),
    }
    
    generation_stats = {
        'n_requested': args.n_generate,
        'n_generated_raw': before_filter,
        'n_after_filter': len(candidates),
    }
    
    scoring_config = {
        'esm2_weight': scorer.esm2_weight,
        'esmfold_weight': scorer.esmfold_weight,
        'rules_weight': scorer.rule_weight,
        'use_esm2': use_esm2,
        'use_esmfold': use_esmfold,
    }
    
    # Build enhanced summary
    summary_data = build_enhanced_summary(
        all_candidates=candidates,
        selected=selected,
        lead=lead_candidate,
        hallmark_pool=[],  # Would need to expose from generator
        target_families=args.target_families,
        generation_stats=generation_stats,
        filter_stats=filter_stats,
        scoring_config=scoring_config,
        selection_stats=selection_stats
    )
    
    # Save
    df = to_dataframe(final_selected)
    
    datetime_str = datetime.now().strftime('%Y%m%d_%H%M%S')
    version_str = f"v{VERSION.replace('.', '_')}"  # v8.3 -> v8_3
    n_generated = len(candidates)
    
    # Naming: [date]_[code]_[n_generated]_[leadID]
    # Example: 20260120_003719_v8_3_n20000_M69
    base_name = f"{datetime_str}_{version_str}_n{n_generated}_{seq_name}"
    folder_name = base_name
    
    if args.output_dir:
        output_dir = args.output_dir
    else:
        results_base = "results/designer_runs"
        output_dir = os.path.join(results_base, folder_name) if os.path.exists(results_base) else folder_name
    
    save_results_enhanced(
        df, selected, candidates, lead_candidate,
        output_dir, base_name, datetime_str,
        summary_data
    )
    
    # Also save all candidates CSV
    df_all = to_dataframe(candidates)
    all_csv = os.path.join(output_dir, f"{base_name}_all.csv")
    df_all.to_csv(all_csv, index=False)
    print(f"Saved: {all_csv}")
    
    # =========================================================================
    # AUTOMATIC MSA GENERATION (v8.3)
    # =========================================================================
    # Run csv2msa_antpack_v5.py on the ranked CSV
    ranked_csv = os.path.join(output_dir, f"{base_name}.csv")
    
    # Find csv2msa script (check common locations)
    msa_script = None
    script_name = "csv2msa_antpack_v5.py"
    search_paths = [
        os.path.dirname(os.path.abspath(__file__)),  # Same directory as this script
        os.getcwd(),  # Current working directory
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),  # Parent directory
    ]
    
    for search_dir in search_paths:
        candidate_path = os.path.join(search_dir, script_name)
        if os.path.exists(candidate_path):
            msa_script = candidate_path
            break
    
    if msa_script and os.path.exists(ranked_csv):
        print(f"\n{'='*70}")
        print("GENERATING MSA...")
        print(f"{'='*70}")
        print(f"Script: {msa_script}")
        print(f"Input: {ranked_csv}")
        
        try:
            import subprocess
            # Run with "n" as input to the initial question
            result = subprocess.run(
                [sys.executable, msa_script, "--table", ranked_csv, "--scheme", "imgt"],
                input="n\n",
                text=True,
                capture_output=True,
                timeout=300  # 5 minute timeout
            )
            
            if result.returncode == 0:
                print("MSA generation completed successfully!")
                if result.stdout:
                    # Print last few lines of output
                    lines = result.stdout.strip().split('\n')
                    for line in lines[-10:]:
                        print(f"  {line}")
            else:
                print(f"MSA generation failed with return code {result.returncode}")
                if result.stderr:
                    print(f"Error: {result.stderr[:500]}")
        except subprocess.TimeoutExpired:
            print("MSA generation timed out after 5 minutes")
        except Exception as e:
            print(f"MSA generation failed: {e}")
    elif not msa_script:
        print(f"\nNote: {script_name} not found. Skipping MSA generation.")
        print(f"To generate MSA manually, run:")
        print(f"  python {script_name} --table {ranked_csv} --scheme imgt")
    
    print(f"\n{'='*70}")
    print("DONE!")
    print(f"{'='*70}")

if __name__ == '__main__':
    main()