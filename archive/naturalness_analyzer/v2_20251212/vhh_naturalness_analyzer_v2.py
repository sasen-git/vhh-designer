#!/usr/bin/env python3
"""
VHH CDR-Framework Naturalness Analyzer
======================================

Uses camelid epistasis data as a PRIOR on CDR3-framework compatibility,
NOT as a binary fold/don't fold gate.

Key insight: For humanized scaffolds, foldability is known from experiments.
The question is: "Given this framework, how natural is this CDR3?"

Two modes:
1. Full sequence input: Extracts CDRs and frameworks automatically
2. Framework + CDR input: Uses provided framework regions with grafted CDRs

Output:
- Naturalness score (lower = more natural)
- Per-feature z-scores (length, charge, aromatics)
- Framework annotation (classical VHH, humanized, VH-like)
- Risk assessment based on CDR-framework compatibility
- Suggested optimizations

Usage:
  # Full sequences
  python vhh_naturalness_analyzer.py --input sequences.csv --epistasis epistasis_v2_full.pkl
  
  # Framework library + CDRs
  python vhh_naturalness_analyzer.py --input designs.csv --frameworks VHH_frameworks.xlsx --epistasis epistasis_v2_full.pkl
"""

import os
import re
import json
import pickle
import argparse
import numpy as np
import pandas as pd
from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Optional, Any
from collections import Counter, defaultdict

# ============================================================================
# CONSTANTS
# ============================================================================

VALID_AA = set('ACDEFGHIKLMNPQRSTVWY')
AROMATIC = set('FWY')
BASIC = set('KRH')
ACIDIC = set('DE')
HYDROPHOBIC = set('AILMFVW')
POLAR = set('STNQ')

HYDROPHOBICITY = {
    'I': 4.5, 'V': 4.2, 'L': 3.8, 'F': 2.8, 'C': 2.5, 'M': 1.9, 'A': 1.8,
    'G': -0.4, 'T': -0.7, 'S': -0.8, 'W': -0.9, 'Y': -1.3, 'P': -1.6,
    'H': -3.2, 'E': -3.5, 'Q': -3.5, 'D': -3.5, 'N': -3.5, 'K': -3.9, 'R': -4.5,
}

# ============================================================================
# DATA CLASSES
# ============================================================================

@dataclass
class FrameworkProfile:
    """Profile of a framework's Vernier positions."""
    name: str
    fr1: str
    fr2: str
    fr3: str
    fr4: str
    
    # Hallmark positions (IMGT numbering)
    pos42: str = '-'  # FR2[1] - Y/F in VHH, V in VH
    pos49: str = '-'  # FR2[8] - E/Q in VHH, G in VH
    pos50: str = '-'  # FR2[9] - R in VHH, L in VH
    pos52: str = '-'  # FR2[11] - G/L/F in VHH, W in VH
    
    # Full Vernier pattern
    vernier_pattern: Dict[str, str] = field(default_factory=dict)
    
    # Classification
    framework_type: str = 'Unknown'  # 'Classical_VHH', 'Humanized', 'VH_like', 'Mixed'
    hallmark_signature: str = ''
    
    def extract_positions(self):
        """Extract hallmark and Vernier positions from FR sequences."""
        fr2 = self.fr2
        fr3 = self.fr3
        
        # Hallmark positions from FR2
        if len(fr2) > 1:
            self.pos42 = fr2[1]
        if len(fr2) > 8:
            self.pos49 = fr2[8]
        if len(fr2) > 9:
            self.pos50 = fr2[9]
        if len(fr2) > 11:
            self.pos52 = fr2[11]
        
        self.hallmark_signature = f"{self.pos42}-{self.pos49}-{self.pos50}-{self.pos52}"
        
        # Full Vernier pattern
        self.vernier_pattern = {}
        
        # FR2 positions
        for i, pos_name in [(1, 'FR2_2'), (8, 'FR2_9'), (9, 'FR2_10'), 
                            (11, 'FR2_12'), (13, 'FR2_14')]:
            if len(fr2) > i:
                self.vernier_pattern[pos_name] = fr2[i]
        
        # FR3 positions
        for i, pos_name in [(5, 'FR3_6'), (7, 'FR3_8'), (11, 'FR3_12'), 
                            (15, 'FR3_16'), (21, 'FR3_22')]:
            if len(fr3) > i:
                self.vernier_pattern[pos_name] = fr3[i]
        
        # Classify framework type
        self._classify_framework()
    
    def _classify_framework(self):
        """Classify the framework type based on hallmark positions."""
        pos42_vhh = self.pos42 in ['F', 'Y']
        pos49_vhh = self.pos49 in ['E', 'Q']
        pos50_vhh = self.pos50 == 'R'
        pos52_vhh = self.pos52 in ['G', 'L', 'F', 'A']
        
        pos42_vh = self.pos42 == 'V'
        pos49_vh = self.pos49 == 'G'
        pos50_vh = self.pos50 == 'L'
        pos52_vh = self.pos52 == 'W'
        
        # Classical VHH: all four hallmarks are VHH-like
        if pos42_vhh and pos49_vhh and pos50_vhh and pos52_vhh:
            self.framework_type = 'Classical_VHH'
        
        # True VH: V at 42 AND W at 52
        elif pos42_vh and pos52_vh:
            self.framework_type = 'VH'
        
        # Humanized VHH: F/Y at 42, G/L/F at 52, but G at 49 or L at 50
        elif pos42_vhh and pos52_vhh and (pos49_vh or pos50_vh):
            self.framework_type = 'Humanized_VHH'
        
        # VH-like: V at 42 but not W at 52
        elif pos42_vh:
            self.framework_type = 'VH_like'
        
        # Mixed/Other
        else:
            self.framework_type = 'Mixed'


@dataclass
class CDRFeatures:
    """Features of a CDR3 sequence."""
    sequence: str
    length: int = 0
    charge: float = 0.0
    aromatic_count: int = 0
    hydrophobicity: float = 0.0
    cysteine_count: int = 0
    glycine_count: int = 0
    proline_count: int = 0
    extraction_warning: str = ''
    
    def compute(self):
        """Compute all features from the sequence."""
        seq = self.sequence.upper()
        self.length = len(seq)
        self.charge = sum(1 for aa in seq if aa in BASIC) - sum(1 for aa in seq if aa in ACIDIC)
        self.aromatic_count = sum(1 for aa in seq if aa in AROMATIC)
        self.cysteine_count = seq.count('C')
        self.glycine_count = seq.count('G')
        self.proline_count = seq.count('P')
        
        if self.length > 0:
            self.hydrophobicity = sum(HYDROPHOBICITY.get(aa, 0) for aa in seq) / self.length
        
        # Validation checks
        if 'WFRQ' in seq or 'WYRQ' in seq or 'WVRQ' in seq:
            self.extraction_warning = "CDR3 contains FR2 anchor (WFRQ/WYRQ/WVRQ) - likely extraction error!"
        elif self.length > 30:
            self.extraction_warning = f"CDR3 unusually long ({self.length} aa) - verify extraction"
        elif self.length < 5:
            self.extraction_warning = f"CDR3 unusually short ({self.length} aa) - verify extraction"


@dataclass
class NaturalnessResult:
    """Result of naturalness assessment."""
    # Core scores
    naturalness_score: float  # Lower = more natural (sum of squared z-scores)
    z_length: float
    z_charge: float
    z_aromatic: float
    
    # Classification
    naturalness_class: str  # 'very_natural', 'natural', 'acceptable', 'unusual', 'outlier'
    
    # Best matching cluster info
    best_cluster_n: int
    best_cluster_pattern: Dict[str, str]
    best_cluster_cdr3_mean: float
    best_cluster_cdr3_std: float
    pattern_match_score: int  # How many Vernier positions match
    
    # Risk assessment
    risk_level: str  # 'low', 'medium', 'high'
    risk_factors: List[str] = field(default_factory=list)
    
    # Suggestions
    suggested_mutations: List[Tuple[str, str, str, str]] = field(default_factory=list)


@dataclass
class AnalysisResult:
    """Complete analysis result for one design."""
    id: str
    framework_name: str
    cdr3_sequence: str
    
    # Framework info
    framework: FrameworkProfile
    cdr3_features: CDRFeatures
    
    # Naturalness assessment
    naturalness: NaturalnessResult
    
    # Summary
    overall_assessment: str
    recommendation: str


# ============================================================================
# VERNIER MODEL
# ============================================================================

class VernierModel:
    """Model for assessing CDR3 naturalness based on camelid Vernier patterns."""
    
    def __init__(self, clusters: Dict, compensation_stats: Dict):
        self.clusters = clusters
        self.compensation_stats = compensation_stats
        self.cluster_list = self._build_cluster_list()
        
        # Pre-compute global statistics for fallback
        self._compute_global_stats()
    
    def _build_cluster_list(self) -> List[Dict]:
        """Convert cluster dict to list with parsed patterns."""
        result = []
        
        for pattern_str, data in self.clusters.items():
            if not isinstance(data, dict):
                continue
            
            pattern = data.get('pattern', {})
            cdr3_stats = data.get('cdr3_length', {})
            charge_stats = data.get('cdr3_charge', {})
            arom_stats = data.get('cdr3_aromatic', {})
            
            n = cdr3_stats.get('n', 0) if isinstance(cdr3_stats, dict) else 0
            
            if n < 50:  # Skip tiny clusters
                continue
            
            result.append({
                'pattern': pattern,
                'n': n,
                'length_mean': cdr3_stats.get('mean', 15) if isinstance(cdr3_stats, dict) else 15,
                'length_std': max(cdr3_stats.get('std', 5), 1) if isinstance(cdr3_stats, dict) else 5,
                'charge_mean': charge_stats.get('mean', 0) if isinstance(charge_stats, dict) else 0,
                'charge_std': max(charge_stats.get('std', 3), 0.5) if isinstance(charge_stats, dict) else 3,
                'aromatic_mean': arom_stats.get('mean', 2) if isinstance(arom_stats, dict) else 2,
                'aromatic_std': max(arom_stats.get('std', 1.5), 0.5) if isinstance(arom_stats, dict) else 1.5,
            })
        
        return result
    
    def _compute_global_stats(self):
        """Compute global CDR3 statistics across all clusters."""
        lengths = []
        charges = []
        aroms = []
        
        for cl in self.cluster_list:
            for _ in range(min(cl['n'], 1000)):  # Weight by n, cap at 1000
                lengths.append(cl['length_mean'])
                charges.append(cl['charge_mean'])
                aroms.append(cl['aromatic_mean'])
        
        if lengths:
            self.global_length_mean = np.mean(lengths)
            self.global_length_std = np.std(lengths)
            self.global_charge_mean = np.mean(charges)
            self.global_charge_std = np.std(charges)
            self.global_aromatic_mean = np.mean(aroms)
            self.global_aromatic_std = np.std(aroms)
        else:
            # Fallback defaults
            self.global_length_mean = 15
            self.global_length_std = 5
            self.global_charge_mean = 0
            self.global_charge_std = 3
            self.global_aromatic_mean = 2
            self.global_aromatic_std = 1.5
    
    @classmethod
    def from_epistasis_file(cls, filepath: str) -> 'VernierModel':
        """Load model from epistasis pickle file."""
        with open(filepath, 'rb') as f:
            data = pickle.load(f)
        
        clusters = data.get('analysis_2_vernier_clusters', {})
        comp_stats = data.get('analysis_1_compensation', {}).get('stats', {})
        
        return cls(clusters, comp_stats)
    
    def score_naturalness(self, framework: FrameworkProfile, 
                          cdr3: CDRFeatures,
                          min_cluster_n: int = 100) -> NaturalnessResult:
        """
        Score how natural this CDR3 is for the given framework.
        
        Uses the Vernier pattern to find similar natural clusters,
        then computes z-scores for CDR3 features.
        """
        design_pattern = framework.vernier_pattern
        
        # Find clusters with similar Vernier patterns
        candidates = []
        
        for cluster in self.cluster_list:
            if cluster['n'] < min_cluster_n:
                continue
            
            # Count how many Vernier positions match
            matches = 0
            total_positions = 0
            
            for pos in ['FR2_2', 'FR2_9', 'FR2_10', 'FR2_12', 'FR2_14',
                        'FR3_6', 'FR3_8', 'FR3_12', 'FR3_16']:
                if pos in design_pattern and pos in cluster['pattern']:
                    total_positions += 1
                    if design_pattern[pos] == cluster['pattern'][pos]:
                        matches += 1
            
            if total_positions > 0:
                match_score = matches
                candidates.append((cluster, match_score))
        
        # Sort by match score (descending) and n (descending)
        candidates.sort(key=lambda x: (-x[1], -x[0]['n']))
        
        # If no good matches, use global statistics
        if not candidates or candidates[0][1] < 2:
            # Fall back to global statistics
            z_len = (cdr3.length - self.global_length_mean) / self.global_length_std
            z_chg = (cdr3.charge - self.global_charge_mean) / self.global_charge_std
            z_aro = (cdr3.aromatic_count - self.global_aromatic_mean) / self.global_aromatic_std
            
            score = z_len**2 + z_chg**2 + 0.5 * z_aro**2
            
            return NaturalnessResult(
                naturalness_score=score,
                z_length=z_len,
                z_charge=z_chg,
                z_aromatic=z_aro,
                naturalness_class=self._classify_score(score),
                best_cluster_n=0,
                best_cluster_pattern={},
                best_cluster_cdr3_mean=self.global_length_mean,
                best_cluster_cdr3_std=self.global_length_std,
                pattern_match_score=0,
                risk_level=self._assess_risk(score, z_len, z_chg, z_aro),
                risk_factors=["No matching Vernier clusters - using global statistics"],
                suggested_mutations=[]
            )
        
        # Use top matching clusters (weighted by n and match score)
        # Take top 5 and compute weighted average
        top_clusters = candidates[:5]
        
        weights = []
        z_lengths = []
        z_charges = []
        z_aromatics = []
        
        for cluster, match_score in top_clusters:
            # Weight by match score and sqrt(n)
            w = match_score * np.sqrt(cluster['n'])
            weights.append(w)
            
            z_len = (cdr3.length - cluster['length_mean']) / cluster['length_std']
            z_chg = (cdr3.charge - cluster['charge_mean']) / cluster['charge_std']
            z_aro = (cdr3.aromatic_count - cluster['aromatic_mean']) / cluster['aromatic_std']
            
            z_lengths.append(z_len)
            z_charges.append(z_chg)
            z_aromatics.append(z_aro)
        
        # Weighted average of z-scores
        total_weight = sum(weights)
        z_len = sum(z * w for z, w in zip(z_lengths, weights)) / total_weight
        z_chg = sum(z * w for z, w in zip(z_charges, weights)) / total_weight
        z_aro = sum(z * w for z, w in zip(z_aromatics, weights)) / total_weight
        
        score = z_len**2 + z_chg**2 + 0.5 * z_aro**2
        
        best_cluster, best_match = top_clusters[0]
        
        # Generate risk factors
        risk_factors = []
        if abs(z_len) > 2:
            if z_len > 0:
                risk_factors.append(f"CDR3 longer than typical ({cdr3.length} vs {best_cluster['length_mean']:.0f} aa)")
            else:
                risk_factors.append(f"CDR3 shorter than typical ({cdr3.length} vs {best_cluster['length_mean']:.0f} aa)")
        
        if abs(z_chg) > 2:
            if z_chg > 0:
                risk_factors.append(f"CDR3 more positive than typical (charge {cdr3.charge:+.0f})")
            else:
                risk_factors.append(f"CDR3 more negative than typical (charge {cdr3.charge:+.0f})")
        
        if abs(z_aro) > 2:
            risk_factors.append(f"Unusual aromatic content ({cdr3.aromatic_count} aromatics)")
        
        if framework.framework_type == 'VH':
            risk_factors.append("VH framework - may require VL pairing")
        
        # Generate mutation suggestions
        suggestions = self._suggest_mutations(design_pattern, best_cluster['pattern'], 
                                              framework.framework_type)
        
        return NaturalnessResult(
            naturalness_score=score,
            z_length=z_len,
            z_charge=z_chg,
            z_aromatic=z_aro,
            naturalness_class=self._classify_score(score),
            best_cluster_n=best_cluster['n'],
            best_cluster_pattern=best_cluster['pattern'],
            best_cluster_cdr3_mean=best_cluster['length_mean'],
            best_cluster_cdr3_std=best_cluster['length_std'],
            pattern_match_score=best_match,
            risk_level=self._assess_risk(score, z_len, z_chg, z_aro),
            risk_factors=risk_factors,
            suggested_mutations=suggestions
        )
    
    def _classify_score(self, score: float) -> str:
        """Classify the naturalness score."""
        if score < 1:
            return 'very_natural'
        elif score < 3:
            return 'natural'
        elif score < 6:
            return 'acceptable'
        elif score < 12:
            return 'unusual'
        else:
            return 'outlier'
    
    def _assess_risk(self, score: float, z_len: float, z_chg: float, z_aro: float) -> str:
        """Assess overall risk level."""
        if score < 3 and all(abs(z) < 2 for z in [z_len, z_chg, z_aro]):
            return 'low'
        elif score < 8 and all(abs(z) < 3 for z in [z_len, z_chg, z_aro]):
            return 'medium'
        else:
            return 'high'
    
    def _suggest_mutations(self, design_pattern: Dict, cluster_pattern: Dict,
                           framework_type: str) -> List[Tuple[str, str, str, str]]:
        """
        Suggest framework mutations to improve naturalness.
        
        Only suggests changes to FR3 positions for humanized frameworks
        (FR2 hallmarks are fixed by design).
        """
        suggestions = []
        
        # For humanized frameworks, don't suggest changing FR2 hallmarks
        skip_positions = set()
        if framework_type in ['Humanized_VHH', 'VH', 'VH_like']:
            skip_positions = {'FR2_2', 'FR2_9', 'FR2_10', 'FR2_12'}
        
        for pos in ['FR2_14', 'FR3_6', 'FR3_8', 'FR3_12', 'FR3_16', 'FR3_22']:
            if pos in skip_positions:
                continue
            
            if pos in design_pattern and pos in cluster_pattern:
                design_aa = design_pattern[pos]
                cluster_aa = cluster_pattern[pos]
                
                if design_aa != cluster_aa:
                    reason = f"Most similar natural clusters use {cluster_aa} at this position"
                    suggestions.append((pos, design_aa, cluster_aa, reason))
        
        return suggestions[:3]  # Top 3 suggestions


# ============================================================================
# SEQUENCE EXTRACTION
# ============================================================================

def extract_cdr3_from_sequence(sequence: str) -> str:
    """
    Extract CDR3 from a full VHH sequence.
    
    CRITICAL: Must find CDR3 AFTER FR2, not just any C...W pattern.
    The C in CDR1 + W in FR2(WFRQ) can be mistaken for CDR3 if not careful.
    """
    seq = sequence.upper().replace('-', '').replace('.', '')
    
    # Step 1: Find FR2 anchor (WFRQ or W.RQ pattern)
    fr2_match = re.search(r'W[YFVILA]RQ', seq[25:60])
    if not fr2_match:
        fr2_match = re.search(r'W.RQ', seq[25:60])
    
    if fr2_match:
        # FR2 starts at this position
        fr2_start = 25 + fr2_match.start()
        # CDR3 must be AFTER FR2 + CDR2 + FR3, typically position 90+
        search_start = fr2_start + 50  # Skip FR2 (14) + CDR2 (~17) + part of FR3 (~20)
    else:
        search_start = 80  # Fallback
    
    # Step 2: Look for CDR3 pattern C...WG AFTER FR2
    # CDR3 typically starts around position 95-105
    cdr3_region = seq[search_start:]
    
    # Find C...WG pattern in the CDR3 region
    cdr3_match = re.search(r'C([A-Z]{3,35})W[GS]', cdr3_region)
    if cdr3_match:
        return 'C' + cdr3_match.group(1)
    
    # Fallback: look for conserved C around position 95-110
    for c_pos in range(max(search_start, 90), min(len(seq) - 5, 115)):
        if seq[c_pos] == 'C':
            # Look for WG within 5-35 residues
            w_pos = seq.find('W', c_pos + 5, c_pos + 40)
            if w_pos > c_pos and seq[w_pos:w_pos+2] in ['WG', 'WS']:
                return seq[c_pos:w_pos]
    
    # Last fallback: find last C before WG at end
    wg_match = re.search(r'W[GS][A-Z]*$', seq[-30:])
    if wg_match:
        wg_pos = len(seq) - 30 + wg_match.start()
        # Find C before WG
        c_pos = seq.rfind('C', wg_pos - 35, wg_pos)
        if c_pos > 0:
            return seq[c_pos:wg_pos]
    
    return ''


def extract_frameworks_from_sequence(sequence: str) -> Optional[FrameworkProfile]:
    """
    Extract framework regions from a full VHH sequence.
    
    Layout: FR1 - CDR1 - FR2 - CDR2 - FR3 - CDR3 - FR4
    """
    seq = sequence.upper().replace('-', '').replace('.', '')
    
    if len(seq) < 100:
        return None
    
    # Find FR2 anchor (WFRQ or W.RQ)
    fr2_match = re.search(r'W[YFVILA]RQ', seq[25:60])
    if not fr2_match:
        fr2_match = re.search(r'W.RQ', seq[25:60])
    
    if not fr2_match:
        return None
    
    fr2_start = 25 + fr2_match.start()
    fr2_end = fr2_start + 14
    
    # Find CDR3: Must be AFTER FR2+CDR2+FR3 (typically position 90+)
    search_start = fr2_end + 40  # Skip CDR2 (~17) + part of FR3 (~25)
    
    # Look for C...WG pattern in CDR3 region
    cdr3_region = seq[search_start:]
    cdr3_match = re.search(r'C([A-Z]{3,35})W[GS]', cdr3_region)
    
    if not cdr3_match:
        # Fallback: look for WG near end and C before it
        wg_match = re.search(r'W[GS]QG', seq[-20:])
        if wg_match:
            wg_pos = len(seq) - 20 + wg_match.start()
            c_pos = seq.rfind('C', wg_pos - 35, wg_pos)
            if c_pos > 0:
                cdr3_start = c_pos
                cdr3_end = wg_pos
            else:
                return None
        else:
            return None
    else:
        cdr3_start = search_start + cdr3_match.start()
        cdr3_end = search_start + cdr3_match.end() - 2  # Exclude WG
    
    # FR3 is between CDR2 and CDR3
    # CDR2 starts ~15-17 residues after FR2
    cdr2_start = fr2_end
    cdr2_approx_len = 17
    fr3_start = cdr2_start + cdr2_approx_len
    fr3_end = cdr3_start
    
    profile = FrameworkProfile(
        name='Extracted',
        fr1=seq[:fr2_start],
        fr2=seq[fr2_start:fr2_end],
        fr3=seq[fr3_start:fr3_end] if fr3_end > fr3_start else '',
        fr4=seq[cdr3_end:] if cdr3_end < len(seq) else ''
    )
    
    profile.extract_positions()
    return profile


# ============================================================================
# MAIN ANALYSIS
# ============================================================================

def analyze_design(design_id: str, 
                   framework: FrameworkProfile,
                   cdr3_seq: str,
                   model: VernierModel) -> AnalysisResult:
    """
    Analyze a single design for CDR3-framework naturalness.
    """
    # Compute CDR3 features
    cdr3 = CDRFeatures(sequence=cdr3_seq)
    cdr3.compute()
    
    # Score naturalness
    naturalness = model.score_naturalness(framework, cdr3)
    
    # Generate overall assessment
    if naturalness.naturalness_class in ['very_natural', 'natural']:
        assessment = f"GOOD: CDR3 is natural for this framework context"
    elif naturalness.naturalness_class == 'acceptable':
        assessment = f"OK: CDR3 is somewhat unusual but within acceptable range"
    elif naturalness.naturalness_class == 'unusual':
        assessment = f"CAUTION: CDR3 is unusual for this framework - higher risk"
    else:
        assessment = f"WARNING: CDR3 is an outlier - may not work well with this framework"
    
    # Add framework context
    if framework.framework_type == 'Humanized_VHH':
        assessment += f" (humanized scaffold: {framework.hallmark_signature})"
    elif framework.framework_type == 'VH':
        assessment += f" (VH framework - typically requires VL pairing)"
    
    # Generate recommendation
    if naturalness.risk_level == 'low':
        recommendation = "Proceed - natural-like CDR3 for this framework"
    elif naturalness.risk_level == 'medium':
        recommendation = "Proceed with monitoring - some deviation from natural patterns"
    else:
        if naturalness.suggested_mutations:
            recommendation = f"Consider FR mutations: {', '.join(m[0]+':'+m[1]+'→'+m[2] for m in naturalness.suggested_mutations)}"
        else:
            recommendation = "High risk - consider CDR3 redesign or alternative framework"
    
    return AnalysisResult(
        id=design_id,
        framework_name=framework.name,
        cdr3_sequence=cdr3_seq,
        framework=framework,
        cdr3_features=cdr3,
        naturalness=naturalness,
        overall_assessment=assessment,
        recommendation=recommendation
    )


def load_frameworks(frameworks_file: str) -> Dict[str, FrameworkProfile]:
    """Load framework library from Excel file."""
    df = pd.read_excel(frameworks_file)
    
    frameworks = {}
    for _, row in df.iterrows():
        name = row.get('Name', row.get('name', 'Unknown'))
        
        profile = FrameworkProfile(
            name=name,
            fr1=str(row.get('FR1', '')),
            fr2=str(row.get('FR2', '')),
            fr3=str(row.get('FR3', '')),
            fr4=str(row.get('FR4', ''))
        )
        profile.extract_positions()
        frameworks[name] = profile
    
    return frameworks


def analyze_batch(input_file: str, 
                  epistasis_file: str,
                  frameworks_file: Optional[str] = None,
                  output_file: Optional[str] = None) -> pd.DataFrame:
    """
    Analyze a batch of designs.
    
    Input modes:
    1. If frameworks_file provided: Input should have 'framework' and 'cdr3' columns
    2. Otherwise: Input should have 'sequence' column (full VHH)
    """
    # Load epistasis model
    print(f"Loading epistasis model from {epistasis_file}...")
    model = VernierModel.from_epistasis_file(epistasis_file)
    print(f"Loaded {len(model.cluster_list)} Vernier clusters")
    
    # Load frameworks if provided
    framework_lib = {}
    if frameworks_file:
        print(f"\nLoading frameworks from {frameworks_file}...")
        framework_lib = load_frameworks(frameworks_file)
        print(f"Loaded {len(framework_lib)} frameworks:")
        for name, fw in framework_lib.items():
            print(f"  {name}: {fw.hallmark_signature} ({fw.framework_type})")
    
    # Load input
    print(f"\nLoading designs from {input_file}...")
    if input_file.endswith('.xlsx') or input_file.endswith('.xls'):
        df = pd.read_excel(input_file)
    else:
        df = pd.read_csv(input_file)
    
    print(f"Found {len(df)} designs")
    print(f"Columns: {list(df.columns)}")
    
    # Determine input mode
    has_framework_col = any(col.lower() in ['framework', 'fw', 'framework_name'] for col in df.columns)
    has_cdr3_col = any(col.lower() in ['cdr3', 'cdr3_seq', 'cdr3_sequence'] for col in df.columns)
    has_sequence_col = any(col.lower() in ['sequence', 'seq', 'aa_sequence'] for col in df.columns)
    
    # Find ID column
    id_col = None
    for col in df.columns:
        if col.lower() in ['id', 'name', 'design_id', 'sequence_id']:
            id_col = col
            break
    if id_col is None:
        df['_id'] = [f"design_{i}" for i in range(len(df))]
        id_col = '_id'
    
    # Analyze each design
    results = []
    
    for idx, row in df.iterrows():
        design_id = str(row[id_col])
        
        try:
            if has_framework_col and has_cdr3_col and framework_lib:
                # Mode 1: Framework library + CDR3
                fw_name = None
                cdr3_seq = None
                
                for col in df.columns:
                    if col.lower() in ['framework', 'fw', 'framework_name']:
                        fw_name = str(row[col])
                    if col.lower() in ['cdr3', 'cdr3_seq', 'cdr3_sequence']:
                        cdr3_seq = str(row[col])
                
                if fw_name not in framework_lib:
                    print(f"  Warning: Unknown framework '{fw_name}' for {design_id}")
                    continue
                
                framework = framework_lib[fw_name]
                
            elif has_sequence_col:
                # Mode 2: Full sequence
                seq_col = None
                for col in df.columns:
                    if col.lower() in ['sequence', 'seq', 'aa_sequence']:
                        seq_col = col
                        break
                
                sequence = str(row[seq_col])
                framework = extract_frameworks_from_sequence(sequence)
                cdr3_seq = extract_cdr3_from_sequence(sequence)
                
                if framework is None:
                    print(f"  Warning: Could not extract framework from {design_id}")
                    continue
            
            else:
                print(f"  Error: Cannot determine input mode for {design_id}")
                continue
            
            # Run analysis
            result = analyze_design(design_id, framework, cdr3_seq, model)
            results.append(result)
            
        except Exception as e:
            print(f"  Error analyzing {design_id}: {e}")
            continue
        
        if (idx + 1) % 100 == 0:
            print(f"  Processed {idx + 1}/{len(df)} designs...")
    
    # Convert to DataFrame
    output_data = []
    for r in results:
        row = {
            'id': r.id,
            'framework_name': r.framework_name,
            'framework_type': r.framework.framework_type,
            'hallmark_signature': r.framework.hallmark_signature,
            
            'cdr3_sequence': r.cdr3_sequence,
            'cdr3_length': r.cdr3_features.length,
            'cdr3_charge': r.cdr3_features.charge,
            'cdr3_aromatic': r.cdr3_features.aromatic_count,
            
            'naturalness_score': r.naturalness.naturalness_score,
            'naturalness_class': r.naturalness.naturalness_class,
            'z_length': r.naturalness.z_length,
            'z_charge': r.naturalness.z_charge,
            'z_aromatic': r.naturalness.z_aromatic,
            
            'risk_level': r.naturalness.risk_level,
            'pattern_match_score': r.naturalness.pattern_match_score,
            'best_cluster_n': r.naturalness.best_cluster_n,
            'best_cluster_cdr3_mean': r.naturalness.best_cluster_cdr3_mean,
            
            'risk_factors': '; '.join(r.naturalness.risk_factors),
            'suggested_mutations': '; '.join(f"{m[0]}:{m[1]}→{m[2]}" for m in r.naturalness.suggested_mutations),
            
            'overall_assessment': r.overall_assessment,
            'recommendation': r.recommendation,
        }
        output_data.append(row)
    
    result_df = pd.DataFrame(output_data)
    
    # Print summary
    print(f"\n{'='*60}")
    print("ANALYSIS SUMMARY")
    print(f"{'='*60}")
    print(f"Total designs analyzed: {len(result_df)}")
    
    print(f"\nNaturalness distribution:")
    for cls, count in result_df['naturalness_class'].value_counts().items():
        pct = 100 * count / len(result_df)
        print(f"  {cls}: {count} ({pct:.1f}%)")
    
    print(f"\nRisk level distribution:")
    for risk, count in result_df['risk_level'].value_counts().items():
        pct = 100 * count / len(result_df)
        print(f"  {risk}: {count} ({pct:.1f}%)")
    
    print(f"\nFramework type distribution:")
    for fw_type, count in result_df['framework_type'].value_counts().items():
        pct = 100 * count / len(result_df)
        print(f"  {fw_type}: {count} ({pct:.1f}%)")
    
    # Save output
    if output_file is None:
        base = os.path.splitext(input_file)[0]
        output_file = f"{base}_naturalness.xlsx"
    
    result_df.to_excel(output_file, index=False)
    print(f"\nResults saved to: {output_file}")
    
    return result_df


# ============================================================================
# CLI
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='VHH CDR-Framework Naturalness Analyzer',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Full sequences
  python vhh_naturalness_analyzer.py -i sequences.csv -e epistasis_v2_full.pkl
  
  # Framework library + CDRs  
  python vhh_naturalness_analyzer.py -i designs.csv -f VHH_frameworks.xlsx -e epistasis_v2_full.pkl

Input formats:
  Mode 1 (with frameworks file):
    - Input: id, framework, cdr3
    - Frameworks: Name, FR1, FR2, FR3, FR4
    
  Mode 2 (full sequences):
    - Input: id, sequence
        """
    )
    
    parser.add_argument('--input', '-i', required=True,
                        help='Input CSV or XLSX file with designs')
    parser.add_argument('--epistasis', '-e', required=True,
                        help='Epistasis pickle file (epistasis_v2_full.pkl)')
    parser.add_argument('--frameworks', '-f', default=None,
                        help='Framework library XLSX (optional)')
    parser.add_argument('--output', '-o', default=None,
                        help='Output XLSX file')
    
    args = parser.parse_args()
    
    analyze_batch(args.input, args.epistasis, args.frameworks, args.output)


if __name__ == '__main__':
    main()
