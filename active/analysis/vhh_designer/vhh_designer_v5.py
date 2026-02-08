#!/usr/bin/env python3
"""
VHH Designer v5 - Protected Framework Rule Application
==============================================

Converts a mouse heavy chain to VHH candidates using proper rule ordering:

  For UNIVERSAL mode (grafting CDRs onto humanized scaffold):
    1. Start with Universal scaffold (FR2 = WFRQAPGQGLEAVA, hallmarks = FGLA)
    2. Apply safe rules (FR3, FR4 positions) - FR2 core is PROTECTED
    3. Hallmark positions (IMGT 42,49,50,52) only changed via explicit patterns
    
  For ORIGINAL mode (keeping input FRs with VHH-izing mutations):
    PASS 1: HALLMARK mutations (IMGT 42,49,50,52 → VHH patterns like FERG, YERL)
    PASS 2: VERNIER archetype (family-specific FR pattern for stability)
    PASS 3: COMPENSATION/EPISTASIS rules (CDR-conditioned fine-tuning)

Family Classification (based on IMGT positions):
  - F_C2: pos42=F, pos49=E, pos50=R, 2 cysteines
  - Y_C2: pos42=Y, pos49=E, pos50=R, 2 cysteines  
  - F_C4: pos42=F, pos49=E, pos50=R, 4 cysteines
  - Y_C4: pos42=Y, pos49=E, pos50=R, 4 cysteines
  - VH_like: pos50=L (human VH-like)
  - Non_classical: other combinations

CRITICAL: Substring indices (FR2_4, FR3_6) are 1-based and NOT equivalent to IMGT positions.
This script maps them correctly and protects FR2 core in universal mode.

Usage:
  python vhh_designer_v5.py -i M69.fasta \\
      --compensation correlation_results_v3_compensation.pkl \\
      --epistasis epistasis_v2_full.pkl \\
      --mode original
"""

import os
import re
import pickle
import argparse
import numpy as np
import pandas as pd
from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Optional, Set
from collections import defaultdict
from datetime import datetime

try:
    from antpack import SingleChainAnnotator
    ANTPACK_AVAILABLE = True
except ImportError:
    ANTPACK_AVAILABLE = False
    print("Warning: AntPack not available. Install: pip install antpack")

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

# IMGT numbering for universal scaffold FR2
# FR2: W-F-R-Q-A-P-G-Q-G-L-E-A-V-A
# IMGT: 39-40-41-42-43-44-45-46-47-48-49-50-51-52
UNIVERSAL_SCAFFOLD_IMGT = {
    'FR2': {
        0: 39, 1: 40, 2: 41, 3: 42, 4: 43, 5: 44, 6: 45, 7: 46,
        8: 47, 9: 48, 10: 49, 11: 50, 12: 51, 13: 52
    },
    # FR2 index -> IMGT position
    'FR2_to_IMGT': {0: 39, 1: 40, 2: 41, 3: 42, 4: 43, 5: 44, 6: 45, 7: 46,
                    8: 47, 9: 48, 10: 49, 11: 50, 12: 51, 13: 52},
    # Substrate index -> IMGT (for epistasis rules like FR2_4)
    'FR2_substr_to_IMGT': {1: 39, 2: 40, 3: 41, 4: 42, 5: 43, 6: 44, 7: 45, 8: 46,
                           9: 47, 10: 48, 11: 49, 12: 50, 13: 51, 14: 52},
}

# Protected positions in universal scaffold - DO NOT MUTATE except for hallmarks
# These are the FR2 core motif (WFRQAPGQGLEAVA) and key structural positions
PROTECTED_IMGT_POSITIONS = {
    # FR2 core motif - only mutate via explicit hallmark patterns
    39, 40, 41,  # W-F-R (start of FR2)
    43, 44, 45, 46, 47, 48,  # A-P-G-Q-G-L (middle of FR2)
    # Hallmark positions are handled separately - 42, 49, 50, 52
}

# Hallmark/interface positions - only mutate with explicit hallmark patterns
HALLMARK_IMGT_POSITIONS = {42, 49, 50, 52}

# VHH Hallmark mutations defined by IMGT positions (42, 49, 50, 52)
VHH_HALLMARKS_IMGT = {
    'FERG': {42: 'F', 49: 'E', 50: 'R', 52: 'G'},  # Classic llama
    'YERL': {42: 'Y', 49: 'E', 50: 'R', 52: 'L'},  # Classic alpaca
    'FERF': {42: 'F', 49: 'E', 50: 'R', 52: 'F'},  # Variant
    'FERA': {42: 'F', 49: 'E', 50: 'R', 52: 'A'},  # Variant
    'FGLA': {42: 'F', 49: 'G', 50: 'L', 52: 'A'},  # Humanized (Universal scaffold)
}

# Mapping from substring index (1-based, as used in epistasis) to IMGT position
# This is critical for correctly applying rules to scaffolds
SUBSTR_TO_IMGT = {
    'FR2': {
        # FR2_1 through FR2_14 (1-based) -> IMGT 39-52
        1: 39, 2: 40, 3: 41, 4: 42, 5: 43, 6: 44, 7: 45, 8: 46,
        9: 47, 10: 48, 11: 49, 12: 50, 13: 51, 14: 52,
    },
    'FR3': {
        # FR3 positions - approximate mapping (depends on CDR2 length)
        # Standard VHH FR3 starts at IMGT 66
        1: 66, 2: 67, 3: 68, 4: 69, 5: 70, 6: 71, 7: 72, 8: 73,
        9: 74, 10: 75, 11: 76, 12: 77, 13: 78, 14: 79, 15: 80, 16: 81,
        17: 82, 18: 82, 19: 83, 20: 84, 21: 85, 22: 86, 23: 87, 24: 88,
        25: 89, 26: 90, 27: 91, 28: 92, 29: 93, 30: 94, 31: 95, 32: 96,
    },
}

AA_CATEGORIES = {
    'hydrophobic': set('AILMFVW'), 'polar': set('STNQ'),
    'basic': set('KRH'), 'acidic': set('DE'),
    'aromatic': set('FWY'), 'special': set('CGP'), 'reactive': set('CY'),
}

def categorize_aa(aa):
    for cat, aas in AA_CATEGORIES.items():
        if aa.upper() in aas:
            return cat
    return 'other'

# ============================================================
# DATA CLASSES
# ============================================================

@dataclass
class CDRSet:
    cdr1: str
    cdr2: str
    cdr3: str
    fr1: str = ''
    fr2: str = ''
    fr3: str = ''
    fr4: str = ''
    imgt_numbering: dict = None
    
    def get_aa_at_imgt(self, imgt_pos):
        """Get amino acid at a specific IMGT position."""
        if self.imgt_numbering and imgt_pos in self.imgt_numbering:
            return self.imgt_numbering[imgt_pos]['aa']
        return None
    
    def get_region_index_at_imgt(self, imgt_pos):
        """Get (region, index) for a specific IMGT position."""
        if self.imgt_numbering and imgt_pos in self.imgt_numbering:
            info = self.imgt_numbering[imgt_pos]
            return info['region'], info['index']
        return None, None
    
    def get_full_sequence(self):
        """Reconstruct full sequence."""
        return self.fr1 + self.cdr1 + self.fr2 + self.cdr2 + self.fr3 + self.cdr3 + self.fr4

@dataclass
class Mutation:
    position: str
    from_aa: str
    to_aa: str
    rule: str
    confidence: float
    source: str  # 'compensation', 'epistasis', 'vernier', 'hallmark'
    family: str = ''
    condition: str = ''
    
    def __str__(self):
        return f"{self.position}:{self.from_aa}->{self.to_aa}"

@dataclass 
class VHHCandidate:
    id: str
    rank: int
    sequence: str
    framework_source: str
    cdr1: str
    cdr2: str
    cdr3: str
    fr1: str
    fr2: str
    fr3: str
    fr4: str
    mutations: List[Mutation]
    family: str = ''
    naturalness_score: float = 0.0
    confidence: float = 0.0
    cdr3_preserved: bool = True
    strategy: str = ''
    warnings: List[str] = field(default_factory=list)

# ============================================================
# CDR EXTRACTION
# ============================================================

_annotator = None
def get_annotator():
    global _annotator
    if _annotator is None and ANTPACK_AVAILABLE:
        _annotator = SingleChainAnnotator(chains=['H'], scheme='imgt')
    return _annotator

def extract_cdrs(sequence):
    """Extract CDRs and FRs using AntPack with IMGT numbering."""
    if not ANTPACK_AVAILABLE:
        print("ERROR: AntPack required")
        return None
    
    annotator = get_annotator()
    seq = sequence.upper().replace('-', '').replace('.', '').replace(' ', '')
    
    try:
        numbering, pid, chain, err = annotator.analyze_seq(seq)
        if numbering is None or len(numbering) == 0:
            print(f"  ANARCI failed: {err}")
            return None
        
        trimmed_seq, trimmed_nums, _, _ = annotator.trim_alignment(seq, (numbering, pid, chain, err))
        labels = annotator.assign_cdr_labels(trimmed_nums, chain)
        
        label_map = {
            'fmwk1': 'FR1', 'fmwk2': 'FR2', 'fmwk3': 'FR3', 'fmwk4': 'FR4',
            'cdr1': 'CDR1', 'cdr2': 'CDR2', 'cdr3': 'CDR3',
            'fwh1': 'FR1', 'fwh2': 'FR2', 'fwh3': 'FR3', 'fwh4': 'FR4',
            'cdrh1': 'CDR1', 'cdrh2': 'CDR2', 'cdrh3': 'CDR3',
        }
        
        regions = {'FR1': '', 'CDR1': '', 'FR2': '', 'CDR2': '', 'FR3': '', 'CDR3': '', 'FR4': ''}
        region_indices = {'FR1': 0, 'CDR1': 0, 'FR2': 0, 'CDR2': 0, 'FR3': 0, 'CDR3': 0, 'FR4': 0}
        imgt_numbering = {}
        
        for i, (aa, label, imgt_pos) in enumerate(zip(trimmed_seq, labels, trimmed_nums)):
            normalized = label_map.get(label, label_map.get(label.lower(), None))
            if normalized and normalized in regions:
                if isinstance(imgt_pos, tuple):
                    imgt_num = imgt_pos[0]
                    insertion = imgt_pos[1].strip() if len(imgt_pos) > 1 else ''
                else:
                    imgt_num = imgt_pos
                    insertion = ''
                
                # Convert to int if possible
                try:
                    imgt_num_int = int(imgt_num)
                except (ValueError, TypeError):
                    imgt_num_int = None
                
                # Store with multiple key formats for lookup flexibility
                info = {
                    'region': normalized,
                    'index': region_indices[normalized],
                    'aa': aa
                }
                
                # Store with integer key (primary)
                if imgt_num_int is not None and not insertion:
                    imgt_numbering[imgt_num_int] = info
                
                # Also store with string key
                if insertion:
                    imgt_numbering[f"{imgt_num}{insertion}"] = info
                else:
                    imgt_numbering[str(imgt_num)] = info
                
                regions[normalized] += aa
                region_indices[normalized] += 1
        
        return CDRSet(
            cdr1=regions['CDR1'], cdr2=regions['CDR2'], cdr3=regions['CDR3'],
            fr1=regions['FR1'], fr2=regions['FR2'], fr3=regions['FR3'], fr4=regions['FR4'],
            imgt_numbering=imgt_numbering
        )
    except Exception as e:
        print(f"  CDR extraction error: {e}")
        return None

# ============================================================
# FAMILY CLASSIFICATION
# ============================================================

def classify_family(cdrs: CDRSet, full_sequence: str = '') -> str:
    """
    Classify VHH family based on IMGT hallmark positions.
    Falls back to FR2 pattern matching if IMGT lookup fails.
    
    Returns one of: F_C2, Y_C2, F_C4, Y_C4, VH_like, Non_classical, Classical_other, Unknown
    """
    pos42 = cdrs.get_aa_at_imgt(42)
    pos49 = cdrs.get_aa_at_imgt(49)
    pos50 = cdrs.get_aa_at_imgt(50)
    pos52 = cdrs.get_aa_at_imgt(52)
    
    # Count cysteines
    if full_sequence:
        n_cys = full_sequence.count('C')
    else:
        n_cys = cdrs.get_full_sequence().count('C')
    
    # Fallback: if IMGT lookup failed, try to extract from FR2 directly
    # FR2 positions: IMGT 39-46 is approximately FR2[0:8] for standard VHH
    # Position 42 is ~FR2[3], position 49 is ~FR2[10], position 50 is ~FR2[11], position 52 is ~FR2[13]
    if pos42 is None and len(cdrs.fr2) >= 14:
        # Standard mapping for VHH FR2:
        # W-Y/F-R-Q-A-P-G-Q/K-G-L/R-E-A/R/F-V-A
        # 39-40-41-42-43-44-45-46-47-48-49-50-51-52-53
        # FR2 index: 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14
        pos42 = cdrs.fr2[3] if len(cdrs.fr2) > 3 else None
        pos49 = cdrs.fr2[10] if len(cdrs.fr2) > 10 else None
        pos50 = cdrs.fr2[11] if len(cdrs.fr2) > 11 else None
        pos52 = cdrs.fr2[13] if len(cdrs.fr2) > 13 else None
    
    # If still no data, use pattern matching on FR2
    if pos42 is None and cdrs.fr2:
        # Look for classic VHH patterns in FR2
        fr2 = cdrs.fr2.upper()
        
        # Classic llama/alpaca: WF/YRQ or WF/YRQA pattern
        if 'FRQ' in fr2:
            pos42 = 'F'
        elif 'YRQ' in fr2:
            pos42 = 'Y'
        
        # Check for E (pos49) and R (pos50) - classic VHH markers
        if 'ER' in fr2 or 'GRER' in fr2 or 'GLEA' in fr2:
            # This looks like: G-L/R-E-A or similar VHH pattern
            pos49 = 'E'
            pos50 = 'R'
        elif 'GL' in fr2:
            # Human-like: G-L (pos49=G, pos50=L)
            pos50 = 'L'
    
    # Set defaults for missing values
    pos42 = pos42 or '-'
    pos49 = pos49 or '-'
    pos50 = pos50 or '-'
    pos52 = pos52 or '-'
    
    # VH-like: pos50 = L (human-like)
    if pos50 == 'L':
        return 'VH_like'
    
    # Classical VHH criteria
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
    
    # W at position 52
    if pos52 == 'W':
        return 'VHH_W52'
    
    # Check if it's likely a mouse/human VH (not camelid)
    # Mouse VH typically has: pos42=W, pos49=G, pos50=L
    if pos50 != 'R' and pos42 not in ['F', 'Y']:
        return 'VH_like'  # Treat as VH-like for humanization purposes
    
    return 'Non_classical'

def get_cdr3_category(cdr3: str) -> dict:
    """Get CDR3 categories for epistasis rule matching."""
    length = len(cdr3)
    charge = sum(1 for a in cdr3 if a in 'KRH') - sum(1 for a in cdr3 if a in 'DE')
    n_cys = cdr3.count('C')
    
    # Length categories
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
    
    # Charge categories
    if charge < -2:
        charge_cat = 'very_neg'
    elif charge < -0.5:
        charge_cat = 'negative'
    elif charge < 1:
        charge_cat = 'neutral'
    elif charge < 3:
        charge_cat = 'positive'
    else:
        charge_cat = 'very_pos'
    
    # Cysteine category
    cys_cat = 'with_cys' if n_cys > 0 else 'no_cys'
    
    return {
        'cdr3_len': len_cat,
        'cdr3_charge': charge_cat,
        'cdr3_cys': cys_cat,
        'length': length,
        'charge': charge,
        'n_cys': n_cys,
    }

# ============================================================
# RULE ENGINES
# ============================================================

class CompensationRulesEngine:
    """
    PASS 1: CDR features → FR residue rules.
    Source: correlation_results_v3_compensation.pkl
    """
    
    def __init__(self, filepath: str):
        with open(filepath, 'rb') as f:
            self.data = pickle.load(f)
        
        self.multi_position_rules = self.data.get('multi_position_rules', [])
        self.vernier_archetypes = self.data.get('vernier_archetypes', {})
        
        print(f"  Compensation: {len(self.multi_position_rules)} rules, {len(self.vernier_archetypes)} families")
    
    def get_cdr_features(self, cdrs: CDRSet) -> dict:
        """Extract CDR position features for rule matching."""
        features = {}
        for i, aa in enumerate(cdrs.cdr1):
            features[f'cdr1[{i}]'] = aa
        for i in range(1, min(4, len(cdrs.cdr1)+1)):
            features[f'cdr1[{-i}]'] = cdrs.cdr1[-i] if i <= len(cdrs.cdr1) else ''
        
        for i, aa in enumerate(cdrs.cdr2):
            features[f'cdr2[{i}]'] = aa
        for i in range(1, min(4, len(cdrs.cdr2)+1)):
            features[f'cdr2[{-i}]'] = cdrs.cdr2[-i] if i <= len(cdrs.cdr2) else ''
        
        for i, aa in enumerate(cdrs.cdr3):
            features[f'cdr3[{i}]'] = aa
        for i in range(1, min(4, len(cdrs.cdr3)+1)):
            features[f'cdr3[{-i}]'] = cdrs.cdr3[-i] if i <= len(cdrs.cdr3) else ''
        
        return features
    
    def find_applicable_rules(self, cdrs: CDRSet, family: str, min_confidence: float = 75.0) -> List[dict]:
        """Find compensation rules that match CDR features and family."""
        features = self.get_cdr_features(cdrs)
        applicable = []
        
        for rule in self.multi_position_rules:
            rule_family = rule.get('family', '')
            
            # Filter by family if specified
            if rule_family and rule_family != family:
                # Also allow rules for "catch-all" families
                if rule_family not in ['Other_VHH', 'Unknown']:
                    continue
            
            for pred in rule.get('top_predictors', []):
                conf = pred.get('conf', 0)
                if conf < min_confidence:
                    continue
                
                # Parse rule like 'cdr1[2]=M'
                match = re.match(r'(cdr\d+)\[(-?\d+)\]=([A-Z])', pred.get('rule', ''))
                if match:
                    cdr, pos, expected = match.groups()
                    feature_key = f'{cdr}[{int(pos)}]'
                    
                    if features.get(feature_key) == expected:
                        applicable.append({
                            'fw_position': rule['fw_position'],
                            'suggested_aa': rule['fw_residue'],
                            'confidence': conf,
                            'rule': pred['rule'],
                            'family': rule_family,
                            'source': 'compensation',
                        })
        
        # Deduplicate: keep highest confidence per position
        by_pos = {}
        for r in applicable:
            pos = r['fw_position']
            if pos not in by_pos or r['confidence'] > by_pos[pos]['confidence']:
                by_pos[pos] = r
        
        return sorted(by_pos.values(), key=lambda x: -x['confidence'])


class EpistasisRulesEngine:
    """
    PASS 2: FR-FR interactions given CDR context.
    Source: epistasis_v2_full.pkl → analysis_4_higher_order_rules
    
    Memory-optimized: only loads rules and vernier clusters, discards the rest.
    """
    
    def __init__(self, filepath: str):
        # Load only what we need to avoid memory issues with large files
        print(f"  Loading epistasis from {filepath}...")
        with open(filepath, 'rb') as f:
            data = pickle.load(f)
        
        # Extract only what we need
        self.higher_order_rules = data.get('analysis_4_higher_order_rules', [])
        self.vernier_clusters = data.get('analysis_2_vernier_clusters', {})
        
        # Explicitly delete the large data to free memory
        del data
        import gc
        gc.collect()
        
        print(f"  Epistasis: {len(self.higher_order_rules)} rules, {len(self.vernier_clusters)} clusters")
    
    def find_applicable_rules(self, cdrs: CDRSet, family: str, current_frs: dict, 
                               min_confidence: float = 70.0) -> List[dict]:
        """
        Find epistasis rules matching family, CDR features, and current FR state.
        
        current_frs: {'FR1': 'ABC...', 'FR2': 'DEF...', ...}
        """
        cdr3_cats = get_cdr3_category(cdrs.cdr3)
        applicable = []
        
        for rule in self.higher_order_rules:
            # Filter by family
            rule_family = rule.get('family', '')
            if rule_family and rule_family != family:
                continue
            
            # Filter by confidence
            conf = rule.get('confidence', 0) * 100  # Convert to percentage
            if conf < min_confidence:
                continue
            
            # Parse condition: "cdr3_len=short AND FR2_12=C"
            condition = rule.get('condition', '')
            if not self._check_condition(condition, cdr3_cats, current_frs):
                continue
            
            # Parse result: "FR3_6=V"
            result = rule.get('result', '')
            match = re.match(r'(FR\d)_(\d+)=([A-Z])', result)
            if not match:
                # Try IMGT format
                match = re.match(r'IMGT(\d+)=([A-Z])', result)
                if match:
                    imgt_num, aa = match.groups()
                    fw_position = f'IMGT{imgt_num}'
                else:
                    continue
            else:
                region, idx, aa = match.groups()
                fw_position = f'{region}[{idx}]'
            
            applicable.append({
                'fw_position': fw_position,
                'suggested_aa': aa,
                'confidence': conf,
                'rule': condition,
                'family': rule_family,
                'source': 'epistasis',
                'condition': condition,
            })
        
        # Deduplicate: keep highest confidence per position
        by_pos = {}
        for r in applicable:
            pos = r['fw_position']
            if pos not in by_pos or r['confidence'] > by_pos[pos]['confidence']:
                by_pos[pos] = r
        
        return sorted(by_pos.values(), key=lambda x: -x['confidence'])
    
    def _check_condition(self, condition: str, cdr3_cats: dict, current_frs: dict) -> bool:
        """Check if condition matches current state."""
        if not condition:
            return True
        
        parts = condition.split(' AND ')
        for part in parts:
            part = part.strip()
            if '=' not in part:
                continue
            
            key, value = part.split('=', 1)
            key = key.strip()
            value = value.strip()
            
            # CDR3 category check
            if key in cdr3_cats:
                if cdr3_cats[key] != value:
                    return False
            # FR position check (e.g., FR2_12=C)
            elif key.startswith('FR'):
                match = re.match(r'(FR\d)_(\d+)', key)
                if match:
                    region, idx = match.groups()
                    idx = int(idx)
                    fr_seq = current_frs.get(region, '')
                    if idx < len(fr_seq):
                        if fr_seq[idx] != value:
                            return False
                    else:
                        return False
        
        return True


class EpistasisRulesEngineFromData:
    """
    Same as EpistasisRulesEngine but initialized from already-loaded data.
    Used to avoid loading the large file twice.
    """
    
    def __init__(self, higher_order_rules: list, vernier_clusters: dict):
        self.higher_order_rules = higher_order_rules
        self.vernier_clusters = vernier_clusters
        print(f"  Epistasis: {len(self.higher_order_rules)} rules, {len(self.vernier_clusters)} clusters")
    
    def find_applicable_rules(self, cdrs: CDRSet, family: str, current_frs: dict, 
                               min_confidence: float = 70.0) -> List[dict]:
        """
        Find epistasis rules matching family, CDR features, and current FR state.
        """
        cdr3_cats = get_cdr3_category(cdrs.cdr3)
        applicable = []
        
        for rule in self.higher_order_rules:
            # Filter by family
            rule_family = rule.get('family', '')
            if rule_family and rule_family != family:
                continue
            
            # Filter by confidence
            conf = rule.get('confidence', 0) * 100  # Convert to percentage
            if conf < min_confidence:
                continue
            
            # Parse condition
            condition = rule.get('condition', '')
            if not self._check_condition(condition, cdr3_cats, current_frs):
                continue
            
            # Parse result
            result = rule.get('result', '')
            match = re.match(r'(FR\d)_(\d+)=([A-Z])', result)
            if not match:
                match = re.match(r'IMGT(\d+)=([A-Z])', result)
                if match:
                    imgt_num, aa = match.groups()
                    fw_position = f'IMGT{imgt_num}'
                else:
                    continue
            else:
                region, idx, aa = match.groups()
                fw_position = f'{region}[{idx}]'
            
            applicable.append({
                'fw_position': fw_position,
                'suggested_aa': aa,
                'confidence': conf,
                'rule': condition,
                'family': rule_family,
                'source': 'epistasis',
                'condition': condition,
            })
        
        # Deduplicate
        by_pos = {}
        for r in applicable:
            pos = r['fw_position']
            if pos not in by_pos or r['confidence'] > by_pos[pos]['confidence']:
                by_pos[pos] = r
        
        return sorted(by_pos.values(), key=lambda x: -x['confidence'])
    
    def _check_condition(self, condition: str, cdr3_cats: dict, current_frs: dict) -> bool:
        """Check if condition matches current state."""
        if not condition:
            return True
        
        parts = condition.split(' AND ')
        for part in parts:
            part = part.strip()
            if '=' not in part:
                continue
            
            key, value = part.split('=', 1)
            key = key.strip()
            value = value.strip()
            
            if key in cdr3_cats:
                if cdr3_cats[key] != value:
                    return False
            elif key.startswith('FR'):
                match = re.match(r'(FR\d)_(\d+)', key)
                if match:
                    region, idx = match.groups()
                    idx = int(idx)
                    fr_seq = current_frs.get(region, '')
                    if idx < len(fr_seq):
                        if fr_seq[idx] != value:
                            return False
                    else:
                        return False
        
        return True


class VernierEngine:
    """
    PASS 3: Family-specific Vernier archetype patterns.
    """
    
    def __init__(self, compensation_data: dict, epistasis_data: dict):
        # Prefer epistasis vernier (from full analysis) over compensation
        self.archetypes = epistasis_data.get('analysis_2_vernier_clusters', {})
        self.comp_archetypes = compensation_data.get('vernier_archetypes', {})
        
        # Merge: epistasis first, then compensation for missing families
        for family, rules in self.comp_archetypes.items():
            if family not in self.archetypes:
                self.archetypes[family] = rules
    
    def get_family_pattern(self, family: str) -> List[dict]:
        """Get Vernier pattern for a family."""
        # Find best matching cluster for family
        family_clusters = []
        
        for cluster_key, cluster_data in self.archetypes.items():
            if not isinstance(cluster_data, dict):
                continue
            
            families = cluster_data.get('families', {})
            if family in families:
                support = families[family]
                pattern = cluster_data.get('pattern', {})
                family_clusters.append({
                    'support': support,
                    'pattern': pattern,
                    'cluster_key': cluster_key,
                })
        
        if not family_clusters:
            return []
        
        # Use cluster with highest support for this family
        best = max(family_clusters, key=lambda x: x['support'])
        pattern = best['pattern']
        
        # Convert pattern to rules
        rules = []
        for pos_name, aa in pattern.items():
            if aa and aa != '-':
                rules.append({
                    'fw_position': pos_name,
                    'suggested_aa': aa,
                    'confidence': 80.0,  # Default vernier confidence
                    'rule': f'Vernier_{family}',
                    'family': family,
                    'source': 'vernier',
                })
        
        return rules


class NaturalnessScorer:
    """Score candidates based on CDR3 statistics from database."""
    
    def __init__(self, epistasis_data: dict):
        self.clusters = epistasis_data.get('analysis_2_vernier_clusters', {})
        
        lengths, charges = [], []
        for cluster in self.clusters.values():
            if isinstance(cluster, dict):
                cdr3_stats = cluster.get('cdr3_length', {})
                n = cdr3_stats.get('n', 0) if isinstance(cdr3_stats, dict) else 0
                if n >= 100:
                    lengths.extend([cdr3_stats.get('mean', 15)] * min(n, 500))
                    charge_stats = cluster.get('cdr3_charge', {})
                    if isinstance(charge_stats, dict):
                        charges.extend([charge_stats.get('mean', 0)] * min(n, 500))
        
        self.len_mean = np.mean(lengths) if lengths else 15
        self.len_std = np.std(lengths) if lengths else 5
        self.chg_mean = np.mean(charges) if charges else 0
        self.chg_std = np.std(charges) if charges else 3
    
    def score(self, cdr3: str) -> float:
        z_len = abs(len(cdr3) - self.len_mean) / max(self.len_std, 1)
        chg = sum(1 for a in cdr3 if a in 'KRH') - sum(1 for a in cdr3 if a in 'DE')
        z_chg = abs(chg - self.chg_mean) / max(self.chg_std, 0.5)
        return round(max(0, 100 - (z_len + z_chg) / 2 * 25), 1)

# ============================================================
# MAIN DESIGNER
# ============================================================

class VHHDesignerV5:
    """
    Three-pass VHH designer:
    1. Compensation rules (CDR → FR)
    2. Epistasis rules (FR-FR given CDR context)
    3. Vernier archetype (family pattern)
    
    Memory-optimized to handle large epistasis files.
    """
    
    def __init__(self, compensation_file: str, epistasis_file: str):
        print("Loading rule engines...")
        import gc
        
        # Load compensation data (small file, ~1MB)
        with open(compensation_file, 'rb') as f:
            self.comp_data = pickle.load(f)
        
        # Load epistasis data, extract what we need, then discard
        print(f"  Loading epistasis file...")
        with open(epistasis_file, 'rb') as f:
            epi_data = pickle.load(f)
        
        # Extract only what's needed for VernierEngine and NaturalnessScorer
        vernier_clusters = epi_data.get('analysis_2_vernier_clusters', {})
        higher_order_rules = epi_data.get('analysis_4_higher_order_rules', [])
        
        # Create a minimal epi_data dict for VernierEngine and NaturalnessScorer
        self.epi_data_minimal = {
            'analysis_2_vernier_clusters': vernier_clusters,
            'analysis_4_higher_order_rules': higher_order_rules,
        }
        
        # Delete the full data to free memory
        del epi_data
        gc.collect()
        
        # Initialize engines (comp_engine will re-load comp file, but it's small)
        self.comp_engine = CompensationRulesEngine(compensation_file)
        
        # Create epistasis engine from already-loaded data
        self.epi_engine = EpistasisRulesEngineFromData(higher_order_rules, vernier_clusters)
        
        self.vernier_engine = VernierEngine(self.comp_data, self.epi_data_minimal)
        self.scorer = NaturalnessScorer(self.epi_data_minimal)
        
        self.scaffold = UNIVERSAL_SCAFFOLD.copy()
        
        # Final cleanup
        gc.collect()
        print("  Engines loaded (memory optimized)")
    
    def design(self, cdrs: CDRSet, n_candidates: int = 92, mode: str = 'both',
               input_sequence: str = None) -> List[VHHCandidate]:
        """
        Design VHH candidates using three-pass rule application.
        
        For Universal mode: FR2 core is protected, only safe FR3/FR4 rules applied
        For Original mode: Hallmarks → Vernier → Compensation/Epistasis
        """
        # Classify input into family
        family = classify_family(cdrs, input_sequence)
        print(f"\nClassified family: {family}")
        
        # Get CDR3 properties for display
        cdr3_cats = get_cdr3_category(cdrs.cdr3)
        print(f"CDR3: length={cdr3_cats['length']} ({cdr3_cats['cdr3_len']}), "
              f"charge={cdr3_cats['charge']} ({cdr3_cats['cdr3_charge']}), "
              f"cys={cdr3_cats['n_cys']} ({cdr3_cats['cdr3_cys']})")
        
        # === Gather all rules (for display and later filtering) ===
        print("\n--- Rule Collection ---")
        
        # Compensation rules (CDR → FR)
        comp_rules = self.comp_engine.find_applicable_rules(cdrs, family, min_confidence=75.0)
        print(f"  Compensation: {len(comp_rules)} rules for {family}")
        for r in comp_rules[:3]:
            print(f"    {r['fw_position']} -> {r['suggested_aa']} ({r['confidence']:.1f}%) [{r['rule']}]")
        
        # Epistasis rules (FR-FR given CDR context)
        # Note: These are checked against INPUT FRs, but will be re-filtered for each mode
        current_frs = {'FR1': cdrs.fr1, 'FR2': cdrs.fr2, 'FR3': cdrs.fr3, 'FR4': cdrs.fr4}
        epi_rules = self.epi_engine.find_applicable_rules(cdrs, family, current_frs, min_confidence=70.0)
        print(f"  Epistasis: {len(epi_rules)} rules for {family}")
        for r in epi_rules[:3]:
            cond = r.get('condition', '')[:35]
            print(f"    {r['fw_position']} -> {r['suggested_aa']} ({r['confidence']:.1f}%) [{cond}]")
        
        # Vernier archetype (family pattern)
        vernier_rules = self.vernier_engine.get_family_pattern(family)
        print(f"  Vernier: {len(vernier_rules)} positions for {family}")
        for r in vernier_rules[:3]:
            print(f"    {r['fw_position']} -> {r['suggested_aa']}")
        
        # === Merge rules with priority ===
        # Priority: epistasis > compensation > vernier (higher confidence wins ties)
        all_rules = self._merge_rules(comp_rules, epi_rules, vernier_rules)
        print(f"\nMerged: {len(all_rules)} unique position rules")
        
        # === Generate candidates ===
        if mode == 'both':
            n_univ, n_orig = n_candidates // 2, n_candidates - n_candidates // 2
        elif mode == 'universal':
            n_univ, n_orig = n_candidates, 0
        else:
            n_univ, n_orig = 0, n_candidates
        
        candidates = []
        
        # Add original input
        if input_sequence:
            orig_lead = VHHCandidate(
                id="Input_Original", rank=0, sequence=input_sequence,
                framework_source='input', family=family,
                cdr1=cdrs.cdr1, cdr2=cdrs.cdr2, cdr3=cdrs.cdr3,
                fr1=cdrs.fr1, fr2=cdrs.fr2, fr3=cdrs.fr3, fr4=cdrs.fr4,
                mutations=[], strategy="original_input"
            )
            candidates.append(orig_lead)
        
        if n_univ > 0:
            print(f"\nGenerating {n_univ} Universal candidates...")
            candidates.extend(self._design_universal(cdrs, family, all_rules, n_univ))
        
        if n_orig > 0 and cdrs.fr1 and cdrs.fr2:
            print(f"\nGenerating {n_orig} Original candidates...")
            candidates.extend(self._design_original(cdrs, family, all_rules, n_orig))
        elif n_orig > 0:
            print("No original FRs, using Universal")
            candidates.extend(self._design_universal(cdrs, family, all_rules, n_orig))
        
        # Score and rank
        for c in candidates:
            c.naturalness_score = self.scorer.score(c.cdr3)
            c.family = family
            if c.mutations:
                c.confidence = 0.4 * c.naturalness_score + 0.6 * np.mean([m.confidence for m in c.mutations])
            else:
                c.confidence = c.naturalness_score * 0.85
        
        # Sort by confidence, pin leads at top
        candidates.sort(key=lambda x: -x.confidence)
        input_orig = [c for c in candidates if c.id == "Input_Original"]
        univ_graft = [c for c in candidates if c.id == "Univ_Graft"]
        others = [c for c in candidates if c.id not in ["Input_Original", "Univ_Graft"]]
        candidates = input_orig + univ_graft + others
        
        for i, c in enumerate(candidates):
            c.rank = i + 1
        
        return candidates[:n_candidates]
    
    def _merge_rules(self, comp_rules: List[dict], epi_rules: List[dict], 
                     vernier_rules: List[dict]) -> List[dict]:
        """
        Merge rules with conflict resolution.
        Higher confidence wins. For ties: epistasis > compensation > vernier.
        """
        merged = {}
        
        # Add in order of priority (lowest to highest)
        for r in vernier_rules:
            pos = r['fw_position']
            merged[pos] = r
        
        for r in comp_rules:
            pos = r['fw_position']
            if pos not in merged or r['confidence'] >= merged[pos]['confidence']:
                merged[pos] = r
        
        for r in epi_rules:
            pos = r['fw_position']
            if pos not in merged or r['confidence'] >= merged[pos]['confidence']:
                merged[pos] = r
        
        return list(merged.values())
    
    def _design_universal(self, cdrs: CDRSet, family: str, all_rules: List[dict], n: int) -> List[VHHCandidate]:
        """
        Generate candidates using Universal scaffold.
        
        Protected positions (FR2 core motif) are NOT mutated by rules.
        Only hallmark positions (IMGT 42, 49, 50, 52) can be changed via hallmark patterns.
        """
        candidates = []
        used = set()
        
        # Show scaffold hallmarks
        print(f"\n  Scaffold hallmarks (IMGT 42,49,50,52): F-G-L-A (humanized baseline)")
        print(f"  Protected FR2 core: WFRQAPGQGLEAVA (positions 39-52 except hallmarks)")
        
        # Pure graft (no mutations)
        candidates.append(self._make("Univ_Graft", cdrs, family, [], "universal_graft", 'universal'))
        
        # Filter rules to exclude protected positions in universal mode
        safe_rules = []
        for r in all_rules:
            mut = self._rule2mut(r, self.scaffold, mode='universal')
            if mut:  # Only keep rules that produce valid mutations
                safe_rules.append(r)
        
        print(f"  {len(all_rules)} total rules → {len(safe_rules)} safe for universal scaffold")
        
        # Single rules (only safe ones)
        for i, r in enumerate(safe_rules[:min(10, n//4)]):
            m = self._rule2mut(r, self.scaffold, mode='universal')
            if m:
                candidates.append(self._make(f"Univ_{r['source']}_{i+1}", cdrs, family, [m], 
                                              f"single_{r['source']}", 'universal'))
        
        # All safe rules combined
        all_muts = [self._rule2mut(r, self.scaffold, mode='universal') for r in safe_rules]
        all_muts = [m for m in all_muts if m]
        if all_muts:
            candidates.append(self._make("Univ_AllRules", cdrs, family, all_muts, "all_rules", 'universal'))
        
        # Pairs of high-confidence safe rules
        high = [r for r in safe_rules if r['confidence'] >= 85]
        for i in range(min(len(high), 8)):
            for j in range(i+1, min(len(high), 8)):
                if high[i]['fw_position'] != high[j]['fw_position']:
                    m1 = self._rule2mut(high[i], self.scaffold, mode='universal')
                    m2 = self._rule2mut(high[j], self.scaffold, mode='universal')
                    if m1 and m2:
                        key = tuple(sorted([str(m1), str(m2)]))
                        if key not in used:
                            used.add(key)
                            candidates.append(self._make(f"Univ_Pair_{len(candidates)}", cdrs, family,
                                                         [m1, m2], "pair", 'universal'))
                if len(candidates) >= n:
                    break
            if len(candidates) >= n:
                break
        
        # Random diverse combinations (only from safe rules)
        import random
        attempts = 0
        max_attempts = n * 5
        while len(candidates) < n and len(safe_rules) >= 2 and attempts < max_attempts:
            attempts += 1
            picked = random.sample(safe_rules, min(3, len(safe_rules)))
            if len({r['fw_position'] for r in picked}) == len(picked):
                muts = [self._rule2mut(r, self.scaffold, mode='universal') for r in picked]
                muts = [m for m in muts if m]
                if muts:
                    key = tuple(sorted([str(m) for m in muts]))
                    if key not in used:
                        used.add(key)
                        candidates.append(self._make(f"Univ_Div_{len(candidates)}", cdrs, family,
                                                     muts, "diverse", 'universal'))
        
        return candidates[:n]
    
    def _design_original(self, cdrs: CDRSet, family: str, all_rules: List[dict], n: int) -> List[VHHCandidate]:
        """Generate candidates using original FRs + hallmark mutations."""
        import random
        candidates = []
        used = set()
        orig = {'FR1': cdrs.fr1, 'FR2': cdrs.fr2, 'FR3': cdrs.fr3, 'FR4': cdrs.fr4}
        
        # Get all applicable mutations upfront
        all_muts = [self._rule2mut(r, orig) for r in all_rules]
        all_muts = [m for m in all_muts if m]
        
        # Hallmark variants for each hallmark pattern
        for hname, hpos in VHH_HALLMARKS_IMGT.items():
            muts = [self._imgt2mut(imgt_pos, aa, cdrs, f"hallmark_{hname}", 95.0)
                    for imgt_pos, aa in hpos.items()]
            muts = [m for m in muts if m]
            if muts:
                candidates.append(self._make_orig(f"Orig_{hname}", cdrs, family, muts,
                                                   f"hallmark_{hname}", orig))
        
        # If no hallmark mutations work (IMGT lookup failed), use rule-only candidates
        base_muts = [self._imgt2mut(imgt_pos, aa, cdrs, "FERG", 95.0)
                     for imgt_pos, aa in VHH_HALLMARKS_IMGT['FERG'].items()]
        base_muts = [m for m in base_muts if m]
        
        # Fallback: if no hallmark mutations, just use rules directly
        if not base_muts:
            print("  Warning: IMGT hallmark positions not found, using rule-only candidates")
            base_muts = []  # Empty base, rules only
        
        # === Strategy 1: Single rule additions ===
        for i, r in enumerate(all_rules):
            if len(candidates) >= n:
                break
            m = self._rule2mut(r, orig)
            if m:
                cand_muts = base_muts + [m] if base_muts else [m]
                key = tuple(sorted([str(x) for x in cand_muts]))
                if key not in used:
                    used.add(key)
                    prefix = "Orig_FERG_" if base_muts else "Orig_"
                    candidates.append(self._make_orig(f"{prefix}{r['source']}_{i+1}", cdrs, family,
                                                       cand_muts, f"single_{r['source']}", orig))
        
        # === Strategy 2: All rules combined ===
        if all_muts:
            cand_muts = base_muts + all_muts
            key = tuple(sorted([str(x) for x in cand_muts]))
            if key not in used:
                used.add(key)
                candidates.append(self._make_orig("Orig_AllRules", cdrs, family,
                                                   cand_muts, "all_rules", orig))
        
        # === Strategy 3: Pairs of rules ===
        for i in range(len(all_rules)):
            if len(candidates) >= n:
                break
            for j in range(i+1, len(all_rules)):
                if len(candidates) >= n:
                    break
                if all_rules[i]['fw_position'] != all_rules[j]['fw_position']:
                    m1, m2 = self._rule2mut(all_rules[i], orig), self._rule2mut(all_rules[j], orig)
                    if m1 and m2:
                        cand_muts = base_muts + [m1, m2]
                        key = tuple(sorted([str(x) for x in cand_muts]))
                        if key not in used:
                            used.add(key)
                            candidates.append(self._make_orig(f"Orig_Pair_{len(candidates)}", cdrs, family,
                                                               cand_muts, "pair", orig))
        
        # === Strategy 4: Triplets of rules ===
        for i in range(len(all_rules)):
            if len(candidates) >= n:
                break
            for j in range(i+1, len(all_rules)):
                if len(candidates) >= n:
                    break
                for k in range(j+1, len(all_rules)):
                    if len(candidates) >= n:
                        break
                    positions = {all_rules[i]['fw_position'], all_rules[j]['fw_position'], all_rules[k]['fw_position']}
                    if len(positions) == 3:  # All different positions
                        m1 = self._rule2mut(all_rules[i], orig)
                        m2 = self._rule2mut(all_rules[j], orig)
                        m3 = self._rule2mut(all_rules[k], orig)
                        if m1 and m2 and m3:
                            cand_muts = base_muts + [m1, m2, m3]
                            key = tuple(sorted([str(x) for x in cand_muts]))
                            if key not in used:
                                used.add(key)
                                candidates.append(self._make_orig(f"Orig_Triple_{len(candidates)}", cdrs, family,
                                                                   cand_muts, "triple", orig))
        
        # === Strategy 5: Different hallmark bases with rules ===
        for hname in ['YERL', 'FERF', 'FERA', 'FGLA']:
            if len(candidates) >= n:
                break
            if hname not in VHH_HALLMARKS_IMGT:
                continue
            hmuts = [self._imgt2mut(imgt_pos, aa, cdrs, hname, 95.0)
                     for imgt_pos, aa in VHH_HALLMARKS_IMGT[hname].items()]
            hmuts = [m for m in hmuts if m]
            
            # Add each rule individually
            for r in all_rules:
                if len(candidates) >= n:
                    break
                m = self._rule2mut(r, orig)
                if m:
                    cand_muts = hmuts + [m]
                    key = tuple(sorted([str(x) for x in cand_muts]))
                    if key not in used:
                        used.add(key)
                        candidates.append(self._make_orig(f"Orig_{hname}_{r['source']}_{len(candidates)}", cdrs, family,
                                                          cand_muts, f"{hname}+{r['source']}", orig))
        
        # === Strategy 6: Random diverse combinations ===
        attempts = 0
        max_attempts = n * 5
        while len(candidates) < n and attempts < max_attempts:
            attempts += 1
            # Pick random subset of rules (1-4 rules)
            n_pick = random.randint(1, min(4, len(all_rules)))
            if len(all_rules) < n_pick:
                continue
            picked = random.sample(all_rules, n_pick)
            
            # Check all different positions
            positions = {r['fw_position'] for r in picked}
            if len(positions) != len(picked):
                continue
            
            muts = [self._rule2mut(r, orig) for r in picked]
            muts = [m for m in muts if m]
            if not muts:
                continue
            
            # Random hallmark base (or none)
            hname = random.choice(['FERG', 'YERL', 'FERF', 'FERA', None])
            if hname:
                hmuts = [self._imgt2mut(imgt_pos, aa, cdrs, hname, 95.0)
                         for imgt_pos, aa in VHH_HALLMARKS_IMGT.get(hname, {}).items()]
                hmuts = [m for m in hmuts if m]
            else:
                hmuts = []
            
            cand_muts = hmuts + muts
            key = tuple(sorted([str(x) for x in cand_muts]))
            if key not in used:
                used.add(key)
                prefix = f"Orig_{hname}_" if hname and hmuts else "Orig_"
                candidates.append(self._make_orig(f"{prefix}Div_{len(candidates)}", cdrs, family,
                                                   cand_muts, "diverse", orig))
        
        return candidates[:n]
    
    def _rule2mut(self, rule: dict, frs: dict, mode: str = 'original') -> Optional[Mutation]:
        """
        Convert rule to mutation given current FRs.
        
        For universal mode, protects FR2 core motif except hallmark positions.
        """
        pos = rule['fw_position']
        
        # Handle FR[idx] format
        match = re.match(r'(FR\d)\[(\d+)\]', pos)
        if match:
            region, idx = match.groups()
            idx = int(idx)
            fr_seq = frs.get(region, '')
            if idx >= len(fr_seq):
                return None
            
            # Check protected positions for universal mode
            if mode == 'universal' and region == 'FR2':
                # Convert 0-based index to IMGT position
                imgt_pos = SUBSTR_TO_IMGT.get('FR2', {}).get(idx + 1)  # +1 because SUBSTR_TO_IMGT is 1-based
                if imgt_pos and imgt_pos in PROTECTED_IMGT_POSITIONS:
                    return None  # Skip protected positions
            
            cur = fr_seq[idx]
            if cur == rule['suggested_aa']:
                return None
            return Mutation(pos, cur, rule['suggested_aa'], rule.get('rule', ''),
                           rule.get('confidence', 75), rule.get('source', 'unknown'),
                           rule.get('family', ''), rule.get('condition', ''))
        
        # Handle FR_idx format (from epistasis, 1-based)
        match = re.match(r'(FR\d)_(\d+)', pos)
        if match:
            region, idx = match.groups()
            idx = int(idx)  # 1-based index from epistasis
            
            # Check protected positions for universal mode
            if mode == 'universal' and region == 'FR2':
                imgt_pos = SUBSTR_TO_IMGT.get('FR2', {}).get(idx)
                if imgt_pos and imgt_pos in PROTECTED_IMGT_POSITIONS:
                    return None  # Skip protected positions
            
            # Convert to 0-based for sequence access
            seq_idx = idx - 1 if idx > 0 else idx
            fr_seq = frs.get(region, '')
            if seq_idx >= len(fr_seq) or seq_idx < 0:
                return None
            cur = fr_seq[seq_idx]
            if cur == rule['suggested_aa']:
                return None
            pos_str = f'{region}[{seq_idx}]'  # Store as 0-based
            return Mutation(pos_str, cur, rule['suggested_aa'], rule.get('rule', ''),
                           rule.get('confidence', 75), rule.get('source', 'unknown'),
                           rule.get('family', ''), rule.get('condition', ''))
        
        return None
    
    def _imgt2mut(self, imgt_pos: int, to_aa: str, cdrs: CDRSet, 
                  rule: str, conf: float) -> Optional[Mutation]:
        """Create mutation from IMGT position."""
        if not cdrs.imgt_numbering or imgt_pos not in cdrs.imgt_numbering:
            return None
        
        info = cdrs.imgt_numbering[imgt_pos]
        cur_aa = info['aa']
        
        if cur_aa == to_aa:
            return None
        
        pos_str = f"IMGT{imgt_pos}"
        return Mutation(pos_str, cur_aa, to_aa, rule, conf, 'hallmark')
    
    def _make(self, id: str, cdrs: CDRSet, family: str, mutations: List[Mutation],
              strategy: str, fw_source: str) -> VHHCandidate:
        """Build candidate using Universal scaffold."""
        frs = {
            'FR1': list(self.scaffold['FR1']),
            'FR2': list(self.scaffold['FR2']),
            'FR3': list(self.scaffold['FR3']),
            'FR4': list(self.scaffold['FR4']),
        }
        
        self._apply_mutations(frs, mutations, cdrs)
        
        fr1 = ''.join(frs['FR1'])
        fr2 = ''.join(frs['FR2'])
        fr3 = ''.join(frs['FR3'])
        fr4 = ''.join(frs['FR4'])
        seq = fr1 + cdrs.cdr1 + fr2 + cdrs.cdr2 + fr3 + cdrs.cdr3 + fr4
        
        return VHHCandidate(id, 0, seq, fw_source, cdrs.cdr1, cdrs.cdr2, cdrs.cdr3,
                           fr1, fr2, fr3, fr4, mutations, family=family, strategy=strategy)
    
    def _make_orig(self, id: str, cdrs: CDRSet, family: str, mutations: List[Mutation],
                   strategy: str, orig_frs: dict) -> VHHCandidate:
        """Build candidate using original FRs."""
        frs = {
            'FR1': list(orig_frs['FR1']),
            'FR2': list(orig_frs['FR2']),
            'FR3': list(orig_frs['FR3']),
            'FR4': list(orig_frs['FR4']),
        }
        
        self._apply_mutations(frs, mutations, cdrs)
        
        fr1 = ''.join(frs['FR1'])
        fr2 = ''.join(frs['FR2'])
        fr3 = ''.join(frs['FR3'])
        fr4 = ''.join(frs['FR4'])
        seq = fr1 + cdrs.cdr1 + fr2 + cdrs.cdr2 + fr3 + cdrs.cdr3 + fr4
        
        return VHHCandidate(id, 0, seq, 'original', cdrs.cdr1, cdrs.cdr2, cdrs.cdr3,
                           fr1, fr2, fr3, fr4, mutations, family=family, strategy=strategy)
    
    def _apply_mutations(self, frs: dict, mutations: List[Mutation], cdrs: CDRSet):
        """Apply mutations to FRs dict (mutates in place)."""
        for m in mutations:
            # FR[idx] format
            match = re.match(r'(FR\d)\[(\d+)\]', m.position)
            if match:
                region, idx = match.groups()
                idx = int(idx)
                if region in frs and idx < len(frs[region]):
                    frs[region][idx] = m.to_aa
                continue
            
            # IMGT format
            match = re.match(r'IMGT(\d+)', m.position)
            if match and cdrs.imgt_numbering:
                imgt_pos = int(match.group(1))
                if imgt_pos in cdrs.imgt_numbering:
                    info = cdrs.imgt_numbering[imgt_pos]
                    region = info['region']
                    idx = info['index']
                    if region in frs and idx < len(frs[region]):
                        frs[region][idx] = m.to_aa

# ============================================================
# OUTPUT FUNCTIONS
# ============================================================

def format_mutations(mutations: List[Mutation]) -> str:
    if not mutations:
        return "None"
    return "; ".join(f"{m.position}:{m.from_aa}->{m.to_aa} [{m.source}, {m.confidence:.0f}%]" 
                     for m in mutations)

def to_dataframe(candidates: List[VHHCandidate], input_seq: str, cdrs: CDRSet) -> pd.DataFrame:
    return pd.DataFrame([{
        'rank': c.rank,
        'id': c.id,
        'sequence': c.sequence,
        'length': len(c.sequence),
        'family': c.family,
        'framework_source': c.framework_source,
        'CDR1': c.cdr1,
        'CDR2': c.cdr2,
        'CDR3': c.cdr3,
        'CDR3_preserved': c.cdr3 == cdrs.cdr3,
        'FR1': c.fr1,
        'FR2': c.fr2,
        'FR3': c.fr3,
        'FR4': c.fr4,
        'confidence': round(c.confidence, 1),
        'naturalness_score': round(c.naturalness_score, 1),
        'n_mutations': len(c.mutations),
        'mutations': format_mutations(c.mutations),
        'rules_applied': "; ".join(m.rule for m in c.mutations) if c.mutations else "None",
        'strategy': c.strategy,
        'input_sequence': input_seq,
    } for c in candidates])

def save_results(df: pd.DataFrame, output_dir: str, name: str, mode: str, datetime_str: str):
    """Save results to folder with xlsx and MSA-ready csv."""
    os.makedirs(output_dir, exist_ok=True)
    
    base_name = f"{datetime_str}_{name}_{mode}"
    xlsx_path = os.path.join(output_dir, f"{base_name}.xlsx")
    csv_path = os.path.join(output_dir, f"{base_name}_sequences.csv")
    
    # Save Excel
    with pd.ExcelWriter(xlsx_path, engine='openpyxl') as writer:
        df.to_excel(writer, sheet_name='Candidates', index=False)
        
        # Summary sheet
        summary = pd.DataFrame({
            'Metric': ['Total candidates', 'Universal', 'Original', 'Input', 
                       'Mean confidence', 'Mean naturalness', 'Family'],
            'Value': [len(df), (df['framework_source'] == 'universal').sum(),
                     (df['framework_source'] == 'original').sum(),
                     (df['framework_source'] == 'input').sum(),
                     f"{df['confidence'].mean():.1f}",
                     f"{df['naturalness_score'].mean():.1f}",
                     df['family'].iloc[0] if len(df) > 0 else 'N/A']
        })
        summary.to_excel(writer, sheet_name='Summary', index=False)
    
    print(f"  Saved: {xlsx_path}")
    
    # Save MSA-ready CSV
    msa_df = df[['id', 'sequence']].copy()
    msa_df.to_csv(csv_path, index=False)
    print(f"  Saved: {csv_path}")

# ============================================================
# MAIN
# ============================================================

def main():
    p = argparse.ArgumentParser(description='VHH Designer v5 - Protected Framework Rule Application')
    p.add_argument('--sequence', '-s', help='Input sequence')
    p.add_argument('--input', '-i', help='Input FASTA file')
    p.add_argument('--compensation', '-c', required=True, 
                   help='Compensation rules file (correlation_results_v3_compensation.pkl)')
    p.add_argument('--epistasis', '-e', required=True,
                   help='Epistasis file (epistasis_v2_full.pkl)')
    p.add_argument('--output-dir', '-o', help='Output directory')
    p.add_argument('--name', help='Sequence name (auto-detected from FASTA)')
    p.add_argument('--n-candidates', '-n', type=int, default=92)
    p.add_argument('--mode', '-m', choices=['both', 'universal', 'original'], default='both')
    args = p.parse_args()
    
    # Get input sequence
    seq_name = args.name
    if args.sequence:
        input_seq = args.sequence.upper().replace(' ', '').replace('\n', '')
        if not seq_name:
            seq_name = "input"
    elif args.input:
        with open(args.input) as f:
            content = f.read()
        if args.input.endswith(('.fasta', '.fa')):
            lines = content.strip().split('\n')
            header = [l for l in lines if l.startswith('>')]
            if header and not seq_name:
                seq_name = header[0][1:].split()[0].replace('/', '_').replace('\\', '_')
            input_seq = ''.join(l for l in lines if not l.startswith('>')).upper()
        else:
            input_seq = content.strip().split('\n')[0].upper()
        if not seq_name:
            seq_name = args.input.rsplit('/', 1)[-1].rsplit('.', 1)[0]
    else:
        p.error("Need --sequence or --input")
        return
    
    seq_name = re.sub(r'[^\w\-]', '_', seq_name)[:50]
    
    print("=" * 70)
    print("VHH DESIGNER v5 - Protected Framework Rule Application")
    print("=" * 70)
    print(f"Sequence: {seq_name}")
    print(f"Mode: {args.mode.upper()}")
    print(f"Input ({len(input_seq)} aa): {input_seq[:50]}...")
    
    # Extract CDRs
    cdrs = extract_cdrs(input_seq)
    if not cdrs:
        print("ERROR: CDR extraction failed")
        return
    
    print(f"\nCDR1: {cdrs.cdr1}")
    print(f"CDR2: {cdrs.cdr2}")
    print(f"CDR3: {cdrs.cdr3}")
    print(f"FR1:  {cdrs.fr1[:20]}... ({len(cdrs.fr1)} aa)" if cdrs.fr1 else "FR1:  (empty)")
    print(f"FR2:  {cdrs.fr2} ({len(cdrs.fr2)} aa)" if cdrs.fr2 else "FR2:  (empty)")
    print(f"FR3:  {cdrs.fr3[:20]}... ({len(cdrs.fr3)} aa)" if cdrs.fr3 else "FR3:  (empty)")
    print(f"FR4:  {cdrs.fr4} ({len(cdrs.fr4)} aa)" if cdrs.fr4 else "FR4:  (empty)")
    
    # Show hallmarks with diagnostic info
    hallmark_aa = []
    imgt_lookup_worked = False
    for imgt_pos in [42, 49, 50, 52]:
        aa = cdrs.get_aa_at_imgt(imgt_pos)
        if aa:
            imgt_lookup_worked = True
        hallmark_aa.append(aa if aa else '?')
    
    if imgt_lookup_worked:
        print(f"Input hallmarks (IMGT 42,49,50,52): {'-'.join(hallmark_aa)}")
    else:
        # Show FR2 fallback positions
        print(f"Input hallmarks (IMGT lookup failed, using FR2 fallback):")
        if len(cdrs.fr2) >= 14:
            fr2_42 = cdrs.fr2[3] if len(cdrs.fr2) > 3 else '?'
            fr2_49 = cdrs.fr2[10] if len(cdrs.fr2) > 10 else '?'
            fr2_50 = cdrs.fr2[11] if len(cdrs.fr2) > 11 else '?'
            fr2_52 = cdrs.fr2[13] if len(cdrs.fr2) > 13 else '?'
            print(f"  FR2 positions: {fr2_42}-{fr2_49}-{fr2_50}-{fr2_52} (pos 3,10,11,13)")
        # Debug: show IMGT keys we have
        if cdrs.imgt_numbering:
            keys_sample = sorted([k for k in cdrs.imgt_numbering.keys() if isinstance(k, int)])[:10]
            print(f"  IMGT keys present: {keys_sample}... (total: {len(cdrs.imgt_numbering)})")
    
    # Show scaffold hallmarks for reference (only for universal/both modes)
    if args.mode in ['universal', 'both']:
        print(f"Scaffold hallmarks (Universal): F-G-L-A (humanized baseline)")
        print(f"  Note: FR2 core (WFRQAPGQGLEAVA) is PROTECTED - only FR3/FR4 rules applied")
    
    # Design
    designer = VHHDesignerV5(args.compensation, args.epistasis)
    candidates = designer.design(cdrs, args.n_candidates, args.mode, input_sequence=input_seq)
    df = to_dataframe(candidates, input_seq, cdrs)
    
    print(f"\n{'='*70}")
    print(f"RESULTS")
    print(f"{'='*70}")
    print(f"Total: {len(df)} | Family: {df['family'].iloc[0] if len(df) > 0 else 'N/A'}")
    print(f"  Input: {(df['framework_source']=='input').sum()}")
    print(f"  Universal: {(df['framework_source']=='universal').sum()}")
    print(f"  Original: {(df['framework_source']=='original').sum()}")
    
    print("\nTop 10:")
    for _, r in df.head(10).iterrows():
        print(f"  {r['rank']:2d}. {r['id'][:40]:40s} [{r['framework_source'][:4]}] conf={r['confidence']:.1f}")
    
    # Save
    datetime_str = datetime.now().strftime('%Y%m%d_%H%M%S')
    script_name = "vhh_designer_v5"
    folder_name = f"{datetime_str}_{script_name}_{seq_name}_{args.mode}"
    
    if args.output_dir:
        output_dir = args.output_dir
    else:
        results_base = "results/analysis_runs"
        if os.path.exists(results_base):
            output_dir = os.path.join(results_base, folder_name)
        else:
            output_dir = folder_name
    
    save_results(df, output_dir, f"{script_name}_{seq_name}", args.mode, datetime_str)
    print("\nDone!")

if __name__ == '__main__':
    main()