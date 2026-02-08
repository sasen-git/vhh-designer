#!/usr/bin/env python3
"""
VHH Designer v7.2 - Compatible with v7.7 Analysis Outputs
==========================================================

Updated to work directly with the JSON outputs from vhh_analysis_unified_v7.7.py:
  - analysis_rules_v7.json (all rules: compensation + triplet)
  - analysis_vernier_archetypes_v7.json (family patterns)
  - analysis_vernier_clusters_v7.json (for naturalness scoring)

For UNIVERSAL mode (grafting CDRs onto humanized scaffold):
  1. Start with Universal scaffold (FR2 = WFRQAPGQGLEAVA, hallmarks = FGLA)
  2. Apply rules from JSON gated by hallmarks/CDR features
  3. FR2 core is PROTECTED (only hallmark positions can change)
    
For ORIGINAL mode (keeping input FRs with VHH-izing mutations):
  PASS 1: HALLMARK mutations (IMGT 42,49,50,52 → target pattern like FERG, YERL)
  PASS 2: Apply rules based on conditions (hallmarks, cdr3_len, cdr3_charge, etc.)

Family Classification (based on IMGT positions):
  - F_C2: pos42=F, 2 cysteines
  - F_C4: pos42=F, 4+ cysteines
  - Y_C2: pos42=Y, 2 cysteines  
  - Y_C4: pos42=Y, 4+ cysteines
  - VH_like: pos50=L (human VH-like)
  - VHH_W52: pos52=W
  - Other_VHH: other combinations

Usage:
  python vhh_designer_v7_2.py -i mouse_hc.fasta \\
      --rules models/analysis/v7.7_*/analysis_rules_v7.json \\
      --archetypes models/analysis/v7.7_*/analysis_vernier_archetypes_v7.json \\
      --target-hallmarks FERG \\
      --mode both

Author: Claude (Anthropic)
Date: January 2026
"""

import os
import re
import json
import argparse
import numpy as np
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

# IMGT position mappings for each FR region
# FR1: positions 1-26, FR2: 39-55 (hallmarks at 42,49,50,52), FR3: 66-104, FR4: 118-128
IMGT_REGIONS = {
    'FR1': (1, 26),
    'FR2': (39, 55),
    'FR3': (66, 104),
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
    'VGLW': {42: 'V', 49: 'G', 50: 'L', 52: 'W'},  # VH-like
    'FGLA': {42: 'F', 49: 'G', 50: 'L', 52: 'A'},  # Humanized
}

# FR2 index to IMGT mapping (0-based index)
# FR2: W-F-R-Q-A-P-G-Q-G-L-E-A-V-A
# IMGT: 39-40-41-42-43-44-45-46-47-48-49-50-51-52-53-54 (but 41 doesn't exist, so offset)
# Actually: IMGT 39-55 maps to FR2 positions
FR2_IDX_TO_IMGT = {
    0: 39, 1: 40, 2: 41, 3: 42, 4: 43, 5: 44, 6: 45, 7: 46,
    8: 47, 9: 48, 10: 49, 11: 50, 12: 51, 13: 52, 14: 53, 15: 54, 16: 55
}

# ============================================================
# DATA CLASSES
# ============================================================

@dataclass
class CDRSet:
    cdr1: str
    cdr2: str
    cdr3: str
    fr1: str = ""
    fr2: str = ""
    fr3: str = ""
    fr4: str = ""

@dataclass
class Mutation:
    position: str      # e.g., "IMGT66" or "FR3_1"
    imgt_num: int      # IMGT position number
    original: str      # Original AA
    mutant: str        # New AA
    source: str        # Rule source
    confidence: float  # Confidence score
    
    def __str__(self):
        return f"{self.position}:{self.original}→{self.mutant}"

@dataclass
class VHHCandidate:
    id: str
    rank: int
    sequence: str
    framework_source: str  # 'universal', 'original', 'input'
    family: str
    cdr1: str
    cdr2: str
    cdr3: str
    fr1: str
    fr2: str
    fr3: str
    fr4: str
    mutations: List[Mutation] = field(default_factory=list)
    strategy: str = ""
    is_lead: bool = False
    confidence: float = 0.0
    naturalness_score: float = 0.0

# ============================================================
# CDR EXTRACTION
# ============================================================

def extract_cdrs(sequence: str) -> Optional[CDRSet]:
    """Extract CDRs and FRs using AntPack."""
    if not ANTPACK_AVAILABLE:
        print("ERROR: AntPack required for CDR extraction")
        return None
    
    try:
        annotator = SingleChainAnnotator()
        numbering, chain_type, e_val = annotator.analyze_seq(sequence)
        
        if chain_type not in ['H', 'K', 'L']:
            print(f"Warning: Chain type {chain_type} detected")
        
        # Extract regions by IMGT numbering
        regions = {'FR1': [], 'CDR1': [], 'FR2': [], 'CDR2': [], 
                   'FR3': [], 'CDR3': [], 'FR4': []}
        
        imgt_ranges = {
            'FR1': (1, 26), 'CDR1': (27, 38), 'FR2': (39, 55), 'CDR2': (56, 65),
            'FR3': (66, 104), 'CDR3': (105, 117), 'FR4': (118, 128)
        }
        
        for i, (pos, aa) in enumerate(zip(numbering, sequence)):
            if pos == '-' or aa == '-':
                continue
            try:
                pos_num = int(pos.rstrip('ABCDEFGHIJ'))
                for region, (start, end) in imgt_ranges.items():
                    if start <= pos_num <= end:
                        regions[region].append(aa)
                        break
            except ValueError:
                continue
        
        return CDRSet(
            cdr1=''.join(regions['CDR1']),
            cdr2=''.join(regions['CDR2']),
            cdr3=''.join(regions['CDR3']),
            fr1=''.join(regions['FR1']),
            fr2=''.join(regions['FR2']),
            fr3=''.join(regions['FR3']),
            fr4=''.join(regions['FR4']),
        )
    except Exception as e:
        print(f"CDR extraction error: {e}")
        return None

# ============================================================
# FAMILY CLASSIFICATION
# ============================================================

def classify_family(cdrs: CDRSet, full_sequence: str = None) -> str:
    """Classify VHH family based on hallmarks and cysteine count."""
    # Get hallmark positions from FR2
    fr2 = cdrs.fr2 if cdrs.fr2 else ""
    
    # FR2 IMGT mapping: index 3=IMGT42, index 10=IMGT49, index 11=IMGT50, index 13=IMGT52
    p42 = fr2[3] if len(fr2) > 3 else "-"
    p49 = fr2[10] if len(fr2) > 10 else "-"
    p50 = fr2[11] if len(fr2) > 11 else "-"
    p52 = fr2[13] if len(fr2) > 13 else "-"
    
    # VH-like check
    if p50 == "L":
        return "VH_like"
    
    # W52 check
    if p52 == "W":
        return "VHH_W52"
    
    # Count cysteines
    seq = full_sequence or (cdrs.fr1 + cdrs.cdr1 + cdrs.fr2 + cdrs.cdr2 + cdrs.fr3 + cdrs.cdr3 + cdrs.fr4)
    cys_count = seq.count("C")
    c_label = "C4" if cys_count >= 4 else "C2"
    
    # F or Y family
    if p42 == "F":
        return f"F_{c_label}"
    elif p42 == "Y":
        return f"Y_{c_label}"
    
    return "Other_VHH"

def get_hallmarks_from_fr2(fr2: str) -> str:
    """Extract hallmark pattern from FR2."""
    if not fr2 or len(fr2) < 14:
        return "????"
    # FR2 indices: 3=IMGT42, 10=IMGT49, 11=IMGT50, 13=IMGT52
    p42 = fr2[3] if len(fr2) > 3 else "?"
    p49 = fr2[10] if len(fr2) > 10 else "?"
    p50 = fr2[11] if len(fr2) > 11 else "?"
    p52 = fr2[13] if len(fr2) > 13 else "?"
    return f"{p42}{p49}{p50}{p52}"

# ============================================================
# CDR3 FEATURES
# ============================================================

def get_cdr3_features(cdr3: str) -> Dict[str, str]:
    """Extract CDR3 features for rule matching."""
    length = len(cdr3)
    
    # Length category
    if length <= 8:
        cdr3_len = "short"
    elif length <= 14:
        cdr3_len = "medium"
    else:
        cdr3_len = "long"
    
    # Charge
    pos_charge = sum(1 for aa in cdr3 if aa in "RKH")
    neg_charge = sum(1 for aa in cdr3 if aa in "DE")
    net_charge = pos_charge - neg_charge
    
    if net_charge >= 2:
        cdr3_charge = "positive"
    elif net_charge <= -2:
        cdr3_charge = "negative"
    else:
        cdr3_charge = "neutral"
    
    # Cysteine
    n_cys = cdr3.count("C")
    if n_cys == 0:
        cdr3_cys = "none"
    elif n_cys == 1:
        cdr3_cys = "one"
    else:
        cdr3_cys = "multiple"
    
    # Terminal residues
    cdr3_end = cdr3[-3:] if len(cdr3) >= 3 else cdr3
    cdr3_last = cdr3[-1] if cdr3 else ""
    
    return {
        "cdr3_len": cdr3_len,
        "cdr3_charge": cdr3_charge,
        "cdr3_cys": cdr3_cys,
        "cdr3_end": cdr3_end,
        "cdr3[-1]": cdr3_last,
        "cdr3[-2]": cdr3[-2] if len(cdr3) >= 2 else "",
        "cdr3[-3]": cdr3[-3] if len(cdr3) >= 3 else "",
    }

# ============================================================
# RULES ENGINE (v7.7 format)
# ============================================================

class RulesEngineV7:
    """
    Rules engine for v7.7 JSON format.
    
    Handles rules with conditions like:
    - hallmarks=FERG
    - cdr3_len=short
    - cdr3_charge=neutral
    - cdr3[-1]=Y
    - family=F_C2 AND cdr3_end=EKF (triplet rules)
    """
    
    def __init__(self, rules_file: str):
        print(f"  Loading rules from {rules_file}...")
        with open(rules_file, 'r') as f:
            self.rules = json.load(f)
        
        # Separate by type
        self.compensation_rules = [r for r in self.rules if r.get('rule_type') == 'cdr_to_imgt']
        self.triplet_rules = [r for r in self.rules if r.get('rule_type') == 'triplet']
        
        # Index compensation rules by condition type
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
        
        print(f"    {len(self.compensation_rules)} compensation rules")
        print(f"    {len(self.triplet_rules)} triplet rules")
        print(f"    {len(self.rules_by_hallmarks)} hallmark patterns")
    
    def find_applicable_rules(self, hallmarks: str, family: str, cdr3_features: Dict[str, str],
                               min_support: int = 1000, min_confidence: float = 0.7) -> List[dict]:
        """
        Find all applicable rules for the given context.
        
        Returns rules sorted by confidence (highest first).
        """
        applicable = []
        
        # 1. Hallmark-based rules
        for rule in self.rules_by_hallmarks.get((hallmarks, family), []):
            if rule.get('support', 0) >= min_support and rule.get('confidence', 0) >= min_confidence:
                applicable.append(rule)
        
        # 2. CDR3 feature-based rules
        for (cond, rule_family), rules in self.rules_by_condition.items():
            if rule_family != family:
                continue
            
            # Parse condition and check if it matches
            if self._check_condition(cond, cdr3_features):
                for rule in rules:
                    if rule.get('support', 0) >= min_support and rule.get('confidence', 0) >= min_confidence:
                        applicable.append(rule)
        
        # Sort by confidence (highest first), then by support
        applicable.sort(key=lambda r: (-r.get('confidence', 0), -r.get('support', 0)))
        
        # Deduplicate: keep best rule per position
        by_position = {}
        for rule in applicable:
            pos = rule.get('position_num', 0)
            if pos not in by_position:
                by_position[pos] = rule
        
        return list(by_position.values())
    
    def _check_condition(self, condition: str, cdr3_features: Dict[str, str]) -> bool:
        """Check if a condition matches the CDR3 features."""
        if '=' not in condition:
            return False
        
        key, value = condition.split('=', 1)
        key = key.strip()
        value = value.strip()
        
        return cdr3_features.get(key) == value
    
    def find_triplet_rules(self, family: str, cdr3: str, min_support: int = 100) -> List[dict]:
        """Find triplet rules for CDR3→FR4 junction."""
        cdr3_end = cdr3[-3:] if len(cdr3) >= 3 else cdr3
        applicable = []
        
        for rule in self.triplet_rules:
            if rule.get('family') != family:
                continue
            if rule.get('support', 0) < min_support:
                continue
            
            # Check if cdr3_end matches
            cond = rule.get('condition', '')
            if f"cdr3_end={cdr3_end}" in cond:
                applicable.append(rule)
        
        return sorted(applicable, key=lambda r: -r.get('confidence', 0))

# ============================================================
# VERNIER/ARCHETYPE ENGINE
# ============================================================

class ArchetypeEngine:
    """Engine for family-specific vernier patterns."""
    
    def __init__(self, archetypes_file: str):
        print(f"  Loading archetypes from {archetypes_file}...")
        with open(archetypes_file, 'r') as f:
            self.archetypes = json.load(f)
        
        # Index by family
        self.by_family = {a['family']: a for a in self.archetypes}
        print(f"    {len(self.archetypes)} family archetypes")
    
    def get_family_pattern(self, family: str) -> Dict[int, str]:
        """Get consensus pattern for family as IMGT position -> AA."""
        arch = self.by_family.get(family, {})
        positions = arch.get('positions', {})
        
        pattern = {}
        for pos_str, data in positions.items():
            try:
                imgt_num = int(pos_str.replace('IMGT', ''))
                pattern[imgt_num] = data.get('consensus', '')
            except ValueError:
                continue
        
        return pattern

# ============================================================
# NATURALNESS SCORER
# ============================================================

class NaturalnessScorer:
    """Simple CDR3 naturalness scorer."""
    
    def __init__(self):
        # Preferred terminal residues for CDR3
        self.preferred_cterm = set("YWFDV")
        self.preferred_nterm = set("CARQ")
    
    def score(self, cdr3: str) -> float:
        if not cdr3:
            return 0.0
        
        score = 50.0  # Base score
        
        # C-terminal preferences
        if cdr3[-1] in self.preferred_cterm:
            score += 20.0
        
        # N-terminal preferences
        if cdr3[0] in self.preferred_nterm:
            score += 10.0
        
        # Length penalty for very short/long
        if len(cdr3) < 6:
            score -= 10.0
        elif len(cdr3) > 20:
            score -= 5.0
        
        # Cysteine check
        if cdr3.count('C') == 0:
            score += 5.0
        elif cdr3.count('C') == 1:
            score -= 5.0  # Odd number of cysteines is usually bad
        
        return min(100.0, max(0.0, score))

# ============================================================
# VHH DESIGNER
# ============================================================

class VHHDesignerV7_2:
    """
    VHH Designer using v7.7 analysis outputs.
    """
    
    def __init__(self, rules_file: str, archetypes_file: str):
        print("\nInitializing VHH Designer v7.2...")
        self.rules_engine = RulesEngineV7(rules_file)
        self.archetype_engine = ArchetypeEngine(archetypes_file)
        self.scorer = NaturalnessScorer()
        self.scaffold = UNIVERSAL_SCAFFOLD
    
    def design(self, cdrs: CDRSet, n_candidates: int = 92, mode: str = 'both',
               input_sequence: str = None, target_hallmarks: str = None) -> List[VHHCandidate]:
        """
        Design VHH candidates.
        
        Args:
            cdrs: Extracted CDR/FR regions
            n_candidates: Number of candidates per mode
            mode: 'universal', 'original', or 'both'
            input_sequence: Original input sequence
            target_hallmarks: Target hallmark pattern (e.g., 'FERG')
        """
        # Classify family
        family = classify_family(cdrs, input_sequence)
        print(f"\nClassified family: {family}")
        
        # Get input hallmarks
        input_hallmarks = get_hallmarks_from_fr2(cdrs.fr2)
        print(f"Input hallmarks (IMGT 42,49,50,52): {input_hallmarks}")
        
        # Determine target hallmarks
        if not target_hallmarks:
            # Auto-select based on family
            if family.startswith('F_'):
                target_hallmarks = 'FERG'
            elif family.startswith('Y_'):
                target_hallmarks = 'YERL'
            else:
                target_hallmarks = 'FERG'
        print(f"Target hallmarks: {target_hallmarks}")
        
        # Get CDR3 features
        cdr3_features = get_cdr3_features(cdrs.cdr3)
        print(f"CDR3: length={cdr3_features['cdr3_len']}, charge={cdr3_features['cdr3_charge']}, "
              f"cys={cdr3_features['cdr3_cys']}, end={cdr3_features['cdr3_end']}")
        
        # Get applicable rules
        print("\n--- Rule Collection ---")
        rules = self.rules_engine.find_applicable_rules(
            target_hallmarks, family, cdr3_features,
            min_support=1000, min_confidence=0.7
        )
        print(f"  Applicable compensation rules: {len(rules)}")
        for r in rules[:5]:
            print(f"    IMGT{r['position_num']} → {r['suggested_aa']} | {r['condition'][:30]} | "
                  f"conf={r['confidence']:.2f} lift={r['lift']:.2f}")
        
        # Get triplet rules
        triplet_rules = self.rules_engine.find_triplet_rules(family, cdrs.cdr3, min_support=100)
        print(f"  Triplet rules: {len(triplet_rules)}")
        if triplet_rules:
            for r in triplet_rules[:3]:
                print(f"    FR4 → {r['suggested_aa']} | {r['condition'][:40]}")
        
        # Get archetype pattern
        archetype = self.archetype_engine.get_family_pattern(family)
        print(f"  Archetype positions: {len(archetype)}")
        
        # Generate candidates
        if mode == 'both':
            n_univ = n_candidates // 2
            n_orig = n_candidates - n_univ
        elif mode == 'universal':
            n_univ, n_orig = n_candidates, 0
        else:
            n_univ, n_orig = 0, n_candidates
        
        candidates = []
        
        # Add original input as lead
        if input_sequence:
            orig_lead = VHHCandidate(
                id="Input_Original", rank=0, sequence=input_sequence,
                framework_source='input', family=family,
                cdr1=cdrs.cdr1, cdr2=cdrs.cdr2, cdr3=cdrs.cdr3,
                fr1=cdrs.fr1, fr2=cdrs.fr2, fr3=cdrs.fr3, fr4=cdrs.fr4,
                mutations=[], strategy="original_input", is_lead=True
            )
            candidates.append(orig_lead)
        
        if n_univ > 0:
            print(f"\nGenerating {n_univ} Universal candidates...")
            candidates.extend(self._design_universal(cdrs, family, rules, archetype, 
                                                      triplet_rules, n_univ, target_hallmarks))
        
        if n_orig > 0 and cdrs.fr1 and cdrs.fr2:
            print(f"\nGenerating {n_orig} Original candidates...")
            candidates.extend(self._design_original(cdrs, family, rules, archetype,
                                                     triplet_rules, n_orig, target_hallmarks))
        elif n_orig > 0:
            print("No original FRs available, using Universal for all")
            candidates.extend(self._design_universal(cdrs, family, rules, archetype,
                                                      triplet_rules, n_orig, target_hallmarks))
        
        # Score and rank
        for c in candidates:
            c.naturalness_score = self.scorer.score(c.cdr3)
            if c.mutations:
                c.confidence = 0.4 * c.naturalness_score + 0.6 * np.mean([m.confidence * 100 for m in c.mutations])
            else:
                c.confidence = c.naturalness_score
        
        # Sort by confidence (leads first)
        candidates.sort(key=lambda c: (-int(c.is_lead), -c.confidence))
        for i, c in enumerate(candidates):
            c.rank = i
        
        return candidates
    
    def _design_universal(self, cdrs: CDRSet, family: str, rules: List[dict],
                           archetype: Dict[int, str], triplet_rules: List[dict],
                           n: int, target_hallmarks: str) -> List[VHHCandidate]:
        """Generate candidates using Universal scaffold."""
        candidates = []
        used = set()
        
        # Pure graft (LEAD for universal)
        graft = self._make_candidate("Univ_Graft", cdrs, family, [], "universal_graft", 'universal')
        graft.is_lead = True
        candidates.append(graft)
        
        # Get mutations from rules
        all_muts = []
        for rule in rules:
            mut = self._rule_to_mutation(rule, self.scaffold, 'universal')
            if mut:
                all_muts.append(mut)
        
        print(f"  {len(all_muts)} applicable mutations for universal scaffold")
        
        # Single mutations
        for i, mut in enumerate(all_muts[:min(15, n//3)]):
            key = str(mut)
            if key not in used:
                used.add(key)
                candidates.append(self._make_candidate(
                    f"Univ_Single_{i+1}", cdrs, family, [mut], f"single_{mut.source}", 'universal'
                ))
        
        # All mutations combined
        if all_muts:
            candidates.append(self._make_candidate(
                "Univ_AllRules", cdrs, family, all_muts, "all_rules", 'universal'
            ))
        
        # Pairs of high-confidence mutations
        high_conf_muts = [m for m in all_muts if m.confidence >= 0.8]
        for i in range(min(len(high_conf_muts), 8)):
            for j in range(i+1, min(len(high_conf_muts), 8)):
                if high_conf_muts[i].imgt_num != high_conf_muts[j].imgt_num:
                    key = tuple(sorted([str(high_conf_muts[i]), str(high_conf_muts[j])]))
                    if key not in used:
                        used.add(key)
                        candidates.append(self._make_candidate(
                            f"Univ_Pair_{len(candidates)}", cdrs, family,
                            [high_conf_muts[i], high_conf_muts[j]], "pair", 'universal'
                        ))
                if len(candidates) >= n:
                    break
            if len(candidates) >= n:
                break
        
        # Random combinations to fill
        import random
        attempts = 0
        while len(candidates) < n and len(all_muts) >= 2 and attempts < n * 3:
            attempts += 1
            picked = random.sample(all_muts, min(3, len(all_muts)))
            if len({m.imgt_num for m in picked}) == len(picked):
                key = tuple(sorted([str(m) for m in picked]))
                if key not in used:
                    used.add(key)
                    candidates.append(self._make_candidate(
                        f"Univ_Combo_{len(candidates)}", cdrs, family,
                        picked, "combo", 'universal'
                    ))
        
        return candidates[:n]
    
    def _design_original(self, cdrs: CDRSet, family: str, rules: List[dict],
                          archetype: Dict[int, str], triplet_rules: List[dict],
                          n: int, target_hallmarks: str) -> List[VHHCandidate]:
        """Generate candidates keeping original FRs with VHH-izing mutations."""
        candidates = []
        used = set()
        orig = {'FR1': cdrs.fr1, 'FR2': cdrs.fr2, 'FR3': cdrs.fr3, 'FR4': cdrs.fr4}
        
        # Get hallmark mutations
        hallmark_muts = []
        if target_hallmarks in VHH_HALLMARKS_IMGT:
            for imgt_pos, aa in VHH_HALLMARKS_IMGT[target_hallmarks].items():
                mut = self._make_hallmark_mutation(imgt_pos, aa, cdrs, target_hallmarks)
                if mut:
                    hallmark_muts.append(mut)
        
        # Hallmark-only variant (LEAD for original)
        if hallmark_muts:
            lead = self._make_candidate_orig(
                f"Orig_{target_hallmarks}", cdrs, family, hallmark_muts,
                f"hallmark_{target_hallmarks}", orig
            )
            lead.is_lead = True
            candidates.append(lead)
        
        # Get rule-based mutations
        all_rule_muts = []
        for rule in rules:
            mut = self._rule_to_mutation(rule, orig, 'original')
            if mut:
                all_rule_muts.append(mut)
        
        print(f"  {len(hallmark_muts)} hallmark mutations, {len(all_rule_muts)} rule mutations")
        
        # Hallmarks + single rules
        for i, rule_mut in enumerate(all_rule_muts[:min(15, n//3)]):
            muts = hallmark_muts + [rule_mut]
            key = tuple(sorted([str(m) for m in muts]))
            if key not in used:
                used.add(key)
                candidates.append(self._make_candidate_orig(
                    f"Orig_{target_hallmarks}_{i+1}", cdrs, family, muts,
                    f"hallmark_plus_rule", orig
                ))
        
        # Hallmarks + all rules
        if all_rule_muts:
            muts = hallmark_muts + all_rule_muts
            candidates.append(self._make_candidate_orig(
                f"Orig_AllRules", cdrs, family, muts, "all_rules", orig
            ))
        
        # Different hallmark patterns
        for hname in ['YERL', 'FERF', 'FERA', 'FQRL']:
            if len(candidates) >= n:
                break
            if hname == target_hallmarks or hname not in VHH_HALLMARKS_IMGT:
                continue
            
            h_muts = []
            for imgt_pos, aa in VHH_HALLMARKS_IMGT[hname].items():
                mut = self._make_hallmark_mutation(imgt_pos, aa, cdrs, hname)
                if mut:
                    h_muts.append(mut)
            
            if h_muts:
                key = tuple(sorted([str(m) for m in h_muts]))
                if key not in used:
                    used.add(key)
                    candidates.append(self._make_candidate_orig(
                        f"Orig_{hname}", cdrs, family, h_muts, f"hallmark_{hname}", orig
                    ))
        
        # Pairs of rule mutations + hallmarks
        for i in range(min(len(all_rule_muts), 8)):
            for j in range(i+1, min(len(all_rule_muts), 8)):
                if len(candidates) >= n:
                    break
                if all_rule_muts[i].imgt_num != all_rule_muts[j].imgt_num:
                    muts = hallmark_muts + [all_rule_muts[i], all_rule_muts[j]]
                    key = tuple(sorted([str(m) for m in muts]))
                    if key not in used:
                        used.add(key)
                        candidates.append(self._make_candidate_orig(
                            f"Orig_Pair_{len(candidates)}", cdrs, family, muts, "pair", orig
                        ))
        
        # Fill with random combos
        import random
        attempts = 0
        while len(candidates) < n and len(all_rule_muts) >= 2 and attempts < n * 3:
            attempts += 1
            picked = random.sample(all_rule_muts, min(3, len(all_rule_muts)))
            if len({m.imgt_num for m in picked}) == len(picked):
                muts = hallmark_muts + picked
                key = tuple(sorted([str(m) for m in muts]))
                if key not in used:
                    used.add(key)
                    candidates.append(self._make_candidate_orig(
                        f"Orig_Combo_{len(candidates)}", cdrs, family, muts, "combo", orig
                    ))
        
        return candidates[:n]
    
    def _rule_to_mutation(self, rule: dict, frs: dict, mode: str) -> Optional[Mutation]:
        """Convert a rule to a mutation."""
        imgt_num = rule.get('position_num', 0)
        suggested_aa = rule.get('suggested_aa', '')
        
        if not imgt_num or not suggested_aa:
            return None
        
        # Determine FR region and position
        region, idx = self._imgt_to_fr(imgt_num)
        if not region:
            return None
        
        # Skip hallmark positions in universal mode (they're protected)
        if mode == 'universal' and imgt_num in HALLMARK_POSITIONS:
            return None
        
        # Get current AA
        fr_seq = frs.get(region, '')
        if idx >= len(fr_seq):
            return None
        
        original_aa = fr_seq[idx]
        
        # Skip if no change needed
        if original_aa == suggested_aa:
            return None
        
        return Mutation(
            position=f"IMGT{imgt_num}",
            imgt_num=imgt_num,
            original=original_aa,
            mutant=suggested_aa,
            source=rule.get('condition', '')[:20],
            confidence=rule.get('confidence', 0.5)
        )
    
    def _make_hallmark_mutation(self, imgt_pos: int, aa: str, cdrs: CDRSet, source: str) -> Optional[Mutation]:
        """Create a hallmark mutation."""
        # FR2 contains hallmark positions
        fr2 = cdrs.fr2
        if not fr2:
            return None
        
        # Map IMGT to FR2 index
        imgt_to_idx = {42: 3, 49: 10, 50: 11, 52: 13}
        idx = imgt_to_idx.get(imgt_pos)
        
        if idx is None or idx >= len(fr2):
            return None
        
        original = fr2[idx]
        if original == aa:
            return None
        
        return Mutation(
            position=f"IMGT{imgt_pos}",
            imgt_num=imgt_pos,
            original=original,
            mutant=aa,
            source=source,
            confidence=0.95
        )
    
    def _imgt_to_fr(self, imgt_num: int) -> Tuple[Optional[str], int]:
        """Convert IMGT position to FR region and 0-based index."""
        if 1 <= imgt_num <= 26:
            return 'FR1', imgt_num - 1
        elif 39 <= imgt_num <= 55:
            return 'FR2', imgt_num - 39
        elif 66 <= imgt_num <= 104:
            return 'FR3', imgt_num - 66
        elif 118 <= imgt_num <= 128:
            return 'FR4', imgt_num - 118
        return None, 0
    
    def _make_candidate(self, name: str, cdrs: CDRSet, family: str,
                        mutations: List[Mutation], strategy: str, source: str) -> VHHCandidate:
        """Create a candidate using universal scaffold."""
        # Start with scaffold
        fr1 = self.scaffold['FR1']
        fr2 = self.scaffold['FR2']
        fr3 = self.scaffold['FR3']
        fr4 = self.scaffold['FR4']
        
        # Apply mutations
        fr1, fr2, fr3, fr4 = self._apply_mutations(fr1, fr2, fr3, fr4, mutations)
        
        # Build sequence
        sequence = fr1 + cdrs.cdr1 + fr2 + cdrs.cdr2 + fr3 + cdrs.cdr3 + fr4
        
        return VHHCandidate(
            id=name, rank=0, sequence=sequence,
            framework_source=source, family=family,
            cdr1=cdrs.cdr1, cdr2=cdrs.cdr2, cdr3=cdrs.cdr3,
            fr1=fr1, fr2=fr2, fr3=fr3, fr4=fr4,
            mutations=mutations, strategy=strategy
        )
    
    def _make_candidate_orig(self, name: str, cdrs: CDRSet, family: str,
                              mutations: List[Mutation], strategy: str, orig: dict) -> VHHCandidate:
        """Create a candidate using original FRs."""
        fr1 = orig['FR1']
        fr2 = orig['FR2']
        fr3 = orig['FR3']
        fr4 = orig['FR4']
        
        # Apply mutations
        fr1, fr2, fr3, fr4 = self._apply_mutations(fr1, fr2, fr3, fr4, mutations)
        
        # Build sequence
        sequence = fr1 + cdrs.cdr1 + fr2 + cdrs.cdr2 + fr3 + cdrs.cdr3 + fr4
        
        return VHHCandidate(
            id=name, rank=0, sequence=sequence,
            framework_source='original', family=family,
            cdr1=cdrs.cdr1, cdr2=cdrs.cdr2, cdr3=cdrs.cdr3,
            fr1=fr1, fr2=fr2, fr3=fr3, fr4=fr4,
            mutations=mutations, strategy=strategy
        )
    
    def _apply_mutations(self, fr1: str, fr2: str, fr3: str, fr4: str,
                          mutations: List[Mutation]) -> Tuple[str, str, str, str]:
        """Apply mutations to FR sequences."""
        fr1 = list(fr1)
        fr2 = list(fr2)
        fr3 = list(fr3)
        fr4 = list(fr4)
        
        for mut in mutations:
            region, idx = self._imgt_to_fr(mut.imgt_num)
            if region == 'FR1' and idx < len(fr1):
                fr1[idx] = mut.mutant
            elif region == 'FR2' and idx < len(fr2):
                fr2[idx] = mut.mutant
            elif region == 'FR3' and idx < len(fr3):
                fr3[idx] = mut.mutant
            elif region == 'FR4' and idx < len(fr4):
                fr4[idx] = mut.mutant
        
        return ''.join(fr1), ''.join(fr2), ''.join(fr3), ''.join(fr4)

# ============================================================
# OUTPUT
# ============================================================

def to_dataframe(candidates: List[VHHCandidate], input_seq: str, cdrs: CDRSet):
    """Convert candidates to DataFrame."""
    import pandas as pd
    
    rows = []
    for c in candidates:
        # Count mutations from input
        if c.framework_source == 'input':
            n_mutations = 0
        else:
            n_mutations = len(c.mutations)
        
        rows.append({
            'rank': c.rank,
            'id': c.id,
            'sequence': c.sequence,
            'framework_source': c.framework_source,
            'family': c.family,
            'strategy': c.strategy,
            'n_mutations': n_mutations,
            'mutations': ';'.join(str(m) for m in c.mutations),
            'confidence': round(c.confidence, 1),
            'naturalness': round(c.naturalness_score, 1),
            'lead': c.is_lead,
            'cdr1': c.cdr1,
            'cdr2': c.cdr2,
            'cdr3': c.cdr3,
            'fr1': c.fr1,
            'fr2': c.fr2,
            'fr3': c.fr3,
            'fr4': c.fr4,
        })
    
    return pd.DataFrame(rows)

def save_results(df, output_dir: str, base_name: str, mode: str, datetime_str: str):
    """Save results to files."""
    os.makedirs(output_dir, exist_ok=True)
    
    # CSV
    csv_path = os.path.join(output_dir, f"{base_name}_{mode}.csv")
    df.to_csv(csv_path, index=False)
    print(f"\nSaved: {csv_path}")
    
    # FASTA
    fasta_path = os.path.join(output_dir, f"{base_name}_{mode}.fasta")
    with open(fasta_path, 'w') as f:
        for _, row in df.iterrows():
            header = f">{row['id']}|{row['framework_source']}|{row['strategy']}|conf={row['confidence']}"
            f.write(f"{header}\n{row['sequence']}\n")
    print(f"Saved: {fasta_path}")
    
    # Summary
    summary = {
        'datetime': datetime_str,
        'mode': mode,
        'n_candidates': len(df),
        'n_universal': int((df['framework_source'] == 'universal').sum()),
        'n_original': int((df['framework_source'] == 'original').sum()),
        'n_input': int((df['framework_source'] == 'input').sum()),
        'n_leads': int(df['lead'].sum()),
        'family': df['family'].iloc[0] if len(df) > 0 else '',
    }
    
    summary_path = os.path.join(output_dir, f"{base_name}_{mode}_summary.json")
    with open(summary_path, 'w') as f:
        json.dump(summary, f, indent=2)
    print(f"Saved: {summary_path}")

# ============================================================
# MAIN
# ============================================================

def main():
    p = argparse.ArgumentParser(
        description="VHH Designer v7.2 - Compatible with v7.7 Analysis Outputs",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage
  python vhh_designer_v7_2.py -i mouse_hc.fasta \\
      --rules models/analysis/v7.7_*/analysis_rules_v7.json \\
      --archetypes models/analysis/v7.7_*/analysis_vernier_archetypes_v7.json

  # With specific hallmarks
  python vhh_designer_v7_2.py -i input.fasta \\
      --rules analysis_rules_v7.json \\
      --archetypes analysis_vernier_archetypes_v7.json \\
      --target-hallmarks FERG \\
      --mode both
"""
    )
    
    p.add_argument('--sequence', '-s', help='Input sequence')
    p.add_argument('--input', '-i', help='Input FASTA file')
    p.add_argument('--rules', '-r', required=True,
                   help='analysis_rules_v7.json from v7.7 analysis')
    p.add_argument('--archetypes', '-a', required=True,
                   help='analysis_vernier_archetypes_v7.json from v7.7 analysis')
    p.add_argument('--target-hallmarks', '-t',
                   help='Target hallmark pattern (e.g., FERG, YERL)')
    p.add_argument('--output-dir', '-o', help='Output directory')
    p.add_argument('--name', help='Sequence name (auto-detected from FASTA)')
    p.add_argument('--n-candidates', '-n', type=int, default=92,
                   help='Number of candidates per mode (default: 92)')
    p.add_argument('--mode', '-m', choices=['both', 'universal', 'original'], default='both',
                   help='Design mode (default: both)')
    
    args = p.parse_args()
    
    # Get input sequence
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
    print("VHH DESIGNER v7.2")
    print("=" * 70)
    print(f"Sequence: {seq_name}")
    print(f"Mode: {args.mode.upper()}")
    print(f"Target hallmarks: {args.target_hallmarks or 'auto'}")
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
    
    # Show input hallmarks
    if cdrs.fr2:
        hallmarks = get_hallmarks_from_fr2(cdrs.fr2)
        print(f"\nInput hallmarks: {hallmarks}")
    
    # Design
    designer = VHHDesignerV7_2(args.rules, args.archetypes)
    candidates = designer.design(
        cdrs, args.n_candidates, args.mode,
        input_sequence=input_seq,
        target_hallmarks=args.target_hallmarks
    )
    
    df = to_dataframe(candidates, input_seq, cdrs)
    
    print(f"\n{'='*70}")
    print(f"RESULTS")
    print(f"{'='*70}")
    print(f"Total: {len(df)} | Leads: {df['lead'].sum()} | Family: {df['family'].iloc[0]}")
    print(f"  Input: {(df['framework_source']=='input').sum()}")
    print(f"  Universal: {(df['framework_source']=='universal').sum()}")
    print(f"  Original: {(df['framework_source']=='original').sum()}")
    
    print("\nLeads:")
    for _, r in df[df['lead']].iterrows():
        print(f"  {r['rank']:2d}. {r['id'][:40]:40s} [{r['framework_source'][:4]}] conf={r['confidence']:.1f}")
    
    # Save
    datetime_str = datetime.now().strftime('%Y%m%d_%H%M%S')
    script_name = "vhh_designer_v7_2"
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
