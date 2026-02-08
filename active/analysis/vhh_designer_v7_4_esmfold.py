#!/usr/bin/env python3
"""
VHH Designer v7.4 - ESMFold Structure Prediction & Comprehensive Scoring
=========================================================================

This version includes:
1. Candidate generation (N candidates, default 500)
2. ESM2 language model scoring (pseudo-perplexity)
3. ESMFold structure prediction (pLDDT confidence)
4. Multi-family probabilistic rule compliance scoring
5. Comprehensive CSV output with all metrics
6. Lead candidates preserved and kept at top

Output CSV includes:
  - generation_order: Original sequence number (1, 2, 3...)
  - detected_family: Family classification
  - vernier_matches: Number of vernier position agreements
  - rules_passed / rules_total: Rule compliance counts
  - framework_identity_pct: % of original framework retained
  - construction_method: How the candidate was built
  - esm2_loss, plddt_mean, plddt_cdr3, combined_score, etc.

Usage:
  python vhh_designer_v7_4_esmfold.py -i input.fasta \\
      --rules analysis_rules_v7.json \\
      --archetypes analysis_vernier_archetypes_v7.json \\
      --n-generate 500 \\
      --n-select 92

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

# ESMFold availability check
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
# CONSTANTS
# ============================================================

UNIVERSAL_SCAFFOLD = {
    'name': 'Universal Humanized VHH',
    'FR1': 'QVQLVESGGGLVQPGGSLRLSCAASG',
    'FR2': 'WFRQAPGQGLEAVA',
    'FR3': 'YYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYC',
    'FR4': 'WGQGTLVTVSS',
}

IMGT_REGIONS = {
    'FR1': (1, 26),
    'FR2': (39, 55),
    'FR3': (66, 104),
    'FR4': (118, 128),
}

HALLMARK_POSITIONS = [42, 49, 50, 52]

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

# Vernier positions (FR3)
VERNIER_POSITIONS = {66, 67, 68, 69, 71, 76, 78, 82, 87, 89, 91, 94}

# All vernier including FR2 hallmarks
ALL_VERNIER_POSITIONS = {42, 49, 50, 52, 66, 67, 68, 69, 71, 76, 78, 82, 87, 89, 91, 94}

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
    position: str
    imgt_num: int
    original: str
    mutant: str
    source: str
    confidence: float
    
    def __str__(self):
        return f"{self.position}:{self.original}→{self.mutant}"

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
    mutations: List[Mutation] = field(default_factory=list)
    strategy: str = ""
    is_lead: bool = False
    generation_order: int = 0  # Original sequence number
    construction_method: str = ""  # Detailed construction description
    target_family: str = ""  # Which family this was designed for
    framework_identity_pct: float = 100.0  # % of original framework retained
    scoring: ScoringResult = field(default_factory=ScoringResult)

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
# FAMILY CLASSIFICATION & FEATURES
# ============================================================

def classify_family(cdrs: CDRSet, full_sequence: str = None) -> str:
    """Classify VHH family based on hallmarks and cysteine count."""
    fr2 = cdrs.fr2 if cdrs.fr2 else ""
    
    p42 = fr2[3] if len(fr2) > 3 else "-"
    p49 = fr2[10] if len(fr2) > 10 else "-"
    p50 = fr2[11] if len(fr2) > 11 else "-"
    p52 = fr2[13] if len(fr2) > 13 else "-"
    
    if p50 == "L":
        return "VH_like"
    if p52 == "W":
        return "VHH_W52"
    
    seq = full_sequence or (cdrs.fr1 + cdrs.cdr1 + cdrs.fr2 + cdrs.cdr2 + cdrs.fr3 + cdrs.cdr3 + cdrs.fr4)
    cys_count = seq.count("C")
    c_label = "C4" if cys_count >= 4 else "C2"
    
    if p42 == "F":
        return f"F_{c_label}"
    elif p42 == "Y":
        return f"Y_{c_label}"
    
    return "Other_VHH"

def get_hallmarks_from_fr2(fr2: str) -> str:
    """Extract hallmark pattern from FR2."""
    if not fr2 or len(fr2) < 14:
        return "????"
    p42 = fr2[3] if len(fr2) > 3 else "?"
    p49 = fr2[10] if len(fr2) > 10 else "?"
    p50 = fr2[11] if len(fr2) > 11 else "?"
    p52 = fr2[13] if len(fr2) > 13 else "?"
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

def calculate_framework_identity(candidate_frs: dict, original_frs: dict) -> float:
    """Calculate % identity between candidate and original frameworks."""
    total_len = 0
    matches = 0
    
    for region in ['FR1', 'FR2', 'FR3', 'FR4']:
        cand = candidate_frs.get(region, '')
        orig = original_frs.get(region, '')
        
        if not orig:
            continue
        
        min_len = min(len(cand), len(orig))
        total_len += len(orig)
        
        for i in range(min_len):
            if cand[i] == orig[i]:
                matches += 1
    
    if total_len == 0:
        return 100.0
    
    return (matches / total_len) * 100.0

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
        """Load model (lazy loading)."""
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
        """Score multiple sequences. Returns list of (loss, perplexity) tuples."""
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
                    print(f"    ESM2 error on sequence: {e}")
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
        """Load ESMFold model (lazy loading)."""
        if self.model is not None:
            return True
        
        if not ESMFOLD_AVAILABLE or not TORCH_AVAILABLE:
            print("  ESMFold not available")
            return False
        
        print("  Loading ESMFold model (this may take a moment)...")
        try:
            self.model = esm.pretrained.esmfold_v1()
            self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
            self.model = self.model.to(self.device)
            self.model.eval()
            
            # Set chunk size for memory efficiency
            self.model.set_chunk_size(128)
            
            print(f"    Device: {self.device}")
            return True
        except Exception as e:
            print(f"    Failed to load ESMFold: {e}")
            return False
    
    def predict_plddt(self, sequence: str) -> Dict[str, float]:
        """
        Predict pLDDT scores for a sequence.
        
        Returns dict with plddt_mean, plddt_median, plddt_min, and per-residue info.
        """
        if not self.load():
            return {
                'plddt_mean': 0.0,
                'plddt_median': 0.0,
                'plddt_min': 0.0,
                'plddt_scores': []
            }
        
        try:
            with torch.no_grad():
                output = self.model.infer_pdb(sequence)
            
            # Parse pLDDT from PDB output (B-factor column)
            plddt_scores = []
            for line in output.split('\n'):
                if line.startswith('ATOM'):
                    try:
                        # B-factor is columns 61-66 in PDB format
                        bfactor = float(line[60:66].strip())
                        plddt_scores.append(bfactor)
                    except:
                        continue
            
            if not plddt_scores:
                return {
                    'plddt_mean': 0.0,
                    'plddt_median': 0.0,
                    'plddt_min': 0.0,
                    'plddt_scores': []
                }
            
            # Get per-residue pLDDT (average over atoms per residue)
            # Group by residue
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
            
            # Don't forget last residue
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
            return {
                'plddt_mean': 0.0,
                'plddt_median': 0.0,
                'plddt_min': 0.0,
                'plddt_scores': []
            }
    
    def predict_batch(self, sequences: List[str], 
                      region_ranges: List[Dict[str, Tuple[int, int]]] = None) -> List[Dict[str, float]]:
        """
        Predict pLDDT for multiple sequences.
        
        Args:
            sequences: List of sequences
            region_ranges: List of dicts with region start/end indices for each sequence
                          e.g., [{'cdr1': (26, 34), 'cdr2': (50, 58), 'cdr3': (95, 110), 'framework': [(0,26), (34,50), (58,95), (110,121)]}]
        
        Returns list of pLDDT result dicts.
        """
        results = []
        
        iterator = enumerate(sequences)
        if TQDM_AVAILABLE:
            iterator = tqdm(list(iterator), desc="ESMFold prediction", leave=False)
        
        for i, seq in iterator:
            result = self.predict_plddt(seq)
            
            # Calculate region-specific pLDDT if ranges provided
            if region_ranges and i < len(region_ranges) and result['plddt_scores']:
                ranges = region_ranges[i]
                plddt_scores = result['plddt_scores']
                
                # CDR pLDDT
                for region in ['cdr1', 'cdr2', 'cdr3']:
                    if region in ranges:
                        start, end = ranges[region]
                        end = min(end, len(plddt_scores))
                        if start < end:
                            region_scores = plddt_scores[start:end]
                            result[f'plddt_{region}'] = np.mean(region_scores) if region_scores else 0.0
                        else:
                            result[f'plddt_{region}'] = 0.0
                    else:
                        result[f'plddt_{region}'] = 0.0
                
                # Framework pLDDT
                if 'framework_ranges' in ranges:
                    fw_scores = []
                    for start, end in ranges['framework_ranges']:
                        end = min(end, len(plddt_scores))
                        if start < end:
                            fw_scores.extend(plddt_scores[start:end])
                    result['plddt_framework'] = np.mean(fw_scores) if fw_scores else 0.0
                else:
                    result['plddt_framework'] = result['plddt_mean']
            
            results.append(result)
        
        return results

# ============================================================
# MULTI-FAMILY PROBABILISTIC SCORER
# ============================================================

class MultiFamilyProbabilisticScorer:
    """Score sequences using multi-family probabilistic rule compliance."""
    
    def __init__(self, rules_file: str, archetypes_file: str):
        print("  Loading rules engine...")
        with open(rules_file, 'r') as f:
            self.rules = json.load(f)
        
        print("  Loading archetypes...")
        with open(archetypes_file, 'r') as f:
            archetypes = json.load(f)
        
        # Index rules by family
        self.rules_by_family = defaultdict(list)
        for rule in self.rules:
            family = rule.get('family', '')
            if family:
                self.rules_by_family[family].append(rule)
        
        # Index archetypes by family
        self.archetypes = {a['family']: a for a in archetypes}
        self.families = list(self.archetypes.keys())
        
        print(f"    {len(self.rules)} rules across {len(self.families)} families")
    
    def count_vernier_matches(self, candidate: VHHCandidate, family: str) -> Tuple[int, int]:
        """Count how many vernier positions match the family archetype."""
        archetype = self.archetypes.get(family, {})
        positions = archetype.get('positions', {})
        
        matches = 0
        total = 0
        
        for pos_str, data in positions.items():
            try:
                imgt_num = int(pos_str.replace('IMGT', ''))
                if imgt_num not in ALL_VERNIER_POSITIONS:
                    continue
                
                consensus_aa = data.get('consensus', '')
                actual_aa = self._get_aa_at_imgt(candidate, imgt_num)
                
                if consensus_aa and actual_aa:
                    total += 1
                    if actual_aa == consensus_aa:
                        matches += 1
            except ValueError:
                continue
        
        return matches, total
    
    def compute_family_probability(self, candidate: VHHCandidate, cdr3_features: Dict[str, str]) -> Dict[str, float]:
        """Compute P(family | sequence) based on vernier matches + CDR compatibility."""
        scores = {}
        
        for family, archetype in self.archetypes.items():
            if family == 'VH_like':
                continue
            
            score = 0.0
            total_weight = 0.0
            
            # 1. Vernier position matches
            positions = archetype.get('positions', {})
            for pos_str, data in positions.items():
                try:
                    imgt_num = int(pos_str.replace('IMGT', ''))
                    consensus_aa = data.get('consensus', '')
                    actual_aa = self._get_aa_at_imgt(candidate, imgt_num)
                    
                    if actual_aa and consensus_aa:
                        weight = 2.0 if imgt_num in VERNIER_POSITIONS else 1.0
                        total_weight += weight
                        if actual_aa == consensus_aa:
                            score += weight
                except ValueError:
                    continue
            
            # 2. CDR3 length compatibility
            cdr3_len = cdr3_features.get('cdr3_len', 'medium')
            total_weight += 1.0
            if family.startswith('F_'):
                score += 1.0 if cdr3_len in ['short', 'medium'] else 0.5
            elif family.startswith('Y_'):
                score += 1.0 if cdr3_len in ['medium', 'long'] else 0.5
            else:
                score += 0.7
            
            # 3. Hallmark compatibility
            hallmarks = get_hallmarks_from_fr2(candidate.fr2)
            total_weight += 2.0
            if family.startswith('F_') and hallmarks[0] == 'F':
                score += 2.0
            elif family.startswith('Y_') and hallmarks[0] == 'Y':
                score += 2.0
            elif family == 'Other_VHH':
                score += 1.0
            
            if total_weight > 0:
                scores[family] = score / total_weight
            else:
                scores[family] = 0.0
        
        # Normalize
        total = sum(scores.values())
        if total > 0:
            scores = {k: v / total for k, v in scores.items()}
        
        return scores
    
    def compute_rule_compliance(self, candidate: VHHCandidate, family: str, 
                                 cdr3_features: Dict[str, str],
                                 min_confidence: float = 0.5) -> Tuple[float, int, int, int, List[str]]:
        """
        Compute rule compliance for a specific family.
        
        Returns (compliance_score, rules_passed, rules_total, rules_applicable, violations).
        """
        rules = self.rules_by_family.get(family, [])
        if not rules:
            return 1.0, 0, 0, 0, []
        
        hallmarks = get_hallmarks_from_fr2(candidate.fr2)
        
        satisfied_weight = 0.0
        total_weight = 0.0
        rules_passed = 0
        rules_applicable = 0
        violations = []
        
        for rule in rules:
            confidence = rule.get('confidence', 0)
            if confidence < min_confidence:
                continue
            
            # Check if rule condition applies
            condition = rule.get('condition', '')
            if not self._condition_applies(condition, hallmarks, cdr3_features):
                continue
            
            rules_applicable += 1
            
            # Check if rule is satisfied
            position_num = rule.get('position_num', 0)
            suggested_aa = rule.get('suggested_aa', '')
            
            if position_num and suggested_aa:
                actual_aa = self._get_aa_at_imgt(candidate, position_num)
                
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
        """Score a candidate using multi-family probabilistic approach."""
        result = ScoringResult()
        
        # Family probabilities
        family_probs = self.compute_family_probability(candidate, cdr3_features)
        result.family_probabilities = family_probs
        
        # Get top family for vernier counting
        top_family = max(family_probs.items(), key=lambda x: x[1])[0] if family_probs else 'F_C2'
        result.vernier_matches, result.vernier_total = self.count_vernier_matches(candidate, top_family)
        
        # Rule compliance for relevant families
        total_rules_passed = 0
        total_rules = 0
        total_applicable = 0
        all_violations = []
        
        for family, prob in family_probs.items():
            if prob < family_threshold:
                continue
            
            compliance, passed, total, applicable, violations = self.compute_rule_compliance(
                candidate, family, cdr3_features, rule_min_confidence
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
        
        return result
    
    def _condition_applies(self, condition: str, hallmarks: str, cdr3_features: Dict[str, str]) -> bool:
        """Check if a rule condition applies."""
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
    
    def _get_aa_at_imgt(self, candidate: VHHCandidate, imgt_num: int) -> str:
        """Get amino acid at IMGT position from candidate."""
        if 1 <= imgt_num <= 26:
            idx = imgt_num - 1
            return candidate.fr1[idx] if idx < len(candidate.fr1) else ""
        elif 39 <= imgt_num <= 55:
            idx = imgt_num - 39
            return candidate.fr2[idx] if idx < len(candidate.fr2) else ""
        elif 66 <= imgt_num <= 104:
            idx = imgt_num - 66
            return candidate.fr3[idx] if idx < len(candidate.fr3) else ""
        elif 118 <= imgt_num <= 128:
            idx = imgt_num - 118
            return candidate.fr4[idx] if idx < len(candidate.fr4) else ""
        return ""

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
                 esmfold_weight: float = 0.3,
                 rule_weight: float = 0.5):
        
        self.multi_family_scorer = MultiFamilyProbabilisticScorer(rules_file, archetypes_file)
        
        self.use_esm2 = use_esm2 and ESM2_AVAILABLE and TORCH_AVAILABLE
        self.use_esmfold = use_esmfold and ESMFOLD_AVAILABLE and TORCH_AVAILABLE
        
        self.esm2_weight = esm2_weight
        self.esmfold_weight = esmfold_weight
        self.rule_weight = rule_weight
        
        # Normalize weights
        total_weight = 0
        if self.use_esm2:
            total_weight += esm2_weight
        if self.use_esmfold:
            total_weight += esmfold_weight
        total_weight += rule_weight
        
        if total_weight > 0:
            self.esm2_weight = esm2_weight / total_weight
            self.esmfold_weight = esmfold_weight / total_weight
            self.rule_weight = rule_weight / total_weight
        
        if self.use_esm2:
            self.esm2_scorer = ESM2Scorer(esm_model)
        else:
            self.esm2_scorer = None
        
        if self.use_esmfold:
            self.esmfold_predictor = ESMFoldPredictor()
        else:
            self.esmfold_predictor = None
        
        print(f"  Scoring weights: ESM2={self.esm2_weight:.2f}, ESMFold={self.esmfold_weight:.2f}, Rules={self.rule_weight:.2f}")
    
    def score_candidates(self, candidates: List[VHHCandidate], cdr3_features: Dict[str, str],
                         original_frs: Dict[str, str] = None) -> List[VHHCandidate]:
        """Score all candidates and update their scoring results."""
        print(f"\nScoring {len(candidates)} candidates...")
        
        # Calculate framework identity for all candidates
        if original_frs:
            for candidate in candidates:
                cand_frs = {
                    'FR1': candidate.fr1,
                    'FR2': candidate.fr2,
                    'FR3': candidate.fr3,
                    'FR4': candidate.fr4
                }
                candidate.framework_identity_pct = calculate_framework_identity(cand_frs, original_frs)
        
        # Step 1: Multi-family probabilistic scoring
        print("  Computing rule compliance scores...")
        for candidate in tqdm(candidates, desc="Rule compliance", disable=not TQDM_AVAILABLE):
            result = self.multi_family_scorer.score(candidate, cdr3_features)
            candidate.scoring = result
        
        # Step 2: ESM2 scoring
        if self.use_esm2:
            print("  Computing ESM2 scores...")
            sequences = [c.sequence for c in candidates]
            esm2_scores = self.esm2_scorer.score_batch(sequences)
            
            for candidate, (loss, perplexity) in zip(candidates, esm2_scores):
                candidate.scoring.esm2_loss = loss
                candidate.scoring.esm2_perplexity = perplexity
        
        # Step 3: ESMFold structure prediction
        if self.use_esmfold:
            print("  Computing ESMFold pLDDT scores...")
            sequences = [c.sequence for c in candidates]
            
            # Build region ranges for each candidate
            region_ranges = []
            for c in candidates:
                fr1_len = len(c.fr1)
                cdr1_len = len(c.cdr1)
                fr2_len = len(c.fr2)
                cdr2_len = len(c.cdr2)
                fr3_len = len(c.fr3)
                cdr3_len = len(c.cdr3)
                fr4_len = len(c.fr4)
                
                cdr1_start = fr1_len
                cdr1_end = cdr1_start + cdr1_len
                fr2_end = cdr1_end + fr2_len
                cdr2_start = fr2_end
                cdr2_end = cdr2_start + cdr2_len
                fr3_end = cdr2_end + fr3_len
                cdr3_start = fr3_end
                cdr3_end = cdr3_start + cdr3_len
                
                region_ranges.append({
                    'cdr1': (cdr1_start, cdr1_end),
                    'cdr2': (cdr2_start, cdr2_end),
                    'cdr3': (cdr3_start, cdr3_end),
                    'framework_ranges': [(0, fr1_len), (cdr1_end, fr2_end), (cdr2_end, fr3_end), (cdr3_end, cdr3_end + fr4_len)]
                })
            
            plddt_results = self.esmfold_predictor.predict_batch(sequences, region_ranges)
            
            for candidate, plddt_result in zip(candidates, plddt_results):
                candidate.scoring.plddt_mean = plddt_result.get('plddt_mean', 0.0)
                candidate.scoring.plddt_median = plddt_result.get('plddt_median', 0.0)
                candidate.scoring.plddt_min = plddt_result.get('plddt_min', 0.0)
                candidate.scoring.plddt_cdr1 = plddt_result.get('plddt_cdr1', 0.0)
                candidate.scoring.plddt_cdr2 = plddt_result.get('plddt_cdr2', 0.0)
                candidate.scoring.plddt_cdr3 = plddt_result.get('plddt_cdr3', 0.0)
                candidate.scoring.plddt_framework = plddt_result.get('plddt_framework', 0.0)
        
        # Step 4: Compute combined scores
        print("  Computing combined scores...")
        
        # Normalize ESM2 scores (lower loss = better, so invert)
        if self.use_esm2:
            esm2_losses = [c.scoring.esm2_loss for c in candidates if c.scoring.esm2_loss > 0]
            if esm2_losses:
                esm2_min, esm2_max = min(esm2_losses), max(esm2_losses)
                esm2_range = esm2_max - esm2_min if esm2_max > esm2_min else 1.0
            else:
                esm2_min, esm2_range = 0, 1
        
        # Normalize pLDDT scores (0-100, higher = better)
        if self.use_esmfold:
            plddt_max = 100.0
        
        for candidate in candidates:
            # Skip scoring for leads - they stay at top regardless
            if candidate.is_lead:
                candidate.scoring.combined_score = 999.0  # Ensures leads stay at top
                continue
            
            # Rule compliance score (already 0-1)
            rule_score = candidate.scoring.weighted_naturalness
            
            # ESM2 score (normalize and invert)
            if self.use_esm2 and candidate.scoring.esm2_loss > 0:
                esm2_normalized = 1.0 - (candidate.scoring.esm2_loss - esm2_min) / esm2_range
            else:
                esm2_normalized = 0.5
            
            # ESMFold pLDDT score (normalize to 0-1)
            if self.use_esmfold and candidate.scoring.plddt_mean > 0:
                plddt_normalized = candidate.scoring.plddt_mean / plddt_max
            else:
                plddt_normalized = 0.5
            
            # Combined score
            combined = (
                self.esm2_weight * esm2_normalized +
                self.esmfold_weight * plddt_normalized +
                self.rule_weight * rule_score
            )
            
            candidate.scoring.combined_score = combined
        
        return candidates

# ============================================================
# CANDIDATE GENERATOR
# ============================================================

class CandidateGenerator:
    """Generate VHH candidates using rules and archetypes."""
    
    def __init__(self, rules_file: str, archetypes_file: str):
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
        self.scaffold = UNIVERSAL_SCAFFOLD
        
        print(f"  {len(self.compensation_rules)} compensation rules")
        print(f"  {len(self.triplet_rules)} triplet rules")
        print(f"  {len(self.archetypes)} family archetypes")
    
    def generate(self, cdrs: CDRSet, n_candidates: int, 
                 target_hallmarks: str = 'FERG',
                 target_families: List[str] = None,
                 mode: str = 'both',
                 input_sequence: str = None) -> List[VHHCandidate]:
        """Generate N candidates with generation_order tracking."""
        candidates = []
        generation_order = 0
        
        if not target_families:
            target_families = ['F_C2', 'Y_C2', 'Other_VHH']
        
        cdr3_features = get_cdr3_features(cdrs.cdr3)
        
        print(f"\nGenerating {n_candidates} candidates...")
        print(f"  Target families: {target_families}")
        print(f"  Target hallmarks: {target_hallmarks}")
        
        # Add original input as LEAD (untouched)
        input_family = classify_family(cdrs)
        input_seq = input_sequence or (cdrs.fr1 + cdrs.cdr1 + cdrs.fr2 + cdrs.cdr2 + cdrs.fr3 + cdrs.cdr3 + cdrs.fr4)
        
        generation_order += 1
        lead = VHHCandidate(
            id="Lead_Input_Original",
            rank=0,
            sequence=input_seq,
            framework_source='input',
            family=input_family,
            cdr1=cdrs.cdr1, cdr2=cdrs.cdr2, cdr3=cdrs.cdr3,
            fr1=cdrs.fr1, fr2=cdrs.fr2, fr3=cdrs.fr3, fr4=cdrs.fr4,
            mutations=[],
            strategy="original_input",
            is_lead=True,
            generation_order=generation_order,
            construction_method="Original input sequence (untouched)",
            target_family=input_family,
            framework_identity_pct=100.0
        )
        candidates.append(lead)
        
        # Calculate candidates per family
        n_remaining = n_candidates - 1
        n_per_family = n_remaining // len(target_families)
        
        for i, family in enumerate(target_families):
            n_for_family = n_per_family
            if i == len(target_families) - 1:
                n_for_family = n_remaining - (len(candidates) - 1)
            
            if n_for_family <= 0:
                continue
            
            # Choose hallmarks
            if family.startswith('Y_'):
                fam_hallmarks = 'YERL'
            else:
                fam_hallmarks = target_hallmarks
            
            # Get rules and archetype
            rules = self._get_rules(fam_hallmarks, family, cdr3_features)
            archetype = self._get_archetype(family)
            triplet_rules = self._get_triplet_rules(family, cdrs.cdr3)
            
            print(f"  {family}: {n_for_family} candidates, {len(rules)} rules")
            
            # Generate
            if mode in ['universal', 'both']:
                n_univ = n_for_family // 2 if mode == 'both' else n_for_family
                new_candidates, generation_order = self._generate_universal(
                    cdrs, family, rules, archetype, triplet_rules, n_univ, 
                    fam_hallmarks, generation_order
                )
                candidates.extend(new_candidates)
            
            if mode in ['original', 'both'] and cdrs.fr1:
                n_orig = n_for_family - (n_for_family // 2 if mode == 'both' else 0)
                new_candidates, generation_order = self._generate_original(
                    cdrs, family, rules, archetype, triplet_rules, n_orig,
                    fam_hallmarks, generation_order
                )
                candidates.extend(new_candidates)
        
        return candidates
    
    def _get_rules(self, hallmarks: str, family: str, cdr3_features: Dict[str, str],
                   min_support: int = 500, min_confidence: float = 0.6) -> List[dict]:
        """Get applicable rules."""
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
    
    def _generate_universal(self, cdrs: CDRSet, family: str, rules: List[dict],
                             archetype: Dict[int, str], triplet_rules: List[dict],
                             n: int, target_hallmarks: str, 
                             generation_order: int) -> Tuple[List[VHHCandidate], int]:
        """Generate candidates using universal scaffold."""
        candidates = []
        used = set()
        
        hallmark_set = set(HALLMARK_POSITIONS)
        
        # Get mutations
        vernier_muts, other_arch_muts = self._get_archetype_mutations(archetype, 'universal')
        rule_muts = self._get_rule_mutations(rules, 'universal')
        triplet_mut = self._get_triplet_mutation(triplet_rules, self.scaffold['FR4'])
        
        # Strategy combinations with descriptions
        strategies = [
            ("graft", [], "Universal scaffold + CDR graft only"),
            ("vernier", vernier_muts, f"Universal + {len(vernier_muts)} vernier archetype mutations"),
            ("full_arch", vernier_muts + other_arch_muts, f"Universal + full {family} archetype ({len(vernier_muts + other_arch_muts)} mutations)"),
            ("vernier_rules", vernier_muts + rule_muts, f"Universal + vernier + {len(rule_muts)} CDR3-context rules"),
            ("full_optimized", vernier_muts + other_arch_muts + rule_muts, f"Universal + full archetype + all rules"),
        ]
        
        if triplet_mut:
            strategies.append(("vernier_triplet", vernier_muts + [triplet_mut], f"Universal + vernier + FR4 junction ({triplet_mut.mutant})"))
            strategies.append(("full_triplet", vernier_muts + other_arch_muts + rule_muts + [triplet_mut], "Universal + full optimization + FR4 junction"))
        
        for strategy_name, mutations, description in strategies:
            if len(candidates) >= n:
                break
            
            key = tuple(sorted([str(m) for m in mutations]))
            if key in used:
                continue
            used.add(key)
            
            generation_order += 1
            candidate = self._make_universal_candidate(
                f"Univ_{family}_{strategy_name}_{generation_order}", cdrs, family, 
                mutations, strategy_name, generation_order, description, family
            )
            candidates.append(candidate)
        
        # Fill with random combinations
        import random
        all_muts = vernier_muts + other_arch_muts + rule_muts
        if triplet_mut:
            all_muts.append(triplet_mut)
        
        attempts = 0
        while len(candidates) < n and len(all_muts) >= 2 and attempts < n * 3:
            attempts += 1
            n_pick = random.randint(2, min(5, len(all_muts)))
            picked = random.sample(all_muts, n_pick)
            
            if len({m.imgt_num for m in picked}) == len(picked):
                key = tuple(sorted([str(m) for m in picked]))
                if key not in used:
                    used.add(key)
                    generation_order += 1
                    
                    description = f"Universal + random combo of {len(picked)} mutations"
                    candidate = self._make_universal_candidate(
                        f"Univ_{family}_combo_{generation_order}", cdrs, family, 
                        picked, "random_combo", generation_order, description, family
                    )
                    candidates.append(candidate)
        
        return candidates[:n], generation_order
    
    def _generate_original(self, cdrs: CDRSet, family: str, rules: List[dict],
                            archetype: Dict[int, str], triplet_rules: List[dict],
                            n: int, target_hallmarks: str,
                            generation_order: int) -> Tuple[List[VHHCandidate], int]:
        """Generate candidates using original FRs with VHH-izing mutations."""
        candidates = []
        used = set()
        
        orig = {'FR1': cdrs.fr1, 'FR2': cdrs.fr2, 'FR3': cdrs.fr3, 'FR4': cdrs.fr4}
        
        # Get hallmark mutations (required for VHH)
        hallmark_muts = self._get_hallmark_mutations(target_hallmarks, cdrs)
        n_hallmark = len(hallmark_muts)
        
        # Get other mutations
        vernier_muts, other_arch_muts = self._get_archetype_mutations(archetype, 'original', orig)
        rule_muts = self._get_rule_mutations(rules, 'original', orig)
        triplet_mut = self._get_triplet_mutation(triplet_rules, cdrs.fr4)
        
        # Strategies with descriptions
        strategies = [
            ("hallmark_only", hallmark_muts, f"Original FRs + {n_hallmark} hallmark mutations ({target_hallmarks})"),
            ("vernier", hallmark_muts + vernier_muts, f"Original + hallmarks + {len(vernier_muts)} vernier mutations"),
            ("full_arch", hallmark_muts + vernier_muts + other_arch_muts, f"Original + hallmarks + full {family} archetype"),
            ("vernier_rules", hallmark_muts + vernier_muts + rule_muts, f"Original + hallmarks + vernier + {len(rule_muts)} rules"),
            ("full_optimized", hallmark_muts + vernier_muts + other_arch_muts + rule_muts, "Original + hallmarks + full optimization"),
        ]
        
        if triplet_mut:
            strategies.append(("vernier_triplet", hallmark_muts + vernier_muts + [triplet_mut], f"Original + hallmarks + vernier + FR4 junction"))
            strategies.append(("full_triplet", hallmark_muts + vernier_muts + other_arch_muts + rule_muts + [triplet_mut], "Original + full optimization + FR4 junction"))
        
        for strategy_name, mutations, description in strategies:
            if len(candidates) >= n:
                break
            
            key = tuple(sorted([str(m) for m in mutations]))
            if key in used:
                continue
            used.add(key)
            
            generation_order += 1
            candidate = self._make_original_candidate(
                f"Orig_{family}_{strategy_name}_{generation_order}", cdrs, family,
                mutations, strategy_name, orig, generation_order, description, family
            )
            candidates.append(candidate)
        
        # Fill with combinations
        import random
        all_optional = vernier_muts + other_arch_muts + rule_muts
        if triplet_mut:
            all_optional.append(triplet_mut)
        
        attempts = 0
        while len(candidates) < n and len(all_optional) >= 2 and attempts < n * 3:
            attempts += 1
            n_pick = random.randint(2, min(5, len(all_optional)))
            picked = random.sample(all_optional, n_pick)
            
            if len({m.imgt_num for m in picked}) == len(picked):
                muts = hallmark_muts + picked
                key = tuple(sorted([str(m) for m in muts]))
                if key not in used:
                    used.add(key)
                    generation_order += 1
                    
                    description = f"Original + hallmarks + random combo of {len(picked)} mutations"
                    candidate = self._make_original_candidate(
                        f"Orig_{family}_combo_{generation_order}", cdrs, family,
                        muts, "random_combo", orig, generation_order, description, family
                    )
                    candidates.append(candidate)
        
        return candidates[:n], generation_order
    
    def _get_archetype_mutations(self, archetype: Dict[int, str], mode: str,
                                   orig: dict = None) -> Tuple[List[Mutation], List[Mutation]]:
        vernier_muts = []
        other_muts = []
        hallmark_set = set(HALLMARK_POSITIONS)
        
        for imgt_num, consensus_aa in archetype.items():
            if imgt_num in hallmark_set:
                continue
            
            region, idx = self._imgt_to_fr(imgt_num)
            if not region:
                continue
            
            if mode == 'universal':
                fr_seq = self.scaffold.get(region, '')
            else:
                fr_seq = orig.get(region, '')
            
            if idx >= len(fr_seq):
                continue
            
            original_aa = fr_seq[idx]
            if original_aa == consensus_aa:
                continue
            
            mut = Mutation(
                position=f"IMGT{imgt_num}",
                imgt_num=imgt_num,
                original=original_aa,
                mutant=consensus_aa,
                source="archetype",
                confidence=0.7
            )
            
            if imgt_num in VERNIER_POSITIONS:
                vernier_muts.append(mut)
            else:
                other_muts.append(mut)
        
        return vernier_muts, other_muts
    
    def _get_rule_mutations(self, rules: List[dict], mode: str, orig: dict = None) -> List[Mutation]:
        mutations = []
        hallmark_set = set(HALLMARK_POSITIONS)
        
        for rule in rules:
            imgt_num = rule.get('position_num', 0)
            suggested_aa = rule.get('suggested_aa', '')
            
            if not imgt_num or not suggested_aa:
                continue
            if imgt_num in hallmark_set:
                continue
            
            region, idx = self._imgt_to_fr(imgt_num)
            if not region:
                continue
            
            if mode == 'universal':
                fr_seq = self.scaffold.get(region, '')
            else:
                fr_seq = orig.get(region, '')
            
            if idx >= len(fr_seq):
                continue
            
            original_aa = fr_seq[idx]
            if original_aa == suggested_aa:
                continue
            
            mutations.append(Mutation(
                position=f"IMGT{imgt_num}",
                imgt_num=imgt_num,
                original=original_aa,
                mutant=suggested_aa,
                source=rule.get('condition', '')[:20],
                confidence=rule.get('confidence', 0.5)
            ))
        
        return mutations
    
    def _get_hallmark_mutations(self, target_hallmarks: str, cdrs: CDRSet) -> List[Mutation]:
        mutations = []
        
        if target_hallmarks not in VHH_HALLMARKS_IMGT:
            return mutations
        
        fr2 = cdrs.fr2
        if not fr2:
            return mutations
        
        imgt_to_idx = {42: 3, 49: 10, 50: 11, 52: 13}
        
        for imgt_pos, target_aa in VHH_HALLMARKS_IMGT[target_hallmarks].items():
            idx = imgt_to_idx.get(imgt_pos)
            if idx is None or idx >= len(fr2):
                continue
            
            original = fr2[idx]
            if original == target_aa:
                continue
            
            mutations.append(Mutation(
                position=f"IMGT{imgt_pos}",
                imgt_num=imgt_pos,
                original=original,
                mutant=target_aa,
                source=target_hallmarks,
                confidence=0.95
            ))
        
        return mutations
    
    def _get_triplet_mutation(self, triplet_rules: List[dict], fr4: str) -> Optional[Mutation]:
        if not triplet_rules:
            return None
        
        best = triplet_rules[0]
        fr4_motif = best.get('suggested_aa', '')
        
        if len(fr4_motif) != 3 or len(fr4) < 3:
            return None
        
        if fr4[:3] == fr4_motif:
            return None
        
        return Mutation(
            position="IMGT118-120",
            imgt_num=118,
            original=fr4[:3],
            mutant=fr4_motif,
            source=f"triplet",
            confidence=best.get('confidence', 0.8)
        )
    
    def _imgt_to_fr(self, imgt_num: int) -> Tuple[Optional[str], int]:
        if 1 <= imgt_num <= 26:
            return 'FR1', imgt_num - 1
        elif 39 <= imgt_num <= 55:
            return 'FR2', imgt_num - 39
        elif 66 <= imgt_num <= 104:
            return 'FR3', imgt_num - 66
        elif 118 <= imgt_num <= 128:
            return 'FR4', imgt_num - 118
        return None, 0
    
    def _make_universal_candidate(self, name: str, cdrs: CDRSet, family: str,
                                    mutations: List[Mutation], strategy: str,
                                    gen_order: int, description: str, 
                                    target_family: str) -> VHHCandidate:
        fr1 = list(self.scaffold['FR1'])
        fr2 = list(self.scaffold['FR2'])
        fr3 = list(self.scaffold['FR3'])
        fr4 = list(self.scaffold['FR4'])
        
        self._apply_mutations(fr1, fr2, fr3, fr4, mutations)
        
        fr1_str, fr2_str = ''.join(fr1), ''.join(fr2)
        fr3_str, fr4_str = ''.join(fr3), ''.join(fr4)
        
        sequence = fr1_str + cdrs.cdr1 + fr2_str + cdrs.cdr2 + fr3_str + cdrs.cdr3 + fr4_str
        
        # Classify the resulting family
        temp_cdrs = CDRSet(cdrs.cdr1, cdrs.cdr2, cdrs.cdr3, fr1_str, fr2_str, fr3_str, fr4_str)
        detected_family = classify_family(temp_cdrs, sequence)
        
        return VHHCandidate(
            id=name, rank=0, sequence=sequence,
            framework_source='universal', family=detected_family,
            cdr1=cdrs.cdr1, cdr2=cdrs.cdr2, cdr3=cdrs.cdr3,
            fr1=fr1_str, fr2=fr2_str, fr3=fr3_str, fr4=fr4_str,
            mutations=mutations, strategy=strategy,
            generation_order=gen_order,
            construction_method=description,
            target_family=target_family,
            framework_identity_pct=0.0  # Will be calculated later vs original
        )
    
    def _make_original_candidate(self, name: str, cdrs: CDRSet, family: str,
                                   mutations: List[Mutation], strategy: str, 
                                   orig: dict, gen_order: int, description: str,
                                   target_family: str) -> VHHCandidate:
        fr1 = list(orig['FR1'])
        fr2 = list(orig['FR2'])
        fr3 = list(orig['FR3'])
        fr4 = list(orig['FR4'])
        
        self._apply_mutations(fr1, fr2, fr3, fr4, mutations)
        
        fr1_str, fr2_str = ''.join(fr1), ''.join(fr2)
        fr3_str, fr4_str = ''.join(fr3), ''.join(fr4)
        
        sequence = fr1_str + cdrs.cdr1 + fr2_str + cdrs.cdr2 + fr3_str + cdrs.cdr3 + fr4_str
        
        # Classify the resulting family
        temp_cdrs = CDRSet(cdrs.cdr1, cdrs.cdr2, cdrs.cdr3, fr1_str, fr2_str, fr3_str, fr4_str)
        detected_family = classify_family(temp_cdrs, sequence)
        
        return VHHCandidate(
            id=name, rank=0, sequence=sequence,
            framework_source='original', family=detected_family,
            cdr1=cdrs.cdr1, cdr2=cdrs.cdr2, cdr3=cdrs.cdr3,
            fr1=fr1_str, fr2=fr2_str, fr3=fr3_str, fr4=fr4_str,
            mutations=mutations, strategy=strategy,
            generation_order=gen_order,
            construction_method=description,
            target_family=target_family
        )
    
    def _apply_mutations(self, fr1: list, fr2: list, fr3: list, fr4: list, mutations: List[Mutation]):
        for mut in mutations:
            if mut.position == "IMGT118-120" and len(mut.mutant) == 3:
                for i, aa in enumerate(mut.mutant):
                    if i < len(fr4):
                        fr4[i] = aa
                continue
            
            region, idx = self._imgt_to_fr(mut.imgt_num)
            if region == 'FR1' and idx < len(fr1):
                fr1[idx] = mut.mutant
            elif region == 'FR2' and idx < len(fr2):
                fr2[idx] = mut.mutant
            elif region == 'FR3' and idx < len(fr3):
                fr3[idx] = mut.mutant
            elif region == 'FR4' and idx < len(fr4):
                fr4[idx] = mut.mutant

# ============================================================
# OUTPUT
# ============================================================

def to_dataframe(candidates: List[VHHCandidate]) -> 'pd.DataFrame':
    """Convert candidates to comprehensive DataFrame."""
    if not PANDAS_AVAILABLE:
        print("Warning: pandas not available")
        return None
    
    rows = []
    for c in candidates:
        # Get top family
        top_family = ""
        top_family_prob = 0.0
        if c.scoring.family_probabilities:
            top_family, top_family_prob = max(c.scoring.family_probabilities.items(), key=lambda x: x[1])
        
        row = {
            # Identification
            'rank': c.rank,
            'generation_order': c.generation_order,
            'id': c.id,
            'is_lead': c.is_lead,
            
            # Family info
            'detected_family': c.family,
            'target_family': c.target_family,
            'top_prob_family': top_family,
            'top_family_prob': round(top_family_prob, 3),
            
            # Construction info
            'framework_source': c.framework_source,
            'strategy': c.strategy,
            'construction_method': c.construction_method,
            'n_mutations': len(c.mutations),
            'mutations': ';'.join(str(m) for m in c.mutations),
            
            # Framework metrics
            'framework_identity_pct': round(c.framework_identity_pct, 1),
            'vernier_matches': c.scoring.vernier_matches,
            'vernier_total': c.scoring.vernier_total,
            
            # Rule compliance
            'rules_passed': c.scoring.rules_passed,
            'rules_applicable': c.scoring.rules_applicable,
            'rules_total': c.scoring.rules_total,
            'weighted_naturalness': round(c.scoring.weighted_naturalness, 4),
            
            # ESM2 scores
            'esm2_loss': round(c.scoring.esm2_loss, 4),
            'esm2_perplexity': round(c.scoring.esm2_perplexity, 2),
            
            # ESMFold pLDDT scores
            'plddt_mean': round(c.scoring.plddt_mean, 1),
            'plddt_median': round(c.scoring.plddt_median, 1),
            'plddt_min': round(c.scoring.plddt_min, 1),
            'plddt_cdr1': round(c.scoring.plddt_cdr1, 1),
            'plddt_cdr2': round(c.scoring.plddt_cdr2, 1),
            'plddt_cdr3': round(c.scoring.plddt_cdr3, 1),
            'plddt_framework': round(c.scoring.plddt_framework, 1),
            
            # Combined score
            'combined_score': round(c.scoring.combined_score, 4) if not c.is_lead else 999.0,
            
            # Violations
            'top_violations': '; '.join(c.scoring.top_violations[:3]),
            
            # Sequence regions
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
            'hallmarks': get_hallmarks_from_fr2(c.fr2),
        }
        rows.append(row)
    
    return pd.DataFrame(rows)

def save_results(df, output_dir: str, base_name: str, datetime_str: str):
    """Save results to files."""
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
            header = f">{row['id']}|{row['framework_source']}|{score_str}|plddt={row['plddt_mean']:.1f}"
            f.write(f"{header}\n{row['sequence']}\n")
    print(f"Saved: {fasta_path}")
    
    # Summary JSON
    summary = {
        'datetime': datetime_str,
        'n_candidates': len(df),
        'n_leads': int(df['is_lead'].sum()),
        'n_universal': int((df['framework_source'] == 'universal').sum()),
        'n_original': int((df['framework_source'] == 'original').sum()),
        'n_input': int((df['framework_source'] == 'input').sum()),
        'families_detected': df['detected_family'].value_counts().to_dict(),
        'score_stats': {
            'combined': {
                'min': float(df[~df['is_lead']]['combined_score'].min()) if len(df[~df['is_lead']]) > 0 else 0,
                'max': float(df[~df['is_lead']]['combined_score'].max()) if len(df[~df['is_lead']]) > 0 else 0,
                'mean': float(df[~df['is_lead']]['combined_score'].mean()) if len(df[~df['is_lead']]) > 0 else 0,
            },
            'plddt_mean': {
                'min': float(df['plddt_mean'].min()),
                'max': float(df['plddt_mean'].max()),
                'mean': float(df['plddt_mean'].mean()),
            },
            'rules_passed': {
                'min': int(df['rules_passed'].min()),
                'max': int(df['rules_passed'].max()),
                'mean': float(df['rules_passed'].mean()),
            }
        },
        'top_10': df.head(10)[['rank', 'id', 'combined_score', 'plddt_mean', 'rules_passed']].to_dict('records')
    }
    
    summary_path = os.path.join(output_dir, f"{base_name}_summary.json")
    with open(summary_path, 'w') as f:
        json.dump(summary, f, indent=2)
    print(f"Saved: {summary_path}")

# ============================================================
# MAIN
# ============================================================

def main():
    p = argparse.ArgumentParser(
        description="VHH Designer v7.4 - ESMFold Structure Prediction & Comprehensive Scoring",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Full pipeline with ESMFold
  python vhh_designer_v7_4_esmfold.py -i input.fasta \\
      --rules analysis_rules_v7.json \\
      --archetypes analysis_vernier_archetypes_v7.json \\
      --n-generate 500 \\
      --n-select 92

  # Without structure prediction (faster)
  python vhh_designer_v7_4_esmfold.py -i input.fasta \\
      --rules analysis_rules_v7.json \\
      --archetypes analysis_vernier_archetypes_v7.json \\
      --no-esmfold
"""
    )
    
    p.add_argument('--sequence', '-s', help='Input sequence')
    p.add_argument('--input', '-i', help='Input FASTA file')
    p.add_argument('--rules', '-r', required=True, help='analysis_rules_v7.json')
    p.add_argument('--archetypes', '-a', required=True, help='analysis_vernier_archetypes_v7.json')
    
    # Generation options
    p.add_argument('--n-generate', type=int, default=500, help='Candidates to generate (default: 500)')
    p.add_argument('--n-select', type=int, default=92, help='Candidates to select (default: 92)')
    p.add_argument('--target-hallmarks', '-t', default='FERG', help='Target hallmarks (default: FERG)')
    p.add_argument('--target-families', nargs='+', default=['F_C2', 'Y_C2', 'Other_VHH'])
    p.add_argument('--mode', '-m', choices=['both', 'universal', 'original'], default='both')
    
    # Scoring options
    p.add_argument('--use-esm2', action='store_true', default=True)
    p.add_argument('--no-esm2', action='store_true', help='Disable ESM2 scoring')
    p.add_argument('--use-esmfold', action='store_true', default=True)
    p.add_argument('--no-esmfold', action='store_true', help='Disable ESMFold (faster)')
    p.add_argument('--esm-model', default='facebook/esm2_t6_8M_UR50D')
    
    # Weights
    p.add_argument('--esm2-weight', type=float, default=0.2, help='ESM2 weight (default: 0.2)')
    p.add_argument('--esmfold-weight', type=float, default=0.3, help='ESMFold weight (default: 0.3)')
    p.add_argument('--rule-weight', type=float, default=0.5, help='Rule compliance weight (default: 0.5)')
    
    # Output
    p.add_argument('--output-dir', '-o', help='Output directory')
    p.add_argument('--name', help='Sequence name')
    
    args = p.parse_args()
    
    use_esm2 = args.use_esm2 and not args.no_esm2
    use_esmfold = args.use_esmfold and not args.no_esmfold
    
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
    print("VHH DESIGNER v7.4 - ESMFold Integration")
    print("=" * 70)
    print(f"Sequence: {seq_name}")
    print(f"Input ({len(input_seq)} aa): {input_seq[:50]}...")
    print(f"\nGeneration: {args.n_generate} → Select top {args.n_select}")
    print(f"Mode: {args.mode}")
    print(f"\nScoring:")
    print(f"  ESM2: {'enabled' if use_esm2 else 'disabled'}")
    print(f"  ESMFold: {'enabled' if use_esmfold else 'disabled'}")
    print(f"  Weights: ESM2={args.esm2_weight}, ESMFold={args.esmfold_weight}, Rules={args.rule_weight}")
    
    # Extract CDRs
    print(f"\n{'='*70}")
    print("CDR EXTRACTION")
    print(f"{'='*70}")
    
    cdrs = extract_cdrs(input_seq)
    if not cdrs:
        print("ERROR: CDR extraction failed")
        return
    
    print(f"CDR1: {cdrs.cdr1}")
    print(f"CDR2: {cdrs.cdr2}")
    print(f"CDR3: {cdrs.cdr3} ({len(cdrs.cdr3)} aa)")
    print(f"Hallmarks: {get_hallmarks_from_fr2(cdrs.fr2)}")
    
    original_frs = {'FR1': cdrs.fr1, 'FR2': cdrs.fr2, 'FR3': cdrs.fr3, 'FR4': cdrs.fr4}
    cdr3_features = get_cdr3_features(cdrs.cdr3)
    
    # Generate candidates
    print(f"\n{'='*70}")
    print("CANDIDATE GENERATION")
    print(f"{'='*70}")
    
    generator = CandidateGenerator(args.rules, args.archetypes)
    candidates = generator.generate(
        cdrs, args.n_generate,
        target_hallmarks=args.target_hallmarks,
        target_families=args.target_families,
        mode=args.mode,
        input_sequence=input_seq
    )
    
    print(f"\nGenerated {len(candidates)} candidates")
    
    # Score candidates
    print(f"\n{'='*70}")
    print("SCORING")
    print(f"{'='*70}")
    
    scorer = CombinedScorer(
        args.rules, args.archetypes,
        esm_model=args.esm_model,
        use_esm2=use_esm2,
        use_esmfold=use_esmfold,
        esm2_weight=args.esm2_weight,
        esmfold_weight=args.esmfold_weight,
        rule_weight=args.rule_weight
    )
    
    candidates = scorer.score_candidates(candidates, cdr3_features, original_frs)
    
    # Rank: Leads stay at top (score=999), others by combined score
    print(f"\n{'='*70}")
    print("RANKING & SELECTION")
    print(f"{'='*70}")
    
    candidates.sort(key=lambda c: (-c.scoring.combined_score,))
    
    # Update ranks (leads get rank 0)
    lead_rank = 0
    other_rank = 1
    for c in candidates:
        if c.is_lead:
            c.rank = lead_rank
            lead_rank += 1
        else:
            c.rank = other_rank
            other_rank += 1
    
    # Re-sort by rank
    candidates.sort(key=lambda c: c.rank)
    
    # Select top N
    selected = candidates[:args.n_select]
    
    print(f"\nSelected top {len(selected)} candidates")
    print(f"\nTop 15:")
    print("-" * 100)
    print(f"{'Rank':>4} {'Gen#':>4} {'ID':40} {'Score':>7} {'pLDDT':>6} {'Rules':>8} {'FW%':>5} {'Family':10}")
    print("-" * 100)
    
    for c in selected[:15]:
        score_str = "LEAD" if c.is_lead else f"{c.scoring.combined_score:.3f}"
        rules_str = f"{c.scoring.rules_passed}/{c.scoring.rules_applicable}"
        print(f"{c.rank:4d} {c.generation_order:4d} {c.id[:40]:40} {score_str:>7} {c.scoring.plddt_mean:6.1f} {rules_str:>8} {c.framework_identity_pct:5.1f} {c.family:10}")
    
    # Score distribution
    non_lead_scores = [c.scoring.combined_score for c in candidates if not c.is_lead]
    if non_lead_scores:
        print(f"\nScore distribution (non-leads):")
        print(f"  Min: {min(non_lead_scores):.3f}")
        print(f"  Max: {max(non_lead_scores):.3f}")
        print(f"  Mean: {np.mean(non_lead_scores):.3f}")
        print(f"  Top 10%: {np.percentile(non_lead_scores, 90):.3f}")
    
    # Save results
    df = to_dataframe(selected)
    
    datetime_str = datetime.now().strftime('%Y%m%d_%H%M%S')
    script_name = "vhh_designer_v7_4"
    folder_name = f"{datetime_str}_{script_name}_{seq_name}"
    
    if args.output_dir:
        output_dir = args.output_dir
    else:
        results_base = "results/designer_runs"
        output_dir = os.path.join(results_base, folder_name) if os.path.exists(results_base) else folder_name
    
    save_results(df, output_dir, f"{script_name}_{seq_name}", datetime_str)
    
    # Save all candidates
    df_all = to_dataframe(candidates)
    all_csv = os.path.join(output_dir, f"{script_name}_{seq_name}_all_{len(candidates)}.csv")
    df_all.to_csv(all_csv, index=False)
    print(f"Saved: {all_csv}")
    
    print(f"\n{'='*70}")
    print("DONE!")
    print(f"{'='*70}")

if __name__ == '__main__':
    main()
