#!/usr/bin/env python3
"""
VHH Designer v7.5 - Fixed AntPack Parsing & IMGT Position Mapping
==================================================================

FIXES from v7.4:
1. Correct AntPack 4-tuple parsing + trim_alignment() for insertion safety
2. IMGT position-based lookup (no more hardcoded FR2 indices)
3. Returns {imgt_pos: aa} dict for safe position access throughout
4. Proper handling of insertions (52A, 52B, etc.)
5. Region extraction aligned to IMGT numbering, not raw sequence indices

This version includes:
- Candidate generation (N candidates, default 500)
- ESM2 language model scoring (pseudo-perplexity)  
- ESMFold structure prediction (pLDDT confidence)
- Multi-family probabilistic rule compliance scoring
- Comprehensive CSV output with all metrics
- Lead candidates preserved and kept at top

Usage:
  python vhh_designer_v7_5_fixed.py -i input.fasta \\
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
    
    def get_aa(self, imgt_pos: Any) -> str:
        """Get AA at IMGT position."""
        return self.imgt_positions.get(imgt_pos, self.imgt_positions.get(str(imgt_pos), ''))

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
    
    Uses positions dict instead of hardcoded string indices.
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
    
    # Count cysteines in full sequence
    if full_sequence:
        cys_count = full_sequence.count('C')
    else:
        # Count from positions
        cys_count = sum(1 for aa in positions.values() if aa == 'C')
    
    c_label = "C4" if cys_count >= 4 else "C2"
    
    # F or Y family based on position 42
    if p42 == 'F':
        return f"F_{c_label}"
    elif p42 == 'Y':
        return f"Y_{c_label}"
    
    return "Other_VHH"

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
        self._load_attempted = False  # Track if we've already tried loading
        self._load_failed = False     # Track if loading failed
    
    def load(self):
        # If already loaded successfully, return True
        if self.model is not None:
            return True
        
        # If we already tried and failed, don't try again
        if self._load_failed:
            return False
        
        # If we haven't tried yet, mark that we're attempting
        if self._load_attempted:
            return False
        self._load_attempted = True
        
        if not ESMFOLD_AVAILABLE or not TORCH_AVAILABLE:
            print("  ESMFold not available (missing dependencies)")
            self._load_failed = True
            return False
        
        print("  Loading ESMFold model (this may take a moment)...")
        try:
            self.model = esm.pretrained.esmfold_v1()
            self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
            self.model = self.model.to(self.device)
            self.model.eval()
            self.model.set_chunk_size(128)
            print(f"    ESMFold loaded successfully. Device: {self.device}")
            return True
        except Exception as e:
            print(f"    Failed to load ESMFold: {e}")
            print(f"    Try: pip install omegaconf fair-esm")
            self._load_failed = True
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
        # Try to load once at the start - if it fails, return empty results for all
        if not self.load():
            print("    Skipping ESMFold predictions (model not available)")
            return [{'plddt_mean': 0.0, 'plddt_median': 0.0, 'plddt_min': 0.0, 'plddt_scores': [],
                     'plddt_cdr1': 0.0, 'plddt_cdr2': 0.0, 'plddt_cdr3': 0.0, 'plddt_framework': 0.0} 
                    for _ in sequences]
        
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
    
    def __init__(self, rules_file: str, archetypes_file: str):
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
                 esmfold_weight: float = 0.3,
                 rule_weight: float = 0.5):
        
        self.multi_family_scorer = MultiFamilyProbabilisticScorer(rules_file, archetypes_file)
        
        self.use_esm2 = use_esm2 and ESM2_AVAILABLE and TORCH_AVAILABLE
        self.use_esmfold = use_esmfold and ESMFOLD_AVAILABLE and TORCH_AVAILABLE
        
        # Normalize weights
        total_weight = rule_weight
        if self.use_esm2:
            total_weight += esm2_weight
        if self.use_esmfold:
            total_weight += esmfold_weight
        
        self.esm2_weight = esm2_weight / total_weight if total_weight > 0 else 0
        self.esmfold_weight = esmfold_weight / total_weight if total_weight > 0 else 0
        self.rule_weight = rule_weight / total_weight if total_weight > 0 else 1
        
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
                 target_families: List[str] = None,
                 mode: str = 'both') -> List[VHHCandidate]:
        """Generate N candidates with proper IMGT position tracking."""
        candidates = []
        generation_order = 0
        
        if not target_families:
            target_families = ['F_C2', 'Y_C2', 'Other_VHH']
        
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
        
        # Calculate candidates per family
        n_remaining = n_candidates - 1
        n_per_family = n_remaining // len(target_families)
        
        for i, family in enumerate(target_families):
            n_for_family = n_per_family
            if i == len(target_families) - 1:
                n_for_family = n_remaining - (len(candidates) - 1)
            
            if n_for_family <= 0:
                continue
            
            fam_hallmarks = 'YERL' if family.startswith('Y_') else target_hallmarks
            
            rules = self._get_rules(fam_hallmarks, family, cdr3_features)
            archetype = self._get_archetype(family)
            triplet_rules = self._get_triplet_rules(family, cdrs.cdr3)
            
            print(f"  {family}: {n_for_family} candidates, {len(rules)} rules")
            
            if mode in ['universal', 'both']:
                n_univ = n_for_family // 2 if mode == 'both' else n_for_family
                new_cands, generation_order = self._generate_universal(
                    cdrs, family, rules, archetype, triplet_rules, n_univ,
                    fam_hallmarks, generation_order
                )
                candidates.extend(new_cands)
            
            if mode in ['original', 'both'] and cdrs.fr1:
                n_orig = n_for_family - (n_for_family // 2 if mode == 'both' else 0)
                new_cands, generation_order = self._generate_original(
                    cdrs, family, rules, archetype, triplet_rules, n_orig,
                    fam_hallmarks, generation_order, original_positions
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
    
    def _generate_universal(self, cdrs: CDRSet, family: str, rules: List[dict],
                             archetype: Dict[int, str], triplet_rules: List[dict],
                             n: int, target_hallmarks: str,
                             generation_order: int) -> Tuple[List[VHHCandidate], int]:
        """Generate using universal scaffold with position-based mutations."""
        candidates = []
        used = set()
        
        hallmark_set = set(HALLMARK_POSITIONS)
        
        # Get mutations using IMGT positions
        vernier_muts = []
        other_arch_muts = []
        
        for imgt_num, consensus_aa in archetype.items():
            if imgt_num in hallmark_set:
                continue
            
            original_aa = self.scaffold_positions.get(imgt_num, '')
            if not original_aa or original_aa == consensus_aa:
                continue
            
            mut = Mutation(
                position=f"IMGT{imgt_num}",
                imgt_num=imgt_num,
                original=original_aa,
                mutant=consensus_aa,
                source="archetype",
                confidence=0.7
            )
            
            if imgt_num in VERNIER_POSITIONS_FR3:
                vernier_muts.append(mut)
            else:
                other_arch_muts.append(mut)
        
        # Rule mutations
        rule_muts = []
        for rule in rules:
            imgt_num = rule.get('position_num', 0)
            suggested_aa = rule.get('suggested_aa', '')
            
            if not imgt_num or not suggested_aa or imgt_num in hallmark_set:
                continue
            
            original_aa = self.scaffold_positions.get(imgt_num, '')
            if not original_aa or original_aa == suggested_aa:
                continue
            
            rule_muts.append(Mutation(
                position=f"IMGT{imgt_num}",
                imgt_num=imgt_num,
                original=original_aa,
                mutant=suggested_aa,
                source=rule.get('condition', '')[:20],
                confidence=rule.get('confidence', 0.5)
            ))
        
        # Triplet mutation
        triplet_mut = None
        if triplet_rules:
            best = triplet_rules[0]
            fr4_motif = best.get('suggested_aa', '')
            if len(fr4_motif) == 3:
                orig_fr4 = UNIVERSAL_SCAFFOLD['FR4'][:3]
                if orig_fr4 != fr4_motif:
                    triplet_mut = Mutation(
                        position="IMGT118-120",
                        imgt_num="118-120",
                        original=orig_fr4,
                        mutant=fr4_motif,
                        source="triplet",
                        confidence=best.get('confidence', 0.8)
                    )
        
        # Strategy combinations
        strategies = [
            ("graft", [], "Universal scaffold + CDR graft only"),
            ("vernier", vernier_muts, f"Universal + {len(vernier_muts)} vernier mutations"),
            ("full_arch", vernier_muts + other_arch_muts, f"Universal + full {family} archetype"),
            ("vernier_rules", vernier_muts + rule_muts, f"Universal + vernier + rules"),
            ("full_optimized", vernier_muts + other_arch_muts + rule_muts, "Universal + full optimization"),
        ]
        
        if triplet_mut:
            strategies.append(("vernier_triplet", vernier_muts + [triplet_mut], "Universal + vernier + FR4 junction"))
        
        for strategy_name, mutations, description in strategies:
            if len(candidates) >= n:
                break
            
            key = tuple(sorted([str(m) for m in mutations]))
            if key in used:
                continue
            used.add(key)
            
            generation_order += 1
            candidate = self._build_universal_candidate(
                f"Univ_{family}_{strategy_name}_{generation_order}",
                cdrs, family, mutations, strategy_name, generation_order, description, family
            )
            candidates.append(candidate)
        
        # Random combinations
        import random
        all_muts = vernier_muts + other_arch_muts + rule_muts
        if triplet_mut:
            all_muts.append(triplet_mut)
        
        attempts = 0
        while len(candidates) < n and len(all_muts) >= 2 and attempts < n * 3:
            attempts += 1
            n_pick = random.randint(2, min(5, len(all_muts)))
            picked = random.sample(all_muts, n_pick)
            
            # Check no duplicate positions
            pos_set = set()
            valid = True
            for m in picked:
                if m.imgt_num in pos_set:
                    valid = False
                    break
                pos_set.add(m.imgt_num)
            
            if valid:
                key = tuple(sorted([str(m) for m in picked]))
                if key not in used:
                    used.add(key)
                    generation_order += 1
                    candidate = self._build_universal_candidate(
                        f"Univ_{family}_combo_{generation_order}",
                        cdrs, family, picked, "random_combo", generation_order,
                        f"Universal + random combo of {len(picked)} mutations", family
                    )
                    candidates.append(candidate)
        
        return candidates[:n], generation_order
    
    def _generate_original(self, cdrs: CDRSet, family: str, rules: List[dict],
                            archetype: Dict[int, str], triplet_rules: List[dict],
                            n: int, target_hallmarks: str, generation_order: int,
                            original_positions: Dict[Any, str]) -> Tuple[List[VHHCandidate], int]:
        """Generate using original FRs with position-based mutations."""
        candidates = []
        used = set()
        
        # Hallmark mutations (required for VHH conversion)
        hallmark_muts = []
        if target_hallmarks in VHH_HALLMARKS_IMGT:
            for imgt_pos, target_aa in VHH_HALLMARKS_IMGT[target_hallmarks].items():
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
        
        # Archetype mutations
        vernier_muts = []
        other_arch_muts = []
        
        for imgt_num, consensus_aa in archetype.items():
            if imgt_num in HALLMARK_POSITIONS:
                continue
            
            original_aa = original_positions.get(imgt_num, '')
            if not original_aa or original_aa == consensus_aa:
                continue
            
            mut = Mutation(
                position=f"IMGT{imgt_num}",
                imgt_num=imgt_num,
                original=original_aa,
                mutant=consensus_aa,
                source="archetype",
                confidence=0.7
            )
            
            if imgt_num in VERNIER_POSITIONS_FR3:
                vernier_muts.append(mut)
            else:
                other_arch_muts.append(mut)
        
        # Rule mutations
        rule_muts = []
        for rule in rules:
            imgt_num = rule.get('position_num', 0)
            suggested_aa = rule.get('suggested_aa', '')
            
            if not imgt_num or not suggested_aa or imgt_num in HALLMARK_POSITIONS:
                continue
            
            original_aa = original_positions.get(imgt_num, '')
            if not original_aa or original_aa == suggested_aa:
                continue
            
            rule_muts.append(Mutation(
                position=f"IMGT{imgt_num}",
                imgt_num=imgt_num,
                original=original_aa,
                mutant=suggested_aa,
                source=rule.get('condition', '')[:20],
                confidence=rule.get('confidence', 0.5)
            ))
        
        # Triplet mutation
        triplet_mut = None
        if triplet_rules:
            best = triplet_rules[0]
            fr4_motif = best.get('suggested_aa', '')
            if len(fr4_motif) == 3 and len(cdrs.fr4) >= 3:
                if cdrs.fr4[:3] != fr4_motif:
                    triplet_mut = Mutation(
                        position="IMGT118-120",
                        imgt_num="118-120",
                        original=cdrs.fr4[:3],
                        mutant=fr4_motif,
                        source="triplet",
                        confidence=best.get('confidence', 0.8)
                    )
        
        # Strategies (all include hallmarks)
        n_hallmark = len(hallmark_muts)
        strategies = [
            ("hallmark_only", hallmark_muts, f"Original FRs + {n_hallmark} hallmark mutations ({target_hallmarks})"),
            ("vernier", hallmark_muts + vernier_muts, f"Original + hallmarks + vernier"),
            ("full_arch", hallmark_muts + vernier_muts + other_arch_muts, f"Original + hallmarks + full archetype"),
            ("vernier_rules", hallmark_muts + vernier_muts + rule_muts, "Original + hallmarks + vernier + rules"),
            ("full_optimized", hallmark_muts + vernier_muts + other_arch_muts + rule_muts, "Original + full optimization"),
        ]
        
        if triplet_mut:
            strategies.append(("vernier_triplet", hallmark_muts + vernier_muts + [triplet_mut], "Original + vernier + FR4 junction"))
        
        for strategy_name, mutations, description in strategies:
            if len(candidates) >= n:
                break
            
            key = tuple(sorted([str(m) for m in mutations]))
            if key in used:
                continue
            used.add(key)
            
            generation_order += 1
            candidate = self._build_original_candidate(
                f"Orig_{family}_{strategy_name}_{generation_order}",
                cdrs, family, mutations, strategy_name, generation_order,
                description, family, original_positions
            )
            candidates.append(candidate)
        
        # Random combinations
        import random
        all_optional = vernier_muts + other_arch_muts + rule_muts
        if triplet_mut:
            all_optional.append(triplet_mut)
        
        attempts = 0
        while len(candidates) < n and len(all_optional) >= 2 and attempts < n * 3:
            attempts += 1
            n_pick = random.randint(2, min(5, len(all_optional)))
            picked = random.sample(all_optional, n_pick)
            
            pos_set = set()
            valid = True
            for m in picked:
                if m.imgt_num in pos_set:
                    valid = False
                    break
                pos_set.add(m.imgt_num)
            
            if valid:
                muts = hallmark_muts + picked
                key = tuple(sorted([str(m) for m in muts]))
                if key not in used:
                    used.add(key)
                    generation_order += 1
                    candidate = self._build_original_candidate(
                        f"Orig_{family}_combo_{generation_order}",
                        cdrs, family, muts, "random_combo", generation_order,
                        f"Original + hallmarks + random combo of {len(picked)} mutations",
                        family, original_positions
                    )
                    candidates.append(candidate)
        
        return candidates[:n], generation_order
    
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
        
        return VHHCandidate(
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
        
        return VHHCandidate(
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

# ============================================================
# OUTPUT
# ============================================================

def to_dataframe(candidates: List[VHHCandidate]) -> 'pd.DataFrame':
    """Convert candidates to comprehensive DataFrame."""
    if not PANDAS_AVAILABLE:
        return None
    
    rows = []
    for c in candidates:
        top_family = ""
        top_family_prob = 0.0
        if c.scoring.family_probabilities:
            top_family, top_family_prob = max(c.scoring.family_probabilities.items(), key=lambda x: x[1])
        
        # Get hallmarks from positions
        hallmarks = get_hallmarks_from_positions(c.imgt_positions)
        
        row = {
            'rank': c.rank,
            'generation_order': c.generation_order,
            'id': c.id,
            'is_lead': c.is_lead,
            
            'detected_family': c.family,
            'target_family': c.target_family,
            'top_prob_family': top_family,
            'top_family_prob': round(top_family_prob, 3),
            'hallmarks': hallmarks,
            
            'framework_source': c.framework_source,
            'strategy': c.strategy,
            'construction_method': c.construction_method,
            'n_mutations': len(c.mutations),
            'mutations': ';'.join(str(m) for m in c.mutations),
            
            'framework_identity_pct': round(c.framework_identity_pct, 1),
            'vernier_matches': c.scoring.vernier_matches,
            'vernier_total': c.scoring.vernier_total,
            
            'rules_passed': c.scoring.rules_passed,
            'rules_applicable': c.scoring.rules_applicable,
            'rules_total': c.scoring.rules_total,
            'weighted_naturalness': round(c.scoring.weighted_naturalness, 4),
            
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

def save_results(df, output_dir: str, base_name: str, datetime_str: str):
    """Save results to files."""
    os.makedirs(output_dir, exist_ok=True)
    
    csv_path = os.path.join(output_dir, f"{base_name}.csv")
    df.to_csv(csv_path, index=False)
    print(f"\nSaved: {csv_path}")
    
    fasta_path = os.path.join(output_dir, f"{base_name}.fasta")
    with open(fasta_path, 'w') as f:
        for _, row in df.iterrows():
            score_str = "LEAD" if row['is_lead'] else f"score={row['combined_score']:.3f}"
            header = f">{row['id']}|{row['framework_source']}|{score_str}|plddt={row['plddt_mean']:.1f}"
            f.write(f"{header}\n{row['sequence']}\n")
    print(f"Saved: {fasta_path}")
    
    summary = {
        'datetime': datetime_str,
        'version': '7.5',
        'n_candidates': len(df),
        'n_leads': int(df['is_lead'].sum()),
        'n_universal': int((df['framework_source'] == 'universal').sum()),
        'n_original': int((df['framework_source'] == 'original').sum()),
        'families_detected': df['detected_family'].value_counts().to_dict(),
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
    
    p.add_argument('--n-generate', type=int, default=500)
    p.add_argument('--n-select', type=int, default=92)
    p.add_argument('--target-hallmarks', '-t', default='FERG')
    p.add_argument('--target-families', nargs='+', default=['F_C2', 'Y_C2', 'Other_VHH'])
    p.add_argument('--mode', '-m', choices=['both', 'universal', 'original'], default='both')
    
    p.add_argument('--use-esm2', action='store_true', default=True)
    p.add_argument('--no-esm2', action='store_true')
    p.add_argument('--use-esmfold', action='store_true', default=True)
    p.add_argument('--no-esmfold', action='store_true')
    p.add_argument('--esm-model', default='facebook/esm2_t6_8M_UR50D')
    
    p.add_argument('--esm2-weight', type=float, default=0.2)
    p.add_argument('--esmfold-weight', type=float, default=0.3)
    p.add_argument('--rule-weight', type=float, default=0.5)
    
    p.add_argument('--output-dir', '-o')
    p.add_argument('--name')
    
    args = p.parse_args()
    
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
    print("VHH DESIGNER v7.5 - Fixed IMGT Position Mapping")
    print("=" * 70)
    print(f"Sequence: {seq_name}")
    print(f"Input ({len(input_seq)} aa): {input_seq[:50]}...")
    print(f"\nGeneration: {args.n_generate} → Select top {args.n_select}")
    print(f"Mode: {args.mode}")
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
    
    generator = CandidateGenerator(args.rules, args.archetypes)
    candidates = generator.generate(
        cdrs, args.n_generate,
        target_hallmarks=args.target_hallmarks,
        target_families=args.target_families,
        mode=args.mode
    )
    
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
    
    # Rank
    print(f"\n{'='*70}")
    print("RANKING & SELECTION")
    print(f"{'='*70}")
    
    candidates.sort(key=lambda c: (-c.scoring.combined_score,))
    
    lead_rank = 0
    other_rank = 1
    for c in candidates:
        if c.is_lead:
            c.rank = lead_rank
            lead_rank += 1
        else:
            c.rank = other_rank
            other_rank += 1
    
    candidates.sort(key=lambda c: c.rank)
    selected = candidates[:args.n_select]
    
    print(f"\nSelected top {len(selected)} candidates")
    print(f"\nTop 15:")
    print("-" * 110)
    print(f"{'Rank':>4} {'Gen#':>4} {'ID':40} {'Score':>7} {'pLDDT':>6} {'Rules':>8} {'Vernier':>8} {'FW%':>5} {'Hallmarks':>9}")
    print("-" * 110)
    
    for c in selected[:15]:
        score_str = "LEAD" if c.is_lead else f"{c.scoring.combined_score:.3f}"
        rules_str = f"{c.scoring.rules_passed}/{c.scoring.rules_applicable}"
        vernier_str = f"{c.scoring.vernier_matches}/{c.scoring.vernier_total}"
        hallmarks = get_hallmarks_from_positions(c.imgt_positions)
        print(f"{c.rank:4d} {c.generation_order:4d} {c.id[:40]:40} {score_str:>7} {c.scoring.plddt_mean:6.1f} {rules_str:>8} {vernier_str:>8} {c.framework_identity_pct:5.1f} {hallmarks:>9}")
    
    # Save
    df = to_dataframe(selected)
    
    datetime_str = datetime.now().strftime('%Y%m%d_%H%M%S')
    script_name = "vhh_designer_v7_5"
    folder_name = f"{datetime_str}_{script_name}_{seq_name}"
    
    if args.output_dir:
        output_dir = args.output_dir
    else:
        results_base = "results/designer_runs"
        output_dir = os.path.join(results_base, folder_name) if os.path.exists(results_base) else folder_name
    
    save_results(df, output_dir, f"{script_name}_{seq_name}", datetime_str)
    
    df_all = to_dataframe(candidates)
    all_csv = os.path.join(output_dir, f"{script_name}_{seq_name}_all_{len(candidates)}.csv")
    df_all.to_csv(all_csv, index=False)
    print(f"Saved: {all_csv}")
    
    print(f"\n{'='*70}")
    print("DONE!")
    print(f"{'='*70}")

if __name__ == '__main__':
    main()
