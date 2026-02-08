#!/usr/bin/env python3
"""
VHH Analysis / Compensation Unified v7 (CSV + IMGT)
====================================================

Single-script engine that:
  - Streams an IMGT-annotated VHH CSV
  - Applies truncation and position exclusions
  - Infers family and hallmarks (if missing)
  - Extracts CDR3 features and conditions
  - Learns CDR→FR position rules (support / confidence / lift)
  - Identifies Vernier archetypes per family
  - Computes mutual information (MI) between positions on a subsample
  - Computes family-aware conservation + correlations (CorrelationAnalyzer)

OUTPUTS:
  - analysis_summary_v7.json        : Metadata, family counts, stats
  - analysis_rules_v7.json          : All compensation rules (unified schema)
  - analysis_vernier_archetypes_v7.json : Per-family vernier consensus
  - analysis_mi_pairs_v7.json       : Top MI position pairs
  - analysis_correlations_v7.json   : Per-family conservation + CDR3↔FR correlations

USAGE:
  python vhh_analysis_unified_v7.py \\
      --csv data/.../vhh_annotated_imgt.csv \\
      -o models/analysis/v7 \\
      --positions vernier \\
      --conditions minimal \\
      --filter-truncated \\
      --min-support 5000 \\
      --min-confidence 0.7 \\
      --min-lift 1.15 \\
      --mi-max-seqs 2000000
"""

import os
import csv
import re
import json
import math
import pickle
import argparse
from collections import defaultdict, Counter
from datetime import datetime
from typing import Dict, List, Optional, Tuple, Set

try:
    import numpy as np
    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False

# ============================================================
# CONSTANTS / IMGT / FAMILY LOGIC
# ============================================================

IMGT_REGIONS = {
    "FR1": (1, 26),
    "CDR1": (27, 38),
    "FR2": (39, 55),
    "CDR2": (56, 65),
    "FR3": (66, 104),
    "CDR3": (105, 117),
    "FR4": (118, 128),
}

HALLMARK_POSITIONS = [42, 49, 50, 52]
VERNIER_POSITIONS = [42, 49, 50, 52, 66, 67, 68, 69, 71, 76, 78, 82, 87, 89, 91, 94]

# Positions to EXCLUDE from rule learning
TRUNCATION_POSITIONS = {1, 2, 3, 126, 127, 128}  # N/C-terminal artifacts
IMGT_GAP_POSITIONS = {10, 73}  # Standard VHH gaps
EXCLUDE_FROM_RULES: Set[int] = TRUNCATION_POSITIONS | IMGT_GAP_POSITIONS

# Truncation filter thresholds
FR1_CORE = list(range(4, 27))   # 23 positions
FR4_CORE = list(range(118, 126))  # 8 positions
FR1_MIN_FILLED = 20
FR4_MIN_FILLED = 6

VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")
AA_POS = set("KRH")
AA_NEG = set("DE")
AA_AROMATIC = set("FWY")

HYDROPHOBICITY = {
    "I": 4.5, "V": 4.2, "L": 3.8, "F": 2.8, "C": 2.5, "M": 1.9, "A": 1.8,
    "G": -0.4, "T": -0.7, "S": -0.8, "W": -0.9, "Y": -1.3, "P": -1.6,
    "H": -3.2, "E": -3.5, "Q": -3.5, "D": -3.5, "N": -3.5, "K": -3.9, "R": -4.5,
}
CHARGE = {aa: (1 if aa in AA_POS else (-1 if aa in AA_NEG else 0)) for aa in VALID_AA}

# ============================================================
# UTILITY FUNCTIONS
# ============================================================

def norm_col(c: str) -> str:
    return re.sub(r"[^0-9a-zA-Z_]+", "_", str(c).strip().lower())

def norm_cell(x: Optional[str]) -> str:
    if x is None:
        return "-"
    x = str(x).strip()
    if x == "" or x.lower() in {"na", "nan", "none"}:
        return "-"
    return x

def passes_truncation_filter(imgt: Dict[int, str]) -> bool:
    """Check FR1/FR4 core completeness."""
    fr1_filled = sum(1 for p in FR1_CORE if imgt.get(p, "-") not in ("-", "X", "*", ""))
    fr4_filled = sum(1 for p in FR4_CORE if imgt.get(p, "-") not in ("-", "X", "*", ""))
    return fr1_filled >= FR1_MIN_FILLED and fr4_filled >= FR4_MIN_FILLED

def infer_hallmarks(imgt: Dict[int, str]) -> str:
    return "".join(imgt.get(p, "-") for p in HALLMARK_POSITIONS)

def infer_family(imgt: Dict[int, str], aa_v_full: str) -> str:
    """Classify VHH family based on hallmarks and cysteine count."""
    p42, p50, p52 = imgt.get(42, "-"), imgt.get(50, "-"), imgt.get(52, "-")
    
    if p50 == "L":
        return "VH_like"
    if p52 == "W":
        return "VHH_W52"
    
    cys_total = (aa_v_full or "").count("C")
    c_label = "C4" if cys_total >= 4 else "C2"
    
    if p42 in {"F", "Y"}:
        return f"{p42}_{c_label}"
    return "Other_VHH"

# ============================================================
# CDR FEATURES / BINNING
# ============================================================

def get_cdr_features(cdr: str) -> Dict[str, float]:
    if not cdr:
        return {"length": 0, "charge": 0, "hydrophobicity": 0.0, "n_cys": 0, "aromatic_count": 0}
    length = len(cdr)
    return {
        "length": length,
        "charge": sum(CHARGE.get(aa, 0) for aa in cdr),
        "hydrophobicity": round(sum(HYDROPHOBICITY.get(aa, 0) for aa in cdr) / max(length, 1), 3),
        "n_cys": cdr.count("C"),
        "aromatic_count": sum(1 for aa in cdr if aa in AA_AROMATIC),
    }

def categorize_cdr3(cdr3: str) -> Dict:
    f = get_cdr_features(cdr3)
    L, q, n = f["length"], f["charge"], f["n_cys"]
    
    len_cat = "short" if L <= 8 else ("medium" if L <= 14 else "long")
    charge_cat = "negative" if q <= -2 else ("positive" if q >= 2 else "neutral")
    cys_cat = "no_cys" if n == 0 else ("one_cys" if n == 1 else "multi_cys")
    
    return {
        "cdr3_len": len_cat,
        "cdr3_charge": charge_cat,
        "cdr3_cys": cys_cat,
        "length": f["length"],
        "charge": f["charge"],
        "hydrophobicity": f["hydrophobicity"],
        "n_cys": f["n_cys"],
        "aromatic_count": f["aromatic_count"],
    }

# ============================================================
# POSITION SELECTION
# ============================================================

def choose_positions(mode: str, available: Set[int], exclude_truncation: bool = True) -> List[int]:
    if mode == "vernier":
        ps = set(VERNIER_POSITIONS)
    elif mode == "fr3fr4":
        ps = set(range(66, 105)) | set(range(118, 129)) | set(HALLMARK_POSITIONS)
    elif mode == "all":
        ps = set()
        for region, (s, e) in IMGT_REGIONS.items():
            if region.startswith("FR"):
                ps.update(range(s, e + 1))
    else:
        raise ValueError(f"Unknown mode: {mode}")
    
    ps &= available
    if exclude_truncation:
        ps -= EXCLUDE_FROM_RULES
    return sorted(ps)

# ============================================================
# CONDITION BUILDING
# ============================================================

def build_conditions(cdr3_info: Dict, cdr3_seq: str, family: str, hallmarks: str, mode: str) -> List[str]:
    """Build list of condition tokens."""
    conds = [
        f"family={family}",
        f"hallmarks={hallmarks}",
        f"cdr3_len={cdr3_info['cdr3_len']}",
        f"cdr3_charge={cdr3_info['cdr3_charge']}",
        f"cdr3_cys={cdr3_info['cdr3_cys']}",
    ]
    
    if cdr3_seq:
        if len(cdr3_seq) >= 1:
            conds.append(f"cdr3[-1]={cdr3_seq[-1]}")
        if len(cdr3_seq) >= 2 and mode == "full":
            conds.append(f"cdr3[-2]={cdr3_seq[-2]}")
        if len(cdr3_seq) >= 3 and mode == "full":
            conds.append(f"cdr3[-3]={cdr3_seq[-3]}")
    
    return conds

# ============================================================
# COMPENSATION ANALYZER (Rule Mining)
# ============================================================

class CompensationAnalyzer:
    """Memory-efficient streaming compensation analyzer with family stratification."""
    
    def __init__(self, positions: List[int]):
        self.positions = positions
        # (family, condition, position) -> Counter(aa)
        self.counts: Dict[Tuple[str, str, int], Counter] = defaultdict(Counter)
        # (family, condition) -> total count (for n_total_for_condition)
        self.cond_totals: Dict[Tuple[str, str], int] = defaultdict(int)
        # (family, position) -> Counter(aa) for family-conditional baseline
        self.baseline_by_family: Dict[Tuple[str, int], Counter] = defaultdict(Counter)
        # position -> Counter(aa) for global baseline (fallback)
        self.baseline_global: Dict[int, Counter] = defaultdict(Counter)
        # Statistics
        self.family_counts: Counter = Counter()
        self.hallmark_counts: Counter = Counter()
        self.total_used = 0
    
    def update(self, imgt: Dict[int, str], cdr3_seq: str, family: str, 
               hallmarks: str, conditions_mode: str):
        cdr3_info = categorize_cdr3(cdr3_seq)
        conds = build_conditions(cdr3_info, cdr3_seq, family, hallmarks, conditions_mode)
        
        self.family_counts[family] += 1
        self.hallmark_counts[hallmarks] += 1
        
        # Update condition totals
        for cond in conds:
            self.cond_totals[(family, cond)] += 1
        
        for pos in self.positions:
            aa = imgt.get(pos, "-")
            if aa not in VALID_AA:
                continue
            # Track both family-specific and global baselines
            self.baseline_by_family[(family, pos)][aa] += 1
            self.baseline_global[pos][aa] += 1
            for cond in conds:
                self.counts[(family, cond, pos)][aa] += 1
        
        self.total_used += 1
    
    def extract_rules(self, min_condition_obs: int, min_aa_support: int,
                      min_confidence: float, min_lift: float, 
                      filter_definitional: bool = True) -> List[Dict]:
        """Extract rules with enhanced schema including condition_tokens and n_total_for_condition.
        
        Args:
            min_condition_obs: Minimum observations for condition+position to be considered
            min_aa_support: Minimum support for the specific AA (can be lower than min_condition_obs)
            min_confidence: Minimum P(AA | condition)
            min_lift: Minimum lift vs family baseline
            filter_definitional: Remove rules where hallmark positions define the condition
        """
        rules = []
        hallmark_nums = set(HALLMARK_POSITIONS)
        
        # Precompute family-conditional baseline probabilities
        baseline_prob_by_family: Dict[Tuple[str, int], Dict[str, float]] = {}
        for (family, pos), aa_counts in self.baseline_by_family.items():
            total = sum(aa_counts.values())
            if total > 0:
                baseline_prob_by_family[(family, pos)] = {aa: c / total for aa, c in aa_counts.items()}
        
        for (family, cond, pos), aa_counts in self.counts.items():
            cond_total = sum(aa_counts.values())
            # Gate on condition observations (not AA support)
            if cond_total < min_condition_obs:
                continue
            
            # Get family-specific baseline for this position
            base_probs = baseline_prob_by_family.get((family, pos), {})
            if not base_probs:
                continue
            
            # Filter definitional rules
            if filter_definitional:
                if (cond.startswith("hallmarks=") or cond.startswith("family=")) and pos in hallmark_nums:
                    continue
            
            # Get n_total_for_condition (across all positions)
            n_total_for_cond = self.cond_totals.get((family, cond), cond_total)
            
            # Get top 3 AAs for this condition
            sorted_aas = sorted(aa_counts.items(), key=lambda kv: -kv[1])
            
            for rank, (aa, count) in enumerate(sorted_aas[:3]):
                # Gate on AA support (separate threshold, typically lower)
                if count < min_aa_support:
                    continue
                
                conf = count / cond_total
                base_p = base_probs.get(aa, 0.0)
                if base_p <= 0:
                    continue
                
                # Lift is now family-conditional (apples-to-apples comparison)
                lift = conf / base_p
                
                if conf < min_confidence or lift < min_lift:
                    continue
                
                # confidence = P(AA | condition, valid position)
                # confidence_all = P(AA | condition) - accounts for gaps/invalid at this position
                conf_all = count / n_total_for_cond if n_total_for_cond > 0 else 0.0
                
                # Build runner-ups list (other top AAs)
                alt_aas = []
                for alt_rank, (alt_aa, alt_count) in enumerate(sorted_aas[:3]):
                    if alt_aa == aa:
                        continue
                    alt_conf = alt_count / cond_total
                    alt_base = base_probs.get(alt_aa, 0.0)
                    if alt_base > 0 and alt_count >= min_aa_support // 2:
                        alt_aas.append({
                            "aa": alt_aa,
                            "confidence": round(alt_conf, 4),
                            "support": int(alt_count),
                            "lift": round(alt_conf / alt_base, 4) if alt_base > 0 else 0.0,
                        })
                
                rules.append({
                    "condition": cond,
                    "condition_tokens": [cond],
                    "position": f"IMGT{pos}",
                    "position_num": pos,
                    "suggested_aa": aa,
                    "confidence": round(conf, 4),
                    "confidence_all": round(conf_all, 4),
                    "support": int(count),
                    "lift": round(lift, 4),
                    "baseline_prob": round(base_p, 4),
                    "n_total_for_condition": int(n_total_for_cond),
                    "alt_aas": alt_aas,  # Runner-ups
                    "rank": rank + 1,  # 1 = top, 2 = runner-up, 3 = second runner-up
                    "source": "vhh_analysis_v7",
                    "rule_type": "cdr_to_imgt",
                    "family": family,
                })
        
        rules.sort(key=lambda r: (-r["support"], -r["confidence"], -r["lift"]))
        return rules
    
    def vernier_archetypes(self, min_support: int = 10000) -> List[Dict]:
        """Extract vernier consensus per family."""
        fam_pos_counts: Dict[str, Dict[int, Counter]] = defaultdict(lambda: defaultdict(Counter))
        
        for (family, cond, pos), aa_counts in self.counts.items():
            if pos not in VERNIER_POSITIONS:
                continue
            if not cond.startswith("family="):
                continue
            for aa, c in aa_counts.items():
                fam_pos_counts[family][pos][aa] += c
        
        archetypes = []
        for family, pos_map in fam_pos_counts.items():
            if self.family_counts[family] < min_support:
                continue
            
            arch = {"family": family, "count": self.family_counts[family], "positions": {}}
            for pos in sorted(VERNIER_POSITIONS):
                if pos in pos_map and pos_map[pos]:
                    aa_counts = pos_map[pos]
                    total = sum(aa_counts.values())
                    aa, c = aa_counts.most_common(1)[0]
                    arch["positions"][f"IMGT{pos}"] = {
                        "consensus": aa,
                        "frequency": round(c / total, 4),
                        "count": c,
                    }
            archetypes.append(arch)
        
        archetypes.sort(key=lambda a: -a["count"])
        return archetypes

# ============================================================
# CORRELATION ANALYZER (Family-aware conservation + CDR3↔FR correlations)
# ============================================================

class CorrAccumulator:
    """Streaming Pearson correlation accumulator."""
    def __init__(self):
        self.n = 0
        self.sum_x = 0.0
        self.sum_y = 0.0
        self.sum_x2 = 0.0
        self.sum_y2 = 0.0
        self.sum_xy = 0.0
    
    def add(self, x: float, y: float):
        self.n += 1
        self.sum_x += x
        self.sum_y += y
        self.sum_x2 += x * x
        self.sum_y2 += y * y
        self.sum_xy += x * y
    
    def corr(self) -> Optional[float]:
        if self.n <= 1:
            return None
        n = float(self.n)
        num = self.sum_xy - (self.sum_x * self.sum_y / n)
        den_x = self.sum_x2 - (self.sum_x * self.sum_x / n)
        den_y = self.sum_y2 - (self.sum_y * self.sum_y / n)
        if den_x <= 0 or den_y <= 0:
            return None
        r = num / math.sqrt(den_x * den_y)
        return float(max(min(r, 1.0), -1.0))


class WelfordAccumulator:
    """Streaming mean/variance using Welford's algorithm."""
    def __init__(self):
        self.n = 0
        self.mean = 0.0
        self.M2 = 0.0
    
    def add(self, x: float):
        self.n += 1
        delta = x - self.mean
        self.mean += delta / self.n
        delta2 = x - self.mean
        self.M2 += delta * delta2
    
    def stats(self) -> Dict:
        if self.n < 2:
            return {"n": self.n, "mean": self.mean, "std": 0.0}
        variance = self.M2 / (self.n - 1)
        return {"n": self.n, "mean": round(self.mean, 4), "std": round(math.sqrt(variance), 4)}


class CorrelationAnalyzer:
    """
    Family-aware conservation & correlations between CDR3 features and FR properties.
    
    Tracks per position:
    - Conservation: major_aa, major_freq, entropy
    - Correlations: cdr3_len/charge/hydro vs is_hydrophobic/positive/negative at FR positions
    
    Also tracks CDR3 stats per family for legacy pickle output.
    """
    
    def __init__(self, positions: List[int]):
        self.positions = positions
        # family -> pos -> Counter(aa)
        self.family_pos_aa_counts: Dict[str, Dict[int, Counter]] = defaultdict(lambda: defaultdict(Counter))
        # family -> total sequences
        self.family_totals: Counter = Counter()
        # family -> pos -> corr_name -> CorrAccumulator
        self.family_pos_corr: Dict[str, Dict[int, Dict[str, CorrAccumulator]]] = defaultdict(
            lambda: defaultdict(dict)
        )
        # CDR3 stats per family (for legacy pickle)
        self.family_cdr3_length: Dict[str, WelfordAccumulator] = defaultdict(WelfordAccumulator)
        self.family_cdr3_charge: Dict[str, WelfordAccumulator] = defaultdict(WelfordAccumulator)
        self.family_cdr3_aromatic: Dict[str, WelfordAccumulator] = defaultdict(WelfordAccumulator)
    
    def _get_corr_acc(self, family: str, pos: int, name: str) -> CorrAccumulator:
        d = self.family_pos_corr[family][pos]
        if name not in d:
            d[name] = CorrAccumulator()
        return d[name]
    
    def add(self, imgt: Dict[int, str], family: str, cdr3_info: Dict):
        self.family_totals[family] += 1
        
        L = float(cdr3_info.get("length", 0))
        q = float(cdr3_info.get("charge", 0))
        h = float(cdr3_info.get("hydrophobicity", 0))
        arom = float(cdr3_info.get("aromatic_count", 0))
        
        # Track CDR3 stats per family
        self.family_cdr3_length[family].add(L)
        self.family_cdr3_charge[family].add(q)
        self.family_cdr3_aromatic[family].add(arom)
        
        for pos in self.positions:
            aa = imgt.get(pos, "-")
            self.family_pos_aa_counts[family][pos][aa] += 1
            
            # Skip gaps for correlation accumulation (avoids bias from variable gap rates)
            if aa not in VALID_AA:
                continue
            
            # Boolean FR properties (only for valid AAs)
            is_hydro = 1.0 if aa in HYDROPHOBICITY and HYDROPHOBICITY.get(aa, 0) > 0 else 0.0
            is_pos = 1.0 if aa in AA_POS else 0.0
            is_neg = 1.0 if aa in AA_NEG else 0.0
            
            # Track correlations (gaps excluded)
            self._get_corr_acc(family, pos, "cdr3_len_vs_is_hydrophobic").add(L, is_hydro)
            self._get_corr_acc(family, pos, "cdr3_len_vs_is_positive").add(L, is_pos)
            self._get_corr_acc(family, pos, "cdr3_charge_vs_is_positive").add(q, is_pos)
            self._get_corr_acc(family, pos, "cdr3_charge_vs_is_negative").add(q, is_neg)
            self._get_corr_acc(family, pos, "cdr3_hydro_vs_is_hydrophobic").add(h, is_hydro)
    
    def get_family_cdr3_stats(self, family: str) -> Dict:
        """Get CDR3 stats for a family (for legacy pickle)."""
        return {
            "cdr3_length": self.family_cdr3_length[family].stats(),
            "cdr3_charge": self.family_cdr3_charge[family].stats(),
            "cdr3_aromatic": self.family_cdr3_aromatic[family].stats(),
        }
    
    def summarize(self) -> Dict:
        """
        Output structure:
        {
          "by_family": {
            "F_C2": {
              "n_sequences": 123456,
              "positions": [
                {
                  "position_num": 71,
                  "n_obs": 120000,
                  "major_aa": "V",
                  "major_freq": 0.72,
                  "entropy": 1.23,
                  "correlations": {
                    "cdr3_len_vs_is_hydrophobic": 0.31,
                    ...
                  }
                },
                ...
              ]
            },
            ...
          }
        }
        """
        out = {"by_family": {}}
        
        for family, pos_map in self.family_pos_aa_counts.items():
            fam_entry = {
                "family": family,
                "n_sequences": int(self.family_totals.get(family, 0)),
                "positions": [],
            }
            
            for pos in sorted(pos_map.keys()):
                aa_counts = pos_map[pos]
                n_obs_total = sum(aa_counts.values())
                if n_obs_total == 0:
                    continue
                
                # Filter out gaps for conservation metrics
                valid_counts = {a: c for a, c in aa_counts.items() if a in VALID_AA}
                n_obs_valid = sum(valid_counts.values())
                
                # Skip all-gap positions
                if n_obs_valid == 0:
                    continue
                
                # Major AA and conservation (excluding gaps)
                major_aa, major_count = max(valid_counts.items(), key=lambda x: x[1])
                major_freq = major_count / n_obs_valid
                
                # Gap rate for reference
                gap_rate = 1.0 - (n_obs_valid / n_obs_total) if n_obs_total > 0 else 0.0
                
                # Shannon entropy (bits) - computed on valid AAs only
                entropy = 0.0
                for a, c in valid_counts.items():
                    p = c / n_obs_valid
                    if p > 0:
                        entropy -= p * math.log(p, 2)
                
                # Correlations
                corr_dict = {}
                corr_map = self.family_pos_corr.get(family, {}).get(pos, {})
                for name, acc in corr_map.items():
                    r = acc.corr()
                    if r is not None:
                        corr_dict[name] = round(r, 4)
                
                fam_entry["positions"].append({
                    "position_num": pos,
                    "n_obs": int(n_obs_valid),  # Valid AAs only
                    "n_obs_total": int(n_obs_total),  # Including gaps
                    "gap_rate": round(gap_rate, 4),
                    "major_aa": major_aa,
                    "major_freq": round(major_freq, 4),
                    "entropy": round(entropy, 4),
                    "correlations": corr_dict,
                })
            
            fam_entry["positions"].sort(key=lambda x: x["position_num"])
            out["by_family"][family] = fam_entry
        
        return out

# ============================================================
# MUTUAL INFORMATION ANALYZER
# ============================================================

class MIAnalyzer:
    """Compute pairwise MI between positions on a subsample (vectorized)."""
    
    def __init__(self, positions: List[int], max_seqs: int = 2_000_000):
        self.positions = positions
        self.max_seqs = max_seqs
        self.n_added = 0
        self.enabled = HAS_NUMPY and max_seqs > 0
        
        if self.enabled:
            self._pos_index = {p: i for i, p in enumerate(self.positions)}
            self._aa_to_idx = {aa: i for i, aa in enumerate(sorted(VALID_AA | {"-"}))}
            self._data = []
    
    def add(self, imgt: Dict[int, str]):
        if not self.enabled or self.n_added >= self.max_seqs:
            return
        
        row = []
        for p in self.positions:
            aa = imgt.get(p, "-")
            if aa not in self._aa_to_idx:
                aa = "-"
            row.append(self._aa_to_idx[aa])
        
        self._data.append(row)
        self.n_added += 1
    
    def compute_mi(self, top_k: int = 500) -> List[Dict]:
        if not self.enabled or self.n_added == 0:
            return []
        
        data = np.array(self._data, dtype=np.int16)
        N, P = data.shape
        
        max_aa = len(self._aa_to_idx)
        
        # Compute marginals
        marginals = []
        for j in range(P):
            counts = np.bincount(data[:, j], minlength=max_aa).astype(float)
            probs = counts / N
            marginals.append(probs)
        
        results = []
        for i in range(P):
            for j in range(i + 1, P):
                x = data[:, i]
                y = data[:, j]
                
                # Vectorized joint distribution using np.add.at
                joint = np.zeros((max_aa, max_aa), dtype=float)
                np.add.at(joint, (x, y), 1.0)
                joint /= N
                
                px = marginals[i]
                py = marginals[j]
                
                # Vectorized MI calculation
                outer = np.outer(px, py) + 1e-12
                # Only compute where joint > 0
                mask = joint > 0
                mi_contrib = np.where(mask, joint * np.log2(joint / outer + 1e-12), 0.0)
                mi = mi_contrib.sum()
                
                if mi > 0:
                    results.append({
                        "pos1": self.positions[i],
                        "pos2": self.positions[j],
                        "mi": round(mi, 4),
                        "source": "vhh_analysis_v7",
                    })
        
        results.sort(key=lambda r: -r["mi"])
        return results[:top_k] if top_k else results

# ============================================================
# TRIPLET ANALYZER
# ============================================================

class TripletAnalyzer:
    """CDR3/FR4 junction triplet patterns."""
    
    def __init__(self):
        self.counts: Dict[str, Dict[str, Counter]] = defaultdict(lambda: defaultdict(Counter))
    
    def update(self, family: str, cdr3: str, fr4_seq: str):
        if len(cdr3) >= 3 and len(fr4_seq) >= 3:
            lhs, rhs = cdr3[-3:], fr4_seq[:3]
            if all(c in VALID_AA for c in lhs + rhs):
                self.counts[family][lhs][rhs] += 1
    
    def extract_rules(self, min_support: int, min_confidence: float, 
                      min_lift: float) -> List[Dict]:
        rules = []
        
        for family, lhs_dict in self.counts.items():
            # Baseline for this family
            all_rhs = Counter()
            for rhs_counts in lhs_dict.values():
                all_rhs.update(rhs_counts)
            total_all = sum(all_rhs.values())
            if total_all < min_support:
                continue
            
            for lhs, rhs_counts in lhs_dict.items():
                total_lhs = sum(rhs_counts.values())
                if total_lhs < min_support // 5:
                    continue
                
                rhs, count = max(rhs_counts.items(), key=lambda kv: kv[1])
                conf = count / total_lhs
                if conf < min_confidence:
                    continue
                
                base_prob = all_rhs.get(rhs, 0) / total_all
                lift = (conf / base_prob) if base_prob > 0 else 999.99
                if lift < min_lift:
                    continue
                
                rules.append({
                    "condition": f"family={family} AND cdr3_end={lhs}",
                    "condition_tokens": [f"family={family}", f"cdr3_end={lhs}"],
                    "position": "TRIPLET:cdr3_fr4",
                    "position_num": 0,
                    "suggested_aa": rhs,
                    "confidence": round(conf, 4),
                    "confidence_all": round(conf, 4),  # Same as confidence for triplets (no gaps)
                    "support": int(count),
                    "lift": round(lift, 4),
                    "baseline_prob": round(base_prob, 4),
                    "n_total_for_condition": int(total_lhs),
                    "source": "vhh_analysis_v7",
                    "rule_type": "triplet",
                    "family": family,
                })
        
        rules.sort(key=lambda r: (-r["confidence"], -r["support"]))
        return rules

# ============================================================
# MAIN
# ============================================================

def main():
    ap = argparse.ArgumentParser(
        description="VHH Analysis / Compensation Unified v7",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Single CSV
  python vhh_analysis_unified_v7.py --csv data.csv -o output/

  # Multiple CSVs
  python vhh_analysis_unified_v7.py --csv data1.csv data2.csv data3.csv -o output/
  
  # Glob pattern (shell expansion)
  python vhh_analysis_unified_v7.py --csv data/*.csv -o output/
  
  # With legacy pickles for existing scripts
  python vhh_analysis_unified_v7.py --csv data/*.csv -o output/ --emit-legacy-pickles
"""
    )
    ap.add_argument("--csv", required=True, nargs="+",
                    help="IMGT-annotated VHH CSV file(s) - accepts multiple files or glob patterns")
    ap.add_argument("-o", "--output", required=True, help="Output directory")
    ap.add_argument("--positions", choices=["all", "vernier", "fr3fr4"], default="vernier")
    ap.add_argument("--conditions", choices=["minimal", "full"], default="minimal")
    ap.add_argument("--filter-truncated", action="store_true")
    ap.add_argument("--keep-definitional", action="store_true")
    # Separate support thresholds (ChatGPT recommendation)
    ap.add_argument("--min-condition-obs", type=int, default=5000,
                    help="Minimum observations for condition+position (default: 5000)")
    ap.add_argument("--min-aa-support", type=int, default=1000,
                    help="Minimum support for specific AA (default: 1000, can be lower than min-condition-obs)")
    ap.add_argument("--min-confidence", type=float, default=0.7)
    ap.add_argument("--min-lift", type=float, default=1.15)
    ap.add_argument("--mi-max-seqs", type=int, default=2000000)
    ap.add_argument("--mi-top-k", type=int, default=500)
    ap.add_argument("--print-every", type=int, default=200000)
    ap.add_argument("--checkpoint-interval", type=int, default=1000000,
                    help="Save checkpoint every N rows (0 to disable)")
    ap.add_argument("--resume", action="store_true",
                    help="Resume from latest checkpoint if available")
    # MI is off by default (visualization only)
    ap.add_argument("--enable-mi", action="store_true",
                    help="Enable MI computation (off by default - visualization only)")
    # Legacy pickle output for existing scripts
    ap.add_argument("--emit-legacy-pickles", action="store_true",
                    help="Also emit legacy pickle files for designer/naturalness scripts")
    args = ap.parse_args()
    
    import time
    import glob as glob_module
    total_start = time.time()
    
    # Expand glob patterns in CSV list
    csv_files = []
    for pattern in args.csv:
        expanded = glob_module.glob(pattern)
        if expanded:
            csv_files.extend(expanded)
        elif os.path.exists(pattern):
            csv_files.append(pattern)
        else:
            print(f"Warning: {pattern} not found, skipping")
    
    if not csv_files:
        print("ERROR: No CSV files found")
        return
    
    csv_files = sorted(set(csv_files))  # Dedupe and sort
    
    os.makedirs(args.output, exist_ok=True)
    checkpoint_path = os.path.join(args.output, "checkpoint.pkl")
    
    print("=" * 70)
    print("VHH ANALYSIS / COMPENSATION UNIFIED v7")
    print("=" * 70)
    print(f"Input files: {len(csv_files)}")
    for f in csv_files[:5]:
        print(f"  - {f}")
    if len(csv_files) > 5:
        print(f"  ... and {len(csv_files) - 5} more")
    print(f"Output: {args.output}")
    print(f"Positions: {args.positions}")
    print(f"Filter truncated: {args.filter_truncated}")
    print(f"Min condition obs: {args.min_condition_obs:,}")
    print(f"Min AA support: {args.min_aa_support:,}")
    print(f"Checkpoint interval: {args.checkpoint_interval:,} rows")
    print(f"Skip MI: {not args.enable_mi}")
    print()
    
    # Read header from first CSV (assume all CSVs have same schema)
    with open(csv_files[0], "r", newline="") as fh:
        header_raw = next(csv.reader(fh))
    header = [norm_col(h) for h in header_raw]
    idx = {h: i for i, h in enumerate(header)}
    
    # Find IMGT columns
    available_imgt: Set[int] = set()
    pos_idx: Dict[int, int] = {}
    for h in header:
        if h.startswith("imgt_"):
            try:
                p = int(h.split("_")[1])
                available_imgt.add(p)
                pos_idx[p] = idx[h]
            except:
                pass
    
    if not available_imgt:
        raise SystemExit("No imgt_* columns found in CSV header.")
    
    target_positions = choose_positions(args.positions, available_imgt, exclude_truncation=True)
    
    print(f"IMGT positions available: {len(available_imgt)}")
    print(f"Target positions: {len(target_positions)}")
    print(f"Excluded: {sorted(EXCLUDE_FROM_RULES & available_imgt)}")
    print()
    
    # Column indices
    aa_idx = idx.get("aa_v_full")
    cdr1_idx, cdr2_idx, cdr3_idx = idx.get("cdr1"), idx.get("cdr2"), idx.get("cdr3")
    fam_idx, hm_idx = idx.get("family"), idx.get("hallmarks")
    
    if cdr3_idx is None:
        raise SystemExit("Missing required column: cdr3")
    
    # Count total rows across all files for ETA
    print("Counting rows for ETA calculation...")
    t0 = time.time()
    total_rows = 0
    for csv_file in csv_files:
        with open(csv_file, "r") as f:
            total_rows += sum(1 for _ in f) - 1  # Subtract header
    print(f"Total rows across {len(csv_files)} file(s): {total_rows:,} (counted in {time.time()-t0:.1f}s)")
    print()
    
    # Initialize analyzers (conditionally)
    comp_analyzer = CompensationAnalyzer(target_positions)
    corr_analyzer = CorrelationAnalyzer(target_positions)
    mi_analyzer = MIAnalyzer(target_positions, max_seqs=args.mi_max_seqs) if args.enable_mi else None
    triplet_analyzer = TripletAnalyzer()
    
    total = filtered = used = 0
    start_row = 0
    
    # Resume from checkpoint if requested
    if args.resume and os.path.exists(checkpoint_path):
        print(f"Loading checkpoint from {checkpoint_path}...")
        try:
            with open(checkpoint_path, "rb") as f:
                ckpt = pickle.load(f)
            comp_analyzer = ckpt["comp_analyzer"]
            corr_analyzer = ckpt.get("corr_analyzer", corr_analyzer)
            mi_analyzer = ckpt.get("mi_analyzer") if args.enable_mi else None
            triplet_analyzer = ckpt["triplet_analyzer"]
            total = ckpt["total"]
            filtered = ckpt["filtered"]
            used = ckpt["used"]
            start_row = ckpt["total"]  # Skip this many rows
            print(f"Resumed from row {start_row:,} (used={used:,}, filtered={filtered:,})")
        except Exception as e:
            print(f"WARNING: Could not load checkpoint: {e}")
            print("Starting from scratch...")
            start_row = 0
    
    def save_checkpoint():
        """Save current state to checkpoint file."""
        ckpt = {
            "comp_analyzer": comp_analyzer,
            "corr_analyzer": corr_analyzer,
            "mi_analyzer": mi_analyzer,
            "triplet_analyzer": triplet_analyzer,
            "total": total,
            "filtered": filtered,
            "used": used,
        }
        tmp_path = checkpoint_path + ".tmp"
        with open(tmp_path, "wb") as f:
            pickle.dump(ckpt, f)
        os.replace(tmp_path, checkpoint_path)  # Atomic replace
        print(f"  [Checkpoint saved at row {total:,}]")
    
    def format_time(seconds):
        """Format seconds as HH:MM:SS."""
        h, rem = divmod(int(seconds), 3600)
        m, s = divmod(rem, 60)
        return f"{h:02d}:{m:02d}:{s:02d}"
    
    last_checkpoint = start_row
    loop_start = time.time()
    rows_to_process = total_rows - start_row
    
    print("=" * 70)
    print("PROCESSING DATA")
    print("=" * 70)
    
    # Process all CSV files
    global_row = 0  # Track position across all files
    for csv_file in csv_files:
        print(f"  Processing: {csv_file}")
        with open(csv_file, "r", newline="") as fh:
            reader = csv.reader(fh)
            next(reader)  # Skip header
            
            for row in reader:
                total += 1
                global_row += 1
                
                # Skip rows if resuming
                if total <= start_row:
                    continue
                
                imgt = {p: norm_cell(row[i] if i < len(row) else "-") for p, i in pos_idx.items()}
                
                if args.filter_truncated and not passes_truncation_filter(imgt):
                    filtered += 1
                    continue
                
                aa_v_full = row[aa_idx] if aa_idx is not None and aa_idx < len(row) else ""
                cdr3 = row[cdr3_idx] if cdr3_idx < len(row) else ""
                
                # Fix: use 'is not None' to handle index 0 correctly
                hallmarks = row[hm_idx] if (hm_idx is not None and hm_idx < len(row) and row[hm_idx]) else infer_hallmarks(imgt)
                family = row[fam_idx] if (fam_idx is not None and fam_idx < len(row) and row[fam_idx]) else infer_family(imgt, aa_v_full)
                
                cdr3_info = categorize_cdr3(cdr3)
                
                # Update analyzers
                comp_analyzer.update(imgt, cdr3, family, hallmarks, args.conditions)
                corr_analyzer.add(imgt, family, cdr3_info)
                if mi_analyzer:
                    mi_analyzer.add(imgt)
                
                # Triplet analysis
                fr4_seq = ''.join(imgt.get(p, '') for p in range(118, 126))
                triplet_analyzer.update(family, cdr3, fr4_seq)
                
                used += 1
                
                if total % args.print_every == 0:
                    elapsed = time.time() - loop_start
                    rows_done = total - start_row
                    rate = rows_done / elapsed if elapsed > 0 else 0
                    remaining = (rows_to_process - rows_done) / rate if rate > 0 else 0
                    pct = 100 * total / total_rows
                    mi_str = f"MI: {mi_analyzer.n_added:,}" if mi_analyzer else "MI: skip"
                    print(f"  [{pct:5.1f}%] Rows: {total:,} | Used: {used:,} | {mi_str} | "
                          f"Elapsed: {format_time(elapsed)} | ETA: {format_time(remaining)}")
                
                # Save checkpoint periodically
                if args.checkpoint_interval > 0 and (total - last_checkpoint) >= args.checkpoint_interval:
                    save_checkpoint()
                    last_checkpoint = total
    
    # Final checkpoint before extraction
    if args.checkpoint_interval > 0:
        save_checkpoint()
    
    loop_elapsed = time.time() - loop_start
    print()
    print(f"Data processing complete in {format_time(loop_elapsed)}")
    
    print()
    print("=" * 70)
    print("EXTRACTING RESULTS")
    print("=" * 70)
    print(f"Total rows: {total:,}")
    print(f"Filtered (truncation): {filtered:,}")
    print(f"Used: {used:,}")
    mi_seqs = mi_analyzer.n_added if mi_analyzer else 0
    print(f"MI sequences: {mi_seqs:,}" + (" (skipped)" if not args.enable_mi else ""))
    print()
    
    # Extract compensation rules (with separate thresholds)
    t0 = time.time()
    comp_rules = comp_analyzer.extract_rules(
        min_condition_obs=args.min_condition_obs,
        min_aa_support=args.min_aa_support,
        min_confidence=args.min_confidence,
        min_lift=args.min_lift,
        filter_definitional=not args.keep_definitional,
    )
    print(f"Compensation rules: {len(comp_rules):,} ({time.time()-t0:.1f}s)")
    
    # Extract triplet rules
    t0 = time.time()
    triplet_rules = triplet_analyzer.extract_rules(
        min_support=args.min_aa_support,  # Use AA support threshold
        min_confidence=args.min_confidence,
        min_lift=args.min_lift,
    )
    print(f"Triplet rules: {len(triplet_rules):,} ({time.time()-t0:.1f}s)")
    
    # Extract vernier archetypes
    t0 = time.time()
    archetypes = comp_analyzer.vernier_archetypes(min_support=max(args.min_condition_obs, 10000))
    print(f"Vernier archetypes: {len(archetypes):,} ({time.time()-t0:.1f}s)")
    
    # Extract MI pairs (if enabled)
    mi_pairs = []
    if mi_analyzer and args.enable_mi:
        t0 = time.time()
        mi_pairs = mi_analyzer.compute_mi(top_k=args.mi_top_k)
        print(f"MI pairs: {len(mi_pairs):,} ({time.time()-t0:.1f}s)")
    else:
        print("MI pairs: skipped")
    
    # Extract correlations (conservation + CDR3↔FR relationships)
    t0 = time.time()
    corr_summary = corr_analyzer.summarize()
    print(f"Correlation families: {len(corr_summary.get('by_family', {})):,} ({time.time()-t0:.1f}s)")
    
    # Combine all rules
    all_rules = comp_rules + triplet_rules
    
    # Build summary
    summary = {
        "schema_version": "v1",
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "metadata": {
            "input_csv": csv_files,  # Now a list
            "n_csv_files": len(csv_files),
            "total_rows": total,
            "filtered_truncation": filtered,
            "used": used,
            "mi_sequences": mi_seqs,
            "skip_mi": not args.enable_mi,
            "positions_mode": args.positions,
            "conditions_mode": args.conditions,
            "target_positions": target_positions,
            "excluded_positions": sorted(EXCLUDE_FROM_RULES),
        },
        "params": {
            "min_condition_obs": args.min_condition_obs,
            "min_aa_support": args.min_aa_support,
            "min_confidence": args.min_confidence,
            "min_lift": args.min_lift,
            "mi_max_seqs": args.mi_max_seqs,
            "mi_top_k": args.mi_top_k,
        },
        "family_counts": dict(comp_analyzer.family_counts),
        "hallmark_counts": dict(comp_analyzer.hallmark_counts.most_common(50)),
        "n_compensation_rules": len(comp_rules),
        "n_triplet_rules": len(triplet_rules),
        "n_archetypes": len(archetypes),
        "n_mi_pairs": len(mi_pairs),
        "n_corr_families": len(corr_summary.get("by_family", {})),
    }
    
    # Write outputs
    print()
    print("=" * 70)
    print("WRITING OUTPUT FILES")
    print("=" * 70)
    
    t0 = time.time()
    with open(os.path.join(args.output, "analysis_summary_v7.json"), "w") as f:
        json.dump(summary, f, indent=2)
    print(f"  analysis_summary_v7.json ({time.time()-t0:.1f}s)")
    
    t0 = time.time()
    with open(os.path.join(args.output, "analysis_rules_v7.json"), "w") as f:
        json.dump(all_rules, f, indent=2)
    print(f"  analysis_rules_v7.json - {len(all_rules):,} rules ({time.time()-t0:.1f}s)")
    
    t0 = time.time()
    with open(os.path.join(args.output, "analysis_vernier_archetypes_v7.json"), "w") as f:
        json.dump(archetypes, f, indent=2)
    print(f"  analysis_vernier_archetypes_v7.json ({time.time()-t0:.1f}s)")
    
    if mi_pairs and args.enable_mi:
        t0 = time.time()
        with open(os.path.join(args.output, "analysis_mi_pairs_v7.json"), "w") as f:
            json.dump(mi_pairs, f, indent=2)
        print(f"  analysis_mi_pairs_v7.json ({time.time()-t0:.1f}s)")
    
    t0 = time.time()
    with open(os.path.join(args.output, "analysis_correlations_v7.json"), "w") as f:
        json.dump(corr_summary, f, indent=2)
    print(f"  analysis_correlations_v7.json ({time.time()-t0:.1f}s)")
    
    # Legacy pickle export for designer/naturalness scripts
    if args.emit_legacy_pickles:
        print()
        print("  [Legacy pickle export]")
        t0 = time.time()
        
        def imgt_to_fw_position(imgt_num: int) -> str:
            """Convert IMGT position number to legacy FRx_y format."""
            if 1 <= imgt_num <= 26:
                return f"FR1_{imgt_num}"
            elif 39 <= imgt_num <= 55:
                return f"FR2_{imgt_num - 38}"  # FR2 starts at IMGT39
            elif 66 <= imgt_num <= 104:
                return f"FR3_{imgt_num - 65}"  # FR3 starts at IMGT66
            elif 118 <= imgt_num <= 128:
                return f"FR4_{imgt_num - 117}"  # FR4 starts at IMGT118
            else:
                return f"IMGT{imgt_num}"
        
        # 1. Compensation pickle (multi_position_rules, vernier_archetypes)
        # Transform rules to legacy format
        legacy_multi_rules = []
        for r in all_rules:
            if r.get("rule_type") == "cdr_to_imgt":
                legacy_multi_rules.append({
                    "condition": r["condition"],
                    "fw_position": imgt_to_fw_position(r["position_num"]),
                    "suggested_aa": r["suggested_aa"],
                    "confidence": r["confidence"] * 100,  # Legacy uses percentages
                    "support": r["support"],
                    "lift": r["lift"],
                    "family": r["family"],
                })
        
        # Transform archetypes to legacy dict format
        legacy_archetypes = {}
        for arch in archetypes:
            family = arch["family"]
            pattern = {}
            for pos_str, data in arch.get("positions", {}).items():
                pattern[pos_str] = data["consensus"]
            legacy_archetypes[family] = {
                "pattern": pattern,
                "count": arch["count"],
            }
        
        compensation_pkl = {
            "multi_position_rules": legacy_multi_rules,
            "vernier_archetypes": legacy_archetypes,
        }
        pkl_path = os.path.join(args.output, "analysis_compensation_legacy.pkl")
        with open(pkl_path, "wb") as f:
            pickle.dump(compensation_pkl, f)
        print(f"    analysis_compensation_legacy.pkl ({time.time()-t0:.1f}s)")
        
        # 2. Epistasis pickle (analysis_2_vernier_clusters, analysis_4_higher_order_rules)
        t0 = time.time()
        
        # Build vernier clusters from correlation data
        # Structure expected: {pattern_str: {'pattern': {...}, 'cdr3_length': {'n':, 'mean':, 'std':}, ...}}
        legacy_clusters = {}
        
        # Group by family vernier pattern  
        for family, fam_data in corr_summary.get("by_family", {}).items():
            # Build pattern string from vernier positions
            pattern = {}
            for pos_data in fam_data.get("positions", []):
                pos_num = pos_data.get("position_num")
                if pos_num in VERNIER_POSITIONS:
                    pattern[f"IMGT{pos_num}"] = pos_data.get("major_aa", "-")
            
            if not pattern:
                continue
            
            pattern_str = "_".join(f"{k}={v}" for k, v in sorted(pattern.items()))
            
            # Get actual CDR3 stats from CorrelationAnalyzer
            cdr3_stats = corr_analyzer.get_family_cdr3_stats(family)
            
            legacy_clusters[pattern_str] = {
                "pattern": pattern,
                "cdr3_length": cdr3_stats["cdr3_length"],
                "cdr3_charge": cdr3_stats["cdr3_charge"],
                "cdr3_aromatic": cdr3_stats["cdr3_aromatic"],
            }
        
        # Higher-order rules - transform compound conditions
        legacy_higher_order = []
        for r in all_rules:
            if len(r.get("condition_tokens", [])) >= 2 or "AND" in r.get("condition", ""):
                legacy_higher_order.append({
                    "condition": r["condition"],
                    "result": f"{r['position']}={r['suggested_aa']}",
                    "confidence": r["confidence"],
                    "support": r["support"],
                    "family": r["family"],
                })
        
        epistasis_pkl = {
            "analysis_2_vernier_clusters": legacy_clusters,
            "analysis_4_higher_order_rules": legacy_higher_order,
        }
        pkl_path = os.path.join(args.output, "analysis_epistasis_legacy.pkl")
        with open(pkl_path, "wb") as f:
            pickle.dump(epistasis_pkl, f)
        print(f"    analysis_epistasis_legacy.pkl ({time.time()-t0:.1f}s)")
    
    # Print top results
    print()
    print("=" * 70)
    print("TOP 15 RULES")
    print("=" * 70)
    for r in comp_rules[:15]:
        print(f"  {r['position']:8s} → {r['suggested_aa']} | {r['condition'][:35]:35s} | "
              f"conf={r['confidence']:.3f} lift={r['lift']:.2f} n={r['support']:,} [{r['family']}]")
    
    print()
    print("=" * 70)
    print("FAMILY DISTRIBUTION")
    print("=" * 70)
    for fam, count in comp_analyzer.family_counts.most_common():
        pct = 100 * count / used if used else 0
        print(f"  {fam}: {count:,} ({pct:.1f}%)")
    
    if mi_pairs and args.enable_mi:
        print()
        print("=" * 70)
        print("TOP 10 MI PAIRS")
        print("=" * 70)
        for p in mi_pairs[:10]:
            print(f"  IMGT{p['pos1']} <-> IMGT{p['pos2']}: MI = {p['mi']:.4f}")
    
    total_elapsed = time.time() - total_start
    print()
    print("=" * 70)
    print(f"COMPLETE - Total time: {format_time(total_elapsed)}")
    print("=" * 70)

if __name__ == "__main__":
    main()
