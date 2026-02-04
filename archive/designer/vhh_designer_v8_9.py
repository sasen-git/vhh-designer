#!/usr/bin/env python3
"""
VHH Designer v8.9 - Restructured Track System with Capped Layers
================================================================

MAJOR UPDATES from v8.8.3:

1. SIMPLIFIED YQRL CONTROL STRUCTURE (v8.9):
   - 0a: GRAFT (0 mutations) - pure CDR graft
   - 0b: Hallmarks + IMGT2 (5 mutations) - ALWAYS include IMGT2 with hallmarks
   - 0c: Hallmarks + IMGT2 + ALL verniers (~20 mutations)
   - 0d+: Progressive framework in PAIRS (2 at a time), all positions >80% consensus
   - Removed separate "hallmarks only (4 mut)" - always include IMGT2

2. ALL HALLMARKS GET GRAFT FIRST (v8.9):
   - 0a: GRAFT (0 mutations) - always first for ALL hallmarks
   - 0b: Hallmarks + IMGT2 (5 mutations)
   - 0c: + ALL verniers (only if hallmark avg vernier consensus ≥60%)

3. TRACK 1 EXPANDED (v8.9):
   - Now includes ALL verniers, not just "load-bearing"
   - Expanded from max 12 to max 20 probes per hallmark
   - Tests each vernier position individually

4. TRACK 2 REDEFINED: "MINIMAL STRUCTURED PERTURBATIONS" (v8.9):
   - Base: Hallmarks + IMGT2 (IMGT2 now protected like hallmarks)
   - Add: Exactly ONE 2-3 position motif (weighted by family support)
   - Add: 0-K extra vernier mutations (K_EXTRA_VERNIERS_T2 = 2, capped)
   - NO framework tier sprinkling - keeps it interpretable

5. TRACK 3 REDEFINED: "MOTIF ATTRIBUTION" (v8.9):
   - Base: Hallmarks + IMGT2
   - Add: Exactly ONE 3-5 position motif (weighted)
   - Add: 0-K extra vernier mutations (K_EXTRA_VERNIERS_T3 = 2, max 3)
   - NO framework tier sprinkling
   - Goal: Interpretable - can analyze "motif X worked/didn't work"

6. TRACK 4 REDEFINED: "PRODUCTION OPTIMIZER" with CAPPED LAYERS (v8.9):
   - Full probabilistic layering for production candidates
   - Vernier layer: Bernoulli + cap (K_VERNIER_MAX = 8)
   - Framework tiers applied in order with individual caps:
     - Tier S (≥95%): K_S_MAX = 10
     - Tier A (90-95%): K_A_MAX = 6
     - Tier B (85-90%): K_B_MAX = 3
     - Tier C (75-85%): K_C_MAX = 2 (only when desperate)
   - CDR-conditional rules: capped at K_CDR_RULES_MAX = 2
   - Adaptive expansion when uniqueness drops

7. UNIVERSAL SCAFFOLD CONTROLS (v8.9):
   - 0a: Complete CDR graft into Universal (0 mutations)
   - 0b: Hallmarks + IMGT2 (5 mutations)
   - 0c: Hallmarks + IMGT2 + ALL verniers

8. IMGT2 NOW PROTECTED LIKE HALLMARKS (v8.9):
   - PROTECTED_POSITIONS = {2, 42, 49, 50, 52}
   - Motifs cannot override these positions

KEY DISTINCTION:
- Track 3 = "Motif attribution": one motif + minimal extras → interpretable
- Track 4 = "Production optimizer": full probabilistic layering → coverage

TRACK STRUCTURE (v8.9):

  ORIGINAL FRAMEWORK:
  Track 0: Controls (exempt)
    - 0a: GRAFT (0 mutations) - ALL hallmarks get this
    - 0b: Hallmarks + IMGT2 (5 mutations)
    - 0c: + ALL verniers (~20 mutations) - if avg vernier ≥60%
    - 0d+: YQRL only: + Progressive FW pairs (>80% consensus)
  Track 1: Single vernier probes - ALL verniers (RANKED)
  Track 2: Minimal structured perturbations - 1 motif + 0-2 verniers (RANKED)
  Track 3: Motif attribution - 1 motif (3-5 pos) + 0-3 verniers (RANKED)
  Track 4: Production optimizer - capped probabilistic layers (RANKED)

  UNIVERSAL FRAMEWORK:
  Track 0: Controls (exempt)
    - 0a: Complete CDR graft (0 mutations)
    - 0b: Hallmarks + IMGT2 (5 mutations)
    - 0c: + ALL verniers (~20 mutations)
  Track 1-4: Similar to Original but CDR-conditional focused

Usage:
  python vhh_designer_v8_9.py -i input.fasta \\
      --rules analysis_rules_v7_all_positions.json \\
      --archetypes analysis_vernier_archetypes_v7.json \\
      --hallmark-db comprehensive_subfamily_analysis_imgt.xlsx \\
      --target-hallmarks AUTO \\
      --n-generate 200000 \\
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
# V8.6: HELPER FUNCTIONS FOR IMPROVED GENERATION
# ============================================================
import math

# Type aliases for clarity
BucketKey = Tuple[str, str, str]  # (hallmark, scaffold_type, family)

def _extract_hallmarks(positions: Dict) -> str:
    """Extract hallmark string from IMGT positions (local helper)."""
    p42 = positions.get(42, positions.get('42', '?'))
    p49 = positions.get(49, positions.get('49', '?'))
    p50 = positions.get(50, positions.get('50', '?'))
    p52 = positions.get(52, positions.get('52', '?'))
    return f"{p42}{p49}{p50}{p52}"

def bucket_key(cand) -> BucketKey:
    """Get bucket key for a candidate: (hallmark, scaffold_type, family)."""
    hallmarks = getattr(cand, 'hallmarks', None)
    if not hallmarks:
        hallmarks = _extract_hallmarks(getattr(cand, 'imgt_positions', {}))
    scaffold = getattr(cand, 'scaffold_type', 'original')
    family = getattr(cand, 'family', '')
    return (hallmarks, scaffold, family)

def track_number(cand) -> int:
    """
    Map design_track string to track number.
    Track 0 = controls (exempt), Tracks 1-4 = ranked lanes.
    """
    dt = getattr(cand, 'design_track', '')
    track_str = getattr(cand, 'track', '')
    
    # Use explicit track field if available
    if track_str:
        if '0' in track_str:
            return 0
        elif '1' in track_str:
            return 1
        elif '2' in track_str:
            return 2
        elif '3' in track_str:
            return 3
        elif '4' in track_str:
            return 4
    
    # Fallback to design_track mapping
    mapping = {
        "minimal_hallmark": 0,
        "single_vernier": 1,
        "paired_vernier": 2,
        "triplet_vernier": 3,
        "motif_combo": 3,
        "optimized": 4,
        "lead": -1,
    }
    return mapping.get(dt, 4)


def temp_scale_prob(p: float, temperature: float = 1.0) -> float:
    """
    Temperature-scaled probability using logit transform.
    
    temperature=1.0 -> unchanged
    temperature<1.0 -> sharper (more aggressive toward 0/1)
    temperature>1.0 -> softer (closer to 0.5)
    """
    p = min(max(p, 1e-6), 1 - 1e-6)
    logit = math.log(p / (1 - p))
    logit_scaled = logit / max(temperature, 1e-6)
    return 1 / (1 + math.exp(-logit_scaled))



def weighted_sample_without_replacement(items, weights, k, rng):
    """Sample up to k distinct items without replacement, proportional to weights.

    items: list[T]
    weights: list[float], same length as items
    k: number of items to sample
    rng: random.Random

    Returns a list of selected items (length <= k).
    """
    if k <= 0 or not items:
        return []
    pool = []
    for item, w in zip(items, weights):
        try:
            w = float(w)
        except Exception:
            continue
        if w > 0:
            pool.append([item, w])
    if not pool:
        return []
    out = []
    for _ in range(min(k, len(pool))):
        total = sum(w for _, w in pool)
        if total <= 0:
            break
        r = rng.random() * total
        c = 0.0
        for i, (item, w) in enumerate(pool):
            c += w
            if r <= c:
                out.append(item)
                pool.pop(i)
                break
    return out


def sprinkle_positions_bernoulli(
    positions: Dict[int, str],
    sprinkle_map: Dict[int, Tuple[str, float]],
    temperature: float = 1.0,
    p_floor: float = 0.0,
    p_cap: float = 1.0,
    rng: random.Random = None,
) -> Dict[int, str]:
    """
    Apply mutations with independent Bernoulli probability per position.
    
    Args:
        positions: Current sequence map {imgt_pos: aa}
        sprinkle_map: {imgt_pos: (target_aa, freq)} where freq is 0..1
        temperature: Scaling factor for probabilities
        p_floor: Minimum probability
        p_cap: Maximum probability
        rng: Random number generator
        
    Returns:
        New positions dict with mutations applied probabilistically
    """
    if rng is None:
        rng = random.Random()
    
    out = dict(positions)
    for pos, (target_aa, freq) in sprinkle_map.items():
        if pos not in out:
            continue
        if out[pos] == target_aa:
            continue
        p = temp_scale_prob(freq, temperature)
        p = min(max(p, p_floor), p_cap)
        if rng.random() < p:
            out[pos] = target_aa
    return out


def weighted_choice(options: List, weights: List[float], rng: random.Random = None) -> Any:
    """
    Weighted random choice using actual frequencies.
    Replaces 50/50 randomization with frequency-based sampling.
    """
    if rng is None:
        rng = random.Random()
    
    s = sum(weights)
    if s <= 0:
        return rng.choice(options)
    r = rng.random() * s
    acc = 0.0
    for opt, w in zip(options, weights):
        acc += w
        if r <= acc:
            return opt
    return options[-1]


def build_motif_pool_from_archetypes(
    archetypes: Dict[str, Any],
    min_len: int = 2,
    max_len: int = 5
) -> List[List[Tuple[int, str]]]:
    """
    Build motif pool from archetype signatures.
    Each motif is a list of (position, amino_acid) tuples.
    
    Args:
        archetypes: Dict mapping family -> archetype data with 'positions' key
        min_len: Minimum motif size
        max_len: Maximum motif size
        
    Returns:
        List of motifs, each motif is [(pos, aa), ...]
    """
    motifs = []
    
    for family, arch_data in archetypes.items():
        positions = arch_data.get('positions', {})
        if not positions:
            continue
        
        # Build list of (pos, consensus_aa) pairs
        pairs = []
        for pos_key, pos_data in positions.items():
            try:
                pos = int(pos_key.replace('IMGT', ''))
            except (ValueError, AttributeError):
                continue
            
            cons_aa = pos_data.get('consensus', '')
            if not cons_aa:
                continue
            
            pairs.append((pos, cons_aa))
        
        if len(pairs) < min_len:
            continue
        
        # Sort by position for consistent ordering
        pairs.sort(key=lambda x: x[0])
        
        # Create sliding window motifs
        for L in range(min_len, min(max_len + 1, len(pairs) + 1)):
            for start in range(len(pairs) - L + 1):
                motif = pairs[start:start + L]
                if motif not in motifs:
                    motifs.append(motif)
    
    return motifs


def sample_motif(
    motifs: List[List[Tuple[int, str]]],
    min_len: int,
    max_len: int,
    rng: random.Random = None
) -> Optional[List[Tuple[int, str]]]:
    """Sample a motif of specified length range."""
    if rng is None:
        rng = random.Random()
    
    eligible = [m for m in motifs if min_len <= len(m) <= max_len]
    return rng.choice(eligible) if eligible else None


def apply_motif(
    positions: Dict[int, str],
    motif: List[Tuple[int, str]]
) -> Dict[int, str]:
    """Apply a motif to positions dict."""
    out = dict(positions)
    for pos, aa in motif:
        if pos in out:
            out[pos] = aa
    return out


def dedup_within_track(candidates: List) -> List:
    """
    Deduplicate sequences within each design_track.
    Preserves same sequence across different tracks as separate entries.
    """
    by_track = defaultdict(list)
    for c in candidates:
        t = track_number(c)
        by_track[t].append(c)
    
    out = []
    for t, arr in by_track.items():
        seen = set()
        for c in arr:
            seq = c.sequence
            if seq in seen:
                continue
            seen.add(seq)
            out.append(c)
    return out


def label_collisions_keep_all(candidates: List) -> List:
    """Keep all candidates but label collision groups."""
    seq_to_ids = defaultdict(list)
    for c in candidates:
        seq_to_ids[c.sequence].append(c.id)
    
    for c in candidates:
        colliders = seq_to_ids[c.sequence]
        c.dup_group_size = len(colliders)
        c.dup_group_ids = ",".join(colliders[:10])
    return candidates


# Default per-bucket track quotas (v8.6)
# Can be overridden via command line
DEFAULT_TRACK_QUOTAS = {
    0: 2,   # Track 0: minimal hallmark controls (exempt)
    1: 5,   # Track 1: single vernier (ranked)
    2: 6,   # Track 2: motifs 2-3 positions (ranked)
    3: 6,   # Track 3: motifs 3-5 positions (ranked)
    4: 6,   # Track 4: optimized/expansive (ranked)
}
DEFAULT_BUCKET_TARGET = 25


def select_with_track_quotas(
    candidates: List,
    n_select_total: int,
    per_bucket_track_quota: Dict[int, int] = None,
    bucket_target: int = None,
    score_attr: str = "combined_score",
) -> List:
    """
    Select candidates with per-bucket, per-track quotas.
    
    v8.8: Controls (ranking_exempt=True) are ALWAYS included first.
    
    Enforces:
    1. ALL controls (ranking_exempt=True) are included first
    2. Remaining slots filled from ranked candidates (Tracks 1-4)
    3. Per-bucket quotas apply only to ranked candidates
    4. Final output is exactly n_select_total
    
    Args:
        candidates: All candidates to select from
        n_select_total: Exact final count (including controls)
        per_bucket_track_quota: {track_num: max_per_bucket} for ranked tracks
        bucket_target: Max candidates per bucket
        score_attr: Attribute to sort by (higher = better)
        
    Returns:
        Selected candidates: [all controls] + [ranked candidates], exactly n_select_total
    """
    if per_bucket_track_quota is None:
        per_bucket_track_quota = DEFAULT_TRACK_QUOTAS
    if bucket_target is None:
        bucket_target = DEFAULT_BUCKET_TARGET
    
    # =========================================================================
    # STEP 1: SEPARATE CONTROLS FROM RANKED
    # =========================================================================
    controls = [c for c in candidates if getattr(c, 'ranking_exempt', False)]
    ranked = [c for c in candidates if not getattr(c, 'ranking_exempt', False)]
    
    print(f"\n  Selection breakdown:")
    print(f"    Controls (ranking_exempt): {len(controls)}")
    print(f"    Ranked candidates: {len(ranked)}")
    
    # =========================================================================
    # STEP 2: ALWAYS INCLUDE ALL CONTROLS
    # =========================================================================
    selected: List = list(controls)
    n_controls = len(selected)
    
    # Check if controls alone exceed budget
    if n_controls >= n_select_total:
        print(f"    Warning: {n_controls} controls exceed n_select={n_select_total}, truncating")
        selected = selected[:n_select_total]
        return selected
    
    # =========================================================================
    # STEP 3: FILL REMAINING SLOTS FROM RANKED CANDIDATES
    # =========================================================================
    n_ranked_slots = n_select_total - n_controls
    print(f"    Slots for ranked: {n_ranked_slots}")
    
    # Group ranked candidates by bucket then track
    by_bucket_track: Dict[BucketKey, Dict[int, List]] = defaultdict(lambda: defaultdict(list))
    
    for c in ranked:
        b = bucket_key(c)
        t = track_number(c)
        if t < 0:  # Skip lead/special
            continue
        by_bucket_track[b][t].append(c)
    
    # Sort each lane by score (lower = better for combined_score)
    for b, tracks in by_bucket_track.items():
        for t, arr in tracks.items():
            arr.sort(key=lambda x: getattr(x.scoring, score_attr, 0))  # ascending = best first
    
    ranked_selected: List = []
    
    # =========================================================================
    # PHASE 1: ENFORCE GLOBAL MINIMUM PER TRACK (v8.8.1)
    # =========================================================================
    # Ensures each track gets fair representation regardless of per-bucket quotas
    global_min_per_track = {
        1: max(10, n_ranked_slots // 10),  # At least 10 or 10% of ranked slots
        2: max(15, n_ranked_slots // 8),   # At least 15 or 12.5% of ranked slots
        3: max(15, n_ranked_slots // 8),   # At least 15 or 12.5% of ranked slots
        4: max(20, n_ranked_slots // 5),   # At least 20 or 20% of ranked slots
    }
    
    # Collect best candidates from each track globally
    by_track_global: Dict[int, List] = defaultdict(list)
    for b, tracks in by_bucket_track.items():
        for t, arr in tracks.items():
            by_track_global[t].extend(arr)
    
    for t in by_track_global:
        by_track_global[t].sort(key=lambda x: getattr(x.scoring, score_attr, 0))
    
    # First take minimum per track
    track_selected_ids = set()
    for t in [1, 2, 3, 4]:
        min_needed = global_min_per_track.get(t, 0)
        available = [c for c in by_track_global.get(t, []) if id(c) not in track_selected_ids]
        take = available[:min_needed]
        for c in take:
            track_selected_ids.add(id(c))
        ranked_selected.extend(take)
        if len(take) < min_needed:
            print(f"    Warning: Track {t} only has {len(take)} candidates (wanted {min_needed})")
    
    # =========================================================================
    # PHASE 2: ENFORCE MINIMUM ORIGINAL SCAFFOLD (v8.8.1)
    # =========================================================================
    # Ensure original scaffold gets fair representation
    min_original = max(10, n_ranked_slots // 5)  # At least 10 or 20% of ranked
    
    original_in_selected = [c for c in ranked_selected if getattr(c, 'scaffold_type', 'original') == 'original']
    all_original_available = [c for c in ranked if getattr(c, 'scaffold_type', 'original') == 'original']
    print(f"    Original scaffold: {len(original_in_selected)} in selected, {len(all_original_available)} total available, need min {min_original}")
    
    if len(original_in_selected) < min_original:
        # Add more original scaffold candidates
        all_original = [c for c in ranked if getattr(c, 'scaffold_type', 'original') == 'original' and id(c) not in track_selected_ids]
        all_original.sort(key=lambda x: getattr(x.scoring, score_attr, 0))
        need = min_original - len(original_in_selected)
        print(f"    Adding {min(need, len(all_original))} more original scaffold candidates")
        for c in all_original[:need]:
            track_selected_ids.add(id(c))
            ranked_selected.append(c)
    
    print(f"    After global minimums: {len(ranked_selected)} (Track mins + scaffold balance)")
    
    # =========================================================================
    # PHASE 3: FILL REMAINING SLOTS BY BUCKET QUOTAS
    # =========================================================================
    # First pass: take quotas per bucket per track (Tracks 1-4 only)
    for b, tracks in by_bucket_track.items():
        bucket_sel: List = []
        
        for t in [1, 2, 3, 4]:  # Skip Track 0, only ranked tracks
            q = per_bucket_track_quota.get(t, 0)
            if q <= 0:
                continue
            # Only take candidates not already selected
            available = [c for c in tracks.get(t, []) if id(c) not in track_selected_ids]
            take = available[:q]
            bucket_sel.extend(take)
            for c in take:
                track_selected_ids.add(id(c))
        
        # If bucket short, fill from any track by score
        max_per_bucket = bucket_target - per_bucket_track_quota.get(0, 0)  # Exclude Track 0 quota
        if len(bucket_sel) < max_per_bucket:
            leftovers: List = []
            for t, arr in tracks.items():
                if t == 0:
                    continue  # Skip controls
                leftovers.extend([x for x in arr if id(x) not in track_selected_ids])
            leftovers.sort(key=lambda x: getattr(x.scoring, score_attr, 0))
            take = leftovers[:max(0, max_per_bucket - len(bucket_sel))]
            bucket_sel.extend(take)
            for c in take:
                track_selected_ids.add(id(c))
        
        ranked_selected.extend(bucket_sel)
    
    # Deduplicate (in case any got added twice)
    seen_ids = set()
    deduped = []
    for c in ranked_selected:
        if id(c) not in seen_ids:
            seen_ids.add(id(c))
            deduped.append(c)
    ranked_selected = deduped
    
    # Sort ranked by score (best first)
    ranked_selected.sort(key=lambda x: getattr(x.scoring, score_attr, 0))
    
    # Clamp to remaining slots
    if len(ranked_selected) > n_ranked_slots:
        ranked_selected = ranked_selected[:n_ranked_slots]
    
    # If too few, fill from remaining best candidates globally
    if len(ranked_selected) < n_ranked_slots:
        selected_ids = {id(x) for x in ranked_selected}
        remaining = [c for c in ranked if id(c) not in selected_ids]
        remaining.sort(key=lambda x: getattr(x.scoring, score_attr, 0))
        need = n_ranked_slots - len(ranked_selected)
        ranked_selected.extend(remaining[:need])
    
    selected.extend(ranked_selected)
    
    # =========================================================================
    # STEP 4: VERIFY COUNT
    # =========================================================================
    if len(selected) != n_select_total:
        raise RuntimeError(
            f"Selection could not reach requested n_select_total={n_select_total}. "
            f"Got {len(selected)} ({n_controls} controls + {len(ranked_selected)} ranked). "
            f"Consider increasing --n-generate or loosening filters."
        )
    
    print(f"    Final: {n_controls} controls + {len(ranked_selected)} ranked = {len(selected)} total")

    return selected


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
HALLMARK_POSITIONS = {42, 49, 50, 52}
HALLMARK_POSITIONS_LIST = [42, 49, 50, 52]  # stable ordering for logs/UI

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

# Vernier positions (CORRECTED v8.8.3)
# Source: 2_Vernier_Consensus sheet from comprehensive_subfamily_analysis_imgt.xlsx
# These are positions that affect CDR conformation
VERNIER_POSITIONS_FR1 = {2, 4}
VERNIER_POSITIONS_FR2 = {41, 47, 52}   # 52 is both hallmark AND vernier
VERNIER_POSITIONS_FR3 = {66, 67, 68, 69, 71, 76, 78, 82, 87, 89, 91, 94}

# Full vernier set (17 positions) - CORRECTED v8.8.3
ALL_VERNIER_POSITIONS_FULL = VERNIER_POSITIONS_FR1 | VERNIER_POSITIONS_FR2 | VERNIER_POSITIONS_FR3

# For backward compatibility (used in some legacy code paths)
ALL_VERNIER_POSITIONS = ALL_VERNIER_POSITIONS_FULL

# Positions to EXCLUDE from framework sprinkling (hallmarks + all verniers)
EXCLUDE_FROM_FRAMEWORK = HALLMARK_POSITIONS | ALL_VERNIER_POSITIONS_FULL  # 20 positions

# For "Control 0c: ALL verniers" - don't double-apply hallmarks
VERNIERS_MINUS_HALLMARKS = ALL_VERNIER_POSITIONS_FULL - HALLMARK_POSITIONS  # 16 positions

# Positions that are extremely conserved in true VHH datasets and should almost
# never be mutated when building "very likely" candidates.
# (From vernier_conservation_report.md: >=95% conserved)
INVARIANT_VHH_POSITIONS = {4, 41, 47, 76, 89, 94}

# ============================================================
# V8.0 CONSTANTS
# ============================================================

VERSION = "8.9"

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

# Sprinkling probabilities by tier (v8.9: used in Track 4 with caps)
TIER_S_PROB = 0.50   # v8.9: reduced from 0.85, caps are the safety net
TIER_A_PROB = 0.25   # v8.9: reduced from 0.60
TIER_B_PROB = 0.10   # v8.9: reduced from 0.25
TIER_C_PROB = 0.05   # v8.9: only used when desperate for uniqueness

# ============================================================
# V8.9: TRACK CAPS AND PROTECTED POSITIONS
# ============================================================

# Protected positions - motifs and random sprinkling cannot override these
# v8.9: Now includes IMGT2 alongside hallmarks
PROTECTED_POSITIONS = {2, 42, 49, 50, 52}  # IMGT2 + hallmarks

# Track 1: Single vernier probes
MAX_SINGLE_PROBES = 20  # v8.9: Expanded from 12 to include all verniers

# Track 2: Minimal structured perturbations (v8.9)
K_EXTRA_VERNIERS_T2 = 2  # Max extra verniers beyond motif

# Track 3: Motif attribution (v8.9)
K_EXTRA_VERNIERS_T3 = 2  # Default extra verniers (max 3)
K_EXTRA_VERNIERS_T3_MAX = 3  # Hard cap

# Track 4: Production optimizer caps (v8.9)
K_VERNIER_MAX = 8      # Max vernier mutations per candidate
K_S_MAX = 10           # Max Tier S framework mutations
K_A_MAX = 6            # Max Tier A framework mutations
K_B_MAX = 3            # Max Tier B framework mutations  
K_C_MAX = 2            # Max Tier C (only when desperate)
K_CDR_RULES_MAX = 2    # Max CDR-conditional rule mutations

# YQRL progressive framework: positions >80% consensus, applied in PAIRS
# v8.9: Include all positions with consensus ≥80%
YQRL_FRAMEWORK_PAIRS = [
    # Pairs selected based on covariance/co-occurrence or structural proximity
    # Format: [(pos1, aa1, freq1), (pos2, aa2, freq2)]
    [(8, 'G', 0.997), (9, 'G', 0.997)],       # FR1 adjacent, both GG
    [(3, 'Q', 0.996), (5, 'V', 0.995)],       # FR1 N-term region
    [(6, 'E', 0.994), (16, 'G', 0.994)],      # FR1
    [(75, 'R', 0.994), (121, 'G', 0.991)],    # FR3/FR4 anchors
    [(23, 'C', 0.990), (104, 'C', 0.958)],    # Conserved cysteines (structurally coupled)
    [(43, 'R', 0.988), (44, 'Q', 0.955)],     # FR2 RQ motif
    [(124, 'V', 0.988), (126, 'V', 0.988)],   # FR4 VxV motif
    [(18, 'S', 0.984), (19, 'L', 0.981)],     # FR1 SL
    [(21, 'L', 0.981), (22, 'S', 0.945)],     # FR1 LS
    [(102, 'Y', 0.983), (98, 'D', 0.975)],    # FR3 pre-CDR3
    [(119, 'G', 0.978), (127, 'S', 0.968)],   # FR4
    [(13, 'V', 0.969), (11, 'G', 0.963)],     # FR1
    [(46, 'P', 0.961), (125, 'T', 0.959)],    # FR2/FR4
    [(74, 'G', 0.956), (122, 'T', 0.946)],    # FR3/FR4
    [(14, 'Q', 0.942), (7, 'S', 0.937)],      # FR1
    [(51, 'E', 0.937), (26, 'S', 0.935)],     # FR2/FR1
    [(81, 'D', 0.928), (118, 'W', 0.926)],    # FR3/FR4
    [(99, 'T', 0.924), (79, 'S', 0.918)],     # FR3
    [(17, 'G', 0.898), (70, 'S', 0.898)],     # FR1/FR3
    [(48, 'K', 0.889), (53, 'V', 0.889)],     # FR2
    [(90, 'Q', 0.886), (97, 'E', 0.883)],     # FR3
    [(100, 'A', 0.881), (103, 'Y', 0.873)],   # FR3
    [(77, 'T', 0.872), (96, 'P', 0.869)],     # FR3
    [(20, 'R', 0.868), (123, 'Q', 0.865)],    # FR1/FR4
    [(72, 'K', 0.856), (45, 'A', 0.852)],     # FR3/FR2
    [(84, 'K', 0.850), (95, 'K', 0.850)],     # FR3 KxxK
    [(86, 'T', 0.847), (88, 'Y', 0.842)],     # FR3
    [(120, 'Q', 0.807)],                       # Singleton - last position >80%
]

# Flatten for total count (for logging)
YQRL_FW_POSITIONS_OVER_80 = set()
for pair in YQRL_FRAMEWORK_PAIRS:
    for pos, aa, freq in pair:
        YQRL_FW_POSITIONS_OVER_80.add(pos)

# How many progressive framework steps for YQRL
YQRL_MAX_FW_STEPS = len(YQRL_FRAMEWORK_PAIRS)  # One step per pair

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

# ============================================================
# V8.5: HIGH-CONFIDENCE EPISTATIC PAIRS
# ============================================================
# When one position in a pair is mutated, the partner MUST also be mutated
# to maintain structural coupling. Only enforced for high-confidence pairs (>=80%)
#
# Format: {pos1: [(partner_pos, my_aa, partner_aa, confidence), ...]}
# This allows lookup: "if I mutate pos1 to aa1, I must also mutate pos2 to aa2"

EPISTATIC_PAIRS_HIGH_CONF = {
    # Core vernier pairs (>=85% co-occurrence)
    68: [(69, 'A', 'D', 0.85), (71, 'A', 'V', 0.80)],
    69: [(68, 'D', 'A', 0.85), (71, 'D', 'V', 0.82)],
    71: [(69, 'V', 'D', 0.82), (76, 'V', 'F', 0.90), (68, 'V', 'A', 0.80)],
    76: [(71, 'F', 'V', 0.90), (78, 'F', 'I', 0.88)],
    78: [(76, 'I', 'F', 0.88), (89, 'I', 'L', 0.83)],
    89: [(78, 'L', 'I', 0.83), (94, 'L', 'R', 0.85), (91, 'L', 'M', 0.80)],
    91: [(89, 'M', 'L', 0.80), (87, 'M', 'V', 0.79)],
    87: [(91, 'V', 'M', 0.79)],
    94: [(89, 'R', 'L', 0.85)],
}

# Minimum confidence to enforce epistatic coupling
EPISTATIC_MIN_CONFIDENCE = 0.80


def get_epistatic_partners(pos: int, mutant_aa: str, min_conf: float = EPISTATIC_MIN_CONFIDENCE) -> List[Tuple[int, str, float]]:
    """
    Get required partner mutations for an epistatic position.
    
    Args:
        pos: IMGT position being mutated
        mutant_aa: Amino acid we're mutating TO
        min_conf: Minimum confidence to enforce coupling
        
    Returns:
        List of (partner_pos, partner_aa, confidence) that MUST be applied
    """
    partners = []
    if pos in EPISTATIC_PAIRS_HIGH_CONF:
        for partner_pos, my_aa, partner_aa, conf in EPISTATIC_PAIRS_HIGH_CONF[pos]:
            if my_aa == mutant_aa and conf >= min_conf:
                partners.append((partner_pos, partner_aa, conf))
    return partners


def apply_with_epistatic_enforcement(
    muts: List, 
    used_positions: Set[int],
    new_mut,
    original_positions: Dict[int, str],
    cdr3_has_cysteine: bool = False,
    min_conf: float = EPISTATIC_MIN_CONFIDENCE
) -> bool:
    """
    Apply a mutation WITH epistatic partner enforcement.
    
    If the mutation triggers epistatic coupling, partner mutations are also added.
    Returns True if mutation (and partners) were successfully applied, False if blocked.
    
    Safety: Will not add cysteine unless CDR3 has cysteine.
    """
    pos = new_mut.imgt_num
    mutant_aa = new_mut.mutant
    
    # Check epistatic partners
    partners = get_epistatic_partners(pos, mutant_aa, min_conf)
    
    # Verify all partners can be added
    partner_muts = []
    for partner_pos, partner_aa, conf in partners:
        if partner_pos in used_positions:
            # Partner already mutated - check if compatible
            existing = [m for m in muts if m.imgt_num == partner_pos]
            if existing and existing[0].mutant != partner_aa:
                # Conflict! Cannot apply this mutation
                return False
            # Already correct, no action needed
            continue
        
        # Safety: no unpaired cysteine
        if partner_aa == 'C' and not cdr3_has_cysteine:
            return False  # Block the whole mutation
        
        orig_aa = original_positions.get(partner_pos, '')
        if orig_aa and orig_aa != partner_aa:
            partner_muts.append(Mutation(
                position=f"IMGT{partner_pos}",
                imgt_num=partner_pos,
                original=orig_aa,
                mutant=partner_aa,
                source=f"epistatic_partner_of_{pos}",
                confidence=conf
            ))
    
    # Apply primary mutation
    muts.append(new_mut)
    used_positions.add(pos)
    
    # Apply partner mutations
    for pm in partner_muts:
        muts.append(pm)
        used_positions.add(pm.imgt_num)
    
    return True

# ============================================================
# V8.4: CDR-CONDITIONAL FRAMEWORK MUTATIONS
# ============================================================

# Positions that were INTENTIONALLY humanized in Vincke Universal scaffold
# These should be RESPECTED in Tracks 0-2, only OVERRIDDEN in Track 3+
HUMANIZED_POSITIONS = {87, 95, 96, 123}

# R50 consensus mutations (from R50 VHH subfamily analysis)
# These are applied in Track 4 to test if "VHH-ifying" helps
R50_CONSENSUS_MUTATIONS = {
    48: ('Q', 'K', 0.87),   # Q→K in R50 VHH consensus
    66: ('Y', 'S', 0.59),   # Y→S in DERL-type R50 (optional, low confidence)
}

# YQRL is the most constrained hallmark - only controls, no optimization
RIGID_HALLMARKS = {'YQRL'}

# Framework positions (for CDR-conditional rule matching)
FRAMEWORK_POSITIONS = set(range(1, 27)) | set(range(40, 56)) | set(range(66, 105)) | set(range(118, 129))


class CDRConditionalRules:
    """
    Load and match CDR-conditional framework rules from analysis_rules_v7_all_positions.json.
    
    Rules are keyed by CDR3 features like:
    - cdr3[-1]=S (C-terminal residue)
    - cdr3_len=short
    - cdr3_charge=negative
    - cdr3_cys=multi_cys
    """
    
    def __init__(self, rules_path: str = None):
        self.rules = []
        self.rules_by_condition = defaultdict(list)
        
        if rules_path:
            self.load_rules(rules_path)
    
    def load_rules(self, rules_path: str):
        """Load rules from JSON file."""
        try:
            with open(rules_path) as f:
                all_rules = json.load(f)
            
            # Filter to framework positions only
            for rule in all_rules:
                pos_str = rule.get('position', '')
                if isinstance(pos_str, str) and pos_str.startswith('IMGT'):
                    try:
                        pos_num = int(pos_str.replace('IMGT', ''))
                        if pos_num in FRAMEWORK_POSITIONS:
                            rule['position_num'] = pos_num
                            self.rules.append(rule)
                            
                            # Index by condition
                            cond = rule.get('condition', '')
                            self.rules_by_condition[cond].append(rule)
                    except ValueError:
                        continue
            
            print(f"  Loaded {len(self.rules)} CDR-conditional framework rules")
            
        except Exception as e:
            print(f"  Warning: Could not load CDR-conditional rules: {e}")
    
    def get_matching_rules(self, cdr3_features: Dict[str, str], 
                           min_confidence: float = 0.70,
                           skip_cysteine: bool = True,
                           cdr3_has_cysteine: bool = False) -> Dict[int, Dict]:
        """
        Get framework mutations that apply to this CDR3.
        
        Returns: {position: {'aa': target_aa, 'confidence': conf, 'condition': cond}}
        """
        # Build list of conditions this CDR3 matches
        matching_conditions = []
        
        # C-terminal residue
        last_aa = cdr3_features.get('cdr3[-1]', '')
        if last_aa:
            matching_conditions.append(f"cdr3[-1]={last_aa}")
        
        # Length class
        len_class = cdr3_features.get('cdr3_len', '')
        if len_class == 'short':
            matching_conditions.append("cdr3_len=short")
        
        # Charge
        charge = cdr3_features.get('cdr3_charge', '')
        if charge == 'negative':
            matching_conditions.append("cdr3_charge=negative")
        elif charge == 'positive':
            matching_conditions.append("cdr3_charge=positive")
        
        # Cysteine content
        cys = cdr3_features.get('cdr3_cys', '')
        if cys == 'multiple':
            matching_conditions.append("cdr3_cys=multi_cys")
        
        # Collect matching rules
        best_by_position = {}
        
        for cond in matching_conditions:
            for rule in self.rules_by_condition.get(cond, []):
                pos = rule.get('position_num')
                aa = rule.get('suggested_aa', '')
                conf = rule.get('confidence', 0)
                
                if conf < min_confidence:
                    continue
                
                # CRITICAL: Skip cysteine mutations unless CDR3 can pair with it
                if skip_cysteine and aa == 'C' and not cdr3_has_cysteine:
                    continue
                
                # Keep highest confidence rule per position
                if pos not in best_by_position or conf > best_by_position[pos]['confidence']:
                    best_by_position[pos] = {
                        'aa': aa,
                        'confidence': conf,
                        'condition': cond,
                        'support': rule.get('support', 0)
                    }
        
        return best_by_position
    
    def get_rules_by_confidence_tier(self, cdr3_features: Dict[str, str],
                                      universal_positions: Dict[int, str],
                                      cdr3_has_cysteine: bool = False) -> Dict[str, List]:
        """
        Get CDR-conditional rules organized by confidence tier.
        
        Returns:
        {
            'high': [(pos, univ_aa, rule_aa, conf), ...],   # >=90%
            'medium': [(pos, univ_aa, rule_aa, conf), ...], # 80-90%
            'medium_humanized': [...],                       # 80-90% at humanized positions
        }
        """
        all_rules = self.get_matching_rules(cdr3_features, 
                                            min_confidence=0.80,
                                            skip_cysteine=True,
                                            cdr3_has_cysteine=cdr3_has_cysteine)
        
        tiers = {
            'high': [],
            'medium': [],
            'medium_humanized': [],
        }
        
        for pos, rule_data in all_rules.items():
            univ_aa = universal_positions.get(pos, '')
            rule_aa = rule_data['aa']
            conf = rule_data['confidence']
            
            # Skip if Universal already has this AA
            if univ_aa == rule_aa:
                continue
            
            entry = (pos, univ_aa, rule_aa, conf, rule_data['condition'])
            
            if conf >= 0.90:
                tiers['high'].append(entry)
            elif pos in HUMANIZED_POSITIONS:
                tiers['medium_humanized'].append(entry)
            else:
                tiers['medium'].append(entry)
        
        # Sort each tier by confidence descending
        for tier in tiers.values():
            tier.sort(key=lambda x: -x[3])
        
        return tiers


# ============================================================
# V8.5: UNIVERSAL FRAMEWORK CONSENSUS TIERS
# ============================================================
# Cross-hallmark consensus for framework positions (excluding already-handled)
# Derived from 12M VHH database - positions where >=85% of hallmarks agree

FRAMEWORK_CONSENSUS_TIERS = {
    # TIER S: >=95% frequency - Near universal, safe to apply early (Track 2)
    'S': {
        3: ('Q', 0.996),   # FR1
        4: ('L', 0.996),   # FR1
        5: ('V', 0.995),   # FR1
        6: ('E', 0.994),   # FR1
        8: ('G', 0.997),   # FR1
        9: ('G', 0.997),   # FR1
        11: ('G', 0.963),  # FR1
        13: ('V', 0.969),  # FR1
        16: ('G', 0.994),  # FR1
        18: ('S', 0.984),  # FR1
        19: ('L', 0.981),  # FR1
        21: ('L', 0.981),  # FR1
        23: ('C', 0.990),  # FR1
        41: ('W', 0.981),  # FR2
        43: ('R', 0.988),  # FR2
        44: ('Q', 0.955),  # FR2
        46: ('P', 0.961),  # FR2
        47: ('G', 0.975),  # FR2
        74: ('G', 0.956),  # FR3
        75: ('R', 0.994),  # FR3
        98: ('D', 0.975),  # FR3
        102: ('Y', 0.983), # FR3
        104: ('C', 0.958), # FR3
        119: ('G', 0.978), # FR4
        121: ('G', 0.991), # FR4
        124: ('V', 0.988), # FR4
        125: ('T', 0.959), # FR4
        126: ('V', 0.988), # FR4
        127: ('S', 0.968), # FR4
    },
    # TIER A: 90-95% frequency - Very high confidence (Track 3)
    'A': {
        7: ('S', 0.937),   # FR1
        14: ('Q', 0.942),  # FR1
        22: ('S', 0.945),  # FR1
        26: ('S', 0.935),  # FR1
        51: ('E', 0.937),  # FR2
        67: ('Y', 0.923),  # FR3
        79: ('S', 0.918),  # FR3
        81: ('D', 0.928),  # FR3
        99: ('T', 0.924),  # FR3
        118: ('W', 0.926), # FR4
        122: ('T', 0.946), # FR4
    },
    # TIER B: 85-90% frequency - High confidence (Track 4 base)
    'B': {
        17: ('G', 0.898),  # FR1
        20: ('R', 0.868),  # FR1
        45: ('A', 0.852),  # FR2
        48: ('K', 0.889),  # FR2
        53: ('V', 0.889),  # FR2
        70: ('S', 0.898),  # FR3
        72: ('K', 0.856),  # FR3
        77: ('T', 0.872),  # FR3
        84: ('K', 0.850),  # FR3
        90: ('Q', 0.886),  # FR3
        95: ('K', 0.850),  # FR3
        96: ('P', 0.869),  # FR3
        97: ('E', 0.883),  # FR3
        100: ('A', 0.881), # FR3
        103: ('Y', 0.873), # FR3
        123: ('Q', 0.865), # FR4
    },
    # TIER C: 75-85% frequency - Medium confidence (Track 4 expansive)
    'C': {
        86: ('T', 0.847),  # FR3
        88: ('Y', 0.842),  # FR3
        120: ('Q', 0.807), # FR4
    },
}

# All framework consensus positions (for easy lookup)
ALL_FRAMEWORK_CONSENSUS_POSITIONS = set()
for tier_data in FRAMEWORK_CONSENSUS_TIERS.values():
    ALL_FRAMEWORK_CONSENSUS_POSITIONS.update(tier_data.keys())


# ============================================================
# V8.8.3: YQRL PROGRESSIVE FRAMEWORK SPRINKLING
# ============================================================

# YQRL-specific vernier consensus (from 2_Vernier_Consensus sheet)
# Format: {position: (consensus_aa, frequency)}
YQRL_VERNIER_MUTATIONS = {
    2:  ('V', 0.880),   # FR1 - universal I→V
    4:  ('L', 0.998),   # FR1
    41: ('W', 0.994),   # FR2
    47: ('G', 0.985),   # FR2
    52: ('L', 1.000),   # FR2 - also hallmark (already applied)
    66: ('N', 0.710),   # FR3
    67: ('Y', 0.960),   # FR3
    68: ('A', 0.855),   # FR3
    69: ('D', 0.907),   # FR3
    71: ('V', 0.935),   # FR3
    76: ('F', 0.974),   # FR3
    78: ('I', 0.937),   # FR3
    82: ('N', 0.870),   # FR3
    87: ('V', 0.767),   # FR3
    89: ('L', 0.982),   # FR3
    91: ('M', 0.967),   # FR3
    94: ('L', 0.974),   # FR3
}

# Verniers to apply in YQRL Control 0c (excludes 52 which is hallmark)
YQRL_VERNIERS_TO_APPLY = {k: v for k, v in YQRL_VERNIER_MUTATIONS.items() if k not in HALLMARK_POSITIONS}

# Framework positions ranked by consensus frequency for YQRL progressive sprinkling
# EXCLUDES: hallmarks (42, 49, 50, 52) and ALL verniers (2, 4, 41, 47, 52, 66-94)
# Use for progressive FW sprinkling after Control 0c
YQRL_FRAMEWORK_MUTATIONS_RANKED = [
    (8, 'G', 0.9970),    # FR1, Tier S
    (9, 'G', 0.9970),    # FR1, Tier S
    (3, 'Q', 0.9960),    # FR1, Tier S
    (5, 'V', 0.9950),    # FR1, Tier S
    (6, 'E', 0.9940),    # FR1, Tier S
    (16, 'G', 0.9940),   # FR1, Tier S
    (75, 'R', 0.9940),   # FR3, Tier S
    (121, 'G', 0.9910),  # FR4, Tier S
    (23, 'C', 0.9900),   # FR1, Tier S
    (43, 'R', 0.9880),   # FR2, Tier S
    (124, 'V', 0.9880),  # FR4, Tier S
    (126, 'V', 0.9880),  # FR4, Tier S
    (18, 'S', 0.9840),   # FR1, Tier S
    (102, 'Y', 0.9830),  # FR3, Tier S
    (19, 'L', 0.9810),   # FR1, Tier S
    (21, 'L', 0.9810),   # FR1, Tier S
    (119, 'G', 0.9780),  # FR4, Tier S
    (98, 'D', 0.9750),   # FR3, Tier S
    (13, 'V', 0.9690),   # FR1, Tier S
    (127, 'S', 0.9680),  # FR4, Tier S
    (11, 'G', 0.9630),   # FR1, Tier S
    (46, 'P', 0.9610),   # FR2, Tier S
    (125, 'T', 0.9590),  # FR4, Tier S
    (104, 'C', 0.9580),  # FR3, Tier S
    (74, 'G', 0.9560),   # FR3, Tier S
    (44, 'Q', 0.9550),   # FR2, Tier S
    (122, 'T', 0.9460),  # FR4, Tier A
    (22, 'S', 0.9450),   # FR1, Tier A
    (14, 'Q', 0.9420),   # FR1, Tier A
    (7, 'S', 0.9370),    # FR1, Tier A
    (51, 'E', 0.9370),   # FR2, Tier A
    (26, 'S', 0.9350),   # FR1, Tier A
    (81, 'D', 0.9280),   # FR3, Tier A
    (118, 'W', 0.9260),  # FR4, Tier A
    (99, 'T', 0.9240),   # FR3, Tier A
    (79, 'S', 0.9180),   # FR3, Tier A
    (17, 'G', 0.8980),   # FR1, Tier B
    (70, 'S', 0.8980),   # FR3, Tier B
    (48, 'K', 0.8890),   # FR2, Tier B
    (53, 'V', 0.8890),   # FR2, Tier B
    (90, 'Q', 0.8860),   # FR3, Tier B
    (97, 'E', 0.8830),   # FR3, Tier B
    (100, 'A', 0.8810),  # FR3, Tier B
    (103, 'Y', 0.8730),  # FR3, Tier B
    (77, 'T', 0.8720),   # FR3, Tier B
    (96, 'P', 0.8690),   # FR3, Tier B
    (20, 'R', 0.8680),   # FR1, Tier B
    (123, 'Q', 0.8650),  # FR4, Tier B
    (72, 'K', 0.8560),   # FR3, Tier B
    (45, 'A', 0.8520),   # FR2, Tier B
    (84, 'K', 0.8500),   # FR3, Tier B
    (95, 'K', 0.8500),   # FR3, Tier B
    (86, 'T', 0.8470),   # FR3, Tier C
    (88, 'Y', 0.8420),   # FR3, Tier C
    (120, 'Q', 0.8070),  # FR4, Tier C
]

# Guardrail: Verify no excluded positions in framework list
_bad_positions = [pos for (pos, aa, freq) in YQRL_FRAMEWORK_MUTATIONS_RANKED if pos in EXCLUDE_FROM_FRAMEWORK]
assert not _bad_positions, f"YQRL_FRAMEWORK_MUTATIONS_RANKED contains excluded positions: {_bad_positions}"

# How many progressive framework steps to generate for YQRL
YQRL_MAX_FW_STEPS = 10  # Generate FW+1 through FW+10


def get_yqrl_control_mutations(control_level: str, original_positions: dict) -> list:
    """
    Get mutations to apply for each YQRL control level.
    
    Args:
        control_level: '0a_graft', '0b_hallmarks', '0c_verniers', or 'fw+N' where N is 1-55
        original_positions: Dict of {imgt_pos: amino_acid} for the input sequence
    
    Returns:
        List of Mutation objects to apply
    """
    mutations = []
    
    if control_level == '0a_graft':
        return []  # No mutations for pure graft
    
    # Hallmarks (42→Y, 49→Q, 50→R, 52→L) + IMGT2 (→V)
    hallmark_targets = {42: 'Y', 49: 'Q', 50: 'R', 52: 'L', 2: 'V'}
    
    for pos, target_aa in hallmark_targets.items():
        orig_aa = original_positions.get(pos, '')
        if orig_aa and orig_aa != target_aa:
            mutations.append(Mutation(
                position=f"IMGT{pos}",
                imgt_num=pos,
                original=orig_aa,
                mutant=target_aa,
                source="yqrl_hallmark" if pos in HALLMARK_POSITIONS else "yqrl_imgt2",
                confidence=1.0 if pos in HALLMARK_POSITIONS else 0.88
            ))
    
    if control_level == '0b_hallmarks':
        return mutations
    
    # Add all verniers (except 52 which is already in hallmarks)
    for pos, (target_aa, freq) in YQRL_VERNIERS_TO_APPLY.items():
        orig_aa = original_positions.get(pos, '')
        if orig_aa and orig_aa != target_aa:
            mutations.append(Mutation(
                position=f"IMGT{pos}",
                imgt_num=pos,
                original=orig_aa,
                mutant=target_aa,
                source="yqrl_vernier",
                confidence=freq
            ))
    
    if control_level == '0c_verniers':
        return mutations
    
    # Progressive framework sprinkling
    if control_level.startswith('fw+'):
        n = int(control_level[3:])
        for pos, target_aa, freq in YQRL_FRAMEWORK_MUTATIONS_RANKED[:n]:
            orig_aa = original_positions.get(pos, '')
            if orig_aa and orig_aa != target_aa:
                mutations.append(Mutation(
                    position=f"IMGT{pos}",
                    imgt_num=pos,
                    original=orig_aa,
                    mutant=target_aa,
                    source=f"yqrl_framework_step{n}",
                    confidence=freq
                ))
        return mutations
    
    raise ValueError(f"Unknown YQRL control level: {control_level}")


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
    track: str = ""                   # v8.5: Explicit track label ("Track 0", "Track 1", etc.)
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

        # v8.6: Build motif pools from archetype consensus signatures
        # Track 2 uses motifs of length 2–3; Track 3 uses motifs of length 3–5
        self.motif_pool_2_3 = build_motif_pool_from_archetypes(self.archetypes, min_len=2, max_len=3)
        self.motif_pool_3_5 = build_motif_pool_from_archetypes(self.archetypes, min_len=3, max_len=5)

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
    
    def __init__(self, rules_file: str, archetypes_file: str, hallmark_db_path: str = None,
                 sprinkle_temperature: float = 1.0):
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

        # v8.4: Load CDR-conditional rules
        # These rules are in the same file but have different condition patterns
        self.cdr_conditional_rules = CDRConditionalRules(rules_file)
        
        # Build scaffold position dict
        self.scaffold_positions = number_scaffold_sequence(
            UNIVERSAL_SCAFFOLD['FR1'], '',  # No CDR in scaffold
            UNIVERSAL_SCAFFOLD['FR2'], '',
            UNIVERSAL_SCAFFOLD['FR3'], '',
            UNIVERSAL_SCAFFOLD['FR4']
        )
        
        # v8.6: Temperature for Bernoulli sprinkling
        self.sprinkle_temperature = sprinkle_temperature
        
        # v8.6: Build motif pool from archetypes for Track 2-3
        self.motif_pool = build_motif_pool_from_archetypes(self.archetypes, min_len=2, max_len=5)
        
        # v8.8: Build SEPARATE pools for Track 2 (2-3 positions) and Track 3 (3-5 positions)
        self.motif_pool_2_3 = build_motif_pool_from_archetypes(self.archetypes, min_len=2, max_len=3)
        self.motif_pool_3_5 = build_motif_pool_from_archetypes(self.archetypes, min_len=3, max_len=5)
        
        print(f"  {len(self.compensation_rules)} compensation rules")
        print(f"  {len(self.triplet_rules)} triplet rules")
        print(f"  {len(self.archetypes)} family archetypes")
        print(f"  v8.4: CDR-conditional rules loaded for Universal tracks")
        print(f"  v8.6: {len(self.motif_pool)} archetype motifs (2-5 positions)")
        print(f"  v8.8: {len(self.motif_pool_2_3)} motifs for Track 2, {len(self.motif_pool_3_5)} motifs for Track 3")
        print(f"  v8.6: Sprinkling temperature: {sprinkle_temperature}")
    
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
        
        # v8.0: Always generates BOTH original framework and Universal scaffold variants
        mode = "BOTH"
        print(f"  Mode: {mode} (original framework + Universal scaffold)")

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
        
        # =================================================================
        # v8.8.2: SEPARATE UNIVERSAL AND ORIGINAL GENERATION
        # =================================================================
        # Universal scaffold: ONE fixed hallmark (FERG default) + FW variation
        # Original scaffold: Distributed across ALL target hallmarks
        # =================================================================
        
        # Pick single hallmark for Universal scaffold
        # Priority: FERG > first F-type hallmark > first hallmark in pool
        universal_hallmark = 'FERG'
        if hallmark_pool:
            if 'FERG' in hallmark_pool:
                universal_hallmark = 'FERG'
            else:
                f_type = [h for h in hallmark_pool if h.startswith('F')]
                if f_type:
                    universal_hallmark = f_type[0]
                else:
                    universal_hallmark = hallmark_pool[0]
        
        print(f"\n  Universal scaffold hallmark: {universal_hallmark} (fixed)")
        
        # Calculate allocation
        n_remaining = n_candidates - 1 - len(universal_controls)
        univ_frac = 0.30  # 30% to Universal
        n_univ_total = int(round(n_remaining * univ_frac))
        n_orig_total = n_remaining - n_univ_total
        
        print(f"  Allocation: {n_univ_total} Universal, {n_orig_total} Original")
        
        # =================================================================
        # GENERATE UNIVERSAL SCAFFOLD CANDIDATES (single hallmark)
        # =================================================================
        print(f"\n  Generating Universal scaffold candidates ({universal_hallmark})...")
        
        # Use first family for Universal (family affects rules lookup)
        univ_family = target_families[0] if target_families else 'F_C2'
        rules = self._get_rules(universal_hallmark, univ_family, cdr3_features)
        archetype = self._get_archetype(univ_family)
        triplet_rules = self._get_triplet_rules(univ_family, cdrs.cdr3)
        
        univ_cands, generation_order = self._generate_universal(
            cdrs, univ_family, rules, archetype, triplet_rules, n_univ_total,
            universal_hallmark, generation_order, cdr3_features
        )
        candidates.extend(univ_cands)
        print(f"    Generated {len(univ_cands)} Universal candidates (all {universal_hallmark})")
        
        # =================================================================
        # GENERATE ORIGINAL SCAFFOLD CANDIDATES (distributed hallmarks)
        # =================================================================
        print(f"\n  Generating Original scaffold candidates (distributed across hallmarks)...")
        
        # Build generation plan for Original only
        n_per_family = n_orig_total // len(target_families)
        generation_plan = []  # List of (family, hallmark, n_candidates)
        
        for i, family in enumerate(target_families):
            n_for_family = n_per_family
            if i == len(target_families) - 1:
                n_for_family = n_orig_total - sum(g[2] for g in generation_plan)
            
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
        
        # Print generation plan (Original scaffold only)
        print(f"  Original scaffold generation plan ({len(generation_plan)} hallmark-family combos):")
        for family, hm, n in generation_plan:
            print(f"    {family} + {hm}: {n} candidates")
        
        # Track which hallmarks have already had controls generated (for deduplication)
        # Controls are hallmark-specific, not family-specific, so we only generate them once
        hallmarks_with_controls = set()
        
        # Execute generation plan (Original scaffold only - Universal already generated above)
        for family, fam_hallmarks, n_for_combo in generation_plan:
            rules = self._get_rules(fam_hallmarks, family, cdr3_features)
            archetype = self._get_archetype(family)
            triplet_rules = self._get_triplet_rules(family, cdrs.cdr3)
            
            print(f"  Original {family}/{fam_hallmarks}: {n_for_combo} candidates, {len(rules)} rules")
            
            # v8.8.2: Generate Original framework variants only
            orig_generated = 0
            if n_for_combo > 0 and cdrs.fr1:
                new_cands, generation_order = self._generate_original(
                    cdrs, family, rules, archetype, triplet_rules, n_for_combo,
                    fam_hallmarks, generation_order, original_positions,
                    hallmarks_with_controls  # Pass the deduplication set
                )
                candidates.extend(new_cands)
                orig_generated = len(new_cands)
            elif n_for_combo > 0 and not cdrs.fr1:
                print(f"    Warning: Skipping original scaffold - FR1 not available")
            
            print(f"    Generated: {orig_generated} original")
        
        # v8.8.2: Show scaffold breakdown before returning
        orig_total = sum(1 for c in candidates if getattr(c, 'scaffold_type', 'original') == 'original')
        univ_total = sum(1 for c in candidates if getattr(c, 'scaffold_type', 'original') == 'universal')
        print(f"\n  Total generated: {orig_total} original, {univ_total} universal")
        
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
        Generate Universal scaffold Track 0 controls (v8.9).
        
        v8.9 STRUCTURE:
        0a) Complete CDR graft into Universal (0 mutations)
        0b) Hallmarks + IMGT2 (5 mutations)
        0c) Hallmarks + IMGT2 + ALL verniers (~20 mutations)
        
        Universal framework = consensus VHH scaffold
        Base sequence: Universal framework with input CDRs grafted in
        """
        candidates = []
        used = set()
        
        # Build Universal framework positions with CDRs grafted in
        universal_positions = self._build_universal_positions_with_cdrs(cdrs)
        
        # --- Control 0a: Complete CDR graft (0 mutations) ---
        # Base: Universal framework + CDRs, no changes to framework
        # Tests: Do CDRs work in a known-good VHH framework without any conversion?
        generation_order += 1
        candidate = self._build_universal_candidate(
            f"Universal_graft_{generation_order}",
            cdrs, family, [], "universal_graft", generation_order,
            "Track0: Universal CDR graft (0 mut) - IMMUTABLE", family
        )
        candidate.design_track = "minimal_hallmark"
        candidate.track = "Track 0a (Universal)"
        candidate.ranking_exempt = True
        candidate.track_info = "Universal_graft_0mut"
        candidate.is_immutable = True
        candidate.track_subtype = "graft"
        candidate.scaffold_type = "universal"
        candidates.append(candidate)
        
        # Prepare IMGT2 mutation (if Universal doesn't already have V)
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
        
        # Select top hallmarks to test (max 4)
        test_hallmarks = hallmark_pool[:4] if len(hallmark_pool) > 4 else hallmark_pool
        
        # Always include YQRL if available
        if 'YQRL' in hallmark_pool and 'YQRL' not in test_hallmarks:
            test_hallmarks = ['YQRL'] + test_hallmarks[:3]
        
        for hallmark in test_hallmarks:
            # Get hallmark mutations
            hallmark_muts = self._get_hallmark_mutations_for_universal(hallmark, universal_positions)
            
            # --- Control 0b: Universal + Hallmarks + IMGT2 (5 mutations) ---
            muts = list(imgt2_mut) + hallmark_muts
            key = ("universal", hallmark, "hallmarks_imgt2")
            if key not in used:
                used.add(key)
                generation_order += 1
                candidate = self._build_universal_candidate(
                    f"Universal_{hallmark}_5mut_{generation_order}",
                    cdrs, family, muts, f"universal_{hallmark}", generation_order,
                    f"Track0: Universal + {hallmark} + IMGT2 (5 mut) - IMMUTABLE", family
                )
                candidate.design_track = "minimal_hallmark"
                candidate.track = "Track 0b (Universal)"
                candidate.ranking_exempt = True
                candidate.track_info = f"Universal_{hallmark}_5mut"
                candidate.is_immutable = True
                candidate.track_subtype = "hallmarks_imgt2"
                candidate.scaffold_type = "universal"
                candidates.append(candidate)
            
            # --- Control 0c: Universal + Hallmarks + IMGT2 + ALL verniers ---
            # v8.9: Include for all hallmarks tested (not just high-consensus)
            vernier_muts = self._get_vernier_mutations_for_universal(hallmark, universal_positions)
            all_muts = list(imgt2_mut) + hallmark_muts + vernier_muts
            
            key = ("universal", hallmark, "full_vernier")
            if key not in used:
                used.add(key)
                generation_order += 1
                candidate = self._build_universal_candidate(
                    f"Universal_{hallmark}_full_{generation_order}",
                    cdrs, family, all_muts, f"universal_{hallmark}_full", generation_order,
                    f"Track0: Universal + {hallmark} + IMGT2 + verniers - IMMUTABLE", family
                )
                candidate.design_track = "minimal_hallmark"
                candidate.track = "Track 0c (Universal)"
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
                             generation_order: int,
                             cdr3_features: Dict[str, str] = None) -> Tuple[List[VHHCandidate], int]:
        """
        Generate using universal scaffold with CDR-CONDITIONAL framework tracks (v8.4).
        
        Track structure for Universal:
        - Track 0: Control (hallmarks only - CDR graft into Universal)
        - Track 1: High confidence CDR-conditional rules (>=90%)
        - Track 2: High + Medium confidence (respect humanization)
        - Track 3: High + Medium (override humanization at 87, 95, 96, 123)
        - Track 4: CDR-conditional + R50 consensus (positions 48, 66)
        
        For YQRL (rigid hallmark): Only Track 0 controls, skip optimization tracks.
        """
        import random
        rng = random.Random(42 + generation_order)
        
        candidates = []
        used = set()
        
        # Check if this is a rigid hallmark (YQRL)
        is_rigid = target_hallmarks in RIGID_HALLMARKS
        
        # Get CDR3 features if not provided
        if cdr3_features is None:
            cdr3_features = get_cdr3_features(cdrs.cdr3)
        
        # Check if CDR3 has cysteine (for safety check)
        cdr3_has_cysteine = cdrs.cdr3.count('C') > 0
        
        hallmark_set = set(HALLMARK_POSITIONS)
        
        # v8.8.1: Get vernier profile early for Track 1 probes
        vernier_profile = {}
        if self.hallmark_db:
            vernier_profile = self.hallmark_db.vernier_profile(target_hallmarks)
        
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
        # GET CDR-CONDITIONAL RULES (v8.4)
        # =====================================================================
        cdr_rules_by_tier = {'high': [], 'medium': [], 'medium_humanized': []}
        
        if hasattr(self, 'cdr_conditional_rules') and self.cdr_conditional_rules:
            cdr_rules_by_tier = self.cdr_conditional_rules.get_rules_by_confidence_tier(
                cdr3_features, 
                self.scaffold_positions,
                cdr3_has_cysteine=cdr3_has_cysteine
            )
            
            print(f"    CDR-conditional rules for {target_hallmarks}:")
            print(f"      High (>=90%): {len(cdr_rules_by_tier['high'])} positions")
            print(f"      Medium (80-90%): {len(cdr_rules_by_tier['medium'])} positions")
            print(f"      Medium humanized: {len(cdr_rules_by_tier['medium_humanized'])} positions")
        
        # =====================================================================
        # TRACK 0: CONTROL (hallmarks only - immutable)
        # =====================================================================
        if hallmark_muts:
            key = tuple(sorted([str(m) for m in hallmark_muts]))
            if key not in used:
                used.add(key)
                generation_order += 1
                candidate = self._build_universal_candidate(
                    f"Univ_{target_hallmarks}_T0_control_{generation_order}",
                    cdrs, family, hallmark_muts, "universal_hallmarks", generation_order,
                    f"Track0: Universal + {target_hallmarks} (CDR graft only)", family
                )
                candidate.design_track = "minimal_hallmark"
                candidate.track = "Track 0 (Universal)"
                candidate.ranking_exempt = True
                candidate.track_info = f"Universal_{target_hallmarks}_T0"
                candidate.is_immutable = True
                candidate.track_subtype = "cdr_graft"
                candidates.append(candidate)
        
        # For YQRL (rigid hallmark): STOP HERE - only controls
        if is_rigid:
            print(f"    {target_hallmarks} is rigid hallmark - skipping optimization tracks")
            return candidates, generation_order
        
        # =====================================================================
        # TRACK 1: VERNIER PROBES FOR UNIVERSAL (v8.8.1)
        # =====================================================================
        # Generate single vernier probes using consensus from hallmark profile
        # =====================================================================
        
        MAX_UNIV_PROBES = 8  # v8.8.1: Add vernier probes to Universal
        
        probe_verniers = []
        for pos in [67, 68, 69, 71, 76, 78, 82, 87, 89, 91, 94]:
            profile_data = vernier_profile.get(pos, {})
            consensus_aa = profile_data.get('cons', '')
            consensus_freq = float(profile_data.get('freq', 0))
            scaffold_aa = self.scaffold_positions.get(pos, '')
            
            if consensus_aa and scaffold_aa and consensus_aa != scaffold_aa and consensus_freq >= 0.50:
                probe_verniers.append((pos, scaffold_aa, consensus_aa, consensus_freq))
        
        # Sort by consensus frequency
        probe_verniers.sort(key=lambda x: -x[3])
        
        univ_probe_count = 0
        for pos, orig, cons, freq in probe_verniers[:MAX_UNIV_PROBES]:
            if univ_probe_count >= MAX_UNIV_PROBES:
                break
            
            muts = list(hallmark_muts)
            muts.append(Mutation(
                position=f"IMGT{pos}",
                imgt_num=pos,
                original=orig,
                mutant=cons,
                source=f"vernier_probe_{target_hallmarks}",
                confidence=freq
            ))
            
            key = tuple(sorted([str(m) for m in muts]))
            if key not in used:
                used.add(key)
                generation_order += 1
                univ_probe_count += 1
                candidate = self._build_universal_candidate(
                    f"Univ_{target_hallmarks}_probe_IMGT{pos}_{generation_order}",
                    cdrs, family, muts, "single_vernier", generation_order,
                    f"Track1: Universal vernier probe IMGT{pos}", family
                )
                candidate.design_track = "single_vernier"
                candidate.track = "Track 1 (Universal)"
                candidate.ranking_exempt = False
                candidate.track_info = f"IMGT{pos}_{orig}>{cons}"
                candidates.append(candidate)
        
        # =====================================================================
        # TRACK 1: HIGH CONFIDENCE CDR-CONDITIONAL (>=90%)
        # =====================================================================
        if cdr_rules_by_tier['high']:
            muts = list(hallmark_muts)
            used_positions = {m.imgt_num for m in muts}
            
            for pos, univ_aa, rule_aa, conf, cond in cdr_rules_by_tier['high']:
                if pos not in used_positions:
                    # SAFETY: Never add unpaired cysteine
                    if rule_aa == 'C' and not cdr3_has_cysteine:
                        continue
                    muts.append(Mutation(
                        position=f"IMGT{pos}",
                        imgt_num=pos,
                        original=univ_aa,
                        mutant=rule_aa,
                        source=f"cdr_conditional_{cond}",
                        confidence=conf
                    ))
                    used_positions.add(pos)
            
            if len(muts) > len(hallmark_muts):
                key = tuple(sorted([str(m) for m in muts]))
                if key not in used:
                    used.add(key)
                    generation_order += 1
                    candidate = self._build_universal_candidate(
                        f"Univ_{target_hallmarks}_T1_high_{generation_order}",
                        cdrs, family, muts, "cdr_conditional_high", generation_order,
                        f"Track1: Universal + high conf CDR-conditional ({len(muts)-len(hallmark_muts)} extra)", family
                    )
                    candidate.design_track = "single_vernier"
                    candidate.track = "Track 1 (Universal)"
                    candidate.ranking_exempt = False  # v8.6: Track 1 is ranked
                    candidate.track_info = f"Universal_{target_hallmarks}_T1_high"
                    candidate.track_subtype = "cdr_conditional"
                    candidates.append(candidate)
        
        # =====================================================================
        # TRACK 2: HIGH + MEDIUM (respect humanization)
        # =====================================================================
        track2_rules = cdr_rules_by_tier['high'] + cdr_rules_by_tier['medium']
        if track2_rules:
            muts = list(hallmark_muts)
            used_positions = {m.imgt_num for m in muts}
            
            for pos, univ_aa, rule_aa, conf, cond in track2_rules:
                if pos not in used_positions:
                    if rule_aa == 'C' and not cdr3_has_cysteine:
                        continue
                    muts.append(Mutation(
                        position=f"IMGT{pos}",
                        imgt_num=pos,
                        original=univ_aa,
                        mutant=rule_aa,
                        source=f"cdr_conditional_{cond}",
                        confidence=conf
                    ))
                    used_positions.add(pos)
            
            if len(muts) > len(hallmark_muts):
                key = tuple(sorted([str(m) for m in muts]))
                if key not in used:
                    used.add(key)
                    generation_order += 1
                    candidate = self._build_universal_candidate(
                        f"Univ_{target_hallmarks}_T2_medium_{generation_order}",
                        cdrs, family, muts, "cdr_conditional_medium", generation_order,
                        f"Track2: Universal + high+med CDR-conditional ({len(muts)-len(hallmark_muts)} extra)", family
                    )
                    candidate.design_track = "paired_vernier"
                    candidate.track = "Track 2 (Universal)"
                    candidate.ranking_exempt = False  # v8.6: Track 2 is ranked
                    candidate.track_info = f"Universal_{target_hallmarks}_T2_medium"
                    candidate.track_subtype = "cdr_conditional"
                    candidates.append(candidate)
        
        # =====================================================================
        # TRACK 3: HIGH + MEDIUM (override humanization)
        # =====================================================================
        track3_rules = cdr_rules_by_tier['high'] + cdr_rules_by_tier['medium'] + cdr_rules_by_tier['medium_humanized']
        if track3_rules:
            muts = list(hallmark_muts)
            used_positions = {m.imgt_num for m in muts}
            
            for pos, univ_aa, rule_aa, conf, cond in track3_rules:
                if pos not in used_positions:
                    if rule_aa == 'C' and not cdr3_has_cysteine:
                        continue
                    muts.append(Mutation(
                        position=f"IMGT{pos}",
                        imgt_num=pos,
                        original=univ_aa,
                        mutant=rule_aa,
                        source=f"cdr_conditional_{cond}",
                        confidence=conf
                    ))
                    used_positions.add(pos)
            
            if len(muts) > len(hallmark_muts):
                key = tuple(sorted([str(m) for m in muts]))
                if key not in used:
                    used.add(key)
                    generation_order += 1
                    candidate = self._build_universal_candidate(
                        f"Univ_{target_hallmarks}_T3_override_{generation_order}",
                        cdrs, family, muts, "cdr_conditional_override", generation_order,
                        f"Track3: Universal + override humanization ({len(muts)-len(hallmark_muts)} extra)", family
                    )
                    candidate.design_track = "triplet_vernier"
                    candidate.track = "Track 3 (Universal)"
                    candidate.ranking_exempt = False  # v8.6: Track 3 is ranked
                    candidate.track_info = f"Universal_{target_hallmarks}_T3_override"
                    candidate.track_subtype = "cdr_conditional"
                    candidates.append(candidate)
        
        # =====================================================================
        # v8.8: SCALED MOTIF-BASED TRACK 2/3 FOR UNIVERSAL
        # =====================================================================
        # Add motif sampling to Universal scaffold (same as original scaffold)
        # This generates many more Track 2/3 candidates to fill quotas
        # =====================================================================
        
        # Build vernier consensus mutations from profile (already loaded at top)
        vernier_consensus_muts = []
        for pos in [67, 68, 69, 71, 76, 78, 82, 87, 89, 91, 94]:
            profile_data = vernier_profile.get(pos, {})
            consensus_aa = profile_data.get('cons', '')
            consensus_freq = float(profile_data.get('freq', 0))
            scaffold_aa = self.scaffold_positions.get(pos, '')
            
            if consensus_aa and scaffold_aa and consensus_aa != scaffold_aa:
                vernier_consensus_muts.append(Mutation(
                    position=f"IMGT{pos}",
                    imgt_num=pos,
                    original=scaffold_aa,
                    mutant=consensus_aa,
                    source=f"vernier_consensus_{target_hallmarks}",
                    confidence=consensus_freq
                ))
        
        # v8.8: Scale Track 2/3 generation
        track2_target = max(30, n // 20)  # 5% of n, min 30
        track3_target = max(30, n // 20)
        SPRINKLE_VARIANTS = 3
        
        print(f"    Track 2/3 (Universal) targets: Track2={track2_target}, Track3={track3_target}")
        
        # Track 2: Motifs 2-3 positions
        motifs_2_3 = getattr(self, "motif_pool_2_3", [])
        track2_count = 0
        track2_attempts = 0
        max_track2_attempts = track2_target * 10
        
        while track2_count < track2_target and track2_attempts < max_track2_attempts:
            track2_attempts += 1
            motif = sample_motif(motifs_2_3, min_len=2, max_len=3, rng=rng)
            if not motif:
                break
            
            for _ in range(SPRINKLE_VARIANTS):
                if track2_count >= track2_target:
                    break
                
                muts = list(hallmark_muts)
                used_positions = {m.imgt_num for m in muts}
                
                # Apply motif (v8.8.2: skip hallmark positions)
                for pos, aa in motif:
                    if pos in hallmark_set:
                        continue  # Never override hallmark positions
                    scaffold_aa = self.scaffold_positions.get(pos, "")
                    if scaffold_aa and scaffold_aa != aa and pos not in used_positions:
                        muts.append(Mutation(
                            position=f"IMGT{pos}", imgt_num=pos,
                            original=scaffold_aa, mutant=aa,
                            source="archetype_motif", confidence=1.0
                        ))
                        used_positions.add(pos)
                
                # Add 0-K extra verniers (capped at K_EXTRA_VERNIERS_T2)
                # Track 2 should remain motif-attributable: avoid unbounded Bernoulli across all verniers.
                temp = getattr(self, "sprinkle_temperature", 1.0)
                extra_pool = [vm for vm in vernier_consensus_muts if vm.imgt_num not in used_positions]
                k_extra = rng.randint(0, K_EXTRA_VERNIERS_T2)  # 0..K
                if extra_pool and k_extra > 0:
                    weights = [temp_scale_prob(float(vm.confidence), temp) for vm in extra_pool]
                    picked = weighted_sample_without_replacement(extra_pool, weights, k_extra, rng)
                    for vm in picked:
                        muts.append(vm)
                        used_positions.add(vm.imgt_num)
                
                key = tuple(sorted([str(m) for m in muts]))
                if key not in used and len(muts) > len(hallmark_muts):
                    used.add(key)
                    generation_order += 1
                    track2_count += 1
                    candidate = self._build_universal_candidate(
                        f"Univ_{target_hallmarks}_motif2_{generation_order}",
                        cdrs, family, muts, "paired_vernier", generation_order,
                        f"Track2: Universal motif (2-3)", family
                    )
                    candidate.design_track = "paired_vernier"
                    candidate.track = "Track 2 (Universal)"
                    candidate.ranking_exempt = False
                    candidate.track_info = "motif2-3:" + ",".join([f"IMGT{p}{aa}" for p, aa in motif])
                    candidates.append(candidate)
        
        print(f"      Track 2 (Universal): generated {track2_count} candidates")
        
        # Track 3: Motifs 3-5 positions
        motifs_3_5 = getattr(self, "motif_pool_3_5", [])
        track3_count = 0
        track3_attempts = 0
        max_track3_attempts = track3_target * 10
        
        while track3_count < track3_target and track3_attempts < max_track3_attempts:
            track3_attempts += 1
            motif = sample_motif(motifs_3_5, min_len=3, max_len=5, rng=rng)
            if not motif:
                break
            
            for _ in range(SPRINKLE_VARIANTS):
                if track3_count >= track3_target:
                    break
                
                muts = list(hallmark_muts)
                used_positions = {m.imgt_num for m in muts}
                
                # Apply motif (v8.8.2: skip hallmark positions)
                for pos, aa in motif:
                    if pos in hallmark_set:
                        continue  # Never override hallmark positions
                    scaffold_aa = self.scaffold_positions.get(pos, "")
                    if scaffold_aa and scaffold_aa != aa and pos not in used_positions:
                        muts.append(Mutation(
                            position=f"IMGT{pos}", imgt_num=pos,
                            original=scaffold_aa, mutant=aa,
                            source="archetype_motif", confidence=1.0
                        ))
                        used_positions.add(pos)
                
                # Sprinkling
                temp = getattr(self, "sprinkle_temperature", 1.0)
                for vm in vernier_consensus_muts:
                    if vm.imgt_num in used_positions:
                        continue
                    prob = temp_scale_prob(float(vm.confidence) * 0.5, temp)
                    if rng.random() < prob:
                        muts.append(vm)
                        used_positions.add(vm.imgt_num)
                
                key = tuple(sorted([str(m) for m in muts]))
                if key not in used and len(muts) > len(hallmark_muts):
                    used.add(key)
                    generation_order += 1
                    track3_count += 1
                    candidate = self._build_universal_candidate(
                        f"Univ_{target_hallmarks}_motif3_{generation_order}",
                        cdrs, family, muts, "motif_triplet", generation_order,
                        f"Track3: Universal motif (3-5)", family
                    )
                    candidate.design_track = "motif_triplet"
                    candidate.track = "Track 3 (Universal)"
                    candidate.ranking_exempt = False
                    candidate.track_info = "motif3-5:" + ",".join([f"IMGT{p}{aa}" for p, aa in motif])
                    candidates.append(candidate)
        
        print(f"      Track 3 (Universal): generated {track3_count} candidates")
        
        # =====================================================================
        # v8.5: TRACK 4 - EXPANSIVE OPTIMIZATION FOR UNIVERSAL
        # =====================================================================
        # Track 4 now includes:
        # - All CDR-conditional rules
        # - R50 consensus mutations
        # - Framework consensus tiers (S, A, B, C)
        # - Scaling complexity when mutation space is exhausted
        # =====================================================================
        
        print(f"    Track 4: EXPANSIVE optimization for Universal...")
        
        # Build framework consensus mutation pools
        tier_s_muts = []
        tier_a_muts = []
        tier_b_muts = []
        tier_c_muts = []
        
        for pos, (cons_aa, freq) in FRAMEWORK_CONSENSUS_TIERS['S'].items():
            scaffold_aa = self.scaffold_positions.get(pos, '')
            if scaffold_aa and scaffold_aa != cons_aa:
                if cons_aa == 'C' and not cdr3_has_cysteine:
                    continue
                tier_s_muts.append(Mutation(
                    position=f"IMGT{pos}", imgt_num=pos,
                    original=scaffold_aa, mutant=cons_aa,
                    source="fw_consensus_tier_S", confidence=freq
                ))
        
        for pos, (cons_aa, freq) in FRAMEWORK_CONSENSUS_TIERS['A'].items():
            scaffold_aa = self.scaffold_positions.get(pos, '')
            if scaffold_aa and scaffold_aa != cons_aa:
                if cons_aa == 'C' and not cdr3_has_cysteine:
                    continue
                tier_a_muts.append(Mutation(
                    position=f"IMGT{pos}", imgt_num=pos,
                    original=scaffold_aa, mutant=cons_aa,
                    source="fw_consensus_tier_A", confidence=freq
                ))
        
        for pos, (cons_aa, freq) in FRAMEWORK_CONSENSUS_TIERS['B'].items():
            scaffold_aa = self.scaffold_positions.get(pos, '')
            if scaffold_aa and scaffold_aa != cons_aa:
                if cons_aa == 'C' and not cdr3_has_cysteine:
                    continue
                tier_b_muts.append(Mutation(
                    position=f"IMGT{pos}", imgt_num=pos,
                    original=scaffold_aa, mutant=cons_aa,
                    source="fw_consensus_tier_B", confidence=freq
                ))
        
        for pos, (cons_aa, freq) in FRAMEWORK_CONSENSUS_TIERS['C'].items():
            scaffold_aa = self.scaffold_positions.get(pos, '')
            if scaffold_aa and scaffold_aa != cons_aa:
                if cons_aa == 'C' and not cdr3_has_cysteine:
                    continue
                tier_c_muts.append(Mutation(
                    position=f"IMGT{pos}", imgt_num=pos,
                    original=scaffold_aa, mutant=cons_aa,
                    source="fw_consensus_tier_C", confidence=freq
                ))
        
        # Build R50 consensus mutations
        r50_muts = []
        for pos, (univ_aa, r50_aa, conf) in R50_CONSENSUS_MUTATIONS.items():
            scaffold_aa = self.scaffold_positions.get(pos, '')
            if scaffold_aa and scaffold_aa != r50_aa:
                if r50_aa == 'C' and not cdr3_has_cysteine:
                    continue
                r50_muts.append(Mutation(
                    position=f"IMGT{pos}", imgt_num=pos,
                    original=scaffold_aa, mutant=r50_aa,
                    source="r50_consensus", confidence=conf
                ))
        
        # Convert CDR rules to mutation objects
        cdr_muts = []
        for pos, univ_aa, rule_aa, conf, cond in track3_rules:
            if rule_aa == 'C' and not cdr3_has_cysteine:
                continue
            cdr_muts.append(Mutation(
                position=f"IMGT{pos}", imgt_num=pos,
                original=univ_aa, mutant=rule_aa,
                source=f"cdr_conditional_{cond}", confidence=conf
            ))
        
        # Base mutations for Track 4
        base_muts = list(hallmark_muts)
        base_positions = {m.imgt_num for m in base_muts}
        
        # Expansive generation with scaling complexity
        n_remaining = n - len(candidates)
        max_attempts = max(2000, n_remaining * 50)
        attempts = 0
        base_prob = 0.4
        prob_increment = 0.02
        stagnation_threshold = 100
        last_count = len(candidates)
        stagnation_counter = 0
        
        while len(candidates) < n and attempts < max_attempts:
            attempts += 1
            
            # Check for stagnation and increase complexity
            if len(candidates) == last_count:
                stagnation_counter += 1
                if stagnation_counter >= stagnation_threshold:
                    base_prob = min(0.95, base_prob + prob_increment)
                    stagnation_counter = 0
            else:
                last_count = len(candidates)
                stagnation_counter = 0
            
            muts = list(base_muts)
            used_positions = set(base_positions)
            
            # v8.6: Temperature for probability scaling
            temp = getattr(self, 'sprinkle_temperature', 1.0)
            
            # Layer 1: CDR-conditional rules (highest priority)
            for cm in cdr_muts:
                if cm.imgt_num not in used_positions:
                    prob = temp_scale_prob(cm.confidence * base_prob * 1.2, temp)
                    if rng.random() < prob:
                        muts.append(cm)
                        used_positions.add(cm.imgt_num)
            
            # Layer 2: R50 consensus
            for rm in r50_muts:
                if rm.imgt_num not in used_positions:
                    prob = temp_scale_prob(rm.confidence * base_prob, temp)
                    if rng.random() < prob:
                        muts.append(rm)
                        used_positions.add(rm.imgt_num)
            
            # Layer 3: Tier S framework (high probability)
            for sm in tier_s_muts:
                if sm.imgt_num not in used_positions:
                    prob = temp_scale_prob(sm.confidence * base_prob * 1.3, temp)
                    if rng.random() < prob:
                        muts.append(sm)
                        used_positions.add(sm.imgt_num)
            
            # Layer 4: Tier A framework
            for am in tier_a_muts:
                if am.imgt_num not in used_positions:
                    prob = temp_scale_prob(am.confidence * base_prob * 1.1, temp)
                    if rng.random() < prob:
                        muts.append(am)
                        used_positions.add(am.imgt_num)
            
            # Layer 5: Tier B framework
            for bm in tier_b_muts:
                if bm.imgt_num not in used_positions:
                    prob = temp_scale_prob(bm.confidence * base_prob * 0.9, temp)
                    if rng.random() < prob:
                        muts.append(bm)
                        used_positions.add(bm.imgt_num)
            
            # Layer 6: Tier C framework (lower probability)
            for cm in tier_c_muts:
                if cm.imgt_num not in used_positions:
                    prob = temp_scale_prob(cm.confidence * base_prob * 0.7, temp)
                    if rng.random() < prob:
                        muts.append(cm)
                        used_positions.add(cm.imgt_num)
            
            key = tuple(sorted([str(m) for m in muts]))
            if key not in used:
                used.add(key)
                generation_order += 1
                n_cdr = len([m for m in muts if 'cdr_conditional' in m.source])
                n_fw = len([m for m in muts if 'fw_consensus' in m.source])
                n_r50 = len([m for m in muts if 'r50_consensus' in m.source])
                
                candidate = self._build_universal_candidate(
                    f"Univ_{target_hallmarks}_T4_exp_{generation_order}",
                    cdrs, family, muts, "expansive_universal", generation_order,
                    f"Track4: Universal expansive ({n_cdr}cdr + {n_fw}fw + {n_r50}r50)", family
                )
                candidate.design_track = "optimized"
                candidate.track = "Track 4 (Universal)"
                candidate.ranking_exempt = False  # RANKED globally
                candidate.track_info = f"Universal_{target_hallmarks}_T4_exp"
                candidate.track_subtype = "expansive"
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
        
        # v8.5: Get CDR3 features for CDR-conditional rules
        cdr3_features = get_cdr3_features(cdrs.cdr3)
        cdr3_has_cysteine = cdrs.cdr3.count('C') > 0
        
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
            print(f"    Position 66 (variable): {pos66_data[1]}→{pos66_data[2]} ({pos66_data[3]:.0%}) - will use weighted sampling")
        
        # =================================================================
        # TRACK 0: CONTROLS (always generated, ranking exempt)
        # =================================================================
        # v8.9 CONTROL STRUCTURE (ALL HALLMARKS):
        #   0a: GRAFT (0 mutations) - pure CDR graft, ALL hallmarks get this
        #   0b: Hallmarks + IMGT2 (5 mutations) - the minimal VHH conversion
        #   0c: + ALL verniers (~20 mutations) - only if avg vernier ≥60%
        #   0d+: YQRL only: Progressive FW pairs (>80% consensus)
        #
        # IMPORTANT: Track 0 candidates are IMMUTABLE - no post-processing.
        # Base sequence: Original input framework (M69) with CDRs grafted
        # =================================================================
        
        is_yqrl = (target_hallmarks == 'YQRL')
        
        # =================================================================
        # Track 0a: GRAFT (ALL HALLMARKS) - 0 mutations
        # =================================================================
        # Base: Original framework + CDRs (no mutations)
        # Tests: Do the CDRs work in completely original context?
        # This is the most conservative control - zero framework changes.
        # =================================================================
        generation_order += 1
        key = (f"GRAFT_{target_hallmarks}",)  # Unique key for graft
        if key not in used:
            used.add(key)
            candidate = self._build_original_candidate(
                f"Orig_{family}_{target_hallmarks}_graft_{generation_order}",
                cdrs, family, [], "minimal_hallmark", generation_order,
                f"Track0: {target_hallmarks} GRAFT (0 mut) - pure CDR graft - IMMUTABLE", family, original_positions
            )
            candidate.design_track = "minimal_hallmark"
            candidate.track = "Track 0a (GRAFT)"
            candidate.ranking_exempt = True
            candidate.track_info = f"{target_hallmarks}_graft_0mut"
            candidate.applied_comp_rules = []
            candidate.is_immutable = True
            candidate.track_subtype = "graft"
            candidates.append(candidate)
            print(f"    Generated: {target_hallmarks} GRAFT control (0 mutations)")
        
        # =================================================================
        # Track 0b: Hallmarks + IMGT2 (5 mutations) - ALL HALLMARKS
        # =================================================================
        # Base: Original framework + CDRs
        # Mutations: IMGT42, IMGT49, IMGT50, IMGT52 (hallmarks) + IMGT2→V
        # Tests: Does minimal VHH conversion enable function?
        # v8.9: Always include IMGT2 with hallmarks (no separate 4-mut track)
        # =================================================================
        generation_order += 1
        key = tuple(sorted([str(m) for m in hallmark_muts]))
        if key not in used:
            used.add(key)
            candidate = self._build_original_candidate(
                f"Orig_{family}_{target_hallmarks}_5mut_{generation_order}",
                cdrs, family, hallmark_muts, "minimal_hallmark", generation_order,
                f"Track0: Hallmarks + IMGT2 (5 mut) - IMMUTABLE", family, original_positions
            )
            candidate.design_track = "minimal_hallmark"
            candidate.track = "Track 0b"
            candidate.ranking_exempt = True
            candidate.track_info = f"{target_hallmarks}_5mut"
            candidate.applied_comp_rules = []
            candidate.is_immutable = True
            candidate.track_subtype = "hallmarks_imgt2"
            candidates.append(candidate)
            print(f"    Generated: {target_hallmarks} hallmarks + IMGT2 control (5 mutations)")
        
        # =================================================================
        # Track 0c: + ALL VERNIERS (~20 mutations)
        # =================================================================
        # Base: Original framework + CDRs
        # Mutations: Hallmarks + IMGT2 + all 16 non-hallmark verniers
        # Only generated if hallmark has average vernier consensus ≥60%
        # v8.9: Threshold reduced from 75% to 60%
        # =================================================================
        
        # Check if this hallmark has sufficient average vernier consensus
        vernier_freqs = [vm.confidence for vm in vernier_consensus_muts]
        is_sufficient_consensus = False
        if vernier_freqs:
            avg_vernier_freq = sum(vernier_freqs) / len(vernier_freqs)
            is_sufficient_consensus = avg_vernier_freq >= 0.60  # v8.9: 60% threshold
            avg_vernier_freq = sum(vernier_freqs) / len(vernier_freqs)
            is_high_consensus = avg_vernier_freq >= 0.75
        
        # Generate full vernier control if sufficient consensus
        if is_yqrl or is_sufficient_consensus:
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
                    f"Orig_{family}_{target_hallmarks}_full_vernier_{generation_order}",
                    cdrs, family, all_vernier_muts, "minimal_hallmark", generation_order,
                    f"Track0: Hallmarks + IMGT2 + ALL verniers ({target_hallmarks}) - IMMUTABLE", family, original_positions
                )
                candidate.design_track = "minimal_hallmark"
                candidate.track = "Track 0c"
                candidate.ranking_exempt = True
                candidate.track_info = f"{target_hallmarks}_full_vernier"
                candidate.applied_comp_rules = []
                candidate.is_immutable = True
                candidate.track_subtype = "full_vernier"
                candidates.append(candidate)
                print(f"    Generated: {target_hallmarks} full vernier control ({len(all_vernier_muts)} mutations)")
        
        # =================================================================
        # v8.9: YQRL PROGRESSIVE FRAMEWORK IN PAIRS (Track 0d+)
        # =================================================================
        # After full vernier control, add framework positions in PAIRS
        # based on covariance/structural coupling. All positions >80% consensus.
        # =================================================================
        if is_yqrl:
            print(f"    Generating YQRL progressive framework controls (FW pairs, {len(YQRL_FRAMEWORK_PAIRS)} steps)...")
            
            # Start with all vernier mutations from 0c
            base_vernier_muts = list(all_vernier_muts) if 'all_vernier_muts' in dir() and all_vernier_muts else list(hallmark_muts)
            if not all_vernier_muts:
                # Rebuild vernier mutations if needed
                base_vernier_muts = list(hallmark_muts)
                for vm in vernier_consensus_muts:
                    if vm.confidence >= 0.50:
                        base_vernier_muts.append(vm)
                if pos66_data:
                    pos, orig, cons, freq = pos66_data[:4]
                    if cons != orig:
                        base_vernier_muts.append(Mutation(
                            position=f"IMGT{pos}",
                            imgt_num=pos,
                            original=orig,
                            mutant=cons,
                            source=f"vernier_consensus_{target_hallmarks}",
                            confidence=freq
                        ))
            
            # Track which positions are already mutated
            mutated_positions = {m.imgt_num for m in base_vernier_muts}
            
            # Cumulative FW mutations - each step adds a PAIR
            cumulative_fw_muts = []
            
            # Generate progressive FW steps - each adds a PAIR (or singleton for last)
            for fw_step, pair in enumerate(YQRL_FRAMEWORK_PAIRS, 1):
                # Add mutations from this pair
                pair_muts = []
                pair_positions = []
                for pos, target_aa, freq in pair:
                    if pos in mutated_positions:
                        continue  # Already mutated by vernier or hallmark
                    orig_aa = original_positions.get(pos, '')
                    if orig_aa and orig_aa != target_aa:
                        pair_muts.append(Mutation(
                            position=f"IMGT{pos}",
                            imgt_num=pos,
                            original=orig_aa,
                            mutant=target_aa,
                            source=f"yqrl_fw_pair{fw_step}",
                            confidence=freq
                        ))
                        pair_positions.append(pos)
                
                if not pair_muts:
                    continue  # All positions in this pair already correct
                
                # Add this pair's mutations to cumulative
                cumulative_fw_muts.extend(pair_muts)
                for m in pair_muts:
                    mutated_positions.add(m.imgt_num)
                
                # Combine vernier + cumulative framework mutations
                all_muts = list(base_vernier_muts) + cumulative_fw_muts
                
                key = tuple(sorted([str(m) for m in all_muts]))
                if key not in used:
                    used.add(key)
                    generation_order += 1
                    
                    # Track name based on step
                    if fw_step <= 26:
                        track_letter = chr(ord('d') + fw_step - 1)  # 0d, 0e, 0f, ...
                    else:
                        track_letter = f"d{fw_step}"  # For many steps
                    
                    # Description of this pair
                    pair_desc = "+".join([f"IMGT{p}" for p in pair_positions])
                    
                    candidate = self._build_original_candidate(
                        f"Orig_{family}_YQRL_FW{fw_step}_{generation_order}",
                        cdrs, family, all_muts, "minimal_hallmark", generation_order,
                        f"Track0: YQRL FW pair {fw_step} ({pair_desc}) - IMMUTABLE", family, original_positions
                    )
                    candidate.design_track = "minimal_hallmark"
                    candidate.track = f"Track 0{track_letter}"
                    candidate.ranking_exempt = True
                    candidate.track_info = f"YQRL_FW{fw_step}_{pair_desc}"
                    candidate.applied_comp_rules = []
                    candidate.is_immutable = True
                    candidate.track_subtype = f"fw_pair_{fw_step}"
                    candidates.append(candidate)
            
            print(f"    Generated YQRL progressive framework controls ({len(cumulative_fw_muts)} total FW mutations)")
        
        # =================================================================
        # TRACK 1: SINGLE VERNIER PROBES (ranked)
        # =================================================================
        # v8.9: Include ALL verniers, not just "load-bearing"
        # Test individual verniers to identify critical positions
        # Base: Hallmarks + IMGT2 + ONE vernier
        # =================================================================
        
        # v8.9: Use ALL vernier consensus mutations (sorted by confidence)
        all_probe_verniers = sorted(vernier_consensus_muts, key=lambda m: -m.confidence)
        
        # Also include pos66 if available
        if pos66_data:
            pos, orig, cons, freq, alt, alt_freq = pos66_data
            if cons != orig:
                pos66_mut = Mutation(
                    position=f"IMGT{pos}",
                    imgt_num=pos,
                    original=orig,
                    mutant=cons,
                    source=f"vernier_consensus_{target_hallmarks}",
                    confidence=freq
                )
                # Add if not already present
                if not any(v.imgt_num == 66 for v in all_probe_verniers):
                    all_probe_verniers.append(pos66_mut)
        
        print(f"    Track 1: {len(all_probe_verniers)} vernier positions available for probing")
        
        # Generate probes - v8.9: no sprinkling variants, just clean probes
        probe_count = 0
        for vm in all_probe_verniers:
            if probe_count >= MAX_SINGLE_PROBES:
                break
            
            # Base = hallmarks + IMGT2 + this one vernier
            muts = list(hallmark_muts) + [vm]
            
            key = tuple(sorted([str(m) for m in muts]))
            if key not in used:
                used.add(key)
                generation_order += 1
                probe_count += 1
                candidate = self._build_original_candidate(
                    f"Orig_{family}_probe_IMGT{vm.imgt_num}_{generation_order}",
                    cdrs, family, muts, "single_vernier", generation_order,
                    f"Track1: Hallmarks + IMGT2 + IMGT{vm.imgt_num} probe", family, original_positions
                )
                candidate.design_track = "single_vernier"
                candidate.track = "Track 1"
                candidate.ranking_exempt = False  # Ranked
                candidate.track_info = f"IMGT{vm.imgt_num}_{vm.original}>{vm.mutant}"
                candidate.applied_comp_rules = []  # No compensation for probes
                candidates.append(candidate)
        
        print(f"      Generated {probe_count} Track 1 probes")
        
        # =================================================================
        # TRACK 2: MINIMAL STRUCTURED PERTURBATIONS (v8.9)
        # =================================================================
        # Goal: Test archetype motifs with minimal additional changes
        # Base: Hallmarks + IMGT2
        # Add: Exactly ONE 2-3 position motif (weighted by family support)
        # Add: 0-K extra vernier mutations (K_EXTRA_VERNIERS_T2 = 2, capped)
        # NO framework tier sprinkling - keeps it interpretable
        # =================================================================
        
        # v8.9: Protected positions - motifs cannot override
        protected_positions = PROTECTED_POSITIONS  # {2, 42, 49, 50, 52}
        
        # Scale track targets
        track2_target = max(30, n // 20)  # 5% of n, min 30
        track3_target = max(30, n // 20)  # 5% of n, min 30
        
        print(f"    Track 2/3 targets: Track2={track2_target}, Track3={track3_target}")
        
        motifs_2_3 = getattr(self, "motif_pool_2_3", [])
        if not motifs_2_3:
            motifs_2_3 = []
        
        # Build pool of extra verniers (excluding protected positions and motif positions)
        # Sorted by confidence for weighted sampling
        extra_vernier_pool = [vm for vm in vernier_consensus_muts 
                              if vm.imgt_num not in protected_positions and vm.confidence >= 0.50]
        extra_vernier_pool.sort(key=lambda m: -m.confidence)
        
        track2_count = 0
        track2_attempts = 0
        max_track2_attempts = track2_target * 10
        
        while track2_count < track2_target and track2_attempts < max_track2_attempts:
            track2_attempts += 1
            
            # Sample ONE motif (weighted by support)
            motif = sample_motif(motifs_2_3, min_len=2, max_len=3, rng=rng)
            if not motif:
                break
            
            # Filter out protected positions from motif
            motif = [(pos, aa) for pos, aa in motif if pos not in protected_positions]
            if len(motif) < 2:
                continue  # Motif too small after filtering
            
            # Base = hallmarks + IMGT2
            muts = list(hallmark_muts)
            used_positions = {m.imgt_num for m in muts}
            
            # Apply motif deterministically
            motif_applied = []
            for pos, aa in motif:
                orig = original_positions.get(pos, "")
                if orig and orig != aa and pos not in used_positions:
                    muts.append(Mutation(
                        position=f"IMGT{pos}", imgt_num=pos,
                        original=orig, mutant=aa,
                        source="archetype_motif_t2",
                        confidence=1.0
                    ))
                    used_positions.add(pos)
                    motif_applied.append((pos, aa))
            
            if not motif_applied:
                continue  # No mutations from this motif
            
            # Add 0-K extra verniers (capped at K_EXTRA_VERNIERS_T2)
            # Sample weighted by confidence
            available_verniers = [vm for vm in extra_vernier_pool if vm.imgt_num not in used_positions]
            k_extra = rng.randint(0, K_EXTRA_VERNIERS_T2)  # 0..K
            
            if available_verniers and k_extra > 0:
                # Weighted sampling by confidence (WITHOUT replacement)
                weights = [vm.confidence for vm in available_verniers]
                picked = weighted_sample_without_replacement(
                    available_verniers, weights, min(k_extra, len(available_verniers)), rng
                )
                for vm in picked:
                    muts.append(vm)
                    used_positions.add(vm.imgt_num)

            
            key = tuple(sorted([str(m) for m in muts]))
            if key not in used and len(muts) > len(hallmark_muts):
                used.add(key)
                generation_order += 1
                track2_count += 1
                
                motif_desc = ",".join([f"IMGT{p}{aa}" for p, aa in motif_applied])
                n_extra = len(muts) - len(hallmark_muts) - len(motif_applied)
                
                candidate = self._build_original_candidate(
                    f"Orig_{family}_t2motif_{generation_order}",
                    cdrs, family, muts, "minimal_structured", generation_order,
                    f"Track2: Motif({motif_desc}) + {n_extra} extra verniers", family, original_positions
                )
                candidate.design_track = "minimal_structured"
                candidate.track = "Track 2"
                candidate.ranking_exempt = False
                candidate.track_info = f"motif:{motif_desc}+{n_extra}v"
                candidate.applied_comp_rules = []
                candidates.append(candidate)
        
        print(f"      Track 2: generated {track2_count} candidates (minimal structured)")
        
        # =================================================================
        # TRACK 3: MOTIF ATTRIBUTION (v8.9)
        # =================================================================
        # Goal: Test 3-5 position archetype motifs for attributable results
        # Base: Hallmarks + IMGT2
        # Add: Exactly ONE 3-5 position motif (weighted by support)
        # Add: 0-K extra verniers (K_EXTRA_VERNIERS_T3 = 2, max 3)
        # NO framework tier sprinkling - interpretability is key
        # =================================================================
        
        motifs_3_5 = getattr(self, "motif_pool_3_5", [])
        if not motifs_3_5:
            motifs_3_5 = []
        
        track3_count = 0
        track3_attempts = 0
        max_track3_attempts = track3_target * 10
        
        while track3_count < track3_target and track3_attempts < max_track3_attempts:
            track3_attempts += 1
            
            # Sample ONE motif (3-5 positions, weighted)
            motif = sample_motif(motifs_3_5, min_len=3, max_len=5, rng=rng)
            if not motif:
                break
            
            # Filter out protected positions from motif
            motif = [(pos, aa) for pos, aa in motif if pos not in protected_positions]
            if len(motif) < 3:
                continue  # Motif too small after filtering
            
            # Base = hallmarks + IMGT2
            muts = list(hallmark_muts)
            used_positions = {m.imgt_num for m in muts}
            
            # Apply motif deterministically (motif positions are LOCKED)
            motif_applied = []
            for pos, aa in motif:
                orig = original_positions.get(pos, "")
                if orig and orig != aa and pos not in used_positions:
                    muts.append(Mutation(
                        position=f"IMGT{pos}", imgt_num=pos,
                        original=orig, mutant=aa,
                        source="archetype_motif_t3",
                        confidence=1.0
                    ))
                    used_positions.add(pos)
                    motif_applied.append((pos, aa))
            
            if len(motif_applied) < 3:
                continue  # Need at least 3 positions from motif
            
            # Add 0-K extra verniers (capped at K_EXTRA_VERNIERS_T3_MAX)
            available_verniers = [vm for vm in extra_vernier_pool if vm.imgt_num not in used_positions]
            k_extra = rng.randint(0, K_EXTRA_VERNIERS_T3)  # Default 0-2
            k_extra = min(k_extra, K_EXTRA_VERNIERS_T3_MAX)  # Hard cap at 3
            
            if available_verniers and k_extra > 0:
                # Weighted sampling by confidence (WITHOUT replacement)
                weights = [vm.confidence for vm in available_verniers]
                picked = weighted_sample_without_replacement(
                    available_verniers, weights, min(k_extra, len(available_verniers)), rng
                )
                for vm in picked:
                    muts.append(vm)
                    used_positions.add(vm.imgt_num)

            
            key = tuple(sorted([str(m) for m in muts]))
            if key not in used:
                used.add(key)
                generation_order += 1
                track3_count += 1
                
                motif_desc = ",".join([f"IMGT{p}{aa}" for p, aa in motif_applied])
                n_extra = len(muts) - len(hallmark_muts) - len(motif_applied)
                
                candidate = self._build_original_candidate(
                    f"Orig_{family}_t3motif_{generation_order}",
                    cdrs, family, muts, "motif_attribution", generation_order,
                    f"Track3: Motif({motif_desc}) + {n_extra} verniers", family, original_positions
                )
                candidate.design_track = "motif_attribution"
                candidate.track = "Track 3"
                candidate.ranking_exempt = False
                candidate.track_info = f"motif:{motif_desc}+{n_extra}v"
                candidate.applied_comp_rules = []
                candidates.append(candidate)
        
        print(f"      Track 3: generated {track3_count} candidates (motif attribution)")

        # Count controls generated so far
        n_controls = len([c for c in candidates if c.ranking_exempt])

        # =================================================================
        # TRACK 4: PRODUCTION OPTIMIZER (v8.9)
        # =================================================================
        # Goal: Generate bulk candidates with high success probability
        # Coverage matters more than attribution
        #
        # v8.9: CAPPED LAYERS (not unbounded Bernoulli)
        # 1. Base: Hallmarks + IMGT2 (always)
        # 2. Vernier layer: Bernoulli + cap (K_VERNIER_MAX = 8)
        # 3. Framework tiers: Bernoulli + per-tier caps
        #    - Tier S: K_S_MAX = 10
        #    - Tier A: K_A_MAX = 6
        #    - Tier B: K_B_MAX = 3
        #    - Tier C: K_C_MAX = 2 (only when desperate)
        # 4. CDR-conditional rules: K_CDR_RULES_MAX = 2
        # 5. Adaptive expansion when uniqueness drops
        # =================================================================
        
        print(f"    Track 4: PRODUCTION OPTIMIZER...")
        
        # Calculate allocation
        n_remaining = n - n_controls
        print(f"    Track allocation: T0-T3={n_controls}, T4(production)={n_remaining}")
        
        # Build all available mutation pools (excluding protected positions)
        
        tier_s_muts = []
        for pos, (cons_aa, freq) in FRAMEWORK_CONSENSUS_TIERS.get('S', {}).items():
            if pos in protected_positions:
                continue
            orig_aa = original_positions.get(pos, '')
            if orig_aa and orig_aa != cons_aa:
                if cons_aa == 'C' and not cdr3_has_cysteine:
                    continue
                tier_s_muts.append(Mutation(
                    position=f"IMGT{pos}",
                    imgt_num=pos,
                    original=orig_aa,
                    mutant=cons_aa,
                    source="fw_consensus_tier_S",
                    confidence=freq
                ))
        tier_s_muts.sort(key=lambda m: -m.confidence)

        tier_a_muts = []
        for pos, (cons_aa, freq) in FRAMEWORK_CONSENSUS_TIERS.get('A', {}).items():
            if pos in protected_positions:
                continue
            orig_aa = original_positions.get(pos, '')
            if orig_aa and orig_aa != cons_aa:
                if cons_aa == 'C' and not cdr3_has_cysteine:
                    continue
                tier_a_muts.append(Mutation(
                    position=f"IMGT{pos}",
                    imgt_num=pos,
                    original=orig_aa,
                    mutant=cons_aa,
                    source="fw_consensus_tier_A",
                    confidence=freq
                ))
        tier_a_muts.sort(key=lambda m: -m.confidence)

        tier_b_muts = []
        for pos, (cons_aa, freq) in FRAMEWORK_CONSENSUS_TIERS['B'].items():
            if pos in protected_positions:
                continue
            orig_aa = original_positions.get(pos, '')
            if orig_aa and orig_aa != cons_aa:
                if cons_aa == 'C' and not cdr3_has_cysteine:
                    continue
                tier_b_muts.append(Mutation(
                    position=f"IMGT{pos}",
                    imgt_num=pos,
                    original=orig_aa,
                    mutant=cons_aa,
                    source="fw_consensus_tier_B",
                    confidence=freq
                ))
        tier_b_muts.sort(key=lambda m: -m.confidence)
        
        tier_c_muts = []
        for pos, (cons_aa, freq) in FRAMEWORK_CONSENSUS_TIERS['C'].items():
            if pos in protected_positions:
                continue
            orig_aa = original_positions.get(pos, '')
            if orig_aa and orig_aa != cons_aa:
                if cons_aa == 'C' and not cdr3_has_cysteine:
                    continue
                tier_c_muts.append(Mutation(
                    position=f"IMGT{pos}",
                    imgt_num=pos,
                    original=orig_aa,
                    mutant=cons_aa,
                    source="fw_consensus_tier_C",
                    confidence=freq
                ))
        tier_c_muts.sort(key=lambda m: -m.confidence)
        
        # Vernier pool for Track 4 (excluding protected positions)
        all_verniers = [vm for vm in vernier_consensus_muts if vm.imgt_num not in protected_positions]
        all_verniers.sort(key=lambda m: -m.confidence)
        
        # Position 66 handling - weighted choice between consensus and alternative
        pos66_weighted = None
        if pos66_data:
            pos, orig, cons, cons_freq, alt, alt_freq = pos66_data
            if cons != orig:
                pos66_weighted = {
                    'pos': pos,
                    'orig': orig,
                    'options': [cons, alt] if alt and alt_freq >= 0.05 else [cons],
                    'weights': [cons_freq, alt_freq] if alt and alt_freq >= 0.05 else [cons_freq],
                }
        
        # CDR-conditional rules (if available)
        cdr_rule_muts = list(rule_muts) if 'rule_muts' in dir() and rule_muts else []
        
        # Adaptive expansion state
        max_attempts = max(2000, n_remaining * 50)
        attempts = 0
        stagnation_threshold = 100
        last_count = len(candidates)
        stagnation_counter = 0
        
        # v8.9: Current caps (start at defaults, expand when struggling)
        current_k_vernier = K_VERNIER_MAX
        current_k_s = K_S_MAX
        current_k_a = K_A_MAX
        current_k_b = K_B_MAX
        current_k_c = 0  # Start with no Tier C
        expansion_level = 0
        
        while len(candidates) < n and attempts < max_attempts:
            attempts += 1
            
            # Check for stagnation and trigger adaptive expansion
            if len(candidates) == last_count:
                stagnation_counter += 1
                if stagnation_counter >= stagnation_threshold:
                    expansion_level += 1
                    stagnation_counter = 0
                    
                    # Expansion order (as per spec):
                    # 1. Increase K_A_max
                    # 2. Increase K_B_max
                    # 3. Increase K_S_max
                    # 4. Enable Tier C with low cap
                    if expansion_level == 1:
                        current_k_a = min(current_k_a + 2, len(tier_a_muts))
                        print(f"      Expansion L1: K_A_max → {current_k_a}")
                    elif expansion_level == 2:
                        current_k_b = min(current_k_b + 2, len(tier_b_muts))
                        print(f"      Expansion L2: K_B_max → {current_k_b}")
                    elif expansion_level == 3:
                        current_k_s = min(current_k_s + 3, len(tier_s_muts))
                        print(f"      Expansion L3: K_S_max → {current_k_s}")
                    elif expansion_level == 4:
                        current_k_c = K_C_MAX  # Enable Tier C
                        print(f"      Expansion L4: Enabling Tier C (K_C={current_k_c})")
                    elif expansion_level >= 5:
                        current_k_vernier = min(current_k_vernier + 2, len(all_verniers))
                        print(f"      Expansion L5+: K_vernier → {current_k_vernier}")
            else:
                last_count = len(candidates)
                stagnation_counter = 0
            
            # Build candidate with capped layers
            muts = list(hallmark_muts)  # Base: hallmarks + IMGT2
            used_positions = {m.imgt_num for m in muts}
            
            temp = getattr(self, 'sprinkle_temperature', 1.0)
            
            # ========== LAYER 1: VERNIERS (Bernoulli + cap) ==========
            vernier_selected = []
            for vm in all_verniers:
                if vm.imgt_num in used_positions:
                    continue
                # Probability based on confidence
                prob = temp_scale_prob(vm.confidence, temp) * 0.6  # Base 60% scaling
                if rng.random() < prob:
                    vernier_selected.append(vm)
            
            # Cap verniers
            if len(vernier_selected) > current_k_vernier:
                vernier_selected.sort(key=lambda m: -m.confidence)
                vernier_selected = vernier_selected[:current_k_vernier]
            
            for vm in vernier_selected:
                muts.append(vm)
                used_positions.add(vm.imgt_num)
            
            # Position 66 special handling (weighted choice)
            if pos66_weighted and pos66_weighted['pos'] not in used_positions:
                if rng.random() < 0.7:  # 70% chance to include pos66
                    chosen_aa = weighted_choice(
                        pos66_weighted['options'],
                        pos66_weighted['weights'],
                        rng
                    )
                    if chosen_aa != pos66_weighted['orig']:
                        muts.append(Mutation(
                            position=f"IMGT{pos66_weighted['pos']}",
                            imgt_num=pos66_weighted['pos'],
                            original=pos66_weighted['orig'],
                            mutant=chosen_aa,
                            source=f"vernier_pos66_weighted",
                            confidence=max(pos66_weighted['weights'])
                        ))
                        used_positions.add(pos66_weighted['pos'])
            
            # ========== LAYER 2: TIER S FRAMEWORK (Bernoulli + cap) ==========
            tier_s_selected = []
            for sm in tier_s_muts:
                if sm.imgt_num in used_positions:
                    continue
                prob = temp_scale_prob(sm.confidence, temp) * TIER_S_PROB
                if rng.random() < prob:
                    tier_s_selected.append(sm)
            
            if len(tier_s_selected) > current_k_s:
                tier_s_selected.sort(key=lambda m: -m.confidence)
                tier_s_selected = tier_s_selected[:current_k_s]
            
            for sm in tier_s_selected:
                muts.append(sm)
                used_positions.add(sm.imgt_num)
            
            # ========== LAYER 3: TIER A FRAMEWORK (Bernoulli + cap) ==========
            tier_a_selected = []
            for am in tier_a_muts:
                if am.imgt_num in used_positions:
                    continue
                prob = temp_scale_prob(am.confidence, temp) * TIER_A_PROB
                if rng.random() < prob:
                    tier_a_selected.append(am)
            
            if len(tier_a_selected) > current_k_a:
                tier_a_selected.sort(key=lambda m: -m.confidence)
                tier_a_selected = tier_a_selected[:current_k_a]
            
            for am in tier_a_selected:
                muts.append(am)
                used_positions.add(am.imgt_num)
            
            # ========== LAYER 4: TIER B FRAMEWORK (Bernoulli + cap) ==========
            tier_b_selected = []
            for bm in tier_b_muts:
                if bm.imgt_num in used_positions:
                    continue
                prob = temp_scale_prob(bm.confidence, temp) * TIER_B_PROB
                if rng.random() < prob:
                    tier_b_selected.append(bm)
            
            if len(tier_b_selected) > current_k_b:
                tier_b_selected.sort(key=lambda m: -m.confidence)
                tier_b_selected = tier_b_selected[:current_k_b]
            
            for bm in tier_b_selected:
                muts.append(bm)
                used_positions.add(bm.imgt_num)
            
            # ========== LAYER 5: TIER C FRAMEWORK (only when enabled) ==========
            if current_k_c > 0:
                tier_c_selected = []
                for cm in tier_c_muts:
                    if cm.imgt_num in used_positions:
                        continue
                    prob = temp_scale_prob(cm.confidence, temp) * TIER_C_PROB
                    if rng.random() < prob:
                        tier_c_selected.append(cm)
                
                if len(tier_c_selected) > current_k_c:
                    tier_c_selected.sort(key=lambda m: -m.confidence)
                    tier_c_selected = tier_c_selected[:current_k_c]
                
                for cm in tier_c_selected:
                    muts.append(cm)
                    used_positions.add(cm.imgt_num)
            
            # ========== LAYER 6: CDR-CONDITIONAL RULES (capped) ==========
            cdr_selected = []
            for crm in cdr_rule_muts:
                if crm.imgt_num in used_positions:
                    continue
                prob = temp_scale_prob(crm.confidence, temp) * 0.3
                if rng.random() < prob:
                    cdr_selected.append(crm)
            
            if len(cdr_selected) > K_CDR_RULES_MAX:
                cdr_selected.sort(key=lambda m: -m.confidence)
                cdr_selected = cdr_selected[:K_CDR_RULES_MAX]
            
            for crm in cdr_selected:
                muts.append(crm)
                used_positions.add(crm.imgt_num)
            
            # ========== DEDUP AND ADD CANDIDATE ==========
            key = tuple(sorted([str(m) for m in muts]))
            if key not in used:
                used.add(key)
                generation_order += 1
                
                n_vern = len([m for m in muts if 'vernier' in m.source.lower()])
                n_fw = len([m for m in muts if 'fw_consensus' in m.source])
                n_cdr = len([m for m in muts if 'cdr' in m.source.lower() or 'rule' in m.source.lower()])
                
                candidate = self._build_original_candidate(
                    f"Orig_{family}_T4_prod_{generation_order}",
                    cdrs, family, muts, "production_optimized", generation_order,
                    f"Track4: Production ({n_vern}v + {n_fw}fw + {n_cdr}cdr)", 
                    family, original_positions
                )
                candidate.design_track = "production_optimizer"
                candidate.track = "Track 4"
                candidate.ranking_exempt = False  # RANKED
                candidate.track_info = f"T4_prod_{n_vern}v_{n_fw}fw"
                candidate.track_subtype = "production"
                candidate.applied_comp_rules = [m.source for m in cdr_selected]
                candidates.append(candidate)
        
        print(f"    Generated {len(candidates)} total candidates (attempts={attempts}, expansion_level={expansion_level})")
        
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
        
        # v8.8.2: Explicitly set scaffold_type (don't rely on dataclass default)
        cand.scaffold_type = "original"
        
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
            
            # Design track fields (v7.12/7.14/v8.1/v8.5)
            'design_track': c.design_track,
            'track': getattr(c, 'track', ''),  # v8.5: Explicit track label
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
        'rank_display': getattr(candidate, 'rank_display', str(candidate.rank)),
        'is_lead': candidate.is_lead,
        'is_control': getattr(candidate, 'ranking_exempt', False),
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
            # Use rank_display for clearer output (controls have brackets)
            rank_str = row.get('rank_display', str(row['rank']))
            ctrl_flag = "_CTRL" if row.get('is_control', False) else ""
            header = f">{rank_str}{ctrl_flag}|{row['id']}|{row['framework_source']}|{score_str}|plddt={row['plddt_mean']:.1f}"
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
    p.add_argument('--seed', type=int, default=42,
               help='Random seed for reproducibility (python + numpy)')

    
    p.add_argument('--n-generate', type=int, default=500)
    p.add_argument('--n-select', type=int, default=198, 
                   help='v8.6: TOTAL final output count INCLUDING controls + ranked')
    p.add_argument('--target-hallmarks', '-t', default='FERG')
    p.add_argument('--target-families', nargs='+', default=['F_C2', 'Y_C2', 'Other_VHH'])
    # v8.0: mode removed - always generates BOTH original and Universal
    
    # v8.6: Track quota controls
    p.add_argument('--bucket-target', type=int, default=25,
                   help='v8.6: Max candidates per bucket (hallmark × scaffold × family)')
    p.add_argument('--track0-quota', type=int, default=2, help='v8.6: Track 0 quota per bucket')
    p.add_argument('--track1-quota', type=int, default=5, help='v8.6: Track 1 quota per bucket')
    p.add_argument('--track2-quota', type=int, default=6, help='v8.6: Track 2 quota per bucket')
    p.add_argument('--track3-quota', type=int, default=6, help='v8.6: Track 3 quota per bucket')
    p.add_argument('--track4-quota', type=int, default=6, help='v8.6: Track 4 quota per bucket')
    p.add_argument('--dedup-mode', choices=['within_track', 'global', 'none'], default='within_track',
                   help='v8.6: Deduplication strategy')
    p.add_argument('--sprinkle-temp', type=float, default=1.0,
                   help='v8.6: Temperature for Bernoulli sprinkling (1.0=unchanged, <1=sharper)')
    
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

    random.seed(args.seed)
    np.random.seed(args.seed)
    
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
    print(f"VHH DESIGNER v{VERSION} - Archetype Motifs & Per-Position Sprinkling")
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
    
    generator = CandidateGenerator(
        args.rules, args.archetypes, 
        hallmark_db_path=args.hallmark_db,
        sprinkle_temperature=args.sprinkle_temp  # v8.6
    )
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
    lead_candidate = next((c for c in candidates if getattr(c, "is_lead", False)), None)
    
    # IMPORTANT: remove lead from the pool completely - main will add it back at the end
    non_lead = [c for c in candidates if not getattr(c, "is_lead", False)]
    
    # Track removal reasons for debugging
    removal_reasons = defaultdict(int)
    filtered = []
    for c in non_lead:
        a42 = c.imgt_positions.get(42, '')
        a49 = c.imgt_positions.get(49, '')
        a50 = c.imgt_positions.get(50, '')
        a52 = c.imgt_positions.get(52, '')
        
        if a50 != 'R':
            removal_reasons['pos50_not_R'] += 1
            continue
        if a52 == 'W':
            removal_reasons['pos52_W_VHlike'] += 1
            continue
        if a42 not in ('F', 'Y'):
            removal_reasons[f'pos42_{a42}_not_FY'] += 1
            continue
        # v8.8.2: Position 49 must be E, Q, or K (not G which is VH-like)
        if a49 not in ('E', 'Q', 'K'):
            removal_reasons[f'pos49_{a49}_not_EQK'] += 1
            continue
        
        filtered.append(c)
    
    # DO NOT re-add lead to candidates - main will add it at the end
    # Lead is kept separate and will be prepended to final_selected
    candidates = filtered
    print(f"  True VHH filter: {before_filter} → {len(candidates)} (removed {before_filter - len(candidates)})")
    
    # v8.8.2: Show scaffold breakdown after filter
    orig_count = sum(1 for c in candidates if getattr(c, 'scaffold_type', 'original') == 'original')
    univ_count = sum(1 for c in candidates if getattr(c, 'scaffold_type', 'original') == 'universal')
    print(f"  After filter scaffold breakdown: {orig_count} original, {univ_count} universal")
    
    if lead_candidate:
        print(f"  Lead candidate kept separately (not in selection pool)")
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
    
    # =========================================================================
    # v8.6: DEDUPLICATION (configurable)
    # =========================================================================
    if args.dedup_mode == 'within_track':
        print(f"\nv8.6: Deduplicating within track (preserves track intent)...")
        before_dedup = len(candidates)
        candidates = dedup_within_track(candidates)
        print(f"  Before: {before_dedup}, After: {len(candidates)}, Removed: {before_dedup - len(candidates)}")
    elif args.dedup_mode == 'global':
        print(f"\nGlobal deduplication (legacy mode)...")
        # Keep legacy global dedup for backward compatibility
        seq_to_best = {}
        for c in candidates:
            if c.sequence not in seq_to_best:
                seq_to_best[c.sequence] = c
            elif c.scoring.combined_score > seq_to_best[c.sequence].scoring.combined_score:
                seq_to_best[c.sequence] = c
        before_dedup = len(candidates)
        candidates = list(seq_to_best.values())
        print(f"  Before: {before_dedup}, After: {len(candidates)}, Removed: {before_dedup - len(candidates)}")
    else:
        print(f"\nv8.6: Skipping deduplication (dedup_mode=none)")
    
    # =========================================================================
    # v8.6: CYSTEINE VALIDATION 
    # =========================================================================
    candidates, cysteine_stats = filter_invalid_cysteines(candidates, verbose=True)
    
    # =========================================================================
    # v8.6: TRACK QUOTA SELECTION
    # =========================================================================
    print(f"\nv8.6: Applying per-track per-bucket selection (target: {args.n_select} total)")
    
    # Build track quotas from args
    track_quotas = {
        0: args.track0_quota,
        1: args.track1_quota,
        2: args.track2_quota,
        3: args.track3_quota,
        4: args.track4_quota,
    }
    print(f"  Track quotas per bucket: {track_quotas}")
    print(f"  Bucket target: {args.bucket_target}")
    
    # Add lead candidate first (always included)
    if lead_candidate:
        lead_candidate.rank = 0
        lead_candidate.is_lead = True
    
    # Use new selection function
    selected = select_with_track_quotas(
        candidates,
        n_select_total=args.n_select - (1 if lead_candidate else 0),  # Reserve 1 slot for lead
        per_bucket_track_quota=track_quotas,
        bucket_target=args.bucket_target,
        score_attr="combined_score"
    )
    
    if lead_candidate is not None:
        selected = [c for c in selected if not getattr(c, "is_lead", False)]

    # Assign ranks: controls first (with bracket notation), then ranked
    # Controls: [1], [2], [3]... (exempt from competitive ranking)
    # Ranked: 1, 2, 3... (competitive ranking by score)
    controls_in_selected = [c for c in selected if getattr(c, 'ranking_exempt', False)]
    ranked_in_selected = [c for c in selected if not getattr(c, 'ranking_exempt', False)]
    
    # Assign ranks to controls: they're exempt, so use bracket notation
    for i, c in enumerate(controls_in_selected):
        c.rank = i + 1
        c.rank_display = f"[{i + 1}]"  # Bracket indicates control/exempt
    
    # Assign ranks to ranked candidates
    for i, c in enumerate(ranked_in_selected):
        c.rank = len(controls_in_selected) + i + 1
        c.rank_display = str(c.rank)  # No brackets for ranked
    
    # Rebuild selected in order: controls first, then ranked
    selected = controls_in_selected + ranked_in_selected
    
    # Build final list: lead + controls + ranked
    if lead_candidate:
        lead_candidate.rank_display = "LEAD"
        final_selected = [lead_candidate] + selected
    else:
        final_selected = selected
    
    # Build selection stats for summary
    selection_stats = {
        'total_selected': len(final_selected),
        'n_controls': len(controls_in_selected),
        'n_ranked': len(ranked_in_selected),
        'track_quotas': track_quotas,
        'bucket_target': args.bucket_target,
        'dedup_mode': args.dedup_mode,
        'by_track': defaultdict(int),
        'by_bucket': defaultdict(int),
        'cysteine_stats': cysteine_stats,
    }
    
    for c in selected:
        t = track_number(c)
        selection_stats['by_track'][f'Track {t}'] += 1
        b = bucket_key(c)
        selection_stats['by_bucket'][str(b)] += 1
    
    selection_stats['by_track'] = dict(selection_stats['by_track'])
    selection_stats['by_bucket'] = dict(selection_stats['by_bucket'])
    
    # Build breakdown for ranked candidates only
    ranked_by_track = defaultdict(int)
    ranked_by_hallmark = defaultdict(int)
    ranked_by_scaffold = defaultdict(int)
    
    for c in ranked_in_selected:
        t = track_number(c)
        ranked_by_track[f'Track {t}'] += 1
        hm = get_hallmarks_from_positions(c.imgt_positions)
        ranked_by_hallmark[hm] += 1
        scaffold = getattr(c, 'scaffold_type', 'original')
        ranked_by_scaffold[scaffold] += 1
    
    # Build breakdown for controls
    controls_by_hallmark = defaultdict(int)
    for c in controls_in_selected:
        hm = get_hallmarks_from_positions(c.imgt_positions)
        controls_by_hallmark[hm] += 1
    
    print(f"\n{'='*70}")
    print("SELECTION SUMMARY")
    print(f"{'='*70}")
    print(f"\nTotal: {len(final_selected)} = 1 lead + {len(controls_in_selected)} controls + {len(ranked_in_selected)} ranked")
    
    print(f"\n  Controls ({len(controls_in_selected)} total):")
    print(f"    By hallmark: {dict(controls_by_hallmark)}")
    
    print(f"\n  Ranked ({len(ranked_in_selected)} total):")
    print(f"    By track: {dict(sorted(ranked_by_track.items()))}")
    print(f"    By hallmark: {dict(sorted(ranked_by_hallmark.items()))}")
    print(f"    By scaffold: {dict(ranked_by_scaffold)}")
    
    print(f"\nTop 15:")
    print("-" * 125)
    print(f"{'Rank':>6} {'Track':>18} {'Gen#':>4} {'ID':35} {'Score':>7} {'pLDDT':>6} {'Rules':>8} {'FW%':>5} {'Hallmarks':>9}")
    print("-" * 125)
    
    for c in final_selected[:15]:
        score_str = "LEAD" if c.is_lead else f"{c.scoring.combined_score:.3f}"
        rules_str = f"{c.scoring.rules_passed}/{c.scoring.rules_applicable}"
        hallmarks = get_hallmarks_from_positions(c.imgt_positions)
        track_str = getattr(c, 'track', c.design_track)[:18]
        rank_str = getattr(c, 'rank_display', str(c.rank))
        print(f"{rank_str:>6} {track_str:>18} {c.generation_order:4d} {c.id[:35]:35} {score_str:>7} {c.scoring.plddt_mean:6.1f} {rules_str:>8} {c.framework_identity_pct:5.1f} {hallmarks:>9}")
    
    # Add ranked breakdown to selection_stats for JSON summary
    selection_stats['ranked_by_track'] = dict(sorted(ranked_by_track.items()))
    selection_stats['ranked_by_hallmark'] = dict(sorted(ranked_by_hallmark.items()))
    selection_stats['ranked_by_scaffold'] = dict(ranked_by_scaffold)
    selection_stats['controls_by_hallmark'] = dict(controls_by_hallmark)
    
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