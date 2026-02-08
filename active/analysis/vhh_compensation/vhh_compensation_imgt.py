#!/usr/bin/env python3
"""
VHH Compensation Analysis - IMGT Positions
===========================================

Analyzes CDR→FR correlations using true IMGT positions from the annotated database.
Outputs rules in JSON format compatible with vhh_designer_v7.

Input: vhh_full_annotated_v3.csv (12M sequences with FR/CDR columns)
Output: compensation_imgt_rules.json

Key difference from old compensation analysis:
- OLD: FR2[3], FR3[6] (substring positions, dataset-specific)
- NEW: IMGT42, IMGT71 (universal IMGT positions)

Usage:
    python vhh_compensation_imgt.py \\
        --csv data/databases/annotated/vhh_full_annotated_v3.csv \\
        --output models/compensation/imgt_v1 \\
        --chunk-size 500000

Author: Claude (Anthropic)
Date: January 2026
"""

import os
import sys
import json
import argparse
import numpy as np
import pandas as pd
from collections import defaultdict, Counter
from typing import Dict, List, Tuple, Optional
from datetime import datetime
import gc

# ============================================================
# IMGT POSITION MAPPING
# ============================================================

# IMGT boundaries for VHH
IMGT_BOUNDARIES = {
    'FR1': (1, 26),
    'CDR1': (27, 38),
    'FR2': (39, 55),
    'CDR2': (56, 65),
    'FR3': (66, 104),
    'CDR3': (105, 117),
    'FR4': (118, 128),
}

# Key positions for analysis
VERNIER_POSITIONS = [42, 49, 50, 52, 66, 67, 68, 69, 71, 76, 78, 82, 87, 89, 91, 94]
HALLMARK_POSITIONS = [42, 49, 50, 52]

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

def residue_category(aa: str) -> str:
    if aa in 'AILMFVW': return 'hydrophobic'
    elif aa in 'STNQ': return 'polar'
    elif aa in 'DE': return 'acidic'
    elif aa in 'KRH': return 'basic'
    elif aa in 'GP': return 'special'
    elif aa in 'CY': return 'reactive'
    else: return 'other'

def get_cdr_features(cdr: str) -> dict:
    """Extract features from a CDR sequence."""
    if not cdr:
        return {'length': 0, 'charge': 0, 'hydrophobicity': 0}
    
    length = len(cdr)
    charge = sum(CHARGE.get(aa, 0) for aa in cdr)
    hydro = sum(HYDROPHOBICITY.get(aa, 0) for aa in cdr) / max(length, 1)
    
    return {
        'length': length,
        'charge': charge,
        'hydrophobicity': round(hydro, 2),
        'n_cys': cdr.count('C'),
    }

def categorize_cdr3(cdr3: str) -> dict:
    """Categorize CDR3 into bins."""
    features = get_cdr_features(cdr3)
    
    # Length categories
    length = features['length']
    if length <= 8:
        len_cat = 'short'
    elif length <= 14:
        len_cat = 'medium'
    else:
        len_cat = 'long'
    
    # Charge categories
    charge = features['charge']
    if charge <= -2:
        charge_cat = 'negative'
    elif charge >= 2:
        charge_cat = 'positive'
    else:
        charge_cat = 'neutral'
    
    # Cysteine categories
    n_cys = features['n_cys']
    if n_cys == 0:
        cys_cat = 'no_cys'
    elif n_cys == 1:
        cys_cat = 'one_cys'
    else:
        cys_cat = 'multi_cys'
    
    return {
        'cdr3_len': len_cat,
        'cdr3_charge': charge_cat,
        'cdr3_cys': cys_cat,
        **features,
    }

# ============================================================
# BUILD IMGT POSITION DICT
# ============================================================

def build_imgt_dict(fr1: str, cdr1: str, fr2: str, cdr2: str, 
                    fr3: str, cdr3: str, fr4: str) -> dict:
    """Build IMGT position -> amino acid mapping."""
    imgt = {}
    
    regions = [
        ('FR1', fr1, 1),
        ('CDR1', cdr1, 27),
        ('FR2', fr2, 39),
        ('CDR2', cdr2, 56),
        ('FR3', fr3, 66),
        ('CDR3', cdr3, 105),
        ('FR4', fr4, 118),
    ]
    
    for region_name, seq, start_pos in regions:
        if seq:
            for i, aa in enumerate(seq):
                imgt[start_pos + i] = aa
    
    return imgt

# ============================================================
# STREAMING ANALYSIS
# ============================================================

class StreamingCompensationAnalyzer:
    """
    Memory-efficient CDR→FR correlation analysis.
    
    For each FR position (IMGT), tracks:
    - Which amino acid appears there
    - Conditioned on CDR features (length, charge, specific residues)
    """
    
    def __init__(self):
        # Main storage: {imgt_pos: {condition: {aa: count}}}
        # condition examples: "cdr3_len=short", "cdr1[3]=S", "family=F_C2"
        self.counts = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
        
        # Track totals for confidence calculation
        self.condition_totals = defaultdict(lambda: defaultdict(int))
        
        # Track position totals
        self.position_totals = defaultdict(lambda: defaultdict(int))
        
        self.total_sequences = 0
    
    def update(self, imgt_dict: dict, cdr1: str, cdr2: str, cdr3: str, 
               family: str, hallmarks: str):
        """Update counts with one sequence."""
        self.total_sequences += 1
        
        # Get CDR3 categories
        cdr3_cats = categorize_cdr3(cdr3)
        
        # Build conditions
        conditions = [
            f"cdr3_len={cdr3_cats['cdr3_len']}",
            f"cdr3_charge={cdr3_cats['cdr3_charge']}",
            f"cdr3_cys={cdr3_cats['cdr3_cys']}",
            f"family={family}",
            f"hallmarks={hallmarks}",
        ]
        
        # Add specific CDR residue conditions (key positions only)
        for i, aa in enumerate(cdr1[:5]):  # First 5 of CDR1
            conditions.append(f"cdr1[{i}]={aa}")
        for i, aa in enumerate(cdr2[:3]):  # First 3 of CDR2
            conditions.append(f"cdr2[{i}]={aa}")
        for i, aa in enumerate(cdr3[:3]):  # First 3 of CDR3
            conditions.append(f"cdr3[{i}]={aa}")
        if len(cdr3) >= 3:
            for i, aa in enumerate(cdr3[-3:]):  # Last 3 of CDR3
                conditions.append(f"cdr3[-{3-i}]={aa}")
        
        # For each FR position, update counts
        for imgt_pos, aa in imgt_dict.items():
            if not isinstance(imgt_pos, int):
                continue
            
            # Only FR positions
            if not any(start <= imgt_pos <= end for region, (start, end) in IMGT_BOUNDARIES.items() 
                      if region.startswith('FR')):
                continue
            
            # Update position total
            self.position_totals[imgt_pos][aa] += 1
            
            # Update conditional counts
            for cond in conditions:
                self.counts[imgt_pos][cond][aa] += 1
                self.condition_totals[imgt_pos][cond] += 1
    
    def extract_rules(self, min_support: int = 1000, min_confidence: float = 0.7,
                      min_lift: float = 1.2) -> List[dict]:
        """
        Extract rules where a condition predicts an FR residue.
        
        Rule format:
        {
            "condition": "cdr3_len=short",
            "position": "IMGT71",
            "suggested_aa": "V",
            "confidence": 0.85,
            "support": 45000,
            "lift": 1.5,
            "source": "compensation_imgt_v1"
        }
        """
        rules = []
        
        for imgt_pos, cond_dict in self.counts.items():
            # Get baseline (unconditional) distribution for this position
            baseline = self.position_totals[imgt_pos]
            baseline_total = sum(baseline.values())
            
            if baseline_total < min_support:
                continue
            
            for condition, aa_counts in cond_dict.items():
                cond_total = self.condition_totals[imgt_pos][condition]
                
                if cond_total < min_support:
                    continue
                
                for aa, count in aa_counts.items():
                    confidence = count / cond_total
                    
                    if confidence < min_confidence:
                        continue
                    
                    if count < min_support:
                        continue
                    
                    # Calculate lift (how much better than baseline)
                    baseline_prob = baseline.get(aa, 0) / baseline_total
                    if baseline_prob > 0:
                        lift = confidence / baseline_prob
                    else:
                        lift = float('inf')
                    
                    if lift < min_lift:
                        continue
                    
                    rules.append({
                        'condition': condition,
                        'position': f'IMGT{imgt_pos}',
                        'suggested_aa': aa,
                        'confidence': round(confidence, 3),
                        'support': count,
                        'lift': round(lift, 2),
                        'baseline_prob': round(baseline_prob, 3),
                        'source': 'compensation_imgt_v1',
                    })
        
        # Sort by support (highest first)
        rules.sort(key=lambda x: -x['support'])
        
        return rules
    
    def get_vernier_archetypes(self, min_support: int = 10000) -> dict:
        """
        Get dominant residues at Vernier positions by family/hallmarks.
        """
        archetypes = defaultdict(dict)
        
        for imgt_pos in VERNIER_POSITIONS:
            if imgt_pos not in self.counts:
                continue
            
            for condition, aa_counts in self.counts[imgt_pos].items():
                # Only family and hallmark conditions
                if not (condition.startswith('family=') or condition.startswith('hallmarks=')):
                    continue
                
                total = sum(aa_counts.values())
                if total < min_support:
                    continue
                
                # Find dominant residue
                dominant_aa, dominant_count = max(aa_counts.items(), key=lambda x: x[1])
                
                archetypes[condition][f'IMGT{imgt_pos}'] = {
                    'aa': dominant_aa,
                    'count': dominant_count,
                    'frequency': round(dominant_count / total, 3),
                }
        
        return dict(archetypes)

# ============================================================
# MAIN PROCESSING
# ============================================================

def process_csv_streaming(csv_path: str, output_dir: str, chunk_size: int = 500000):
    """Process the annotated CSV in streaming fashion."""
    
    print(f"Processing: {csv_path}")
    print(f"Chunk size: {chunk_size:,}")
    
    os.makedirs(output_dir, exist_ok=True)
    
    analyzer = StreamingCompensationAnalyzer()
    
    # Process in chunks
    chunk_num = 0
    total_processed = 0
    
    for chunk in pd.read_csv(csv_path, chunksize=chunk_size, low_memory=False):
        chunk_num += 1
        chunk_start = (chunk_num - 1) * chunk_size
        
        print(f"\nChunk {chunk_num}: rows {chunk_start:,} - {chunk_start + len(chunk):,}")
        
        for idx, row in chunk.iterrows():
            try:
                # Get regions
                fr1 = str(row.get('fr1', '')) if pd.notna(row.get('fr1')) else ''
                cdr1 = str(row.get('cdr1', '')) if pd.notna(row.get('cdr1')) else ''
                fr2 = str(row.get('fr2', '')) if pd.notna(row.get('fr2')) else ''
                cdr2 = str(row.get('cdr2', '')) if pd.notna(row.get('cdr2')) else ''
                fr3 = str(row.get('fr3', '')) if pd.notna(row.get('fr3')) else ''
                cdr3 = str(row.get('cdr3', '')) if pd.notna(row.get('cdr3')) else ''
                fr4 = str(row.get('fr4', '')) if pd.notna(row.get('fr4')) else ''
                
                family = str(row.get('family', 'Unknown')) if pd.notna(row.get('family')) else 'Unknown'
                hallmarks = str(row.get('hallmarks', '----')) if pd.notna(row.get('hallmarks')) else '----'
                
                # Skip if missing key regions
                if not cdr3 or not fr3:
                    continue
                
                # Build IMGT dict
                imgt_dict = build_imgt_dict(fr1, cdr1, fr2, cdr2, fr3, cdr3, fr4)
                
                # Update analyzer
                analyzer.update(imgt_dict, cdr1, cdr2, cdr3, family, hallmarks)
                
                total_processed += 1
                
            except Exception as e:
                continue
        
        print(f"  Processed: {total_processed:,} sequences")
        gc.collect()
    
    print(f"\n{'='*60}")
    print(f"ANALYSIS COMPLETE")
    print(f"{'='*60}")
    print(f"Total sequences: {analyzer.total_sequences:,}")
    
    # Extract rules
    print("\nExtracting rules...")
    rules = analyzer.extract_rules(min_support=5000, min_confidence=0.70, min_lift=1.1)
    print(f"  Found {len(rules)} rules (support >= 5000, confidence >= 70%, lift >= 1.1)")
    
    # Get vernier archetypes
    print("\nExtracting Vernier archetypes...")
    archetypes = analyzer.get_vernier_archetypes(min_support=10000)
    print(f"  Found archetypes for {len(archetypes)} conditions")
    
    # Save rules as JSON
    rules_path = os.path.join(output_dir, 'compensation_imgt_rules.json')
    with open(rules_path, 'w') as f:
        json.dump(rules, f, indent=2)
    print(f"\nSaved: {rules_path}")
    
    # Save archetypes as JSON
    archetypes_path = os.path.join(output_dir, 'vernier_archetypes_imgt.json')
    with open(archetypes_path, 'w') as f:
        json.dump(archetypes, f, indent=2)
    print(f"Saved: {archetypes_path}")
    
    # Print summary
    print(f"\n{'='*60}")
    print("TOP RULES BY SUPPORT")
    print(f"{'='*60}")
    for r in rules[:20]:
        print(f"  {r['position']:8s} -> {r['suggested_aa']} | {r['condition']:30s} | "
              f"conf={r['confidence']:.2f} lift={r['lift']:.1f} n={r['support']:,}")
    
    # Summary by position
    print(f"\n{'='*60}")
    print("RULES BY IMGT POSITION")
    print(f"{'='*60}")
    by_pos = defaultdict(list)
    for r in rules:
        by_pos[r['position']].append(r)
    
    for pos in sorted(by_pos.keys(), key=lambda x: int(x.replace('IMGT', ''))):
        pos_rules = by_pos[pos]
        print(f"  {pos}: {len(pos_rules)} rules")
        for r in pos_rules[:3]:
            print(f"    -> {r['suggested_aa']} ({r['condition'][:25]}...) conf={r['confidence']:.2f}")
    
    # Summary by condition type
    print(f"\n{'='*60}")
    print("RULES BY CONDITION TYPE")
    print(f"{'='*60}")
    by_type = defaultdict(int)
    for r in rules:
        cond = r['condition']
        if '=' in cond:
            cond_type = cond.split('=')[0]
            by_type[cond_type] += 1
    
    for cond_type, count in sorted(by_type.items(), key=lambda x: -x[1]):
        print(f"  {cond_type}: {count} rules")
    
    return rules, archetypes


def main():
    parser = argparse.ArgumentParser(description='VHH Compensation Analysis - IMGT Positions')
    parser.add_argument('--csv', required=True, help='Annotated CSV file (vhh_full_annotated_v3.csv)')
    parser.add_argument('--output', '-o', default='models/compensation/imgt_v1', help='Output directory')
    parser.add_argument('--chunk-size', type=int, default=500000, help='Chunk size for streaming')
    
    args = parser.parse_args()
    
    print("=" * 70)
    print("VHH COMPENSATION ANALYSIS - IMGT POSITIONS")
    print("=" * 70)
    print(f"Input: {args.csv}")
    print(f"Output: {args.output}")
    print(f"Chunk size: {args.chunk_size:,}")
    print()
    
    start_time = datetime.now()
    
    rules, archetypes = process_csv_streaming(args.csv, args.output, args.chunk_size)
    
    elapsed = datetime.now() - start_time
    print(f"\nTotal time: {elapsed}")
    print("Done!")


if __name__ == '__main__':
    main()
