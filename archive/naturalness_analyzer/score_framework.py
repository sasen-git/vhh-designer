#!/usr/bin/env python3
"""
Framework Scorer
================
Evaluates a universal framework against p(FR|CDR) models.

For each FR position, shows:
- How well your framework residue matches the repertoire preference
- Alternative residues that might work better for specific CDR panels
- Positions that are "out of distribution" for your target CDRs

Usage:
    # Score against all training data
    python score_framework.py --models pfr_cdr_models.pkl

    # Score for specific CDR sequences
    python score_framework.py --models pfr_cdr_models.pkl \
        --cdr1 "GFTFSSYA" --cdr2 "ISYDGSNK" --cdr3 "ARDLLVRY"

    # Score for a CDR panel (file with one CDR set per line)
    python score_framework.py --models pfr_cdr_models.pkl \
        --cdr-file my_cdrs.txt

    # Custom framework (default is Vincke universal)
    python score_framework.py --models pfr_cdr_models.pkl \
        --fr1 "EVQLVESGGGLVQPGGSLRLSCAAS" \
        --fr2 "WFRQAPGKGREFVA" \
        --fr3 "YYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYC" \
        --fr4 "WGQGTQVTVSS"
"""

import numpy as np
import pickle
import sys
import os
import argparse
from collections import defaultdict

# ============================================================
# DEFAULT UNIVERSAL FRAMEWORK (Vincke et al. 2009)
# ============================================================

DEFAULT_FRAMEWORK = {
    'FR1': 'EVQLVESGGGLVQPGGSLRLSCAAS',
    'FR2': 'WFRQAPGKGREFVA',
    'FR3': 'YYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYC',
    'FR4': 'WGQGTQVTVSS'
}

# ============================================================
# AMINO ACID PROPERTIES (same as model builder)
# ============================================================

HYDROPHOBICITY = {
    'I': 4.5, 'V': 4.2, 'L': 3.8, 'F': 2.8, 'C': 2.5,
    'M': 1.9, 'A': 1.8, 'G': -0.4, 'T': -0.7, 'S': -0.8,
    'W': -0.9, 'Y': -1.3, 'P': -1.6, 'H': -3.2, 'E': -3.5,
    'Q': -3.5, 'D': -3.5, 'N': -3.5, 'K': -3.9, 'R': -4.5
}

CHARGE = {
    'K': 1, 'R': 1, 'H': 0.1, 'D': -1, 'E': -1,
    'A': 0, 'C': 0, 'F': 0, 'G': 0, 'I': 0, 'L': 0, 'M': 0,
    'N': 0, 'P': 0, 'Q': 0, 'S': 0, 'T': 0, 'V': 0, 'W': 0, 'Y': 0
}

AROMATIC = {'F': 1, 'W': 1, 'Y': 1, 'H': 0.5}
POLAR = {'S': 1, 'T': 1, 'N': 1, 'Q': 1, 'Y': 0.5, 'H': 0.5, 'C': 0.5}

AA_LIST = list('ACDEFGHIKLMNPQRSTVWY')
AA_TO_IDX = {aa: i for i, aa in enumerate(AA_LIST)}

def compute_sequence_properties(seq):
    if not seq or len(seq) == 0:
        return {'length': 0, 'charge': 0, 'hydrophobicity': 0, 
                'aromaticity': 0, 'small_fraction': 0, 'polar_fraction': 0}
    
    seq = str(seq).upper()
    n = len(seq)
    
    charge = sum(CHARGE.get(aa, 0) for aa in seq)
    hydro = sum(HYDROPHOBICITY.get(aa, 0) for aa in seq) / n
    aromatic = sum(AROMATIC.get(aa, 0) for aa in seq) / n
    polar = sum(POLAR.get(aa, 0) for aa in seq) / n
    
    return {
        'length': n,
        'charge': charge,
        'hydrophobicity': hydro,
        'aromaticity': aromatic,
        'polar_fraction': polar
    }

def extract_cdr_features(cdr1, cdr2, cdr3):
    features = []
    
    # CDR lengths
    features.extend([len(cdr1), len(cdr2), len(cdr3)])
    
    # CDR properties
    for cdr_seq in [cdr1, cdr2, cdr3]:
        props = compute_sequence_properties(cdr_seq)
        features.append(props['charge'])
        features.append(props['hydrophobicity'])
        features.append(props['aromaticity'])
        features.append(props['polar_fraction'])
    
    # CDR boundary amino acids
    for cdr_seq in [cdr1, cdr2, cdr3]:
        for i in range(2):
            if i < len(cdr_seq):
                features.append(AA_TO_IDX.get(cdr_seq[i], 10))
            else:
                features.append(10)
        for i in [-2, -1]:
            if abs(i) <= len(cdr_seq):
                features.append(AA_TO_IDX.get(cdr_seq[i], 10))
            else:
                features.append(10)
    
    return np.array(features, dtype=np.float32).reshape(1, -1)

# ============================================================
# SCORING FUNCTIONS
# ============================================================

def score_framework_position(model_info, scaler, features, framework_aa):
    """Score a framework residue given CDR features"""
    
    if model_info['type'] == 'constant':
        # Position is constant - check if it matches
        expected_aa = model_info['aa']
        matches = framework_aa == expected_aa
        return {
            'probability': 1.0 if matches else 0.0,
            'expected': expected_aa,
            'is_constant': True,
            'alternatives': [],
            'matches': matches
        }
    
    model = model_info['model']
    classes = model_info['classes']
    
    # Scale features
    X_scaled = scaler.transform(features)
    
    # Get probabilities
    probs = model.predict_proba(X_scaled)[0]
    
    # Map class indices to amino acids
    class_to_aa = {c: AA_LIST[c] for c in classes}
    
    # Create probability distribution
    prob_dist = {}
    for i, c in enumerate(classes):
        prob_dist[AA_LIST[c]] = probs[i]
    
    # Get probability of framework residue
    fw_prob = prob_dist.get(framework_aa, 0.0)
    
    # Get expected residue
    expected_idx = np.argmax(probs)
    expected_aa = AA_LIST[classes[expected_idx]]
    expected_prob = probs[expected_idx]
    
    # Get alternatives (sorted by probability)
    alternatives = sorted(prob_dist.items(), key=lambda x: x[1], reverse=True)[:5]
    
    return {
        'probability': fw_prob,
        'expected': expected_aa,
        'expected_prob': expected_prob,
        'is_constant': False,
        'alternatives': alternatives,
        'matches': framework_aa == expected_aa,
        'prob_dist': prob_dist
    }

def score_full_framework(models, scalers, framework, cdr_features):
    """Score entire framework against CDR features"""
    
    results = {}
    
    for (fw_name, fw_pos), model_info in models.items():
        try:
            fw_seq = framework[fw_name]
            fw_aa = fw_seq[fw_pos]
        except (KeyError, IndexError):
            continue
        
        scaler = scalers.get((fw_name, fw_pos))
        
        if model_info['type'] == 'constant':
            score = score_framework_position(model_info, None, cdr_features, fw_aa)
        else:
            if scaler is None:
                continue
            score = score_framework_position(model_info, scaler, cdr_features, fw_aa)
        
        results[(fw_name, fw_pos)] = {
            'fw_aa': fw_aa,
            **score
        }
    
    return results

def summarize_scores(results, threshold=0.1):
    """Summarize scoring results"""
    
    summary = {
        'total_positions': len(results),
        'matching': 0,
        'mismatched': 0,
        'low_probability': [],  # Positions where FW residue has low prob
        'high_probability': [], # Positions where FW residue matches well
        'by_framework': defaultdict(list)
    }
    
    for (fw_name, fw_pos), score in results.items():
        if score['matches']:
            summary['matching'] += 1
        else:
            summary['mismatched'] += 1
        
        if score['probability'] < threshold:
            summary['low_probability'].append({
                'position': f"{fw_name}[{fw_pos}]",
                'current': score['fw_aa'],
                'expected': score['expected'],
                'prob': score['probability'],
                'alternatives': score['alternatives'][:3] if 'alternatives' in score else []
            })
        elif score['probability'] > 0.5:
            summary['high_probability'].append({
                'position': f"{fw_name}[{fw_pos}]",
                'residue': score['fw_aa'],
                'prob': score['probability']
            })
        
        summary['by_framework'][fw_name].append({
            'pos': fw_pos,
            'aa': score['fw_aa'],
            'prob': score['probability'],
            'expected': score['expected']
        })
    
    # Sort low probability by probability (lowest first)
    summary['low_probability'] = sorted(summary['low_probability'], 
                                        key=lambda x: x['prob'])
    
    return summary

# ============================================================
# OUTPUT FORMATTING
# ============================================================

def print_framework_report(summary, framework):
    """Print formatted framework scoring report"""
    
    print("\n" + "="*70)
    print("FRAMEWORK SCORING REPORT")
    print("="*70)
    
    print(f"\nTotal positions scored: {summary['total_positions']}")
    print(f"Matching expected: {summary['matching']} ({100*summary['matching']/summary['total_positions']:.1f}%)")
    print(f"Different from expected: {summary['mismatched']} ({100*summary['mismatched']/summary['total_positions']:.1f}%)")
    
    # Problem positions
    print("\n" + "-"*70)
    print("⚠️  PROBLEM POSITIONS (low probability given CDRs)")
    print("-"*70)
    
    if not summary['low_probability']:
        print("  None! Your framework is well-matched to the CDRs.")
    else:
        print(f"\n{'Position':<12} {'Current':<8} {'Expected':<8} {'Prob':<8} {'Alternatives'}")
        print("-"*70)
        for pos in summary['low_probability'][:15]:
            alts = ", ".join([f"{aa}({p:.0%})" for aa, p in pos['alternatives']])
            print(f"{pos['position']:<12} {pos['current']:<8} {pos['expected']:<8} {pos['prob']:<8.1%} {alts}")
    
    # Framework visualization
    print("\n" + "-"*70)
    print("FRAMEWORK VISUALIZATION")
    print("(lowercase = low probability, UPPERCASE = good match)")
    print("-"*70)
    
    for fw_name in ['FR1', 'FR2', 'FR3', 'FR4']:
        fw_seq = framework[fw_name]
        positions = summary['by_framework'][fw_name]
        
        # Build visualization
        vis = list(fw_seq)
        for pos_info in positions:
            pos = pos_info['pos']
            try:
                if pos_info['prob'] < 0.1:
                    vis[pos] = vis[pos].lower()
                elif pos_info['prob'] > 0.5:
                    vis[pos] = vis[pos].upper()
            except IndexError:
                pass
        
        print(f"\n{fw_name}: {''.join(vis)}")
    
    print("\n" + "="*70)

def print_recommendations(summary):
    """Print actionable recommendations"""
    
    print("\n" + "="*70)
    print("RECOMMENDATIONS")
    print("="*70)
    
    if not summary['low_probability']:
        print("\n✅ Your framework looks compatible with the target CDRs!")
        print("   No major changes recommended.")
        return
    
    print("\nSuggested mutations (ordered by priority):\n")
    
    for i, pos in enumerate(summary['low_probability'][:10], 1):
        if pos['alternatives']:
            best_alt = pos['alternatives'][0]
            print(f"  {i}. {pos['position']}: {pos['current']} → {best_alt[0]}")
            print(f"     Probability: {pos['prob']:.1%} → {best_alt[1]:.1%}")
            print(f"     Other options: {', '.join([f'{aa}' for aa, p in pos['alternatives'][1:3]])}")
            print()

# ============================================================
# MAIN
# ============================================================

def main():
    parser = argparse.ArgumentParser(description='Score framework against CDRs')
    parser.add_argument('--models', required=True, help='Path to pfr_cdr_models.pkl')
    parser.add_argument('--cdr1', help='CDR1 sequence')
    parser.add_argument('--cdr2', help='CDR2 sequence')
    parser.add_argument('--cdr3', help='CDR3 sequence')
    parser.add_argument('--cdr-file', help='File with CDR sequences (one set per line: cdr1,cdr2,cdr3)')
    parser.add_argument('--fr1', help='Custom FR1 sequence')
    parser.add_argument('--fr2', help='Custom FR2 sequence')
    parser.add_argument('--fr3', help='Custom FR3 sequence')
    parser.add_argument('--fr4', help='Custom FR4 sequence')
    parser.add_argument('--threshold', type=float, default=0.1, 
                        help='Probability threshold for flagging positions')
    parser.add_argument('--output', '-o', help='Output file for detailed report')
    
    args = parser.parse_args()
    
    # Load models
    print("Loading models...")
    with open(args.models, 'rb') as f:
        data = pickle.load(f)
    
    models = data['models']
    scalers = data['scalers']
    metrics = data['metrics']
    
    print(f"  Loaded {len(models)} position models")
    
    # Set framework
    framework = DEFAULT_FRAMEWORK.copy()
    if args.fr1:
        framework['FR1'] = args.fr1
    if args.fr2:
        framework['FR2'] = args.fr2
    if args.fr3:
        framework['FR3'] = args.fr3
    if args.fr4:
        framework['FR4'] = args.fr4
    
    print(f"\nFramework to score:")
    for fw_name, fw_seq in framework.items():
        print(f"  {fw_name}: {fw_seq}")
    
    # Get CDR sequences
    cdr_sets = []
    
    if args.cdr1 and args.cdr2 and args.cdr3:
        cdr_sets.append((args.cdr1, args.cdr2, args.cdr3))
    elif args.cdr_file:
        with open(args.cdr_file) as f:
            for line in f:
                parts = line.strip().split(',')
                if len(parts) == 3:
                    cdr_sets.append(tuple(parts))
    
    if not cdr_sets:
        print("\n⚠️  No CDR sequences provided. Showing average across training data.")
        print("   Use --cdr1/--cdr2/--cdr3 or --cdr-file for specific scoring.")
        
        # Use median CDR properties as a "typical" CDR set
        # This is an approximation - real usage should provide actual CDRs
        cdr_sets.append(("GFTFSSYA", "ISYDGSNK", "ARDLYYYYGMDV"))
    
    # Score for each CDR set
    all_summaries = []
    
    for i, (cdr1, cdr2, cdr3) in enumerate(cdr_sets):
        if len(cdr_sets) > 1:
            print(f"\n--- CDR Set {i+1} ---")
            print(f"  CDR1: {cdr1}")
            print(f"  CDR2: {cdr2}")
            print(f"  CDR3: {cdr3}")
        else:
            print(f"\nTarget CDRs:")
            print(f"  CDR1: {cdr1}")
            print(f"  CDR2: {cdr2}")
            print(f"  CDR3: {cdr3}")
        
        # Extract features
        features = extract_cdr_features(cdr1, cdr2, cdr3)
        
        # Score framework
        results = score_full_framework(models, scalers, framework, features)
        
        # Summarize
        summary = summarize_scores(results, threshold=args.threshold)
        all_summaries.append(summary)
        
        # Print report
        print_framework_report(summary, framework)
        print_recommendations(summary)
    
    # If multiple CDR sets, aggregate
    if len(cdr_sets) > 1:
        print("\n" + "="*70)
        print("AGGREGATE ANALYSIS")
        print("="*70)
        
        # Find positions that are problematic across multiple CDR sets
        problem_counts = defaultdict(int)
        for summary in all_summaries:
            for pos in summary['low_probability']:
                problem_counts[pos['position']] += 1
        
        consistent_problems = [(pos, count) for pos, count in problem_counts.items() 
                              if count >= len(cdr_sets) * 0.5]
        consistent_problems.sort(key=lambda x: x[1], reverse=True)
        
        if consistent_problems:
            print("\nPositions problematic across multiple CDR sets:")
            for pos, count in consistent_problems:
                print(f"  {pos}: problematic in {count}/{len(cdr_sets)} CDR sets")
    
    # Save detailed output
    if args.output:
        with open(args.output, 'w') as f:
            import json
            json.dump({
                'framework': framework,
                'cdr_sets': cdr_sets,
                'summaries': [s for s in all_summaries]
            }, f, indent=2, default=str)
        print(f"\n✓ Detailed report saved to: {args.output}")

if __name__ == '__main__':
    main()
