#!/usr/bin/env python3
"""
Conditional p(FR | CDR) Model Builder
=====================================
Trains predictive models for each framework position given CDR features.

For each FR position, answers:
- "Given these CDRs, what amino acid does nature choose here?"
- "How 'weird' is residue X at this position for these CDRs?"

Features used:
- CDR lengths (CDR1, CDR2, CDR3)
- CDR physicochemical properties (charge, hydrophobicity, aromaticity)
- CDR boundary amino acids (first/last 3 of each CDR)

Usage:
    # Build models from NPZ files
    python build_pfr_cdr_models.py /path/to/shards/*.npz --output ./models
    
    # Then use score_framework.py to evaluate your framework

Output:
    - pfr_cdr_models.pkl: Trained models for each FR position
    - model_summary.txt: Feature importances and model quality metrics
"""

import numpy as np
from collections import defaultdict, Counter
import pickle
import sys
import os
import time
import argparse
from glob import glob
import warnings
warnings.filterwarnings('ignore')

# ============================================================
# AMINO ACID PROPERTIES
# ============================================================

# Hydrophobicity (Kyte-Doolittle scale)
HYDROPHOBICITY = {
    'I': 4.5, 'V': 4.2, 'L': 3.8, 'F': 2.8, 'C': 2.5,
    'M': 1.9, 'A': 1.8, 'G': -0.4, 'T': -0.7, 'S': -0.8,
    'W': -0.9, 'Y': -1.3, 'P': -1.6, 'H': -3.2, 'E': -3.5,
    'Q': -3.5, 'D': -3.5, 'N': -3.5, 'K': -3.9, 'R': -4.5
}

# Charge at pH 7
CHARGE = {
    'K': 1, 'R': 1, 'H': 0.1,  # Positive
    'D': -1, 'E': -1,          # Negative
    'A': 0, 'C': 0, 'F': 0, 'G': 0, 'I': 0, 'L': 0, 'M': 0,
    'N': 0, 'P': 0, 'Q': 0, 'S': 0, 'T': 0, 'V': 0, 'W': 0, 'Y': 0
}

# Aromatic residues
AROMATIC = {'F': 1, 'W': 1, 'Y': 1, 'H': 0.5}

# Small residues
SMALL = {'G': 1, 'A': 1, 'S': 1, 'C': 1, 'T': 0.5, 'P': 0.5}

# Polar residues
POLAR = {'S': 1, 'T': 1, 'N': 1, 'Q': 1, 'Y': 0.5, 'H': 0.5, 'C': 0.5}

# All amino acids
AA_LIST = list('ACDEFGHIKLMNPQRSTVWY')
AA_TO_IDX = {aa: i for i, aa in enumerate(AA_LIST)}

def compute_sequence_properties(seq):
    """Compute physicochemical properties of a sequence"""
    if not seq or len(seq) == 0:
        return {'length': 0, 'charge': 0, 'hydrophobicity': 0, 
                'aromaticity': 0, 'small_fraction': 0, 'polar_fraction': 0}
    
    seq = str(seq).upper()
    n = len(seq)
    
    charge = sum(CHARGE.get(aa, 0) for aa in seq)
    hydro = sum(HYDROPHOBICITY.get(aa, 0) for aa in seq) / n
    aromatic = sum(AROMATIC.get(aa, 0) for aa in seq) / n
    small = sum(SMALL.get(aa, 0) for aa in seq) / n
    polar = sum(POLAR.get(aa, 0) for aa in seq) / n
    
    return {
        'length': n,
        'charge': charge,
        'hydrophobicity': hydro,
        'aromaticity': aromatic,
        'small_fraction': small,
        'polar_fraction': polar
    }

def aa_to_onehot(aa):
    """Convert amino acid to one-hot vector"""
    vec = [0] * 20
    if aa in AA_TO_IDX:
        vec[AA_TO_IDX[aa]] = 1
    return vec

# ============================================================
# FEATURE EXTRACTION
# ============================================================

def extract_cdr_features(cdr1, cdr2, cdr3):
    """Extract feature vector from CDR sequences"""
    
    features = []
    feature_names = []
    
    # CDR lengths
    features.extend([len(cdr1), len(cdr2), len(cdr3)])
    feature_names.extend(['cdr1_len', 'cdr2_len', 'cdr3_len'])
    
    # CDR properties
    for cdr_name, cdr_seq in [('cdr1', cdr1), ('cdr2', cdr2), ('cdr3', cdr3)]:
        props = compute_sequence_properties(cdr_seq)
        features.append(props['charge'])
        features.append(props['hydrophobicity'])
        features.append(props['aromaticity'])
        features.append(props['polar_fraction'])
        feature_names.extend([
            f'{cdr_name}_charge', f'{cdr_name}_hydro', 
            f'{cdr_name}_aromatic', f'{cdr_name}_polar'
        ])
    
    # CDR boundary amino acids (simplified - just use indices)
    # First and last 2 AAs of each CDR
    for cdr_name, cdr_seq in [('cdr1', cdr1), ('cdr2', cdr2), ('cdr3', cdr3)]:
        # First 2 positions
        for i in range(2):
            if i < len(cdr_seq):
                aa = cdr_seq[i]
                features.append(AA_TO_IDX.get(aa, 10))  # Use index, default to middle
            else:
                features.append(10)
            feature_names.append(f'{cdr_name}_pos{i}')
        
        # Last 2 positions
        for i in [-2, -1]:
            if abs(i) <= len(cdr_seq):
                aa = cdr_seq[i]
                features.append(AA_TO_IDX.get(aa, 10))
            else:
                features.append(10)
            feature_names.append(f'{cdr_name}_pos{i}')
    
    return np.array(features, dtype=np.float32), feature_names

# ============================================================
# DATA COLLECTION
# ============================================================

def extract_frameworks(full_seq, cdr1, cdr2, cdr3):
    """Extract framework regions"""
    if not full_seq or not cdr1 or not cdr2 or not cdr3:
        return None
    
    full_seq, cdr1, cdr2, cdr3 = str(full_seq), str(cdr1), str(cdr2), str(cdr3)
    
    cdr1_start = full_seq.find(cdr1)
    cdr2_start = full_seq.find(cdr2)
    cdr3_start = full_seq.find(cdr3)
    
    if cdr1_start == -1 or cdr2_start == -1 or cdr3_start == -1:
        return None
    if not (cdr1_start < cdr2_start < cdr3_start):
        return None
    
    fr1 = full_seq[:cdr1_start]
    fr2 = full_seq[cdr1_start + len(cdr1):cdr2_start]
    fr3 = full_seq[cdr2_start + len(cdr2):cdr3_start]
    fr4 = full_seq[cdr3_start + len(cdr3):]
    
    if len(fr1) < 15 or len(fr2) < 10 or len(fr3) < 30 or len(fr4) < 8:
        return None
    if len(fr1) > 35 or len(fr2) > 20 or len(fr3) > 45 or len(fr4) > 15:
        return None
    
    return {
        'FR1': fr1, 'FR2': fr2, 'FR3': fr3, 'FR4': fr4,
        'cdr1': cdr1, 'cdr2': cdr2, 'cdr3': cdr3
    }

def collect_training_data(npz_files, max_samples=None, verbose=True):
    """Collect training data from NPZ files, sampling evenly across shards"""
    
    # Target positions for modeling (most important FR positions)
    # Using relative positions from start/end
    target_positions = {
        'FR1': list(range(-10, 0)),  # Last 10 positions
        'FR2': list(range(15)),       # First 15 positions (full FR2)
        'FR3': list(range(10)) + list(range(-10, 0)),  # First 10 + Last 10
        'FR4': list(range(11))        # First 11 positions (full FR4)
    }
    
    # Data storage: {(fw_name, fw_pos): {'X': [], 'y': []}}
    position_data = defaultdict(lambda: {'X': [], 'y': []})
    
    n_collected = 0
    feature_names = None
    
    # Calculate samples per shard for even distribution
    if max_samples and len(npz_files) > 0:
        samples_per_shard = max_samples // len(npz_files)
        if verbose:
            print(f"  Sampling ~{samples_per_shard:,} per shard ({max_samples:,} total across {len(npz_files)} shards)")
    else:
        samples_per_shard = None
        if verbose:
            print(f"  Using ALL samples from {len(npz_files)} shards")
    
    for npz_path in npz_files:
        if verbose:
            print(f"  Loading: {os.path.basename(npz_path)}", end=" ", flush=True)
        
        try:
            data = np.load(npz_path, allow_pickle=True)
            full_seqs = data['aa_v_full']
            cdr1s = data['cdr1']
            cdr2s = data['cdr2']
            cdr3s = data['cdr3']
        except Exception as e:
            print(f"ERROR: {e}")
            continue
        
        n_seqs = len(full_seqs)
        
        # Determine indices to sample from this shard
        if samples_per_shard and samples_per_shard < n_seqs:
            indices = np.random.choice(n_seqs, size=samples_per_shard, replace=False)
        else:
            indices = range(n_seqs)
        
        n_valid = 0
        for i in indices:
            fw = extract_frameworks(full_seqs[i], cdr1s[i], cdr2s[i], cdr3s[i])
            if fw is None:
                continue
            
            # Extract features
            features, feature_names = extract_cdr_features(
                fw['cdr1'], fw['cdr2'], fw['cdr3']
            )
            
            # Collect targets for each FR position
            for fw_name, positions in target_positions.items():
                fw_seq = fw[fw_name]
                for pos in positions:
                    try:
                        aa = fw_seq[pos]
                        if aa not in AA_TO_IDX:
                            continue
                        
                        key = (fw_name, pos)
                        position_data[key]['X'].append(features)
                        position_data[key]['y'].append(AA_TO_IDX[aa])
                    except IndexError:
                        continue
            
            n_valid += 1
            n_collected += 1
        
        if verbose:
            print(f"→ {n_valid:,} sampled")
    
    # Convert to numpy arrays
    for key in position_data:
        position_data[key]['X'] = np.array(position_data[key]['X'])
        position_data[key]['y'] = np.array(position_data[key]['y'])
    
    return dict(position_data), feature_names, n_collected

# ============================================================
# MODEL TRAINING
# ============================================================

def train_models(position_data, feature_names, verbose=True):
    """Train logistic regression model for each FR position"""
    
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler
    from sklearn.model_selection import cross_val_score
    
    models = {}
    scalers = {}
    metrics = {}
    
    positions = sorted(position_data.keys())
    
    if verbose:
        print(f"\nTraining models for {len(positions)} FR positions...")
    
    for i, (fw_name, fw_pos) in enumerate(positions):
        X = position_data[(fw_name, fw_pos)]['X']
        y = position_data[(fw_name, fw_pos)]['y']
        
        if len(X) < 100:
            continue
        
        # Check class distribution
        unique_classes = np.unique(y)
        if len(unique_classes) < 2:
            # Only one amino acid at this position - store as constant
            models[(fw_name, fw_pos)] = {
                'type': 'constant',
                'value': unique_classes[0],
                'aa': AA_LIST[unique_classes[0]]
            }
            metrics[(fw_name, fw_pos)] = {
                'accuracy': 1.0,
                'n_samples': len(X),
                'n_classes': 1,
                'top_class': AA_LIST[unique_classes[0]],
                'top_class_freq': 1.0
            }
            continue
        
        # Scale features
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X)
        scalers[(fw_name, fw_pos)] = scaler
        
        # Train logistic regression
        # Use class_weight='balanced' to handle imbalanced classes
        model = LogisticRegression(
            max_iter=500,
            multi_class='multinomial',
            solver='lbfgs',
            class_weight='balanced',
            n_jobs=-1
        )
        
        try:
            model.fit(X_scaled, y)
            
            # Compute accuracy
            train_acc = model.score(X_scaled, y)
            
            # Get class distribution
            class_counts = Counter(y)
            top_class = max(class_counts, key=class_counts.get)
            top_class_freq = class_counts[top_class] / len(y)
            
            models[(fw_name, fw_pos)] = {
                'type': 'logistic',
                'model': model,
                'classes': model.classes_
            }
            
            metrics[(fw_name, fw_pos)] = {
                'accuracy': train_acc,
                'n_samples': len(X),
                'n_classes': len(unique_classes),
                'top_class': AA_LIST[top_class],
                'top_class_freq': top_class_freq,
                'improvement_over_baseline': train_acc - top_class_freq
            }
            
        except Exception as e:
            if verbose:
                print(f"  Warning: Failed to train {fw_name}[{fw_pos}]: {e}")
            continue
        
        if verbose and (i + 1) % 20 == 0:
            print(f"  Trained {i+1}/{len(positions)} models...")
    
    return models, scalers, metrics

# ============================================================
# FEATURE IMPORTANCE
# ============================================================

def compute_feature_importance(models, scalers, feature_names):
    """Extract feature importances from trained models"""
    
    importances = {}
    
    for (fw_name, fw_pos), model_info in models.items():
        if model_info['type'] != 'logistic':
            continue
        
        model = model_info['model']
        
        # Get coefficients (absolute mean across classes)
        coefs = np.abs(model.coef_).mean(axis=0)
        
        # Map to feature names
        importance_dict = {name: coefs[i] for i, name in enumerate(feature_names)}
        
        # Sort by importance
        sorted_importance = sorted(importance_dict.items(), key=lambda x: x[1], reverse=True)
        
        importances[(fw_name, fw_pos)] = sorted_importance
    
    return importances

# ============================================================
# SUMMARY STATISTICS
# ============================================================

def compute_summary_stats(metrics, importances, feature_names):
    """Compute summary statistics across all positions"""
    
    # Average accuracy by framework region
    fw_stats = defaultdict(list)
    for (fw_name, fw_pos), m in metrics.items():
        fw_stats[fw_name].append(m['accuracy'])
    
    fw_summary = {fw: np.mean(accs) for fw, accs in fw_stats.items()}
    
    # Most important features overall
    feature_importance_sum = defaultdict(float)
    feature_importance_count = defaultdict(int)
    
    for (fw_name, fw_pos), imp_list in importances.items():
        for feat_name, imp_value in imp_list:
            feature_importance_sum[feat_name] += imp_value
            feature_importance_count[feat_name] += 1
    
    avg_importance = {
        name: feature_importance_sum[name] / feature_importance_count[name]
        for name in feature_importance_sum
    }
    
    top_features = sorted(avg_importance.items(), key=lambda x: x[1], reverse=True)
    
    # Positions with highest model improvement (most CDR-dependent)
    improvement_ranking = sorted(
        [(k, v['improvement_over_baseline']) for k, v in metrics.items() 
         if 'improvement_over_baseline' in v],
        key=lambda x: x[1],
        reverse=True
    )
    
    return {
        'fw_avg_accuracy': fw_summary,
        'top_features': top_features[:15],
        'most_cdr_dependent_positions': improvement_ranking[:20],
        'total_positions_modeled': len(metrics)
    }

# ============================================================
# MAIN
# ============================================================

def main():
    parser = argparse.ArgumentParser(description='Build p(FR|CDR) models')
    parser.add_argument('npz_files', nargs='+', help='NPZ files to process')
    parser.add_argument('--output', '-o', default='.', help='Output directory')
    parser.add_argument('--max-samples', type=int, default=None, 
                        help='Max training samples (default: use all)')
    
    args = parser.parse_args()
    
    # Expand wildcards
    npz_files = []
    for pattern in args.npz_files:
        if '*' in pattern:
            npz_files.extend(glob(pattern))
        else:
            npz_files.append(pattern)
    
    npz_files = [f for f in npz_files if f.endswith('.npz') and os.path.exists(f)]
    
    if not npz_files:
        print("ERROR: No valid NPZ files found")
        sys.exit(1)
    
    print("="*70)
    print("CONDITIONAL p(FR | CDR) MODEL BUILDER")
    print("="*70)
    print(f"\nFound {len(npz_files)} NPZ files")
    print(f"Max training samples: {'ALL' if args.max_samples is None else f'{args.max_samples:,}'}\n")
    
    # Collect data
    print("Step 1: Collecting training data...")
    position_data, feature_names, n_samples = collect_training_data(
        npz_files, max_samples=args.max_samples
    )
    print(f"\n  Collected {n_samples:,} samples")
    print(f"  Positions with data: {len(position_data)}")
    print(f"  Features: {len(feature_names)}")
    
    # Train models
    print("\nStep 2: Training models...")
    models, scalers, metrics = train_models(position_data, feature_names)
    print(f"  Trained {len(models)} models")
    
    # Compute importances
    print("\nStep 3: Computing feature importances...")
    importances = compute_feature_importance(models, scalers, feature_names)
    
    # Summary stats
    summary = compute_summary_stats(metrics, importances, feature_names)
    
    # Save everything
    os.makedirs(args.output, exist_ok=True)
    
    model_path = os.path.join(args.output, 'pfr_cdr_models.pkl')
    with open(model_path, 'wb') as f:
        pickle.dump({
            'models': models,
            'scalers': scalers,
            'metrics': metrics,
            'feature_names': feature_names,
            'importances': importances,
            'summary': summary
        }, f)
    print(f"\n✓ Models saved: {model_path}")
    
    # Write summary
    summary_path = os.path.join(args.output, 'model_summary.txt')
    with open(summary_path, 'w') as f:
        f.write("="*70 + "\n")
        f.write("p(FR | CDR) MODEL SUMMARY\n")
        f.write("="*70 + "\n\n")
        
        f.write(f"Total samples: {n_samples:,}\n")
        f.write(f"Positions modeled: {summary['total_positions_modeled']}\n\n")
        
        f.write("AVERAGE ACCURACY BY FRAMEWORK REGION\n")
        f.write("-"*40 + "\n")
        for fw, acc in sorted(summary['fw_avg_accuracy'].items()):
            f.write(f"  {fw}: {acc:.1%}\n")
        
        f.write("\n\nTOP 15 PREDICTIVE FEATURES\n")
        f.write("-"*40 + "\n")
        for feat, imp in summary['top_features']:
            f.write(f"  {feat}: {imp:.4f}\n")
        
        f.write("\n\nMOST CDR-DEPENDENT FR POSITIONS\n")
        f.write("(Highest improvement over baseline)\n")
        f.write("-"*40 + "\n")
        for (fw, pos), improvement in summary['most_cdr_dependent_positions']:
            m = metrics[(fw, pos)]
            f.write(f"  {fw}[{pos}]: +{improvement:.1%} over baseline "
                   f"(acc={m['accuracy']:.1%}, baseline={m['top_class_freq']:.1%})\n")
    
    print(f"  Summary: {summary_path}")
    
    # Print key results
    print("\n" + "="*70)
    print("KEY RESULTS")
    print("="*70)
    
    print("\nAverage accuracy by framework:")
    for fw, acc in sorted(summary['fw_avg_accuracy'].items()):
        print(f"  {fw}: {acc:.1%}")
    
    print("\nTop predictive features:")
    for feat, imp in summary['top_features'][:10]:
        print(f"  {feat}: {imp:.4f}")
    
    print("\nMost CDR-dependent FR positions:")
    for (fw, pos), improvement in summary['most_cdr_dependent_positions'][:10]:
        print(f"  {fw}[{pos}]: +{improvement:.1%}")
    
    print(f"\n{'='*70}")
    print("NEXT: Use score_framework.py to evaluate your universal framework")
    print(f"{'='*70}")

if __name__ == '__main__':
    main()
