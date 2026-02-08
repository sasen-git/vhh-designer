#!/usr/bin/env python3
"""
Build models from existing shards WITHOUT recomputing correlations.
Much faster - just samples 100k sequences for model training.

Usage:
    pip install scikit-learn
    python build_models_only.py VHH_shards/*.npz --output ./results/2_models
"""

import argparse
import numpy as np
import pickle
import gc
import sys
import time
from pathlib import Path

# ============================================================
# CONFIGURATION
# ============================================================

AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY")
AA_TO_IDX = {aa: i for i, aa in enumerate(AMINO_ACIDS)}

CDR_POSITIONS = {
    'cdr1': list(range(6)),
    'cdr2': list(range(6)),
    'cdr3': [0, 1, 2, -3, -2, -1],
}

FR_POSITIONS = {
    'FR1': list(range(25)),
    'FR2': list(range(14)),
    'FR3': list(range(32)),
    'FR4': list(range(11)),
}

# ============================================================
# PROGRESS BAR
# ============================================================

class ProgressBar:
    def __init__(self, total, width=40, prefix=''):
        self.total = total
        self.width = width
        self.prefix = prefix
        self.current = 0
        self.start_time = time.time()
        
    def update(self, n=1):
        self.current += n
        self._display()
    
    def _display(self):
        if self.total == 0:
            return
        pct = self.current / self.total
        filled = int(self.width * pct)
        bar = '█' * filled + '░' * (self.width - filled)
        elapsed = time.time() - self.start_time
        if self.current > 0:
            eta = (elapsed / self.current) * (self.total - self.current)
            eta_str = f"ETA {eta:.0f}s"
        else:
            eta_str = "ETA --"
        sys.stdout.write(f'\r  {self.prefix}|{bar}| {self.current:,}/{self.total:,} ({pct*100:.1f}%) {eta_str}   ')
        sys.stdout.flush()
    
    def finish(self):
        elapsed = time.time() - self.start_time
        sys.stdout.write(f'\r  {self.prefix}|{"█" * self.width}| {self.total:,}/{self.total:,} (100%) Done in {elapsed:.1f}s\n')
        sys.stdout.flush()

# ============================================================
# HELPERS
# ============================================================

def get_position_value(seq, pos):
    if seq is None or len(seq) == 0:
        return None
    try:
        if pos < 0:
            if len(seq) >= abs(pos):
                return seq[pos]
        else:
            if len(seq) > pos:
                return seq[pos]
    except:
        pass
    return None

def extract_frameworks(full_seq, cdr1, cdr2, cdr3):
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
    
    if len(fr1) < 10 or len(fr2) < 10 or len(fr3) < 20 or len(fr4) < 5:
        return None
    if len(fr1) > 40 or len(fr2) > 25 or len(fr3) > 50 or len(fr4) > 20:
        return None
    
    return {'FR1': fr1, 'FR2': fr2, 'FR3': fr3, 'FR4': fr4,
            'cdr1': cdr1, 'cdr2': cdr2, 'cdr3': cdr3}

def extract_features(seq_data):
    features = {}
    for cdr in ['cdr1', 'cdr2', 'cdr3']:
        seq = seq_data.get(cdr)
        if not seq:
            return None
        features[f'{cdr}_len'] = len(seq)
        for pos in CDR_POSITIONS[cdr]:
            aa = get_position_value(seq, pos)
            features[f'{cdr}_{pos}'] = AA_TO_IDX.get(aa, -1)
    return features

def extract_targets(seq_data):
    targets = {}
    for fw in ['FR1', 'FR2', 'FR3', 'FR4']:
        seq = seq_data.get(fw)
        if not seq:
            return None
        for pos in FR_POSITIONS[fw]:
            aa = get_position_value(seq, pos)
            targets[f'{fw}_{pos}'] = AA_TO_IDX.get(aa, -1)
    return targets

# ============================================================
# MAIN
# ============================================================

def main():
    parser = argparse.ArgumentParser(description='Build models from shards (skip correlations)')
    parser.add_argument('npz_files', nargs='+', help='NPZ shard files')
    parser.add_argument('--output', '-o', required=True, help='Output directory')
    parser.add_argument('--max-samples', type=int, default=100000, help='Max samples for training')
    
    args = parser.parse_args()
    
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    npz_files = sorted(args.npz_files)
    max_samples = args.max_samples
    
    print()
    print("=" * 60)
    print("  Building Models (skipping correlation computation)")
    print("=" * 60)
    print(f"  Shards: {len(npz_files)}")
    print(f"  Max samples: {max_samples:,}")
    print("=" * 60)
    print()
    
    # Reservoir sampling across all shards
    samples = []
    n_seen = 0
    
    overall_start = time.time()
    
    for shard_idx, npz_path in enumerate(npz_files):
        shard_name = Path(npz_path).name
        print(f"  [{shard_idx+1}/{len(npz_files)}] {shard_name}")
        
        shard_start = time.time()
        
        try:
            data = np.load(npz_path, allow_pickle=True)
            
            if 'aa_v_full' not in data:
                print(f"    Skipping - no aa_v_full")
                continue
                
            full_seqs = data['aa_v_full']
            cdr1s = data['cdr1']
            cdr2s = data['cdr2']
            cdr3s = data['cdr3']
            
            n_total = len(full_seqs)
            n_sampled = 0
            
            for i in range(n_total):
                fw = extract_frameworks(full_seqs[i], cdr1s[i], cdr2s[i], cdr3s[i])
                if fw is None:
                    continue
                
                features = extract_features(fw)
                targets = extract_targets(fw)
                
                if features is None or targets is None:
                    continue
                
                n_seen += 1
                
                # Reservoir sampling
                if len(samples) < max_samples:
                    samples.append((features, targets))
                    n_sampled += 1
                else:
                    j = np.random.randint(0, n_seen)
                    if j < max_samples:
                        samples[j] = (features, targets)
                        n_sampled += 1
            
            data.close()
            
            shard_time = time.time() - shard_start
            print(f"    ✓ Sampled from {n_total:,} seqs | Time: {shard_time:.1f}s | Reservoir: {len(samples):,}")
            
        except Exception as e:
            print(f"    Error: {e}")
        
        gc.collect()
    
    print()
    print(f"  Total sequences seen: {n_seen:,}")
    print(f"  Samples collected: {len(samples):,}")
    print()
    
    # Build models
    print("-" * 60)
    print("  Training models...")
    print("-" * 60)
    
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler
    
    # Convert to arrays
    feature_names = sorted(samples[0][0].keys())
    target_names = sorted(samples[0][1].keys())
    
    X = np.array([[s[0][f] for f in feature_names] for s in samples])
    
    # Filter valid samples
    valid_mask = (X >= 0).all(axis=1)
    X = X[valid_mask]
    valid_samples = [s for s, v in zip(samples, valid_mask) if v]
    
    print(f"  Valid samples: {len(valid_samples):,}")
    
    # Scale
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    # Train models
    models = {}
    n_targets = len(target_names)
    
    pbar = ProgressBar(n_targets)
    
    for target_name in target_names:
        y = np.array([s[1][target_name] for s in valid_samples])
        
        valid_target = y >= 0
        if valid_target.sum() < 100:
            pbar.update(1)
            continue
        
        X_train = X_scaled[valid_target]
        y_train = y[valid_target]
        
        n_classes = len(np.unique(y_train))
        if n_classes < 2:
            pbar.update(1)
            continue
        
        try:
            model = LogisticRegression(max_iter=500, solver='lbfgs', n_jobs=1, random_state=42)
            model.fit(X_train, y_train)
            
            models[target_name] = {
                'model': model,
                'scaler': scaler,
                'feature_names': feature_names,
                'classes': model.classes_,
            }
        except:
            pass
        
        pbar.update(1)
    
    pbar.finish()
    
    # Save
    result = {
        'models': models,
        'feature_names': feature_names,
        'n_samples': len(valid_samples),
    }
    
    with open(output_dir / 'pfr_cdr_models.pkl', 'wb') as f:
        pickle.dump(result, f)
    
    total_time = time.time() - overall_start
    
    print()
    print("=" * 60)
    print("  COMPLETE")
    print("=" * 60)
    print(f"  Models built: {len(models)}")
    print(f"  Saved to: {output_dir / 'pfr_cdr_models.pkl'}")
    print(f"  Total time: {total_time:.1f}s")
    print("=" * 60)
    print()

if __name__ == '__main__':
    main()
