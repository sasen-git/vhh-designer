#!/usr/bin/env python3
"""
LOW MEMORY VERSION - CDR-Framework Analysis Pipeline

Processes shards ONE AT A TIME, accumulating statistics without
holding all sequences in memory. Suitable for large datasets on
memory-constrained systems (e.g., WSL with limited RAM).

Usage:
    python run_full_analysis_lowmem.py /path/to/shards/*.npz --output ./results
    python run_full_analysis_lowmem.py /path/to/shards/*.npz --output ./results --max-shards 10
"""

import argparse
import numpy as np
import pickle
import gc
import os
import sys
import time
from pathlib import Path
from collections import defaultdict
from datetime import datetime

# ============================================================
# CONFIGURATION
# ============================================================

AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY")
AA_TO_IDX = {aa: i for i, aa in enumerate(AMINO_ACIDS)}
N_AA = len(AMINO_ACIDS)

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

CDR_FR_PAIRS = [
    ('cdr1', 'FR1'), ('cdr1', 'FR2'),
    ('cdr2', 'FR2'), ('cdr2', 'FR3'),
    ('cdr3', 'FR3'), ('cdr3', 'FR4'),
]

VERNIER_POSITIONS = {
    'FR1': [24],
    'FR2': [1, 11, 12],
    'FR3': [2, 23, 24, 29, 30],
    'FR4': [1, 3],
}

# ============================================================
# PROGRESS BAR
# ============================================================

class ProgressBar:
    """Simple progress bar that updates in place."""
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
        
        sys.stdout.write(f'\r    {self.prefix}|{bar}| {self.current:,}/{self.total:,} ({pct*100:.1f}%) {eta_str}   ')
        sys.stdout.flush()
    
    def finish(self):
        elapsed = time.time() - self.start_time
        sys.stdout.write(f'\r    {self.prefix}|{"█" * self.width}| {self.total:,}/{self.total:,} (100%) Done in {elapsed:.1f}s\n')
        sys.stdout.flush()

# ============================================================
# HELPER FUNCTIONS
# ============================================================

def get_cdr_position_value(seq, pos):
    """Get AA at position (handles negative indexing safely)."""
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

def get_fw_position_value(seq, pos):
    """Get AA at framework position."""
    if seq is None or len(seq) <= pos:
        return None
    return seq[pos]

def extract_frameworks(full_seq, cdr1, cdr2, cdr3):
    """Extract framework regions by finding CDRs in full sequence."""
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
    
    # Sanity checks on framework lengths
    if len(fr1) < 10 or len(fr2) < 10 or len(fr3) < 20 or len(fr4) < 5:
        return None
    if len(fr1) > 40 or len(fr2) > 25 or len(fr3) > 50 or len(fr4) > 20:
        return None
    
    return {'FR1': fr1, 'FR2': fr2, 'FR3': fr3, 'FR4': fr4,
            'cdr1': cdr1, 'cdr2': cdr2, 'cdr3': cdr3}

def process_shard_streaming(npz_path, corr_accumulator, model_accumulator=None, show_progress=True):
    """Process a single shard, streaming sequences and updating accumulators."""
    try:
        data = np.load(npz_path, allow_pickle=True)
        
        # Check for the expected format with separate arrays
        if 'aa_v_full' in data and 'cdr1' in data and 'cdr2' in data and 'cdr3' in data:
            full_seqs = data['aa_v_full']
            cdr1s = data['cdr1']
            cdr2s = data['cdr2']
            cdr3s = data['cdr3']
            
            n_total = len(full_seqs)
            n_valid = 0
            
            if show_progress:
                pbar = ProgressBar(n_total, prefix='')
            
            for i in range(n_total):
                fw = extract_frameworks(full_seqs[i], cdr1s[i], cdr2s[i], cdr3s[i])
                if fw:
                    corr_accumulator.add_sequence(fw)
                    if model_accumulator:
                        model_accumulator.add_sequence(fw)
                    n_valid += 1
                
                if show_progress and i % 10000 == 0:
                    pbar.update(10000 if i > 0 else 0)
            
            if show_progress:
                pbar.current = n_total
                pbar.finish()
            
            data.close()
            return n_valid, n_total
        
        # Try alternative format with annotations dict
        n_valid = 0
        n_total = 0
        
        if 'annotations' in data:
            annotations = data['annotations']
        elif 'arr_0' in data:
            annotations = data['arr_0']
        else:
            keys = list(data.keys())
            if keys:
                annotations = data[keys[0]]
            else:
                data.close()
                return 0, 0
        
        if hasattr(annotations, 'item'):
            annotations = annotations.item()
        
        if isinstance(annotations, dict):
            annotations = [annotations]
        
        n_total = len(annotations)
        if show_progress:
            pbar = ProgressBar(n_total, prefix='')
        
        for i, ann in enumerate(annotations):
            if not isinstance(ann, dict):
                continue
            
            regions = ann.get('regions', {})
            if not regions:
                continue
            
            seq_data = {
                'cdr1': regions.get('CDR1', {}).get('sequence'),
                'cdr2': regions.get('CDR2', {}).get('sequence'),
                'cdr3': regions.get('CDR3', {}).get('sequence'),
                'FR1': regions.get('FR1', {}).get('sequence'),
                'FR2': regions.get('FR2', {}).get('sequence'),
                'FR3': regions.get('FR3', {}).get('sequence'),
                'FR4': regions.get('FR4', {}).get('sequence'),
            }
            
            if seq_data['cdr3'] and seq_data['FR4']:
                corr_accumulator.add_sequence(seq_data)
                if model_accumulator:
                    model_accumulator.add_sequence(seq_data)
                n_valid += 1
            
            if show_progress and i % 10000 == 0:
                pbar.update(10000 if i > 0 else 0)
        
        if show_progress:
            pbar.current = n_total
            pbar.finish()
                
        data.close()
        return n_valid, n_total
        
    except Exception as e:
        print(f"\n  Warning: Error processing {npz_path}: {e}")
        return 0, 0

# ============================================================
# STREAMING STATISTICS ACCUMULATORS
# ============================================================

class StreamingCorrelationAccumulator:
    """Accumulates correlation statistics without storing all sequences."""
    
    def __init__(self):
        # Single position counts: counts[cdr][cdr_pos][fw][fw_pos][cdr_aa][fw_aa]
        self.joint_counts = defaultdict(lambda: defaultdict(lambda: defaultdict(
            lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(int))))))
        
        # Marginal counts
        self.cdr_counts = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
        self.fw_counts = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
        
        # Triplet counts for junctions
        self.triplet_counts = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
        
        self.n_sequences = 0
    
    def add_sequence(self, seq_data):
        """Add one sequence's statistics."""
        self.n_sequences += 1
        
        # Process each CDR-FW pair
        for cdr, fw in CDR_FR_PAIRS:
            cdr_seq = seq_data.get(cdr)
            fw_seq = seq_data.get(fw)
            
            if not cdr_seq or not fw_seq:
                continue
            
            # Single position correlations
            for cdr_pos in CDR_POSITIONS[cdr]:
                cdr_aa = get_cdr_position_value(cdr_seq, cdr_pos)
                if cdr_aa not in AA_TO_IDX:
                    continue
                
                self.cdr_counts[cdr][cdr_pos][cdr_aa] += 1
                
                for fw_pos in FR_POSITIONS[fw]:
                    fw_aa = get_fw_position_value(fw_seq, fw_pos)
                    if fw_aa not in AA_TO_IDX:
                        continue
                    
                    self.joint_counts[cdr][cdr_pos][fw][fw_pos][cdr_aa][fw_aa] += 1
            
            # Framework marginals
            for fw_pos in FR_POSITIONS[fw]:
                fw_aa = get_fw_position_value(fw_seq, fw_pos)
                if fw_aa in AA_TO_IDX:
                    self.fw_counts[fw][fw_pos][fw_aa] += 1
        
        # CDR3-FR4 triplet junction
        cdr3 = seq_data.get('cdr3')
        fr4 = seq_data.get('FR4')
        if cdr3 and fr4 and len(cdr3) >= 3 and len(fr4) >= 3:
            triplet = cdr3[-3:]
            fr4_start = fr4[:3]
            if all(aa in AA_TO_IDX for aa in triplet) and all(aa in AA_TO_IDX for aa in fr4_start):
                self.triplet_counts['cdr3_fr4'][triplet][fr4_start] += 1
    
    def compute_results(self, min_count=50, min_confidence=90.0):
        """Compute final correlation results from accumulated counts."""
        results = {
            'n_sequences': self.n_sequences,
            'universal_rules': [],
            'triplet_rules': [],
            'position_stats': {},
        }
        
        # Compute single-position rules
        for cdr in self.joint_counts:
            for cdr_pos in self.joint_counts[cdr]:
                for fw in self.joint_counts[cdr][cdr_pos]:
                    for fw_pos in self.joint_counts[cdr][cdr_pos][fw]:
                        joint = self.joint_counts[cdr][cdr_pos][fw][fw_pos]
                        
                        for cdr_aa in joint:
                            total = sum(joint[cdr_aa].values())
                            if total < min_count:
                                continue
                            
                            # Find most common FW residue
                            best_fw_aa = max(joint[cdr_aa], key=joint[cdr_aa].get)
                            best_count = joint[cdr_aa][best_fw_aa]
                            confidence = 100.0 * best_count / total
                            
                            if confidence >= min_confidence:
                                is_vernier = fw_pos in VERNIER_POSITIONS.get(fw, [])
                                results['universal_rules'].append({
                                    'cdr': cdr,
                                    'cdr_pos': cdr_pos,
                                    'cdr_aa': cdr_aa,
                                    'fw': fw,
                                    'fw_pos': fw_pos,
                                    'predicted_fw_aa': best_fw_aa,
                                    'confidence': confidence,
                                    'n': total,
                                    'is_vernier': is_vernier,
                                })
        
        # Compute triplet rules
        for junction in self.triplet_counts:
            for triplet in self.triplet_counts[junction]:
                total = sum(self.triplet_counts[junction][triplet].values())
                if total < min_count:
                    continue
                
                best_fr = max(self.triplet_counts[junction][triplet], 
                             key=self.triplet_counts[junction][triplet].get)
                best_count = self.triplet_counts[junction][triplet][best_fr]
                confidence = 100.0 * best_count / total
                
                if confidence >= min_confidence:
                    results['triplet_rules'].append({
                        'junction': junction,
                        'cdr_triplet': triplet,
                        'fw_triplet': best_fr,
                        'confidence': confidence,
                        'n': total,
                    })
        
        # Sort rules by confidence
        results['universal_rules'].sort(key=lambda x: x['confidence'], reverse=True)
        results['triplet_rules'].sort(key=lambda x: x['confidence'], reverse=True)
        
        return results

# ============================================================
# MODEL BUILDING (Streaming version)
# ============================================================

class StreamingModelDataAccumulator:
    """Accumulates data for model building without storing all sequences."""
    
    def __init__(self, max_samples=100000):
        self.max_samples = max_samples
        self.samples = []  # Will store (features, targets) tuples
        self.n_seen = 0
        
    def add_sequence(self, seq_data):
        """Reservoir sampling to keep max_samples representative sequences."""
        self.n_seen += 1
        
        # Extract features
        features = self._extract_features(seq_data)
        if features is None:
            return
            
        targets = self._extract_targets(seq_data)
        if targets is None:
            return
        
        # Reservoir sampling
        if len(self.samples) < self.max_samples:
            self.samples.append((features, targets))
        else:
            # Randomly replace with decreasing probability
            j = np.random.randint(0, self.n_seen)
            if j < self.max_samples:
                self.samples[j] = (features, targets)
    
    def _extract_features(self, seq_data):
        """Extract CDR features for model."""
        features = {}
        
        for cdr in ['cdr1', 'cdr2', 'cdr3']:
            seq = seq_data.get(cdr)
            if not seq:
                return None
            
            features[f'{cdr}_len'] = len(seq)
            
            # Terminal residues
            for pos in CDR_POSITIONS[cdr]:
                aa = get_cdr_position_value(seq, pos)
                features[f'{cdr}_{pos}'] = AA_TO_IDX.get(aa, -1)
        
        return features
    
    def _extract_targets(self, seq_data):
        """Extract framework targets for model."""
        targets = {}
        
        for fw in ['FR1', 'FR2', 'FR3', 'FR4']:
            seq = seq_data.get(fw)
            if not seq:
                return None
            
            for pos in FR_POSITIONS[fw]:
                aa = get_fw_position_value(seq, pos)
                targets[f'{fw}_{pos}'] = AA_TO_IDX.get(aa, -1)
        
        return targets
    
    def build_models(self):
        """Build logistic regression models from accumulated samples."""
        from sklearn.linear_model import LogisticRegression
        from sklearn.preprocessing import StandardScaler
        
        if len(self.samples) < 100:
            print(f"  Warning: Only {len(self.samples)} samples, skipping model building")
            return None
        
        print(f"  Building models from {len(self.samples):,} samples...")
        
        # Convert to arrays
        feature_names = sorted(self.samples[0][0].keys())
        target_names = sorted(self.samples[0][1].keys())
        
        X = np.array([[s[0][f] for f in feature_names] for s in self.samples])
        
        # Remove samples with invalid features
        valid_mask = (X >= 0).all(axis=1)
        X = X[valid_mask]
        
        valid_samples = [s for s, v in zip(self.samples, valid_mask) if v]
        
        print(f"  {len(valid_samples):,} valid samples after filtering")
        
        # Scale features
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X)
        
        # Build model for each framework position
        models = {}
        n_targets = len(target_names)
        
        print(f"  Training {n_targets} position models...")
        pbar = ProgressBar(n_targets, prefix='')
        
        for i, target_name in enumerate(target_names):
            y = np.array([s[1][target_name] for s in valid_samples])
            
            # Skip if target has invalid values
            valid_target = y >= 0
            if valid_target.sum() < 100:
                pbar.update(1)
                continue
            
            X_train = X_scaled[valid_target]
            y_train = y[valid_target]
            
            # Skip if not enough classes
            n_classes = len(np.unique(y_train))
            if n_classes < 2:
                pbar.update(1)
                continue
            
            try:
                model = LogisticRegression(
                    max_iter=500,
                    solver='lbfgs',
                    n_jobs=1,  # Reduce memory
                    random_state=42
                )
                model.fit(X_train, y_train)
                
                models[target_name] = {
                    'model': model,
                    'scaler': scaler,
                    'feature_names': feature_names,
                    'classes': model.classes_,
                }
            except Exception as e:
                pass  # Skip failed models silently
            
            pbar.update(1)
        
        pbar.finish()
        
        return {
            'models': models,
            'feature_names': feature_names,
            'n_samples': len(valid_samples),
        }

# ============================================================
# MAIN PIPELINE
# ============================================================

def format_time(seconds):
    """Format seconds into human readable string."""
    if seconds < 60:
        return f"{seconds:.1f}s"
    elif seconds < 3600:
        return f"{seconds/60:.1f}m"
    else:
        return f"{seconds/3600:.1f}h"

def main():
    parser = argparse.ArgumentParser(description='Low-memory CDR-Framework Analysis')
    parser.add_argument('npz_files', nargs='+', help='NPZ shard files')
    parser.add_argument('--output', '-o', required=True, help='Output directory')
    parser.add_argument('--max-shards', type=int, default=None, help='Max shards to process')
    parser.add_argument('--max-samples-model', type=int, default=100000, 
                        help='Max samples for model building (default: 100000)')
    parser.add_argument('--skip-models', action='store_true', help='Skip model building')
    parser.add_argument('--min-count', type=int, default=50, help='Min count for rules')
    parser.add_argument('--min-confidence', type=float, default=90.0, help='Min confidence %')
    
    args = parser.parse_args()
    
    # Setup output directory
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    npz_files = sorted(args.npz_files)
    if args.max_shards:
        npz_files = npz_files[:args.max_shards]
    
    print()
    print(f"{'='*70}")
    print(f"  LOW-MEMORY CDR-Framework Analysis Pipeline")
    print(f"{'='*70}")
    print(f"  Shards to process: {len(npz_files)}")
    print(f"  Output directory:  {output_dir}")
    print(f"  Model samples:     {args.max_samples_model:,}")
    print(f"{'='*70}")
    print()
    
    overall_start = time.time()
    shard_times = []
    
    # Initialize accumulators
    corr_accumulator = StreamingCorrelationAccumulator()
    model_accumulator = None if args.skip_models else StreamingModelDataAccumulator(args.max_samples_model)
    
    # Process shards one at a time
    print("STEP 1: Processing shards (streaming mode)")
    print("-" * 70)
    
    total_valid = 0
    total_seqs = 0
    
    for i, npz_file in enumerate(npz_files):
        shard_start = time.time()
        shard_name = Path(npz_file).name
        
        print(f"\n  [{i+1}/{len(npz_files)}] {shard_name}")
        
        n_valid, n_total = process_shard_streaming(
            npz_file, corr_accumulator, model_accumulator, show_progress=True
        )
        
        shard_time = time.time() - shard_start
        shard_times.append(shard_time)
        
        total_valid += n_valid
        total_seqs += n_total
        
        # Calculate ETA
        avg_time = sum(shard_times) / len(shard_times)
        remaining = len(npz_files) - (i + 1)
        eta = avg_time * remaining
        
        print(f"    ✓ {n_valid:,}/{n_total:,} valid | Time: {format_time(shard_time)} | ETA: {format_time(eta)}")
        
        # Force garbage collection after each shard
        gc.collect()
    
    step1_time = time.time() - overall_start
    
    print(f"\n  Total: {total_valid:,} valid sequences from {total_seqs:,}")
    print(f"  Step 1 time: {format_time(step1_time)}")
    
    # Compute correlation results
    print()
    print("-" * 70)
    print("STEP 2: Computing correlation statistics")
    print("-" * 70)
    
    step2_start = time.time()
    corr_results = corr_accumulator.compute_results(
        min_count=args.min_count, 
        min_confidence=args.min_confidence
    )
    step2_time = time.time() - step2_start
    
    # Save correlation results
    corr_dir = output_dir / "1_correlations"
    corr_dir.mkdir(exist_ok=True)
    
    with open(corr_dir / "correlation_results.pkl", 'wb') as f:
        pickle.dump(corr_results, f)
    
    # Write summary
    with open(corr_dir / "correlation_summary.txt", 'w') as f:
        f.write(f"CDR-Framework Correlation Analysis (Low-Memory)\n")
        f.write(f"=" * 50 + "\n\n")
        f.write(f"Total sequences: {total_valid:,}\n")
        f.write(f"Universal rules found: {len(corr_results['universal_rules'])}\n")
        f.write(f"Triplet rules found: {len(corr_results['triplet_rules'])}\n\n")
        
        f.write(f"Top Universal Rules:\n")
        f.write(f"-" * 40 + "\n")
        for rule in corr_results['universal_rules'][:20]:
            v = "[V]" if rule['is_vernier'] else ""
            f.write(f"  {rule['cdr']}[{rule['cdr_pos']}]={rule['cdr_aa']} → "
                   f"{rule['fw']}[{rule['fw_pos']}]={rule['predicted_fw_aa']} "
                   f"({rule['confidence']:.1f}%, n={rule['n']}) {v}\n")
        
        f.write(f"\nTop Triplet Rules:\n")
        f.write(f"-" * 40 + "\n")
        for rule in corr_results['triplet_rules'][:10]:
            f.write(f"  {rule['cdr_triplet']} → {rule['fw_triplet']} "
                   f"({rule['confidence']:.1f}%, n={rule['n']})\n")
    
    print(f"  ✓ Found {len(corr_results['universal_rules']):,} universal rules")
    print(f"  ✓ Found {len(corr_results['triplet_rules']):,} triplet rules")
    print(f"  ✓ Saved to {corr_dir}")
    print(f"  Step 2 time: {format_time(step2_time)}")
    
    # Free correlation accumulator memory
    del corr_accumulator
    gc.collect()
    
    # Build models if requested
    if model_accumulator and not args.skip_models:
        print()
        print("-" * 70)
        print("STEP 3: Building predictive models")
        print("-" * 70)
        
        step3_start = time.time()
        model_dir = output_dir / "2_models"
        model_dir.mkdir(exist_ok=True)
        
        model_results = model_accumulator.build_models()
        step3_time = time.time() - step3_start
        
        if model_results:
            with open(model_dir / "pfr_cdr_models.pkl", 'wb') as f:
                pickle.dump(model_results, f)
            print(f"  ✓ Built {len(model_results['models'])} position models")
            print(f"  ✓ Saved to {model_dir}")
        
        print(f"  Step 3 time: {format_time(step3_time)}")
        
        del model_accumulator
        gc.collect()
    
    # Final summary
    total_time = time.time() - overall_start
    
    print()
    print(f"{'='*70}")
    print(f"  COMPLETE")
    print(f"{'='*70}")
    print(f"  Total sequences:  {total_valid:,}")
    print(f"  Universal rules:  {len(corr_results['universal_rules']):,}")
    print(f"  Triplet rules:    {len(corr_results['triplet_rules']):,}")
    print(f"  Total time:       {format_time(total_time)}")
    print(f"  Output:           {output_dir}")
    print(f"{'='*70}")
    print()
    
    # Write final report
    with open(output_dir / "ANALYSIS_REPORT.md", 'w') as f:
        f.write(f"# CDR-Framework Analysis Report\n\n")
        f.write(f"**Generated:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        f.write(f"## Summary\n\n")
        f.write(f"- **Sequences analyzed:** {total_valid:,}\n")
        f.write(f"- **Shards processed:** {len(npz_files)}\n")
        f.write(f"- **Processing time:** {format_time(total_time)}\n")
        f.write(f"- **Universal rules:** {len(corr_results['universal_rules']):,}\n")
        f.write(f"- **Triplet rules:** {len(corr_results['triplet_rules']):,}\n\n")
        
        f.write(f"## Top Rules\n\n")
        f.write(f"### Single Position (Top 10)\n\n")
        f.write(f"| CDR Position | AA | → | FR Position | AA | Confidence | N |\n")
        f.write(f"|--------------|----|----|-------------|----|-----------:|--:|\n")
        for rule in corr_results['universal_rules'][:10]:
            f.write(f"| {rule['cdr']}[{rule['cdr_pos']}] | {rule['cdr_aa']} | → | "
                   f"{rule['fw']}[{rule['fw_pos']}] | {rule['predicted_fw_aa']} | "
                   f"{rule['confidence']:.1f}% | {rule['n']:,} |\n")
        
        f.write(f"\n### Junction Triplets (Top 10)\n\n")
        f.write(f"| CDR3 End | → | FR4 Start | Confidence | N |\n")
        f.write(f"|----------|---|-----------|----------:|--:|\n")
        for rule in corr_results['triplet_rules'][:10]:
            f.write(f"| {rule['cdr_triplet']} | → | {rule['fw_triplet']} | "
                   f"{rule['confidence']:.1f}% | {rule['n']:,} |\n")

if __name__ == '__main__':
    main()