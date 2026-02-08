#!/usr/bin/env python3
"""
Visualize p(FR|CDR) Model Results
=================================
Creates visualizations from the trained models showing:
1. Feature importance by FR position
2. Most CDR-dependent positions
3. CDR property influences on FR

Usage:
    python visualize_pfr_models.py pfr_cdr_models.pkl --output ./figures
"""

import pickle
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import argparse
import os
from collections import defaultdict

plt.style.use('seaborn-v0_8-whitegrid')

AA_LIST = list('ACDEFGHIKLMNPQRSTVWY')

def load_models(pkl_path):
    with open(pkl_path, 'rb') as f:
        return pickle.load(f)

def create_feature_importance_heatmap(data, output_dir):
    """Create heatmap of feature importances across FR positions"""
    
    importances = data['importances']
    feature_names = data['feature_names']
    
    # Get all positions
    positions = sorted(importances.keys())
    
    # Create matrix
    matrix = np.zeros((len(feature_names), len(positions)))
    
    for j, pos in enumerate(positions):
        imp_dict = dict(importances[pos])
        for i, feat in enumerate(feature_names):
            matrix[i, j] = imp_dict.get(feat, 0)
    
    # Create figure
    fig, ax = plt.subplots(figsize=(20, 10))
    
    # Position labels
    pos_labels = [f"{p[0]}[{p[1]}]" for p in positions]
    
    # Plot heatmap
    im = ax.imshow(matrix, aspect='auto', cmap='YlOrRd')
    
    ax.set_xticks(range(len(positions)))
    ax.set_xticklabels(pos_labels, rotation=90, fontsize=8)
    ax.set_yticks(range(len(feature_names)))
    ax.set_yticklabels(feature_names, fontsize=9)
    
    ax.set_xlabel('Framework Position')
    ax.set_ylabel('CDR Feature')
    ax.set_title('Feature Importance by Framework Position', fontsize=14, fontweight='bold')
    
    plt.colorbar(im, ax=ax, label='Importance')
    
    plt.tight_layout()
    path = os.path.join(output_dir, 'feature_importance_heatmap.png')
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {path}")

def create_position_summary(data, output_dir):
    """Create summary of most important FR positions"""
    
    metrics = data['metrics']
    
    # Separate by framework
    fw_data = defaultdict(list)
    for (fw_name, fw_pos), m in metrics.items():
        fw_data[fw_name].append({
            'pos': fw_pos,
            'accuracy': m['accuracy'],
            'baseline': m['top_class_freq'],
            'improvement': m.get('improvement_over_baseline', 0),
            'top_class': m['top_class'],
            'n_classes': m.get('n_classes', 1)
        })
    
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    for ax, fw_name in zip(axes.flatten(), ['FR1', 'FR2', 'FR3', 'FR4']):
        fw = sorted(fw_data[fw_name], key=lambda x: x['pos'])
        
        positions = [d['pos'] for d in fw]
        accuracies = [d['accuracy'] for d in fw]
        baselines = [d['baseline'] for d in fw]
        n_classes = [d['n_classes'] for d in fw]
        
        x = np.arange(len(positions))
        width = 0.35
        
        bars1 = ax.bar(x - width/2, accuracies, width, label='Model Accuracy', color='steelblue')
        bars2 = ax.bar(x + width/2, baselines, width, label='Baseline (most common AA)', color='lightcoral', alpha=0.7)
        
        ax.set_xlabel('Position')
        ax.set_ylabel('Accuracy')
        ax.set_title(f'{fw_name} Position Analysis', fontsize=12, fontweight='bold')
        ax.set_xticks(x)
        ax.set_xticklabels(positions, fontsize=8)
        ax.legend()
        ax.set_ylim(0, 1.1)
        
        # Add diversity indicator
        for i, (bar, nc) in enumerate(zip(bars2, n_classes)):
            if nc > 1:
                ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02,
                       f'{nc}', ha='center', fontsize=7, color='darkred')
    
    plt.suptitle('Framework Position Conservation and Predictability\n(Number above bar = # of amino acids observed)',
                 fontsize=14, fontweight='bold')
    plt.tight_layout()
    
    path = os.path.join(output_dir, 'position_summary.png')
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {path}")

def create_top_features_plot(data, output_dir):
    """Create bar plot of top predictive features"""
    
    summary = data['summary']
    top_features = summary['top_features']
    
    fig, ax = plt.subplots(figsize=(12, 8))
    
    features = [f[0] for f in top_features]
    importances = [f[1] for f in top_features]
    
    # Color by feature type
    colors = []
    for f in features:
        if 'len' in f:
            colors.append('#e74c3c')
        elif 'charge' in f or 'hydro' in f or 'aromatic' in f or 'polar' in f:
            colors.append('#27ae60')
        else:
            colors.append('#3498db')
    
    ax.barh(range(len(features)), importances, color=colors, alpha=0.8)
    ax.set_yticks(range(len(features)))
    ax.set_yticklabels(features)
    ax.set_xlabel('Average Importance Across All FR Positions')
    ax.set_title('Top Predictive CDR Features\n(Red=Length, Green=Property, Blue=Position)', 
                 fontsize=14, fontweight='bold')
    ax.invert_yaxis()
    
    plt.tight_layout()
    path = os.path.join(output_dir, 'top_features.png')
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {path}")

def create_cdr_dependent_positions(data, output_dir):
    """Create visualization of CDR-dependent vs conserved positions"""
    
    metrics = data['metrics']
    
    # Calculate variability for each position
    positions = []
    for (fw_name, fw_pos), m in metrics.items():
        positions.append({
            'fw': fw_name,
            'pos': fw_pos,
            'label': f"{fw_name}[{fw_pos}]",
            'baseline': m['top_class_freq'],
            'n_classes': m.get('n_classes', 1),
            'top_aa': m['top_class']
        })
    
    # Sort by baseline (conserved positions first)
    positions = sorted(positions, key=lambda x: x['baseline'], reverse=True)
    
    fig, ax = plt.subplots(figsize=(16, 8))
    
    labels = [p['label'] for p in positions]
    baselines = [p['baseline'] for p in positions]
    n_classes = [p['n_classes'] for p in positions]
    
    # Color by conservation
    colors = ['#27ae60' if b > 0.95 else '#f1c40f' if b > 0.8 else '#e74c3c' for b in baselines]
    
    bars = ax.bar(range(len(labels)), baselines, color=colors, alpha=0.8)
    
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=90, fontsize=8)
    ax.set_ylabel('Conservation (Most Common AA Frequency)')
    ax.set_title('Framework Position Conservation\n(Green: >95%, Yellow: 80-95%, Red: <80%)',
                 fontsize=14, fontweight='bold')
    
    ax.axhline(y=0.95, color='green', linestyle='--', alpha=0.5)
    ax.axhline(y=0.80, color='orange', linestyle='--', alpha=0.5)
    
    plt.tight_layout()
    path = os.path.join(output_dir, 'position_conservation.png')
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {path}")

def create_all_visualizations(pkl_path, output_dir):
    """Create all visualizations"""
    
    os.makedirs(output_dir, exist_ok=True)
    
    print("Loading models...")
    data = load_models(pkl_path)
    
    print("\nGenerating visualizations...")
    create_feature_importance_heatmap(data, output_dir)
    create_position_summary(data, output_dir)
    create_top_features_plot(data, output_dir)
    create_cdr_dependent_positions(data, output_dir)
    
    print("\n✓ All visualizations complete!")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('pkl_file', help='Path to pfr_cdr_models.pkl')
    parser.add_argument('--output', '-o', default='./figures', help='Output directory')
    
    args = parser.parse_args()
    create_all_visualizations(args.pkl_file, args.output)

if __name__ == '__main__':
    main()
