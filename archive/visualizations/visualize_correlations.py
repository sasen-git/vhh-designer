#!/usr/bin/env python3
"""
Visualization Script for CDR-Framework Correlations
===================================================
Creates comprehensive visualizations from correlation_results.pkl

Usage:
    python visualize_correlations.py correlation_results.pkl
    
    # Or specify output directory:
    python visualize_correlations.py correlation_results.pkl --output ./figures
"""

import pickle
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
import seaborn as sns
from collections import defaultdict
import argparse
import os

# Set style
plt.style.use('seaborn-v0_8-whitegrid')

def load_results(pkl_path):
    """Load correlation results"""
    with open(pkl_path, 'rb') as f:
        data = pickle.load(f)
    return data

def create_single_aa_heatmaps(stats, output_dir):
    """Create heatmaps for single AA correlations"""
    
    associations = stats['single_aa_associations']
    
    # Group by CDR-FW pair
    cdr_fw_data = defaultdict(list)
    for assoc in associations:
        key = (assoc['cdr'], assoc['fw'])
        cdr_fw_data[key].append(assoc)
    
    # Create figure with subplots for each CDR-FW combination
    fig, axes = plt.subplots(3, 4, figsize=(24, 18))
    
    cdrs = ['cdr1', 'cdr2', 'cdr3']
    fws = ['FR1', 'FR2', 'FR3', 'FR4']
    
    for i, cdr in enumerate(cdrs):
        for j, fw in enumerate(fws):
            ax = axes[i, j]
            
            data = cdr_fw_data.get((cdr, fw), [])
            
            if not data:
                ax.text(0.5, 0.5, 'No data', ha='center', va='center')
                ax.set_title(f'{cdr.upper()} → {fw}')
                continue
            
            # Get unique positions
            cdr_positions = sorted(set(d['cdr_pos'] for d in data))
            fw_positions = sorted(set(d['fw_pos'] for d in data))
            
            # Create matrix of max confidence for each position pair
            matrix = np.zeros((len(cdr_positions), len(fw_positions)))
            
            for d in data:
                ci = cdr_positions.index(d['cdr_pos'])
                fi = fw_positions.index(d['fw_pos'])
                matrix[ci, fi] = max(matrix[ci, fi], d['confidence'])
            
            # Plot
            sns.heatmap(matrix, ax=ax, cmap='YlOrRd', vmin=50, vmax=100,
                       xticklabels=fw_positions, yticklabels=cdr_positions,
                       cbar_kws={'label': 'Max Confidence %'})
            
            ax.set_title(f'{cdr.upper()} → {fw}', fontsize=11, fontweight='bold')
            ax.set_xlabel(f'{fw} position')
            ax.set_ylabel(f'{cdr.upper()} position')
    
    plt.suptitle('Single Amino Acid Position Correlations\n(Max confidence per position pair)', 
                 fontsize=14, fontweight='bold')
    plt.tight_layout()
    
    path = os.path.join(output_dir, 'single_aa_heatmaps.png')
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {path}")

def create_triplet_summary(stats, output_dir):
    """Create triplet association summary visualization"""
    
    associations = stats['triplet_associations']
    
    # Group by junction type
    junction_data = defaultdict(list)
    for assoc in associations:
        key = (assoc['junction_from'], assoc['junction_to'])
        junction_data[key].append(assoc)
    
    fig, axes = plt.subplots(2, 3, figsize=(20, 12))
    axes = axes.flatten()
    
    junction_order = [
        ('CDR1_start', 'FR1_end'),
        ('CDR1_end', 'FR2_start'),
        ('CDR2_start', 'FR2_end'),
        ('CDR2_end', 'FR3_start'),
        ('CDR3_start', 'FR3_end'),
        ('CDR3_end', 'FR4_start'),
    ]
    
    colors = ['#e74c3c', '#e67e22', '#f1c40f', '#27ae60', '#3498db', '#9b59b6']
    
    for idx, (junction, ax) in enumerate(zip(junction_order, axes)):
        data = junction_data.get(junction, [])
        
        if not data:
            ax.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax.transAxes)
            ax.set_title(f'{junction[0]} → {junction[1]}')
            continue
        
        # Sort by confidence and take top 15
        data_sorted = sorted(data, key=lambda x: x['confidence'], reverse=True)[:15]
        
        labels = [f"{d['cdr_triplet']}→{d['predicted_fw_motif']}" for d in data_sorted]
        confidences = [d['confidence'] for d in data_sorted]
        ns = [d['n'] for d in data_sorted]
        
        bars = ax.barh(range(len(labels)), confidences, color=colors[idx], alpha=0.7)
        ax.set_yticks(range(len(labels)))
        ax.set_yticklabels(labels, fontsize=8)
        ax.set_xlabel('Confidence %')
        ax.set_xlim(0, 105)
        ax.axvline(x=90, color='gray', linestyle='--', alpha=0.5)
        
        # Add sample sizes
        for bar, n in zip(bars, ns):
            ax.text(bar.get_width() + 1, bar.get_y() + bar.get_height()/2,
                   f'n={n:,}', va='center', fontsize=7)
        
        ax.set_title(f'{junction[0]} → {junction[1]}', fontsize=11, fontweight='bold')
        ax.invert_yaxis()
    
    plt.suptitle('Triplet Motif Associations by Junction Type\n(Top 15 per junction)', 
                 fontsize=14, fontweight='bold')
    plt.tight_layout()
    
    path = os.path.join(output_dir, 'triplet_associations.png')
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {path}")

def create_length_correlations(stats, output_dir):
    """Create length correlation visualization"""
    
    correlations = stats['length_correlations']
    
    fig, axes = plt.subplots(3, 4, figsize=(20, 15))
    
    cdrs = ['cdr1', 'cdr2', 'cdr3']
    fws = ['FR1', 'FR2', 'FR3', 'FR4']
    
    for i, cdr in enumerate(cdrs):
        for j, fw in enumerate(fws):
            ax = axes[i, j]
            
            # Get data for this CDR-FW pair
            data = [c for c in correlations if c['cdr'] == cdr and c['fw'] == fw]
            
            if not data:
                ax.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax.transAxes)
                ax.set_title(f'{cdr.upper()} vs {fw}')
                continue
            
            # Create scatter plot
            cdr_lens = [d['cdr_len'] for d in data]
            fw_lens = [d['predicted_fw_len'] for d in data]
            confidences = [d['confidence'] for d in data]
            ns = [d['n'] for d in data]
            
            # Size by sample count
            sizes = [min(n/10, 500) for n in ns]
            
            scatter = ax.scatter(cdr_lens, fw_lens, c=confidences, s=sizes,
                               cmap='RdYlGn', vmin=50, vmax=100, alpha=0.7)
            
            ax.set_xlabel(f'{cdr.upper()} length')
            ax.set_ylabel(f'{fw} length')
            ax.set_title(f'{cdr.upper()} vs {fw}', fontsize=11, fontweight='bold')
            
            # Calculate correlation
            if len(cdr_lens) > 2:
                from scipy import stats as sp_stats
                r, p = sp_stats.spearmanr(cdr_lens, fw_lens)
                ax.text(0.05, 0.95, f'r={r:.2f}', transform=ax.transAxes,
                       fontsize=10, fontweight='bold', va='top')
    
    # Add colorbar
    fig.subplots_adjust(right=0.92)
    cbar_ax = fig.add_axes([0.94, 0.15, 0.02, 0.7])
    cbar = fig.colorbar(scatter, cax=cbar_ax)
    cbar.set_label('Confidence %')
    
    plt.suptitle('CDR-Framework Length Correlations\n(Size = sample count, Color = confidence)', 
                 fontsize=14, fontweight='bold')
    
    path = os.path.join(output_dir, 'length_correlations.png')
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {path}")

def create_top_rules_summary(stats, output_dir):
    """Create summary of top rules with biological interpretation"""
    
    fig = plt.figure(figsize=(20, 16))
    
    # Panel 1: Top single AA rules
    ax1 = fig.add_subplot(2, 2, 1)
    
    top_aa = sorted(stats['single_aa_associations'], 
                   key=lambda x: x['confidence'], reverse=True)[:20]
    
    labels = [f"{a['cdr']}[{a['cdr_pos']}]={a['cdr_aa']}→{a['fw']}[{a['fw_pos']}]={a['predicted_fw_aa']}" 
              for a in top_aa]
    confidences = [a['confidence'] for a in top_aa]
    
    colors = ['#e74c3c' if 'cdr1' in l else '#27ae60' if 'cdr2' in l else '#3498db' for l in labels]
    
    bars = ax1.barh(range(len(labels)), confidences, color=colors, alpha=0.7)
    ax1.set_yticks(range(len(labels)))
    ax1.set_yticklabels(labels, fontsize=9)
    ax1.set_xlabel('Confidence %')
    ax1.set_xlim(80, 102)
    ax1.set_title('Top 20 Single AA Rules', fontsize=12, fontweight='bold')
    ax1.invert_yaxis()
    
    # Panel 2: Top triplet rules
    ax2 = fig.add_subplot(2, 2, 2)
    
    top_trip = sorted(stats['triplet_associations'],
                     key=lambda x: x['confidence'], reverse=True)[:20]
    
    labels = [f"{t['cdr_triplet']}→{t['predicted_fw_motif']}" for t in top_trip]
    confidences = [t['confidence'] for t in top_trip]
    junctions = [t['junction_from'].replace('_', ' ') for t in top_trip]
    
    bars = ax2.barh(range(len(labels)), confidences, color='coral', alpha=0.7)
    ax2.set_yticks(range(len(labels)))
    ax2.set_yticklabels(labels, fontsize=9)
    
    # Add junction labels
    for bar, junc in zip(bars, junctions):
        ax2.text(bar.get_width() + 0.5, bar.get_y() + bar.get_height()/2,
                junc, va='center', fontsize=7, style='italic')
    
    ax2.set_xlabel('Confidence %')
    ax2.set_xlim(80, 110)
    ax2.set_title('Top 20 Triplet Rules', fontsize=12, fontweight='bold')
    ax2.invert_yaxis()
    
    # Panel 3: Distribution of confidences
    ax3 = fig.add_subplot(2, 2, 3)
    
    aa_conf = [a['confidence'] for a in stats['single_aa_associations']]
    trip_conf = [t['confidence'] for t in stats['triplet_associations']]
    
    ax3.hist(aa_conf, bins=20, alpha=0.6, label=f'Single AA (n={len(aa_conf):,})', color='steelblue')
    ax3.hist(trip_conf, bins=20, alpha=0.6, label=f'Triplet (n={len(trip_conf):,})', color='coral')
    ax3.axvline(x=90, color='red', linestyle='--', label='90% threshold')
    ax3.axvline(x=95, color='green', linestyle='--', label='95% threshold')
    ax3.set_xlabel('Confidence %')
    ax3.set_ylabel('Count')
    ax3.set_title('Distribution of Rule Confidence', fontsize=12, fontweight='bold')
    ax3.legend()
    
    # Panel 4: Summary statistics
    ax4 = fig.add_subplot(2, 2, 4)
    ax4.axis('off')
    
    n_seq = stats['summary']['n_sequences']
    n_aa = stats['summary']['n_single_aa_associations']
    n_trip = stats['summary']['n_triplet_associations']
    
    # Count high confidence rules
    n_aa_90 = len([a for a in stats['single_aa_associations'] if a['confidence'] >= 90])
    n_aa_95 = len([a for a in stats['single_aa_associations'] if a['confidence'] >= 95])
    n_trip_90 = len([t for t in stats['triplet_associations'] if t['confidence'] >= 90])
    n_trip_95 = len([t for t in stats['triplet_associations'] if t['confidence'] >= 95])
    
    summary_text = f"""
╔══════════════════════════════════════════════════════════════╗
║                    ANALYSIS SUMMARY                          ║
╠══════════════════════════════════════════════════════════════╣
║  Sequences Analyzed: {n_seq:>15,}                         ║
║                                                              ║
║  SINGLE AA ASSOCIATIONS                                      ║
║    Total rules: {n_aa:>20,}                         ║
║    ≥90% confidence: {n_aa_90:>16,}                         ║
║    ≥95% confidence: {n_aa_95:>16,}                         ║
║                                                              ║
║  TRIPLET ASSOCIATIONS                                        ║
║    Total rules: {n_trip:>20,}                         ║
║    ≥90% confidence: {n_trip_90:>16,}                         ║
║    ≥95% confidence: {n_trip_95:>16,}                         ║
╚══════════════════════════════════════════════════════════════╝

KEY BIOLOGICAL RULES DISCOVERED:
"""
    
    # Add top biological rules
    for a in top_aa[:5]:
        summary_text += f"\n• {a['cdr']}[{a['cdr_pos']}]={a['cdr_aa']} → {a['fw']}[{a['fw_pos']}]={a['predicted_fw_aa']} ({a['confidence']:.1f}%)"
    
    ax4.text(0.05, 0.95, summary_text, transform=ax4.transAxes,
            fontsize=10, verticalalignment='top', fontfamily='monospace',
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.9))
    
    plt.suptitle('CDR-Framework Correlation Analysis Summary', fontsize=14, fontweight='bold')
    plt.tight_layout()
    
    path = os.path.join(output_dir, 'analysis_summary.png')
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {path}")

def create_visualizations(stats, output_dir):
    """Create all visualizations"""
    
    os.makedirs(output_dir, exist_ok=True)
    
    print("\nGenerating visualizations...")
    create_single_aa_heatmaps(stats, output_dir)
    create_triplet_summary(stats, output_dir)
    create_length_correlations(stats, output_dir)
    create_top_rules_summary(stats, output_dir)
    print("\n✓ All visualizations complete!")

def main():
    parser = argparse.ArgumentParser(description='Visualize CDR-Framework correlations')
    parser.add_argument('pkl_file', help='Path to correlation_results.pkl')
    parser.add_argument('--output', '-o', default='./figures', help='Output directory')
    
    args = parser.parse_args()
    
    print("Loading results...")
    data = load_results(args.pkl_file)
    stats = data['statistics']
    
    print(f"Loaded: {stats['summary']['n_sequences']:,} sequences")
    print(f"  Single AA associations: {stats['summary']['n_single_aa_associations']:,}")
    print(f"  Triplet associations: {stats['summary']['n_triplet_associations']:,}")
    
    create_visualizations(stats, args.output)

if __name__ == '__main__':
    main()
