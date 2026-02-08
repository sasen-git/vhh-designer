#!/usr/bin/env python3
"""
Visualize VHH Epistasis Analysis Results
"""

import pickle
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from collections import defaultdict

# Load results
with open('epistasis_overnight_full.pkl', 'rb') as f:
    results = pickle.load(f)

# Create figure directory
import os
os.makedirs('epistasis_figures', exist_ok=True)

# =============================================================================
# 1. MI HEATMAP - Which positions co-evolve?
# =============================================================================
print("Creating MI heatmap...")

mi_data = results['analysis_5_mutual_information']
matrix = np.array(mi_data['matrix'])
positions = mi_data['positions']

# Create shorter labels
short_labels = []
for p in positions:
    if 'FR2_' in p:
        short_labels.append(f"F2.{p.split('_')[1]}")
    elif 'FR3_' in p:
        short_labels.append(f"F3.{p.split('_')[1]}")
    elif 'FR4_' in p:
        short_labels.append(f"F4.{p.split('_')[1]}")
    elif 'CDR1' in p:
        short_labels.append(f"C1.{p.split('_')[1][:3]}")
    elif 'CDR2' in p:
        short_labels.append(f"C2.{p.split('_')[1][:3]}")
    elif 'CDR3' in p:
        short_labels.append(f"C3.{p.split('_')[1][:3]}")
    else:
        short_labels.append(p[:6])

fig, ax = plt.subplots(figsize=(14, 12))
im = ax.imshow(matrix, cmap='YlOrRd', aspect='auto')

ax.set_xticks(range(len(short_labels)))
ax.set_yticks(range(len(short_labels)))
ax.set_xticklabels(short_labels, rotation=90, fontsize=8)
ax.set_yticklabels(short_labels, fontsize=8)

plt.colorbar(im, label='Mutual Information (bits)')
ax.set_title('Position Co-evolution in VHH Sequences\n(Higher MI = Stronger Coupling)', fontsize=14)

# Add region separators
fr2_end = sum(1 for p in positions if 'FR2_' in p)
fr3_end = fr2_end + sum(1 for p in positions if 'FR3_' in p)
fr4_end = fr3_end + sum(1 for p in positions if 'FR4_' in p)

for pos in [fr2_end, fr3_end, fr4_end]:
    ax.axhline(pos - 0.5, color='white', linewidth=2)
    ax.axvline(pos - 0.5, color='white', linewidth=2)

plt.tight_layout()
plt.savefig('epistasis_figures/1_mi_heatmap.png', dpi=150)
plt.close()
print("  Saved: epistasis_figures/1_mi_heatmap.png")

# =============================================================================
# 2. FR-CDR COUPLING NETWORK - What talks to what?
# =============================================================================
print("Creating FR-CDR coupling diagram...")

fig, ax = plt.subplots(figsize=(12, 8))

# Get top FR-CDR pairs
pairs = mi_data['top_pairs']
fr_cdr_pairs = [p for p in pairs if p['type'] == 'FR-CDR'][:15]

# Position coordinates for visualization
fr_positions = {'FR2_2': 0, 'FR2_9': 1, 'FR2_10': 2, 'FR2_12': 3, 'FR2_14': 4,
                'FR3_1': 5, 'FR3_6': 6, 'FR3_20': 7, 'FR4_3': 8}
cdr_positions = {'CDR1_first': 0, 'CDR1_2nd': 1, 'CDR1_last': 2, 'CDR1_2ndlast': 3,
                 'CDR2_first': 4, 'CDR2_2nd': 5, 'CDR2_last': 6, 'CDR2_2ndlast': 7,
                 'CDR3_first': 8, 'CDR3_2nd': 9, 'CDR3_last': 10, 'CDR3_2ndlast': 11, 'CDR3_3rdlast': 12}

# Draw FR positions (top row)
fr_y = 0.8
for name, x_idx in fr_positions.items():
    x = 0.1 + x_idx * 0.09
    color = '#2196F3' if 'FR2' in name else ('#4CAF50' if 'FR3' in name else '#FF9800')
    ax.scatter(x, fr_y, s=500, c=color, zorder=3)
    ax.text(x, fr_y + 0.08, name.replace('_', '\n'), ha='center', fontsize=8)

# Draw CDR positions (bottom row)  
cdr_y = 0.2
for name, x_idx in cdr_positions.items():
    x = 0.05 + x_idx * 0.07
    color = '#E91E63' if 'CDR1' in name else ('#9C27B0' if 'CDR2' in name else '#F44336')
    ax.scatter(x, cdr_y, s=500, c=color, zorder=3)
    label = name.replace('CDR', 'C').replace('_', '\n')
    ax.text(x, cdr_y - 0.1, label, ha='center', fontsize=7)

# Draw connections
for pair in fr_cdr_pairs:
    p1, p2 = pair['pos1'], pair['pos2']
    mi = pair['mi']
    
    # Figure out which is FR and which is CDR
    if 'FR' in p1:
        fr_pos, cdr_pos = p1, p2
    else:
        fr_pos, cdr_pos = p2, p1
    
    if fr_pos in fr_positions and cdr_pos in cdr_positions:
        x1 = 0.1 + fr_positions[fr_pos] * 0.09
        x2 = 0.05 + cdr_positions[cdr_pos] * 0.07
        
        # Line thickness based on MI
        linewidth = mi * 10
        alpha = min(0.9, mi * 3)
        
        ax.plot([x1, x2], [fr_y, cdr_y], 'k-', linewidth=linewidth, alpha=alpha)
        
        # Add MI value
        mid_x, mid_y = (x1 + x2) / 2, (fr_y + cdr_y) / 2
        ax.text(mid_x, mid_y, f'{mi:.2f}', fontsize=7, ha='center', 
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))

# Legend
legend_elements = [
    mpatches.Patch(color='#2196F3', label='FR2'),
    mpatches.Patch(color='#4CAF50', label='FR3'),
    mpatches.Patch(color='#FF9800', label='FR4'),
    mpatches.Patch(color='#E91E63', label='CDR1'),
    mpatches.Patch(color='#9C27B0', label='CDR2'),
    mpatches.Patch(color='#F44336', label='CDR3'),
]
ax.legend(handles=legend_elements, loc='upper right')

ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.axis('off')
ax.set_title('Framework ↔ CDR Coupling\n(Line thickness = Mutual Information strength)', fontsize=14)

plt.tight_layout()
plt.savefig('epistasis_figures/2_fr_cdr_coupling.png', dpi=150)
plt.close()
print("  Saved: epistasis_figures/2_fr_cdr_coupling.png")

# =============================================================================
# 3. VERNIER CLUSTER DISTRIBUTION
# =============================================================================
print("Creating Vernier cluster chart...")

clusters = results['analysis_2_vernier_clusters']

# Re-classify clusters properly
family_counts = defaultdict(int)
family_clusters = defaultdict(list)

for pattern_str, info in clusters.items():
    n = info['n']
    p = info['pattern']
    
    if p['FR2_12'] == 'L':
        fam = 'VH-like (L50)'
    elif p['FR2_12'] == 'W':
        fam = 'VHH-W52'
    elif p['FR2_12'] in 'GFA' and p['FR3_6'] == 'K':
        fam = 'VHH-standard (K71)'
    elif p['FR2_12'] in 'GFA':
        fam = 'VHH-other'
    else:
        fam = 'Other'
    
    family_counts[fam] += n
    family_clusters[fam].append((p, n))

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Pie chart
ax = axes[0]
labels = list(family_counts.keys())
sizes = list(family_counts.values())
colors = ['#2196F3', '#E91E63', '#4CAF50', '#FF9800', '#9E9E9E']
explode = [0.02] * len(labels)

ax.pie(sizes, labels=labels, autopct='%1.1f%%', colors=colors[:len(labels)],
       explode=explode, startangle=90)
ax.set_title(f'VHH Family Distribution\n(n={sum(sizes):,} sequences)', fontsize=12)

# Bar chart of top patterns
ax = axes[1]
top_clusters = sorted(clusters.items(), key=lambda x: -x[1]['n'])[:10]

patterns = []
counts = []
colors_bar = []

for pattern_str, info in top_clusters:
    p = info['pattern']
    # Create readable pattern string
    pat_str = f"F2.12={p['FR2_12']}, F3.6={p['FR3_6']}"
    patterns.append(pat_str)
    counts.append(info['n'])
    
    # Color by family
    if p['FR2_12'] == 'L':
        colors_bar.append('#2196F3')
    elif p['FR2_12'] == 'W':
        colors_bar.append('#E91E63')
    elif p['FR3_6'] == 'K':
        colors_bar.append('#4CAF50')
    else:
        colors_bar.append('#FF9800')

bars = ax.barh(range(len(patterns)), counts, color=colors_bar)
ax.set_yticks(range(len(patterns)))
ax.set_yticklabels(patterns, fontsize=9)
ax.set_xlabel('Number of Sequences')
ax.set_title('Top 10 Vernier Patterns', fontsize=12)
ax.invert_yaxis()

# Add count labels
for i, (bar, count) in enumerate(zip(bars, counts)):
    ax.text(bar.get_width() + 20000, bar.get_y() + bar.get_height()/2, 
            f'{count:,}', va='center', fontsize=9)

plt.tight_layout()
plt.savefig('epistasis_figures/3_vernier_clusters.png', dpi=150)
plt.close()
print("  Saved: epistasis_figures/3_vernier_clusters.png")

# =============================================================================
# 4. COMPENSATION RULES - What compensates for what?
# =============================================================================
print("Creating compensation rules chart...")

rules = results['analysis_1_compensation']['rules']

# Filter to substantial effects with good sample size
good_rules = [r for r in rules if r['n1'] >= 1000 and r['n2'] >= 1000]

# Separate by feature type
length_rules = [r for r in good_rules if 'length' in r['feature']][:20]
charge_rules = [r for r in good_rules if 'charge' in r['feature']][:20]

fig, axes = plt.subplots(1, 2, figsize=(14, 8))

# CDR3 Length compensation
ax = axes[0]
if length_rules:
    positions = [f"{r['position']}\n{r['res1']}→{r['res2']}" for r in length_rules]
    effects = [r['effect'] for r in length_rules]
    colors = ['#4CAF50' if e > 0 else '#F44336' for e in effects]
    
    bars = ax.barh(range(len(positions)), effects, color=colors)
    ax.set_yticks(range(len(positions)))
    ax.set_yticklabels(positions, fontsize=8)
    ax.axvline(0, color='black', linewidth=0.5)
    ax.set_xlabel('CDR3 Length Difference (amino acids)')
    ax.set_title('CDR3 Length Compensation\n(Green = longer CDR3, Red = shorter)', fontsize=11)
    ax.invert_yaxis()

# Charge compensation  
ax = axes[1]
if charge_rules:
    positions = [f"{r['position']}\n{r['res1']}→{r['res2']}" for r in charge_rules]
    effects = [r['effect'] for r in charge_rules]
    colors = ['#2196F3' if e > 0 else '#FF9800' for e in effects]
    
    bars = ax.barh(range(len(positions)), effects, color=colors)
    ax.set_yticks(range(len(positions)))
    ax.set_yticklabels(positions, fontsize=8)
    ax.axvline(0, color='black', linewidth=0.5)
    ax.set_xlabel('CDR3 Charge Difference')
    ax.set_title('CDR3 Charge Compensation\n(Blue = more positive, Orange = more negative)', fontsize=11)
    ax.invert_yaxis()

plt.tight_layout()
plt.savefig('epistasis_figures/4_compensation_rules.png', dpi=150)
plt.close()
print("  Saved: epistasis_figures/4_compensation_rules.png")

# =============================================================================
# 5. KEY INSIGHT SUMMARY
# =============================================================================
print("Creating summary figure...")

fig, ax = plt.subplots(figsize=(12, 8))
ax.axis('off')

summary_text = """
VHH EPISTASIS ANALYSIS SUMMARY
═══════════════════════════════════════════════════════════════════

DATASET: 10.76 million VHH sequences

FAMILY DISTRIBUTION (corrected):
  • VH-like (L50):      1.64M  (15%)  - Conventional VH framework
  • VHH-W52:            1.87M  (17%)  - Tryptophan at position 52
  • VHH-standard (K71):  502K   (5%)  - Classic camelid with Lys71
  • VHH-other:          5.77M  (54%)  - VHH variants

KEY FINDINGS:

1. STRONGEST FR↔CDR COUPLINGS (Mutual Information):
   • FR2_12 (pos 50) ↔ CDR1_2nd:  MI = 0.38  ← The hallmark position!
   • FR4_3 ↔ CDR3_3rdlast:        MI = 0.28  ← FR4 matters for CDR3!
   • FR2_12 ↔ CDR2_2nd:           MI = 0.26  ← Position 50 is central

2. COMPENSATION PATTERNS:
   • FR3_17: K vs L → CDR3 length +10 aa difference
   • FR3_30: E vs F → CDR3 length +9.7 aa difference
   These framework positions predict CDR3 length!

3. HIGHER-ORDER RULES (>95% confidence):
   • short CDR3 + FR2_12=F → FR3_6=V (96%)
   • short CDR3 + FR3_6=K → FR3_12=S (95%)
   Framework positions are interdependent with CDR properties

IMPLICATIONS FOR HUMANIZATION:
   Position 50 (FR2_12) is the master switch - it correlates with:
   - CDR1 composition
   - CDR2 composition  
   - Other FR2 positions (MI > 0.5 with FR2_2, FR2_9, FR2_10)
   
   When changing pos 50, expect to need CDR adjustments!
"""

ax.text(0.05, 0.95, summary_text, transform=ax.transAxes, fontsize=10,
        verticalalignment='top', fontfamily='monospace',
        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

plt.savefig('epistasis_figures/5_summary.png', dpi=150, bbox_inches='tight')
plt.close()
print("  Saved: epistasis_figures/5_summary.png")

print("\n" + "="*60)
print("All figures saved to epistasis_figures/")
print("="*60)
