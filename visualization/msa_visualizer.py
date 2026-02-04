#!/usr/bin/env python3
"""
align_vs_lead_clear3_antpack_legend_v7.py

v7 Changes:
- Switched from Clustal Omega to MAFFT --localpair for better CDR alignment
- MAFFT's local pairwise alignment keeps motifs together (e.g., DYD aligns with DYD)
- Falls back to Clustal Omega if MAFFT is not installed

Install MAFFT: sudo apt install mafft

Updated to handle sorting even when no CDRs are detected.
Falls back to framework-only similarity scoring.
"""

import os, sys, shutil, subprocess, argparse, csv
from collections import Counter, defaultdict
import math

os.environ.setdefault("MPLBACKEND", "Agg")  # headless
from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from pymsaviz import MsaViz
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.patches import Rectangle
import numpy as np
from matplotlib.backends.backend_agg import FigureCanvasAgg
import matplotlib.colors as mcolors

# --- Per-AA palette (Zappo-like) ---
AA_PALETTE = {
    "A":"#C8C8C8","R":"#145AFF","N":"#00DCDC","D":"#E60A0A","C":"#E6E600",
    "Q":"#00DCDC","E":"#E60A0A","G":"#EBEBEB","H":"#8282D2","I":"#0F820F",
    "L":"#0F820F","K":"#145AFF","M":"#E6E600","F":"#3232AA","P":"#DC9682",
    "S":"#FA9600","T":"#FA9600","W":"#B45AB4","Y":"#3232AA","V":"#0F820F"
}

# --- Group palette (6 distinct, bright, no grey) ---
GROUP_PALETTE = {
    "hydrophobic": "#FF7F00",  # A V I L M (orange)
    "polar"      : "#4DAF4A",  # S T N Q C (green)
    "aromatic"   : "#984EA3",  # F W Y (purple)
    "positive"   : "#E41A1C",  # K R H (red)
    "negative"   : "#377EB8",  # D E (blue)
    "special"    : "#FFD92F",  # G P (yellow)
}

# Choose ONE primary group per AA (anchor for similarity vs lead)
PRIMARY_GROUP = {
    "A":"hydrophobic","V":"hydrophobic","I":"hydrophobic","L":"hydrophobic","M":"hydrophobic",
    "F":"aromatic","W":"aromatic","Y":"aromatic",
    "C":"polar","S":"polar","T":"polar","N":"polar","Q":"polar",
    "K":"positive","R":"positive","H":"positive",
    "D":"negative","E":"negative",
    "G":"special","P":"special",
}

def lighten_hex(hex_color: str, amt: float) -> str:
    hex_color = hex_color.lstrip("#")
    r = int(hex_color[0:2], 16); g = int(hex_color[2:4], 16); b = int(hex_color[4:6], 16)
    r = int(r + (255 - r) * amt); g = int(g + (255 - g) * amt); b = int(b + (255 - b) * amt)
    return f"#{r:02x}{g:02x}{b:02x}"

def similar_group(a: str, b: str) -> bool:
    ga = PRIMARY_GROUP.get(a.upper()); gb = PRIMARY_GROUP.get(b.upper())
    return ga is not None and gb is not None and ga == gb

def run_clustalo(in_fa: str, out_sto: str):
    """Legacy Clustal Omega alignment - kept for fallback"""
    exe = shutil.which("clustalo") or sys.exit("ERROR: clustalo not found. sudo apt update && sudo apt install -y clustalo")
    cmd = [exe, "-i", in_fa, "-o", out_sto, "--outfmt=st", "--force"]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stdout); print(r.stderr, file=sys.stderr)
        sys.exit(f"ERROR: clustalo failed ({r.returncode})")
    if not os.path.exists(out_sto) or os.path.getsize(out_sto) == 0:
        sys.exit("ERROR: empty alignment.")

def run_mafft(in_fa: str, out_fa: str):
    """
    MAFFT with --localpair for better CDR alignment.
    This maximizes local matches and keeps motifs together.
    """
    exe = shutil.which("mafft")
    if not exe:
        print("[WARN] mafft not found, falling back to clustalo. Install with: sudo apt install mafft", file=sys.stderr)
        # Fall back to clustalo, but output to fasta format
        out_sto = out_fa.replace('.fa', '.sto')
        run_clustalo(in_fa, out_sto)
        # Convert stockholm to fasta
        msa = AlignIO.read(out_sto, "stockholm")
        AlignIO.write(msa, out_fa, "fasta")
        return
    
    # --localpair: All pairwise local alignments, more accurate for variable regions like CDRs
    # --maxiterate 1000: More refinement iterations
    cmd = [exe, "--localpair", "--maxiterate", "1000", "--quiet", in_fa]
    
    print(f"[INFO] Running MAFFT with --localpair for better CDR alignment...", file=sys.stderr)
    r = subprocess.run(cmd, capture_output=True, text=True)
    
    if r.returncode != 0:
        print(f"[WARN] MAFFT failed: {r.stderr}", file=sys.stderr)
        print("[WARN] Falling back to clustalo...", file=sys.stderr)
        out_sto = out_fa.replace('.fa', '.sto')
        run_clustalo(in_fa, out_sto)
        msa = AlignIO.read(out_sto, "stockholm")
        AlignIO.write(msa, out_fa, "fasta")
        return
    
    # Write MAFFT output (comes on stdout)
    with open(out_fa, 'w') as f:
        f.write(r.stdout)
    
    if not os.path.exists(out_fa) or os.path.getsize(out_fa) == 0:
        sys.exit("ERROR: empty alignment from MAFFT.")

def lead_reference_per_column(msa, lead_rows):
    L = msa.get_alignment_length()
    ref = []
    
    # If no leads, use consensus
    if not lead_rows:
        for j in range(L):
            col = [msa[i].seq[j] for i in range(len(msa))]
            col = [aa for aa in col if aa not in "-.*?"]
            ref.append(Counter(col).most_common(1)[0][0] if col else None)
        return ref
    
    # Original lead-based logic
    for j in range(L):
        col = [msa[r].seq[j] for r in lead_rows]
        col = [aa for aa in col if aa not in "-.*?"]
        ref.append(Counter(col).most_common(1)[0][0] if col else None)
    return ref

def strict_similarity(a: str, b: str) -> bool:
    """Strict similarity = same PRIMARY_GROUP as lead residue"""
    ga = PRIMARY_GROUP.get(a.upper()); gb = PRIMARY_GROUP.get(b.upper())
    return ga is not None and gb is not None and ga == gb

# ---------- AntPack CSV integration ----------
def read_antpack_cdr_csv(csv_path):
    rows = []
    with open(csv_path, newline="") as f:
        r = csv.reader(f)
        header = next(r)
        idx = {h.lower().strip(): i for i, h in enumerate(header)}
        for line in r:
            row = {k: (line[v] if v < len(line) else "") for k, v in idx.items()}
            rows.append(row)
    return rows

def map_cdr_seq_to_alignment_cols(aligned_seq, cdr_seq):
    if not cdr_seq:
        return None
    AA = set("ACDEFGHIKLMNPQRSTVWY")

    s_aln = str(aligned_seq).upper()
    ungapped = "".join(ch for ch in s_aln if ch in AA)

    cdr = "".join(ch for ch in (cdr_seq or "").upper() if ch in AA)
    if not cdr:
        return None

    start = ungapped.find(cdr)
    if start == -1:
        return None

    pos2col = {}
    p = 0
    for j, ch in enumerate(s_aln):
        if ch in AA:
            p += 1
            pos2col[p] = j + 1

    a = pos2col.get(start + 1)
    b = pos2col.get(start + len(cdr))
    return (a, b) if a and b and a <= b else None

def find_best_cdr_mapping(msa, cdr_rows, lead_rows):
    """
    Return (best_seq_idx, best_row_idx, mapped, hits)
    where mapped = {'CDR1': (a,b), ...} in 1-based alignment columns.
    """
    def get(row, k):
        return row.get(k, row.get(k.lower(), ""))

    def score_row_on_seq(aligned_seq, row):
        mapped = {}
        hits = 0
        for tag, key in (("CDR1","cdr1_seq"), ("CDR2","cdr2_seq"), ("CDR3","cdr3_seq")):
            s = (get(row, key) or "").strip().upper()
            if not s:
                continue
            rng = map_cdr_seq_to_alignment_cols(aligned_seq, s)
            if rng:
                mapped[tag] = rng
                hits += 1
        return hits, mapped

    best_hits = -1
    best = (None, None, {}, 0)

    for seq_idx, rec in enumerate(msa):
        for row_idx, row in enumerate(cdr_rows):
            hits, mapped = score_row_on_seq(msa[seq_idx].seq, row)
            if hits > best_hits:
                best_hits = hits
                best = (seq_idx, row_idx, mapped, hits)
            elif hits == best_hits and hits > 0:
                if lead_rows and seq_idx == lead_rows[0] and best[0] != lead_rows[0]:
                    best = (seq_idx, row_idx, mapped, hits)

    return best

def _split_range_across_wrap(a, b, W):
    if not W or W <= 0:
        return [(a, b)]
    out = []
    i = a
    while i <= b:
        k = (i - 1) // W
        row_end = (k + 1) * W
        j = min(b, row_end)
        a_loc = i - k * W
        b_loc = j - k * W
        out.append((a_loc, b_loc))
        i = j + 1
    return out

def split_mapped_for_wrap(mapped, wrap_len):
    try:
        W = int(wrap_len)
    except Exception:
        W = None
    if not W or W <= 0:
        return mapped

    out = {}
    for name, (a, b) in mapped.items():
        segs = _split_range_across_wrap(a, b, W)
        if len(segs) == 1:
            out[name] = segs[0]
        else:
            for idx, (aa, bb) in enumerate(segs, 1):
                out[f"{name}_part{idx}"] = (aa, bb)
    return out

def draw_cdr_boxes(ax, mapped_ranges, fill="#7b777d", alpha=0.25, z=8,
                   edge=None, outline_lw=1.2, linestyle="-",
                   cell_left=-1.0, cell_right=0):
    """Draw CDR boxes spanning whole residue cells."""
    if not mapped_ranges:
        return
        
    if edge is None:
        edge = fill

    y0, y1 = ax.get_ylim()
    if y0 > y1:
        y0, y1 = y1, y0
    pad_y = 0.0

    for name, (a, b) in mapped_ranges.items():
        left  = (a + cell_left)
        right = (b + cell_right)
        x = left
        w = right - left
        h = (y1 - y0) - 2 * pad_y

        ax.add_patch(Rectangle((x, y0 + pad_y), w, h,
                               facecolor=fill, alpha=alpha, lw=0, zorder=z))
        ax.add_patch(Rectangle((x, y0 + pad_y), w, h,
                               fill=False, edgecolor=edge,
                               linewidth=outline_lw, linestyle=linestyle,
                               zorder=z+1))

def compose_with_top_legend(fig_main, pale_amt=0.4):
    """
    Compose a new figure with a boxed legend at the top and the alignment below.
    """
    import matplotlib.colors as mcolors
    ax1 = fig_main.get_axes()[0]

    # Create colors with spacing between pairs
    colors = []
    labels = []
    group_data = [("hydrophobic", "AVILM"), ("polar", "STNQC"),
                ("aromatic", "FWY"), ("positive", "KRH"),
                ("negative", "DE"), ("special", "GP")]

    for i, (group_name, letters) in enumerate(group_data):
        colors.extend([GROUP_PALETTE[group_name], lighten_hex(GROUP_PALETTE[group_name], pale_amt)])
        labels.append(f"{group_name}: {letters}")

        if i < len(group_data) - 1:
            colors.append('white')

    cmap = mcolors.ListedColormap(colors)
    norm = mcolors.BoundaryNorm(range(len(colors) + 1), len(colors))

    dummy_data = np.arange(len(colors)).reshape(1, -1)
    im1 = ax1.imshow(dummy_data, cmap=cmap, norm=norm, aspect='auto')

    pos = ax1.get_position()
    
    # Scale all spacing values to maintain constant absolute size
    reference_height = 5.0  # Adjust this to match your 1-2 sequence figure height
    fig_height = fig_main.get_figheight()
    scale_factor = reference_height / fig_height
    
    # Gap between alignment and legend
    legend_gap = 0.15 * scale_factor
    # Legend height
    legend_height = 0.15 * scale_factor
    # Gap between legend and title
    title_gap = 0.25 * scale_factor  # Adjust this value to control title-legend spacing
    
    cax1 = fig_main.add_axes([pos.x0, pos.y1 + legend_gap, 
                              pos.width, legend_height])
    cb1 = fig_main.colorbar(im1, cax=cax1, orientation="horizontal", 
                           ticklocation='top', drawedges=False)

    tick_positions = []
    for i in range(len(group_data)):
        pair_start = i * 3
        tick_positions.append(pair_start + 1)

    cb1.set_ticks(tick_positions)
    cb1.set_ticklabels(labels)
    cb1.outline.set_visible(False)

    fig_main.suptitle("Groups — bright=dissimilar, pale=similar", 
                     x=0.5, y=cax1.get_position().y1 + title_gap)

def bolden_lead_rows(fig, msa, lead_rows):
    if not fig.axes: return
    ax = fig.axes[0]
    wanted_ids = set()
    for i in lead_rows:
        original_id = msa[i].id.split('__')[0]
        wanted_ids.add(original_id)
    
    row_y = []
    had_labels = False

    for t in fig.findobj(matplotlib.text.Text):
        s = t.get_text()
        for wanted_id in wanted_ids:
            if s.startswith(wanted_id):
                y = t.get_position()[1]
                row_y.append(y); had_labels = True
                t.set_fontweight("bold")
                t.set_path_effects([pe.Stroke(linewidth=1.6, foreground="white"), pe.Normal()])
                break

    if not had_labels:
        y0, y1 = ax.get_ylim()
        if y0 > y1: y0, y1 = y1, y0
        for i in lead_rows:
            row_y.append(y1 - (i + 1))

    tol = 0.25
    AA_CHARS = set("ACDEFGHIKLMNPQRSTVWY-.*?")
    for t in fig.findobj(matplotlib.text.Text):
        s = t.get_text()
        if len(s) == 1 and s in AA_CHARS:
            y = t.get_position()[1]
            if any(abs(y - yy) < tol for yy in row_y):
                t.set_fontweight("bold")
                t.set_path_effects([pe.Stroke(linewidth=0.9, foreground="white"), pe.Normal()])

def get_cdr_columns_from_mapped(mapped_ranges):
    cdr_cols = set()
    for name, (a, b) in mapped_ranges.items():
        if "CDR" in name.upper():
            for col in range(a, b + 1):
                cdr_cols.add(col)
    return cdr_cols

def calculate_similarity_score_corrected(seq_record, ref_chars, ref_record, cdr_cols, weights, cdr3_cols=None):
    """
    Enhanced version with CDR3-specific weighting.
    Now handles case where cdr_cols is empty (no CDRs detected).
    """
    seq_str = str(seq_record.seq)
    ref_str = str(ref_record.seq) if ref_record else ""
    L = len(seq_str)
    
    has_cdrs = bool(cdr_cols)
    
    cdr_dissimilar = 0
    cdr_similar = 0
    cdr3_dissimilar = 0
    cdr3_similar = 0
    fr_dissimilar = 0
    fr_similar = 0
    mutation_positions = []
    
    cdr3_multiplier = weights.get("cdr3_multiplier", 1.0)
    has_cdr3_weight = cdr3_multiplier > 1.0 and cdr3_cols and has_cdrs
    
    # Count mutations by type and location
    for j in range(min(L, len(ref_chars))):
        aa = seq_str[j].upper()
        ref_aa = ref_chars[j]
        
        if aa in "-.*?" or ref_aa is None or aa == ref_aa:
            continue
            
        mutation_positions.append(j + 1)
        is_similar = similar_group(aa, ref_aa)
        
        if has_cdrs:
            # Normal CDR-based scoring
            is_cdr = (j + 1) in cdr_cols
            is_cdr3 = has_cdr3_weight and (j + 1) in cdr3_cols
            
            if is_cdr:
                if is_cdr3:
                    if is_similar:
                        cdr3_similar += 1
                    else:
                        cdr3_dissimilar += 1
                else:
                    if is_similar:
                        cdr_similar += 1
                    else:
                        cdr_dissimilar += 1
            else:
                if is_similar:
                    fr_similar += 1
                else:
                    fr_dissimilar += 1
        else:
            # No CDRs - treat all as framework
            if is_similar:
                fr_similar += 1
            else:
                fr_dissimilar += 1
    
    # CDR length penalty calculation (only if CDRs exist)
    cdr_length_penalty = 0
    if has_cdrs and ref_record:
        seq_cdr_gaps = sum(1 for j in range(L) 
                          if (j + 1) in cdr_cols and seq_str[j] in "-.*?")
        ref_cdr_gaps = sum(1 for j in range(len(ref_str)) 
                          if (j + 1) in cdr_cols and ref_str[j] in "-.*?")
        
        seq_cdr_residues = sum(1 for j in range(L) 
                              if (j + 1) in cdr_cols and seq_str[j] not in "-.*?")
        ref_cdr_residues = sum(1 for j in range(len(ref_str)) 
                              if (j + 1) in cdr_cols and ref_str[j] not in "-.*?")
        
        gap_diff = abs(seq_cdr_gaps - ref_cdr_gaps)
        residue_diff = abs(seq_cdr_residues - ref_cdr_residues)
        cdr_length_penalty = gap_diff + residue_diff
    
    # Overall length penalty (always calculate)
    overall_length_penalty = 0
    if ref_record:
        seq_gaps = sum(1 for ch in seq_str if ch in "-.*?")
        ref_gaps = sum(1 for ch in ref_str if ch in "-.*?")
        overall_length_penalty = abs(seq_gaps - ref_gaps)
    
    # Calculate weighted score
    if has_cdrs:
        # Normal CDR-aware scoring
        score = (weights["cdr_dissimilar"] * cdr_dissimilar +
                 weights["cdr_similar"] * cdr_similar +
                 weights["cdr_dissimilar"] * cdr3_dissimilar * cdr3_multiplier +
                 weights["cdr_similar"] * cdr3_similar * cdr3_multiplier +
                 weights["fr_dissimilar"] * fr_dissimilar +
                 weights["fr_similar"] * fr_similar +
                 weights["cdr_length"] * cdr_length_penalty)
    else:
        # No CDRs - use framework-only scoring with length penalty
        score = (weights["fr_dissimilar"] * fr_dissimilar +
                 weights["fr_similar"] * fr_similar +
                 weights.get("overall_length", 5.0) * overall_length_penalty)
    
    # Tie breakers
    total_mutations = len(mutation_positions)
    earliest_mutation = min(mutation_positions) if mutation_positions else float('inf')
    
    # Build details string
    if has_cdrs:
        if has_cdr3_weight and (cdr3_dissimilar > 0 or cdr3_similar > 0):
            details = (f"CDR3_dis:{cdr3_dissimilar}, CDR3_sim:{cdr3_similar}, "
                      f"CDR_dis:{cdr_dissimilar}, CDR_sim:{cdr_similar}, "
                      f"FR_dis:{fr_dissimilar}, FR_sim:{fr_similar}, CDR_gaps_Δ:{cdr_length_penalty}")
        else:
            details = (f"CDR_dis:{cdr_dissimilar + cdr3_dissimilar}, "
                      f"CDR_sim:{cdr_similar + cdr3_similar}, "
                      f"FR_dis:{fr_dissimilar}, FR_sim:{fr_similar}, CDR_gaps_Δ:{cdr_length_penalty}")
    else:
        # No CDRs detected
        details = (f"Dissimilar:{fr_dissimilar}, Similar:{fr_similar}, "
                  f"Length_Δ:{overall_length_penalty}")
    
    return {
        'score': score,
        'total_mutations': total_mutations,
        'earliest_mutation': earliest_mutation,
        'cdr_dissimilar': cdr_dissimilar,
        'cdr_similar': cdr_similar,
        'cdr3_dissimilar': cdr3_dissimilar,
        'cdr3_similar': cdr3_similar,
        'fr_dissimilar': fr_dissimilar,
        'fr_similar': fr_similar,
        'cdr_length_penalty': cdr_length_penalty,
        'overall_length_penalty': overall_length_penalty,
        'details': details
    }

def embed_scores_in_sequence_ids(msa, score_display_info):
    """
    Create new MSA with scores embedded directly in sequence IDs.
    """
    from Bio.Align import MultipleSeqAlignment
    
    id_to_score = {info['sequence_id']: info for info in score_display_info}
    
    new_records = []
    for record in msa:
        original_id = record.id
        
        if original_id in id_to_score:
            score_info = id_to_score[original_id]
            if score_info['is_lead']:
                new_id = f"{original_id}__LEAD"
            else:
                new_id = f"{original_id}__{score_info['score']:.1f}"
        else:
            new_id = original_id
        
        new_record = SeqRecord(record.seq, id=new_id, description=record.description)
        new_records.append(new_record)
    
    return MultipleSeqAlignment(new_records)

def sort_msa_by_similarity_corrected(msa, lead_rows, ref_chars, cdr_cols, weights, mapped_ranges=None):
    """
    Enhanced version that handles CDR3-specific weighting.
    Now works even when cdr_cols is empty (no CDRs detected).
    """
    ref_record = msa[lead_rows[0]] if lead_rows else msa[0]
    
    has_cdrs = bool(cdr_cols)
    
    # Get CDR3 columns if CDR3 multiplier is set and we have CDRs
    cdr3_cols = None
    if has_cdrs and weights.get("cdr3_multiplier", 1.0) > 1.0 and mapped_ranges:
        cdr3_cols = get_cdr3_columns(mapped_ranges)
        if cdr3_cols:
            print(f"[INFO] Applying {weights['cdr3_multiplier']}x weight to {len(cdr3_cols)} CDR3 positions", file=sys.stderr)
    
    if not has_cdrs:
        print("[INFO] No CDRs detected - using framework-only similarity scoring", file=sys.stderr)
    
    # Calculate scores for all sequences
    seq_scores = []
    for i, record in enumerate(msa):
        is_lead = i in lead_rows
        
        if is_lead:
            score_info = {
                'score': -1000,
                'total_mutations': 0,
                'earliest_mutation': float('inf'),
                'cdr_dissimilar': 0,
                'cdr_similar': 0,
                'cdr3_dissimilar': 0,
                'cdr3_similar': 0,
                'fr_dissimilar': 0,
                'fr_similar': 0,
                'cdr_length_penalty': 0,
                'overall_length_penalty': 0,
                'details': "LEAD"
            }
        else:
            score_info = calculate_similarity_score_corrected(record, ref_chars, ref_record, cdr_cols, weights, cdr3_cols)
        
        seq_scores.append((i, record, score_info, is_lead))
    
    # Sort
    def sort_key(item):
        i, record, score_info, is_lead = item
        return (
            score_info['score'],
            score_info['total_mutations'],
            score_info['earliest_mutation'],
            record.id.lower()
        )
    
    seq_scores.sort(key=sort_key)
    
    # Create new MSA
    sorted_records = [item[1] for item in seq_scores]
    sorted_msa = MultipleSeqAlignment(sorted_records)
    
    # Update lead_rows indices
    new_lead_rows = []
    for new_idx, (old_idx, record, score_info, is_lead) in enumerate(seq_scores):
        if is_lead:
            new_lead_rows.append(new_idx)
    
    # Prepare score info
    score_display_info = []
    for new_idx, (old_idx, record, score_info, is_lead) in enumerate(seq_scores):
        score_display_info.append({
            'sequence_id': record.id,
            'score': score_info['score'],
            'details': score_info['details'],
            'is_lead': is_lead
        })
    
    return sorted_msa, new_lead_rows, score_display_info

def get_cdr3_columns(mapped_ranges):
    """Extract just the CDR3 columns from mapped ranges."""
    cdr3_cols = set()
    for name, (a, b) in mapped_ranges.items():
        if "CDR3" in name.upper():
            for col in range(a, b + 1):
                cdr3_cols.add(col)
    return cdr3_cols

def prompt_for_sorting(has_cdrs=True):
    """
    Enhanced prompt that explains what will be compared based on CDR detection.
    """
    print("\n" + "="*80)
    print("SEQUENCE SORTING OPTIONS")
    print("="*80)
    
    if has_cdrs:
        print("\nWould you like to sort sequences by similarity to the lead/consensus?")
        print("• If 'n': Keep sequences in upload order")
        print("• If 'y': Sort by CDR-aware similarity scoring")
    else:
        print("\n⚠️  No CDR regions were detected in the sequences.")
        print("\nWould you like to sort sequences by overall similarity?")
        print("• If 'n': Keep sequences in upload order")
        print("• If 'y': Sort by overall mutation similarity and length differences")
    print()
    
    # Get basic y/n choice
    while True:
        choice = input("Sort sequences by similarity? (y/n): ").strip().lower()
        if choice in ['y', 'yes']:
            return True, prompt_for_weights(has_cdrs)
        elif choice in ['n', 'no']:
            return False, None
        else:
            print("Please enter 'y' or 'n'")

def prompt_for_weights(has_cdrs=True):
    """
    Enhanced version with different prompts based on CDR detection.
    """
    if has_cdrs:
        print("\nSORTING ALGORITHM (CDR-aware):")
        print("The algorithm prioritizes mutations in this order:")
        print("1. CDR dissimilar mutations (different physicochemical groups)")
        print("2. CDR similar mutations (same physicochemical groups)")
        print("3. Framework dissimilar mutations")
        print("4. Framework similar mutations")
        print("5. CDR length differences (gaps and residue count differences)")
    else:
        print("\nSORTING ALGORITHM (Framework-only):")
        print("Since no CDRs were detected, the algorithm will use:")
        print("1. Overall dissimilar mutations (different physicochemical groups)")
        print("2. Overall similar mutations (same physicochemical groups)")
        print("3. Overall length differences (gaps)")
    print()
    
    if has_cdrs:
        default_weights = {
            "cdr_dissimilar": 16.0,
            "cdr_similar": 3.0,
            "fr_dissimilar": 0.5,
            "fr_similar": 0.1,
            "cdr_length": 100.0,
            "cdr3_multiplier": 1.0
        }
    else:
        default_weights = {
            "cdr_dissimilar": 0.0,  # Not used
            "cdr_similar": 0.0,      # Not used
            "fr_dissimilar": 10.0,   # Main weight
            "fr_similar": 2.0,       # Secondary weight
            "cdr_length": 0.0,       # Not used
            "cdr3_multiplier": 1.0,  # Not used
            "overall_length": 5.0    # New weight for overall length
        }
    
    print("Default weights:")
    if has_cdrs:
        for key, value in default_weights.items():
            if key != "cdr3_multiplier":
                print(f"  {key.replace('_', ' ').title()}: {value}")
    else:
        print(f"  Dissimilar Mutations: {default_weights['fr_dissimilar']}")
        print(f"  Similar Mutations: {default_weights['fr_similar']}")
        print(f"  Overall Length: {default_weights['overall_length']}")
    print()
    
    # Ask about CDR3 emphasis only if we have CDRs
    if has_cdrs:
        print("CDR3 WEIGHTING:")
        print("CDR3 is often the most important region for antibody specificity.")
        cdr3_response = input("Apply extra weight to CDR3 mutations? (y/n, default=n): ").strip().lower()
        
        if cdr3_response in ['y', 'yes']:
            while True:
                multiplier_input = input("CDR3 weight multiplier (e.g., 2.0 = double weight, default=2.0): ").strip()
                if not multiplier_input:
                    default_weights["cdr3_multiplier"] = 2.0
                    break
                try:
                    default_weights["cdr3_multiplier"] = float(multiplier_input)
                    if default_weights["cdr3_multiplier"] < 1.0:
                        print("  Multiplier should be >= 1.0. Using 1.0 (no extra weight).")
                        default_weights["cdr3_multiplier"] = 1.0
                    break
                except ValueError:
                    print("  Please enter a valid number")
            
            print(f"[INFO] CDR3 mutations will be weighted {default_weights['cdr3_multiplier']}x more than other CDR mutations")
    
    print("\nTie-breakers (in order):")
    print("1. Fewer total mutations")
    print("2. Earlier mutation position")
    print("3. Alphabetical by sequence ID")
    print()
    
    use_custom = input("Customize weights? (y/n, default=n): ").strip().lower()
    
    if use_custom in ['y', 'yes']:
        weights = {"cdr3_multiplier": default_weights.get("cdr3_multiplier", 1.0)}
        print("\nEnter custom weights (press Enter to use default):")
        
        for key, default_val in default_weights.items():
            if key == "cdr3_multiplier":
                continue
            # Skip irrelevant weights based on CDR detection
            if not has_cdrs and key in ["cdr_dissimilar", "cdr_similar", "cdr_length"]:
                weights[key] = 0.0
                continue
            if has_cdrs and key == "overall_length":
                continue
                
            while True:
                user_input = input(f"  {key.replace('_', ' ').title()} [{default_val}]: ").strip()
                if not user_input:
                    weights[key] = default_val
                    break
                try:
                    weights[key] = float(user_input)
                    break
                except ValueError:
                    print("    Please enter a valid number")
        
        # Ensure all required keys exist
        for key in default_weights:
            if key not in weights:
                weights[key] = default_weights[key]
        
        return weights
    else:
        return default_weights

def export_sorting_results(msa, score_display_info, output_path, weights):
    """
    Export sequence order and scores to CSV.
    """
    import csv
    
    rows = []
    for idx, info in enumerate(score_display_info):
        seq_id = info['sequence_id']
        
        if info['is_lead']:
            rows.append({
                'rank': idx + 1,
                'sequence_id': seq_id,
                'score': 'LEAD',
                'details': info['details']
            })
        else:
            rows.append({
                'rank': idx + 1,
                'sequence_id': seq_id,
                'score': f"{info['score']:.2f}",
                'details': info['details']
            })
    
    # Write CSV
    with open(output_path, 'w', newline='') as f:
        fieldnames = ['rank', 'sequence_id', 'score', 'details']
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    
    # Also write weights used
    weights_path = output_path.replace('.csv', '_weights.txt')
    with open(weights_path, 'w') as f:
        f.write("Weights used for sorting:\n")
        f.write("=" * 50 + "\n")
        for key, val in weights.items():
            if val > 0:  # Only show weights that were actually used
                f.write(f"{key.replace('_', ' ').title()}: {val}\n")
    
    print(f"\n[INFO] Exported sorting results to: {output_path}")
    print(f"[INFO] Exported weights to: {weights_path}")


# =========================
# IMGT numbering + critical-position markers under consensus
# =========================
# Hallmark IMGT positions (VHH-critical)
HALLMARK_IMGT_POSITIONS = {42, 49, 50, 52}

# Vernier IMGT positions (edit this list to match your definition)
VERNIER_IMGT_POSITIONS = {66, 67, 68, 69, 71, 76, 78, 82, 87, 89, 91, 94}

AA_SET = set("ACDEFGHIKLMNPQRSTVWY")

def _antpack_imgt_numbering(ungapped_seq: str):
    """Return list of IMGT numbers aligned to ungapped_seq indices (0-based), or None."""
    try:
        from antpack.numbering_tools.single_chain_annotator import SingleChainAnnotator
    except Exception:
        return None

    annot = SingleChainAnnotator(chains=["H"], scheme="imgt")
    numbering, pid, chain, err = annot.analyze_seq(ungapped_seq)
    if err or numbering is None:
        return None

    # Prefer trim_alignment if available: gives a per-residue list of IMGT numbers for the trimmed sequence.
    try:
        trimmed_seq, trimmed_nums, _, _ = annot.trim_alignment(ungapped_seq, (numbering, pid, chain, err))
        if trimmed_nums and len(trimmed_seq) == len(trimmed_nums):
            out = []
            for n in trimmed_nums:
                try:
                    if isinstance(n, (tuple, list)):
                        n = n[0]
                    out.append(int(str(n).strip()))
                except Exception:
                    out.append(None)

            # Map trimmed back onto full ungapped_seq by subsequence match (common: only ends trimmed)
            idx = ungapped_seq.find(trimmed_seq)
            if idx != -1:
                full = [None] * len(ungapped_seq)
                for i, n in enumerate(out):
                    full[idx + i] = n
                return full
    except Exception:
        pass

    # Fallback: extract integer-like positions from numbering, if length matches.
    nums = []
    for n in numbering:
        try:
            if isinstance(n, (tuple, list)):
                n = n[0]
            nums.append(int(str(n).strip()))
        except Exception:
            nums.append(None)

    if len(nums) == len(ungapped_seq):
        return nums
    return None


def build_alignmentcol_to_imgt_from_lead(msa: MultipleSeqAlignment, lead_row_idx: int):
    """Map 1-based alignment columns to IMGT numbers using lead row AntPack IMGT numbering."""
    lead_aln = str(msa[lead_row_idx].seq).upper()
    ungapped = "".join(ch for ch in lead_aln if ch in AA_SET)

    imgt_per_res = _antpack_imgt_numbering(ungapped)
    if not imgt_per_res:
        return {}

    col_to_imgt = {}
    residx = 0
    for col, ch in enumerate(lead_aln, start=1):
        if ch in AA_SET:
            residx += 1
            if residx - 1 < len(imgt_per_res):
                n = imgt_per_res[residx - 1]
                if isinstance(n, int):
                    col_to_imgt[col] = n
    return col_to_imgt


def annotate_imgt_and_markers(ax, msa: MultipleSeqAlignment, lead_rows, wrap_len: int,
                              hallmark_positions=HALLMARK_IMGT_POSITIONS,
                              vernier_positions=VERNIER_IMGT_POSITIONS,
                              show_every: int = 10):
    """
    Draw IMGT position numbers below the consensus line, plus:
      H (red) at hallmark positions
      V (blue) at vernier positions

    Note: for best alignment, run with a large wrap_len (e.g. 10000) so wrapping is effectively off.
    """
    if not lead_rows:
        return
    if wrap_len and wrap_len < 200:
        print("[WARN] IMGT annotation is best with large --wrap (e.g. 10000).", file=sys.stderr)

    col_to_imgt = build_alignmentcol_to_imgt_from_lead(msa, lead_rows[0])
    if not col_to_imgt:
        print("[WARN] Could not build IMGT mapping from lead (AntPack missing or numbering failed).", file=sys.stderr)
        return

    y0, y1 = ax.get_ylim()
    y_bottom = max(y0, y1)

    # Two annotation rows below the alignment
    y_nums = y_bottom + 0.70
    y_mark = y_bottom + 1.25

    important = set(hallmark_positions) | set(vernier_positions)

    for col, imgt in col_to_imgt.items():
        # Numbers row: every N, plus always at important sites
        if (show_every and imgt % show_every == 0) or (imgt in important):
            ax.text(col, y_nums, str(imgt),
                    ha="center", va="center", fontsize=5,
                    color="black", clip_on=False)

        # Marker row
        if imgt in hallmark_positions:
            ax.text(col, y_mark, "H",
                    ha="center", va="center", fontsize=6, fontweight="bold",
                    color="red", clip_on=False)
        elif imgt in vernier_positions:
            ax.text(col, y_mark, "V",
                    ha="center", va="center", fontsize=6, fontweight="bold",
                    color="#145AFF", clip_on=False)


def main():
    ap = argparse.ArgumentParser(
        description="Color differences vs lead(s) with scores embedded in sequence IDs."
    )
    ap.add_argument("-i","--in", dest="infile", default='from_excel.fa')
    ap.add_argument("--sto", default="out/alignment.sto")
    ap.add_argument("-o","--out", default="alignment_vs_lead.png")
    ap.add_argument("--wrap", type=int, default=10000)
    ap.add_argument("--dpi", type=int, default=350)
    ap.add_argument("--lead-prefix", default="lead")
    ap.add_argument("--style", choices=["group","aa"], default="group", help="Color style")
    ap.add_argument("--pale", type=float, default=0.70, help="Lighten amount for conservative swaps (0..1)")
    ap.add_argument("--show-consensus", action="store_true", help="Include consensus bar")
    ap.add_argument("--cdr-antpack-csv", dest="cdr_csv", default='', help="AntPack CDR CSV (cdr1_seq/cdr2_seq/cdr3_seq)")
    ap.add_argument("--auto-sort", choices=["y","n","ask"], default="ask", 
                    help="Auto-sort sequences: y=sort, n=keep order, ask=prompt")
    args = ap.parse_args()

    if not os.path.exists(args.infile):
        sys.exit(f"ERROR: can't find {args.infile}")

    # Check if we need to align first or if it's pre-aligned
    try:
        msa = AlignIO.read(args.infile, "fasta")
    except:
        raw = list(SeqIO.parse(args.infile, "fasta"))
        leads_exist = any(r.id.lower().startswith(args.lead_prefix.lower()) for r in raw)
        if not leads_exist and args.lead_prefix != "___no_lead___":
            sys.exit(f"ERROR: No leads found. Name at least one record starting with '{args.lead_prefix}' (no spaces).")

        # Use MAFFT for better CDR alignment (falls back to clustalo if not installed)
        out_aligned = args.sto.replace('.sto', '_aligned.fa')
        run_mafft(args.infile, out_aligned)
        msa = AlignIO.read(out_aligned, "fasta")

    # Find lead rows
    lead_rows = [i for i, rec in enumerate(msa) if rec.id.lower().startswith(args.lead_prefix.lower())]
    
    if not lead_rows and args.lead_prefix != "___no_lead___":
        ids = ", ".join(r.id for r in msa)
        sys.exit(f"ERROR: No lead rows in alignment. Got IDs: {ids}")
    elif not lead_rows:
        print(f"[INFO] No lead sequences found. Using consensus-based coloring.", file=sys.stderr)
    else:
        lead_ids = [msa[i].id for i in lead_rows]
        print(f"[INFO] Found {len(lead_rows)} lead sequence(s): {', '.join(lead_ids)}", file=sys.stderr)

    # Get reference sequence per column
    ref = lead_reference_per_column(msa, lead_rows)
    L = msa.get_alignment_length()

    # Map CDRs if CSV provided
    mapped_ranges = {}
    cdr_cols = set()
    if args.cdr_csv and os.path.exists(args.cdr_csv):
        try:
            rows = read_antpack_cdr_csv(args.cdr_csv)
            
            if rows:
                seq_idx, row_idx, mapped, hits = find_best_cdr_mapping(msa, rows, lead_rows)
                
                if hits > 0 and mapped:
                    reference_id = msa[seq_idx].id
                    print(f"[INFO] Found and mapped {len(mapped)} CDR region(s) from sequence {reference_id}", file=sys.stderr)
                    
                    mapped_ranges = split_mapped_for_wrap(mapped, args.wrap)
                    cdr_cols = get_cdr_columns_from_mapped(mapped_ranges)
                else:
                    print("[INFO] No CDR regions could be mapped. Proceeding without CDR boxes.", file=sys.stderr)
            else:
                print("[WARN] AntPack CSV had no usable rows.", file=sys.stderr)
        except Exception as e:
            print(f"[WARN] Failed to read AntPack CSV: {e}. Proceeding without CDR boxes.", file=sys.stderr)

    has_cdrs = bool(cdr_cols)
    if not has_cdrs:
        print("[INFO] No CDR regions detected - will use framework-only scoring if sorting is enabled.", file=sys.stderr)

    # Sorting logic
    should_sort = False
    weights = None
    score_display_info = []
    
    if args.auto_sort == "y":
        should_sort = True
        if has_cdrs:
            weights = {
                "cdr_dissimilar": 16.0, "cdr_similar": 3.0,
                "fr_dissimilar": 0.5, "fr_similar": 0.1, 
                "cdr_length": 100.0, "cdr3_multiplier": 1.0
            }
        else:
            weights = {
                "cdr_dissimilar": 0.0, "cdr_similar": 0.0,
                "fr_dissimilar": 10.0, "fr_similar": 2.0,
                "cdr_length": 0.0, "cdr3_multiplier": 1.0,
                "overall_length": 5.0
            }
    elif args.auto_sort == "n":
        should_sort = False
    else:  # ask
        if sys.stdin.isatty():
            should_sort, weights = prompt_for_sorting(has_cdrs)
        else:
            print("[INFO] Non-interactive mode: keeping upload order.", file=sys.stderr)
            should_sort = False

    # Apply sorting if requested
    if should_sort and weights:
        mode_desc = "CDR-aware" if has_cdrs else "framework-only"
        print(f"\n[INFO] Sorting {len(msa)} sequences using {mode_desc} similarity to {'lead' if lead_rows else 'consensus'}...", file=sys.stderr)
        msa, lead_rows, score_display_info = sort_msa_by_similarity_corrected(msa, lead_rows, ref, cdr_cols, weights, mapped_ranges)
        ref = lead_reference_per_column(msa, lead_rows)
        
        # Embed scores in sequence IDs
        msa = embed_scores_in_sequence_ids(msa, score_display_info)
        
        # Print detailed scoring info
        print("\n[DEBUG] Top 5 sequence scores:", file=sys.stderr)
        for i, info in enumerate(score_display_info[:5]):
            print(f"  {i+1}. {info['sequence_id']}: {info['score']:.1f} ({info['details']})", file=sys.stderr)

    # Set up MsaViz 
    base = {aa: "#ffffff" for aa in AA_PALETTE}
    mv = MsaViz(msa, start=1, end=None, wrap_length=args.wrap,
                show_consensus=args.show_consensus, color_scheme=None)
    mv.set_custom_color_scheme(base)

    # Color palettes
    if args.style == "aa":
        bright_color = lambda aa: AA_PALETTE.get(aa.upper())
        pale_cache = {aa: lighten_hex(col, args.pale) for aa, col in AA_PALETTE.items()}
        pale_color = lambda aa: pale_cache.get(aa.upper())
    else:  # group
        def bright_color(aa):
            g = PRIMARY_GROUP.get(aa.upper())
            return GROUP_PALETTE.get(g)
        pale_cache = {g: lighten_hex(col, args.pale) for g, col in GROUP_PALETTE.items()}
        def pale_color(aa):
            g = PRIMARY_GROUP.get(aa.upper())
            return pale_cache.get(g)

    # Find columns with differences for markers
    diff_cols = []
    for j in range(L):
        r = ref[j]
        anydiff = False
        for i in range(len(msa)):
            if i in lead_rows: continue
            aa = msa[i].seq[j].upper()
            if aa in "-.*?": continue
            if r is None or aa != r:
                anydiff = True
                break
        if anydiff:
            diff_cols.append(j+1)

    # Per-cell coloring function
    def color_func(row_pos, col_pos, aa, _msa):
        aa = aa.upper()
        if row_pos in lead_rows or aa in "-.*?": 
            return None
        r = ref[col_pos]
        if r is None: 
            return bright_color(aa)
        if aa == r:   
            return None
        return pale_color(aa) if strict_similarity(aa, r) else bright_color(aa)

    mv.set_custom_color_func(color_func)
    if diff_cols: 
        mv.add_markers(diff_cols, marker="o")

    # Render plot
    fig = mv.plotfig()
    ax = fig.axes[0]

    # IMGT numbering + hallmark/vernier markers under consensus
    annotate_imgt_and_markers(ax, msa, lead_rows, args.wrap)

    # Draw CDR boxes (only if we have them)
    if mapped_ranges:
        draw_cdr_boxes(ax, mapped_ranges, fill="#7b777d", alpha=0.25, z=8)
        print(f"[INFO] Drew CDR boxes for {len(mapped_ranges)} regions.", file=sys.stderr)

    # Bold lead rows
    bolden_lead_rows(fig, msa, lead_rows)

    # add the legend strip on top and save
    compose_with_top_legend(fig_main=fig)
    fig.savefig(args.out, dpi=args.dpi)
    plt.close(fig)
    print(f"Done → {args.out}")
    
    print(f"\n[SUCCESS] Saved alignment visualization: {args.out}", file=sys.stderr)
    
    if should_sort:
        mode_desc = "CDR-aware" if has_cdrs else "framework-only"
        print(f"[INFO] Sequences were sorted using {mode_desc} weights:", file=sys.stderr)
        for key, val in weights.items():
            if val > 0:  # Only show weights that were used
                print(f"  {key.replace('_', ' ').title()}: {val}", file=sys.stderr)
        
        # Prompt for export
        if sys.stdin.isatty() and score_display_info:
            export_response = input("\nExport sequence ranking to CSV? (y/n, default=n): ").strip().lower()
            if export_response in ['y', 'yes']:
                export_path = args.out.replace('.png', '_ranking.csv')
                export_sorting_results(msa, score_display_info, export_path, weights)

if __name__ == "__main__":
    main()