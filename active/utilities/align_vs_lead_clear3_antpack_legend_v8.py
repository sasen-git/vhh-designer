#!/usr/bin/env python3
"""
align_vs_lead_clear3_antpack_legend_v8.py

v8 Changes:
- Added IMGT position track below consensus
- Shows H (hallmark) and V (vernier) markers at key VHH positions
- IMGT numbering derived from AntPack

v7 Changes:
- Switched from Clustal Omega to MAFFT --localpair for better CDR alignment
- MAFFT's local pairwise alignment keeps motifs together (e.g., DYD aligns with DYD)
- Falls back to Clustal Omega if MAFFT is not installed
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

# --- IMGT Key Positions for VHH ---
IMGT_HALLMARKS = {42, 49, 50, 52}
IMGT_VERNIERS = {66, 67, 68, 69, 71, 76, 78, 80, 82, 83, 84, 85, 87, 89, 91, 94, 108}

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

def get_imgt_numbering_for_sequence(sequence, scheme="imgt"):
    """
    Use AntPack to get IMGT numbering for a sequence.
    Returns dict mapping 0-based sequence position to IMGT number.
    """
    try:
        from antpack import SingleChainAnnotator
    except ImportError:
        print("[WARN] AntPack not available for IMGT numbering", file=sys.stderr)
        return {}
    
    annot = SingleChainAnnotator(chains=["H"], scheme=scheme)
    try:
        numbering, pid, chain, err = annot.analyze_seq(str(sequence))
        if err or not numbering:
            print(f"[WARN] AntPack returned error or empty numbering: {err}", file=sys.stderr)
            return {}
        
        # Debug: show first few elements to understand format
        if numbering:
            print(f"[DEBUG] AntPack numbering format sample: {numbering[:3]} (type: {type(numbering[0])})", file=sys.stderr)
        
        # numbering format can vary - handle different cases
        pos_to_imgt = {}
        for i, item in enumerate(numbering):
            try:
                if item is None:
                    continue
                
                imgt_num = None
                
                # Tuple/list format: (position, insertion_code) or just (position,)
                if isinstance(item, (tuple, list)):
                    if len(item) >= 1 and item[0] is not None:
                        imgt_num = item[0]
                
                # String format: "1", "27", "27A", etc.
                elif isinstance(item, str):
                    num_str = ''.join(c for c in item if c.isdigit())
                    if num_str:
                        imgt_num = int(num_str)
                
                # Integer/float format
                elif isinstance(item, (int, float)):
                    imgt_num = int(item)
                
                # Store the mapping
                if imgt_num is not None:
                    pos_to_imgt[i] = int(imgt_num)
                    
            except (ValueError, TypeError, IndexError) as e:
                continue
        
        print(f"[INFO] IMGT mapping: {len(pos_to_imgt)} positions mapped", file=sys.stderr)
        return pos_to_imgt
        
    except Exception as e:
        print(f"[WARN] IMGT numbering failed: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return {}

def get_alignment_col_to_imgt(aligned_seq, pos_to_imgt):
    """
    Map alignment columns (1-based) to IMGT positions.
    aligned_seq: string with gaps
    pos_to_imgt: dict from 0-based ungapped position to IMGT number
    Returns dict mapping 1-based alignment column to IMGT number.
    """
    col_to_imgt = {}
    ungapped_pos = 0
    
    for col, aa in enumerate(str(aligned_seq)):
        if aa not in "-.*?":
            if ungapped_pos in pos_to_imgt:
                col_to_imgt[col + 1] = pos_to_imgt[ungapped_pos]
            ungapped_pos += 1
    
    return col_to_imgt

def add_imgt_track(fig, ax, col_to_imgt, wrap_length, total_cols):
    """
    Add IMGT position markers by highlighting hallmark and vernier columns
    with subtle background colors, and add IMGT number labels.
    """
    if not col_to_imgt:
        return
    
    # Get axes limits
    y_min, y_max = ax.get_ylim()
    
    # Collect positions for hallmarks and verniers
    hallmark_cols = []
    vernier_cols = []
    
    for col, imgt_pos in col_to_imgt.items():
        if imgt_pos in IMGT_HALLMARKS:
            hallmark_cols.append((col, imgt_pos))
        elif imgt_pos in IMGT_VERNIERS:
            vernier_cols.append((col, imgt_pos))
    
    # Add vertical spans for hallmark positions (red, more visible)
    for col, imgt_pos in hallmark_cols:
        ax.axvspan(col - 0.5, col + 0.5, alpha=0.15, color='#E41A1C', zorder=0)
    
    # Add vertical spans for vernier positions (blue, subtle)
    for col, imgt_pos in vernier_cols:
        ax.axvspan(col - 0.5, col + 0.5, alpha=0.08, color='#377EB8', zorder=0)
    
    # Add IMGT number labels at key positions (below the alignment)
    # y_max is at the bottom in inverted coordinates
    label_y = y_max + 0.8
    
    # Label hallmarks with IMGT numbers (red, bold)
    for col, imgt_pos in hallmark_cols:
        ax.text(col, label_y, str(imgt_pos), 
                ha='center', va='top', fontsize=5, fontweight='bold',
                color='#E41A1C', rotation=0)
    
    # Label verniers with IMGT numbers (blue, smaller)
    for col, imgt_pos in vernier_cols:
        ax.text(col, label_y + 0.6, str(imgt_pos),
                ha='center', va='top', fontsize=4,
                color='#377EB8', rotation=90)
    
    # Print summary
    print(f"[INFO] Highlighted {len(hallmark_cols)} hallmark (red) and {len(vernier_cols)} vernier (blue) columns", file=sys.stderr)

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
    """
    exe = shutil.which("mafft")
    if not exe:
        print("[WARN] mafft not found, falling back to clustalo.", file=sys.stderr)
        out_sto = out_fa.replace('.fa', '.sto')
        run_clustalo(in_fa, out_sto)
        msa = AlignIO.read(out_sto, "stockholm")
        AlignIO.write(msa, out_fa, "fasta")
        return
    
    cmd = [exe, "--localpair", "--maxiterate", "1000", "--quiet", in_fa]
    print(f"[INFO] Running MAFFT with --localpair...", file=sys.stderr)
    r = subprocess.run(cmd, capture_output=True, text=True)
    
    if r.returncode != 0:
        print(f"[WARN] MAFFT failed, falling back to clustalo...", file=sys.stderr)
        out_sto = out_fa.replace('.fa', '.sto')
        run_clustalo(in_fa, out_sto)
        msa = AlignIO.read(out_sto, "stockholm")
        AlignIO.write(msa, out_fa, "fasta")
        return
    
    with open(out_fa, 'w') as f:
        f.write(r.stdout)
    
    if not os.path.exists(out_fa) or os.path.getsize(out_fa) == 0:
        sys.exit("ERROR: empty alignment from MAFFT.")

def lead_reference_per_column(msa, lead_rows):
    L = msa.get_alignment_length()
    ref = []
    
    if not lead_rows:
        for j in range(L):
            col = [msa[i].seq[j] for i in range(len(msa))]
            col = [aa for aa in col if aa not in "-.*?"]
            ref.append(Counter(col).most_common(1)[0][0] if col else None)
        return ref
    
    for j in range(L):
        col = [msa[r].seq[j] for r in lead_rows]
        col = [aa for aa in col if aa not in "-.*?"]
        ref.append(Counter(col).most_common(1)[0][0] if col else None)
    return ref

def strict_similarity(a: str, b: str) -> bool:
    ga = PRIMARY_GROUP.get(a.upper()); gb = PRIMARY_GROUP.get(b.upper())
    return ga is not None and gb is not None and ga == gb

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
    import matplotlib.colors as mcolors
    ax1 = fig_main.get_axes()[0]

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
    
    reference_height = 5.0
    fig_height = fig_main.get_figheight()
    scale_factor = reference_height / fig_height
    
    legend_gap = 0.15 * scale_factor
    legend_height = 0.15 * scale_factor
    title_gap = 0.25 * scale_factor
    
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

    fig_main.suptitle("Groups — bright=dissimilar, pale=similar  |  H=hallmark V=vernier", 
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
    
    for j in range(min(L, len(ref_chars))):
        aa = seq_str[j].upper()
        ref_aa = ref_chars[j]
        
        if aa in "-.*?" or ref_aa is None or aa == ref_aa:
            continue
            
        mutation_positions.append(j + 1)
        is_similar = similar_group(aa, ref_aa)
        
        if has_cdrs:
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
            if is_similar:
                fr_similar += 1
            else:
                fr_dissimilar += 1
    
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
    
    overall_length_penalty = 0
    if ref_record:
        seq_gaps = sum(1 for ch in seq_str if ch in "-.*?")
        ref_gaps = sum(1 for ch in ref_str if ch in "-.*?")
        overall_length_penalty = abs(seq_gaps - ref_gaps)
    
    if has_cdrs:
        score = (weights["cdr_dissimilar"] * cdr_dissimilar +
                 weights["cdr_similar"] * cdr_similar +
                 weights["cdr_dissimilar"] * cdr3_dissimilar * cdr3_multiplier +
                 weights["cdr_similar"] * cdr3_similar * cdr3_multiplier +
                 weights["fr_dissimilar"] * fr_dissimilar +
                 weights["fr_similar"] * fr_similar +
                 weights["cdr_length"] * cdr_length_penalty)
    else:
        score = (weights["fr_dissimilar"] * fr_dissimilar +
                 weights["fr_similar"] * fr_similar +
                 weights.get("overall_length", 5.0) * overall_length_penalty)
    
    total_mutations = len(mutation_positions)
    earliest_mutation = min(mutation_positions) if mutation_positions else float('inf')
    
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
    ref_record = msa[lead_rows[0]] if lead_rows else msa[0]
    
    has_cdrs = bool(cdr_cols)
    
    cdr3_cols = None
    if has_cdrs and weights.get("cdr3_multiplier", 1.0) > 1.0 and mapped_ranges:
        cdr3_cols = get_cdr3_columns(mapped_ranges)
        if cdr3_cols:
            print(f"[INFO] Applying {weights['cdr3_multiplier']}x weight to {len(cdr3_cols)} CDR3 positions", file=sys.stderr)
    
    if not has_cdrs:
        print("[INFO] No CDRs detected - using framework-only similarity scoring", file=sys.stderr)
    
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
    
    def sort_key(item):
        i, record, score_info, is_lead = item
        return (
            score_info['score'],
            score_info['total_mutations'],
            score_info['earliest_mutation'],
            record.id.lower()
        )
    
    seq_scores.sort(key=sort_key)
    
    sorted_records = [item[1] for item in seq_scores]
    sorted_msa = MultipleSeqAlignment(sorted_records)
    
    new_lead_rows = []
    for new_idx, (old_idx, record, score_info, is_lead) in enumerate(seq_scores):
        if is_lead:
            new_lead_rows.append(new_idx)
    
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
    cdr3_cols = set()
    for name, (a, b) in mapped_ranges.items():
        if "CDR3" in name.upper():
            for col in range(a, b + 1):
                cdr3_cols.add(col)
    return cdr3_cols

def prompt_for_sorting(has_cdrs=True):
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
    
    while True:
        choice = input("Sort sequences by similarity? (y/n): ").strip().lower()
        if choice in ['y', 'yes']:
            return True, prompt_for_weights(has_cdrs)
        elif choice in ['n', 'no']:
            return False, None
        else:
            print("Please enter 'y' or 'n'")

def prompt_for_weights(has_cdrs=True):
    if has_cdrs:
        default_weights = {
            "cdr_dissimilar": 16.0, "cdr_similar": 3.0,
            "fr_dissimilar": 0.5, "fr_similar": 0.1,
            "cdr_length": 100.0, "cdr3_multiplier": 1.0
        }
    else:
        default_weights = {
            "cdr_dissimilar": 0.0, "cdr_similar": 0.0,
            "fr_dissimilar": 10.0, "fr_similar": 2.0,
            "cdr_length": 0.0, "cdr3_multiplier": 1.0,
            "overall_length": 5.0
        }
    
    print("\nUsing default weights. Customize? (y/n, default=n): ", end="")
    use_custom = input().strip().lower()
    
    if use_custom in ['y', 'yes']:
        print("\nEnter custom weights (press Enter to use default):")
        weights = {}
        for key, default_val in default_weights.items():
            while True:
                user_input = input(f"  {key} [{default_val}]: ").strip()
                if not user_input:
                    weights[key] = default_val
                    break
                try:
                    weights[key] = float(user_input)
                    break
                except ValueError:
                    print("    Please enter a valid number")
        return weights
    else:
        return default_weights

def main():
    ap = argparse.ArgumentParser(
        description="Color differences vs lead(s) with IMGT position markers."
    )
    ap.add_argument("-i","--in", dest="infile", default='from_excel.fa')
    ap.add_argument("--sto", default="out/alignment.sto")
    ap.add_argument("-o","--out", default="alignment_vs_lead.png")
    ap.add_argument("--wrap", type=int, default=10000)
    ap.add_argument("--dpi", type=int, default=350)
    ap.add_argument("--lead-prefix", default="lead")
    ap.add_argument("--style", choices=["group","aa"], default="group")
    ap.add_argument("--pale", type=float, default=0.70)
    ap.add_argument("--show-consensus", action="store_true")
    ap.add_argument("--show-imgt", action="store_true", default=True, help="Show IMGT H/V markers")
    ap.add_argument("--cdr-antpack-csv", dest="cdr_csv", default='')
    ap.add_argument("--auto-sort", choices=["y","n","ask"], default="ask")
    args = ap.parse_args()

    if not os.path.exists(args.infile):
        sys.exit(f"ERROR: can't find {args.infile}")

    try:
        msa = AlignIO.read(args.infile, "fasta")
    except:
        raw = list(SeqIO.parse(args.infile, "fasta"))
        leads_exist = any(r.id.lower().startswith(args.lead_prefix.lower()) for r in raw)
        if not leads_exist and args.lead_prefix != "___no_lead___":
            sys.exit(f"ERROR: No leads found.")

        out_aligned = args.sto.replace('.sto', '_aligned.fa')
        run_mafft(args.infile, out_aligned)
        msa = AlignIO.read(out_aligned, "fasta")

    lead_rows = [i for i, rec in enumerate(msa) if rec.id.lower().startswith(args.lead_prefix.lower())]
    
    if not lead_rows and args.lead_prefix != "___no_lead___":
        sys.exit(f"ERROR: No lead rows in alignment.")
    elif not lead_rows:
        print(f"[INFO] No lead sequences found. Using consensus-based coloring.", file=sys.stderr)

    ref = lead_reference_per_column(msa, lead_rows)
    L = msa.get_alignment_length()

    # Get IMGT numbering for lead/first sequence
    col_to_imgt = {}
    if args.show_imgt:
        ref_seq_idx = lead_rows[0] if lead_rows else 0
        ref_seq = str(msa[ref_seq_idx].seq).replace("-", "")
        pos_to_imgt = get_imgt_numbering_for_sequence(ref_seq)
        if pos_to_imgt:
            col_to_imgt = get_alignment_col_to_imgt(msa[ref_seq_idx].seq, pos_to_imgt)
            n_hallmarks = sum(1 for c, i in col_to_imgt.items() if i in IMGT_HALLMARKS)
            n_verniers = sum(1 for c, i in col_to_imgt.items() if i in IMGT_VERNIERS)
            print(f"[INFO] IMGT mapping: {n_hallmarks} hallmark positions, {n_verniers} vernier positions", file=sys.stderr)

    # Map CDRs
    mapped_ranges = {}
    cdr_cols = set()
    if args.cdr_csv and os.path.exists(args.cdr_csv):
        try:
            rows = read_antpack_cdr_csv(args.cdr_csv)
            if rows:
                seq_idx, row_idx, mapped, hits = find_best_cdr_mapping(msa, rows, lead_rows)
                if hits > 0 and mapped:
                    mapped_ranges = split_mapped_for_wrap(mapped, args.wrap)
                    cdr_cols = get_cdr_columns_from_mapped(mapped_ranges)
        except Exception as e:
            print(f"[WARN] Failed to read CDR CSV: {e}", file=sys.stderr)

    has_cdrs = bool(cdr_cols)

    # Sorting
    should_sort = False
    weights = None
    score_display_info = []
    
    if args.auto_sort == "y":
        should_sort = True
        weights = {"cdr_dissimilar": 16.0, "cdr_similar": 3.0,
                   "fr_dissimilar": 0.5, "fr_similar": 0.1, 
                   "cdr_length": 100.0, "cdr3_multiplier": 1.0} if has_cdrs else {
                   "cdr_dissimilar": 0.0, "cdr_similar": 0.0,
                   "fr_dissimilar": 10.0, "fr_similar": 2.0,
                   "cdr_length": 0.0, "cdr3_multiplier": 1.0, "overall_length": 5.0}
    elif args.auto_sort == "n":
        should_sort = False
    else:
        if sys.stdin.isatty():
            should_sort, weights = prompt_for_sorting(has_cdrs)
        else:
            should_sort = False

    if should_sort and weights:
        msa, lead_rows, score_display_info = sort_msa_by_similarity_corrected(msa, lead_rows, ref, cdr_cols, weights, mapped_ranges)
        ref = lead_reference_per_column(msa, lead_rows)
        msa = embed_scores_in_sequence_ids(msa, score_display_info)
        
        # Re-map IMGT columns after sorting
        if args.show_imgt and pos_to_imgt:
            ref_seq_idx = lead_rows[0] if lead_rows else 0
            col_to_imgt = get_alignment_col_to_imgt(msa[ref_seq_idx].seq, pos_to_imgt)

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
    else:
        def bright_color(aa):
            g = PRIMARY_GROUP.get(aa.upper())
            return GROUP_PALETTE.get(g)
        pale_cache = {g: lighten_hex(col, args.pale) for g, col in GROUP_PALETTE.items()}
        def pale_color(aa):
            g = PRIMARY_GROUP.get(aa.upper())
            return pale_cache.get(g)

    # Find difference columns
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

    # Render
    fig = mv.plotfig()
    ax = fig.axes[0]

    # Draw CDR boxes
    if mapped_ranges:
        draw_cdr_boxes(ax, mapped_ranges, fill="#7b777d", alpha=0.25, z=8)

    # Add IMGT H/V markers below the alignment
    if col_to_imgt and args.show_imgt:
        add_imgt_track(fig, ax, col_to_imgt, args.wrap, L)

    bolden_lead_rows(fig, msa, lead_rows)
    compose_with_top_legend(fig_main=fig)
    
    fig.savefig(args.out, dpi=args.dpi, bbox_inches='tight')
    plt.close(fig)
    print(f"Done → {args.out}")

if __name__ == "__main__":
    main()
