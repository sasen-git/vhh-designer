#!/usr/bin/env python3
"""
csv2msa_antpack_v3.py

Updated wrapper that handles sorting even when no CDRs are detected.
Falls back to framework-only similarity scoring.
"""

import os, sys, re, csv, argparse, subprocess, shutil

os.environ.setdefault("MPLBACKEND", "Agg")

# Tag definitions
KNOWN_TAGS = {
    "MSKIK": {"name": "MSKIK", "position": "N-term"},
    "YPYDVPDYAYPYDVPDYAYPYDVPDYA": {"name": "3xHA", "position": "both"},
    "WSHPQFEK": {"name": "Strep-tag", "position": "C-term"},
    "EPEA": {"name": "C-tag", "position": "C-term"},
    "HHHHHH": {"name": "6xHis-tag", "position": "both"},
    "HHHHHHHH": {"name": "8xHis-tag", "position": "both"},
}

def detect_tags_in_sequence(seq):
    """Detect all known tags in a sequence."""
    found_tags = {}
    seq_upper = seq.upper()
    
    for tag_seq, tag_info in KNOWN_TAGS.items():
        positions = []
        tag_name = tag_info["name"]
        
        start = 0
        while True:
            idx = seq_upper.find(tag_seq, start)
            if idx == -1:
                break
            
            at_nterm = (idx < 20)
            at_cterm = (idx > len(seq) - len(tag_seq) - 20)
            
            if tag_info["position"] == "N-term" and at_nterm:
                positions.append(("N-term", idx, idx + len(tag_seq)))
            elif tag_info["position"] == "C-term" and at_cterm:
                positions.append(("C-term", idx, idx + len(tag_seq)))
            elif tag_info["position"] == "both":
                if at_nterm:
                    positions.append(("N-term", idx, idx + len(tag_seq)))
                elif at_cterm:
                    positions.append(("C-term", idx, idx + len(tag_seq)))
            
            start = idx + 1
        
        if positions:
            found_tags[tag_name] = positions
    
    return found_tags

def analyze_tags_in_dataset(rows):
    """Analyze all sequences for tags."""
    tag_analysis = {}
    
    for is_lead, rid, seq in rows:
        found = detect_tags_in_sequence(seq)
        for tag_name, positions in found.items():
            if tag_name not in tag_analysis:
                tag_analysis[tag_name] = []
            tag_analysis[tag_name].append((rid, positions))
    
    return tag_analysis

def remove_tags_from_sequence(seq, tags_to_remove):
    """Remove specified tags from sequence."""
    if not tags_to_remove:
        return seq, []
    
    all_removals = []
    for tag_name, positions in tags_to_remove.items():
        for pos_type, start, end in positions:
            all_removals.append((start, end, tag_name, pos_type))
    
    all_removals.sort(reverse=True)
    
    cleaned_seq = seq
    removed_report = []
    
    for start, end, tag_name, pos_type in all_removals:
        removed_seq = cleaned_seq[start:end]
        cleaned_seq = cleaned_seq[:start] + cleaned_seq[end:]
        removed_report.append(f"{tag_name} at {pos_type}")
    
    return cleaned_seq, removed_report

def prompt_for_tag_removal(tag_analysis, total_sequences):
    """Ask user which tags to remove."""
    if not tag_analysis:
        print("\n[INFO] No common affinity tags detected in sequences.")
        return set()
    
    print("\n" + "="*80)
    print("AFFINITY TAG DETECTION")
    print("="*80)
    print("\nThe following affinity tags were detected in your sequences:\n")
    
    for tag_name, seq_list in sorted(tag_analysis.items()):
        count = len(seq_list)
        percentage = (count / total_sequences) * 100
        print(f"  • {tag_name}: {count}/{total_sequences} sequences ({percentage:.0f}%)")
    
    print("\nWould you like to remove any of these tags before alignment?")
    print("(Removing tags can improve alignment quality for the functional domains)")
    print()
    
    tags_to_remove = set()
    for tag_name, seq_list in sorted(tag_analysis.items()):
        count = len(seq_list)
        response = input(f"Remove {tag_name}? (y/n, default=y): ").strip().lower()
        if response in ['y', 'yes', '']:
            tags_to_remove.add(tag_name)
    
    return tags_to_remove

def apply_tag_removal(rows, tag_analysis, tags_to_remove):
    """Remove selected tags from all sequences and report results."""
    if not tags_to_remove:
        return rows
    
    print("\n" + "="*80)
    print("TAG REMOVAL REPORT")
    print("="*80)
    
    seq_tags = {}
    for tag_name in tags_to_remove:
        if tag_name in tag_analysis:
            for rid, positions in tag_analysis[tag_name]:
                if rid not in seq_tags:
                    seq_tags[rid] = {}
                seq_tags[rid][tag_name] = positions
    
    cleaned_rows = []
    removal_stats = {tag: {"removed": 0, "not_found": []} for tag in tags_to_remove}
    
    for is_lead, rid, seq in rows:
        if rid in seq_tags:
            cleaned_seq, removed_report = remove_tags_from_sequence(seq, seq_tags[rid])
            cleaned_rows.append((is_lead, rid, cleaned_seq))
            
            for tag_name in tags_to_remove:
                if tag_name in seq_tags[rid]:
                    removal_stats[tag_name]["removed"] += 1
                else:
                    removal_stats[tag_name]["not_found"].append(rid)
        else:
            cleaned_rows.append((is_lead, rid, seq))
            for tag_name in tags_to_remove:
                removal_stats[tag_name]["not_found"].append(rid)
    
    total = len(rows)
    for tag_name in sorted(tags_to_remove):
        stats = removal_stats[tag_name]
        removed_count = stats["removed"]
        not_found = stats["not_found"]
        
        print(f"\n{tag_name}:")
        print(f"  ✓ Removed from {removed_count}/{total} sequences")
        
        if not_found:
            if len(not_found) <= 5:
                print(f"  ⚠ Not found in: {', '.join(not_found)}")
            else:
                print(f"  ⚠ Not found in {len(not_found)} sequences: {', '.join(not_found[:3])}, ...")
    
    print("\n" + "="*80)
    return cleaned_rows

def sanitize_id(s: str) -> str:
    s = (s or "").strip()
    s = re.sub(r"\s+", "_", s)
    return re.sub(r"[^A-Za-z0-9_.:+-]", "_", s)

def truthy(x) -> bool:
    return str(x).strip().lower() in {"1","t","true","y","yes"}

def read_threecol_csv(path):
    with open(path, newline="") as f:
        r = csv.DictReader(f)
        cols = {h.lower().strip(): h for h in r.fieldnames}

        def pick(*cands):
            for c in cands:
                if c in cols: return cols[c]
            want = {re.sub(r"[^a-z]","",c) for c in cands}
            for h in r.fieldnames:
                if re.sub(r"[^a-z]","",h.lower()) in want:
                    return h
            return None

        h_lead = pick("lead","is_lead","lead?","lead_tf")
        h_id   = pick("seq id","sequence id","id","name","seqid")
        h_seq  = pick("aa sequence","sequence","aa","protein sequence","seq")
        if not (h_lead and h_id and h_seq):
            sys.exit(f"ERROR: CSV must have columns like [lead, Seq ID, AA sequence]. Found: {r.fieldnames}")

        rows=[]
        for row in r:
            is_lead = truthy(row[h_lead])
            rid     = sanitize_id(row[h_id])
            seq     = re.sub(r"[^A-Za-z]","", str(row[h_seq]).upper())
            if not rid or not seq:
                continue
            rows.append((is_lead, rid, seq))
        return rows

def write_fasta(path, records):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w") as f:
        for rid, seq in records:
            f.write(f">{rid}\n{seq}\n")

def antpack_cdrs_for_leads(leads, scheme="imgt", chains="H"):
    """
    Use AntPack to analyze CDR regions for lead sequences.
    Returns empty CDRs if analysis fails rather than crashing.
    """
    try:
        from antpack import SingleChainAnnotator
    except ImportError:
        print("[WARN] AntPack not available. CDR analysis will be skipped.", file=sys.stderr)
        return [{"id": rid, "cdr1_seq": "", "cdr2_seq": "", "cdr3_seq": ""} for rid, seq in leads]
    
    annot = SingleChainAnnotator(chains=[c.strip() for c in chains.split(",") if c.strip()],
                                 scheme=scheme)
    out_rows=[]
    
    for rid, seq in leads:
        try:
            seq_str = str(seq) if seq is not None else ""
            if not seq_str:
                print(f"[WARN] Empty sequence for {rid}", file=sys.stderr)
                out_rows.append({
                    "id": rid,
                    "cdr1_seq": "",
                    "cdr2_seq": "",
                    "cdr3_seq": "",
                })
                continue
                
            numbering, pid, chain, err = annot.analyze_seq(seq_str)
            
            if err or not numbering:
                print(f"[WARN] AntPack could not analyze sequence {rid}: {err or 'Unknown error'}", file=sys.stderr)
                out_rows.append({
                    "id": rid,
                    "cdr1_seq": "",
                    "cdr2_seq": "",
                    "cdr3_seq": "",
                })
                continue
            
            trimmed_seq, trimmed_nums, _, _ = annot.trim_alignment(seq_str, (numbering, pid, chain, err))
            labels = annot.assign_cdr_labels(trimmed_nums, chain)

            def spans(nums, labs):
                ans={}; cur=None; start=None
                for i,(n,lb) in enumerate(zip(nums,labs)):
                    if lb in ("cdr1","cdr2","cdr3"):
                        if cur is None: cur, start = lb, i
                    else:
                        if cur is not None:
                            ans[cur]=(start,i-1); cur=None
                if cur is not None:
                    ans[cur]=(start,len(labs)-1)
                return ans

            s = spans(trimmed_nums, labels)

            def extract(tag):
                if tag in s:
                    a,b = s[tag]
                    return trimmed_seq[a:b+1]
                return ""

            out_rows.append({
                "id": rid,
                "cdr1_seq": extract("cdr1"),
                "cdr2_seq": extract("cdr2"),
                "cdr3_seq": extract("cdr3"),
            })
            
            cdrs_found = [tag for tag in ["cdr1", "cdr2", "cdr3"] if extract(tag)]
            if cdrs_found:
                print(f"[INFO] AntPack found CDRs for {rid}: {', '.join([c.upper() for c in cdrs_found])}", file=sys.stderr)
            else:
                print(f"[WARN] AntPack analyzed {rid} but found no CDR regions", file=sys.stderr)
                
        except Exception as e:
            print(f"[WARN] AntPack failed to process sequence {rid}: {e}", file=sys.stderr)
            out_rows.append({
                "id": rid,
                "cdr1_seq": "",
                "cdr2_seq": "",
                "cdr3_seq": "",
            })
            
    return out_rows
    
def prompt_for_manual_cdr_input(lead_records):
    """
    Prompt user to manually enter CDR sequences if AntPack failed.
    """
    print("\n" + "="*80)
    print("MANUAL CDR INPUT")
    print("="*80)
    print("\nAntPack was unable to automatically detect CDR regions.")
    print("Would you like to manually enter CDR sequences?")
    print("(This will enable CDR-aware similarity scoring)")
    print()
    
    response = input("Enter CDR sequences manually? (y/n): ").strip().lower()
    if response not in ['y', 'yes']:
        return []
    
    manual_cdr_rows = []
    
    for rid, seq in lead_records:
        print(f"\n--- Enter CDRs for: {rid} ---")
        print(f"Full sequence length: {len(seq)} aa")
        print()
        
        cdr1 = input("CDR1 sequence (or press Enter to skip): ").strip().upper()
        cdr1 = re.sub(r"[^A-Z]", "", cdr1)
        
        cdr2 = input("CDR2 sequence (or press Enter to skip): ").strip().upper()
        cdr2 = re.sub(r"[^A-Z]", "", cdr2)
        
        cdr3 = input("CDR3 sequence (or press Enter to skip): ").strip().upper()
        cdr3 = re.sub(r"[^A-Z]", "", cdr3)
        
        # Validate that CDRs exist in sequence
        validated_cdrs = {"id": rid, "cdr1_seq": "", "cdr2_seq": "", "cdr3_seq": ""}
        
        if cdr1 and cdr1 in seq:
            validated_cdrs["cdr1_seq"] = cdr1
            print(f"  ✓ CDR1 found in sequence")
        elif cdr1:
            print(f"  ✗ CDR1 not found in sequence - skipping")
        
        if cdr2 and cdr2 in seq:
            validated_cdrs["cdr2_seq"] = cdr2
            print(f"  ✓ CDR2 found in sequence")
        elif cdr2:
            print(f"  ✗ CDR2 not found in sequence - skipping")
        
        if cdr3 and cdr3 in seq:
            validated_cdrs["cdr3_seq"] = cdr3
            print(f"  ✓ CDR3 found in sequence")
        elif cdr3:
            print(f"  ✗ CDR3 not found in sequence - skipping")
        
        manual_cdr_rows.append(validated_cdrs)
    
    return manual_cdr_rows

def write_cdr_csv(path, rows):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["id","cdr1_seq","cdr2_seq","cdr3_seq"])
        w.writeheader()
        for r in rows: w.writerow(r)

def main():
    ap = argparse.ArgumentParser(
        description="CSV → FASTA, AntPack CDRs, Clustal Omega MSA, and colored figure with optional sorting (works even without CDRs)"
    )
    ap.add_argument("--table", required=True, help="Input CSV with columns: lead, Seq ID, AA sequence")
    ap.add_argument("--scheme", default="imgt", choices=["imgt","kabat","martin","aho"],
                    help="AntPack numbering scheme (default: imgt)")
    ap.add_argument("--auto-sort", choices=["y","n","ask"], default="ask", 
                    help="Auto-sort by similarity: y=sort, n=keep order, ask=prompt")
    args = ap.parse_args()

    csv_path = args.table
    if not os.path.exists(csv_path):
        sys.exit(f"ERROR: CSV not found: {csv_path}")

    # output folder: <csvbase>_<scriptbase>/
    script_base = os.path.splitext(os.path.basename(sys.argv[0]))[0]
    csv_dir  = os.path.dirname(csv_path) or "."
    csv_base = os.path.splitext(os.path.basename(csv_path))[0]
    out_dir  = os.path.join(csv_dir, f"{csv_base}_{script_base}")
    os.makedirs(out_dir, exist_ok=True)

    # file stems inside the folder
    base     = os.path.join(out_dir, f"{csv_base}_{script_base}")
    fa_path  = f"{base}.fa"
    sto_path = f"{base}.sto"
    cdr_path = f"{base}_cdrs.csv"
    png_path = f"{base}.png"

    # --- read CSV
    rows = read_threecol_csv(csv_path)

    # --- Tag detection and removal
    tag_analysis = analyze_tags_in_dataset(rows)
    if tag_analysis and sys.stdin.isatty():
        tags_to_remove = prompt_for_tag_removal(tag_analysis, len(rows))
        if tags_to_remove:
            rows = apply_tag_removal(rows, tag_analysis, tags_to_remove)
    elif tag_analysis:
        print(f"\n[INFO] Detected tags in dataset but running in non-interactive mode. Keeping tags.", file=sys.stderr)
        for tag_name, seq_list in sorted(tag_analysis.items()):
            count = len(seq_list)
            print(f"  • {tag_name}: {count}/{len(rows)} sequences", file=sys.stderr)
    
    # Check for long sequences and ask about wrapping
    all_seqs = [(rid, seq) for _, rid, seq in rows]
    long_sequences = [(rid, len(seq)) for rid, seq in all_seqs if len(seq) > 200]
    wrap_length = "10000"  # default no wrap
    
    if long_sequences:
        max_len = max(len(seq) for _, seq in all_seqs)
        print(f"\n[INFO] Found {len(long_sequences)} sequence(s) > 200 amino acids (longest: {max_len} aa)")
        if sys.stdin.isatty():
            response = input("Would you like to wrap the alignment? (y/n): ").strip().lower()
            if response in ['y', 'yes']:
                wrap_input = input(f"Enter wrap length (press Enter for 80): ").strip()
                wrap_length = wrap_input if wrap_input else "80"

    # Separate leads and non-leads
    leads = [(rid, seq) for is_lead, rid, seq in rows if is_lead]
    non_leads = [(rid, seq) for is_lead, rid, seq in rows if not is_lead]
    
    # decide lead_prefix for renderer
    if not leads:
        print("[INFO] No leads in CSV → consensus mode", file=sys.stderr)
        lead_prefix = "___no_lead___"
        records = all_seqs
    else:
        # Rename leads to have a unique prefix to avoid ID conflicts
        lead_records = [(f"lead_{rid}", seq) for rid, seq in leads]
        records = lead_records + non_leads
        lead_prefix = "lead_"
        
    # --- write FASTA
    write_fasta(fa_path, records)
    print(f"✓ FASTA → {fa_path}")

    # --- AntPack CDRs (always try, but handle failures gracefully)
    cdr_rows = []
    if leads:
        # Rename leads for CDR analysis to match FASTA
        leads_for_antpack = [(f"lead_{rid}", seq) for rid, seq in leads]
        cdr_rows = antpack_cdrs_for_leads(leads_for_antpack, scheme=args.scheme, chains="H")
    elif records:
        # No leads - try first sequence
        print("[INFO] No leads found. Attempting CDR analysis on first sequence...", file=sys.stderr)
        cdr_rows = antpack_cdrs_for_leads([records[0]], scheme=args.scheme, chains="H")

    # Check if any CDRs were actually found
    has_any_cdrs = False
    if cdr_rows:
        for row in cdr_rows:
            if row.get("cdr1_seq") or row.get("cdr2_seq") or row.get("cdr3_seq"):
                has_any_cdrs = True
                break

    # NEW: If no CDRs found and interactive, offer manual input
    if not has_any_cdrs and sys.stdin.isatty():
        print("[INFO] No CDR regions detected in any sequences.", file=sys.stderr)
        
        # Prepare lead records for manual input - use RENAMED IDs to match FASTA
        if leads:
            leads_for_manual = [(f"lead_{rid}", seq) for rid, seq in leads]
        else:
            leads_for_manual = [records[0]]
        
        manual_cdr_rows = prompt_for_manual_cdr_input(leads_for_manual)
        
        if manual_cdr_rows:
            cdr_rows = manual_cdr_rows  # Replace empty/failed rows with manual input
            
            # Recheck if we now have CDRs
            has_any_cdrs = any(row.get("cdr1_seq") or row.get("cdr2_seq") or row.get("cdr3_seq") 
                            for row in cdr_rows)
            
            if has_any_cdrs:
                print("\n[INFO] Manual CDR input accepted. CDR-aware sorting will be available.", file=sys.stderr)

    # Write CDR CSV AFTER manual input has been collected
    write_cdr_csv(cdr_path, cdr_rows)
    print(f"✓ CDR CSV → {cdr_path}")

    # --- call the renderer
    renderer = shutil.which("python") or sys.executable
    script = os.path.join(os.path.dirname(__file__), "align_vs_lead_clear3_antpack_legend_v7.py")
    if not os.path.exists(script):
        sys.exit(f"ERROR: align_vs_lead_clear3_antpack_legend_v6.py not found next to this wrapper.")

    cmd = [
        renderer, script,
        "-i", fa_path,
        "--sto", sto_path,
        "-o", png_path,
        "--wrap", wrap_length,
        "--dpi", "350",
        "--lead-prefix", lead_prefix,
        "--show-consensus",
        "--auto-sort", args.auto_sort,
    ]
    
    # Always add CDR CSV (renderer will handle empty CDRs gracefully)
    if cdr_rows:
        cmd.extend(["--cdr-antpack-csv", cdr_path])

    print("→ running:", " ".join(cmd))
    r = subprocess.run(cmd)
    if r.returncode != 0:
        sys.exit(f"Renderer failed with exit code {r.returncode}")

    print(f"\n{'='*60}")
    print("ALL OUTPUTS COMPLETED SUCCESSFULLY")
    print(f"{'='*60}")
    print(f"Output directory: {out_dir}")
    print(f" - FASTA:     {fa_path}")
    print(f" - Stockholm: {sto_path}")
    print(f" - CDR CSV:   {cdr_path}")
    if has_any_cdrs:
        print(f"   (CDR regions detected)")
    else:
        print(f"   (No CDR regions detected - framework-only analysis)")
    print(f" - Figure:    {png_path}")
    print(f"{'='*60}")

if __name__ == "__main__":
    main()

    #python csv2msa_antpack_v5.py --table 20251103_M69_V_Saerens_Scfv.csv --scheme imgt