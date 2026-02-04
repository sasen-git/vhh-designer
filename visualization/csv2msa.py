#!/usr/bin/env python3
"""
csv2msa_antpack_v6.py

v6: Uses v8 renderer with IMGT H/V position markers
"""

import os, sys, re, csv, argparse, subprocess, shutil

os.environ.setdefault("MPLBACKEND", "Agg")

KNOWN_TAGS = {
    "MSKIK": {"name": "MSKIK", "position": "N-term"},
    "YPYDVPDYAYPYDVPDYAYPYDVPDYA": {"name": "3xHA", "position": "both"},
    "WSHPQFEK": {"name": "Strep-tag", "position": "C-term"},
    "EPEA": {"name": "C-tag", "position": "C-term"},
    "HHHHHH": {"name": "6xHis-tag", "position": "both"},
    "HHHHHHHH": {"name": "8xHis-tag", "position": "both"},
}

def detect_tags_in_sequence(seq):
    found_tags = {}
    seq_upper = seq.upper()
    for tag_seq, tag_info in KNOWN_TAGS.items():
        positions = []
        tag_name = tag_info["name"]
        start = 0
        while True:
            idx = seq_upper.find(tag_seq, start)
            if idx == -1: break
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
            if not rid or not seq: continue
            rows.append((is_lead, rid, seq))
        return rows

def write_fasta(path, records):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w") as f:
        for rid, seq in records:
            f.write(f">{rid}\n{seq}\n")

def antpack_cdrs_for_leads(leads, scheme="imgt", chains="H"):
    try:
        from antpack import SingleChainAnnotator
    except ImportError:
        return [{"id": rid, "cdr1_seq": "", "cdr2_seq": "", "cdr3_seq": ""} for rid, seq in leads]
    annot = SingleChainAnnotator(chains=[c.strip() for c in chains.split(",") if c.strip()], scheme=scheme)
    out_rows=[]
    for rid, seq in leads:
        try:
            seq_str = str(seq) if seq is not None else ""
            if not seq_str:
                out_rows.append({"id": rid, "cdr1_seq": "", "cdr2_seq": "", "cdr3_seq": ""})
                continue
            numbering, pid, chain, err = annot.analyze_seq(seq_str)
            if err or not numbering:
                out_rows.append({"id": rid, "cdr1_seq": "", "cdr2_seq": "", "cdr3_seq": ""})
                continue
            trimmed_seq, trimmed_nums, _, _ = annot.trim_alignment(seq_str, (numbering, pid, chain, err))
            labels = annot.assign_cdr_labels(trimmed_nums, chain)
            def spans(nums, labs):
                ans={}; cur=None; start=None
                for i,(n,lb) in enumerate(zip(nums,labs)):
                    if lb in ("cdr1","cdr2","cdr3"):
                        if cur is None: cur, start = lb, i
                    else:
                        if cur is not None: ans[cur]=(start,i-1); cur=None
                if cur is not None: ans[cur]=(start,len(labs)-1)
                return ans
            s = spans(trimmed_nums, labels)
            def extract(tag):
                if tag in s:
                    a,b = s[tag]
                    return trimmed_seq[a:b+1]
                return ""
            out_rows.append({"id": rid, "cdr1_seq": extract("cdr1"), "cdr2_seq": extract("cdr2"), "cdr3_seq": extract("cdr3")})
        except Exception as e:
            out_rows.append({"id": rid, "cdr1_seq": "", "cdr2_seq": "", "cdr3_seq": ""})
    return out_rows

def write_cdr_csv(path, rows):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["id","cdr1_seq","cdr2_seq","cdr3_seq"])
        w.writeheader()
        for r in rows: w.writerow(r)

def main():
    ap = argparse.ArgumentParser(description="CSV → FASTA, AntPack CDRs, MSA with IMGT markers")
    ap.add_argument("--table", required=True)
    ap.add_argument("--scheme", default="imgt", choices=["imgt","kabat","martin","aho"])
    ap.add_argument("--auto-sort", choices=["y","n","ask"], default="ask")
    args = ap.parse_args()

    csv_path = args.table
    if not os.path.exists(csv_path):
        sys.exit(f"ERROR: CSV not found: {csv_path}")

    script_base = os.path.splitext(os.path.basename(sys.argv[0]))[0]
    csv_dir  = os.path.dirname(csv_path) or "."
    csv_base = os.path.splitext(os.path.basename(csv_path))[0]
    out_dir  = os.path.join(csv_dir, f"{csv_base}_{script_base}")
    os.makedirs(out_dir, exist_ok=True)

    base     = os.path.join(out_dir, f"{csv_base}_{script_base}")
    fa_path  = f"{base}.fa"
    sto_path = f"{base}.sto"
    cdr_path = f"{base}_cdrs.csv"
    png_path = f"{base}.png"

    rows = read_threecol_csv(csv_path)
    all_seqs = [(rid, seq) for _, rid, seq in rows]
    wrap_length = "10000"

    leads = [(rid, seq) for is_lead, rid, seq in rows if is_lead]
    non_leads = [(rid, seq) for is_lead, rid, seq in rows if not is_lead]
    
    if not leads:
        lead_prefix = "___no_lead___"
        records = all_seqs
    else:
        lead_records = [(f"lead_{rid}", seq) for rid, seq in leads]
        records = lead_records + non_leads
        lead_prefix = "lead_"
        
    write_fasta(fa_path, records)
    print(f"✓ FASTA → {fa_path}")

    cdr_rows = []
    if leads:
        leads_for_antpack = [(f"lead_{rid}", seq) for rid, seq in leads]
        cdr_rows = antpack_cdrs_for_leads(leads_for_antpack, scheme=args.scheme, chains="H")
    elif records:
        cdr_rows = antpack_cdrs_for_leads([records[0]], scheme=args.scheme, chains="H")

    has_any_cdrs = any(row.get("cdr1_seq") or row.get("cdr2_seq") or row.get("cdr3_seq") for row in cdr_rows) if cdr_rows else False
    write_cdr_csv(cdr_path, cdr_rows)
    print(f"✓ CDR CSV → {cdr_path}")

    renderer = shutil.which("python") or sys.executable
    # v8 renderer with IMGT markers
    script = os.path.join(os.path.dirname(__file__), "align_vs_lead_clear3_antpack_legend_v8.py")
    if not os.path.exists(script):
        sys.exit(f"ERROR: align_vs_lead_clear3_antpack_legend_v8.py not found next to this wrapper.")

    cmd = [
        renderer, script,
        "-i", fa_path,
        "--sto", sto_path,
        "-o", png_path,
        "--wrap", wrap_length,
        "--dpi", "350",
        "--lead-prefix", lead_prefix,
        "--show-consensus",
        "--show-imgt",  # Enable IMGT H/V markers
        "--auto-sort", args.auto_sort,
    ]
    
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
    print(f" - Figure:    {png_path}")
    print(f"   (IMGT H=hallmark V=vernier markers enabled)")

if __name__ == "__main__":
    main()
