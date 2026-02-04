#!/usr/bin/env python3
from __future__ import annotations

import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# Where your files live
SEQ_DIR = Path("/home/sasenefrem/KA-Search/data/raw/sequences")

# -----------------------------
# Utilities
# -----------------------------
AA_RE = re.compile(r"[^A-Z]")

def norm_col(s: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", s.strip().lower())

ID_CANDIDATES = {
    "id", "name", "sequencename", "sequence_name", "seqname", "header", "recordid", "record_id"
}
SEQ_CANDIDATES = {
    "seq", "sequence", "aa", "aaseq", "aminoacidsequence", "proteinsequence", "protein", "aminoacids"
}

LINKER_RE = re.compile(r"^[GSA T]{0,30}$".replace(" ", ""))  # small GSAT linkers

@dataclass
class TagRemoval:
    tag: str
    terminal: str  # "N" or "C"
    removed: str

def strip_terminal_tags(seq: str) -> Tuple[str, List[TagRemoval]]:
    """
    Conservative: only strips known tags if they appear flush at N- or C-terminus (optionally with short GSAT linkers).
    """
    removals: List[TagRemoval] = []
    s = seq

    def strip_n(prefix: str, tagname: str) -> bool:
        nonlocal s, removals
        if s.startswith(prefix):
            s = s[len(prefix):]
            removals.append(TagRemoval(tagname, "N", prefix))
            # also remove a short linker immediately after
            m = re.match(r"^[GSAT]{1,20}", s)
            if m and LINKER_RE.match(m.group(0)):
                link = m.group(0)
                s = s[len(link):]
                removals.append(TagRemoval("linker", "N", link))
            return True
        return False

    def strip_c(suffix: str, tagname: str) -> bool:
        nonlocal s, removals
        if s.endswith(suffix):
            s = s[:-len(suffix)]
            removals.append(TagRemoval(tagname, "C", suffix))
            # also remove a short linker immediately before
            # (only if it's clearly a GSAT linker at the end)
            m = re.search(r"[GSAT]{1,20}$", s)
            if m and LINKER_RE.match(m.group(0)):
                link = m.group(0)
                s = s[:-len(link)]
                removals.append(TagRemoval("linker", "C", link))
            return True
        return False

    # Known tags
    HA = "YPYDVPDYA"
    FLAG = "DYKDDDDK"
    STREP2 = "WSHPQFEK"
    CTAG = "EPEA"

    # N-term: SKIK/MSKIK helpers
    changed = True
    while changed:
        changed = False
        if strip_n("MSKIK", "MSKIK"):
            changed = True
        if strip_n("SKIK", "SKIK"):
            changed = True

    # C-term: repeated HA motifs (up to 3), with optional short linkers between them
    # We'll strip from the very end only.
    def strip_ha_repeats() -> bool:
        nonlocal s, removals
        original = s
        count = 0
        while True:
            # allow optional linker before HA at the very end
            m = re.search(r"([GSAT]{0,20})(" + re.escape(HA) + r")$", s)
            if not m:
                break
            link, motif = m.group(1), m.group(2)
            if link and not LINKER_RE.match(link):
                break
            s = s[: len(s) - (len(link) + len(motif))]
            if link:
                removals.append(TagRemoval("linker", "C", link))
            removals.append(TagRemoval("HA", "C", motif))
            count += 1
            if count >= 3:
                break
        return s != original

    strip_ha_repeats()

    # TwinStrep common exact form (one-piece) OR StrepII
    # (Try longer first)
    twin_strep = STREP2 + "GGGSGGGSGGSA" + STREP2
    if strip_c(twin_strep, "TwinStrep"):
        pass
    else:
        strip_c(STREP2, "StrepII")

    # C-tag
    strip_c(CTAG, "C-tag(EPEA)")

    # His tag (6-10 H at C-terminus)
    m_his = re.search(r"(H{6,10})$", s)
    if m_his:
        his = m_his.group(1)
        s = s[:-len(his)]
        removals.append(TagRemoval("His-tag", "C", his))

    # FLAG (either terminal)
    if strip_c(FLAG, "FLAG"):
        pass
    if strip_n(FLAG, "FLAG"):
        pass

    return s, removals

def clean_sequence(raw: str) -> str:
    s = raw.strip().upper()
    s = AA_RE.sub("", s)
    return s

def list_sequence_files() -> List[Path]:
    exts = {".txt", ".fa", ".fasta", ".csv", ".xlsx", ".xls"}
    if not SEQ_DIR.exists():
        print(f"ERROR: Directory not found: {SEQ_DIR}")
        sys.exit(1)
    files = [p for p in SEQ_DIR.iterdir() if p.is_file() and p.suffix.lower() in exts]
    return sorted(files)

# -----------------------------
# Loaders
# -----------------------------
def load_fasta(path: Path) -> List[Tuple[str, str]]:
    records: List[Tuple[str, str]] = []
    cur_id: Optional[str] = None
    cur_seq: List[str] = []
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if cur_id is not None:
                records.append((cur_id, "".join(cur_seq)))
            cur_id = line[1:].strip() or f"record_{len(records)+1}"
            cur_seq = []
        else:
            cur_seq.append(line)
    if cur_id is not None:
        records.append((cur_id, "".join(cur_seq)))
    return records

def load_txt(path: Path) -> List[Tuple[str, str]]:
    text = path.read_text(encoding="utf-8", errors="ignore")
    if ">" in text:
        # treat as FASTA-like
        return load_fasta(path)

    records: List[Tuple[str, str]] = []
    lines = [ln.strip() for ln in text.splitlines() if ln.strip()]
    for i, ln in enumerate(lines, start=1):
        # If line looks like "ID<tab>SEQUENCE" or "ID SEQUENCE"
        parts = re.split(r"[\t,; ]+", ln, maxsplit=1)
        if len(parts) == 2 and len(clean_sequence(parts[1])) > 50:
            rid, seq = parts[0], parts[1]
        else:
            rid, seq = f"seq_{i:04d}", ln
        records.append((rid, seq))
    return records

def paste_sequence_record() -> Tuple[str, str]:
    print("\nPaste an AA sequence (single-line or multi-line).")
    print("Finish by pressing ENTER on an empty line.\n")

    lines = []
    while True:
        try:
            ln = input()
        except EOFError:
            break
        if ln.strip() == "":
            break
        lines.append(ln.strip())

    raw = "".join(lines).strip()
    if not raw:
        raise ValueError("No sequence pasted.")
    return ("PASTED_SEQ", raw)


def detect_columns(cols: List[str]) -> Tuple[Optional[str], Optional[str]]:
    norm_map = {c: norm_col(c) for c in cols}
    id_col = None
    seq_col = None

    for c, n in norm_map.items():
        if n in ID_CANDIDATES:
            id_col = c
            break

    for c, n in norm_map.items():
        if n in SEQ_CANDIDATES:
            seq_col = c
            break

    return id_col, seq_col

def load_table(path: Path) -> List[Tuple[str, str]]:
    import pandas as pd

    if path.suffix.lower() == ".csv":
        df = pd.read_csv(path)
    else:
        df = pd.read_excel(path)

    cols = list(df.columns)
    id_col, seq_col = detect_columns(cols)

    print("\nDetected columns:")
    print("  " + ", ".join([str(c) for c in cols]))
    print(f"\nGuessed ID column:  {id_col}")
    print(f"Guessed Seq column: {seq_col}")

    # Ask for confirmation / manual override
    ok = input("Is this correct? [y/n] ").strip().lower()
    if ok != "y":
        id_col = input("Type the ID column name exactly (or leave blank to auto-generate IDs): ").strip() or None
        seq_col = input("Type the sequence column name exactly: ").strip()
        if seq_col not in df.columns:
            raise ValueError(f"Sequence column not found: {seq_col}")
        if id_col is not None and id_col not in df.columns:
            raise ValueError(f"ID column not found: {id_col}")

    records: List[Tuple[str, str]] = []
    for i, row in df.iterrows():
        rid = str(row[id_col]).strip() if id_col else f"row_{i+1}"
        seq = str(row[seq_col]).strip()
        records.append((rid, seq))
    return records

def load_any(path: Path) -> List[Tuple[str, str]]:
    suf = path.suffix.lower()
    if suf in {".fa", ".fasta"}:
        return load_fasta(path)
    if suf == ".txt":
        return load_txt(path)
    if suf in {".csv", ".xlsx", ".xls"}:
        return load_table(path)
    raise ValueError(f"Unsupported file type: {path}")

# -----------------------------
# Annotation (AntPack + ANARCI)
# -----------------------------
def annotate_imgt_with_antpack(seq: str) -> Dict[str, str]:
    """
    Returns dict with keys:
      full, fwrk1, cdr1, fwrk2, cdr2, fwrk3, cdr3, fwrk4
    """
    from antpack import SingleChainAnnotator

    sc = SingleChainAnnotator(["H", "K", "L"], scheme="imgt")
    numbering, pid, chain, err = sc.analyze_seq(seq)
    if err:
        raise RuntimeError(f"AntPack error: {err} (pid={pid}, chain={chain})")

    trimmed_seq, trimmed_numbering, exstart, exend = sc.trim_alignment(seq, (numbering, pid, chain, err))
    labels = sc.assign_cdr_labels(trimmed_numbering, chain=chain, scheme="imgt")

    regions = {
        "full": trimmed_seq,
        "fmwk1": [],
        "cdr1": [],
        "fmwk2": [],
        "cdr2": [],
        "fmwk3": [],
        "cdr3": [],
        "fmwk4": [],
    }
    for aa, lab in zip(trimmed_seq, labels):
        if lab in regions:
            regions[lab].append(aa)

    return {
        "full": regions["full"],
        "fwrk1": "".join(regions["fmwk1"]),
        "cdr1": "".join(regions["cdr1"]),
        "fwrk2": "".join(regions["fmwk2"]),
        "cdr2": "".join(regions["cdr2"]),
        "fwrk3": "".join(regions["fmwk3"]),
        "cdr3": "".join(regions["cdr3"]),
        "fwrk4": "".join(regions["fmwk4"]),
    }

def try_run_anarci(seq: str) -> Optional[str]:
    """
    Just attempts ANARCI IMGT numbering; returns a short status string if available.
    """
    try:
        from anarci import anarci  # type: ignore
    except Exception:
        return None

    try:
        # based on ANARCI python usage examples/issues
        numbering, alignment_details, hit_tables = anarci(
            sequences=[("query", seq)],
            scheme="imgt",
            output=False,
        )
        # If it got here, it ran. We keep it minimal.
        return "ANARCI: OK (IMGT numbering ran)"
    except Exception as e:
        return f"ANARCI: FAILED ({e})"

# -----------------------------
# Interactive selection
# -----------------------------
def pick_file_or_paste(files: List[Path]) -> Tuple[str, Optional[Path]]:
    """
    Returns ("file", path) or ("paste", None)
    """
    print(f"\nLooking in: {SEQ_DIR}\n")
    for i, p in enumerate(files, start=1):
        print(f"{i:2d}) {p.name}")
    print(" P) Paste a sequence (no file)")

    while True:
        choice = input("\nChoose a file number, or 'P' to paste: ").strip().lower()
        if choice == "p":
            return ("paste", None)
        if choice.isdigit() and 1 <= int(choice) <= len(files):
            return ("file", files[int(choice) - 1])
        print("Invalid choice. Try again.")


def pick_sequences(ids: List[str]) -> List[str]:
    print("\nFirst 20 sequence IDs:")
    show = ids[:20]
    for i, rid in enumerate(show, start=1):
        print(f"{i:2d}) {rid}")

    print("\nPick sequences:")
    print("  - Enter comma-separated numbers (e.g. 1,3,7)")
    print("  - OR type an ID exactly")
    print("  - OR type 'q' to quit")
    s = input("> ").strip()
    if s.lower() == "q":
        return []

    if re.fullmatch(r"[0-9,\s]+", s):
        nums = [int(x) for x in re.split(r"[,\s]+", s) if x.strip()]
        chosen = []
        for n in nums:
            if 1 <= n <= len(show):
                chosen.append(show[n - 1])
        return chosen

    # otherwise treat as ID
    return [s]

def main() -> None:
    files = list_sequence_files()
    if not files:
        print(f"No sequence files found in {SEQ_DIR} (.txt/.fa/.fasta/.csv/.xlsx).")
        # still allow paste even if there are no files
        files = []

    mode, path = pick_file_or_paste(files)

    if mode == "paste":
        try:
            rid, raw = paste_sequence_record()
            records = [(rid, raw)]
            print("\nUsing pasted sequence.\n")
        except Exception as e:
            print(f"ERROR: {e}")
            sys.exit(1)
    else:
        print(f"\nLoading: {path}\n")
        try:
            records = load_any(path)
        except Exception as e:
            print(f"ERROR loading file: {e}")
            sys.exit(1)

   
    # Clean + tag-strip + filter by length
    seqs: Dict[str, str] = {}
    tag_notes: Dict[str, List[TagRemoval]] = {}
    too_short: List[str] = []

    for rid, raw in records:
        seq = clean_sequence(raw)
        if not seq:
            continue

        cleaned, removals = strip_terminal_tags(seq)

        if len(cleaned) <= 100:
            too_short.append(rid)
            continue

        # If duplicate IDs, make unique
        base = rid
        k = 2
        while rid in seqs:
            rid = f"{base}__{k}"
            k += 1

        seqs[rid] = cleaned
        if removals:
            tag_notes[rid] = removals

    if not seqs:
        print("No usable sequences found after cleaning/tag stripping/length filter (>100 aa).")
        sys.exit(1)

    if too_short:
        print(f"Skipped {len(too_short)} sequences with length <= 100 aa (after cleaning/tag stripping).")

    ids = list(seqs.keys())

    chosen_ids = pick_sequences(ids)
    if not chosen_ids:
        print("Done.")
        return

    # If user typed an ID not in the first 20, try to find it
    resolved: List[str] = []
    for cid in chosen_ids:
        if cid in seqs:
            resolved.append(cid)
            continue
        # case-insensitive match
        hits = [rid for rid in ids if rid.lower() == cid.lower()]
        if len(hits) == 1:
            resolved.append(hits[0])
            continue
        # substring match
        hits = [rid for rid in ids if cid.lower() in rid.lower()]
        if len(hits) == 1:
            resolved.append(hits[0])
            continue
        if hits:
            print(f"\nMultiple matches for '{cid}':")
            for h in hits[:20]:
                print(f"  - {h}")
            pick = input("Type the exact ID you want: ").strip()
            if pick in seqs:
                resolved.append(pick)
            else:
                print(f"Could not resolve '{cid}'. Skipping.")
        else:
            print(f"Could not find ID '{cid}'. Skipping.")

    if not resolved:
        print("No valid selections. Done.")
        return

    # Check deps early
    try:
        import antpack  # noqa: F401
    except Exception as e:
        print("\nERROR: antpack is not importable in this environment.")
        print("Install it (and set up licensing if needed):")
        print("  pip install antpack")
        print("  AntPack-setup")
        print(f"Details: {e}")
        sys.exit(1)

    for rid in resolved:
        seq = seqs[rid]
        print("\n" + "=" * 80)
        print(f"ID: {rid}")
        print(f"Length (cleaned): {len(seq)}")

        # Tags removed?
        if rid in tag_notes:
            print("\nTags removed:")
            for r in tag_notes[rid]:
                print(f"  - {r.terminal}-term: {r.tag}  ({r.removed})")
        else:
            print("\nTags removed: none detected")

        # ANARCI sanity check (optional)
        anarci_status = try_run_anarci(seq)
        if anarci_status:
            print(f"\n{anarci_status}")

        # AntPack extraction
        try:
            out = annotate_imgt_with_antpack(seq)
        except Exception as e:
            print(f"\nAntPack annotation failed: {e}")
            continue

        print("\nIMGT regions (AntPack):")
        print(f"\nFULL:\n{out['full']}")
        print(f"\nFR1:\n{out['fwrk1']}")
        print(f"\nCDR1:\n{out['cdr1']}")
        print(f"\nFR2:\n{out['fwrk2']}")
        print(f"\nCDR2:\n{out['cdr2']}")
        print(f"\nFR3:\n{out['fwrk3']}")
        print(f"\nCDR3:\n{out['cdr3']}")
        print(f"\nFR4:\n{out['fwrk4']}")

    print("\nDone.\n")

if __name__ == "__main__":
    main()
