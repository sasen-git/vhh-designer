# models/analysis/

Analysis run outputs from `vhh_analysis_unified_v7.7.py`. Each timestamped folder
is a separate run with different parameters.

## Runs

| Folder | Positions Analyzed | Notes |
|--------|-------------------|-------|
| `v7.4_2026-01-12/` | Standard (vernier only) | Earlier unified version |
| `v7.7_2026-01-13_221057/` | **All positions** | Full FR1-FR4 analysis |
| `v7.7_2026-01-14_010118/` | Vernier only | v7.7, vernier-focused |
| `v7.7_2026-01-14_073437/` | Standard | v7.7, standard positions |

## Files in Each Run

### Used by the Designer (`vhh_designer_v9_1.py`)

| File | Size | Designer arg | What it contains |
|------|------|-------------|------------------|
| `analysis_rules_v7.json` | 1–5 MB | `--rules` | All compensation + epistasis rules. Each rule has a condition (e.g. `hallmarks=FERG`, `cdr3_len=short`, `cdr3[-1]=S`) and suggested mutations with confidence scores. This is the main model output. |
| `analysis_vernier_archetypes_v7.json` | ~13 KB | `--archetypes` | Per-family vernier consensus patterns (e.g. F_C2, Y_C2, VH_like) with consensus AA, frequency, and distribution at each IMGT position. Used for consensus-based mutations. |

### Informational Only (not fed into designer)

| File | Size | What it contains |
|------|------|------------------|
| `analysis_correlations_v7.json` | 40–300 KB | Per-position conservation stats (entropy, major AA frequency) and CDR3↔FR correlation coefficients (e.g. cdr3_len vs hydrophobicity). This is the statistical evidence used to *derive* the rules — the rules file is what the designer actually consumes. |
| `analysis_summary_v7.json` | ~3 KB | Run metadata: family counts, sequence totals, parameters used. |
| `run_config.json` | ~3 KB | Exact CLI parameters and settings for the run. |
| `analysis_vernier_clusters_v7.json` | **~1 GB** | Raw vernier cluster assignments for all sequences. Too large for git — kept locally only. |

### Legacy PKL Files (not in git, kept locally)

| File | Size | What it is |
|------|------|-----------|
| `checkpoint.pkl` | 600–750 MB | Full analyzer state — can resume interrupted runs |
| `analysis_epistasis_legacy.pkl` | 577–615 MB | Old-format epistasis output for backward compatibility |
| `analysis_compensation_legacy.pkl` | 8–74 KB | Old-format compensation output |

These PKL files are gitignored due to size. They exist on the local machine in each run folder.

## How They Were Generated

```bash
python active/analysis/vhh_analysis_unified_v7.7.py \
    --csv data/databases/annotated/vhh_full_annotated_v4_dedup_indiv_imgtfull/*.csv \
    -o models/analysis/v7.7_TIMESTAMP \
    --positions all_positions \    # or "vernier" for vernier-only
    --filter-truncated \
    --mi-max-seqs 2000000
```

Input: ~12M IMGT-annotated VHH sequences from `data/databases/annotated/`

## How They Are Consumed

```bash
python active/analysis/vhh_designer_v9_1.py \
    -i input.fasta \
    --rules models/analysis/v7.7_2026-01-13_221057/analysis_rules_v7_all_positions.json \
    --archetypes models/analysis/v7.7_2026-01-13_221057/analysis_vernier_archetypes_v7.json \
    --n-generate 100000 \
    --n-select 198
```

## Pipeline Summary

```
vhh_analysis_unified_v7.7.py
    ├── analysis_rules_v7.json ─────────┐
    ├── analysis_vernier_archetypes_v7.json ──→ vhh_designer_v9_1.py → designed VHHs
    ├── analysis_correlations_v7.json    (informational — not consumed)
    ├── analysis_summary_v7.json         (informational — not consumed)
    ├── analysis_vernier_clusters_v7.json (1GB, local only)
    ├── checkpoint.pkl                   (local only)
    ├── analysis_epistasis_legacy.pkl    (local only, backward compat)
    └── analysis_compensation_legacy.pkl (local only, backward compat)
```

## Note on "Correlations" vs "Rules"

The correlations file contains raw statistical associations (e.g. "at IMGT position 71,
hydrophobic residues correlate with longer CDR3 loops with r=0.34"). The rules file
contains the actionable output derived from those correlations (e.g. "if cdr3_len=long
and hallmarks=FERG, suggest V at IMGT71 with 87% confidence"). The designer only
needs the rules — the correlations are for human inspection and debugging.
