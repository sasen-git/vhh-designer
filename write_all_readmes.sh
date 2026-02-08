#!/usr/bin/env bash
set -euo pipefail

ROOT="/home/sasenefrem/KA-Search"
A="$ROOT/active/analysis"
AR="$ROOT/archive"
DRY_RUN=true

if [[ "${1:-}" == "--execute" ]]; then
    DRY_RUN=false
    echo "⚠️  EXECUTE MODE"
else
    echo "🔍 DRY RUN (pass --execute to apply)"
fi
echo ""

run_cmd() {
    if $DRY_RUN; then
        echo "  [DRY] $*"
    else
        echo "  [RUN] $*"
        eval "$@"
    fi
}

# ============================================================================
echo "========================================================================"
echo "  1. Archive superseded standalone IMGT scripts"
echo "========================================================================"

run_cmd mkdir -p "'$AR/standalone_imgt'"

for f in vhh_epistasis_v3_imgt.py \
         vhh_compensation_imgt_v6_csvcols.py \
         analyze_cdr_framework_correlations_imgt_csv_v3.py \
         convert_epistasis_to_designer.py; do
    if [ -f "$A/$f" ]; then
        run_cmd mv "'$A/$f'" "'$AR/standalone_imgt/'"
    else
        echo "  [SKIP] $f (not found or already moved)"
    fi
done

echo ""

# ============================================================================
echo "========================================================================"
echo "  2. Write READMEs for ALL folders"
echo "========================================================================"

if $DRY_RUN; then
    echo "  [DRY] Would write READMEs for:"
    echo "    models/README.md"
    echo "    models/epistasis/README.md"
    echo "    models/correlations/README.md"
    echo "    data/README.md"
    echo "    data/databases/README.md"
    echo "    data/raw/README.md"
    echo "    results/README.md"
    echo "    docs/README.md"
    echo "    tools/README.md"
    echo "    logs/README.md"
    echo "    legacy/README.md"
    echo "    active/analysis/README.md"
    echo "    active/database/README.md"
    echo "    active/utilities/README.md"
    echo "    archive/README.md"
else

# ===================== models/ =====================
cat > "$ROOT/models/README.md" << 'EOF'
# models/

Trained model outputs and statistical analysis results. All binary files (PKL, NPZ)
are gitignored. Only JSON config/rule files are tracked.

## Subfolders

| Folder | Contents | Status |
|--------|----------|--------|
| `analysis/` | Unified v7.x analysis run outputs (rules, archetypes, correlations) | **Production** — see analysis/README.md |
| `epistasis/` | Legacy epistasis pipeline outputs (pre-unified) | **Superseded** by analysis/ |
| `correlations/` | Legacy correlation/compensation pipeline outputs (pre-unified) | **Superseded** by analysis/ |
| `Testing-Walking-charge/` | Experimental charge-walking analysis — did not produce useful results | **Abandoned test** |

## Loose Files

| File | What it is |
|------|-----------|
| `family_fw_consensus.json` | Per-family framework consensus sequences (tracked) |
| `hallmark_fw_consensus.json` | Hallmark-grouped framework consensus sequences (tracked) |
| `claude.zip` | Archived Claude-generated scripts from early iterations |

## Current Production Model

The designer (`vhh_designer_v9_1.py`) reads from `analysis/`:

```bash
--rules   models/analysis/v7.7_2026-01-13_221057/analysis_rules_v7_all_positions.json
--archetypes models/analysis/v7.7_2026-01-13_221057/analysis_vernier_archetypes_v7.json
```

See `analysis/README.md` for full details on each JSON file.
EOF
echo "  [RUN] Wrote models/README.md"

# ===================== models/epistasis/ =====================
cat > "$ROOT/models/epistasis/README.md" << 'EOF'
# models/epistasis/

Legacy epistasis analysis outputs from the standalone `vhh_epistasis_overnight_*.py`
scripts (pre-unified pipeline). These used substring-indexed positions (not IMGT)
and PKL output format.

**Superseded by** `models/analysis/` which uses IMGT numbering and JSON output.

## Subfolders

| Folder | Contents |
|--------|----------|
| `current/` | epistasis_v2_full.pkl (36MB, substring positions), imgt_rules.json, correlation_results_v4_epistasis_only.pkl |
| `checkpoints/` | v2 intermediate checkpoints (compensation, clusters, models, rules) |
| `legacy/` | epistasis_overnight_full.pkl (12MB, v1 from Dec 7) |
| `legacy_checkpoints/` | v1 intermediate checkpoints |
| `imgt_v1/` | First attempt at IMGT epistasis — 185-byte checkpoint PKLs (effectively empty/failed runs) |
| `imgt_light_v1/` | Light IMGT attempt — 318-byte PKL (also failed) |

## Note

The `imgt_v1/` and `imgt_light_v1/` folders contain PKLs that are only ~185 bytes,
meaning the runs likely failed or produced no meaningful data. These were the early
attempts that led to building the unified `vhh_analysis_unified_v7.7.py` pipeline,
which successfully completed the IMGT analysis.
EOF
echo "  [RUN] Wrote models/epistasis/README.md"

# ===================== models/correlations/ =====================
cat > "$ROOT/models/correlations/README.md" << 'EOF'
# models/correlations/

Legacy correlation and compensation analysis outputs from standalone scripts.

**Superseded by** `models/analysis/` unified pipeline.

## Files

| File | Size | What it is |
|------|------|-----------|
| `compensation_results.pkl` | 22KB | v1 compensation output (Dec 6) |
| `correlation_results_v3_compensation.pkl` | 194KB | v3 merged compensation + correlation (Jan 2) |
| `correlation_results_v4_merged.pkl` | 533KB | v4 merged with epistasis (Jan 3) |
| `imgt_csv_v3_sanity/correlation_results_v3_summary.json` | tracked | Sanity check summary from IMGT CSV correlation run |
| `imgt_csv_v3_sanity/correlation_rules_v3.json` | tracked | Rules from IMGT CSV correlation sanity check |

## Note

These PKLs use the old schema (`multi_position_rules`, `top_predictors`). The
current designer (v9.1) expects the new JSON schema from `models/analysis/`.
These files are kept for reference but are not used in the current pipeline.
EOF
echo "  [RUN] Wrote models/correlations/README.md"

# ===================== data/ =====================
cat > "$ROOT/data/README.md" << 'EOF'
# data/ (51GB)

All sequence data for the VHH design pipeline. Entirely gitignored.

See also: `VHH_DATABASE_CONSOLIDATION_README.md` in this folder for detailed
provenance and consolidation history.

## Subfolders

| Folder | Size | Contents |
|--------|------|----------|
| `databases/` | ~47GB | Production NPZ shards, annotated CSVs, legacy databases |
| `archive/` | ~3GB | Older annotation versions (v3, pre-dedup) |
| `raw/` | ~1GB | Original source files before processing |

## Data Sources (~12M sequences total)

| Source | Sequences | Description |
|--------|-----------|-------------|
| OAS Camel | 1.35M | Observed Antibody Space — camelid repertoires |
| INDI NGS | 10.7M | INDI Molecular nanobody NGS data (6 shards) |
| Patents/Structures | 19K | Patent sequences, SAbDab/PDB structures, curated VHHs |

## Pipeline

```
raw/ sources → build scripts → databases/shards/*.npz
                             → annotate_all_shards_v8.py
                             → databases/annotated/*.csv (with IMGT positions)
                             → vhh_analysis_unified_v7.7.py → models/
```
EOF
echo "  [RUN] Wrote data/README.md"

# ===================== data/databases/ =====================
cat > "$ROOT/data/databases/README.md" << 'EOF'
# data/databases/

Production databases used by the pipeline.

## Subfolders

| Folder | Contents |
|--------|----------|
| `shards/` | 8 NPZ files (~12M sequences) — primary input for analysis pipeline |
| `production/` | VHH_db_unified_v2.npz (67MB) — consolidated single-file database |
| `annotated/` | IMGT-augmented CSV versions of shards |
| `legacy/` | Older database versions (VHH_db_final.npz, camel_vhh_clean_db.npz, etc.) |

## Shards (primary)

| File | Sequences | Source |
|------|-----------|--------|
| vhh_annotated.npz | 19,409 | Patents, structures, curated |
| vhh_oas_camel.npz | 1,354,537 | OAS Camel repertoire |
| vhh_indi_ngs_000.npz | 2,000,000 | INDI NGS batch 0 |
| vhh_indi_ngs_001.npz | 2,000,000 | INDI NGS batch 1 |
| vhh_indi_ngs_002.npz | 2,000,000 | INDI NGS batch 2 |
| vhh_indi_ngs_003.npz | 2,000,000 | INDI NGS batch 3 |
| vhh_indi_ngs_004.npz | 2,000,000 | INDI NGS batch 4 |
| vhh_indi_ngs_005.npz | 675,894 | INDI NGS batch 5 (remainder) |
| **Total** | **~12,049,840** | |

## Annotated CSV Versions

| Folder | Description |
|--------|-------------|
| `vhh_full_annotated_v4/` | Base annotation (AntPack numbering + family classification) |
| `vhh_full_annotated_v4_dedup/` | Deduplicated version + dedup temp buckets |
| `vhh_full_annotated_v4_dedup_indiv_imgtfull/` | Full IMGT position columns per shard — **primary input for unified analysis** |
EOF
echo "  [RUN] Wrote data/databases/README.md"

# ===================== data/raw/ =====================
cat > "$ROOT/data/raw/README.md" << 'EOF'
# data/raw/

Original source files before processing into the production database.

## Subfolders

| Folder | Contents |
|--------|----------|
| `oas_paper/compressed/` | OAS camelid bulk download (SRR3544217-SRR3544222 CSVs, gzipped) |
| `sequences/` | Input FASTA/XLSX files |

## Key Files

| File | What it is |
|------|-----------|
| `sequences/M69.fasta` | Primary input antibody (lead VH sequence for VHH conversion) |
| `sequences/M69_HC_LC.fasta` | Heavy + light chain version |
| `sequences/INDI_patent_sequences.xlsx` | Patent-derived VHH sequences |
| `sequences/INDI_structure_sequences.xlsx` | Structure-derived VHH sequences |
| `sequences/VHH_annotated_full.xlsx` | Curated annotated VHH collection |
| `oas_paper/compressed/vhh_sequences.csv.gz` | Processed VHH sequences from OAS |
| `oas_paper/compressed/all_nano_structures.zip` | Nanobody structural data |
EOF
echo "  [RUN] Wrote data/raw/README.md"

# ===================== results/ =====================
cat > "$ROOT/results/README.md" << 'EOF'
# results/

All pipeline output — designer runs, NPZ scans, analysis results. Entirely gitignored.

## Subfolders

| Folder | Contents |
|--------|----------|
| `design_runs/` | 60 timestamped designer output directories (v7.5 through v9.1) |
| `analysis_runs/` | Early designer runs (v2-v7) + naturalness analysis outputs |
| `npz_scans/` | NPZ database scan results organized by date |
| `analysis_outputs/` | Loose analysis CSVs/XLSXs moved from root |

## Design Runs

Each `design_runs/YYYYMMDD_HHMMSS_vX_Y_nNNNNN_M69/` folder contains:
- `*.csv` — selected candidates
- `*_all.csv` — all generated candidates
- `*.fasta` — sequences in FASTA format
- `*_summary.json` — run configuration and statistics
- `*_csv2msa_antpack_*/` — MSA alignment subfolder (if post-processed)

## Notable Runs

| Run | Designer | Candidates |
|-----|----------|------------|
| `20260201_234207_v9_1_n100000_t08_M69/` | v9.1 (latest) | 100K generated |
| `20260201_221918_v9_0_n100000_t08_M69/` | v9.0 | 100K generated |
| `20260201_204445_v8_9_4_n100000_t08_M69/` | v8.9.4 | 100K generated |
EOF
echo "  [RUN] Wrote results/README.md"

# ===================== docs/ =====================
cat > "$ROOT/docs/README.md" << 'EOF'
# docs/

Project documentation.

## Files

| File | What it is |
|------|-----------|
| `VHH_DESIGNER_ULTIMATE_GUIDE_V90.md` | Comprehensive guide to the VHH designer pipeline — architecture, usage, biology background |
| `scripts_log.md` | Chronological log of script development and changes |
| `REORGANIZATION_PROPOSAL.md` | Initial repo reorganization plan |
| `CHANGELOG_epistasis_pipeline.md` | Epistasis pipeline version history |
| `CHANGELOG_full_analysis.md` | Full analysis version history |
| `CHANGELOG_naturalness_analyzer.md` | Naturalness analyzer version history |
| `CHANGELOG_npz_scanner.md` | NPZ scanner version history |
EOF
echo "  [RUN] Wrote docs/README.md"

# ===================== tools/ =====================
cat > "$ROOT/tools/README.md" << 'EOF'
# tools/

Maintenance and QC utilities.

## Scripts

| Script | Purpose |
|--------|---------|
| `dedup_annotated_imgt.py` | Deduplicates IMGT-annotated CSV shards by sequence identity |
| `dedupe_csv_headers.py` | Fixes duplicate column headers in CSVs |
| `qc_imgt_coverage_fast.py` | Quick QC check: what % of IMGT positions are filled per shard |
| `diagnose_unknown_files.py` | Identifies and categorizes unknown files in the repo |
| `kasearch_paths.py` | Centralized path management for all KA-Search scripts |
| `kasearch_reorganize.py` | Repo reorganization helper |
| `maintain.py` | General maintenance (cleanup temp files, check consistency) |
| `summarize_sweep.py` | Summarizes parameter sweep results across multiple runs |
| `pathlib_list_v2.py` | File listing with filtering |

## Other

| File | What it is |
|------|-----------|
| `.cache.json` | Cached metadata for maintenance scripts |
EOF
echo "  [RUN] Wrote tools/README.md"

# ===================== logs/ =====================
cat > "$ROOT/logs/README.md" << 'EOF'
# logs/

Runtime logs from pipeline execution. Gitignored.

Contains stdout/stderr captures from long-running analysis and designer jobs.
EOF
echo "  [RUN] Wrote logs/README.md"

# ===================== legacy/ =====================
cat > "$ROOT/legacy/README.md" << 'EOF'
# legacy/

Placeholder for very old/deprecated items. Currently mostly empty (~8KB).
Bulk historical scripts are in `archive/` instead.
EOF
echo "  [RUN] Wrote legacy/README.md"

# ===================== active/analysis/ =====================
cat > "$A/README.md" << 'EOF'
# active/analysis/

Current production scripts for the VHH design pipeline.

## Core Pipeline

| Script | Purpose |
|--------|---------|
| `vhh_designer_v9_1.py` | Main VHH designer — converts VH→VHH using learned rules (JSON input) |
| `vhh_analysis_unified_v7.7.py` | Unified epistasis + compensation + correlation analysis (IMGT, JSON output) |
| `vhh_naturalness_analyzer_v4.py` | Scores sequences for naturalness/evolutionary fitness |
| `annotate_all_shards_v8.py` | Annotates database shards with IMGT positions |

## Supporting Scripts

| Script | Purpose |
|--------|---------|
| `augment_imgt_positions_from_csvs_v2.py` | Augments CSVs with full IMGT position columns |
| `charge_balance_vhh_windows.py` | Charge balance analysis across VHH windows |
| `extract_script_headers.py` | Extracts docstrings/headers from all scripts |
| `patches.py` | Hotfix patches for data issues |
| `run_with_positions.py` | Runs analysis with specific position sets |

## Data Flow

```
NPZ shards (12M sequences)
    ↓
annotate_all_shards_v8.py → IMGT-augmented CSVs
    ↓
vhh_analysis_unified_v7.7.py → JSON rules (epistasis + compensation + correlation)
    ↓
vhh_designer_v9_1.py → designed VHH sequences (CSV + FASTA)
    ↓
../utilities/csv2msa_antpack_v6.py → MSA visualization
vhh_naturalness_analyzer_v4.py → fitness scoring
```

## Version History

All older iterations are in `archive/` organized by component:
- `archive/designer/` — v2 through v9_0
- `archive/analysis_unified/` — v7 through v7.6
- `archive/compensation/` — v1 through v5
- `archive/standalone_imgt/` — pre-unified IMGT scripts (epistasis_v3, compensation_v6, correlations_v3)
EOF
echo "  [RUN] Wrote active/analysis/README.md"

# ===================== active/database/ =====================
cat > "$ROOT/active/database/README.md" << 'EOF'
# active/database/

## Scripts

| Script | Purpose |
|--------|---------|
| `npz_fullscan_v6_integrated.py` | Scans NPZ databases for sequence matches with integrated filtering |

## Database Location

Production databases are in `data/databases/`:
- `shards/` — 8 NPZ files (~12M sequences total)
- `production/` — consolidated VHH_db_unified_v2.npz
- `annotated/` — IMGT-augmented CSV versions
- `legacy/` — older database versions

See `data/VHH_DATABASE_CONSOLIDATION_README.md` for full details.
EOF
echo "  [RUN] Wrote active/database/README.md"

# ===================== active/utilities/ =====================
cat > "$ROOT/active/utilities/README.md" << 'EOF'
# active/utilities/

Post-processing and helper tools.

## MSA & Visualization

| Script | Purpose |
|--------|---------|
| `csv2msa_antpack_v6.py` | Converts designer CSV → MSA with AntPack numbering (latest) |
| `csv2msa_antpack_v5_imgt_v2.py` | MSA converter, IMGT numbering variant |
| `csv2msa_antpack_v5.py` | Earlier MSA converter |
| `align_vs_lead_clear3_antpack_legend_v8.py` | Alignment visualization vs lead sequence (latest) |
| `align_vs_lead_clear3_antpack_legend_v8_imgt.py` | Same, IMGT variant |
| `align_vs_lead_clear3_antpack_legend_v7.py` | Earlier alignment visualization |
| `align_vs_lead_clear3_antpack_legend_v7_imgt_v2.py` | Earlier, IMGT variant |

## Sequence Tools

| Script | Purpose |
|--------|---------|
| `dna_translator.py` | Translates nucleotide → amino acid sequences |
| `pull_cdrs.py` | Extracts CDR regions from annotated sequences |
| `scfv_anarci.py` | ANARCI numbering for scFv sequences |
| `pathlib_list_v2.py` | File listing utility |
EOF
echo "  [RUN] Wrote active/utilities/README.md"

# ===================== archive/ =====================
cat > "$AR/README.md" << 'EOF'
# archive/

Historical versions of all pipeline components, organized by function.

| Folder | Contents |
|--------|----------|
| `designer/` | VHH designer v2–v9_0 (~40 versions) |
| `analysis_unified/` | Unified analysis v7–v7.6 |
| `standalone_imgt/` | Pre-unified IMGT scripts (superseded by unified v7.7) |
| `compensation/` | Compensation analysis v1–v5 |
| `annotate_shards/` | Shard annotation v1–v7 |
| `epistasis_pipeline/` | Epistasis overnight v1, v2 |
| `epistasis_to_imgt/` | IMGT conversion scripts v1–v3 |
| `naturalness_analyzer/` | Naturalness analyzer v1–v3, foldability |
| `npz_scanner/` | NPZ fullscan v2–v9 |
| `full_analysis/` | Full analysis runners v1–v5 |
| `correlation_analysis/` | CDR-framework correlation scripts |
| `database_builders/` | Database build scripts |
| `visualizations/` | Visualization scripts |
| `one_off/` | One-time diagnostic/processing scripts |
| `augment_imgt/` | IMGT augmentation v1 |
| `Backup/` | Data backups (annotated CSVs) |
| `misc/` | Miscellaneous (DXF, zips) |
EOF
echo "  [RUN] Wrote archive/README.md"

fi

echo ""

# ============================================================================
echo "========================================================================"
echo "  3. Git commit + push"
echo "========================================================================"

if $DRY_RUN; then
    echo "  [DRY] Would commit and push"
else
    cd "$ROOT"

    # Force-add READMEs in gitignored directories
    git add -f models/README.md
    git add -f models/epistasis/README.md
    git add -f models/correlations/README.md
    git add -f models/analysis/README.md
    git add -f data/README.md
    git add -f data/databases/README.md
    git add -f data/raw/README.md
    git add -f results/README.md
    git add -f logs/README.md
    git add -f legacy/README.md

    # These aren't gitignored so normal add works
    git add active/analysis/README.md
    git add active/database/README.md
    git add active/utilities/README.md
    git add archive/README.md
    git add docs/README.md
    git add tools/README.md

    # Add archived scripts
    git add -A

    echo ""
    git status --short
    TOTAL=$(git status --short | wc -l)
    echo "  $TOTAL files changed"
    echo ""

    read -p "Commit and push? [y/N] " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        git commit -m "Add READMEs for all folders, archive superseded IMGT scripts

READMEs added for:
  models/ — overview + Testing-Walking-charge (abandoned test), claude.zip
  models/epistasis/ — legacy PKL pipeline (superseded by analysis/)
  models/correlations/ — legacy correlation PKLs (superseded by analysis/)
  data/ — 51GB overview, data sources, pipeline
  data/databases/ — shard index, annotated CSV versions
  data/raw/ — original source files
  results/ — design runs, NPZ scans, analysis outputs
  docs/ — documentation index
  tools/ — maintenance/QC scripts
  logs/ — runtime logs
  legacy/ — placeholder
  active/analysis/ — current production scripts + data flow
  active/database/ — NPZ scanner
  active/utilities/ — MSA converters, alignment tools
  archive/ — historical version index

Archived to archive/standalone_imgt/:
  vhh_epistasis_v3_imgt.py
  vhh_compensation_imgt_v6_csvcols.py
  analyze_cdr_framework_correlations_imgt_csv_v3.py
  convert_epistasis_to_designer.py
(All superseded by vhh_analysis_unified_v7.7.py)"

        git push
        echo "✅ Pushed!"
    fi
fi

echo ""
echo "========================================================================"
if $DRY_RUN; then
    echo "  DRY RUN complete. Re-run with --execute to apply."
else
    echo "  Done!"
fi
echo "========================================================================"
