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
    fi
done

echo ""
echo "  [KEEP] vhh_analysis_unified_v7.7.py"
echo "  [KEEP] vhh_naturalness_analyzer_v4.py"
echo "  [KEEP] vhh_designer_v9_1.py"
echo ""

# ============================================================================
echo "========================================================================"
echo "  2. Write local READMEs"
echo "========================================================================"

if $DRY_RUN; then
    echo "  [DRY] Would write active/analysis/README.md"
    echo "  [DRY] Would write active/database/README.md"
    echo "  [DRY] Would write active/utilities/README.md"
    echo "  [DRY] Would write archive/README.md"
else

# --- active/analysis/README.md ---
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

# --- active/database/README.md ---
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

# --- active/utilities/README.md ---
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

# --- archive/README.md ---
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
    git add -A
    echo ""
    git status --short
    TOTAL=$(git status --short | wc -l)
    echo "  $TOTAL files changed"
    echo ""

    read -p "Commit and push? [y/N] " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        git commit -m "Archive superseded IMGT scripts, add local READMEs

Archived to archive/standalone_imgt/:
  vhh_epistasis_v3_imgt.py
  vhh_compensation_imgt_v6_csvcols.py
  analyze_cdr_framework_correlations_imgt_csv_v3.py
  convert_epistasis_to_designer.py
(All superseded by vhh_analysis_unified_v7.7.py)

Added READMEs:
  active/analysis/README.md  — script table + data flow diagram
  active/database/README.md  — database locations
  active/utilities/README.md — post-processing tools
  archive/README.md          — archive folder index"

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
