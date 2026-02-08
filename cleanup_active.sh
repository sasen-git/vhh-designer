#!/usr/bin/env bash
set -euo pipefail

# ============================================================================
# CLEAN UP active/ — keep only latest versions, archive the rest
#
# Usage:
#   bash cleanup_active.sh              # dry run
#   bash cleanup_active.sh --execute    # do it + commit
# ============================================================================

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
echo "  1. DESIGNER — keep v9_1, archive everything else"
echo "========================================================================"

run_cmd mkdir -p "'$AR/designer'"

# Move the existing subfolder (v2-v5_1, v7) into archive
if [ -d "$A/vhh_designer" ]; then
    for f in "$A"/vhh_designer/*.py; do
        [ -f "$f" ] && run_cmd mv "'$f'" "'$AR/designer/'"
    done
    run_cmd rmdir "'$A/vhh_designer'" 2>/dev/null || true
fi

# Move all loose designer iterations EXCEPT v9_1
for f in "$A"/vhh_designer_v*.py; do
    [ -f "$f" ] || continue
    base=$(basename "$f")
    if [[ "$base" == "vhh_designer_v9_1.py" ]]; then
        echo "  [KEEP] $base"
    else
        run_cmd mv "'$f'" "'$AR/designer/'"
    fi
done

echo ""

# ============================================================================
echo "========================================================================"
echo "  2. ANALYSIS UNIFIED — keep v7.7, archive v7-v7.6"
echo "========================================================================"

run_cmd mkdir -p "'$AR/analysis_unified'"

for f in "$A"/vhh_analysis_unified_v*.py; do
    [ -f "$f" ] || continue
    base=$(basename "$f")
    if [[ "$base" == "vhh_analysis_unified_v7.7.py" ]]; then
        echo "  [KEEP] $base"
    else
        run_cmd mv "'$f'" "'$AR/analysis_unified/'"
    fi
done

echo ""

# ============================================================================
echo "========================================================================"
echo "  3. COMPENSATION — keep v6_csvcols, archive the rest"
echo "========================================================================"

run_cmd mkdir -p "'$AR/compensation'"

# Move existing subfolder contents
if [ -d "$A/vhh_compensation" ]; then
    for f in "$A"/vhh_compensation/*.py; do
        [ -f "$f" ] && run_cmd mv "'$f'" "'$AR/compensation/'"
    done
    run_cmd rmdir "'$A/vhh_compensation'" 2>/dev/null || true
fi

# Move loose compensation files except v6
for f in "$A"/vhh_compensation_imgt_v*.py; do
    [ -f "$f" ] || continue
    base=$(basename "$f")
    if [[ "$base" == "vhh_compensation_imgt_v6_csvcols.py" ]]; then
        echo "  [KEEP] $base"
    else
        run_cmd mv "'$f'" "'$AR/compensation/'"
    fi
done

echo ""

# ============================================================================
echo "========================================================================"
echo "  4. ANNOTATE SHARDS — keep v8, archive v1-v7"
echo "========================================================================"

run_cmd mkdir -p "'$AR/annotate_shards'"

# Move existing subfolder contents
if [ -d "$A/annotate_all_shards" ]; then
    for f in "$A"/annotate_all_shards/*.py; do
        [ -f "$f" ] && run_cmd mv "'$f'" "'$AR/annotate_shards/'"
    done
    run_cmd rmdir "'$A/annotate_all_shards'" 2>/dev/null || true
fi

# v8 stays (it's already loose in active/analysis/)
echo "  [KEEP] annotate_all_shards_v8.py"

echo ""

# ============================================================================
echo "========================================================================"
echo "  5. NATURALNESS ANALYZER — keep v4, archive the rest"
echo "========================================================================"

# Existing archive already has v1, v2 — add v3 and foldability
if [ -d "$A/vhh_naturalness_analyzer" ]; then
    for f in "$A"/vhh_naturalness_analyzer/*.py; do
        [ -f "$f" ] && run_cmd mv "'$f'" "'$AR/naturalness_analyzer/'"
    done
    run_cmd rmdir "'$A/vhh_naturalness_analyzer'" 2>/dev/null || true
fi

echo "  [KEEP] vhh_naturalness_analyzer_v4.py"

echo ""

# ============================================================================
echo "========================================================================"
echo "  6. EPISTASIS — keep v3_imgt, archive older subfolder contents"
echo "========================================================================"

# The vhh_epistasis/ subfolder has overnight_final, overnight_v2, visualize
# These are the legacy (pre-IMGT) pipeline — move to existing archive folder
if [ -d "$A/vhh_epistasis" ]; then
    for f in "$A"/vhh_epistasis/*.py; do
        [ -f "$f" ] || continue
        base=$(basename "$f")
        # Check if already in archive
        if [ ! -f "$AR/epistasis_pipeline/$base" ]; then
            run_cmd mv "'$f'" "'$AR/epistasis_pipeline/'"
        else
            echo "  [SKIP] $base (already in archive)"
            run_cmd rm "'$f'"
        fi
    done
    run_cmd rmdir "'$A/vhh_epistasis'" 2>/dev/null || true
fi

# vhh_epistasis_to_imgt/ — older IMGT conversion iterations
run_cmd mkdir -p "'$AR/epistasis_to_imgt'"
if [ -d "$A/vhh_epistasis_to_imgt" ]; then
    for f in "$A"/vhh_epistasis_to_imgt/*.py; do
        [ -f "$f" ] && run_cmd mv "'$f'" "'$AR/epistasis_to_imgt/'"
    done
    run_cmd rmdir "'$A/vhh_epistasis_to_imgt'" 2>/dev/null || true
fi

echo "  [KEEP] vhh_epistasis_v3_imgt.py"

echo ""

# ============================================================================
echo "========================================================================"
echo "  7. AUGMENT IMGT — keep v2, archive v1"
echo "========================================================================"

run_cmd mkdir -p "'$AR/augment_imgt'"

for f in "$A"/augment_imgt_positions_from_csvs_v*.py; do
    [ -f "$f" ] || continue
    base=$(basename "$f")
    if [[ "$base" == "augment_imgt_positions_from_csvs_v2.py" ]]; then
        echo "  [KEEP] $base"
    else
        run_cmd mv "'$f'" "'$AR/augment_imgt/'"
    fi
done

echo ""

# ============================================================================
echo "========================================================================"
echo "  8. CHANGELOGS — move to docs/"
echo "========================================================================"

for f in "$AR"/*/CHANGELOG.md; do
    [ -f "$f" ] || continue
    folder=$(basename "$(dirname "$f")")
    run_cmd mv "'$f'" "'$ROOT/docs/CHANGELOG_${folder}.md'"
done

echo ""

# ============================================================================
echo "========================================================================"
echo "  SUMMARY — what stays in active/analysis/"
echo "========================================================================"

echo "  Current (latest) scripts:"
echo "    vhh_designer_v9_1.py"
echo "    vhh_analysis_unified_v7.7.py"
echo "    vhh_compensation_imgt_v6_csvcols.py"
echo "    vhh_naturalness_analyzer_v4.py"
echo "    vhh_epistasis_v3_imgt.py"
echo "    annotate_all_shards_v8.py"
echo "    augment_imgt_positions_from_csvs_v2.py"
echo "    analyze_cdr_framework_correlations_imgt_csv_v3.py"
echo "    charge_balance_vhh_windows.py"
echo "    convert_epistasis_to_designer.py"
echo "    extract_script_headers.py"
echo "    patches.py"
echo "    run_with_positions.py"
echo ""
echo "  active/database/"
echo "    npz_fullscan_v6_integrated.py"
echo ""
echo "  active/utilities/"
echo "    dna_translator.py, pathlib_list_v2.py, pull_cdrs.py, scfv_anarci.py"
echo "    csv2msa_antpack_v5.py, csv2msa_antpack_v5_imgt_v2.py, csv2msa_antpack_v6.py"
echo "    align_vs_lead_clear3_antpack_legend_v7.py (+ v7_imgt_v2, v8, v8_imgt)"

echo ""

# ============================================================================
echo "========================================================================"
echo "  9. GIT COMMIT"
echo "========================================================================"

if $DRY_RUN; then
    echo "  [DRY] Would run:"
    echo "    cd $ROOT && git add -A && git commit && git push"
else
    cd "$ROOT"
    git add -A

    echo ""
    echo "--- Changes ---"
    git status --short | head -40
    TOTAL=$(git status --short | wc -l)
    echo "  ... $TOTAL files changed"
    echo ""

    read -p "Commit and push? [y/N] " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        git commit -m "Clean up active/: keep latest versions, archive iterations

Kept in active/analysis/ (current production):
  vhh_designer_v9_1.py
  vhh_analysis_unified_v7.7.py
  vhh_compensation_imgt_v6_csvcols.py
  vhh_naturalness_analyzer_v4.py
  vhh_epistasis_v3_imgt.py
  annotate_all_shards_v8.py
  augment_imgt_positions_from_csvs_v2.py
  + standalone tools (charge_balance, convert_epistasis, patches, etc.)

Archived iterations to:
  archive/designer/           v2-v9_0 (~40 versions)
  archive/analysis_unified/   v7-v7.6
  archive/compensation/       v1-v5
  archive/annotate_shards/    v1-v7
  archive/epistasis_to_imgt/  v1-v3
  archive/augment_imgt/       v1"

        git push
        echo ""
        echo "✅ Pushed!"
    else
        echo "  Aborted."
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
