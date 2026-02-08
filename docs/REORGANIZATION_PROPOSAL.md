# KA-Search Reorganization Proposal
Generated: 2025-12-31

## Current Issues
1. Scripts scattered across root and Archive/One-off with unclear relationships
2. Multiple iterations of same script without version tracking (npz_fullscan_v2-v9)
3. No clear separation between active work and archived experiments
4. PKL files, JSON files, Excel results in separate folders but not linked to generating scripts
5. No changelog documenting what changed between versions
6. Zone.Identifier files (Windows artifacts) cluttering the repository
7. File List folder not integrated into structure

## Automatic Cleanup
The reorganization script will automatically:
- Delete all `*:Zone.Identifier` and `*.Zone.Identifier` files
- Move the `File List/` folder to `active/utilities/file_list/`

## Proposed Structure

```
KA-Search/
│
├── README.md                           # Project overview, quick start
├── CHANGELOG.md                        # High-level project changes
│
├── tools/                              # MAINTENANCE SCRIPTS (stay here permanently)
│   ├── kasearch_reorganize.py          # This reorganization script
│   ├── kasearch_maintain.py            # Ongoing maintenance (version, archive, etc.)
│   ├── find_missing_files.py           # File discovery utility
│   └── diagnose_unknown_files.py       # File analysis utility
│
├── active/                             # CURRENT WORKING SCRIPTS
│   ├── analysis/
│   │   ├── vhh_naturalness_analyzer_v3.py      # Current naturalness tool
│   │   ├── vhh_epistasis_overnight_final.py    # Current epistasis pipeline
│   │   └── visualize_epistasis.py
│   │
│   ├── database/
│   │   ├── npz_fullscan_v6_integrated.py       # Current NPZ scanner
│   │   └── build_camel_vhh_db.py               # Current DB builder
│   │
│   ├── utilities/
│   │   ├── dna_translator.py                   # DNA→AA translation (6-frame)
│   │   ├── pathlib_list.py                     # File listing generator
│   │   └── file_list/                          # Output from pathlib_list.py
│   │       └── file_list_YYYY-MM-DD.txt
│   │
│   └── config/
│       └── default_thresholds.yaml             # Configurable parameters
│
├── archive/                            # HISTORICAL VERSIONS (read-only reference)
│   │
│   ├── epistasis_pipeline/
│   │   ├── CHANGELOG.md                # What changed each version
│   │   ├── v1_20251201_overnight/
│   │   │   ├── vhh_epistasis_overnight.py
│   │   │   └── notes.md
│   │   ├── v2_20251205_position_fix/
│   │   │   ├── vhh_epistasis_overnight_v2.py
│   │   │   └── notes.md
│   │   └── v3_20251210_final/
│   │       ├── vhh_epistasis_overnight_final.py
│   │       └── notes.md
│   │
│   ├── npz_scanner/
│   │   ├── CHANGELOG.md
│   │   ├── v2_20251103_initial/
│   │   ├── v3_20251103_bugfix/
│   │   ├── v5_20251104_/
│   │   ├── v6_20251105_interactive/
│   │   ├── v7_20251106_interactive/
│   │   ├── v8_20251107_vhh/
│   │   └── v9_20251108_vhh_final/
│   │
│   ├── naturalness_analyzer/
│   │   ├── CHANGELOG.md
│   │   ├── v1_20251210_basic/
│   │   ├── v2_20251212_cdr3_fix/
│   │   └── v3_20251213_edge_cases/
│   │
│   ├── full_analysis/
│   │   ├── CHANGELOG.md
│   │   ├── run_full_analysis.py
│   │   ├── run_full_analysis_fast.py
│   │   ├── run_full_analysis_fast2.py
│   │   ├── run_full_analysis_lowmem.py
│   │   └── run_full_analysis_lowmem_v3.py
│   │
│   ├── database_builders/
│   │   ├── CHANGELOG.md
│   │   ├── build_camel_vhh_db.py
│   │   ├── build_camel_vhh_db_one_step.py
│   │   ├── build_models_only.py
│   │   └── build_pfr_cdr_models.py
│   │
│   ├── correlation_analysis/
│   │   ├── analyze_cdr_framework_advanced.py
│   │   ├── analyze_cdr_framework_correlations.py
│   │   └── vhh_global_compensation_analysis.py
│   │
│   ├── one_off/                        # Truly one-time scripts
│   │   ├── debug_npz_format.py
│   │   ├── diagnose_npz.py
│   │   ├── diagnose_oas_aligned.py
│   │   ├── diagnose_oas_data.py
│   │   ├── download_oas_camel.py
│   │   ├── process_sabdab_pdb.py
│   │   ├── process_indi_merge_final.py
│   │   ├── process_vhh_with_antigens.py
│   │   ├── shard_database.py
│   │   └── interactive_comprehensive_cdr_analysis.py
│   │
│   └── visualizations/
│       ├── visualize_correlations.py
│       ├── visualize_pfr_models.py
│       └── figures/
│           └── (move epistasis_figures here)
│
├── data/                               # INPUT DATA
│   ├── raw/
│   │   ├── sequences/
│   │   │   ├── 20251206_Mix_sequences.csv
│   │   │   ├── 20251208_PPO2_185_seq.csv
│   │   │   ├── 20251213_AVIDa_hIL6.csv
│   │   │   └── M69.fasta
│   │   │
│   │   ├── oas_paper/                  # OAS extracted data
│   │   │   ├── Compressed/
│   │   │   ├── INDI/
│   │   │   └── SRR*/
│   │   │
│   │   └── structures/                 # PDB files
│   │       └── all_nano_structures/
│   │
│   ├── databases/                      # NPZ databases
│   │   ├── production/
│   │   │   └── VHH_db_unified_v2.npz   # Current production DB
│   │   ├── shards/
│   │   │   └── (vhh_indi_ngs_*.npz)
│   │   └── legacy/
│   │       ├── VHH_db_final.npz
│   │       ├── VHH_db_unified.npz
│   │       ├── camel_vhh_clean_db.npz
│   │       └── camel_vhh_test_db.npz
│   │
│   └── external/                       # Downloaded/external data
│       ├── OAS-aligned-paper-version-20230111.tar
│       ├── PLAbDab_nano_VHH.csv.gz
│       └── PLAbDab_nano_all.csv.gz
│
├── models/                             # TRAINED MODELS & RESULTS
│   ├── epistasis/
│   │   ├── current/
│   │   │   ├── epistasis_v2_full.pkl   # Current production model
│   │   │   └── epistasis_v2_summary.json
│   │   ├── checkpoints/
│   │   │   ├── epistasis_v2_checkpoint_1_compensation.pkl
│   │   │   ├── epistasis_v2_checkpoint_2_clusters.pkl
│   │   │   ├── epistasis_v2_checkpoint_3_models.pkl
│   │   │   └── epistasis_v2_checkpoint_4_rules.pkl
│   │   └── legacy/
│   │       ├── epistasis_overnight_full.pkl
│   │       └── epistasis_overnight_checkpoints/
│   │
│   ├── correlations/
│   │   ├── correlation_results.pkl
│   │   └── correlation_summary.txt
│   │
│   └── pfr_cdr/
│       └── pfr_cdr_models.pkl
│
├── results/                            # OUTPUT RESULTS
│   ├── analysis_runs/
│   │   ├── 2025-12/
│   │   │   ├── 20251206_Mix_naturalness_FINAL.xlsx
│   │   │   ├── 20251213_AVIDa_hIL6_Results.xlsx
│   │   │   └── 20251215_PPO2_Isotope_controls.xlsx
│   │   └── 2025-11/
│   │       └── (earlier results)
│   │
│   ├── fullscan_runs/
│   │   └── (move runs_fullscan contents)
│   │
│   └── kasearch_runs/
│       └── (move runs/ contents)
│
├── docs/                               # DOCUMENTATION
│   ├── setup.md
│   ├── epistasis_methodology.md
│   ├── naturalness_scoring.md
│   └── database_schema.md
│
└── legacy/                             # UNTOUCHED ARCHIVE (safety backup)
    └── original_archive_20251231/      # Copy of current Archive/ before reorg
```

## Naming Conventions (Going Forward)

### Scripts (version only - no date in filename)
```
{purpose}_v{version}.py
Examples:
  vhh_naturalness_analyzer_v4.py
  npz_fullscan_v10.py
  cdr_extractor_v2.py
```

### Data Files (date prefix)
```
{date}_{project}_{description}.{ext}
Examples:
  20260101_IL6_binding_sequences.csv
  20260115_PPO3_translated.xlsx
```

### Results/Outputs (date prefix)
```
{date}_{input_name}_{analysis_type}.{ext}
Examples:
  20260101_IL6_binding_naturalness.xlsx
  20260115_PPO3_epistasis_scan.csv
```

### Archive Folders (version + date for historical context)
```
archive/{project}/v{version}_{date}/
Examples:
  archive/npz_scanner/v3_20251103/
  archive/epistasis_pipeline/v2_20251205/
```

## Key Changes Summary

| Current Location | New Location | Rationale |
|-----------------|--------------|-----------|
| Root `*.py` files | `active/analysis/` or `active/utilities/` | Group by function |
| `Archive/One-off Py Scripts/npz_fullscan_v*.py` | `archive/npz_scanner/v*/` | Version history |
| `Archive/One-off Py Scripts/run_full_analysis*.py` | `archive/full_analysis/` | Group iterations |
| `PKL Files/` | `models/epistasis/` | Organize by model type |
| `Excel Results/` | `results/analysis_runs/YYYY-MM/` | Organize by date |
| `Archive/NPZ Files/` | `data/databases/` | Separate prod/legacy |
| `Archive/Raw Excel Data/` | `data/raw/sequences/` | Input data |
| `Archive/runs/` | `results/kasearch_runs/` | Output results |
| `File List/` | `active/utilities/file_list/` | Utility output |
| `M69.fasta` | `data/raw/sequences/` | Reference VHH sequence |
| `*:Zone.Identifier` | (deleted) | Windows artifacts |

## Reference Files

### M69.fasta
- Single VHH sequence (119 aa)
- Header: `>M69`
- Sequence preview: `EIQLQQSGAELMKPGASVKISCKASGYTFSSYWIEWIKQRPGHGLEWIGQ...`
- Purpose: Reference/test sequence for validation

## Workflow: Using the Tools

### You Only Need TWO Files

| File | Purpose |
|------|---------|
| `maintain.py` | **THE maintenance tool** - does everything |
| `kasearch_paths.py` | Import this in your scripts for path resolution |

### Daily Usage - Just Run:
```bash
cd /home/sasenefrem/KA-Search
python maintain.py
```

That's it! It will automatically:
1. 🧹 Clean all Zone.Identifier files
2. 📁 Find unorganized files and offer to move them
3. 📝 Detect which scripts changed
4. 🔢 Auto-determine version bump (major 4.0 vs minor 3.1)
5. 📋 Generate changelog entries from your `# CHANGED`, `# FIXED`, `# NEW` comments
6. 💾 Track everything in a cache

### How Version Detection Works

The tool automatically determines if a change is major or minor:

| Change Type | Version Bump | Triggers |
|-------------|--------------|----------|
| **Major** | 3.0 → 4.0 | Functions added/removed, classes added/removed, >20% size change |
| **Minor** | 3.0 → 3.1 | Code logic changed, bug fixes, refactoring |
| **Patch** | 3.1 → 3.1 | Comments only, docstring only (just logs, no rename) |

### Auto-Changelog from Comments

Add these comments in your code and they'll be auto-extracted:
```python
# CHANGED: Updated CDR3 extraction logic
# FIXED: Edge case with short sequences
# NEW: Added support for humanized frameworks
# ADDED: Batch processing mode
# REMOVED: Deprecated legacy function
```

When you run `python maintain.py`, it finds these and uses them as changelog entries!

### Command Options

```bash
python maintain.py              # Interactive - asks before changes
python maintain.py --auto       # Auto-approve all safe changes
python maintain.py --status     # Just show what needs attention, don't change
python maintain.py --help       # Show help
```

### Initial Setup (One Time)

For the FIRST run to reorganize your existing structure:
```bash
# 1. Download maintain.py and kasearch_paths.py to KA-Search root

# 2. Also download kasearch_reorganize.py for the initial big reorg
python kasearch_reorganize.py --dry-run    # Preview
python kasearch_reorganize.py --backup     # Execute

# 3. After initial reorg, just use maintain.py going forward
python maintain.py
```

## Safeguards Against Breaking Scripts

### Option 1: Scan Before Reorganizing (Recommended First Step)
```bash
python scan_path_dependencies.py
```
This scans ALL Python files and reports:
- Which files have hardcoded paths
- Which paths will break after reorganization
- What the new paths should be

Outputs:
- `path_dependencies_report.txt` - Human-readable report
- `path_dependencies.json` - Machine-readable for auto-fixing

### Option 2: Use the Paths Module (Recommended for New Scripts)
Instead of hardcoding paths, import the central paths module:

```python
# OLD WAY (will break):
db_path = "Archive/NPZ Files/VHH_db_unified_v2.npz"
model_path = "PKL Files/epistasis_v2_full.pkl"

# NEW WAY (works before AND after reorganization):
from kasearch_paths import PATHS, get_epistasis_model, get_production_database

db_path = PATHS.database_production          # Auto-resolves
model_path = PATHS.epistasis_model           # Auto-resolves

# Or use convenience functions:
db = get_production_database()               # With fallback search
model = get_epistasis_model()                # With fallback search
```

### Option 3: Create Symlinks (Quick Fix for Old Scripts)
```bash
python kasearch_reorganize.py --backup --symlinks
```
This creates symbolic links from old paths to new paths:
- `Archive/NPZ Files/VHH_db_unified_v2.npz` → `data/databases/production/VHH_db_unified_v2.npz`
- `PKL Files/epistasis_v2_full.pkl` → `models/epistasis/current/epistasis_v2_full.pkl`
- etc.

Old scripts continue working without modification!

### Option 4: Manual Update
Use the `path_dependencies_report.txt` to manually update scripts.

## Recommended Workflow

```bash
# 1. SCAN first to see what would break
python scan_path_dependencies.py
cat path_dependencies_report.txt

# 2. If many scripts would break, use symlinks:
python kasearch_reorganize.py --backup --symlinks

# 3. If few/no scripts would break, just reorganize:
python kasearch_reorganize.py --backup

# 4. For new scripts, always use kasearch_paths.py:
from kasearch_paths import PATHS
```
