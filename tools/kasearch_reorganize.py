#!/usr/bin/env python3
"""
KA-Search Reorganization Script
================================
Moves files from current structure to proposed new organization.

IMPORTANT: Run diagnose_unknown_files.py FIRST to understand unknown files!

Usage:
    python kasearch_reorganize.py --dry-run     # Preview changes
    python kasearch_reorganize.py --execute     # Actually move files
    python kasearch_reorganize.py --backup      # Create backup first, then move

Steps:
1. Creates new directory structure
2. Copies current Archive/ to legacy/ (safety backup)
3. Moves files to new locations
4. Generates CHANGELOGs for script families
5. Creates README files
"""

import os
import re
import sys
import json
import shutil
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Tuple, Optional


# ============================================================================
# CONFIGURATION - ADJUST THESE MAPPINGS AS NEEDED
# ============================================================================

# Root-level scripts -> new locations
ROOT_SCRIPTS = {
    'vhh_epistasis_overnight_final.py': 'active/analysis/',
    'vhh_epistasis_overnight_v2.py': 'archive/epistasis_pipeline/v2_20251205/',
    'vhh_foldability_analyzer.py': 'active/analysis/',
    'visualize_epistasis.py': 'active/analysis/',
    'dna_translator.py': 'active/utilities/',
    'pathlib_list.py': 'active/utilities/',
}

# File List folder -> utilities
FILE_LIST_SCRIPTS = {
    'File List/': 'active/utilities/file_list/',  # Contains pathlib_list.py output
}

# Archive/One-off Py Scripts -> organized locations
ONEOFF_SCRIPTS = {
    # NPZ Scanner family
    'npz_fullscan_v2.py': ('archive/npz_scanner/v2_20251103/', 'v2'),
    'npz_fullscan_v3.py': ('archive/npz_scanner/v3_20251103/', 'v3'),
    'npz_fullscan_v5.py': ('archive/npz_scanner/v5_20251104/', 'v5'),
    'npz_fullscan_v6_interactive.py': ('archive/npz_scanner/v6_20251105/', 'v6'),
    'npz_fullscan_v7_interactive.py': ('archive/npz_scanner/v7_20251106/', 'v7'),
    'npz_fullscan_v8_vhh.py': ('archive/npz_scanner/v8_20251107/', 'v8'),
    'npz_fullscan_v9_vhh.py': ('archive/npz_scanner/v9_20251108/', 'v9'),
    
    # Full analysis family
    'run_full_analysis.py': ('archive/full_analysis/', 'v1'),
    'run_full_analysis_fast.py': ('archive/full_analysis/', 'v2-fast'),
    'run_full_analysis_fast2.py': ('archive/full_analysis/', 'v3-fast2'),
    'run_full_analysis_lowmem.py': ('archive/full_analysis/', 'v4-lowmem'),
    'run_full_analysis_lowmem_v3.py': ('archive/full_analysis/', 'v5-lowmem3'),
    
    # Database builders
    'build_camel_vhh_db.py': ('archive/database_builders/', 'original'),
    'build_camel_vhh_db_one_step.py': ('archive/database_builders/', 'one-step'),
    'build_models_only.py': ('archive/database_builders/', 'models-only'),
    'build_pfr_cdr_models.py': ('archive/database_builders/', 'pfr-models'),
    
    # Correlation analysis
    'analyze_cdr_framework_advanced.py': ('archive/correlation_analysis/', 'advanced'),
    'analyze_cdr_framework_correlations.py': ('archive/correlation_analysis/', 'basic'),
    'vhh_global_compensation_analysis.py': ('archive/correlation_analysis/', 'compensation'),
    
    # Epistasis
    'vhh_full_epistasis_overnight.py': ('archive/epistasis_pipeline/v1_20251201/', 'v1'),
    
    # Diagnostics & One-off (stay in one_off)
    'debug_npz_format.py': ('archive/one_off/', None),
    'diagnose_npz.py': ('archive/one_off/', None),
    'diagnose_oas_aligned.py': ('archive/one_off/', None),
    'diagnose_oas_data.py': ('archive/one_off/', None),
    'download_oas_camel.py': ('archive/one_off/', None),
    'process_sabdab_pdb.py': ('archive/one_off/', None),
    'process_indi_merge_final.py': ('archive/one_off/', None),
    'process_vhh_with_antigens.py': ('archive/one_off/', None),
    'process_camel_vhh_pipeline.py': ('archive/one_off/', None),
    'shard_database.py': ('archive/one_off/', None),
    'interactive_comprehensive_cdr_analysis.py': ('archive/one_off/', None),
    
    # Framework scoring/optimization (naturalness-related)
    'score_framework.py': ('archive/naturalness_analyzer/', 'score'),
    'vhh_framework_optimizer.py': ('archive/naturalness_analyzer/', 'optimizer'),
    
    # Visualization scripts
    'visualize_correlations.py': ('archive/visualizations/', None),
    'visualize_pfr_models.py': ('archive/visualizations/', None),
}

# Data file mappings
DATA_MAPPINGS = {
    'Archive/NPZ Files/VHH_db_unified_v2.npz': 'data/databases/production/',
    'Archive/NPZ Files/VHH_db_final.npz': 'data/databases/legacy/',
    'Archive/NPZ Files/VHH_db_unified.npz': 'data/databases/legacy/',
    'Archive/NPZ Files/camel_vhh_clean_db.npz': 'data/databases/legacy/',
    'Archive/NPZ Files/camel_vhh_test_db.npz': 'data/databases/legacy/',
    
    'Archive/Raw Excel Data/': 'data/raw/sequences/',
    'Archive/VHH_shards/': 'data/databases/shards/',
    'Archive/Databases/': 'data/external/',
    
    'Archive/extracted/oas-paper/INDI/': 'data/raw/oas_paper/INDI/',
    'Archive/extracted/oas-paper/Compressed/': 'data/raw/oas_paper/compressed/',
    'Archive/extracted/oas-paper/all_nano_structures/': 'data/raw/structures/',
    
    'M69.fasta': 'data/raw/sequences/',
}

# Model/PKL mappings
MODEL_MAPPINGS = {
    'PKL Files/epistasis_v2_full.pkl': 'models/epistasis/current/',
    'PKL Files/epistasis_overnight_full.pkl': 'models/epistasis/legacy/',
    'PKL Files/epistasis_v2_checkpoint_1_compensation.pkl': 'models/epistasis/checkpoints/',
    'PKL Files/epistasis_v2_checkpoint_2_clusters.pkl': 'models/epistasis/checkpoints/',
    'PKL Files/epistasis_v2_checkpoint_3_models.pkl': 'models/epistasis/checkpoints/',
    'PKL Files/epistasis_v2_checkpoint_4_rules.pkl': 'models/epistasis/checkpoints/',
    'PKL Files/epistasis_overnight_checkpoint_1_compensation.pkl': 'models/epistasis/legacy_checkpoints/',
    'PKL Files/epistasis_overnight_checkpoint_2_clusters.pkl': 'models/epistasis/legacy_checkpoints/',
    'PKL Files/epistasis_overnight_checkpoint_3_models.pkl': 'models/epistasis/legacy_checkpoints/',
    'PKL Files/epistasis_overnight_checkpoint_4_rules.pkl': 'models/epistasis/legacy_checkpoints/',
    'PKL Files/compensation_results.pkl': 'models/correlations/',
    
    'Model Results/1_correlations/': 'models/correlations/',
    'Model Results/2_models/': 'models/pfr_cdr/',
    
    'JSON_Files/': 'models/epistasis/summaries/',
}

# Results mappings
RESULTS_MAPPINGS = {
    'Excel Results/': 'results/analysis_runs/2025-12/',
    'Archive/runs/': 'results/kasearch_runs/',
    'Archive/runs_fullscan/': 'results/fullscan_runs/',
    'Model Results/20251119_camel/': 'results/npz_scans/2025-11-19/',
    'Model Results/20251120_camel/': 'results/npz_scans/2025-11-20/',
    'Model Results/20251121_camel/': 'results/npz_scans/2025-11-21/',
    'Model Results/20251124_camel/': 'results/npz_scans/2025-11-24/',
}


# ============================================================================
# DIRECTORY STRUCTURE
# ============================================================================

NEW_DIRECTORIES = [
    'active/analysis',
    'active/database',
    'active/utilities',
    'active/utilities/file_list',
    'active/config',
    
    'tools',  # Maintenance scripts live here permanently
    
    'archive/epistasis_pipeline',
    'archive/npz_scanner',
    'archive/naturalness_analyzer',
    'archive/full_analysis',
    'archive/database_builders',
    'archive/correlation_analysis',
    'archive/one_off',
    'archive/visualizations',
    
    'data/raw/sequences',
    'data/raw/oas_paper',
    'data/raw/structures',
    'data/databases/production',
    'data/databases/shards',
    'data/databases/legacy',
    'data/external',
    
    'models/epistasis/current',
    'models/epistasis/checkpoints',
    'models/epistasis/legacy',
    'models/epistasis/legacy_checkpoints',
    'models/epistasis/summaries',
    'models/correlations',
    'models/pfr_cdr',
    
    'results/analysis_runs/2025-12',
    'results/analysis_runs/2025-11',
    'results/kasearch_runs',
    'results/fullscan_runs',
    'results/npz_scans',
    
    'docs',
    'legacy/original_archive_20251231',
]

# Scripts that should STAY in tools/ and not be moved
TOOLS_SCRIPTS = [
    'maintain.py',
    'kasearch_reorganize.py',
    'kasearch_maintain.py',
    'kasearch_paths.py',
    'find_missing_files.py',
    'diagnose_unknown_files.py',
    'scan_path_dependencies.py',
]


# ============================================================================
# FUNCTIONS
# ============================================================================

def create_directories(root: Path, dry_run: bool = True):
    """Create new directory structure."""
    print("\n📁 Creating directory structure...")
    
    for dir_path in NEW_DIRECTORIES:
        full_path = root / dir_path
        if dry_run:
            print(f"  [DRY-RUN] Would create: {dir_path}")
        else:
            full_path.mkdir(parents=True, exist_ok=True)
            print(f"  ✅ Created: {dir_path}")


def backup_archive(root: Path, dry_run: bool = True):
    """Create safety backup of current Archive."""
    src = root / 'Archive'
    dst = root / 'legacy' / 'original_archive_20251231'
    
    if not src.exists():
        print("  ⚠️  No Archive/ directory found")
        return
    
    if dry_run:
        print(f"  [DRY-RUN] Would backup Archive/ to legacy/original_archive_20251231/")
    else:
        if dst.exists():
            print(f"  ⚠️  Backup already exists: {dst}")
            return
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copytree(src, dst)
        print(f"  ✅ Backed up Archive/ to {dst}")


def move_file(src: Path, dst: Path, dry_run: bool = True) -> bool:
    """Move a single file."""
    if not src.exists():
        return False
    
    if dry_run:
        print(f"  [DRY-RUN] {src.name} -> {dst}")
        return True
    else:
        dst.parent.mkdir(parents=True, exist_ok=True)
        if dst.exists():
            print(f"  ⚠️  Target exists, skipping: {dst}")
            return False
        shutil.copy2(src, dst)
        print(f"  ✅ {src.name} -> {dst}")
        return True


def move_directory(src: Path, dst: Path, dry_run: bool = True) -> bool:
    """Move/copy a directory."""
    if not src.exists():
        return False
    
    if dry_run:
        print(f"  [DRY-RUN] {src} -> {dst}")
        return True
    else:
        dst.parent.mkdir(parents=True, exist_ok=True)
        if dst.exists():
            print(f"  ⚠️  Target exists, skipping: {dst}")
            return False
        shutil.copytree(src, dst)
        print(f"  ✅ {src} -> {dst}")
        return True


def clean_zone_identifiers(root: Path, dry_run: bool = True) -> int:
    """Delete all Zone.Identifier files (Windows artifact)."""
    print("\n🧹 Cleaning Zone.Identifier files...")
    
    zone_files = []
    for f in root.rglob('*'):
        if 'Zone.Identifier' in f.name:
            zone_files.append(f)
    
    if not zone_files:
        print("  No Zone.Identifier files found.")
        return 0
    
    deleted = 0
    for zf in zone_files:
        if dry_run:
            print(f"  [DRY-RUN] Would delete: {zf.relative_to(root)}")
            deleted += 1
        else:
            try:
                zf.unlink()
                deleted += 1
            except Exception as e:
                print(f"  ⚠️ Failed: {zf}: {e}")
    
    print(f"  {'Would delete' if dry_run else 'Deleted'}: {deleted} Zone.Identifier files")
    return deleted


def setup_tools_directory(root: Path, dry_run: bool = True):
    """Move/copy maintenance scripts to tools/ directory."""
    print("\n🔧 Setting up tools/ directory...")
    
    tools_dir = root / 'tools'
    
    for script in TOOLS_SCRIPTS:
        # Check if script exists in root
        src = root / script
        dst = tools_dir / script
        
        if src.exists() and not dst.exists():
            if dry_run:
                print(f"  [DRY-RUN] Would move: {script} -> tools/")
            else:
                tools_dir.mkdir(parents=True, exist_ok=True)
                shutil.copy2(src, dst)
                print(f"  ✅ Copied: {script} -> tools/")
        elif dst.exists():
            print(f"  ✓ Already in tools/: {script}")
        elif not src.exists():
            # Check if it's in tools/ already
            if (tools_dir / script).exists():
                print(f"  ✓ Already in tools/: {script}")
            else:
                print(f"  ℹ️ Not found (will be added when you download): {script}")


def add_version_to_filename(filepath: Path, version: int) -> str:
    """Generate versioned filename (no date - just version number)."""
    base = filepath.stem
    # Remove existing version if present
    base = re.sub(r'_v\d+(_\d{8})?$', '', base)
    
    return f"{base}_v{version}.py"


def create_symlinks(root: Path, dry_run: bool = True):
    """
    Create symbolic links from old paths to new paths.
    This allows old scripts to continue working without modification.
    """
    print("\n🔗 Creating backward-compatibility symlinks...")
    
    # Symlink mapping: old_path -> new_path
    symlink_map = {
        'Archive/NPZ Files/VHH_db_unified_v2.npz': 'data/databases/production/VHH_db_unified_v2.npz',
        'PKL Files/epistasis_v2_full.pkl': 'models/epistasis/current/epistasis_v2_full.pkl',
        'PKL Files/epistasis_overnight_full.pkl': 'models/epistasis/legacy/epistasis_overnight_full.pkl',
        'Model Results/1_correlations': 'models/correlations',
        'Model Results/2_models': 'models/pfr_cdr',
        'Archive/VHH_shards': 'data/databases/shards',
        'Archive/Raw Excel Data': 'data/raw/sequences',
        'Excel Results': 'results/analysis_runs',
        'JSON_Files': 'models/epistasis/summaries',
    }
    
    created = 0
    for old_rel, new_rel in symlink_map.items():
        old_path = root / old_rel
        new_path = root / new_rel
        
        # Only create symlink if:
        # 1. New path exists
        # 2. Old path doesn't exist (or is the same as new)
        if not new_path.exists():
            continue
        
        if old_path.exists() and not old_path.is_symlink():
            # Old path still exists as real file/dir - skip
            continue
        
        if dry_run:
            print(f"  [DRY-RUN] Would symlink: {old_rel} -> {new_rel}")
            created += 1
        else:
            try:
                # Create parent directories if needed
                old_path.parent.mkdir(parents=True, exist_ok=True)
                
                # Remove existing symlink if present
                if old_path.is_symlink():
                    old_path.unlink()
                
                # Create relative symlink
                rel_target = os.path.relpath(new_path, old_path.parent)
                old_path.symlink_to(rel_target)
                print(f"  ✅ Symlinked: {old_rel} -> {new_rel}")
                created += 1
            except Exception as e:
                print(f"  ⚠️ Failed to create symlink {old_rel}: {e}")
    
    print(f"  {'Would create' if dry_run else 'Created'}: {created} symlinks")
    return created


def cleanup_old_folders(root: Path):
    """Delete old folders after successful reorganization."""
    print("\n🗑️  Cleaning up old folders...")
    
    # Folders to delete (only if they're now empty or files have been moved)
    old_folders = [
        'Archive/One-off Py Scripts',
        'Archive/NPZ Files', 
        'Archive/Raw Excel Data',
        'Archive/VHH_shards',
        'Archive/runs',
        'Archive/runs_fullscan',
        'Archive/epistasis_figures',
        'Archive/Zip Files',
        'PKL Files',
        'Excel Results',
        'JSON_Files',
        'Model Results',
        'File List',
    ]
    
    deleted = 0
    for folder in old_folders:
        folder_path = root / folder
        if folder_path.exists():
            try:
                shutil.rmtree(folder_path)
                print(f"  ✅ Deleted: {folder}/")
                deleted += 1
            except Exception as e:
                print(f"  ⚠️ Could not delete {folder}/: {e}")
    
    # Check if Archive is now empty (or only has legacy backup)
    archive_path = root / 'Archive'
    if archive_path.exists():
        remaining = list(archive_path.iterdir())
        # Filter out hidden files and legacy
        remaining = [r for r in remaining if not r.name.startswith('.')]
        
        if not remaining:
            try:
                shutil.rmtree(archive_path)
                print(f"  ✅ Deleted empty Archive/")
                deleted += 1
            except:
                pass
        else:
            print(f"\n  ℹ️ Archive/ still has contents: {[r.name for r in remaining[:5]]}")
            print(f"     Review and delete manually if no longer needed.")
    
    print(f"\n  Deleted {deleted} old folders")


def reorganize_scripts(root: Path, dry_run: bool = True, add_versions: bool = False):
    """Move scripts to new locations."""
    print("\n📜 Reorganizing scripts...")
    
    # Root-level scripts (skip tools scripts)
    print("\n  Root scripts:")
    for script, dest in ROOT_SCRIPTS.items():
        if script in TOOLS_SCRIPTS:
            print(f"  ⏭️ Skipping (tools script): {script}")
            continue
            
        src = root / script
        if not src.exists():
            continue
            
        dst = root / dest / script
        move_file(src, dst, dry_run)
    
    # Archive/One-off scripts
    print("\n  Archive/One-off scripts:")
    for script, (dest, version_tag) in ONEOFF_SCRIPTS.items():
        src = root / 'Archive' / 'One-off Py Scripts' / script
        if not src.exists():
            continue
        
        # Optionally add version to filename
        if add_versions and version_tag:
            # Extract version number from tag
            version_match = re.search(r'v?(\d+)', version_tag)
            if version_match:
                version_num = int(version_match.group(1))
                new_name = add_version_to_filename(src, version_num)
                dst = root / dest / new_name
            else:
                dst = root / dest / script
        else:
            dst = root / dest / script
            
        move_file(src, dst, dry_run)


def reorganize_data(root: Path, dry_run: bool = True):
    """Move data files to new locations."""
    print("\n📊 Reorganizing data files...")
    
    for src_pattern, dest in DATA_MAPPINGS.items():
        src = root / src_pattern
        
        if src_pattern.endswith('/'):
            # It's a directory
            if src.exists():
                dst = root / dest
                move_directory(src, dst, dry_run)
        else:
            # It's a file
            dst = root / dest / src.name
            move_file(src, dst, dry_run)


def reorganize_models(root: Path, dry_run: bool = True):
    """Move model files to new locations."""
    print("\n🧠 Reorganizing model files...")
    
    for src_pattern, dest in MODEL_MAPPINGS.items():
        src = root / src_pattern
        
        if src_pattern.endswith('/'):
            # It's a directory
            if src.exists():
                dst = root / dest
                move_directory(src, dst, dry_run)
        else:
            # It's a file
            dst = root / dest / src.name
            move_file(src, dst, dry_run)


def reorganize_results(root: Path, dry_run: bool = True):
    """Move results files to new locations."""
    print("\n📈 Reorganizing results...")
    
    for src_pattern, dest in RESULTS_MAPPINGS.items():
        src = root / src_pattern
        
        if src_pattern.endswith('/'):
            if src.exists():
                dst = root / dest
                move_directory(src, dst, dry_run)
        else:
            dst = root / dest / src.name
            move_file(src, dst, dry_run)


def generate_changelogs(root: Path, dry_run: bool = True):
    """Generate CHANGELOG.md files for script families."""
    print("\n📝 Generating changelogs...")
    
    changelogs = {
        'archive/npz_scanner': """# NPZ Scanner Changelog

## v9 (2025-11-08)
- VHH-specific optimizations
- File: `npz_fullscan_v9_vhh.py`

## v8 (2025-11-07)
- Added VHH support
- File: `npz_fullscan_v8_vhh.py`

## v7 (2025-11-06)
- Interactive mode improvements
- File: `npz_fullscan_v7_interactive.py`

## v6 (2025-11-05)
- Added interactive mode
- File: `npz_fullscan_v6_interactive.py`

## v5 (2025-11-04)
- Bug fixes and optimizations
- File: `npz_fullscan_v5.py`

## v3 (2025-11-03)
- Initial improvements over v2
- File: `npz_fullscan_v3.py`

## v2 (2025-11-03)
- First iteration with IMGT CDR extraction
- File: `npz_fullscan_v2.py`
""",
        'archive/epistasis_pipeline': """# Epistasis Pipeline Changelog

## v3 (2025-12-10) - CURRENT
- Final production version
- Fixed position mapping issues
- Comprehensive output format
- File: `vhh_epistasis_overnight_final.py`

## v2 (2025-12-05)
- Position mapping corrections
- Added Vernier cluster analysis
- File: `vhh_epistasis_overnight_v2.py`

## v1 (2025-12-01)
- Initial overnight analysis script
- File: `vhh_full_epistasis_overnight.py`
""",
        'archive/full_analysis': """# Full Analysis Changelog

## v5 (2025-11-15) - lowmem_v3
- Further memory optimizations
- File: `run_full_analysis_lowmem_v3.py`

## v4 (2025-11-12) - lowmem
- Low memory version for large datasets
- File: `run_full_analysis_lowmem.py`

## v3 (2025-11-10) - fast2
- Additional speed optimizations
- File: `run_full_analysis_fast2.py`

## v2 (2025-11-08) - fast
- Speed optimizations with vectorization
- File: `run_full_analysis_fast.py`

## v1 (2025-11-05) - original
- Initial version
- File: `run_full_analysis.py`
""",
        'archive/naturalness_analyzer': """# Naturalness Analyzer Changelog

## v3 (2025-12-13) - CURRENT
- Fixed CDR3 extraction edge cases
- 100% extraction success on test set
- File: `vhh_naturalness_analyzer_v3.py`

## v2 (2025-12-12)
- CDR3 extraction bug fixes
- File: `vhh_naturalness_analyzer_v2.py`

## v1 (2025-12-10)
- Initial naturalness scoring implementation
- File: `vhh_naturalness_analyzer_v1.py`
""",
    }
    
    for path, content in changelogs.items():
        changelog_path = root / path / 'CHANGELOG.md'
        if dry_run:
            print(f"  [DRY-RUN] Would create: {changelog_path}")
        else:
            changelog_path.parent.mkdir(parents=True, exist_ok=True)
            changelog_path.write_text(content)
            print(f"  ✅ Created: {changelog_path}")


def generate_readme(root: Path, dry_run: bool = True):
    """Generate main README.md."""
    print("\n📖 Generating README...")
    
    readme_content = """# KA-Search VHH Analysis

Nanobody/VHH sequence analysis using KA-Search and epistasis modeling.

## Directory Structure

```
KA-Search/
├── active/          # Current working scripts
├── archive/         # Historical versions (organized by project)
├── data/            # Input data (sequences, databases, structures)
├── models/          # Trained models and checkpoints
├── results/         # Analysis outputs
├── docs/            # Documentation
└── legacy/          # Backup of original Archive/
```

## Quick Start

### Run Naturalness Analysis
```bash
python active/analysis/vhh_naturalness_analyzer_v3.py \\
    --input data/raw/sequences/your_sequences.csv \\
    --epistasis models/epistasis/current/epistasis_v2_full.pkl
```

### Run NPZ Scan
```bash
python active/database/npz_fullscan_v6_integrated.py \\
    --input your_sequences.csv \\
    --database data/databases/production/VHH_db_unified_v2.npz
```

## Key Models

- **Epistasis Model**: `models/epistasis/current/epistasis_v2_full.pkl`
  - Vernier cluster patterns from 11.5M camelid sequences
  - CDR3 length/charge/aromatic statistics per cluster
  
- **Correlation Results**: `models/correlations/correlation_results.pkl`
  - Position-specific CDR-framework coupling rules
  - Triplet co-occurrence statistics

## Documentation

See `docs/` for:
- `epistasis_methodology.md` - How the epistasis model works
- `naturalness_scoring.md` - Naturalness score interpretation
- `database_schema.md` - NPZ database format

## Maintenance

Use `kasearch_maintain.py` to:
- Check naming conventions: `python kasearch_maintain.py check`
- Create new scripts: `python kasearch_maintain.py new my_script`
- Archive old versions: `python kasearch_maintain.py archive script_name`

Reorganized: """ + datetime.now().strftime('%Y-%m-%d')
    
    readme_path = root / 'README.md'
    if dry_run:
        print(f"  [DRY-RUN] Would create: {readme_path}")
    else:
        readme_path.write_text(readme_content)
        print(f"  ✅ Created: {readme_path}")


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        print("\nOptions:")
        print("  --dry-run      Preview changes without moving files")
        print("  --execute      Actually reorganize files")
        print("  --backup       Create backup first, then reorganize")
        print("  --add-versions Also rename files with version numbers (e.g., script_v2.py)")
        print("  --symlinks     Create symlinks from old paths to new (for backward compatibility)")
        print("  --cleanup      Delete old folders (run after successful reorganization)")
        sys.exit(1)
    
    mode = sys.argv[1].lower()
    
    # Check if just doing cleanup
    if mode == '--cleanup':
        root = Path.cwd()
        if not (root / 'active').exists():
            print("Error: No 'active/' directory found. Run reorganization first.")
            sys.exit(1)
        print("=" * 70)
        print("CLEANUP MODE - Deleting old folders")
        print("=" * 70)
        confirm = input("\n⚠️  This will delete old folders (Archive/, PKL Files/, etc.). Continue? [y/n]: ")
        if confirm.lower() == 'y':
            cleanup_old_folders(root)
            print("\n✅ Cleanup complete!")
        else:
            print("Aborted.")
        sys.exit(0)
    
    dry_run = mode == '--dry-run'
    do_backup = mode == '--backup'
    add_versions = '--add-versions' in sys.argv
    create_symlinks_flag = '--symlinks' in sys.argv
    do_cleanup = '--cleanup' in sys.argv  # Can combine with --backup or --execute
    
    root = Path.cwd()
    
    # Verify we're in KA-Search
    if not (root / 'Archive').exists():
        print("Error: No Archive/ directory found.")
        print("Please run this from the KA-Search root directory.")
        sys.exit(1)
    
    print("=" * 70)
    print("KA-SEARCH REORGANIZATION")
    print("=" * 70)
    print(f"Root: {root}")
    print(f"Mode: {'DRY-RUN (preview only)' if dry_run else 'EXECUTE (will move files)'}")
    print(f"Add versions to filenames: {add_versions}")
    print(f"Create backward-compat symlinks: {create_symlinks_flag}")
    print("=" * 70)
    
    if not dry_run:
        confirm = input("\n⚠️  This will reorganize your files. Continue? [y/n]: ")
        if confirm.lower() != 'y':
            print("Aborted.")
            sys.exit(0)
    
    # Step 1: Create directories
    create_directories(root, dry_run)
    
    # Step 2: Setup tools directory (maintenance scripts stay here)
    setup_tools_directory(root, dry_run)
    
    # Step 3: Clean Zone.Identifier files
    clean_zone_identifiers(root, dry_run)
    
    # Step 4: Backup (if requested)
    if do_backup or not dry_run:
        print("\n💾 Creating backup...")
        backup_archive(root, dry_run)
    
    # Step 5: Reorganize
    reorganize_scripts(root, dry_run, add_versions)
    reorganize_data(root, dry_run)
    reorganize_models(root, dry_run)
    reorganize_results(root, dry_run)
    
    # Step 6: Move File List folder
    print("\n📋 Moving File List folder...")
    file_list_src = root / 'File List'
    file_list_dst = root / 'active' / 'utilities' / 'file_list'
    if file_list_src.exists():
        move_directory(file_list_src, file_list_dst, dry_run)
    
    # Step 7: Create symlinks (if requested)
    if create_symlinks_flag:
        create_symlinks(root, dry_run)
    
    # Step 8: Generate documentation
    generate_changelogs(root, dry_run)
    generate_readme(root, dry_run)
    
    # Step 9: Copy kasearch_paths.py to tools
    paths_module_src = root / 'kasearch_paths.py'
    paths_module_dst = root / 'tools' / 'kasearch_paths.py'
    if paths_module_src.exists() and not paths_module_dst.exists():
        if dry_run:
            print(f"\n  [DRY-RUN] Would copy kasearch_paths.py to tools/")
        else:
            shutil.copy2(paths_module_src, paths_module_dst)
            print(f"\n  ✅ Copied kasearch_paths.py to tools/")
    
    # Step 10: Cleanup old folders (if requested)
    if do_cleanup and not dry_run:
        cleanup_old_folders(root)
    elif do_cleanup and dry_run:
        print("\n🗑️  [DRY-RUN] Would delete old folders:")
        old_folders = [
            'Archive/One-off Py Scripts',
            'Archive/NPZ Files', 
            'Archive/Raw Excel Data',
            'Archive/VHH_shards',
            'Archive/runs',
            'Archive/runs_fullscan',
            'PKL Files',
            'Excel Results',
            'JSON_Files',
            'Model Results',
            'File List',
        ]
        for folder in old_folders:
            if (root / folder).exists():
                print(f"    {folder}/")
    
    print("\n" + "=" * 70)
    if dry_run:
        print("DRY-RUN COMPLETE - No files were moved")
        print("Run with --execute or --backup to actually reorganize")
        print("\nOptions:")
        print("  --add-versions  Rename scripts with version numbers")
        print("  --symlinks      Create symlinks for backward compatibility")
        print("  --cleanup       Delete old folders after reorganization")
    else:
        print("REORGANIZATION COMPLETE")
        
        if not do_cleanup:
            print("\n⚠️  Old folders still exist (safe mode).")
            print("   To delete them, run again with --cleanup")
            print("   Or delete manually after verifying the new structure.")
        
        print("\n" + "=" * 70)
        print("GOING FORWARD - Just use maintain.py!")
        print("=" * 70)
        print("""
    python maintain.py          # Auto-maintains everything:
                                #   - Cleans Zone.Identifier files
                                #   - Organizes new files  
                                #   - Auto-versions changed scripts
                                #   - Updates changelog

    python maintain.py --status # Just check status
    python maintain.py --auto   # Auto-approve all changes
""")
        print("For imports in your scripts:")
        print("    from kasearch_paths import PATHS")
    print("=" * 70)


if __name__ == '__main__':
    main()
