#!/usr/bin/env python3
"""
KA-Search Paths Configuration
==============================
Central path resolution for all KA-Search scripts.

IMPORT THIS MODULE instead of hardcoding paths!

Usage in your scripts:
    from kasearch_paths import PATHS, get_path, find_file
    
    # Get database path (works before AND after reorganization)
    db_path = PATHS.database_production
    
    # Or use flexible finder
    epistasis_model = find_file('epistasis_v2_full.pkl')
    
    # Get path with fallback
    shards_dir = get_path('shards', fallback='Archive/VHH_shards')

Benefits:
    - Works before AND after reorganization
    - Auto-detects which structure is in use
    - Single place to update if paths change
    - Graceful fallbacks
"""

import os
from pathlib import Path
from typing import Optional, Union, List


def _find_root() -> Path:
    """Find the KA-Search root directory."""
    # Try current directory first
    cwd = Path.cwd()
    
    # Check for telltale directories
    markers = ['Archive', 'active', 'tools', 'PKL Files', 'data']
    
    # Walk up to find root
    current = cwd
    for _ in range(5):  # Max 5 levels up
        for marker in markers:
            if (current / marker).exists():
                return current
        current = current.parent
    
    # Fallback to cwd
    return cwd


# Auto-detect root
ROOT = _find_root()


class KASearchPaths:
    """
    Central path configuration with automatic fallback.
    
    Checks new paths first, falls back to old paths if not found.
    """
    
    def __init__(self, root: Path = None):
        self.root = root or ROOT
        self._detect_structure()
    
    def _detect_structure(self):
        """Detect if we're using old or new structure."""
        self.is_reorganized = (self.root / 'active').exists()
        self.is_legacy = (self.root / 'Archive').exists()
    
    def _resolve(self, new_path: str, old_path: str) -> Path:
        """Try new path first, fall back to old."""
        new = self.root / new_path
        old = self.root / old_path
        
        if new.exists():
            return new
        elif old.exists():
            return old
        else:
            # Return new path (will be created)
            return new
    
    # ==========================================================================
    # DATABASE PATHS
    # ==========================================================================
    
    @property
    def database_production(self) -> Path:
        """Current production VHH database."""
        return self._resolve(
            'data/databases/production/VHH_db_unified_v2.npz',
            'Archive/NPZ Files/VHH_db_unified_v2.npz'
        )
    
    @property
    def database_legacy(self) -> Path:
        """Legacy database directory."""
        return self._resolve(
            'data/databases/legacy',
            'Archive/NPZ Files'
        )
    
    @property
    def database_shards(self) -> Path:
        """Sharded database directory."""
        return self._resolve(
            'data/databases/shards',
            'Archive/VHH_shards'
        )
    
    @property
    def databases(self) -> Path:
        """All databases directory."""
        return self._resolve(
            'data/databases',
            'Archive/NPZ Files'
        )
    
    # ==========================================================================
    # MODEL PATHS
    # ==========================================================================
    
    @property
    def epistasis_model(self) -> Path:
        """Current epistasis model."""
        return self._resolve(
            'models/epistasis/current/epistasis_v2_full.pkl',
            'PKL Files/epistasis_v2_full.pkl'
        )
    
    @property
    def epistasis_legacy(self) -> Path:
        """Legacy epistasis model."""
        return self._resolve(
            'models/epistasis/legacy/epistasis_overnight_full.pkl',
            'PKL Files/epistasis_overnight_full.pkl'
        )
    
    @property
    def epistasis_checkpoints(self) -> Path:
        """Epistasis checkpoints directory."""
        return self._resolve(
            'models/epistasis/checkpoints',
            'PKL Files'
        )
    
    @property
    def correlation_results(self) -> Path:
        """Correlation analysis results."""
        return self._resolve(
            'models/correlations/correlation_results.pkl',
            'Model Results/1_correlations/correlation_results.pkl'
        )
    
    @property
    def pfr_cdr_models(self) -> Path:
        """PFR-CDR models."""
        return self._resolve(
            'models/pfr_cdr/pfr_cdr_models.pkl',
            'Model Results/2_models/pfr_cdr_models.pkl'
        )
    
    @property
    def models(self) -> Path:
        """All models directory."""
        return self._resolve(
            'models',
            'PKL Files'
        )
    
    # ==========================================================================
    # DATA PATHS
    # ==========================================================================
    
    @property
    def raw_sequences(self) -> Path:
        """Raw input sequences directory."""
        return self._resolve(
            'data/raw/sequences',
            'Archive/Raw Excel Data'
        )
    
    @property
    def raw_structures(self) -> Path:
        """PDB structure files."""
        return self._resolve(
            'data/raw/structures',
            'Archive/extracted/oas-paper/all_nano_structures'
        )
    
    @property
    def external_data(self) -> Path:
        """External/downloaded data."""
        return self._resolve(
            'data/external',
            'Archive/Databases'
        )
    
    # ==========================================================================
    # RESULTS PATHS
    # ==========================================================================
    
    @property
    def results_analysis(self) -> Path:
        """Analysis results directory."""
        return self._resolve(
            'results/analysis_runs',
            'Excel Results'
        )
    
    @property
    def results_kasearch(self) -> Path:
        """KA-Search run results."""
        return self._resolve(
            'results/kasearch_runs',
            'Archive/runs'
        )
    
    @property
    def results_fullscan(self) -> Path:
        """Fullscan results."""
        return self._resolve(
            'results/fullscan_runs',
            'Archive/runs_fullscan'
        )
    
    # ==========================================================================
    # SCRIPT PATHS
    # ==========================================================================
    
    @property
    def active_analysis(self) -> Path:
        """Active analysis scripts."""
        return self._resolve(
            'active/analysis',
            '.'
        )
    
    @property
    def active_utilities(self) -> Path:
        """Utility scripts."""
        return self._resolve(
            'active/utilities',
            '.'
        )
    
    @property
    def tools(self) -> Path:
        """Maintenance tools."""
        return self._resolve(
            'tools',
            '.'
        )


# Global instance
PATHS = KASearchPaths()


def get_path(name: str, fallback: str = None) -> Path:
    """
    Get a path by name with optional fallback.
    
    Args:
        name: Path name (e.g., 'epistasis_model', 'database_production')
        fallback: Fallback path if name not found
        
    Returns:
        Resolved Path object
    """
    # Try as attribute
    if hasattr(PATHS, name):
        return getattr(PATHS, name)
    
    # Try common aliases
    aliases = {
        'db': 'database_production',
        'database': 'database_production',
        'epistasis': 'epistasis_model',
        'shards': 'database_shards',
        'results': 'results_analysis',
        'sequences': 'raw_sequences',
        'models': 'models',
    }
    
    if name in aliases:
        return getattr(PATHS, aliases[name])
    
    # Use fallback
    if fallback:
        return PATHS.root / fallback
    
    raise ValueError(f"Unknown path name: {name}")


def find_file(filename: str, search_dirs: List[str] = None) -> Optional[Path]:
    """
    Find a file by name, searching common directories.
    
    Args:
        filename: Name of file to find
        search_dirs: Optional list of directories to search
        
    Returns:
        Path to file if found, None otherwise
    """
    if search_dirs is None:
        search_dirs = [
            PATHS.models,
            PATHS.databases,
            PATHS.raw_sequences,
            PATHS.results_analysis,
            PATHS.root,
        ]
    
    for search_dir in search_dirs:
        if not isinstance(search_dir, Path):
            search_dir = PATHS.root / search_dir
        
        # Direct check
        if (search_dir / filename).exists():
            return search_dir / filename
        
        # Recursive search
        matches = list(search_dir.rglob(filename))
        if matches:
            return matches[0]
    
    return None


def ensure_dir(path: Union[str, Path]) -> Path:
    """Ensure a directory exists, create if needed."""
    if isinstance(path, str):
        path = PATHS.root / path
    path.mkdir(parents=True, exist_ok=True)
    return path


def get_output_path(prefix: str, extension: str = 'xlsx') -> Path:
    """
    Generate an output path with date prefix.
    
    Args:
        prefix: Descriptive prefix (e.g., 'naturalness_analysis')
        extension: File extension
        
    Returns:
        Path like results/analysis_runs/2025-12/20251231_naturalness_analysis.xlsx
    """
    from datetime import datetime
    
    date_str = datetime.now().strftime('%Y%m%d')
    month_dir = datetime.now().strftime('%Y-%m')
    
    output_dir = ensure_dir(PATHS.results_analysis / month_dir)
    return output_dir / f"{date_str}_{prefix}.{extension}"


# =============================================================================
# CONVENIENCE FUNCTIONS FOR COMMON FILES
# =============================================================================

def get_epistasis_model() -> Path:
    """Get the current epistasis model, with fallback."""
    model = PATHS.epistasis_model
    if model.exists():
        return model
    
    # Try to find it
    found = find_file('epistasis_v2_full.pkl')
    if found:
        return found
    
    # Last resort
    found = find_file('epistasis_overnight_full.pkl')
    if found:
        return found
    
    raise FileNotFoundError("Could not find epistasis model. Expected at: " + str(model))


def get_production_database() -> Path:
    """Get the production VHH database."""
    db = PATHS.database_production
    if db.exists():
        return db
    
    found = find_file('VHH_db_unified_v2.npz')
    if found:
        return found
    
    raise FileNotFoundError("Could not find production database. Expected at: " + str(db))


def get_shard_files() -> List[Path]:
    """Get all shard NPZ files."""
    shard_dir = PATHS.database_shards
    if shard_dir.exists():
        return sorted(shard_dir.glob('*.npz'))
    
    # Fallback
    old_dir = PATHS.root / 'Archive' / 'VHH_shards'
    if old_dir.exists():
        return sorted(old_dir.glob('*.npz'))
    
    return []


# =============================================================================
# PRINT PATH INFO (for debugging)
# =============================================================================

def print_paths():
    """Print all resolved paths for debugging."""
    print("=" * 60)
    print("KA-SEARCH PATHS")
    print(f"Root: {PATHS.root}")
    print(f"Structure: {'Reorganized' if PATHS.is_reorganized else 'Legacy'}")
    print("=" * 60)
    
    attrs = [a for a in dir(PATHS) if not a.startswith('_') and isinstance(getattr(type(PATHS), a, None), property)]
    
    for attr in sorted(attrs):
        path = getattr(PATHS, attr)
        exists = "✅" if path.exists() else "❌"
        print(f"{exists} {attr}: {path}")


if __name__ == '__main__':
    print_paths()
