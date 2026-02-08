#!/usr/bin/env python3
"""
KA-Search Maintain
===========================
One script to rule them all.

Just run:
    python maintain.py

It will automatically:
    1. Clean Zone.Identifier files
    2. Detect new/unorganized files and offer to move them
    3. Delete duplicate files (source files that already exist in destination)
    4. Detect modified scripts and auto-version them
    5. Generate changelog entries from detected changes
    6. Clean up old folders (Archive/, PKL Files/, etc.)
    7. Show project status

Usage:
    python maintain.py              # Full maintenance pass (interactive)
    python maintain.py --auto       # Non-interactive (auto-approve safe changes)
    python maintain.py --status     # Just show status, don't change anything
    python maintain.py --help       # Show this help
"""

import os
import re
import sys
import ast
import json
import shutil
import hashlib
import difflib
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Tuple, Optional, Set
from dataclasses import dataclass, field


# =============================================================================
# CONFIGURATION
# =============================================================================

# Where different file types should live
FILE_DESTINATIONS = {
    # Scripts by pattern
    'epistasis': 'active/analysis/',
    'naturalness': 'active/analysis/',
    'analyzer': 'active/analysis/',
    'visualize': 'active/analysis/',
    'npz_fullscan': 'active/database/',
    'npz_scan': 'active/database/',
    'build_': 'active/database/',
    'translator': 'active/utilities/',
    'diagnose': 'tools/',
    'debug': 'tools/',
    
    # Data files
    '.csv': 'data/raw/sequences/',
    '.xlsx': 'data/raw/sequences/',
    '.tsv': 'data/raw/sequences/',
    '.fasta': 'data/raw/sequences/',
    '.fa': 'data/raw/sequences/',
    '.npz': 'data/databases/',
    '.pkl': 'models/',
    '.json': 'models/',
}

# Scripts that stay in root/tools
TOOL_SCRIPTS = {'maintain.py', 'kasearch_paths.py', 'kasearch_reorganize.py'}

# Cache file for tracking script states
CACHE_FILE = 'tools/.cache.json'

# Directories to ignore
IGNORE_DIRS = {'.git', '__pycache__', '.venv', 'venv', 'node_modules', 'legacy'}


# =============================================================================
# DATA CLASSES
# =============================================================================

@dataclass
class ScriptState:
    """Tracks state of a script for change detection."""
    path: str
    hash: str
    version: float
    functions: List[str]
    classes: List[str]
    line_count: int
    last_seen: str
    docstring_hash: str = ''


@dataclass
class DetectedChange:
    """A detected change in a script."""
    script_path: Path
    change_type: str  # 'major', 'minor', 'patch'
    description: str
    old_version: float
    new_version: float
    details: List[str] = field(default_factory=list)


# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def find_root() -> Path:
    """Find KA-Search root directory."""
    cwd = Path.cwd()
    markers = ['Archive', 'active', 'tools', 'PKL Files', 'data']
    
    current = cwd
    for _ in range(5):
        for marker in markers:
            if (current / marker).exists():
                return current
        current = current.parent
    return cwd


def compute_hash(filepath: Path) -> str:
    """Compute hash of file content."""
    try:
        content = filepath.read_bytes()
        return hashlib.md5(content).hexdigest()[:12]
    except:
        return ''


def parse_version(filename: str) -> float:
    """Extract version number from filename (supports 3, 3.1, 3.2 etc)."""
    match = re.search(r'_v(\d+)(?:\.(\d+))?(?:\.py)?$', filename)
    if match:
        major = int(match.group(1))
        minor = int(match.group(2)) if match.group(2) else 0
        return major + minor / 10
    return 0.0


def format_version(version: float) -> str:
    """Format version number for filename."""
    if version == int(version):
        return f"v{int(version)}"
    else:
        major = int(version)
        minor = int(round((version - major) * 10))
        return f"v{major}.{minor}"


def analyze_python_file(filepath: Path) -> dict:
    """Analyze a Python file for functions, classes, docstring."""
    info = {
        'functions': [],
        'classes': [],
        'docstring': '',
        'docstring_hash': '',
        'line_count': 0,
        'todos': [],
        'changes': [],  # Lines with # CHANGED, # FIXED, # NEW, etc.
    }
    
    try:
        content = filepath.read_text(encoding='utf-8', errors='replace')
        lines = content.split('\n')
        info['line_count'] = len(lines)
        
        # Find comment markers
        for i, line in enumerate(lines):
            line_upper = line.upper()
            for marker in ['# CHANGED', '# FIXED', '# NEW', '# ADDED', '# REMOVED', '# BUG', '# UPDATE']:
                if marker in line_upper:
                    # Extract the comment
                    comment = line.split('#', 1)[1].strip() if '#' in line else line.strip()
                    info['changes'].append(comment)
        
        # Parse AST
        try:
            tree = ast.parse(content)
            info['docstring'] = ast.get_docstring(tree) or ''
            info['docstring_hash'] = hashlib.md5(info['docstring'].encode()).hexdigest()[:8]
            
            for node in ast.iter_child_nodes(tree):
                if isinstance(node, ast.FunctionDef):
                    info['functions'].append(node.name)
                elif isinstance(node, ast.ClassDef):
                    info['classes'].append(node.name)
        except SyntaxError:
            pass
            
    except Exception as e:
        pass
    
    return info


def determine_change_type(old_state: ScriptState, new_info: dict, new_hash: str) -> Tuple[str, List[str]]:
    """
    Determine if change is major, minor, or patch.
    
    Major (X.0): 
        - Functions added/removed
        - Classes added/removed
        - >20% line count change
        
    Minor (X.Y):
        - Logic changes (hash changed but structure same)
        - Significant edits
        
    Patch (same version, just log):
        - Comments only
        - Docstring only
        - Whitespace
    """
    details = []
    
    old_funcs = set(old_state.functions)
    new_funcs = set(new_info['functions'])
    old_classes = set(old_state.classes)
    new_classes = set(new_info['classes'])
    
    added_funcs = new_funcs - old_funcs
    removed_funcs = old_funcs - new_funcs
    added_classes = new_classes - old_classes
    removed_classes = old_classes - new_classes
    
    # Check for major changes
    is_major = False
    
    if added_funcs:
        details.append(f"Added functions: {', '.join(added_funcs)}")
        is_major = True
    if removed_funcs:
        details.append(f"Removed functions: {', '.join(removed_funcs)}")
        is_major = True
    if added_classes:
        details.append(f"Added classes: {', '.join(added_classes)}")
        is_major = True
    if removed_classes:
        details.append(f"Removed classes: {', '.join(removed_classes)}")
        is_major = True
    
    # Line count change
    if old_state.line_count > 0:
        line_change = abs(new_info['line_count'] - old_state.line_count) / old_state.line_count
        if line_change > 0.20:
            details.append(f"Significant size change: {old_state.line_count} → {new_info['line_count']} lines ({line_change*100:.0f}%)")
            is_major = True
    
    if is_major:
        return 'major', details
    
    # Check for minor changes (code changed but structure same)
    if new_hash != old_state.hash:
        # Check if only docstring changed
        if new_info['docstring_hash'] != old_state.docstring_hash:
            details.append("Docstring updated")
        
        # Use any # CHANGED comments as details
        if new_info['changes']:
            details.extend(new_info['changes'][:3])  # Max 3
        
        if not details:
            details.append("Code modifications")
        
        return 'minor', details
    
    return 'patch', details


def bump_version(current: float, change_type: str) -> float:
    """Calculate new version number."""
    major = int(current) if current else 1
    minor = int(round((current - major) * 10)) if current else 0
    
    if change_type == 'major':
        return float(major + 1)
    elif change_type == 'minor':
        return major + (minor + 1) / 10
    else:
        return current  # patch doesn't bump version


# =============================================================================
# CACHE MANAGEMENT
# =============================================================================

def load_cache(root: Path) -> Dict[str, ScriptState]:
    """Load script state cache."""
    cache_path = root / CACHE_FILE
    if not cache_path.exists():
        return {}
    
    try:
        data = json.loads(cache_path.read_text())
        return {
            path: ScriptState(**state) 
            for path, state in data.items()
        }
    except:
        return {}


def save_cache(root: Path, cache: Dict[str, ScriptState]):
    """Save script state cache."""
    cache_path = root / CACHE_FILE
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    data = {
        path: {
            'path': state.path,
            'hash': state.hash,
            'version': state.version,
            'functions': state.functions,
            'classes': state.classes,
            'line_count': state.line_count,
            'last_seen': state.last_seen,
            'docstring_hash': state.docstring_hash,
        }
        for path, state in cache.items()
    }
    cache_path.write_text(json.dumps(data, indent=2))


# =============================================================================
# MAINTENANCE TASKS
# =============================================================================

def clean_zone_identifiers(root: Path) -> int:
    """Delete all Zone.Identifier files."""
    count = 0
    for f in root.rglob('*'):
        if 'Zone.Identifier' in f.name:
            try:
                f.unlink()
                count += 1
            except:
                pass
    return count


def find_unorganized_files(root: Path) -> List[Tuple[Path, str]]:
    """Find files in root that should be organized elsewhere."""
    unorganized = []
    
    for item in root.iterdir():
        if item.is_dir():
            continue
        if item.name.startswith('.'):
            continue
        if item.name in TOOL_SCRIPTS:
            continue
        if item.suffix == '.md':
            continue
            
        # Determine where it should go
        destination = None
        
        # Check by name pattern
        name_lower = item.name.lower()
        for pattern, dest in FILE_DESTINATIONS.items():
            if pattern.startswith('.'):
                if item.suffix.lower() == pattern:
                    destination = dest
                    break
            elif pattern in name_lower:
                destination = dest
                break
        
        if destination:
            unorganized.append((item, destination))
    
    return unorganized


def detect_script_changes(root: Path, cache: Dict[str, ScriptState]) -> List[DetectedChange]:
    """Detect changes in scripts since last run."""
    changes = []
    
    # Find all Python scripts (not in ignored dirs)
    for py_file in root.rglob('*.py'):
        # Skip ignored
        if any(d in str(py_file) for d in IGNORE_DIRS):
            continue
        if py_file.name in TOOL_SCRIPTS:
            continue
        if py_file.name == 'kasearch_paths.py':
            continue
            
        rel_path = str(py_file.relative_to(root))
        current_hash = compute_hash(py_file)
        current_info = analyze_python_file(py_file)
        current_version = parse_version(py_file.name)
        
        if rel_path in cache:
            old_state = cache[rel_path]
            
            # Check if changed
            if current_hash != old_state.hash:
                change_type, details = determine_change_type(old_state, current_info, current_hash)
                
                if change_type != 'patch' or details:
                    new_version = bump_version(old_state.version or current_version, change_type)
                    
                    # Build description from details or comments
                    if current_info['changes']:
                        description = '; '.join(current_info['changes'][:2])
                    elif details:
                        description = '; '.join(details[:2])
                    else:
                        description = f"{change_type.title()} update"
                    
                    changes.append(DetectedChange(
                        script_path=py_file,
                        change_type=change_type,
                        description=description,
                        old_version=old_state.version or current_version,
                        new_version=new_version,
                        details=details,
                    ))
        else:
            # New script - add to cache
            pass
    
    return changes


def update_script_version(filepath: Path, old_version: float, new_version: float) -> Path:
    """Rename script with new version number."""
    old_name = filepath.name
    
    # Remove old version
    base = re.sub(r'_v\d+(\.\d+)?\.py$', '', old_name)
    base = re.sub(r'\.py$', '', base)
    
    new_name = f"{base}_{format_version(new_version)}.py"
    new_path = filepath.parent / new_name
    
    if new_path != filepath:
        filepath.rename(new_path)
        return new_path
    return filepath


def update_changelog(root: Path, script_name: str, version: float, description: str):
    """Add entry to changelog."""
    changelog_path = root / 'CHANGELOG.md'
    today = datetime.now().strftime('%Y-%m-%d')
    
    # Create if doesn't exist
    if not changelog_path.exists():
        changelog_path.write_text("# KA-Search Changelog\n\n")
    
    content = changelog_path.read_text()
    
    # Format version for display
    ver_str = format_version(version).replace('v', '')
    new_entry = f"\n## {script_name} v{ver_str} ({today})\n- {description}\n"
    
    # Insert after header
    if '\n## ' in content:
        first_entry = content.find('\n## ')
        content = content[:first_entry] + new_entry + content[first_entry:]
    else:
        content += new_entry
    
    changelog_path.write_text(content)


def update_cache_for_script(cache: Dict[str, ScriptState], root: Path, filepath: Path):
    """Update cache with current script state."""
    rel_path = str(filepath.relative_to(root))
    info = analyze_python_file(filepath)
    
    cache[rel_path] = ScriptState(
        path=rel_path,
        hash=compute_hash(filepath),
        version=parse_version(filepath.name),
        functions=info['functions'],
        classes=info['classes'],
        line_count=info['line_count'],
        last_seen=datetime.now().isoformat(),
        docstring_hash=info['docstring_hash'],
    )


# =============================================================================
# MAIN MAINTENANCE ROUTINE
# =============================================================================

def run_maintenance(root: Path, auto_mode: bool = False, status_only: bool = False):
    """Run full maintenance pass."""
    print("=" * 60)
    print("🔧 KA-SEARCH MAINTENANCE")
    print("=" * 60)
    print(f"Root: {root}")
    print(f"Mode: {'Status Only' if status_only else 'Auto' if auto_mode else 'Interactive'}")
    print()
    
    # Load cache
    cache = load_cache(root)
    changes_made = False
    
    # -------------------------------------------------------------------------
    # Step 1: Clean Zone.Identifier files (always)
    # -------------------------------------------------------------------------
    if not status_only:
        zone_count = clean_zone_identifiers(root)
        if zone_count > 0:
            print(f"🧹 Cleaned {zone_count} Zone.Identifier files")
            changes_made = True
    
    # -------------------------------------------------------------------------
    # Step 2: Find unorganized files
    # -------------------------------------------------------------------------
    unorganized = find_unorganized_files(root)
    
    if unorganized:
        print(f"\n📁 UNORGANIZED FILES ({len(unorganized)} found)")
        print("-" * 50)
        
        to_move = []
        to_delete = []  # Source files that already exist in destination
        
        for filepath, destination in unorganized:
            dest_path = root / destination / filepath.name
            if dest_path.exists():
                # File already exists in destination - source is duplicate
                to_delete.append((filepath, dest_path))
                print(f"  {filepath.name}")
                print(f"    → already in {destination} (duplicate)")
            else:
                to_move.append((filepath, destination))
                print(f"  {filepath.name}")
                print(f"    → {destination}")
        
        if not status_only:
            # Handle files to move
            if to_move:
                if auto_mode or input("\nMove these files? [y/n]: ").strip().lower() == 'y':
                    for filepath, destination in to_move:
                        dest_dir = root / destination
                        dest_dir.mkdir(parents=True, exist_ok=True)
                        dest_path = dest_dir / filepath.name
                        
                        shutil.move(str(filepath), str(dest_path))
                        print(f"  ✅ Moved: {filepath.name}")
                        changes_made = True
            
            # Handle duplicates - delete source since dest exists
            if to_delete:
                print(f"\n🗑️  DUPLICATE FILES ({len(to_delete)} found)")
                print("-" * 50)
                for src, dst in to_delete:
                    print(f"  {src.name} (in root)")
                    print(f"    = {dst.relative_to(root)} (already exists)")
                
                if auto_mode or input("\nDelete duplicates from root? [y/n]: ").strip().lower() == 'y':
                    for src, dst in to_delete:
                        try:
                            src.unlink()
                            print(f"  ✅ Deleted: {src.name}")
                            changes_made = True
                        except Exception as e:
                            print(f"  ⚠️ Could not delete {src.name}: {e}")
    else:
        print("\n📁 All files organized ✅")
    
    # -------------------------------------------------------------------------
    # Step 3: Detect and handle script changes
    # -------------------------------------------------------------------------
    changes = detect_script_changes(root, cache)
    
    if changes:
        print(f"\n📝 SCRIPT CHANGES DETECTED ({len(changes)} scripts)")
        print("-" * 50)
        
        for change in changes:
            old_ver = format_version(change.old_version) if change.old_version else 'v0'
            new_ver = format_version(change.new_version)
            
            print(f"\n  {change.script_path.name}")
            print(f"    Type: {change.change_type.upper()}")
            print(f"    Version: {old_ver} → {new_ver}")
            print(f"    Changes: {change.description}")
            
            if change.details:
                for detail in change.details[:3]:
                    print(f"      • {detail}")
        
        if not status_only:
            if auto_mode or input("\nApply version updates and log changes? [y/n]: ").strip().lower() == 'y':
                for change in changes:
                    # Update version in filename
                    if change.change_type in ['major', 'minor']:
                        new_path = update_script_version(
                            change.script_path, 
                            change.old_version, 
                            change.new_version
                        )
                        print(f"  ✅ Renamed: {change.script_path.name} → {new_path.name}")
                    else:
                        new_path = change.script_path
                    
                    # Update changelog
                    script_base = re.sub(r'_v\d+(\.\d+)?\.py$', '', change.script_path.name)
                    script_base = re.sub(r'\.py$', '', script_base)
                    update_changelog(root, script_base, change.new_version, change.description)
                    print(f"  ✅ Logged: {change.description[:50]}...")
                    
                    # Update cache
                    update_cache_for_script(cache, root, new_path)
                    changes_made = True
    else:
        print("\n📝 No script changes detected ✅")
    
    # -------------------------------------------------------------------------
    # Step 4: Update cache for all scripts (so we track them going forward)
    # -------------------------------------------------------------------------
    if not status_only:
        for py_file in root.rglob('*.py'):
            if any(d in str(py_file) for d in IGNORE_DIRS):
                continue
            rel_path = str(py_file.relative_to(root))
            if rel_path not in cache:
                update_cache_for_script(cache, root, py_file)
        
        save_cache(root, cache)
    
    # -------------------------------------------------------------------------
    # Step 5: Check for old folders that should be deleted
    # -------------------------------------------------------------------------
    old_folders_to_delete = []
    old_folder_candidates = [
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
    
    for folder in old_folder_candidates:
        folder_path = root / folder
        if folder_path.exists():
            old_folders_to_delete.append(folder_path)
    
    # Also check if Archive itself is mostly empty
    archive_path = root / 'Archive'
    if archive_path.exists():
        # Check what's left
        remaining = [r for r in archive_path.iterdir() if not r.name.startswith('.')]
        # If only legacy backup or empty dirs remain
        if not remaining or all(r.name.startswith('original_') or r.name == 'legacy' for r in remaining):
            old_folders_to_delete.append(archive_path)
    
    if old_folders_to_delete:
        print(f"\n🗑️  OLD FOLDERS ({len(old_folders_to_delete)} found)")
        print("-" * 50)
        for folder in old_folders_to_delete:
            rel = folder.relative_to(root) if folder.is_relative_to(root) else folder
            # Count files inside
            try:
                file_count = sum(1 for _ in folder.rglob('*') if _.is_file())
                print(f"  {rel}/ ({file_count} files)")
            except:
                print(f"  {rel}/")
        
        if not status_only:
            if auto_mode or input("\nDelete old folders? [y/n]: ").strip().lower() == 'y':
                for folder in old_folders_to_delete:
                    try:
                        # First try to delete any remaining files
                        for item in folder.rglob('*'):
                            if item.is_file():
                                try:
                                    item.unlink()
                                except:
                                    pass
                        # Then delete empty directories from bottom up
                        for item in sorted(folder.rglob('*'), reverse=True):
                            if item.is_dir():
                                try:
                                    item.rmdir()
                                except:
                                    pass
                        # Finally delete the folder itself
                        if folder.exists():
                            shutil.rmtree(folder, ignore_errors=True)
                        
                        if not folder.exists():
                            print(f"  ✅ Deleted: {folder.name}/")
                            changes_made = True
                        else:
                            # Check what's left
                            remaining = list(folder.rglob('*'))
                            print(f"  ⚠️ Could not fully delete {folder.name}/ ({len(remaining)} items remain)")
                    except Exception as e:
                        print(f"  ⚠️ Error with {folder.name}/: {type(e).__name__}: {e}")
    
    # -------------------------------------------------------------------------
    # Step 5: Summary
    # -------------------------------------------------------------------------
    print("\n" + "=" * 60)
    if status_only:
        print("STATUS CHECK COMPLETE")
    elif changes_made:
        print("✅ MAINTENANCE COMPLETE - Changes applied")
    else:
        print("✅ MAINTENANCE COMPLETE - No changes needed")
    print("=" * 60)
    
    # Quick stats
    py_count = len(list(root.rglob('*.py')))
    print(f"\n📊 Project Stats:")
    print(f"   Python scripts: {py_count}")
    print(f"   Tracked in cache: {len(cache)}")
    
    if (root / 'CHANGELOG.md').exists():
        changelog_lines = len((root / 'CHANGELOG.md').read_text().split('\n'))
        print(f"   Changelog entries: ~{changelog_lines // 3}")


# =============================================================================
# MAIN
# =============================================================================

def main():
    if '--help' in sys.argv or '-h' in sys.argv:
        print(__doc__)
        sys.exit(0)
    
    root = find_root()
    
    auto_mode = '--auto' in sys.argv
    status_only = '--status' in sys.argv
    
    run_maintenance(root, auto_mode, status_only)


if __name__ == '__main__':
    main()