#!/usr/bin/env python3
"""
Extract Script Headers/Docstrings
==================================
Scans a directory for Python files and extracts the module-level
docstrings and initial comments to create a documentation log.

Usage:
    python extract_script_headers.py /path/to/scripts
    python extract_script_headers.py /path/to/scripts --output my_log.md
    python extract_script_headers.py /path/to/scripts --format csv
"""

import os
import sys
import ast
import re
import argparse
from pathlib import Path
from datetime import datetime
from typing import Optional, Dict, List, Tuple


def extract_docstring_ast(filepath: str) -> Optional[str]:
    """Extract module-level docstring using AST parsing."""
    try:
        with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
            source = f.read()
        
        tree = ast.parse(source)
        docstring = ast.get_docstring(tree)
        return docstring
    except SyntaxError:
        return None
    except Exception as e:
        return None


def extract_header_comments(filepath: str) -> str:
    """Extract initial comment block (# style) if no docstring."""
    comments = []
    try:
        with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
            for line in f:
                stripped = line.strip()
                # Skip shebang
                if stripped.startswith('#!'):
                    continue
                # Skip encoding declarations
                if stripped.startswith('# -*-') or stripped.startswith('# coding'):
                    continue
                # Collect comments
                if stripped.startswith('#'):
                    comments.append(stripped[1:].strip())
                # Stop at first non-comment, non-empty line
                elif stripped and not stripped.startswith('"""') and not stripped.startswith("'''"):
                    break
                # Also stop if we hit a docstring (we'll get that separately)
                elif stripped.startswith('"""') or stripped.startswith("'''"):
                    break
    except Exception:
        pass
    
    return '\n'.join(comments) if comments else ''


def get_file_stats(filepath: str) -> Dict:
    """Get basic file statistics."""
    try:
        stat = os.stat(filepath)
        with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
            lines = sum(1 for _ in f)
        
        # Count functions and classes
        with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
            source = f.read()
        
        try:
            tree = ast.parse(source)
            functions = sum(1 for node in ast.walk(tree) if isinstance(node, ast.FunctionDef))
            classes = sum(1 for node in ast.walk(tree) if isinstance(node, ast.ClassDef))
        except:
            functions = 0
            classes = 0
        
        return {
            'lines': lines,
            'size_kb': round(stat.st_size / 1024, 1),
            'modified': datetime.fromtimestamp(stat.st_mtime).strftime('%Y-%m-%d %H:%M'),
            'functions': functions,
            'classes': classes
        }
    except Exception as e:
        return {'lines': 0, 'size_kb': 0, 'modified': 'unknown', 'functions': 0, 'classes': 0}


def extract_version(docstring: str) -> Optional[str]:
    """Try to extract version info from docstring."""
    if not docstring:
        return None
    
    patterns = [
        r'v(\d+\.?\d*\.?\d*)',
        r'version\s*:?\s*(\d+\.?\d*\.?\d*)',
        r'V(\d+\.?\d*\.?\d*)',
    ]
    
    for pattern in patterns:
        match = re.search(pattern, docstring, re.IGNORECASE)
        if match:
            return match.group(1)
    
    return None


def scan_directory(directory: str, recursive: bool = True) -> List[Tuple[str, str, Dict]]:
    """Scan directory for Python files and extract headers."""
    results = []
    
    path = Path(directory)
    
    if recursive:
        py_files = list(path.rglob('*.py'))
    else:
        py_files = list(path.glob('*.py'))
    
    # Sort by name
    py_files = sorted(py_files, key=lambda x: x.name.lower())
    
    for py_file in py_files:
        filepath = str(py_file)
        
        # Skip __pycache__ and hidden directories
        if '__pycache__' in filepath or '/.' in filepath:
            continue
        
        # Extract docstring
        docstring = extract_docstring_ast(filepath)
        
        # If no docstring, try header comments
        if not docstring:
            docstring = extract_header_comments(filepath)
        
        if not docstring:
            docstring = "(No documentation found)"
        
        # Get file stats
        stats = get_file_stats(filepath)
        
        results.append((filepath, docstring, stats))
    
    return results


def group_by_folder(results: List[Tuple[str, str, Dict]], base_dir: str) -> Dict[str, List[Tuple[str, str, Dict]]]:
    """Group results by their parent folder."""
    from collections import defaultdict
    
    grouped = defaultdict(list)
    
    for filepath, docstring, stats in results:
        rel_path = os.path.relpath(filepath, base_dir)
        folder = os.path.dirname(rel_path)
        if not folder:
            folder = "(root)"
        grouped[folder].append((filepath, docstring, stats))
    
    # Sort folders: root first, then alphabetically
    sorted_grouped = {}
    if "(root)" in grouped:
        sorted_grouped["(root)"] = grouped["(root)"]
    
    for folder in sorted(grouped.keys()):
        if folder != "(root)":
            sorted_grouped[folder] = grouped[folder]
    
    return sorted_grouped


def format_markdown(results: List[Tuple[str, str, Dict]], base_dir: str) -> str:
    """Format results as Markdown, organized by folder."""
    lines = []
    
    lines.append("# Script Documentation Log")
    lines.append(f"\nGenerated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append(f"Directory: `{base_dir}`")
    lines.append(f"Total scripts: {len(results)}")
    
    # Group by folder
    grouped = group_by_folder(results, base_dir)
    
    # Folder summary
    lines.append(f"\nFolders: {len(grouped)}")
    lines.append("")
    for folder, scripts in grouped.items():
        lines.append(f"- `{folder}/` ({len(scripts)} scripts)")
    
    # Summary table by folder
    lines.append("\n---\n")
    lines.append("## Summary by Folder\n")
    
    for folder, scripts in grouped.items():
        lines.append(f"### 📁 `{folder}/`\n")
        lines.append("| Script | Lines | Size | Functions | Classes | Modified |")
        lines.append("|--------|-------|------|-----------|---------|----------|")
        
        for filepath, docstring, stats in scripts:
            name = Path(filepath).name
            lines.append(f"| {name} | {stats['lines']} | {stats['size_kb']}KB | {stats['functions']} | {stats['classes']} | {stats['modified']} |")
        
        lines.append("")
    
    # Detailed sections by folder
    lines.append("---\n")
    lines.append("## Script Details\n")
    
    for folder, scripts in grouped.items():
        lines.append(f"## 📁 `{folder}/`\n")
        
        for filepath, docstring, stats in scripts:
            name = Path(filepath).name
            rel_path = os.path.relpath(filepath, base_dir)
            version = extract_version(docstring)
            
            lines.append(f"### `{name}`")
            if version:
                lines.append(f"**Version:** {version}")
            lines.append(f"**Path:** `{rel_path}`")
            lines.append(f"**Stats:** {stats['lines']} lines | {stats['size_kb']}KB | {stats['functions']} functions | {stats['classes']} classes")
            lines.append(f"**Modified:** {stats['modified']}")
            lines.append("")
            lines.append("**Description:**")
            lines.append("```")
            lines.append(docstring.strip())
            lines.append("```")
            lines.append("")
            lines.append("---\n")
    
    return '\n'.join(lines)


def format_csv(results: List[Tuple[str, str, Dict]], base_dir: str) -> str:
    """Format results as CSV with folder organization."""
    import csv
    import io
    
    output = io.StringIO()
    writer = csv.writer(output)
    
    # Header - added folder column
    writer.writerow(['folder', 'filename', 'path', 'lines', 'size_kb', 'functions', 'classes', 'modified', 'version', 'docstring'])
    
    # Group and sort by folder
    grouped = group_by_folder(results, base_dir)
    
    for folder, scripts in grouped.items():
        for filepath, docstring, stats in scripts:
            name = Path(filepath).name
            rel_path = os.path.relpath(filepath, base_dir)
            version = extract_version(docstring) or ''
            
            # Clean docstring for CSV
            clean_doc = docstring.replace('\n', ' | ').strip()
            
            writer.writerow([
                folder,
                name,
                rel_path,
                stats['lines'],
                stats['size_kb'],
                stats['functions'],
                stats['classes'],
                stats['modified'],
                version,
                clean_doc
            ])
    
    return output.getvalue()


def format_json(results: List[Tuple[str, str, Dict]], base_dir: str) -> str:
    """Format results as JSON, organized by folder."""
    import json
    
    grouped = group_by_folder(results, base_dir)
    
    output = {
        'generated': datetime.now().isoformat(),
        'directory': base_dir,
        'total_scripts': len(results),
        'total_folders': len(grouped),
        'folders': {}
    }
    
    for folder, scripts in grouped.items():
        output['folders'][folder] = {
            'script_count': len(scripts),
            'scripts': []
        }
        
        for filepath, docstring, stats in scripts:
            name = Path(filepath).name
            rel_path = os.path.relpath(filepath, base_dir)
            version = extract_version(docstring)
            
            output['folders'][folder]['scripts'].append({
                'filename': name,
                'path': rel_path,
                'version': version,
                'docstring': docstring.strip(),
                'stats': stats
            })
    
    return json.dumps(output, indent=2)


def main():
    parser = argparse.ArgumentParser(
        description='Extract docstrings and headers from Python scripts',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python extract_script_headers.py ~/KA-Search/
  python extract_script_headers.py ~/KA-Search/ --output scripts_log.md
  python extract_script_headers.py ~/KA-Search/ --format csv --output scripts.csv
  python extract_script_headers.py ~/KA-Search/ --format json --output scripts.json
  python extract_script_headers.py ~/KA-Search/ --no-recursive
        """
    )
    
    parser.add_argument('directory', nargs='?', default='.', 
                        help='Directory to scan (default: current directory)')
    parser.add_argument('--output', '-o', 
                        help='Output file (default: print to stdout)')
    parser.add_argument('--format', '-f', choices=['markdown', 'csv', 'json'], 
                        default='markdown', help='Output format (default: markdown)')
    parser.add_argument('--no-recursive', action='store_true',
                        help='Do not scan subdirectories')
    
    args = parser.parse_args()
    
    # Expand path
    directory = os.path.expanduser(args.directory)
    
    if not os.path.isdir(directory):
        print(f"Error: '{directory}' is not a valid directory", file=sys.stderr)
        sys.exit(1)
    
    print(f"Scanning: {directory}", file=sys.stderr)
    
    # Scan
    results = scan_directory(directory, recursive=not args.no_recursive)
    
    if not results:
        print("No Python files found.", file=sys.stderr)
        sys.exit(0)
    
    print(f"Found {len(results)} Python files", file=sys.stderr)
    
    # Format output
    if args.format == 'markdown':
        output = format_markdown(results, directory)
    elif args.format == 'csv':
        output = format_csv(results, directory)
    elif args.format == 'json':
        output = format_json(results, directory)
    
    # Write or print
    if args.output:
        output_path = os.path.expanduser(args.output)
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write(output)
        print(f"Saved to: {output_path}", file=sys.stderr)
    else:
        print(output)


if __name__ == '__main__':
    main()
