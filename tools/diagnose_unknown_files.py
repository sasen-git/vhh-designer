#!/usr/bin/env python3
"""
Diagnose Unknown Files in KA-Search Directory
==============================================
Extracts docstrings, function names, and key patterns to help categorize files.

Run this on your local machine:
    python diagnose_unknown_files.py /path/to/KA-Search

Output: unknown_files_report.txt with summaries of each file
"""

import os
import re
import ast
import sys
from pathlib import Path
from datetime import datetime

# Files to specifically examine (adjust paths as needed)
UNKNOWN_FILES = [
    "dna_translator.py",
    "pathlib_list.py", 
    "M69.fasta",
    "vhh_framework_optimizer.py",
    "score_framework.py",
    "visualize_epistasis.py",
    "visualize_correlations.py",
    "visualize_pfr_models.py",
    "Archive/One-off Py Scripts/process_sabdab_pdb.py",
    "Archive/One-off Py Scripts/process_vhh_with_antigens.py",
    "Archive/One-off Py Scripts/vhh_global_compensation_analysis.py",
    "Archive/One-off Py Scripts/analyze_cdr_framework_advanced.py",
    "Archive/One-off Py Scripts/analyze_cdr_framework_correlations.py",
    "Archive/One-off Py Scripts/vhh_full_epistasis_overnight.py",
    "Archive/One-off Py Scripts/build_camel_vhh_db.py",
    "Archive/One-off Py Scripts/build_camel_vhh_db_one_step.py",
    "Archive/One-off Py Scripts/build_models_only.py",
    "Archive/One-off Py Scripts/build_pfr_cdr_models.py",
    "Archive/One-off Py Scripts/npz_fullscan_v2.py",
    "Archive/One-off Py Scripts/npz_fullscan_v3.py",
    "Archive/One-off Py Scripts/npz_fullscan_v5.py",
    "Archive/One-off Py Scripts/run_full_analysis.py",
    "Archive/One-off Py Scripts/run_full_analysis_fast.py",
    "Archive/One-off Py Scripts/run_full_analysis_fast2.py",
    "Archive/One-off Py Scripts/run_full_analysis_lowmem.py",
    "Archive/One-off Py Scripts/run_full_analysis_lowmem_v3.py",
]


def extract_python_info(filepath: Path) -> dict:
    """Extract docstring, functions, classes, and imports from a Python file."""
    info = {
        'type': 'python',
        'docstring': None,
        'functions': [],
        'classes': [],
        'imports': [],
        'main_block': False,
        'argparse': False,
        'key_patterns': [],
        'first_50_lines': [],
    }
    
    try:
        content = filepath.read_text(encoding='utf-8', errors='replace')
        lines = content.split('\n')
        info['first_50_lines'] = lines[:50]
        info['total_lines'] = len(lines)
        
        # Try to parse AST
        try:
            tree = ast.parse(content)
            
            # Get module docstring
            info['docstring'] = ast.get_docstring(tree)
            
            # Get top-level functions and classes
            for node in ast.iter_child_nodes(tree):
                if isinstance(node, ast.FunctionDef):
                    info['functions'].append(node.name)
                elif isinstance(node, ast.ClassDef):
                    info['classes'].append(node.name)
                elif isinstance(node, ast.Import):
                    for alias in node.names:
                        info['imports'].append(alias.name)
                elif isinstance(node, ast.ImportFrom):
                    info['imports'].append(f"from {node.module}")
                        
        except SyntaxError:
            # Fall back to regex
            info['docstring'] = "PARSE_ERROR"
            
        # Check for patterns
        if 'if __name__' in content:
            info['main_block'] = True
        if 'argparse' in content:
            info['argparse'] = True
        if 'pickle' in content:
            info['key_patterns'].append('uses_pickle')
        if 'pandas' in content:
            info['key_patterns'].append('uses_pandas')
        if 'numpy' in content or 'np.' in content:
            info['key_patterns'].append('uses_numpy')
        if 'matplotlib' in content or 'plt.' in content:
            info['key_patterns'].append('visualization')
        if 'kasearch' in content.lower():
            info['key_patterns'].append('uses_kasearch')
        if 'anarci' in content.lower():
            info['key_patterns'].append('uses_anarci')
        if '.npz' in content:
            info['key_patterns'].append('npz_files')
        if 'CDR3' in content or 'cdr3' in content:
            info['key_patterns'].append('cdr3_analysis')
        if 'epistasis' in content.lower():
            info['key_patterns'].append('epistasis')
        if 'framework' in content.lower():
            info['key_patterns'].append('framework_analysis')
            
    except Exception as e:
        info['error'] = str(e)
        
    return info


def extract_fasta_info(filepath: Path) -> dict:
    """Extract info from a FASTA file."""
    info = {
        'type': 'fasta',
        'sequences': [],
        'total_sequences': 0,
    }
    
    try:
        content = filepath.read_text(encoding='utf-8', errors='replace')
        lines = content.strip().split('\n')
        
        current_header = None
        current_seq = []
        
        for line in lines:
            if line.startswith('>'):
                if current_header:
                    info['sequences'].append({
                        'header': current_header,
                        'length': len(''.join(current_seq)),
                        'seq_preview': ''.join(current_seq)[:50]
                    })
                current_header = line[1:].strip()
                current_seq = []
            else:
                current_seq.append(line.strip())
                
        # Last sequence
        if current_header:
            info['sequences'].append({
                'header': current_header,
                'length': len(''.join(current_seq)),
                'seq_preview': ''.join(current_seq)[:50]
            })
            
        info['total_sequences'] = len(info['sequences'])
        
    except Exception as e:
        info['error'] = str(e)
        
    return info


def analyze_file(filepath: Path) -> dict:
    """Analyze a file and return its info."""
    if not filepath.exists():
        return {'error': 'FILE_NOT_FOUND', 'path': str(filepath)}
    
    suffix = filepath.suffix.lower()
    
    if suffix == '.py':
        return extract_python_info(filepath)
    elif suffix in ['.fasta', '.fa', '.faa']:
        return extract_fasta_info(filepath)
    else:
        # Generic file info
        try:
            stat = filepath.stat()
            return {
                'type': 'other',
                'suffix': suffix,
                'size_kb': stat.st_size / 1024,
                'modified': datetime.fromtimestamp(stat.st_mtime).isoformat(),
            }
        except Exception as e:
            return {'error': str(e)}


def generate_report(root_dir: Path, output_file: Path):
    """Generate a report of all unknown files."""
    report_lines = [
        "=" * 80,
        "UNKNOWN FILES DIAGNOSTIC REPORT",
        f"Root: {root_dir}",
        f"Generated: {datetime.now().isoformat()}",
        "=" * 80,
        "",
    ]
    
    for rel_path in UNKNOWN_FILES:
        filepath = root_dir / rel_path
        report_lines.append("-" * 80)
        report_lines.append(f"FILE: {rel_path}")
        report_lines.append("-" * 80)
        
        info = analyze_file(filepath)
        
        if 'error' in info:
            report_lines.append(f"  ERROR: {info['error']}")
            report_lines.append("")
            continue
            
        report_lines.append(f"  Type: {info.get('type', 'unknown')}")
        
        if info['type'] == 'python':
            report_lines.append(f"  Lines: {info.get('total_lines', '?')}")
            report_lines.append(f"  Has main block: {info.get('main_block', False)}")
            report_lines.append(f"  Has argparse: {info.get('argparse', False)}")
            
            if info.get('docstring'):
                doc_preview = info['docstring'][:300].replace('\n', '\n    ')
                report_lines.append(f"  Docstring:\n    {doc_preview}...")
            else:
                report_lines.append("  Docstring: None")
                
            if info.get('functions'):
                report_lines.append(f"  Functions: {', '.join(info['functions'][:10])}")
            if info.get('classes'):
                report_lines.append(f"  Classes: {', '.join(info['classes'][:10])}")
            if info.get('key_patterns'):
                report_lines.append(f"  Patterns: {', '.join(info['key_patterns'])}")
                
            # Show first few lines if no docstring
            if not info.get('docstring'):
                report_lines.append("  First 20 lines:")
                for line in info.get('first_50_lines', [])[:20]:
                    report_lines.append(f"    {line}")
                    
        elif info['type'] == 'fasta':
            report_lines.append(f"  Sequences: {info.get('total_sequences', 0)}")
            for seq in info.get('sequences', [])[:5]:
                report_lines.append(f"    > {seq['header']}")
                report_lines.append(f"      Length: {seq['length']} aa")
                report_lines.append(f"      Preview: {seq['seq_preview']}...")
                
        report_lines.append("")
        
    # Write report
    output_file.write_text('\n'.join(report_lines))
    print(f"Report written to: {output_file}")
    

def main():
    if len(sys.argv) < 2:
        print("Usage: python diagnose_unknown_files.py /path/to/KA-Search")
        print("\nThis will examine unknown files and generate a report.")
        sys.exit(1)
        
    root_dir = Path(sys.argv[1])
    if not root_dir.exists():
        print(f"Error: Directory not found: {root_dir}")
        sys.exit(1)
        
    output_file = root_dir / "unknown_files_report.txt"
    generate_report(root_dir, output_file)
    

if __name__ == '__main__':
    main()
