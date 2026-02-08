#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
NPZ FULLSCAN v6 INTERACTIVE - Enhanced Interactive Mode
Full-featured antibody search with interactive setup

FEATURES:
- Run without arguments for fully interactive mode
- Interactive threshold analysis (sampling mode)
- Full exhaustive scan with all v6 improvements
- Automatic recommendation of optimal thresholds

USAGE:
1. Interactive mode (just hit enter):
   python3 npz_fullscan_v6_interactive.py
   
2. Command line mode (as before):
   python3 npz_fullscan_v6_interactive.py --db-root /path/to/db --query-seq "EIQLQQ..."
"""

import os
import re
import time
import glob
import signal
import argparse
import datetime as dt
import sys
import json
import platform
from typing import Dict, List, Tuple, Optional
from pathlib import Path
import random

import numpy as np
import pandas as pd
from tqdm import tqdm

# Optional: ANARCI for query extraction
try:
    from anarci import run_anarci
    _ANARCI_AVAILABLE = True
except Exception:
    run_anarci = None
    _ANARCI_AVAILABLE = False

# ------------------------
# Excel Export
# ------------------------
def export_to_excel(csv_path: str, results_csv_path: str, output_folder: str) -> Optional[str]:
    """
    Export CSV results to a formatted Excel file in the main output folder.
    CSVs are in subfolder, Excel goes in main folder.
    Returns the path to the Excel file if successful, None otherwise.
    """
    try:
        # Create Excel path in MAIN folder (not CSV subfolder)
        base_name = os.path.basename(csv_path).replace('.csv', '')
        excel_path = os.path.join(output_folder, base_name + '.xlsx')
        
        # Read the summary CSV
        with open(csv_path, 'r') as f:
            lines = f.readlines()
            # Find where data starts (after metadata)
            data_start = 0
            for i, line in enumerate(lines):
                if line.startswith('shard_id,'):
                    data_start = i
                    break
        
        # Read summary data
        summary_df = pd.read_csv(csv_path, skiprows=data_start)
        
        # Read results data if it exists
        results_df = None
        if os.path.exists(results_csv_path):
            # Skip comment lines
            with open(results_csv_path, 'r') as f:
                lines = f.readlines()
                data_start = 0
                for i, line in enumerate(lines):
                    if not line.startswith('#'):
                        data_start = i
                        break
            results_df = pd.read_csv(results_csv_path, skiprows=data_start)
        
        # Create Excel writer
        with pd.ExcelWriter(excel_path, engine='openpyxl') as writer:
            # Write summary sheet
            summary_df.to_excel(writer, sheet_name='Summary', index=False)
            
            # Write results sheet if available
            if results_df is not None and not results_df.empty:
                # Split large results into multiple sheets if needed
                max_rows = 1048576 - 2  # Excel row limit minus header
                if len(results_df) > max_rows:
                    for i, start_idx in enumerate(range(0, len(results_df), max_rows)):
                        end_idx = min(start_idx + max_rows, len(results_df))
                        sheet_name = f'Results_{i+1}'
                        results_df.iloc[start_idx:end_idx].to_excel(
                            writer, sheet_name=sheet_name, index=False
                        )
                else:
                    results_df.to_excel(writer, sheet_name='Results', index=False)
            
            # Add metadata sheet
            metadata = {
                'File': [os.path.basename(csv_path)],
                'Created': [dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S')],
                'Total Hits': [len(results_df) if results_df is not None else 0],
                'Summary Rows': [len(summary_df)]
            }
            metadata_df = pd.DataFrame(metadata)
            metadata_df.to_excel(writer, sheet_name='Info', index=False)
            
            # Format the sheets
            for sheet_name in writer.sheets:
                worksheet = writer.sheets[sheet_name]
                # Adjust column widths
                for column_cells in worksheet.columns:
                    length = max(len(str(cell.value or '')) for cell in column_cells)
                    adjusted_width = min(length + 2, 50)
                    worksheet.column_dimensions[column_cells[0].column_letter].width = adjusted_width
        
        return excel_path
    
    except Exception as e:
        print(f"Warning: Could not create Excel file: {e}")
        print("Results are saved in CSV format.")
        return None

# ------------------------
# Constants
# ------------------------
AA_PAD = 124  
AA_GAP = 0

# Amino acid similarity groups for position constraints
SIMILAR_GROUPS = {
    'H': ['H', 'K', 'R'],      # Positive charged
    'K': ['H', 'K', 'R'],      # Positive charged  
    'R': ['H', 'K', 'R'],      # Positive charged
    'D': ['D', 'E'],           # Negative charged
    'E': ['D', 'E'],           # Negative charged
    'S': ['S', 'T'],           # Polar uncharged (small)
    'T': ['S', 'T'],           # Polar uncharged (small)
    'N': ['N', 'Q'],           # Polar uncharged (amide)
    'Q': ['N', 'Q'],           # Polar uncharged (amide)
    'L': ['L', 'I', 'V'],      # Hydrophobic aliphatic
    'I': ['L', 'I', 'V'],      # Hydrophobic aliphatic
    'V': ['L', 'I', 'V'],      # Hydrophobic aliphatic
    'F': ['F', 'Y', 'W'],      # Aromatic
    'Y': ['F', 'Y', 'W'],      # Aromatic
    'W': ['F', 'Y', 'W'],      # Aromatic
    'A': ['A', 'G'],           # Small nonpolar
    'G': ['A', 'G'],           # Small nonpolar
    'C': ['C'],                # Cysteine (special)
    'M': ['M'],                # Methionine (special)
    'P': ['P'],                # Proline (special)
}

# ------------------------
# Timeout Handler
# ------------------------
class _Timeout(Exception):
    pass

def _alarm_handler(signum, frame):
    raise _Timeout()

# ------------------------
# Decoder functions
# ------------------------
def _ints_to_aa(int_row):
    """Ultra-robust decoder for NPZ integer arrays."""
    out = []
    for v in int_row:
        if v is None:
            continue
            
        if isinstance(v, (str, bytes)):
            if isinstance(v, bytes):
                try:
                    v = v.decode('utf-8', errors='ignore')
                except:
                    continue
            v = str(v).strip()
            if len(v) == 1 and v.upper() in 'ACDEFGHIKLMNPQRSTVWY':
                out.append(v.upper())
            continue
        
        try:
            v_str = str(v).strip()
            if len(v_str) == 1 and v_str.upper() in 'ACDEFGHIKLMNPQRSTVWY':
                out.append(v_str.upper())
                continue
        except:
            pass
        
        try:
            iv = int(v)
            if iv == AA_GAP or iv == AA_PAD:
                continue
            elif 65 <= iv <= 90:
                out.append(chr(iv))
            elif 97 <= iv <= 122:
                out.append(chr(iv).upper())
        except (ValueError, TypeError):
            continue
    return "".join(out)

# ------------------------
# CDR/FR Extraction
# ------------------------
def extract_cdr1_from_npz(int_row):
    """Extract CDR-H1 using correct IMGT positions (27-38)."""
    # IMGT CDR-H1 is at positions 27-38
    # Adjusted for NPZ array structure
    cdr1 = _ints_to_aa(int_row[26:39])  # Was [39:46], now corrected
    return cdr1.replace("-", "")

def extract_cdr2_from_npz(int_row):
    """Extract CDR-H2 using correct IMGT positions (56-65)."""
    # IMGT CDR-H2 is at positions 56-65
    # The old code was incorrectly combining two regions
    cdr2 = _ints_to_aa(int_row[55:66])  # Was [66:75]+[85:91], now corrected
    return cdr2.replace("-", "")

def extract_cdr3_from_npz(int_row):
    """Extract CDR-H3 using correct IMGT positions (105-117, after C at 104)."""
    # IMGT CDR-H3 starts at position 105 (after conserved C at 104)
    # Extract broader region and find actual CDR3 boundaries
    cdr3_region = int_row[104:128]  # Was [150:190], now corrected
    cdr3_aa = _ints_to_aa(cdr3_region)
    
    if not cdr3_aa:
        return ""
    
    # CDR3 starts AFTER the conserved C (not including it)
    # If first char is C, skip it (that's position 104)
    if cdr3_aa and cdr3_aa[0] == 'C':
        cdr3_aa = cdr3_aa[1:]
    
    # Find the W that marks the end of CDR3
    w_pos = cdr3_aa.find('W')
    if w_pos > 0:
        return cdr3_aa[:w_pos].replace("-", "")
    
    # If no W found, take reasonable length (up to 13 residues)
    return cdr3_aa[:13].replace("-", "")

def extract_frameworks_from_npz(int_row):
    """Extract framework regions (everything except CDRs)."""
    # Framework regions based on IMGT numbering:
    # FR1: positions 1-26
    # FR2: positions 39-55
    # FR3: positions 66-104
    # FR4: positions 118-128+
    
    fr1 = _ints_to_aa(int_row[0:26])    # Before CDR1
    fr2 = _ints_to_aa(int_row[39:55])   # Between CDR1 and CDR2
    fr3 = _ints_to_aa(int_row[66:104])  # Between CDR2 and CDR3
    fr4 = _ints_to_aa(int_row[118:128]) # After CDR3
    
    # Combine and clean
    frameworks = (fr1 + fr2 + fr3 + fr4).replace("-", "")
    return frameworks

def extract_whole_from_npz(int_row):
    """Extract full sequence from NPZ row."""
    # Extract the complete variable region
    # Based on analysis, the sequence starts earlier than position 16
    full_seq = _ints_to_aa(int_row[0:128])  # Was [16:190], now broader range
    # Remove gaps and clean up
    return full_seq.replace("-", "")

def heuristic_cdrh3(full_seq: str) -> str:
    """Backup CDR3 extraction using heuristics."""
    s = re.sub(r"[^A-Z]", "", full_seq.upper())
    end = None
    for pat in [r'WGQG', r'W.QG', r'WG.G', r'W..G', r'WG..']:
        m = list(re.finditer(pat, s))
        if m:
            end = m[-1].start()
            break
    if end is None:
        w = s.rfind('W', max(0, len(s) - 30))
        if w == -1:
            return ""
        end = w
    c_start = s.rfind('C', max(0, end - 35), end)
    if c_start == -1:
        c_start = s.rfind('C', 0, end)
        if c_start == -1:
            return ""
    if end - c_start < 3:
        return ""
    return s[c_start:end]

# ------------------------
# Identity calculation
# ------------------------
def edit_dist(s1: str, s2: str) -> int:
    """Calculate Levenshtein distance between two strings."""
    if len(s1) < len(s2):
        return edit_dist(s2, s1)
    
    if len(s2) == 0:
        return len(s1)
    
    prev_row = range(len(s2) + 1)
    
    for i, c1 in enumerate(s1):
        curr_row = [i + 1]
        for j, c2 in enumerate(s2):
            insertions = prev_row[j + 1] + 1
            deletions = curr_row[j] + 1
            substitutions = prev_row[j] + (c1 != c2)
            curr_row.append(min(insertions, deletions, substitutions))
        prev_row = curr_row
    
    return prev_row[-1]

def calc_id_by_region(qval: str, dval: str) -> float:
    """Calculate identity for a region using edit distance."""
    if not qval or not dval:
        return 0.0
    max_len = max(len(qval), len(dval))
    if max_len == 0:
        return 0.0
    dist = edit_dist(qval, dval)
    return 1.0 - (dist / max_len)

# ------------------------
# Query CDR extraction via ANARCI
# ------------------------
def extract_cdrs_from_query_sequence(seq: str, scheme: str = "imgt", timeout_s: int = 5) -> Dict[str, str]:
    """Extract CDRs from query using ANARCI."""
    if not _ANARCI_AVAILABLE:
        h3 = heuristic_cdrh3(seq)
        return {"cdr1": "", "cdr2": "", "cdr3": h3}
    
    seq_clean = seq.upper().replace("-", "").replace(".", "")
    
    if hasattr(signal, 'SIGALRM'):
        signal.signal(signal.SIGALRM, _alarm_handler)
        signal.alarm(timeout_s)
    
    try:
        res = run_anarci([("H", seq_clean)], scheme=scheme, allowed_species=None)
        if hasattr(signal, 'SIGALRM'):
            signal.alarm(0)
        
        if not res or not res[1] or not res[1][0] or not res[1][0][0]:
            return {"cdr1": "", "cdr2": "", "cdr3": heuristic_cdrh3(seq_clean)}
        
        numbering = res[1][0][0][0]
        
        def pos_val(pos_tuple):
            pos = pos_tuple[0]
            ins = pos_tuple[1] if len(pos_tuple) > 1 else ""
            if isinstance(ins, str) and ins.strip():
                return float(f"{pos}.{ord(ins.upper())-64:02d}")
            return float(pos)
        
        def grab(start, end):
            out = []
            for entry in numbering:
                if not isinstance(entry, (tuple, list)) or len(entry) != 2:
                    continue
                t, aa = entry
                if not isinstance(t, (tuple, list)) or aa in (None, "-"):
                    continue
                v = pos_val(t)
                if start <= v <= end + 0.09:
                    out.append(aa)
            return "".join(out)
        
        return {
            "cdr1": grab(27, 38),
            "cdr2": grab(56, 65),
            "cdr3": grab(105, 117)
        }
    
    except _Timeout:
        if hasattr(signal, 'SIGALRM'):
            signal.alarm(0)
        return {"cdr1": "", "cdr2": "", "cdr3": heuristic_cdrh3(seq_clean)}
    except Exception:
        if hasattr(signal, 'SIGALRM'):
            signal.alarm(0)
        return {"cdr1": "", "cdr2": "", "cdr3": heuristic_cdrh3(seq_clean)}

# ------------------------
# Interactive Mode Functions
# ------------------------
def get_user_input_interactive():
    """Get all parameters interactively from the user."""
    print("\n" + "="*80)
    print("NPZ FULLSCAN v6 - INTERACTIVE MODE")
    print("="*80)
    print("Welcome! Let's set up your antibody search.\n")
    
    params = {}
    
    # Pre-configured database paths
    DB_BASE_PATHS = [
        os.path.expanduser("~/KA-Search/extracted/oas-paper/oasdb_20230111/Heavy"),
        "/home/sasenefrem/KA-Search/extracted/oas-paper/oasdb_20230111/Heavy",
        "./Heavy",  # Fallback for local directory
    ]
    
    # Available species
    SPECIES = {
        '1': 'Human',
        '2': 'Mouse',
        '3': 'Camel',
        '4': 'Rabbit',
        '5': 'Rat',
        '6': 'Rhesus',
        '7': 'Humanised',
    }
    
    # Find the base path that exists
    db_base = None
    for path in DB_BASE_PATHS:
        if os.path.exists(path):
            db_base = path
            break
    
    # Get database path
    while True:
        print("Step 1: Species Selection")
        print("-" * 40)
        
        if db_base:
            print(f"Database location: {db_base}\n")
            print("Select species:")
            for key, species in SPECIES.items():
                species_path = os.path.join(db_base, species)
                if os.path.exists(species_path):
                    npz_count = len(glob.glob(os.path.join(species_path, "*.npz")))
                    if npz_count > 0:
                        print(f"  {key}. {species:<12} ({npz_count} shards)")
                    else:
                        print(f"  {key}. {species:<12} (no data)")
                else:
                    print(f"  {key}. {species:<12} (not found)")
            print("  8. Custom path")
            
            choice = input("\nYour choice (1-8): ").strip()
            
            if choice in SPECIES:
                species = SPECIES[choice]
                db_root = os.path.join(db_base, species)
                
                if not os.path.exists(db_root):
                    print(f"  ⚠ Path does not exist: {db_root}")
                    continue
                
                npz_files = glob.glob(os.path.join(db_root, "*.npz"))
                if not npz_files:
                    print(f"  ⚠ No NPZ files found for {species}")
                    cont = input("  Continue anyway? (y/n): ").lower()
                    if cont != 'y':
                        continue
                else:
                    print(f"  ✓ Selected {species}: {len(npz_files)} NPZ files")
                
                params['db_root'] = db_root
                params['species'] = species
                break
                
            elif choice == '8':
                # Custom path option
                db_root = input("Enter custom database path: ").strip()
                if not db_root:
                    print("  ⚠ Database path cannot be empty. Please try again.\n")
                    continue
                
                # Expand user path
                db_root = os.path.expanduser(db_root)
                
                if not os.path.exists(db_root):
                    print(f"  ⚠ Path does not exist: {db_root}")
                    retry = input("  Try again? (y/n): ").lower()
                    if retry == 'y':
                        continue
                    
                npz_files = glob.glob(os.path.join(db_root, "*.npz"))
                if not npz_files:
                    print(f"  ⚠ No NPZ files found in: {db_root}")
                    cont = input("  Continue anyway? (y/n): ").lower()
                    if cont != 'y':
                        continue
                else:
                    print(f"  ✓ Found {len(npz_files)} NPZ files")
                
                params['db_root'] = db_root
                break
            else:
                print("  ⚠ Invalid choice. Please select 1-8.")
        else:
            # No default database found, ask for custom path
            print("Default database not found at expected locations.")
            print("Please enter the path to your database.\n")
            
            db_root = input("Enter database path (e.g., /path/to/Heavy/Human): ").strip()
            if not db_root:
                print("  ⚠ Database path cannot be empty. Please try again.\n")
                continue
            
            # Expand user path
            db_root = os.path.expanduser(db_root)
            
            # Check if path exists and has NPZ files
            if not os.path.exists(db_root):
                print(f"  ⚠ Path does not exist: {db_root}")
                create = input("  Would you like to try another path? (y/n): ").lower()
                if create != 'n':
                    continue
            
            npz_files = glob.glob(os.path.join(db_root, "*.npz"))
            if not npz_files:
                print(f"  ⚠ No NPZ files found in: {db_root}")
                cont = input("  Continue anyway? (y/n): ").lower()
                if cont != 'y':
                    continue
            else:
                print(f"  ✓ Found {len(npz_files)} NPZ files")
            
            params['db_root'] = db_root
            break
    
    # Get query sequence
    print("\nStep 2: Query Sequence")
    print("-" * 40)
    
    # Check for saved queries
    saved_queries_dir = os.path.expanduser("~/.antibody_queries")
    saved_files = []
    if os.path.exists(saved_queries_dir):
        saved_files = [f for f in os.listdir(saved_queries_dir) if f.endswith(('.txt', '.fasta', '.fa'))]
    
    if saved_files:
        print("Found saved queries:")
        for i, fname in enumerate(saved_files, 1):
            print(f"  {i}. {fname}")
        print(f"  {len(saved_files)+1}. Enter new sequence")
        print(f"  {len(saved_files)+2}. Load from custom file")
        
        choice = input(f"\nYour choice (1-{len(saved_files)+2}): ").strip()
        
        try:
            choice_num = int(choice)
            if 1 <= choice_num <= len(saved_files):
                # Load saved query
                query_file = os.path.join(saved_queries_dir, saved_files[choice_num-1])
                with open(query_file, 'r') as f:
                    content = f.read()
                    if content.startswith('>'):
                        # FASTA format
                        lines = content.split('\n')
                        query_seq = ''.join([line for line in lines[1:] if not line.startswith('>')])
                    else:
                        query_seq = content
                print(f"  ✓ Loaded {saved_files[choice_num-1]}")
            elif choice_num == len(saved_files) + 1:
                # Enter new sequence
                print("Enter your query sequence (paste and press Enter):")
                query_seq = input().strip()
            elif choice_num == len(saved_files) + 2:
                # Load from custom file
                file_path = input("Enter file path: ").strip()
                file_path = os.path.expanduser(file_path)
                with open(file_path, 'r') as f:
                    content = f.read()
                    if content.startswith('>'):
                        lines = content.split('\n')
                        query_seq = ''.join([line for line in lines[1:] if not line.startswith('>')])
                    else:
                        query_seq = content
                print(f"  ✓ Loaded from {file_path}")
            else:
                raise ValueError("Invalid choice")
        except (ValueError, FileNotFoundError) as e:
            print(f"  Error: {e}")
            print("Enter your query sequence (paste and press Enter):")
            query_seq = input().strip()
    else:
        print("Options:")
        print("  1. Paste sequence directly")
        print("  2. Load from file")
        
        choice = input("\nYour choice (1-2): ").strip()
        
        if choice == '2':
            file_path = input("Enter file path: ").strip()
            file_path = os.path.expanduser(file_path)
            try:
                with open(file_path, 'r') as f:
                    content = f.read()
                    if content.startswith('>'):
                        lines = content.split('\n')
                        query_seq = ''.join([line for line in lines[1:] if not line.startswith('>')])
                    else:
                        query_seq = content
                print(f"  ✓ Loaded from {file_path}")
            except Exception as e:
                print(f"  Error: {e}")
                print("Enter your query sequence (paste and press Enter):")
                query_seq = input().strip()
        else:
            print("Enter your query sequence (paste and press Enter):")
            query_seq = input().strip()
    
    while True:
        if not query_seq:
            print("  ⚠ Query sequence cannot be empty. Please try again.\n")
            query_seq = input("Enter sequence: ").strip()
            continue
        
        # Clean sequence
        query_clean = re.sub(r"[^A-Za-z]", "", query_seq)
        print(f"  ✓ Sequence length: {len(query_clean)} aa")
        
        # Offer to save
        save = input("\nSave this query for future use? (y/n): ").strip().lower()
        if save == 'y':
            os.makedirs(saved_queries_dir, exist_ok=True)
            name = input("Name for this query (e.g., 'my_antibody'): ").strip()
            if not name:
                name = f"query_{len(os.listdir(saved_queries_dir))+1}"
            if not name.endswith(('.txt', '.fasta')):
                name += '.txt'
            save_path = os.path.join(saved_queries_dir, name)
            with open(save_path, 'w') as f:
                f.write(query_clean)
            print(f"  ✓ Saved to {save_path}")
        
        params['query_seq'] = query_seq
        break
    
    # Get search regions
    print("\nStep 3: Search Regions")
    print("-" * 40)
    print("What regions do you want to search?")
    print("  1. Full sequence")
    print("  2. CDRs only (CDR1, CDR2, CDR3)")
    print("  3. Individual CDRs (choose specific)")
    print("  4. Frameworks only")
    print("  5. Custom combination")
    
    while True:
        choice = input("\nYour choice (1-5): ").strip()
        
        if choice == '1':
            params['regions'] = 'full'
            break
        elif choice == '2':
            params['regions'] = 'cdr1,cdr2,cdr3'
            break
        elif choice == '3':
            print("\n  Select CDRs to search (comma-separated):")
            print("    Options: cdr1, cdr2, cdr3")
            cdrs = input("  Enter choices: ").strip().lower()
            params['regions'] = cdrs
            break
        elif choice == '4':
            params['regions'] = 'frameworks'
            break
        elif choice == '5':
            print("\n  Enter custom regions (comma-separated):")
            print("    Options: full, cdr1, cdr2, cdr3, frameworks")
            custom = input("  Enter choices: ").strip().lower()
            params['regions'] = custom
            break
        else:
            print("  ⚠ Invalid choice. Please enter 1-5.")
    
    # Ask about identity thresholds
    print("\nStep 4: Identity Thresholds")
    print("-" * 40)
    print("Do you want to set minimum identity thresholds?")
    print("  1. No thresholds (return all matches)")
    print("  2. Set specific thresholds manually")
    print("  3. Use automatic threshold detection (recommended)")
    
    while True:
        choice = input("\nYour choice (1-3): ").strip()
        
        if choice == '1':
            # No thresholds
            params['use_thresholds'] = False
            break
        elif choice == '2':
            # Manual thresholds
            params['use_thresholds'] = True
            regions_list = params['regions'].split(',')
            for region in regions_list:
                region = region.strip()
                while True:
                    threshold = input(f"  Minimum identity for {region.upper()} (0.0-1.0, or Enter to skip): ").strip()
                    if not threshold:
                        break
                    try:
                        val = float(threshold)
                        if 0.0 <= val <= 1.0:
                            params[f'min_id_{region}'] = val
                            break
                        else:
                            print("    ⚠ Please enter a value between 0.0 and 1.0")
                    except ValueError:
                        print("    ⚠ Please enter a valid number")
            break
        elif choice == '3':
            params['use_thresholds'] = True
            params['run_analysis'] = True
            
            # Ask about detailed test run output
            print("\nAutomatic threshold detection will analyze your database.")
            show_analysis = input("Show detailed analysis? (y/n): ").strip().lower()
            params['show_analysis'] = (show_analysis == 'y')
            
            if params['show_analysis']:
                # Get sample size and target
                # Suggest sample size based on database size
                if 'camel' in params.get('db_root', '').lower():
                    default_sample = 10000  # Much larger for small databases
                    suggested = "(default 10000 for Camel)"
                elif 'rabbit' in params.get('db_root', '').lower() or 'rat' in params.get('db_root', '').lower():
                    default_sample = 5000
                    suggested = "(default 5000 for small database)"
                else:
                    default_sample = 2000  # Larger default for bigger databases too
                    suggested = "(default 2000)"
                
                while True:
                    sample_size = input(f"  Sample size per shard {suggested}: ").strip()
                    if not sample_size:
                        params['sample_size'] = default_sample
                        break
                    try:
                        size = int(sample_size)
                        if size > 0:
                            params['sample_size'] = size
                            break
                        else:
                            print("    ⚠ Sample size must be positive")
                    except ValueError:
                        print("    ⚠ Please enter a valid number")
                
                while True:
                    # Suggest target based on database
                    if 'camel' in params.get('db_root', '').lower():
                        default_target = 1000
                        suggested = "(default 1000 for Camel)"
                    elif 'rabbit' in params.get('db_root', '').lower() or 'rat' in params.get('db_root', '').lower():
                        default_target = 2000
                        suggested = "(default 2000 for small database)"
                    else:
                        default_target = 10000
                        suggested = "(default 10000)"
                    
                    target = input(f"  Target number of hits {suggested}: ").strip()
                    if not target:
                        params['target_hits'] = default_target
                        break
                    try:
                        hits = int(target)
                        if hits > 0:
                            params['target_hits'] = hits
                            break
                        else:
                            print("    ⚠ Target hits must be positive")
                    except ValueError:
                        print("    ⚠ Please enter a valid number")
            else:
                # Use defaults silently - but smart defaults based on database
                if 'camel' in params.get('db_root', '').lower():
                    params['sample_size'] = 10000  # Much larger for Camel
                    params['target_hits'] = 1000   # More realistic target
                    print("  Using defaults: 10,000 samples/shard (Camel database), targeting ~1,000 hits")
                elif 'rabbit' in params.get('db_root', '').lower() or 'rat' in params.get('db_root', '').lower():
                    params['sample_size'] = 5000
                    params['target_hits'] = 2000
                    print("  Using defaults: 5,000 samples/shard, targeting ~2,000 hits")
                else:
                    params['sample_size'] = 2000
                    params['target_hits'] = 10000
                    print("  Using defaults: 2,000 samples/shard, targeting ~10,000 hits")
            break
        else:
            print("  ⚠ Invalid choice. Please enter 1-3.")
    
    # Output options
    print("\nStep 5: Output Options")
    print("-" * 40)
    
    outdir = input("Output directory (default: ./results): ").strip()
    params['outdir'] = outdir if outdir else './results'
    
    tag = input("Optional tag for output files: ").strip()
    if tag:
        params['tag'] = tag
    
    # Advanced options
    print("\nStep 6: Advanced Options")
    print("-" * 40)
    advanced = input("Configure advanced options? (y/n): ").strip().lower()
    
    # Check if searching full sequence
    searching_full = 'full' in params.get('regions', '')
    
    if advanced == 'y':
        # Length filtering - only show for non-full searches
        if not searching_full:
            print("\n  Length Filtering:")
            for cdr in ['cdr1', 'cdr2', 'cdr3']:
                if cdr in params.get('regions', ''):
                    window = input(f"    {cdr.upper()} length window (±aa, or Enter to skip): ").strip()
                    if window:
                        try:
                            val = int(window)
                            if val >= 0:
                                params[f'{cdr}_len_window'] = val
                        except ValueError:
                            pass
        
        # Other options
        print("\n  Other Options:")
        
        scheme = input("    Numbering scheme (imgt/kabat, default imgt): ").strip().lower()
        if scheme in ['imgt', 'kabat']:
            params['scheme'] = scheme
        else:
            params['scheme'] = 'imgt'
        
        # Only show sequence length filter if searching full
        if searching_full:
            seq_filter = input("    Enable sequence length pre-filter? (y/n, default y): ").strip().lower()
            params['disable_seq_length_prefilter'] = (seq_filter == 'n')
        else:
            params['disable_seq_length_prefilter'] = True  # Disable for CDR-only searches
        
        # Only show CDR length filter if NOT searching full
        if not searching_full:
            cdr_filter = input("    Enable CDR length pre-filter? (y/n, default y): ").strip().lower()
            params['disable_cdr_length_prefilter'] = (cdr_filter == 'n')
        else:
            params['disable_cdr_length_prefilter'] = True  # Disable for full searches
    else:
        # Show defaults being used
        print("\n  Using default advanced options:")
        if not searching_full:
            print("    • No CDR length windows (accepting all lengths)")
        print("    • Numbering scheme: IMGT")
        if searching_full:
            print("    • Sequence length pre-filter: ENABLED")
        if not searching_full:
            print("    • CDR length pre-filter: ENABLED")
        params['scheme'] = 'imgt'
        params['disable_seq_length_prefilter'] = not searching_full
        params['disable_cdr_length_prefilter'] = False
    
    return params

def convert_params_to_args(params):
    """Convert interactive params to argparse namespace."""
    args = argparse.Namespace()
    
    # Required params
    args.db_root = params.get('db_root')
    args.query_seq = params.get('query_seq')
    args.regions = params.get('regions', 'full')
    args.outdir = params.get('outdir', './results')
    
    # Optional params
    args.tag = params.get('tag', None)
    args.scheme = params.get('scheme', 'imgt')
    
    # Thresholds
    args.min_id_full = params.get('min_id_full', None)
    args.min_id_cdr1 = params.get('min_id_cdr1', None)
    args.min_id_cdr2 = params.get('min_id_cdr2', None)
    args.min_id_cdr3 = params.get('min_id_cdr3', None)
    args.min_id_frameworks = params.get('min_id_frameworks', None)
    
    # Length windows
    args.cdr1_len_window = params.get('cdr1_len_window', None)
    args.cdr2_len_window = params.get('cdr2_len_window', None)
    args.cdr3_len_window = params.get('cdr3_len_window', None)
    
    # Filters
    args.disable_seq_length_prefilter = params.get('disable_seq_length_prefilter', False)
    args.disable_cdr_length_prefilter = params.get('disable_cdr_length_prefilter', False)
    args.seq_length_threshold = params.get('seq_length_threshold', 0.25)
    
    # Analysis params
    args.sample_size = params.get('sample_size', 1000)
    args.target_hits = params.get('target_hits', 10000)
    
    return args

# ============================================================================
# SAMPLING AND ANALYSIS FUNCTIONS
# ============================================================================

def sample_and_analyze_shard(npz_file: str, query_cdrs: Dict, sample_size: int, 
                            regions: List[str]) -> pd.DataFrame:
    """Sample sequences from a shard and calculate identities."""
    try:
        data = np.load(npz_file, allow_pickle=True)
        
        # Look for the "numberings" key specifically (OAS database format)
        if "numberings" not in data:
            # Fallback to other possible keys
            if "arr_0" in data.files:
                arr = data["arr_0"]
            elif "sequences" in data.files:
                arr = data["sequences"]
            else:
                data.close()
                return pd.DataFrame()
        else:
            arr = data["numberings"]
        
        # Random sample
        n_seqs = len(arr)
        if n_seqs == 0:
            data.close()
            return pd.DataFrame()
            
        if n_seqs <= sample_size:
            indices = list(range(n_seqs))
        else:
            indices = random.sample(range(n_seqs), sample_size)
        
        results = []
        for idx in indices:
            row = arr[idx]
            result = {"shard": os.path.basename(npz_file), "index": idx}
            
            # Calculate identities for each region
            identities = []
            for region in regions:
                if region == "cdr1":
                    db_seq = extract_cdr1_from_npz(row)
                    q_seq = query_cdrs.get("cdr1", "")
                elif region == "cdr2":
                    db_seq = extract_cdr2_from_npz(row)
                    q_seq = query_cdrs.get("cdr2", "")
                elif region == "cdr3":
                    db_seq = extract_cdr3_from_npz(row)
                    q_seq = query_cdrs.get("cdr3", "")
                elif region == "frameworks":
                    db_seq = extract_frameworks_from_npz(row)
                    q_seq = query_cdrs.get("frameworks", "")
                elif region == "full":
                    db_seq = extract_whole_from_npz(row)
                    q_seq = query_cdrs.get("full", "")
                else:
                    continue
                
                if db_seq and q_seq:
                    identity = calc_id_by_region(q_seq, db_seq)
                    result[f"id_{region}"] = identity
                    identities.append(identity)
            
            if identities:
                result["avg_id"] = np.mean(identities)
                results.append(result)
        
        data.close()
        return pd.DataFrame(results)
    
    except Exception as e:
        print(f"Error sampling {npz_file}: {e}")
        return pd.DataFrame()

def run_interactive_analysis(db_root: str, query_seq: str, regions: List[str],
                            sample_size: int = 1000, target_hits: int = 10000) -> Dict:
    """Run interactive analysis to determine optimal thresholds."""
    
    print("\n" + "="*80)
    print("INTERACTIVE THRESHOLD ANALYSIS")
    print("="*80)
    print(f"Sampling {sample_size} sequences per shard...")
    print(f"Target: Find thresholds to get ~{target_hits:,} hits")
    print("="*80 + "\n")
    
    # Extract query CDRs
    q_full = re.sub(r"[^A-Za-z]", "", query_seq.upper())
    query_cdrs = {"full": q_full}
    
    if set(regions) != {"full"}:
        qc = extract_cdrs_from_query_sequence(q_full)
        query_cdrs.update(qc)
    
    if 'frameworks' in regions:
        query_cdrs['frameworks'] = q_full  # Placeholder
    
    # Find shards
    npz_files = sorted(glob.glob(os.path.join(db_root, "*.npz")))
    
    if not npz_files:
        print("ERROR: No NPZ files found!")
        return {}
    
    print(f"Found {len(npz_files)} shards")
    
    # Estimate database size based on number of shards
    # More accurate estimates based on typical OAS database sizes
    if 'camel' in db_root.lower():
        # Camel has 2 shards with ~50k each
        avg_seqs_per_shard = 50_000
    elif 'human' in db_root.lower():
        # Human has 379 shards with varying sizes
        avg_seqs_per_shard = 500_000
    elif 'mouse' in db_root.lower():
        # Mouse has 35 shards
        avg_seqs_per_shard = 300_000
    elif len(npz_files) <= 3:  # Small databases (Rabbit, Rat, Rhesus)
        avg_seqs_per_shard = 100_000
    elif len(npz_files) > 100:  # Large database
        avg_seqs_per_shard = 500_000
    else:  # Medium database
        avg_seqs_per_shard = 200_000
    
    estimated_db_size = len(npz_files) * avg_seqs_per_shard
    print(f"Estimated database size: ~{estimated_db_size:,} sequences")
    
    # Limit sample size to not exceed estimated shard size
    effective_sample_size = min(sample_size, int(avg_seqs_per_shard * 0.5))  # Sample up to 50% of a shard
    if effective_sample_size < sample_size:
        print(f"Adjusted sample size to {effective_sample_size:,} (50% of estimated shard size)")
        sample_size = effective_sample_size
    
    # Limit to first 10 shards for very large databases
    shards_to_sample = min(len(npz_files), 10)
    print(f"Sampling from {shards_to_sample} shards")
    print(f"Total sequences to sample: {shards_to_sample * sample_size:,}\n")
    
    # Sample from shards
    all_samples = []
    # Sample from limited number of shards
    sample_bar = tqdm(npz_files[:shards_to_sample], desc="Sampling shards")
    
    for npz_file in sample_bar:
        df = sample_and_analyze_shard(npz_file, query_cdrs, sample_size, regions)
        if not df.empty:
            all_samples.append(df)
    
    if not all_samples:
        print("ERROR: No samples collected!")
        return {}
    
    combined = pd.concat(all_samples, ignore_index=True)
    
    # Analyze distributions
    print(f"\n{'='*80}")
    print("IDENTITY DISTRIBUTION ANALYSIS")
    print(f"{'='*80}\n")
    
    recommendations = {}
    
    # Use the estimated database size from earlier
    # (already calculated above)
    
    for region in regions:
        col = f"id_{region}"
        if col not in combined.columns:
            continue
        
        values = combined[col].dropna().values  # Remove NaN values
        
        if len(values) == 0:
            print(f"\n{region.upper()} Statistics:")
            print(f"  No valid values extracted for {region.upper()}")
            continue
        
        print(f"\n{region.upper()} Statistics:")
        print(f"  Samples with data: {len(values)}/{len(combined)} ({len(values)*100/len(combined):.1f}%)")
        print(f"  Max:    {values.max():.3f} ({values.max()*100:.1f}%)")
        print(f"  Q95:    {np.percentile(values, 95):.3f} ({np.percentile(values, 95)*100:.1f}%)")
        print(f"  Q90:    {np.percentile(values, 90):.3f} ({np.percentile(values, 90)*100:.1f}%)")
        print(f"  Q75:    {np.percentile(values, 75):.3f} ({np.percentile(values, 75)*100:.1f}%)")
        print(f"  Median: {np.median(values):.3f} ({np.median(values)*100:.1f}%)")
        print(f"  Q25:    {np.percentile(values, 25):.3f} ({np.percentile(values, 25)*100:.1f}%)")
        print(f"  Mean:   {values.mean():.3f} ({values.mean()*100:.1f}%)")
        
        # Estimate threshold for target hits
        total_sampled = len(combined)
        # Use the database size estimate from earlier
        scale_factor = estimated_db_size / total_sampled
        needed_in_sample = min(target_hits / scale_factor, len(values))
        
        if needed_in_sample >= len(values):
            # Not enough good matches to reach target
            print(f"\n  ⚠ Warning: Based on sampling, database may not have {target_hits:,} matches")
            print(f"  at any reasonable threshold for {region.upper()}")
            
            # Suggest a more reasonable threshold (e.g., Q90 or Q95)
            if len(values) > 10:
                recommended_threshold = np.percentile(values, 90)  # Top 10%
                percentile_rank = 90
                estimated_hits = len(values) * scale_factor * 0.1  # Rough estimate
                print(f"\n  SUGGESTED THRESHOLD: {recommended_threshold:.3f} ({recommended_threshold*100:.1f}%)")
                print(f"  This is the 90th percentile (top 10% of matches)")
                print(f"  Estimated hits: ~{int(estimated_hits):,}")
            else:
                # Very few matches found
                recommended_threshold = np.percentile(values, 50)  # Median
                percentile_rank = 50
                estimated_hits = len(values) * scale_factor * 0.5
                print(f"\n  SUGGESTED THRESHOLD: {recommended_threshold:.3f} ({recommended_threshold*100:.1f}%)")
                print(f"  Using median due to limited matches")
                print(f"  Estimated hits: ~{int(estimated_hits):,}")
            
            recommendations[region] = recommended_threshold
        else:
            # Normal case - we have enough samples
            sorted_vals = np.sort(values)[::-1]  # Sort descending
            threshold_idx = int(needed_in_sample)
            recommended_threshold = sorted_vals[threshold_idx]
            
            # Calculate correct percentile
            # If we're at position threshold_idx in a descending sort,
            # the percentile is (1 - threshold_idx/len) * 100
            percentile_rank = (1 - (threshold_idx / len(values))) * 100
            
            print(f"\n  RECOMMENDED THRESHOLD: {recommended_threshold:.3f} ({recommended_threshold*100:.1f}%)")
            print(f"  This is the {percentile_rank:.1f}th percentile (top {100-percentile_rank:.1f}% of matches)")
            print(f"  Estimated to give ~{target_hits:,} hits across full database")
            
            recommendations[region] = recommended_threshold
    
    # Show top matches
    print(f"\n{'='*80}")
    print("TOP 10 MATCHES IN SAMPLE")
    print(f"{'='*80}\n")
    
    if 'avg_id' in combined.columns:
        top10 = combined.nlargest(10, 'avg_id')
        for i, (_, row) in enumerate(top10.iterrows(), 1):
            print(f"#{i:2d} Shard: {row['shard'][:40]:40s} (idx {row['index']:,})")
            for region in regions:
                col = f"id_{region}"
                if col in row:
                    print(f"    {region:10s}: {row[col]*100:5.1f}%")
            print(f"    Average:    {row['avg_id']*100:5.1f}%")
            print()
    
    return recommendations

# ============================================================================
# FULL SCAN FUNCTIONS (from v6)
# ============================================================================

def write_metadata_block(file_handle, args, query_cdrs, start_time, config):
    """Write comprehensive metadata block to CSV file."""
    file_handle.write("# ============================================\n")
    file_handle.write("# NPZ FULLSCAN v6 INTERACTIVE - RUN METADATA\n")
    file_handle.write("# ============================================\n")
    file_handle.write(f"# Run started: {start_time.strftime('%Y-%m-%d %H:%M:%S')}\n")
    file_handle.write(f"# Script: {os.path.basename(__file__)}\n")
    file_handle.write(f"# Command: {' '.join(sys.argv)}\n")
    file_handle.write(f"# Platform: {platform.platform()}\n")
    file_handle.write(f"# Python: {sys.version.split()[0]}\n")
    file_handle.write("#\n")
    file_handle.write("# DATABASE INFORMATION\n")
    file_handle.write(f"# Database path: {args.db_root}\n")
    file_handle.write(f"# Species: {get_species_from_path(args.db_root)}\n")
    file_handle.write("#\n")
    file_handle.write("# QUERY SEQUENCE\n")
    file_handle.write(f"# Length: {len(query_cdrs['full'])} aa\n")
    file_handle.write(f"# Sequence: {query_cdrs['full']}\n")
    file_handle.write(f"# CDR1: {query_cdrs.get('cdr1', 'N/A')} (len={len(query_cdrs.get('cdr1', ''))})\n")
    file_handle.write(f"# CDR2: {query_cdrs.get('cdr2', 'N/A')} (len={len(query_cdrs.get('cdr2', ''))})\n")
    file_handle.write(f"# CDR3: {query_cdrs.get('cdr3', 'N/A')} (len={len(query_cdrs.get('cdr3', ''))})\n")
    file_handle.write("#\n")
    file_handle.write("# SEARCH PARAMETERS\n")
    file_handle.write(f"# Regions searched: {args.regions}\n")
    file_handle.write(f"# Numbering scheme: {args.scheme}\n")
    
    # Write thresholds
    for key, value in config.items():
        if key.startswith('min_id_') and value is not None:
            region = key.replace('min_id_', '').upper()
            file_handle.write(f"# {region} min identity: {value}\n")
    
    if config.get('use_length_filter'):
        file_handle.write("#\n# LENGTH FILTERING\n")
        for cdr, window in config['length_windows'].items():
            file_handle.write(f"# {cdr.upper()} length window: ±{window} aa\n")
    
    if config.get('enable_seq_length_prefilter'):
        file_handle.write("#\n# PRE-FILTERING\n")
        file_handle.write("# Sequence length pre-filter: ENABLED\n")
        file_handle.write(f"# Sequence length threshold: {config.get('seq_length_threshold', 0.25)}\n")
    if config.get('enable_cdr_length_prefilter'):
        file_handle.write("# CDR length pre-filter: ENABLED\n")
    
    if config.get('interactive_analysis_run'):
        file_handle.write("#\n# INTERACTIVE ANALYSIS\n")
        file_handle.write("# Interactive analysis was run: YES\n")
        file_handle.write(f"# Used recommended thresholds: {'YES' if config.get('used_recommendations') else 'NO'}\n")
    
    file_handle.write("# ============================================\n\n")

def get_species_from_path(db_path: str) -> str:
    """Extract species from database path."""
    path_parts = db_path.rstrip('/').split('/')
    for part in reversed(path_parts):
        if part.lower() in ['human', 'mouse', 'rabbit', 'camel', 'llama', 'alpaca', 
                           'rat', 'chicken', 'rhesus', 'cynomolgus']:
            return part.capitalize()
    return "Unknown"

def create_output_structure(outdir: str, species: str, regions: str, 
                           thresholds: Dict, tag: Optional[str] = None) -> Tuple[str, str, str]:
    """Create organized output folder and filename."""
    # Create base output directory
    os.makedirs(outdir, exist_ok=True)
    
    # Create dated folder
    date_str = dt.datetime.now().strftime("%Y%m%d")
    folder_name = f"{date_str}_{species.lower()}"
    output_folder = os.path.join(outdir, folder_name)
    os.makedirs(output_folder, exist_ok=True)
    
    # Create CSVs subfolder
    csv_folder = os.path.join(output_folder, "csvs")
    os.makedirs(csv_folder, exist_ok=True)
    
    # Build filename with parameters
    filename_parts = [
        "npz_scan",
        species.lower(),
        regions.replace(',', '_')
    ]
    
    # Add thresholds to filename
    for region, threshold in thresholds.items():
        if threshold is not None:
            region_name = region.replace('min_id_', '')
            filename_parts.append(f"{region_name}{int(threshold*100)}")
    
    # Add tag if provided
    if tag:
        filename_parts.append(tag)
    
    # Add timestamp
    time_str = dt.datetime.now().strftime("%H%M%S")
    filename_parts.append(time_str)
    
    filename = "_".join(filename_parts) + ".csv"
    csv_path = os.path.join(csv_folder, filename)  # CSVs go in subfolder
    
    return csv_path, output_folder, csv_folder

def process_npz_file_fullscan(npz_path: str, config: Dict, query_cdrs: Dict, 
                             shard_idx: int, total_shards: int) -> Dict:
    """Process a single NPZ file for full scan."""
    shard_name = os.path.basename(npz_path)
    start_time = time.time()
    
    try:
        data = np.load(npz_path, allow_pickle=True)
        
        # Look for the "numberings" key specifically (OAS database format)
        if "numberings" not in data:
            # Fallback to other possible keys
            if "arr_0" in data.files:
                arr = data["arr_0"]
            elif "sequences" in data.files:
                arr = data["sequences"]
            else:
                data.close()
                return {
                    "shard": shard_name,
                    "n_total": 0,
                    "n_hits": 0,
                    "n_skipped": 0,
                    "n_filtered_length": 0,
                    "n_filtered_seq_length": 0,
                    "n_filtered_cdr_length": 0,
                    "max_id": 0.0,
                    "min_id": 0.0,
                    "q75": 0.0,
                    "best_seq": "",
                    "worst_seq": "",
                    "time_s": 0.0,
                    "hits": []  # Add hits list
                }
        else:
            arr = data["numberings"]
        
        n_seqs = len(arr)
        regions = config["regions"].split(",")
        
        # Initialize counters
        n_hits = 0
        n_skipped = 0
        n_filtered_length = 0
        n_filtered_seq_length = 0
        n_filtered_cdr_length = 0
        
        identities = []
        best_id = 0.0
        worst_id = 1.0
        best_seq = ""
        worst_seq = ""
        hits_data = []  # Store actual hit data
        
        # Process sequences with progress bar
        inner_bar = tqdm(range(n_seqs), 
                        desc=f"  Shard {shard_idx}/{total_shards}",
                        ncols=100, position=1, leave=False)
        
        for i in inner_bar:
            row = arr[i]
            
            # Sequence length pre-filter for full search
            if config.get("enable_seq_length_prefilter") and "full" in regions:
                db_full = extract_whole_from_npz(row)
                if not db_full:
                    n_skipped += 1
                    continue
                
                len_ratio = abs(len(db_full) - len(query_cdrs["full"])) / len(query_cdrs["full"])
                if len_ratio > config.get("seq_length_threshold", 0.25):
                    n_filtered_seq_length += 1
                    continue
            
            # CDR length pre-filter
            if config.get("enable_cdr_length_prefilter"):
                skip_cdr = False
                for cdr in ["cdr1", "cdr2", "cdr3"]:
                    if cdr not in regions:
                        continue
                    
                    if cdr == "cdr1":
                        db_cdr = extract_cdr1_from_npz(row)
                    elif cdr == "cdr2":
                        db_cdr = extract_cdr2_from_npz(row)
                    else:
                        db_cdr = extract_cdr3_from_npz(row)
                    
                    q_cdr_len = config["query_lengths"].get(cdr, 0)
                    if q_cdr_len > 0:
                        len_diff = abs(len(db_cdr) - q_cdr_len)
                        if len_diff > 3:  # Hard-coded threshold
                            skip_cdr = True
                            break
                
                if skip_cdr:
                    n_filtered_cdr_length += 1
                    continue
            
            # Calculate identities and extract sequences
            region_identities = []
            region_sequences = {}
            passes_threshold = True
            
            # Always extract full sequence and ALL CDRs for results
            full_seq = extract_whole_from_npz(row)
            
            # Always extract all CDRs for reporting
            db_cdr1 = extract_cdr1_from_npz(row)
            db_cdr2 = extract_cdr2_from_npz(row)
            db_cdr3 = extract_cdr3_from_npz(row)
            
            # Store all CDRs
            all_cdrs = {
                'cdr1': db_cdr1,
                'cdr2': db_cdr2,
                'cdr3': db_cdr3
            }
            
            # Calculate identities for ALL CDRs (for reporting)
            all_identities = {}
            if db_cdr1 and query_cdrs.get('cdr1'):
                all_identities['cdr1'] = calc_id_by_region(query_cdrs['cdr1'], db_cdr1)
            else:
                all_identities['cdr1'] = 0.0
                
            if db_cdr2 and query_cdrs.get('cdr2'):
                all_identities['cdr2'] = calc_id_by_region(query_cdrs['cdr2'], db_cdr2)
            else:
                all_identities['cdr2'] = 0.0
                
            if db_cdr3 and query_cdrs.get('cdr3'):
                all_identities['cdr3'] = calc_id_by_region(query_cdrs['cdr3'], db_cdr3)
            else:
                all_identities['cdr3'] = 0.0
            
            # Calculate full sequence identity if we have it
            if full_seq and query_cdrs.get('full'):
                all_identities['full'] = calc_id_by_region(query_cdrs['full'], full_seq)
            else:
                all_identities['full'] = 0.0
            
            # Now check thresholds only for SEARCHED regions
            for region in regions:
                region = region.strip()
                
                # Get sequences for searched regions
                if region == "cdr1":
                    db_seq = db_cdr1
                    q_seq = query_cdrs.get("cdr1", "")
                elif region == "cdr2":
                    db_seq = db_cdr2
                    q_seq = query_cdrs.get("cdr2", "")
                elif region == "cdr3":
                    db_seq = db_cdr3
                    q_seq = query_cdrs.get("cdr3", "")
                elif region == "frameworks":
                    db_seq = extract_frameworks_from_npz(row)
                    q_seq = query_cdrs.get("frameworks", "")
                    region_sequences['frameworks'] = db_seq
                    if db_seq and q_seq:
                        all_identities['frameworks'] = calc_id_by_region(q_seq, db_seq)
                elif region == "full":
                    db_seq = full_seq
                    q_seq = query_cdrs.get("full", "")
                else:
                    continue
                
                if not db_seq or not q_seq:
                    passes_threshold = False
                    break
                
                # Length filter for searched regions
                if config.get("use_length_filter"):
                    window = config["length_windows"].get(region, None)
                    if window is not None:
                        if abs(len(db_seq) - len(q_seq)) > window:
                            n_filtered_length += 1
                            passes_threshold = False
                            break
                
                # Get identity for this searched region
                identity = all_identities.get(region, 0.0)
                region_identities.append(identity)
                
                # Check threshold for searched regions only
                min_id = config.get(f"min_id_{region}")
                if min_id is not None and identity < min_id:
                    passes_threshold = False
                    break
            
            if passes_threshold and region_identities:
                avg_identity = np.mean(region_identities)  # Average of SEARCHED regions only
                n_hits += 1
                identities.append(avg_identity)
                
                # Store hit data with ALL CDRs and identities
                hit_data = {
                    "shard": shard_name,
                    "index": i,
                    "avg_identity": avg_identity,
                    "full_sequence": full_seq,
                    "full_identity": all_identities.get('full', 0.0),
                    # Always include all CDRs
                    "cdr1_seq": all_cdrs['cdr1'],
                    "cdr1_identity": all_identities['cdr1'],
                    "cdr2_seq": all_cdrs['cdr2'],
                    "cdr2_identity": all_identities['cdr2'],
                    "cdr3_seq": all_cdrs['cdr3'],
                    "cdr3_identity": all_identities['cdr3']
                }
                
                # Add frameworks if it was extracted
                if 'frameworks' in region_sequences:
                    hit_data["frameworks_seq"] = region_sequences['frameworks']
                    hit_data["frameworks_identity"] = all_identities.get('frameworks', 0.0)
                
                hits_data.append(hit_data)
                
                if avg_identity > best_id:
                    best_id = avg_identity
                    if "full" in regions:
                        best_seq = full_seq[:50]
                    else:
                        best_seq = f"CDR3: {region_sequences.get('cdr3', 'N/A')}"
                
                if avg_identity < worst_id:
                    worst_id = avg_identity
                    if "full" in regions:
                        worst_seq = full_seq[:50]
                    else:
                        worst_seq = f"CDR3: {region_sequences.get('cdr3', 'N/A')}"
            
            # Update progress
            if (i + 1) % 1000 == 0:
                inner_bar.set_postfix({"hits": n_hits})
        
        inner_bar.close()
        data.close()
        
        # Calculate statistics
        q75 = np.percentile(identities, 75) if identities else 0.0
        
        return {
            "shard": shard_name,
            "n_total": n_seqs,
            "n_hits": n_hits,
            "n_skipped": n_skipped,
            "n_filtered_length": n_filtered_length,
            "n_filtered_seq_length": n_filtered_seq_length,
            "n_filtered_cdr_length": n_filtered_cdr_length,
            "max_id": best_id,
            "min_id": worst_id if n_hits > 0 else 0.0,
            "q75": q75,
            "best_seq": best_seq,
            "worst_seq": worst_seq,
            "time_s": time.time() - start_time,
            "hits": hits_data  # Return actual hit data
        }
    
    except Exception as e:
        print(f"\nError processing {shard_name}: {e}")
        return {
            "shard": shard_name,
            "n_total": 0,
            "n_hits": 0,
            "n_skipped": 0,
            "n_filtered_length": 0,
            "n_filtered_seq_length": 0,
            "n_filtered_cdr_length": 0,
            "max_id": 0.0,
            "min_id": 0.0,
            "q75": 0.0,
            "best_seq": "",
            "worst_seq": "",
            "time_s": 0.0,
            "hits": []
        }
    """Process a single NPZ file for full scan."""
    shard_name = os.path.basename(npz_path)
    start_time = time.time()
    
    try:
        data = np.load(npz_path, allow_pickle=True)
        
        # Look for the "numberings" key specifically (OAS database format)
        if "numberings" not in data:
            # Fallback to other possible keys
            if "arr_0" in data.files:
                arr = data["arr_0"]
            elif "sequences" in data.files:
                arr = data["sequences"]
            else:
                data.close()
                return {
                    "shard": shard_name,
                    "n_total": 0,
                    "n_hits": 0,
                    "n_skipped": 0,
                    "n_filtered_length": 0,
                    "n_filtered_seq_length": 0,
                    "n_filtered_cdr_length": 0,
                    "max_id": 0.0,
                    "min_id": 0.0,
                    "q75": 0.0,
                    "best_seq": "",
                    "worst_seq": "",
                    "time_s": 0.0
                }
        else:
            arr = data["numberings"]
        
        n_seqs = len(arr)
        regions = config["regions"].split(",")
        
        # Initialize counters
        n_hits = 0
        n_skipped = 0
        n_filtered_length = 0
        n_filtered_seq_length = 0
        n_filtered_cdr_length = 0
        
        identities = []
        best_id = 0.0
        worst_id = 1.0
        best_seq = ""
        worst_seq = ""
        
        # Process sequences with progress bar
        inner_bar = tqdm(range(n_seqs), 
                        desc=f"  Shard {shard_idx}/{total_shards}",
                        ncols=100, position=1, leave=False)
        
        for i in inner_bar:
            row = arr[i]
            
            # Sequence length pre-filter for full search
            if config.get("enable_seq_length_prefilter") and "full" in regions:
                db_full = extract_whole_from_npz(row)
                if not db_full:
                    n_skipped += 1
                    continue
                
                len_ratio = abs(len(db_full) - len(query_cdrs["full"])) / len(query_cdrs["full"])
                if len_ratio > config.get("seq_length_threshold", 0.25):
                    n_filtered_seq_length += 1
                    continue
            
            # CDR length pre-filter
            if config.get("enable_cdr_length_prefilter"):
                skip_cdr = False
                for cdr in ["cdr1", "cdr2", "cdr3"]:
                    if cdr not in regions:
                        continue
                    
                    if cdr == "cdr1":
                        db_cdr = extract_cdr1_from_npz(row)
                    elif cdr == "cdr2":
                        db_cdr = extract_cdr2_from_npz(row)
                    else:
                        db_cdr = extract_cdr3_from_npz(row)
                    
                    q_cdr_len = config["query_lengths"].get(cdr, 0)
                    if q_cdr_len > 0:
                        len_diff = abs(len(db_cdr) - q_cdr_len)
                        if len_diff > 3:  # Hard-coded threshold
                            skip_cdr = True
                            break
                
                if skip_cdr:
                    n_filtered_cdr_length += 1
                    continue
            
            # Calculate identities
            region_identities = []
            passes_threshold = True
            
            for region in regions:
                region = region.strip()
                
                # Extract sequences
                if region == "cdr1":
                    db_seq = extract_cdr1_from_npz(row)
                    q_seq = query_cdrs.get("cdr1", "")
                elif region == "cdr2":
                    db_seq = extract_cdr2_from_npz(row)
                    q_seq = query_cdrs.get("cdr2", "")
                elif region == "cdr3":
                    db_seq = extract_cdr3_from_npz(row)
                    q_seq = query_cdrs.get("cdr3", "")
                elif region == "frameworks":
                    db_seq = extract_frameworks_from_npz(row)
                    q_seq = query_cdrs.get("frameworks", "")
                elif region == "full":
                    db_seq = extract_whole_from_npz(row)
                    q_seq = query_cdrs.get("full", "")
                else:
                    continue
                
                if not db_seq or not q_seq:
                    passes_threshold = False
                    break
                
                # Length filter
                if config.get("use_length_filter"):
                    window = config["length_windows"].get(region, None)
                    if window is not None:
                        if abs(len(db_seq) - len(q_seq)) > window:
                            n_filtered_length += 1
                            passes_threshold = False
                            break
                
                # Calculate identity using edit distance
                identity = calc_id_by_region(q_seq, db_seq)
                region_identities.append(identity)
                
                # Check threshold
                min_id = config.get(f"min_id_{region}")
                if min_id is not None and identity < min_id:
                    passes_threshold = False
                    break
            
            if passes_threshold and region_identities:
                avg_identity = np.mean(region_identities)
                n_hits += 1
                identities.append(avg_identity)
                
                if avg_identity > best_id:
                    best_id = avg_identity
                    if "full" in regions:
                        best_seq = extract_whole_from_npz(row)[:50]
                    else:
                        best_seq = f"CDR3: {extract_cdr3_from_npz(row)}"
                
                if avg_identity < worst_id:
                    worst_id = avg_identity
                    if "full" in regions:
                        worst_seq = extract_whole_from_npz(row)[:50]
                    else:
                        worst_seq = f"CDR3: {extract_cdr3_from_npz(row)}"
            
            # Update progress
            if (i + 1) % 1000 == 0:
                inner_bar.set_postfix({"hits": n_hits})
        
        inner_bar.close()
        data.close()
        
        # Calculate statistics
        q75 = np.percentile(identities, 75) if identities else 0.0
        
        return {
            "shard": shard_name,
            "n_total": n_seqs,
            "n_hits": n_hits,
            "n_skipped": n_skipped,
            "n_filtered_length": n_filtered_length,
            "n_filtered_seq_length": n_filtered_seq_length,
            "n_filtered_cdr_length": n_filtered_cdr_length,
            "max_id": best_id,
            "min_id": worst_id if n_hits > 0 else 0.0,
            "q75": q75,
            "best_seq": best_seq,
            "worst_seq": worst_seq,
            "time_s": time.time() - start_time
        }
    
    except Exception as e:
        print(f"\nError processing {shard_name}: {e}")
        return {
            "shard": shard_name,
            "n_total": 0,
            "n_hits": 0,
            "n_skipped": 0,
            "n_filtered_length": 0,
            "n_filtered_seq_length": 0,
            "n_filtered_cdr_length": 0,
            "max_id": 0.0,
            "min_id": 0.0,
            "q75": 0.0,
            "best_seq": "",
            "worst_seq": "",
            "time_s": 0.0
        }

# ============================================================================
# MAIN FUNCTION
# ============================================================================

def main():
    # Check if running in interactive mode (no arguments)
    if len(sys.argv) == 1:
        # Interactive mode
        params = get_user_input_interactive()
        args = convert_params_to_args(params)
        
        # Handle analysis if requested
        if params.get('run_analysis'):
            regions_list = args.regions.split(',')
            
            # Run analysis based on user preference
            if params.get('show_analysis', True):
                recommendations = run_interactive_analysis(
                    args.db_root, args.query_seq, regions_list,
                    args.sample_size, args.target_hits
                )
            else:
                # Run silently
                print("\nRunning threshold analysis...")
                import io
                import contextlib
                
                f = io.StringIO()
                with contextlib.redirect_stdout(f):
                    recommendations = run_interactive_analysis(
                        args.db_root, args.query_seq, regions_list,
                        args.sample_size, args.target_hits
                    )
                print("Analysis complete.\n")
            
            if recommendations:
                print("\n" + "="*80)
                print("RECOMMENDED THRESHOLDS")
                print("="*80)
                print(f"\nBased on analysis, these thresholds should give ~{args.target_hits:,} hits:")
                for region, threshold in recommendations.items():
                    print(f"  {region.upper()}: {threshold:.3f} ({threshold*100:.1f}%)")
                
                print("\nOptions:")
                print("  1. Use recommended thresholds")
                print("  2. Adjust thresholds manually")
                print("  3. Proceed without thresholds")
                print("  4. Cancel and exit")
                
                while True:
                    apply = input("\nYour choice (1-4): ").strip()
                    
                    if apply == '1':
                        # Apply recommendations
                        for region, threshold in recommendations.items():
                            setattr(args, f"min_id_{region}", threshold)
                        print("✓ Recommended thresholds applied")
                        break
                    elif apply == '2':
                        # Manual adjustment
                        print("\nAdjust thresholds (Enter to keep recommended value):")
                        for region in regions_list:
                            region = region.strip()
                            recommended = recommendations.get(region)
                            if recommended:
                                prompt = f"  {region.upper()} (recommended: {recommended:.3f}): "
                            else:
                                prompt = f"  {region.upper()}: "
                            
                            value = input(prompt).strip()
                            if value:
                                try:
                                    val = float(value)
                                    if 0.0 <= val <= 1.0:
                                        setattr(args, f"min_id_{region}", val)
                                        print(f"    ✓ Set to {val:.3f}")
                                    else:
                                        print("    ⚠ Invalid value, using recommended")
                                        if recommended:
                                            setattr(args, f"min_id_{region}", recommended)
                                except ValueError:
                                    print("    ⚠ Invalid value, using recommended")
                                    if recommended:
                                        setattr(args, f"min_id_{region}", recommended)
                            elif recommended:
                                setattr(args, f"min_id_{region}", recommended)
                                print(f"    ✓ Using recommended {recommended:.3f}")
                        break
                    elif apply == '3':
                        print("✓ Proceeding without thresholds")
                        break
                    elif apply == '4':
                        print("\nSearch cancelled.")
                        return
                    else:
                        print("Please enter 1, 2, 3, or 4")
    
    else:
        # Command line mode
        parser = argparse.ArgumentParser(description="NPZ FULLSCAN v6 INTERACTIVE")
        
        parser.add_argument("--db-root", help="Path to NPZ database (or use --species)")
        parser.add_argument("--species", choices=['human', 'mouse', 'camel', 'rabbit', 'rat', 'rhesus', 'humanised'],
                          help="Select species instead of providing db-root")
        parser.add_argument("--query-seq", help="Query sequence (direct input)")
        parser.add_argument("--query-file", help="File containing query sequence(s)")
        parser.add_argument("--regions", default="full",
                          help="Regions to search (comma-separated: full,cdr1,cdr2,cdr3,frameworks)")
        parser.add_argument("--outdir", default="./results", help="Output directory")
        parser.add_argument("--tag", help="Optional tag for output files")
        
        # Identity thresholds
        parser.add_argument("--min-id-full", type=float, help="Min identity for full sequence")
        parser.add_argument("--min-id-cdr1", type=float, help="Min identity for CDR1")
        parser.add_argument("--min-id-cdr2", type=float, help="Min identity for CDR2")
        parser.add_argument("--min-id-cdr3", type=float, help="Min identity for CDR3")
        parser.add_argument("--min-id-frameworks", type=float, help="Min identity for frameworks")
        
        # Length filters
        parser.add_argument("--cdr1-len-window", type=int, help="CDR1 length tolerance")
        parser.add_argument("--cdr2-len-window", type=int, help="CDR2 length tolerance")
        parser.add_argument("--cdr3-len-window", type=int, help="CDR3 length tolerance")
        
        # Other options
        parser.add_argument("--scheme", default="imgt", choices=["imgt", "kabat"],
                          help="Numbering scheme")
        parser.add_argument("--disable-seq-length-prefilter", action="store_true",
                          help="Disable sequence length pre-filter")
        parser.add_argument("--disable-cdr-length-prefilter", action="store_true",
                          help="Disable CDR length pre-filter")
        parser.add_argument("--seq-length-threshold", type=float, default=0.25,
                          help="Sequence length difference threshold")
        parser.add_argument("--sample-size", type=int, default=1000,
                          help="Sample size for interactive analysis")
        parser.add_argument("--target-hits", type=int, default=10000,
                          help="Target number of hits for threshold detection")
        
        args = parser.parse_args()
        
        # Set smart defaults for sample size and target hits based on species
        if args.species:
            if args.sample_size == 1000:  # Still using default
                if args.species == 'camel':
                    args.sample_size = 10000
                    print(f"Using larger sample size for Camel: {args.sample_size}")
                elif args.species in ['rabbit', 'rat', 'rhesus']:
                    args.sample_size = 5000
                    print(f"Using adjusted sample size for {args.species}: {args.sample_size}")
                else:
                    args.sample_size = 2000
            
            if args.target_hits == 10000:  # Still using default
                if args.species == 'camel':
                    args.target_hits = 1000
                    print(f"Using realistic target for Camel: {args.target_hits} hits")
                elif args.species in ['rabbit', 'rat', 'rhesus']:
                    args.target_hits = 2000
                    print(f"Using adjusted target for {args.species}: {args.target_hits} hits")
        
        # Handle query file or direct sequence
        if args.query_file:
            # Read query from file
            try:
                with open(args.query_file, 'r') as f:
                    content = f.read()
                    # Handle FASTA format
                    if content.startswith('>'):
                        lines = content.split('\n')
                        args.query_seq = ''.join([line for line in lines if not line.startswith('>')])
                    else:
                        args.query_seq = content
                    # Clean the sequence
                    args.query_seq = re.sub(r"[^A-Za-z]", "", args.query_seq)
                print(f"✓ Loaded query from {args.query_file} ({len(args.query_seq)} aa)")
            except Exception as e:
                print(f"Error reading query file: {e}")
                sys.exit(1)
        elif not args.query_seq:
            print("Error: Either --query-seq or --query-file must be specified")
            parser.print_help()
            sys.exit(1)
        
        # Handle species shortcut
        if args.species and not args.db_root:
            # Try to find the database path
            DB_BASE_PATHS = [
                os.path.expanduser("~/KA-Search/extracted/oas-paper/oasdb_20230111/Heavy"),
                "/home/sasenefrem/KA-Search/extracted/oas-paper/oasdb_20230111/Heavy",
                "./Heavy",
            ]
            
            db_base = None
            for path in DB_BASE_PATHS:
                if os.path.exists(path):
                    db_base = path
                    break
            
            if db_base:
                args.db_root = os.path.join(db_base, args.species.capitalize())
                print(f"Using database: {args.db_root}")
            else:
                print("Error: Could not find database base path. Please use --db-root instead.")
                sys.exit(1)
        elif not args.db_root:
            print("Error: Either --db-root or --species must be specified")
            parser.print_help()
            sys.exit(1)
        
        # Ask about interactive analysis
        print("\n" + "="*80)
        print("NPZ FULLSCAN v6 INTERACTIVE")
        print("="*80 + "\n")
        
        print("Do you want to:")
        print("  1. Run threshold analysis first (recommended)")
        print("  2. Skip directly to full scan")
        print("  3. Run analysis only (no full scan)")
        
        while True:
            choice = input("\nYour choice (1/2/3): ").strip()
            if choice in ['1', '2', '3']:
                break
            print("Please enter 1, 2, or 3")
        
        if choice == '1' or choice == '3':
            regions_list = args.regions.split(',')
            
            # Determine if we show analysis or run silently
            show_analysis = (choice == '1')  # Show for option 1, silent for option 3 in cmd mode
            
            if show_analysis:
                recommendations = run_interactive_analysis(
                    args.db_root, args.query_seq, regions_list,
                    args.sample_size, args.target_hits
                )
            else:
                # Run silently for command line mode
                print("\nRunning threshold analysis...")
                import io
                import contextlib
                
                f = io.StringIO()
                with contextlib.redirect_stdout(f):
                    recommendations = run_interactive_analysis(
                        args.db_root, args.query_seq, regions_list,
                        args.sample_size, args.target_hits
                    )
                print("Analysis complete.\n")
            
            if recommendations:
                print("\n" + "="*80)
                print("RECOMMENDED THRESHOLDS")
                print("="*80)
                print(f"\nBased on analysis, these thresholds should give ~{args.target_hits:,} hits:")
                for region, threshold in recommendations.items():
                    print(f"  {region.upper()}: {threshold:.3f} ({threshold*100:.1f}%)")
                
                print("\nOptions:")
                print("  1. Use recommended thresholds")
                print("  2. Adjust thresholds manually")
                print("  3. Cancel and exit")
                
                while True:
                    apply = input("\nYour choice (1-3): ").strip()
                    
                    if apply == '1':
                        # Apply recommendations
                        for region, threshold in recommendations.items():
                            setattr(args, f"min_id_{region}", threshold)
                        print("✓ Recommended thresholds applied")
                        break
                    elif apply == '2':
                        # Manual adjustment
                        print("\nAdjust thresholds (Enter to keep recommended value):")
                        for region in regions_list:
                            region = region.strip()
                            recommended = recommendations.get(region)
                            if recommended:
                                prompt = f"  {region.upper()} (recommended: {recommended:.3f}): "
                            else:
                                prompt = f"  {region.upper()}: "
                            
                            value = input(prompt).strip()
                            if value:
                                try:
                                    val = float(value)
                                    if 0.0 <= val <= 1.0:
                                        setattr(args, f"min_id_{region}", val)
                                        print(f"    ✓ Set to {val:.3f}")
                                    else:
                                        print("    ⚠ Invalid value, using recommended")
                                        if recommended:
                                            setattr(args, f"min_id_{region}", recommended)
                                except ValueError:
                                    print("    ⚠ Invalid value, using recommended")
                                    if recommended:
                                        setattr(args, f"min_id_{region}", recommended)
                            elif recommended:
                                setattr(args, f"min_id_{region}", recommended)
                                print(f"    ✓ Using recommended {recommended:.3f}")
                        break
                    elif apply == '3':
                        print("\nSearch cancelled.")
                        return
                    else:
                        print("Please enter 1, 2, or 3")
            
            if choice == '3':
                print("\nAnalysis complete. Exiting without full scan.")
                return
    
    # Run full scan
    start_time = dt.datetime.now()
    
    # Extract query CDRs
    q_full = re.sub(r"[^A-Za-z]", "", args.query_seq.upper())
    query_cdrs = {"full": q_full}
    
    regions_list = [r.strip() for r in args.regions.split(",")]
    
    # Handle frameworks region
    if "frameworks" in regions_list:
        query_cdrs["frameworks"] = q_full  # Use full sequence for frameworks
    
    # Only try to extract CDRs if we're searching for them
    if set(regions_list) != {"full"} and not (len(regions_list) == 1 and regions_list[0] == "frameworks"):
        qc = extract_cdrs_from_query_sequence(q_full, args.scheme)
        query_cdrs.update(qc)
        
        # Show extracted CDRs
        print("\n" + "="*80)
        print("QUERY CDR EXTRACTION")
        print("="*80)
        print(f"Query sequence length: {len(q_full)} aa")
        print(f"\nExtracted regions from query:")
        
        if query_cdrs.get('cdr1'):
            print(f"  CDR1: {query_cdrs['cdr1']} (len={len(query_cdrs['cdr1'])})")
        else:
            print(f"  CDR1: Could not extract (using heuristic search)")
            
        if query_cdrs.get('cdr2'):
            print(f"  CDR2: {query_cdrs['cdr2']} (len={len(query_cdrs['cdr2'])})")
        else:
            print(f"  CDR2: Could not extract (using heuristic search)")
            
        if query_cdrs.get('cdr3'):
            print(f"  CDR3: {query_cdrs['cdr3']} (len={len(query_cdrs['cdr3'])})")
        else:
            # Try heuristic CDR3 extraction as fallback
            h3 = heuristic_cdrh3(q_full)
            if h3:
                query_cdrs['cdr3'] = h3
                print(f"  CDR3: {h3} (len={len(h3)}) [heuristic]")
            else:
                print(f"  CDR3: Could not extract")
        print("="*80 + "\n")
    else:
        # Full sequence search - CDR extraction not needed
        print("\n" + "="*80)
        print("QUERY SEQUENCE")
        print("="*80)
        print(f"Full sequence length: {len(q_full)} aa")
        print(f"Searching: FULL SEQUENCE")
        print("CDR extraction: Not needed for full sequence search")
        print("="*80 + "\n")
    
    # Build config
    cfg = {
        "scheme": args.scheme,
        "regions": args.regions,
        "min_id_full": args.min_id_full,
        "min_id_cdr1": args.min_id_cdr1,
        "min_id_cdr2": args.min_id_cdr2,
        "min_id_cdr3": args.min_id_cdr3,
        "min_id_frameworks": args.min_id_frameworks,
        "use_length_filter": False,
        "length_windows": {},
        "query_lengths": {
            "cdr1": len(query_cdrs.get("cdr1", "")),
            "cdr2": len(query_cdrs.get("cdr2", "")),
            "cdr3": len(query_cdrs.get("cdr3", ""))
        },
        "enable_seq_length_prefilter": not args.disable_seq_length_prefilter,
        "enable_cdr_length_prefilter": not args.disable_cdr_length_prefilter,
        "seq_length_threshold": args.seq_length_threshold,
    }
    
    # Add length windows
    if args.cdr1_len_window is not None:
        cfg["use_length_filter"] = True
        cfg["length_windows"]["cdr1"] = args.cdr1_len_window
    if args.cdr2_len_window is not None:
        cfg["use_length_filter"] = True
        cfg["length_windows"]["cdr2"] = args.cdr2_len_window
    if args.cdr3_len_window is not None:
        cfg["use_length_filter"] = True
        cfg["length_windows"]["cdr3"] = args.cdr3_len_window
    
    # Find shards
    npz_files = sorted(glob.glob(os.path.join(args.db_root, "*.npz")))
    
    if not npz_files:
        print("ERROR: No NPZ files found!")
        return
    
    print(f"Found {len(npz_files)} shards to process\n")
    
    # Create output structure
    species = get_species_from_path(args.db_root)
    thresholds = {
        f'min_id_{r}': cfg.get(f'min_id_{r}')
        for r in regions_list
    }
    
    csv_path, output_folder, csv_folder = create_output_structure(
        args.outdir, species, args.regions, thresholds, args.tag
    )
    
    print(f"Output folder: {output_folder}")
    print(f"CSV subfolder: {csv_folder}")
    print(f"CSV files: {os.path.basename(csv_path)}")
    print("📝 Note: CSVs update in real-time as hits are found!")
    print("📊 Excel will be created in main folder at the end\n")
    
    # Create results file path
    results_path = csv_path.replace('.csv', '_results.csv')
    
    # Write CSV headers
    with open(csv_path, "w", encoding="utf-8") as fh:
        write_metadata_block(fh, args, query_cdrs, start_time, cfg)
        fh.write("shard_id,total_seqs,hits,skipped,length_filtered,seq_length_filtered,"
                 "cdr_length_filtered,max_id,min_id,q75_id,best_seq_preview,worst_seq_preview,time_s\n")
    
    # Write results file header
    with open(results_path, "w", encoding="utf-8") as fh:
        fh.write("# Detailed hit results\n")
        fh.write(f"# Query CDR1: {query_cdrs.get('cdr1', 'N/A')}\n")
        fh.write(f"# Query CDR2: {query_cdrs.get('cdr2', 'N/A')}\n")
        fh.write(f"# Query CDR3: {query_cdrs.get('cdr3', 'N/A')}\n")
        fh.write(f"# Searched regions: {args.regions}\n")
        fh.write("#\n")
        
        # Always include ALL CDRs in header, regardless of search
        header_parts = ["shard", "index", "avg_identity", "full_identity",
                       "cdr1_identity", "cdr1_seq",
                       "cdr2_identity", "cdr2_seq", 
                       "cdr3_identity", "cdr3_seq"]
        
        # Add frameworks if it was searched
        if "frameworks" in regions_list:
            header_parts.extend(["frameworks_identity", "frameworks_seq"])
        
        header_parts.append("full_sequence")
        fh.write(",".join(header_parts) + "\n")
    
    # Run full scan
    print("="*80)
    print("RUNNING FULL SCAN")
    print("="*80 + "\n")
    
    total_hits = 0
    total_filtered_length = 0
    total_filtered_seq_length = 0
    total_filtered_cdr_length = 0
    global_best_id = 0.0
    global_worst_id = 1.0
    global_best_seq = ""
    global_worst_seq = ""
    
    outer_bar = tqdm(total=len(npz_files), desc="Processing shards", ncols=100, position=0, leave=True)
    
    for i, npz_path in enumerate(npz_files, start=1):
        stats = process_npz_file_fullscan(npz_path, cfg, query_cdrs, i, len(npz_files))
        
        # Update totals
        total_hits += stats["n_hits"]
        total_filtered_length += stats["n_filtered_length"]
        total_filtered_seq_length += stats["n_filtered_seq_length"]
        total_filtered_cdr_length += stats["n_filtered_cdr_length"]
        
        # Track global best/worst
        if stats["max_id"] > global_best_id:
            global_best_id = stats["max_id"]
            global_best_seq = stats["best_seq"]
        if stats["n_hits"] > 0 and stats["min_id"] < global_worst_id:
            global_worst_id = stats["min_id"]
            global_worst_seq = stats["worst_seq"]
        
        # Append summary to CSV
        with open(csv_path, "a", encoding="utf-8") as fh:
            fh.write(
                f"{stats['shard']},{stats['n_total']},{stats['n_hits']},{stats['n_skipped']},"
                f"{stats['n_filtered_length']},{stats['n_filtered_seq_length']},"
                f"{stats['n_filtered_cdr_length']},{stats['max_id']:.6f},{stats['min_id']:.6f},"
                f"{stats['q75']:.6f},\"{stats['best_seq']}\",\"{stats['worst_seq']}\",{stats['time_s']:.3f}\n"
            )
            fh.flush()  # Real-time update!
        
        # Append hit details to results file
        if stats["hits"]:
            with open(results_path, "a", encoding="utf-8") as fh:
                for hit in stats["hits"]:
                    row_parts = [
                        hit["shard"],
                        str(hit["index"]),
                        f"{hit['avg_identity']:.6f}",
                        f"{hit.get('full_identity', 0.0):.6f}"
                    ]
                    
                    # Always add all CDRs (they're always present now)
                    for cdr in ['cdr1', 'cdr2', 'cdr3']:
                        identity_key = f"{cdr}_identity"
                        seq_key = f"{cdr}_seq"
                        
                        row_parts.append(f"{hit.get(identity_key, 0.0):.6f}")
                        
                        seq = hit.get(seq_key, "")
                        if seq:
                            seq = seq.replace('"', '""')
                            row_parts.append(f'"{seq}"')
                        else:
                            row_parts.append('""')
                    
                    # Add frameworks if present
                    if "frameworks_identity" in hit:
                        row_parts.append(f"{hit['frameworks_identity']:.6f}")
                        seq = hit.get("frameworks_seq", "").replace('"', '""')
                        row_parts.append(f'"{seq}"')
                    
                    # Add full sequence
                    full_seq = hit.get("full_sequence", "").replace('"', '""')
                    row_parts.append(f'"{full_seq}"')
                    
                    fh.write(",".join(row_parts) + "\n")
                    fh.flush()  # Force write to disk immediately
        
        # Also flush summary after each shard
        with open(csv_path, "a", encoding="utf-8") as fh:
            fh.flush()
        
        outer_bar.update(1)
    
    outer_bar.close()
    
    # Calculate elapsed time
    end_time = dt.datetime.now()
    elapsed = (end_time - start_time).total_seconds()
    
    # Write summary footer
    with open(csv_path, "a", encoding="utf-8") as fh:
        fh.write("\n# ============================================\n")
        fh.write("# RUN SUMMARY\n")
        fh.write("# ============================================\n")
        fh.write(f"# Run completed: {end_time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        fh.write(f"# Elapsed time: {elapsed:.2f} seconds ({elapsed/60:.2f} minutes)\n")
        fh.write(f"# Total shards processed: {len(npz_files)}\n")
        fh.write(f"# Total hits: {total_hits:,}\n")
        if cfg["use_length_filter"]:
            fh.write(f"# Filtered by CDR length: {total_filtered_length:,}\n")
        if total_filtered_seq_length > 0:
            fh.write(f"# Filtered by sequence length: {total_filtered_seq_length:,}\n")
        if total_filtered_cdr_length > 0:
            fh.write(f"# Filtered by CDR length pre-filter: {total_filtered_cdr_length:,}\n")
        fh.write(f"# Best identity found: {global_best_id:.6f}\n")
        fh.write(f"# Worst identity found: {global_worst_id:.6f}\n")
        fh.write("# ============================================\n")
    
    # Print summary
    print(f"\n✅ Done – all {len(npz_files)} shards processed.")
    print(f"\n{'='*80}")
    print("RUN SUMMARY")
    print(f"{'='*80}")
    print(f"Completed: {end_time.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Elapsed time: {elapsed:.2f} seconds ({elapsed/60:.2f} minutes)")
    print(f"Total hits: {total_hits:,}")
    
    # Export to Excel (in the main dated folder, not CSV subfolder)
    excel_path = None
    if total_hits > 0:
        print("\n📊 Creating Excel file...")
        excel_path = export_to_excel(csv_path, results_path, output_folder)
    
    # Final output summary
    print(f"\n{'='*80}")
    print("FILES SAVED")
    print(f"{'='*80}")
    
    print(f"\n📁 Output folder: {output_folder}")
    
    if excel_path and os.path.exists(excel_path):
        print(f"\n📊 Excel file (main folder):")
        print(f"   {os.path.basename(excel_path)}")
        print(f"   • Summary sheet")
        print(f"   • Results sheet ({total_hits:,} rows)")
        print(f"   • Info sheet")
        print(f"   • Contains everything - ready to open!")
    
    print(f"\n📁 CSV files (csvs/ subfolder):")
    print(f"   1. {os.path.basename(csv_path)}")
    print(f"      • Statistics per shard")
    print(f"      • Processing times")
    print(f"   2. {os.path.basename(results_path)}")
    print(f"      • {total_hits:,} hits with ALL sequences")
    print(f"      • Real-time updates during processing")
    
    print(f"\n✅ Folder structure:")
    print(f"   {output_folder}/")
    print(f"   ├── {os.path.basename(excel_path) if excel_path else 'npz_scan.xlsx'} (Excel with everything)")
    print(f"   └── csvs/")
    print(f"       ├── {os.path.basename(csv_path)}")
    print(f"       └── {os.path.basename(results_path)}")
    
    print(f"\n📊 All {total_hits:,} sequences saved in real-time as found!")
    print("="*80 + "\n")


if __name__ == "__main__":
    main()