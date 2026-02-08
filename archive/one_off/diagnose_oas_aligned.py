#!/usr/bin/env python3
"""
Advanced diagnostic for OAS-aligned NPZ files.
Based on the paper's description of 200 unique positions in canonical alignment.
"""

import numpy as np
import glob
import os
import sys

def inspect_oas_aligned_structure(npz_path):
    """Inspect OAS-aligned NPZ structure with understanding of 200-position alignment."""
    
    print("="*80)
    print(f"OAS-ALIGNED NPZ DIAGNOSTIC")
    print(f"File: {os.path.basename(npz_path)}")
    print("="*80)
    
    # Load the file
    data = np.load(npz_path, allow_pickle=True)
    
    # Check structure
    print("\n1. NPZ STRUCTURE:")
    print("-" * 40)
    for key in data.files:
        arr = data[key]
        print(f"  '{key}': shape={arr.shape}, dtype={arr.dtype}")
    
    # Get the main array
    if 'numberings' in data:
        arr = data['numberings']
    else:
        arr = data.get('arr_0', data[data.files[0]])
    
    print(f"\n2. ALIGNMENT STRUCTURE (200 positions as per paper):")
    print("-" * 40)
    print(f"  Total sequences: {len(arr)}")
    if len(arr) > 0:
        print(f"  Array shape per sequence: {arr[0].shape if hasattr(arr[0], 'shape') else len(arr[0])}")
        print(f"  Expected: 200 positions for canonical alignment")
        actual_positions = len(arr[0])
        print(f"  Actual positions: {actual_positions}")
        if actual_positions == 200:
            print("  ✓ Matches expected 200-position alignment")
        else:
            print(f"  ⚠ Different from expected ({actual_positions} vs 200)")
    
    # Analyze sequences with V-domain understanding
    print(f"\n3. V-DOMAIN ANALYSIS:")
    print("-" * 40)
    print("  Note: OAS-aligned contains ONLY V-domain (no constant region)")
    print("  Expected V-domain length: 110-130 aa")
    
    # Sample sequences and check structure
    sample_size = min(100, len(arr))
    v_domain_lengths = []
    has_cdr3_end = 0
    truncated_looking = 0
    
    for i in range(sample_size):
        row = arr[i]
        
        # Try different extraction methods
        # Method 1: Extract positions 0-128 (standard V-domain)
        v_domain = extract_sequence(row[0:128])
        v_domain_lengths.append(len(v_domain))
        
        # Check for CDR3/FR4 boundary
        if 'WGQG' in v_domain or 'WGKG' in v_domain or 'WGRG' in v_domain:
            has_cdr3_end += 1
        
        # Check if looks truncated (ends without FR4)
        if len(v_domain) < 100 and not any(motif in v_domain for motif in ['WGQG', 'WGKG', 'WGRG']):
            truncated_looking += 1
        
        if i < 3:  # Show first few sequences
            print(f"\n  Sequence {i+1}:")
            print(f"    V-domain (0-128): {len(v_domain)} aa")
            print(f"    Start: {v_domain[:30]}...")
            print(f"    End: ...{v_domain[-30:]}")
            
            # Check CDR3 region specifically
            cdr3_region = extract_sequence(row[104:128])
            print(f"    CDR3 region (104-128): {cdr3_region}")
    
    # Statistics
    print(f"\n4. V-DOMAIN STATISTICS (from {sample_size} sequences):")
    print("-" * 40)
    v_domain_lengths.sort()
    print(f"  Length range: {min(v_domain_lengths)}-{max(v_domain_lengths)} aa")
    print(f"  Mean length: {sum(v_domain_lengths)/len(v_domain_lengths):.1f} aa")
    print(f"  Median length: {v_domain_lengths[len(v_domain_lengths)//2]} aa")
    print(f"  With FR4 motif: {has_cdr3_end}/{sample_size} ({has_cdr3_end*100/sample_size:.1f}%)")
    print(f"  Truncated-looking: {truncated_looking}/{sample_size} ({truncated_looking*100/sample_size:.1f}%)")
    
    # Distribution
    print(f"\n  Length distribution:")
    ranges = [(60,70), (70,80), (80,90), (90,100), (100,110), (110,120), (120,130), (130,140)]
    for start, end in ranges:
        count = sum(1 for l in v_domain_lengths if start <= l < end)
        pct = count * 100 / len(v_domain_lengths)
        bar = "█" * int(pct / 2)
        print(f"    {start:3d}-{end:3d} aa: {bar} {count:3d} ({pct:5.1f}%)")
    
    # Check what's at position 128-200
    print(f"\n5. POSITIONS 128-200 (beyond typical V-domain):")
    print("-" * 40)
    beyond_v = extract_sequence(arr[0][128:200])
    print(f"  Content after position 128: {len(beyond_v)} aa")
    if beyond_v:
        print(f"  Sequence: {beyond_v}")
    else:
        print(f"  Empty or gaps")
    
    # Special check for camel long CDR3s
    print(f"\n6. CAMEL-SPECIFIC ISSUES:")
    print("-" * 40)
    print("  Paper notes: Camel sequences with CDR3 >37 residues have issues")
    
    long_cdr3_count = 0
    for i in range(sample_size):
        row = arr[i]
        cdr3 = extract_cdr3(row)
        if len(cdr3) > 37:
            long_cdr3_count += 1
            if long_cdr3_count == 1:  # Show first example
                print(f"  Example long CDR3: {cdr3} ({len(cdr3)} aa)")
    
    print(f"  Sequences with CDR3 >37 aa: {long_cdr3_count}/{sample_size}")
    
    data.close()
    
    # Conclusion
    print(f"\n7. INTERPRETATION:")
    print("="*80)
    
    if sum(l < 100 for l in v_domain_lengths) > len(v_domain_lengths) * 0.5:
        print("⚠ ISSUE DETECTED: Many sequences appear incomplete")
        print("  - Over 50% of sequences are <100 aa (expected 110-130 aa for V-domain)")
        print("  - This suggests either:")
        print("    1. Sequences are genuinely truncated in the database")
        print("    2. Alignment issues specific to camel sequences")
        print("    3. Different numbering scheme causing extraction problems")
        print("\n  RECOMMENDATION: Contact OAS maintainers about camel data quality")
    else:
        print("✓ Sequences appear complete")
        print("  - Most sequences are in expected V-domain range (110-130 aa)")
    
    return v_domain_lengths

def extract_sequence(int_row):
    """Extract amino acid sequence from integer array."""
    out = []
    for v in int_row:
        if v is None:
            continue
        
        try:
            iv = int(v)
            # Skip gaps and padding (21, 124, 0 are common gap/pad codes)
            if iv in [0, 21, 124]:
                continue
            # ASCII letters
            elif 65 <= iv <= 90:  # A-Z
                out.append(chr(iv))
            elif 97 <= iv <= 122:  # a-z
                out.append(chr(iv).upper())
            # Direct amino acid codes (1-20)
            elif 1 <= iv <= 20:
                aa_map = {
                    1: 'A', 2: 'C', 3: 'D', 4: 'E', 5: 'F',
                    6: 'G', 7: 'H', 8: 'I', 9: 'K', 10: 'L',
                    11: 'M', 12: 'N', 13: 'P', 14: 'Q', 15: 'R',
                    16: 'S', 17: 'T', 18: 'V', 19: 'W', 20: 'Y'
                }
                if iv in aa_map:
                    out.append(aa_map[iv])
        except:
            continue
    
    return "".join(out).replace("-", "")

def extract_cdr3(int_row):
    """Extract CDR3 from aligned sequence."""
    # CDR3 typically at positions 105-117 in IMGT
    cdr3_region = int_row[104:128]
    cdr3_aa = extract_sequence(cdr3_region)
    
    if not cdr3_aa:
        return ""
    
    # Skip conserved C
    if cdr3_aa[0] == 'C':
        cdr3_aa = cdr3_aa[1:]
    
    # Find FR4 start (W)
    w_pos = cdr3_aa.find('W')
    if w_pos > 0:
        return cdr3_aa[:w_pos]
    
    return cdr3_aa[:13]  # Default max CDR3 length

def main():
    if len(sys.argv) > 1:
        npz_path = sys.argv[1]
        if os.path.exists(npz_path):
            inspect_oas_aligned_structure(npz_path)
            return
    
    # Try to find NPZ files
    patterns = [
        "*.npz",
        "Camel/*.npz",
        "../Camel/*.npz",
        "*/Camel/*.npz"
    ]
    
    for pattern in patterns:
        files = glob.glob(pattern)
        if files:
            print(f"Found {len(files)} NPZ files")
            
            # Analyze all files
            all_lengths = []
            for npz_file in files:
                lengths = inspect_oas_aligned_structure(npz_file)
                all_lengths.extend(lengths)
                
                if len(files) > 1 and npz_file != files[-1]:
                    print("\n" + "="*80 + "\n")
            
            # Overall summary
            if len(files) > 1:
                print("\n" + "="*80)
                print("OVERALL SUMMARY:")
                print("="*80)
                print(f"Total sequences analyzed: {len(all_lengths)}")
                print(f"Mean V-domain length: {sum(all_lengths)/len(all_lengths):.1f} aa")
                truncated = sum(1 for l in all_lengths if l < 100)
                print(f"Likely truncated (<100 aa): {truncated}/{len(all_lengths)} ({truncated*100/len(all_lengths):.1f}%)")
            return
    
    print("ERROR: No NPZ files found")
    print("Usage: python3 diagnose_oas.py [npz_file]")
    print("Or run from directory containing NPZ files")

if __name__ == "__main__":
    main()
