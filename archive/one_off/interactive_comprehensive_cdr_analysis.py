#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Interactive Comprehensive CDR Analysis with Threshold Recommendations

This script:
1. Asks you questions about your search preferences
2. Runs a comprehensive analysis sampling sequences from the database
3. Calculates optimal thresholds to achieve your target hit count
4. Provides ready-to-run KA-Search commands
"""

import os
import re
import glob
import signal
import argparse
import numpy as np
import pandas as pd
from typing import Dict, List, Tuple


# Constants
AA_PAD = 124  
AA_GAP = 0

class _Timeout(Exception):
    pass

def _alarm_handler(signum, frame):
    raise _Timeout()


# ============================================================================
# SEQUENCE EXTRACTION FUNCTIONS
# ============================================================================

def _ints_to_aa(int_row):
    """Ultra-robust decoder."""
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


def extract_cdr1_from_npz(int_row):
    return _ints_to_aa(int_row[39:46])


def extract_cdr2_from_npz(int_row):
    part1 = _ints_to_aa(int_row[66:75])
    part2 = _ints_to_aa(int_row[85:91])
    return part1 + part2


def extract_cdr3_from_npz(int_row):
    cdr3_region = int_row[150:190]
    cdr3_aa = _ints_to_aa(cdr3_region)
    if cdr3_aa and cdr3_aa[0] == 'C':
        w_pos = cdr3_aa.find('W')
        if w_pos > 0:
            return cdr3_aa[:w_pos]
    return cdr3_aa


def extract_whole_from_npz(int_row):
    return _ints_to_aa(int_row[16:190])


def heuristic_cdrh3(full_seq: str) -> str:
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


def fast_edit_distance(s1: str, s2: str) -> int:
    len1, len2 = len(s1), len(s2)
    if len1 < len2:
        s1, s2 = s2, s1
        len1, len2 = len2, len1
    if len2 == 0:
        return len1
    prev_row = list(range(len2 + 1))
    for i, c1 in enumerate(s1):
        curr_row = [i + 1]
        for j, c2 in enumerate(s2):
            insertions = prev_row[j + 1] + 1
            deletions = curr_row[j] + 1
            substitutions = prev_row[j] + (c1 != c2)
            curr_row.append(min(insertions, deletions, substitutions))
        prev_row = curr_row
    return prev_row[-1]


def alignment_identity(a: str, b: str) -> float:
    if not a or not b:
        return 0.0
    max_len = max(len(a), len(b))
    if max_len == 0:
        return 0.0
    dist = fast_edit_distance(a, b)
    identity = 1.0 - (dist / max_len)
    return max(0.0, identity)


def extract_query_cdrs(seq: str, scheme: str = "imgt") -> Dict[str, str]:
    seq_clean = seq.upper().replace("-", "").replace(".", "")
    try:
        from anarci import run_anarci
    except ImportError:
        return {"cdr1": "", "cdr2": "", "cdr3": heuristic_cdrh3(seq_clean)}
    
    signal.signal(signal.SIGALRM, _alarm_handler)
    signal.alarm(5)
    
    try:
        res = run_anarci([("H", seq_clean)], scheme=scheme, allowed_species=None)
        signal.alarm(0)
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
            "cdr1": grab(27.0, 38.9),
            "cdr2": grab(56.0, 65.9),
            "cdr3": grab(105.0, 117.9)
        }
    except Exception:
        signal.alarm(0)
        return {"cdr1": "", "cdr2": "", "cdr3": heuristic_cdrh3(seq_clean)}


# ============================================================================
# SHARD ANALYSIS
# ============================================================================

def analyze_shard_comprehensive(npz_path: str, query_data: Dict, sample_size: int, 
                                 shard_num: int, total_shards: int, regions_to_calc: List[str]):
    """Comprehensive per-shard analysis."""
    
    shard_name = os.path.basename(npz_path)
    
    try:
        # Load the shard
        arr = np.load(npz_path, allow_pickle=True)
        
        # Check for numberings key (this is what KA-Search uses)
        if "numberings" not in arr:
            arr.close()
            return None, None
        
        numberings = arr["numberings"]
        total_seqs = len(numberings)
        
        if total_seqs == 0:
            arr.close()
            return None, None
        
        # Random sample
        sample_indices = np.random.choice(total_seqs, min(sample_size, total_seqs), replace=False)
        
        results = []
        
        # Get query data
        q_cdr1 = query_data.get("cdr1", "")
        q_cdr2 = query_data.get("cdr2", "")
        q_cdr3 = query_data.get("cdr3", "")
        q_full = query_data.get("full", "")
        
        for idx in sample_indices:
            row = numberings[idx]
            
            # Extract CDRs and full sequence
            db_cdr1 = extract_cdr1_from_npz(row)
            db_cdr2 = extract_cdr2_from_npz(row)
            db_cdr3 = extract_cdr3_from_npz(row)
            
            # Fallback for CDR3
            if not db_cdr3 or not db_cdr3.startswith('C'):
                full = extract_whole_from_npz(row)
                db_cdr3 = heuristic_cdrh3(full)
            
            db_full = extract_whole_from_npz(row)
            
            # Calculate identities based on what user requested
            result = {
                "shard": shard_name,
                "index": int(idx),
                "cdr1_seq": db_cdr1,
                "cdr2_seq": db_cdr2,
                "cdr3_seq": db_cdr3,
            }
            
            if 'cdr1' in regions_to_calc:
                result["id_cdr1"] = alignment_identity(q_cdr1, db_cdr1) if q_cdr1 and db_cdr1 else 0.0
            if 'cdr2' in regions_to_calc:
                result["id_cdr2"] = alignment_identity(q_cdr2, db_cdr2) if q_cdr2 and db_cdr2 else 0.0
            if 'cdr3' in regions_to_calc:
                result["id_cdr3"] = alignment_identity(q_cdr3, db_cdr3) if q_cdr3 and db_cdr3 else 0.0
            if 'full' in regions_to_calc:
                result["id_full"] = alignment_identity(q_full, db_full) if q_full and db_full else 0.0
            
            # Calculate combined scores if using multiple CDRs
            cdr_ids = []
            if 'cdr1' in regions_to_calc:
                cdr_ids.append(result["id_cdr1"])
            if 'cdr2' in regions_to_calc:
                cdr_ids.append(result["id_cdr2"])
            if 'cdr3' in regions_to_calc:
                cdr_ids.append(result["id_cdr3"])
            
            if len(cdr_ids) > 0:
                result["cdr_avg"] = np.mean(cdr_ids)
                result["cdr_min"] = np.min(cdr_ids)
            
            results.append(result)
        
        arr.close()
        
        df = pd.DataFrame(results)
        
        # Per-shard statistics
        stats = {
            "shard_num": shard_num,
            "shard_name": shard_name,
            "total_seqs": total_seqs,
            "sampled": len(df),
        }
        
        # Add max/mean for each region
        for region in regions_to_calc:
            col = f"id_{region}"
            if col in df.columns:
                stats[f"max_{region}"] = df[col].max()
                stats[f"mean_{region}"] = df[col].mean()
        
        if len(cdr_ids) > 0:
            stats["max_cdr_avg"] = df["cdr_avg"].max()
            stats["max_cdr_min"] = df["cdr_min"].max()
        
        # Count triple matches at different thresholds (if using CDRs)
        if set(['cdr1', 'cdr2', 'cdr3']).issubset(regions_to_calc):
            stats["triple_40"] = len(df[(df["id_cdr1"] >= 0.40) & (df["id_cdr2"] >= 0.40) & (df["id_cdr3"] >= 0.40)])
            stats["triple_35"] = len(df[(df["id_cdr1"] >= 0.35) & (df["id_cdr2"] >= 0.35) & (df["id_cdr3"] >= 0.35)])
            stats["triple_30"] = len(df[(df["id_cdr1"] >= 0.30) & (df["id_cdr2"] >= 0.30) & (df["id_cdr3"] >= 0.30)])
        
        # Progress indicator
        pct = (shard_num / total_shards) * 100
        progress_str = f"[{shard_num:3d}/{total_shards}] {pct:5.1f}% | {shard_name[:40]:40s} | "
        
        if 'cdr1' in regions_to_calc and 'cdr2' in regions_to_calc and 'cdr3' in regions_to_calc:
            progress_str += f"Max: C1={stats['max_cdr1']:.2f} C2={stats['max_cdr2']:.2f} C3={stats['max_cdr3']:.2f}"
        elif 'full' in regions_to_calc:
            progress_str += f"Max Full={stats['max_full']:.2f}"
        
        print(progress_str)
        
        return df, stats
        
    except Exception as e:
        # Silently skip failed shards like the original script
        return None, None


# ============================================================================
# THRESHOLD CALCULATION
# ============================================================================

def calculate_optimal_thresholds(combined_df: pd.DataFrame, stats_df: pd.DataFrame,
                                  target_hits: int, regions: List[str]):
    """Calculate optimal thresholds for target hit count."""
    
    # Calculate scale factor
    total_sampled = stats_df['sampled'].sum()
    num_shards = len(stats_df)
    full_db_size = num_shards * 5_000_000
    scale_factor = full_db_size / total_sampled
    
    needed_in_sample = target_hits / scale_factor
    
    print(f"\n{'='*100}")
    print("THRESHOLD RECOMMENDATIONS FOR YOUR TARGET")
    print(f"{'='*100}")
    print(f"Target hits in full DB:  {target_hits:,}")
    print(f"Scale factor:            {scale_factor:.1f}x")
    print(f"Needed in sample:        ~{needed_in_sample:.1f} sequences\n")
    
    recommendations = []
    
    # If using all 3 CDRs, test both equal and differential thresholds
    if set(['cdr1', 'cdr2', 'cdr3']).issubset(regions):
        print("Testing threshold configurations (ALL 3 CDRs must meet criteria):\n")
        print(f"{'Configuration':<25} {'Sample':<10} {'Full DB Est.':<15} {'Status'}")
        print("-" * 70)
        
        # Test equal thresholds from 50% down to 10%
        equal_configs = []
        for thresh in [0.50, 0.45, 0.40, 0.38, 0.37, 0.35, 0.33, 0.30, 0.28, 0.25, 0.23, 0.20, 0.18, 0.15, 0.12, 0.10]:
            count = len(combined_df[
                (combined_df['id_cdr1'] >= thresh) & 
                (combined_df['id_cdr2'] >= thresh) & 
                (combined_df['id_cdr3'] >= thresh)
            ])
            expected = int(count * scale_factor)
            
            # Determine status
            if expected == 0:
                status = "❌ Zero hits"
            elif expected < target_hits * 0.5:
                status = "⚠️  Few hits"
            elif target_hits * 0.5 <= expected <= target_hits * 5:
                status = "✅ GOOD RANGE"
                equal_configs.append((thresh, thresh, thresh, expected))
            elif expected < target_hits * 10:
                status = "⚠️  Many hits"
            else:
                status = "❌ Too many"
            
            print(f"All CDRs >= {thresh:.0%}    {count:<10} ~{expected:<14,} {status}")
        
        # Calculate ONE BEST differential config based on actual data
        print()
        print("="*70)
        print("CALCULATING OPTIMAL DIFFERENTIAL THRESHOLDS")
        print("="*70)
        
        # Analyze the distribution for each CDR
        max_cdr1 = combined_df['id_cdr1'].max()
        max_cdr2 = combined_df['id_cdr2'].max()
        max_cdr3 = combined_df['id_cdr3'].max()
        
        q75_cdr1 = combined_df['id_cdr1'].quantile(0.75)
        q75_cdr2 = combined_df['id_cdr2'].quantile(0.75)
        q75_cdr3 = combined_df['id_cdr3'].quantile(0.75)
        
        print(f"\nCDR Identity Analysis:")
        print(f"  CDR1: Max={max_cdr1:.2%}, Q75={q75_cdr1:.2%}")
        print(f"  CDR2: Max={max_cdr2:.2%}, Q75={q75_cdr2:.2%}")
        print(f"  CDR3: Max={max_cdr3:.2%}, Q75={q75_cdr3:.2%}")
        
        # Find a base threshold that gets us in the ballpark
        # Start by finding what equal threshold would give us close to target
        base_thresh = None
        for thresh in np.arange(0.10, 0.60, 0.01):
            count = len(combined_df[
                (combined_df['id_cdr1'] >= thresh) & 
                (combined_df['id_cdr2'] >= thresh) & 
                (combined_df['id_cdr3'] >= thresh)
            ])
            expected = int(count * scale_factor)
            if expected > 0 and expected <= target_hits * 10:
                base_thresh = thresh
                break
        
        if base_thresh is None:
            base_thresh = 0.25  # fallback
        
        print(f"\nBase threshold (equal): {base_thresh:.0%}")
        
        # Now create differential based on relative CDR performance
        # CDR3 should be highest (most important)
        # Adjust based on how much "room" each CDR has
        
        # Calculate how much to adjust each CDR relative to base
        # If a CDR has high max, we can be more stringent
        # If a CDR has low max, we need to be more lenient
        
        cdr_maxes = [max_cdr1, max_cdr2, max_cdr3]
        avg_max = np.mean(cdr_maxes)
        
        # Adjustments: if CDR max is above average, we can increase threshold
        # If below average, we need to decrease threshold
        adjustment = 0.03  # 3% adjustment range
        
        diff_cdr1 = base_thresh
        diff_cdr2 = base_thresh
        diff_cdr3 = base_thresh
        
        # CDR3 gets priority (increase if possible)
        if max_cdr3 > avg_max:
            diff_cdr3 = min(base_thresh + adjustment, max_cdr3 * 0.8)  # Don't go above 80% of max
        
        # CDR2 moderate
        if max_cdr2 > avg_max * 0.9:
            diff_cdr2 = min(base_thresh + adjustment * 0.5, max_cdr2 * 0.8)
        elif max_cdr2 < avg_max * 0.9:
            diff_cdr2 = max(base_thresh - adjustment * 0.5, 0.10)
        
        # CDR1 most lenient
        if max_cdr1 < avg_max:
            diff_cdr1 = max(base_thresh - adjustment, 0.10)
        
        # Ensure CDR3 >= CDR2 >= CDR1
        if diff_cdr2 > diff_cdr3:
            diff_cdr2 = diff_cdr3 - 0.01
        if diff_cdr1 > diff_cdr2:
            diff_cdr1 = diff_cdr2 - 0.01
        
        # Round to nearest 1%
        diff_cdr1 = round(diff_cdr1 * 100) / 100
        diff_cdr2 = round(diff_cdr2 * 100) / 100
        diff_cdr3 = round(diff_cdr3 * 100) / 100
        
        # Test this configuration
        diff_count = len(combined_df[
            (combined_df['id_cdr1'] >= diff_cdr1) &
            (combined_df['id_cdr2'] >= diff_cdr2) &
            (combined_df['id_cdr3'] >= diff_cdr3)
        ])
        diff_expected = int(diff_count * scale_factor)
        
        print(f"\nOptimal Differential Configuration:")
        print(f"  CDR1 >= {diff_cdr1:.0%}")
        print(f"  CDR2 >= {diff_cdr2:.0%}")
        print(f"  CDR3 >= {diff_cdr3:.0%}")
        print(f"  Expected hits: ~{diff_expected:,}")
        
        if diff_expected > 0:
            recommendations.append(("Differential", diff_cdr1, diff_cdr2, diff_cdr3, diff_expected))
        
        # Add best equal configs
        for t1, t2, t3, exp in equal_configs[:3]:  # Top 3 equal configs
            recommendations.append((f"Equal {t1:.0%}", t1, t2, t3, exp))
    
    # If using full sequence
    elif 'full' in regions:
        print("Testing full sequence thresholds:\n")
        print(f"{'Threshold':<15} {'Sample':<10} {'Full DB Est.':<15} {'Status'}")
        print("-" * 60)
        
        for thresh in np.arange(0.10, 0.65, 0.02):
            count = len(combined_df[combined_df['id_full'] >= thresh])
            expected = int(count * scale_factor)
            
            if expected == 0:
                status = "❌ Zero"
            elif target_hits * 0.5 <= expected <= target_hits * 5:
                status = "✅ GOOD"
                recommendations.append((f"Full >= {thresh:.0%}", thresh, None, None, expected))
            elif expected < target_hits * 0.5:
                status = "⚠️  Few"
            else:
                status = "⚠️  Many"
            
            print(f"Full >= {thresh:.0%}    {count:<10} ~{expected:<14,} {status}")
    
    # Show top recommendations
    if recommendations:
        print(f"\n{'='*100}")
        print(f"🎯 RECOMMENDED CONFIGURATION (closest to {target_hits:,} hits)")
        print(f"{'='*100}\n")
        
        # Sort by how close to target
        recommendations.sort(key=lambda x: abs(x[4] - target_hits))
        
        # Show best one
        name, t1, t2, t3, expected = recommendations[0]
        
        print(f"Configuration: {name}")
        if t2 is not None and t3 is not None:
            # CDR thresholds
            print(f"  CDR1 >= {t1:.0%}")
            print(f"  CDR2 >= {t2:.0%}")
            print(f"  CDR3 >= {t3:.0%}")
        else:
            # Full sequence threshold
            print(f"  Full sequence >= {t1:.0%}")
        
        print(f"  Expected hits: ~{expected:,}")
        print()
    
    return recommendations


# ============================================================================
# INTERACTIVE QUESTIONS
# ============================================================================

def ask_question(prompt: str, options: List[str]) -> str:
    """Ask user to choose from options."""
    print(f"\n{prompt}")
    for i, opt in enumerate(options, 1):
        print(f"  {i}. {opt}")
    
    while True:
        try:
            choice = input("\nYour choice (number): ").strip()
            idx = int(choice) - 1
            if 0 <= idx < len(options):
                return options[idx]
            print(f"Please enter a number between 1 and {len(options)}")
        except (ValueError, KeyboardInterrupt, EOFError):
            print("Invalid input. Please enter a number.")


def ask_integer(prompt: str, default: int, min_val: int = 1) -> int:
    """Ask for integer input with default."""
    while True:
        try:
            response = input(f"\n{prompt} [default: {default}]: ").strip()
            if not response:
                return default
            val = int(response)
            if val >= min_val:
                return val
            print(f"Please enter a number >= {min_val}")
        except (ValueError, KeyboardInterrupt, EOFError):
            print("Invalid input. Please enter a number.")


def interactive_setup():
    """Ask user questions and return configuration."""
    
    print("="*100)
    print("INTERACTIVE COMPREHENSIVE CDR ANALYSIS")
    print("="*100)
    print("\nThis tool will:")
    print("  1. Ask you questions about your search preferences")
    print("  2. Sample random sequences from the database")
    print("  3. Analyze identity distributions")
    print("  4. Recommend optimal thresholds for your target hit count")
    print("="*100)
    
    # Question 1: Sample size
    sample_size = ask_integer(
        "How many sequences to sample per shard?",
        default=1000,
        min_val=10
    )
    
    # Question 2: Target hits
    target_hits = ask_integer(
        "How many hits do you want to find?",
        default=10000,
        min_val=100
    )
    
    # Question 3: Species
    species = ask_question(
        "Which species?",
        ["Human", "Camel", "Mouse", "Rabbit", "Rat", "Rhesus", "Humanised"]
    )
    
    # Set database path based on species
    if species == "Human":
        db_path = "/home/sasenefrem/KA-Search/extracted/oas-paper/oasdb_20230111/Heavy/Human"
    elif species == "Camel":
        db_path = "/home/sasenefrem/KA-Search/extracted/oas-paper/oasdb_20230111/Heavy/Camel"
    else:
        db_path = f"/home/sasenefrem/KA-Search/extracted/oas-paper/oasdb_20230111/Heavy/{species}"
    
    # Question 4: Optimize for
    optimize = ask_question(
        "Optimize for full sequence or CDRs?",
        ["CDRs only", "Full sequence", "Both (CDRs + Full)"]
    )
    
    # Question 5: Which CDRs?
    regions = []
    if optimize in ["CDRs only", "Both (CDRs + Full)"]:
        cdr_choice = ask_question(
            "Which CDRs?",
            ["All 3 CDRs (CDR1 + CDR2 + CDR3)", "CDR3 only", "CDR2 + CDR3", "Custom"]
        )
        
        if cdr_choice == "All 3 CDRs (CDR1 + CDR2 + CDR3)":
            regions = ['cdr1', 'cdr2', 'cdr3']
        elif cdr_choice == "CDR3 only":
            regions = ['cdr3']
        elif cdr_choice == "CDR2 + CDR3":
            regions = ['cdr2', 'cdr3']
        else:  # Custom
            if input("\nInclude CDR1? (y/n): ").strip().lower() == 'y':
                regions.append('cdr1')
            if input("Include CDR2? (y/n): ").strip().lower() == 'y':
                regions.append('cdr2')
            if input("Include CDR3? (y/n): ").strip().lower() == 'y':
                regions.append('cdr3')
    
    if optimize in ["Full sequence", "Both (CDRs + Full)"]:
        regions.append('full')
    
    # Question 6: Query sequence
    print("\n" + "="*100)
    query_seq = input("Enter your query sequence: ").strip()
    
    return {
        'sample_size': sample_size,
        'target_hits': target_hits,
        'species': species,
        'db_path': db_path,
        'optimize': optimize,
        'regions': regions,
        'query_seq': query_seq
    }


# ============================================================================
# MAIN ANALYSIS
# ============================================================================

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description="Interactive Comprehensive CDR Analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Interactive mode (asks questions):
  python3 interactive_comprehensive_cdr_analysis.py
  
  # Command-line mode (skip questions):
  python3 interactive_comprehensive_cdr_analysis.py \\
    --query-seq "EIQLQQSGAELMKPGASVKISCKASGYTFSSYWIEWIKQRP..." \\
    --species Camel \\
    --sample-size 1000 \\
    --target-hits 10000 \\
    --optimize cdrs \\
    --regions cdr1,cdr2,cdr3
        """
    )
    
    parser.add_argument("--query-seq", help="Query sequence (heavy chain amino acids)")
    parser.add_argument("--species", choices=["Human", "Camel", "Mouse", "Rabbit", "Rat", "Rhesus", "Humanised"],
                        help="Species to search")
    parser.add_argument("--sample-size", type=int, help="Sequences to sample per shard (default: 1000)")
    parser.add_argument("--target-hits", type=int, help="Target number of hits (default: 10000)")
    parser.add_argument("--optimize", choices=["cdrs", "full", "both"],
                        help="Optimize for CDRs, full sequence, or both")
    parser.add_argument("--regions", help="Comma-separated regions (e.g., cdr1,cdr2,cdr3)")
    
    args = parser.parse_args()
    
    # Check if we should run in interactive or command-line mode
    if args.query_seq and args.species:
        # Command-line mode - use provided arguments
        config = {}
        config['query_seq'] = args.query_seq
        config['species'] = args.species
        config['sample_size'] = args.sample_size if args.sample_size else 1000
        config['target_hits'] = args.target_hits if args.target_hits else 10000
        
        # Parse optimize and regions
        if args.optimize:
            if args.optimize == "cdrs":
                config['optimize'] = "CDRs only"
            elif args.optimize == "full":
                config['optimize'] = "Full sequence"
            else:
                config['optimize'] = "Both (CDRs + Full)"
        else:
            config['optimize'] = "CDRs only"
        
        if args.regions:
            config['regions'] = args.regions.split(',')
        else:
            # Default based on optimize
            if config['optimize'] == "CDRs only":
                config['regions'] = ['cdr1', 'cdr2', 'cdr3']
            elif config['optimize'] == "Full sequence":
                config['regions'] = ['full']
            else:
                config['regions'] = ['cdr1', 'cdr2', 'cdr3', 'full']
        
        # Set database path
        config['db_path'] = f"/home/sasenefrem/KA-Search/extracted/oas-paper/oasdb_20230111/Heavy/{config['species']}"
        
        print("="*100)
        print("COMMAND-LINE MODE - Using provided arguments")
        print("="*100)
    else:
        # Interactive mode - ask questions
        config = interactive_setup()
    
    print("\n" + "="*100)
    print("CONFIGURATION SUMMARY")
    print("="*100)
    print(f"Sample size:      {config['sample_size']} sequences per shard")
    print(f"Target hits:      {config['target_hits']:,}")
    print(f"Species:          {config['species']}")
    print(f"Database:         {config['db_path']}")
    print(f"Optimize for:     {config['optimize']}")
    print(f"Regions:          {', '.join(config['regions'])}")
    print(f"Query length:     {len(config['query_seq'])} aa")
    print("="*100)
    
    # Extract query CDRs
    print("\n🔄 Extracting query CDRs...")
    query_full = re.sub(r"[^A-Za-z]", "", config['query_seq'].upper())
    query_cdrs = extract_query_cdrs(query_full)
    query_cdrs['full'] = query_full
    
    print(f"  CDR1: {query_cdrs['cdr1']} (len={len(query_cdrs['cdr1'])})")
    print(f"  CDR2: {query_cdrs['cdr2']} (len={len(query_cdrs['cdr2'])})")
    print(f"  CDR3: {query_cdrs['cdr3']} (len={len(query_cdrs['cdr3'])})")
    print(f"  Full: {len(query_full)} aa")
    
    # Find shards
    print(f"\n🔄 Finding shards in {config['db_path']}...")
    
    # First check if path exists
    if not os.path.exists(config['db_path']):
        print(f"\n❌ ERROR: Database path does not exist!")
        print(f"   Path: {config['db_path']}")
        print(f"\nPlease check:")
        print(f"  1. Is the database installed?")
        print(f"  2. Is the path correct?")
        print(f"  3. Do you have permission to access it?")
        return
    
    npz_files = sorted(glob.glob(os.path.join(config['db_path'], "*.npz")))
    total_shards = len(npz_files)
    
    if total_shards == 0:
        print(f"\n❌ No .npz files found in {config['db_path']}")
        print(f"\nDirectory contents:")
        try:
            contents = os.listdir(config['db_path'])
            for item in contents[:20]:
                print(f"   {item}")
            if len(contents) > 20:
                print(f"   ... and {len(contents) - 20} more items")
        except Exception as e:
            print(f"   Cannot list directory: {e}")
        return
    
    print(f"✅ Found {total_shards} shards")
    if total_shards <= 5:
        print(f"   Shard files:")
        for f in npz_files:
            print(f"     - {os.path.basename(f)}")
    print(f"   Total sequences to sample: {total_shards * config['sample_size']:,}\n")
    
    print("="*100)
    print("PROCESSING SHARDS")
    print("="*100)
    
    all_data = []
    all_stats = []
    failed_count = 0
    
    for i, npz_file in enumerate(npz_files, 1):
        # Show what we're processing (but only for first few and last few to avoid spam)
        if i <= 3 or i >= total_shards - 2 or i % 50 == 0:
            print(f"\nProcessing shard {i}/{total_shards}: {os.path.basename(npz_file)}")
        
        df, stats = analyze_shard_comprehensive(
            npz_file, query_cdrs, config['sample_size'], i, total_shards, config['regions']
        )
        if df is not None:
            all_data.append(df)
            all_stats.append(stats)
        else:
            failed_count += 1
            if failed_count <= 5:  # Show first 5 failures
                print(f"   ⚠️  Shard {i} failed to process")
    
    print(f"\n{'='*100}")
    print(f"SHARD PROCESSING SUMMARY")
    print(f"{'='*100}")
    print(f"Total shards found: {total_shards}")
    print(f"Successfully processed: {len(all_data)}")
    print(f"Failed: {failed_count}")
    
    # Check if we got any data
    if len(all_data) == 0:
        print(f"\n{'='*100}")
        print("❌ ERROR: No shards were successfully processed!")
        print(f"{'='*100}")
        print(f"\nPossible issues:")
        print(f"  1. Database files might be empty or corrupted")
        print(f"  2. Database files might have different structure than expected")
        print(f"  3. You might need to re-download/re-extract the database")
        print(f"\nDatabase path tried: {config['db_path']}")
        print(f"\nTry running this to inspect a shard file:")
        if total_shards > 0:
            print(f"  python3 inspect_shard.py {npz_files[0]}")
        print()
        return
    
    # Combine all data
    combined = pd.concat(all_data, ignore_index=True)
    stats_df = pd.DataFrame(all_stats)
    
    # Save raw data
    output_csv = f"{config['species'].lower()}_comprehensive_results.csv"
    output_stats = f"{config['species'].lower()}_per_shard_stats.csv"
    
    combined.to_csv(output_csv, index=False)
    stats_df.to_csv(output_stats, index=False)
    
    print(f"\n{'='*100}")
    print(f"OVERALL SUMMARY ({len(combined):,} sequences from {total_shards} shards)")
    print(f"{'='*100}\n")
    
    # Overall statistics
    print("MAXIMUM IDENTITIES FOUND:")
    for region in config['regions']:
        col = f'id_{region}'
        if col in combined.columns:
            max_val = combined[col].max()
            print(f"  {region.upper():10s}: {max_val:.4f} ({max_val*100:.2f}%)")
    
    if 'cdr_avg' in combined.columns:
        print(f"  {'CDR Avg':10s}: {combined['cdr_avg'].max():.4f} ({combined['cdr_avg'].max()*100:.2f}%)")
        print(f"  {'CDR Min':10s}: {combined['cdr_min'].max():.4f} ({combined['cdr_min'].max()*100:.2f}%)")
    
    print("\nMEAN IDENTITIES:")
    for region in config['regions']:
        col = f'id_{region}'
        if col in combined.columns:
            mean_val = combined[col].mean()
            print(f"  {region.upper():10s}: {mean_val:.4f} ({mean_val*100:.2f}%)")
    
    # Calculate optimal thresholds (moved before top 20)
    recommendations = calculate_optimal_thresholds(
        combined, stats_df, config['target_hits'], config['regions']
    )
    
    # Show top matches at the end
    if 'cdr_avg' in combined.columns:
        print(f"\n{'='*100}")
        print("TOP 20 SEQUENCES (by average of selected CDRs)")
        print(f"{'='*100}")
        
        top20 = combined.nlargest(20, "cdr_avg")
        
        for i, (_, row) in enumerate(top20.iterrows(), 1):
            print(f"\n#{i:2d} Shard: {row['shard'][:50]:50s} (idx {row['index']:,})")
            if 'id_cdr1' in row:
                print(f"    CDR1: {row['cdr1_seq']:20s} → {row['id_cdr1']*100:5.2f}%")
            if 'id_cdr2' in row:
                print(f"    CDR2: {row['cdr2_seq']:20s} → {row['id_cdr2']*100:5.2f}%")
            if 'id_cdr3' in row:
                print(f"    CDR3: {row['cdr3_seq']:20s} → {row['id_cdr3']*100:5.2f}%")
            if 'cdr_avg' in row:
                print(f"    AVG:  {row['cdr_avg']*100:5.2f}%", end="")
            if 'cdr_min' in row:
                print(f"  |  MIN: {row['cdr_min']*100:5.2f}%", end="")
            if 'id_full' in row:
                print(f"  |  FULL: {row['id_full']*100:5.2f}%")
            else:
                print()
    
    print(f"\n{'='*100}")
    print(f"✅ Analysis complete!")
    print(f"   Raw data:       {output_csv}")
    print(f"   Shard stats:    {output_stats}")
    print(f"{'='*100}\n")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\n\n❌ Analysis interrupted by user")
    except Exception as e:
        print(f"\n\n❌ Error: {e}")
        import traceback
        traceback.print_exc()