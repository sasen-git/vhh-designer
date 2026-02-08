#!/usr/bin/env python3
"""
Convert Epistasis Analysis Output to VHH Designer Format

This script converts the epistasis_v2_full.pkl output into the format expected
by vhh_designer_v2.py, creating properly structured:
- multi_position_rules (from analysis_4_higher_order_rules)
- vernier_archetypes (from analysis_2_vernier_clusters)

Usage:
    python convert_epistasis_to_designer.py \
        --epistasis epistasis_v2_full.pkl \
        --output correlation_results_v3_compensation.pkl \
        [--existing correlation_results_v3.pkl]  # Optional: merge with existing

Author: Claude (Anthropic)
Date: January 2025
"""

import os
import re
import pickle
import argparse
from collections import defaultdict
from typing import Dict, List, Tuple, Any


def parse_result_string(result: str) -> Tuple[str, str, str]:
    """
    Parse epistasis result string like 'FR3_6=V' into components.
    
    Returns: (region, index, amino_acid) or (None, None, None) if parse fails
    """
    match = re.match(r'(FR\d)_(\d+)=(\w)', result)
    if match:
        return match.groups()
    return None, None, None


def parse_condition_string(condition: str) -> Dict[str, Any]:
    """
    Parse epistasis condition string like 'cdr3_len=short AND FR2_12=C'
    or new IMGT format 'cdr3_len_cat=short AND IMGT49=G'
    
    Returns dict with parsed conditions
    """
    parsed = {
        'cdr_features': [],
        'fr_requirements': [],
        'raw': condition
    }
    
    # Split on AND
    parts = [p.strip() for p in condition.split(' AND ')]
    
    for part in parts:
        # CDR length conditions (both old and new format)
        if 'cdr3_len=' in part or 'cdr3_len_cat=' in part:
            parsed['cdr_features'].append(('cdr3_length', part.split('=')[1]))
        elif 'cdr2_len=' in part:
            parsed['cdr_features'].append(('cdr2_length', part.split('=')[1]))
        elif 'cdr1_len=' in part:
            parsed['cdr_features'].append(('cdr1_length', part.split('=')[1]))
        # CDR charge conditions
        elif 'cdr3_charge=' in part or 'cdr3_charge_cat=' in part:
            parsed['cdr_features'].append(('cdr3_charge', part.split('=')[1]))
        # CDR terminal residues
        elif 'cdr3_last=' in part:
            parsed['cdr_features'].append(('cdr3_last', part.split('=')[1]))
        elif 'cdr3_first=' in part:
            parsed['cdr_features'].append(('cdr3_first', part.split('=')[1]))
        # IMGT position requirements (new format)
        elif re.match(r'IMGT\d+=\w', part):
            match = re.match(r'(IMGT\d+)=(\w)', part)
            if match:
                position, aa = match.groups()
                parsed['fr_requirements'].append({
                    'position': position,
                    'required_aa': aa
                })
        # FR position requirements (old format)
        elif re.match(r'FR\d_\d+=\w', part):
            match = re.match(r'(FR\d)_(\d+)=(\w)', part)
            if match:
                region, idx, aa = match.groups()
                parsed['fr_requirements'].append({
                    'position': f'{region}[{idx}]',
                    'required_aa': aa
                })
    
    return parsed


def convert_higher_order_rules(rules: List[Dict]) -> List[Dict]:
    """
    Convert analysis_4_higher_order_rules to multi_position_rules format.
    
    Handles both old format (FR3_6=V) and new IMGT format (IMGT69=V).
    
    Input format (old):
        {
            'family': 'Classical_other',
            'condition': 'cdr3_len=short AND FR2_12=C',
            'result': 'FR3_6=V',
            'support': 3184,
            'confidence': 0.986,
            'count': 3140
        }
    
    Input format (new IMGT):
        {
            'family': 'Classical_other',
            'condition': 'cdr3_len_cat=short AND IMGT49=G',
            'result': 'IMGT69=V',
            'support': 3184,
            'confidence': 0.986,
            'count': 3140
        }
    
    Output format:
        {
            'fw_position': 'IMGT69' or 'FR3[6]',
            'suggested_aa': 'V',
            'confidence': 98.6,
            'condition': {...parsed...},
            'family': 'Classical_other',
            'support': 3184,
            'source': 'epistasis'
        }
    """
    converted = []
    seen = set()  # Avoid duplicates
    
    for rule in rules:
        result_str = rule.get('result', '')
        
        # Try IMGT format first (IMGT69=V)
        imgt_match = re.match(r'(IMGT\d+)=(\w)', result_str)
        if imgt_match:
            fw_position = imgt_match.group(1)
            aa = imgt_match.group(2)
        else:
            # Try old format (FR3_6=V)
            old_match = re.match(r'(FR\d)_(\d+)=(\w)', result_str)
            if old_match:
                region, idx, aa = old_match.groups()
                fw_position = f'{region}[{idx}]'
            else:
                continue
        
        # Parse condition
        condition = parse_condition_string(rule.get('condition', ''))
        
        # Create unique key to avoid duplicates
        key = (fw_position, aa, rule.get('family', ''), rule.get('condition', ''))
        if key in seen:
            continue
        seen.add(key)
        
        converted.append({
            'fw_position': fw_position,
            'suggested_aa': aa,
            'confidence': rule.get('confidence', 0.5) * 100,  # Convert to percentage
            'condition': condition,
            'condition_raw': rule.get('condition', ''),
            'family': rule.get('family', 'all'),
            'support': rule.get('support', 0),
            'count': rule.get('count', 0),
            'source': 'epistasis'
        })
    
    # Sort by confidence (descending) then support (descending)
    converted.sort(key=lambda x: (-x['confidence'], -x['support']))
    
    return converted


def convert_vernier_clusters(clusters: Dict, min_support: int = 10000, top_n_per_family: int = 5) -> Dict[str, List[Dict]]:
    """
    Convert analysis_2_vernier_clusters to vernier_archetypes format.
    
    Handles both old format (FR2_4, FR3_6) and new IMGT format (IMGT42, IMGT69).
    
    Input format (old):
        ('Q', 'E', 'F', 'A', 'V', 'G', 'I', 'N'): {
            'n': 590284,
            'pattern': {
                'FR2_4': 'Q', 'FR2_11': 'E', ...
            },
            'families': {'F_C2': 362842, ...},
        }
    
    Input format (new IMGT):
        ('Q', 'E', 'F', 'A', 'V', 'G', 'I', 'N'): {
            'n': 590284,
            'pattern': {
                'IMGT42': 'Q', 'IMGT47': 'E', ...
            },
            'families': {'F_C2': 362842, ...},
        }
    
    Output format:
        {
            'F_C2': [
                {'fw_position': 'IMGT42', 'suggested_aa': 'Q', 'confidence': 80.0, 'source': 'vernier'},
                ...
            ],
        }
    """
    # Collect clusters by dominant family
    family_clusters = defaultdict(list)
    
    for pattern_tuple, cluster_data in clusters.items():
        n = cluster_data.get('n', 0)
        if n < min_support:
            continue
        
        families = cluster_data.get('families', {})
        if not families:
            continue
        
        # Find dominant family
        dominant_family = max(families, key=families.get)
        dominant_count = families[dominant_family]
        
        # Calculate family purity (what % of this cluster is the dominant family)
        purity = dominant_count / n if n > 0 else 0
        
        family_clusters[dominant_family].append({
            'pattern': cluster_data.get('pattern', {}),
            'n': n,
            'purity': purity,
            'cdr3_stats': cluster_data.get('cdr3_length', {}),
            'cdr3_charge': cluster_data.get('cdr3_charge', {}),
        })
    
    # For each family, get top clusters and extract rules
    archetypes = {}
    
    target_families = ['F_C2', 'Y_C2', 'F_C4', 'Y_C4', 'Non_classical', 'Classical_other', 'VH_like']
    
    for family in target_families:
        if family not in family_clusters:
            continue
        
        # Sort by support * purity
        sorted_clusters = sorted(
            family_clusters[family],
            key=lambda x: x['n'] * x['purity'],
            reverse=True
        )
        
        # Take top N clusters
        top_clusters = sorted_clusters[:top_n_per_family]
        
        # Extract position rules from top cluster (most representative)
        if top_clusters:
            best = top_clusters[0]
            rules = []
            
            for pos_key, aa in best['pattern'].items():
                # Check if already IMGT format
                if pos_key.startswith('IMGT'):
                    fw_position = pos_key
                else:
                    # Convert FR2_4 -> FR2[4] (old format)
                    match = re.match(r'(FR\d)_(\d+)', pos_key)
                    if match:
                        region, idx = match.groups()
                        fw_position = f'{region}[{idx}]'
                    else:
                        continue
                
                rules.append({
                    'fw_position': fw_position,
                    'suggested_aa': aa,
                    'confidence': 80.0,  # Vernier positions are stable
                    'support': best['n'],
                    'purity': best['purity'],
                    'source': 'vernier'
                })
            
            archetypes[family] = rules
            
            # Also store cluster metadata
            archetypes[f'{family}_metadata'] = {
                'top_cluster_support': best['n'],
                'top_cluster_purity': best['purity'],
                'cdr3_mean_length': best['cdr3_stats'].get('mean', 0),
                'cdr3_std_length': best['cdr3_stats'].get('std', 0),
                'n_clusters': len(family_clusters[family]),
            }
    
    return archetypes


def extract_cdr_conditional_rules(rules: List[Dict]) -> List[Dict]:
    """
    Extract rules that have clear CDR conditions (not just FR requirements).
    These are the most useful for CDR-dependent framework selection.
    """
    cdr_conditional = []
    
    for rule in rules:
        condition = rule.get('condition', {})
        if isinstance(condition, str):
            condition = parse_condition_string(condition)
        
        # Check if there are CDR-related conditions
        cdr_features = condition.get('cdr_features', [])
        if cdr_features:
            cdr_conditional.append(rule)
    
    return cdr_conditional


def create_simplified_rules(converted_rules: List[Dict], vernier_archetypes: Dict) -> List[Dict]:
    """
    Create simplified rules that are easy to apply:
    - Filter to high-confidence rules (>70%)
    - Group by position
    - Include both epistasis and vernier sources
    """
    simplified = []
    position_rules = defaultdict(list)
    
    # Add epistasis rules
    for rule in converted_rules:
        conf = rule.get('confidence', 0)
        if conf >= 70:
            pos = rule.get('fw_position', '')
            if pos:
                position_rules[pos].append(rule)
    
    # Add vernier rules
    for family, rules in vernier_archetypes.items():
        if family.endswith('_metadata'):
            continue
        if not isinstance(rules, list):
            continue
        for rule in rules:
            pos = rule.get('fw_position', '')
            if pos:
                rule_with_family = rule.copy()
                rule_with_family['family'] = family
                position_rules[pos].append(rule_with_family)
    
    # Create simplified list
    for pos, rules in position_rules.items():
        # Sort by confidence (with default 0 for missing)
        rules.sort(key=lambda x: -x.get('confidence', 0))
        simplified.extend(rules)
    
    return simplified


def merge_with_existing(existing_path: str, new_rules: Dict) -> Dict:
    """
    Merge new rules with existing correlation results.
    Existing rules take precedence if there are conflicts.
    """
    with open(existing_path, 'rb') as f:
        existing = pickle.load(f)
    
    merged = existing.copy()
    
    # Merge multi_position_rules
    existing_rules = merged.get('multi_position_rules', [])
    new_multi = new_rules.get('multi_position_rules', [])
    
    # Create set of existing rule keys
    existing_keys = set()
    for r in existing_rules:
        key = (r.get('fw_position'), r.get('suggested_aa'))
        existing_keys.add(key)
    
    # Add new rules that don't conflict
    for r in new_multi:
        key = (r.get('fw_position'), r.get('suggested_aa'))
        if key not in existing_keys:
            existing_rules.append(r)
    
    merged['multi_position_rules'] = existing_rules
    
    # Merge vernier_archetypes
    existing_vernier = merged.get('vernier_archetypes', {})
    new_vernier = new_rules.get('vernier_archetypes', {})
    
    for family, rules in new_vernier.items():
        if family not in existing_vernier:
            existing_vernier[family] = rules
    
    merged['vernier_archetypes'] = existing_vernier
    
    # Preserve cdr_conditional_rules from new (epistasis) data
    if 'cdr_conditional_rules' in new_rules:
        merged['cdr_conditional_rules'] = new_rules['cdr_conditional_rules']
    
    # Preserve conversion_metadata from new data
    if 'conversion_metadata' in new_rules:
        merged['conversion_metadata'] = new_rules['conversion_metadata']
    
    return merged


def print_summary(output_data: Dict):
    """Print summary of converted data."""
    print("\n" + "=" * 70)
    print("CONVERSION SUMMARY")
    print("=" * 70)
    
    # Multi-position rules
    rules = output_data.get('multi_position_rules', [])
    print(f"\nMulti-position rules: {len(rules)}")
    
    if rules:
        # Count by region
        region_counts = defaultdict(int)
        for r in rules:
            pos = r.get('fw_position', '')
            region = pos.split('[')[0] if '[' in pos else 'unknown'
            region_counts[region] += 1
        
        print("  By region:")
        for region, count in sorted(region_counts.items()):
            print(f"    {region}: {count}")
        
        # Confidence distribution
        conf_ranges = {'90-100': 0, '80-90': 0, '70-80': 0, '<70': 0}
        for r in rules:
            conf = r.get('confidence', 0)
            if conf >= 90:
                conf_ranges['90-100'] += 1
            elif conf >= 80:
                conf_ranges['80-90'] += 1
            elif conf >= 70:
                conf_ranges['70-80'] += 1
            else:
                conf_ranges['<70'] += 1
        
        print("  By confidence:")
        for range_name, count in conf_ranges.items():
            print(f"    {range_name}%: {count}")
    
    # Vernier archetypes
    vernier = output_data.get('vernier_archetypes', {})
    families = [k for k in vernier.keys() if not k.endswith('_metadata')]
    print(f"\nVernier archetypes: {len(families)} families")
    
    for family in families:
        rules = vernier.get(family, [])
        metadata = vernier.get(f'{family}_metadata', {})
        support = metadata.get('top_cluster_support', 0)
        purity = metadata.get('top_cluster_purity', 0)
        cdr3_mean = metadata.get('cdr3_mean_length', 0)
        print(f"  {family}: {len(rules)} positions, support={support:,}, purity={purity:.1%}, CDR3~{cdr3_mean:.1f}aa")
    
    # CDR-conditional rules
    cdr_rules = output_data.get('cdr_conditional_rules', [])
    print(f"\nCDR-conditional rules: {len(cdr_rules)}")
    
    print("=" * 70)


def main():
    parser = argparse.ArgumentParser(
        description='Convert epistasis analysis output to VHH designer format',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Basic conversion
    python convert_epistasis_to_designer.py \\
        --epistasis epistasis_v2_full.pkl \\
        --output correlation_results_v3_compensation.pkl
    
    # Merge with existing rules
    python convert_epistasis_to_designer.py \\
        --epistasis epistasis_v2_full.pkl \\
        --existing correlation_results_v3.pkl \\
        --output correlation_results_v3_compensation.pkl
    
    # Custom thresholds
    python convert_epistasis_to_designer.py \\
        --epistasis epistasis_v2_full.pkl \\
        --output correlation_results_v3_compensation.pkl \\
        --min-support 5000 \\
        --min-confidence 75
        """
    )
    
    parser.add_argument('--epistasis', '-e', required=True,
                        help='Input epistasis pickle file (epistasis_v2_full.pkl)')
    parser.add_argument('--output', '-o', required=True,
                        help='Output correlation results pickle file')
    parser.add_argument('--existing', '-x',
                        help='Optional: existing correlation results to merge with')
    parser.add_argument('--min-support', type=int, default=10000,
                        help='Minimum support for Vernier clusters (default: 10000)')
    parser.add_argument('--min-confidence', type=float, default=70.0,
                        help='Minimum confidence for rules (default: 70.0)')
    parser.add_argument('--top-clusters', type=int, default=5,
                        help='Top N Vernier clusters per family (default: 5)')
    parser.add_argument('--verbose', '-v', action='store_true',
                        help='Print detailed conversion info')
    
    args = parser.parse_args()
    
    # Load epistasis data
    print(f"Loading epistasis data from: {args.epistasis}")
    with open(args.epistasis, 'rb') as f:
        epistasis_data = pickle.load(f)
    
    print(f"  Keys: {list(epistasis_data.keys())}")
    
    # Convert higher-order rules (handle both v2 and v3 key names)
    print("\nConverting higher-order rules...")
    higher_order = epistasis_data.get('analysis_4_higher_order_rules', 
                   epistasis_data.get('higher_order_rules', []))
    print(f"  Input rules: {len(higher_order)}")
    
    converted_rules = convert_higher_order_rules(higher_order)
    print(f"  Converted rules: {len(converted_rules)}")
    
    # Filter by confidence
    filtered_rules = [r for r in converted_rules if r['confidence'] >= args.min_confidence]
    print(f"  After confidence filter (>={args.min_confidence}%): {len(filtered_rules)}")
    
    # Convert Vernier clusters (handle both v2 and v3 key names)
    print("\nConverting Vernier clusters...")
    vernier_clusters = epistasis_data.get('analysis_2_vernier_clusters',
                       epistasis_data.get('vernier_clusters', {}))
    print(f"  Input clusters: {len(vernier_clusters)}")
    
    vernier_archetypes = convert_vernier_clusters(
        vernier_clusters,
        min_support=args.min_support,
        top_n_per_family=args.top_clusters
    )
    
    families = [k for k in vernier_archetypes.keys() if not k.endswith('_metadata')]
    print(f"  Families extracted: {families}")
    
    # Extract CDR-conditional rules
    print("\nExtracting CDR-conditional rules...")
    cdr_conditional = extract_cdr_conditional_rules(filtered_rules)
    print(f"  CDR-conditional rules: {len(cdr_conditional)}")
    
    # Create output data structure
    output_data = {
        'multi_position_rules': filtered_rules,
        'vernier_archetypes': vernier_archetypes,
        'cdr_conditional_rules': cdr_conditional,
        'conversion_metadata': {
            'source_file': args.epistasis,
            'min_support': args.min_support,
            'min_confidence': args.min_confidence,
            'total_input_rules': len(higher_order),
            'total_input_clusters': len(vernier_clusters),
        }
    }
    
    # Merge with existing if provided
    if args.existing and os.path.exists(args.existing):
        print(f"\nMerging with existing rules from: {args.existing}")
        output_data = merge_with_existing(args.existing, output_data)
    
    # Create simplified rules for easy lookup
    output_data['simplified_rules'] = create_simplified_rules(
        output_data['multi_position_rules'],
        output_data['vernier_archetypes']
    )
    
    # Print summary
    print_summary(output_data)
    
    # Save output
    print(f"\nSaving to: {args.output}")
    with open(args.output, 'wb') as f:
        pickle.dump(output_data, f)
    
    print("Done!")
    
    # Verbose output
    if args.verbose:
        print("\n" + "=" * 70)
        print("SAMPLE RULES")
        print("=" * 70)
        
        print("\nTop 10 multi-position rules:")
        for i, rule in enumerate(output_data['multi_position_rules'][:10]):
            print(f"  {i+1}. {rule['fw_position']} -> {rule['suggested_aa']} "
                  f"(conf={rule['confidence']:.1f}%, family={rule.get('family', 'all')})")
            print(f"      Condition: {rule.get('condition_raw', 'N/A')}")
        
        print("\nVernier archetype positions (F_C2):")
        if 'F_C2' in vernier_archetypes:
            for rule in vernier_archetypes['F_C2']:
                print(f"  {rule['fw_position']} -> {rule['suggested_aa']} "
                      f"(support={rule.get('support', 0):,})")


if __name__ == '__main__':
    main()