#!/usr/bin/env python3
"""
VHH Framework Optimizer v2 - Fixed CDR extraction and scoring

Key fixes from v1:
1. Better CDR extraction (excludes FR1 anchor from CDR1)
2. Higher default confidence threshold (90%)
3. Ignores cdr1[0:3] rules (often contaminated by FR anchor)
4. Added user's custom scaffold
"""

import argparse
import pickle
import sys
import csv
import re
from collections import defaultdict, Counter
from dataclasses import dataclass
from typing import List, Dict, Tuple, Optional

# ============================================================
# SCAFFOLDS - Including user's custom scaffold
# ============================================================

SCAFFOLDS = {
    'humanized': {
        'name': 'Humanized VHH (VH3-based)',
        'FR1': 'EVQLVESGGGLVQPGGSLRLSCAAS',
        'FR2': 'WVRQAPGKGLEWVS',
        'FR3': 'RFTISRDNSKNTLYLQMNSLRAEDTAVYYC',
        'FR4': 'WGQGTLVTVSS',
        'family': 'humanized',
    },
    'custom_hybrid': {
        'name': 'Custom Humanized Hybrid (WFR+GLEWA)',
        'FR1': 'QVQLVESGGGLVQPGGSLRLSC',
        'FR2': 'WFRQAPGQGLEAVA',
        'FR3': 'KGRFTISRDNSKNTLYLQMNSLRAEDTAVYYC',
        'FR4': 'WGQGTLVTVSS',
        'family': 'humanized',
    },
    'VHH1_llama': {
        'name': 'Llama VHH1 (F-E-R-G)',
        'FR1': 'QVQLQESGGGLVQAGGSLRLSCAAS',
        'FR2': 'WFRQAPGKEREGVS',
        'FR3': 'RFTISRDNAKNTVYLQMNSLKPEDTAVYYC',
        'FR4': 'WGQGTQVTVSS',
        'family': 'F_C2',
    },
    'F_C2_alpaca': {
        'name': 'Alpaca IGHV3-3 (F-E-R-F)',
        'FR1': 'QVQLVESGGGLVQPGGSLRLSCAAS',
        'FR2': 'WFRQAPGKEREFVS',
        'FR3': 'RFTISRDSAKNTVYLQMNSLRAEDTAVYYC',
        'FR4': 'WGQGTLVTVSS',
        'family': 'F_C2',
    },
    'Y_C2_alpaca': {
        'name': 'Alpaca IGHV3S53 (Y-E-R-L)',
        'FR1': 'EVQLVESGGGLVQPGGSLRLSCAAS',
        'FR2': 'WYRQAPGKQRELVA',
        'FR3': 'RFTISRDNAKNTVYLQMNSLRAEDTAVYYC',
        'FR4': 'WGQGTLVTVSS',
        'family': 'Y_C2',
    },
}

# ============================================================
# CDR EXTRACTION - Fixed to exclude FR1 anchor
# ============================================================

def extract_cdrs(sequence: str) -> Dict[str, str]:
    """
    Extract CDRs using conserved anchors.
    CDR1 starts AFTER the Cys-Anchor motif, not including it.
    """
    seq = sequence.upper().replace(' ', '').replace('\n', '')
    
    result = {
        'FR1': '', 'CDR1': '', 'FR2': '', 'CDR2': '', 
        'FR3': '', 'CDR3': '', 'FR4': '',
        'valid': False, 'error': ''
    }
    
    if len(seq) < 100:
        result['error'] = 'Too short'
        return result
    
    # Find conserved Cys (~position 22)
    cys1_pos = -1
    for i in range(18, 28):
        if i < len(seq) and seq[i] == 'C':
            cys1_pos = i
            break
    if cys1_pos == -1:
        cys1_pos = 22
    
    # FR1 ends with C (inclusive)
    result['FR1'] = seq[:cys1_pos + 1]
    
    # Find FR2 start (W[YFV]RQ pattern)
    fr2_start = -1
    for pattern in [r'W[YFVIL]RQ', r'W.RQ']:
        m = re.search(pattern, seq[cys1_pos+1:cys1_pos+25])
        if m:
            fr2_start = cys1_pos + 1 + m.start()
            break
    
    if fr2_start == -1:
        # Default: CDR1 is ~8-12 residues after Cys
        fr2_start = cys1_pos + 10
    
    # CDR1: between Cys and FR2-W
    result['CDR1'] = seq[cys1_pos + 1:fr2_start]
    
    # FR2 is ~14 aa
    fr2_end = fr2_start + 14
    result['FR2'] = seq[fr2_start:fr2_end]
    
    # Find FR4 start (WG pattern near end)
    fr4_start = -1
    for pattern in [r'WG[QR]GT', r'WGQG', r'WG']:
        m = re.search(pattern, seq[-20:])
        if m:
            fr4_start = len(seq) - 20 + m.start()
            break
    if fr4_start == -1:
        fr4_start = len(seq) - 11
    
    result['FR4'] = seq[fr4_start:]
    
    # Find Cys before CDR3
    cys2_pos = -1
    for i in range(fr4_start - 25, fr4_start - 3):
        if 0 <= i < len(seq) and seq[i] == 'C':
            cys2_pos = i
            break
    if cys2_pos == -1:
        yyc = seq.rfind('YYC', fr2_end + 20, fr4_start)
        cys2_pos = yyc + 2 if yyc != -1 else fr4_start - 15
    
    result['CDR3'] = seq[cys2_pos + 1:fr4_start]
    
    # CDR2 and FR3
    remaining = cys2_pos - fr2_end + 1
    # CDR2 is typically 8-17 aa, FR3 is 32-38 aa
    cdr2_len = max(5, min(17, remaining - 32))
    
    result['CDR2'] = seq[fr2_end:fr2_end + cdr2_len]
    result['FR3'] = seq[fr2_end + cdr2_len:cys2_pos + 1]
    
    # Validate
    total = sum(len(result[r]) for r in ['FR1', 'CDR1', 'FR2', 'CDR2', 'FR3', 'CDR3', 'FR4'])
    result['valid'] = (total == len(seq))
    
    return result

# ============================================================
# SCORING - Fixed with high confidence and exclusions
# ============================================================

def get_predictions(cdr1: str, cdr2: str, cdr3: str, data: Dict, family: str, 
                   min_conf: float = 90.0, exclude_early_cdr1: bool = True) -> Dict:
    """
    Get framework predictions.
    
    Args:
        exclude_early_cdr1: Skip cdr1[0], cdr1[1], cdr1[2] rules (often contaminated)
    """
    predictions = {}
    cdr_deps = data.get('family_cdr_dependent', {}).get(family, {})
    
    # Extract features
    features = {}
    for cdr_name, cdr_seq in [('cdr1', cdr1), ('cdr2', cdr2), ('cdr3', cdr3)]:
        if not cdr_seq:
            continue
        for pos in [0, 1, 2, -3, -2, -1]:
            try:
                if pos >= 0:
                    if pos < len(cdr_seq):
                        features[f'{cdr_name}[{pos}]'] = cdr_seq[pos]
                else:
                    if abs(pos) <= len(cdr_seq):
                        features[f'{cdr_name}[{pos}]'] = cdr_seq[pos]
            except:
                pass
    
    for fw_pos, info in cdr_deps.items():
        for rule_key, rule_details in info.get('rules', {}).items():
            # Skip low confidence
            if rule_details['confidence'] < min_conf:
                continue
            
            # Skip early CDR1 rules if requested
            if exclude_early_cdr1:
                if 'cdr1[0]=' in rule_key or 'cdr1[1]=' in rule_key or 'cdr1[2]=' in rule_key:
                    continue
            
            try:
                parts = rule_key.split('=')
                cdr_part, required_aa = parts[0], parts[1]
                
                if features.get(cdr_part) == required_aa:
                    if fw_pos not in predictions or rule_details['confidence'] > predictions[fw_pos]['confidence']:
                        predictions[fw_pos] = {
                            'predicted_aa': rule_details['fw_aa'],
                            'confidence': rule_details['confidence'],
                            'rule': rule_key,
                        }
            except:
                continue
    
    return predictions

def guess_family(cdr1: str, cdr2: str, cdr3: str) -> str:
    """Guess germline family from CDR characteristics."""
    cdr3_len = len(cdr3) if cdr3 else 0
    total_cys = sum(s.count('C') for s in [cdr1, cdr2, cdr3] if s)
    
    if total_cys >= 2 or cdr3_len >= 18:
        return 'F_C4'
    elif cdr3_len >= 14:
        return 'F_C2'
    else:
        return 'Y_C2'

def score_framework(fw: Dict[str, str], predictions: Dict) -> Tuple[float, List[Dict]]:
    """Score framework compatibility."""
    matches, mismatches = [], []
    
    for fw_pos, pred in predictions.items():
        try:
            region = fw_pos.split('[')[0]
            idx = int(fw_pos.split('[')[1].rstrip(']'))
            
            if idx < len(fw.get(region, '')):
                actual = fw[region][idx]
                if actual == pred['predicted_aa']:
                    matches.append({'pos': fw_pos, 'conf': pred['confidence']})
                else:
                    mismatches.append({
                        'position': fw_pos,
                        'current': actual,
                        'suggested': pred['predicted_aa'],
                        'confidence': pred['confidence'],
                        'rule': pred['rule'],
                    })
        except:
            continue
    
    if not predictions:
        return 100.0, []
    
    total_w = sum(m['conf'] for m in matches) + sum(m['confidence'] for m in mismatches)
    match_w = sum(m['conf'] for m in matches)
    score = 100.0 * match_w / total_w if total_w > 0 else 100.0
    
    return score, sorted(mismatches, key=lambda x: -x['confidence'])

def apply_mutations(fw: Dict[str, str], mismatches: List[Dict]) -> Dict[str, str]:
    """Apply mutations to framework."""
    result = {k: list(v) for k, v in fw.items()}
    
    for mm in mismatches:
        try:
            region = mm['position'].split('[')[0]
            idx = int(mm['position'].split('[')[1].rstrip(']'))
            if region in result and idx < len(result[region]):
                result[region][idx] = mm['suggested']
        except:
            continue
    
    return {k: ''.join(v) for k, v in result.items()}

def assemble(fw: Dict, cdr1: str, cdr2: str, cdr3: str) -> str:
    """Assemble full sequence."""
    return fw.get('FR1','') + cdr1 + fw.get('FR2','') + cdr2 + fw.get('FR3','') + cdr3 + fw.get('FR4','')

# ============================================================
# MAIN
# ============================================================

@dataclass
class Result:
    id: str
    original_seq: str
    cdr1: str
    cdr2: str
    cdr3: str
    fr1: str
    fr2: str
    fr3: str
    fr4: str
    valid: bool
    family: str
    orig_score: float
    orig_n_issues: int
    best_scaffold: str
    best_score: float
    best_n_mut: int
    best_mutations: str
    best_grafted: str
    human_n_mut: int
    human_mutations: str
    human_grafted: str
    custom_n_mut: int
    custom_mutations: str
    custom_grafted: str
    error: str = ''

def process(seq_id: str, seq: str, data: Dict, min_conf: float) -> Result:
    """Process one sequence."""
    
    regions = extract_cdrs(seq)
    
    if not regions['CDR1'] or not regions['CDR3']:
        return Result(
            id=seq_id, original_seq=seq,
            cdr1='', cdr2='', cdr3='',
            fr1='', fr2='', fr3='', fr4='',
            valid=False, family='',
            orig_score=0, orig_n_issues=0,
            best_scaffold='', best_score=0, best_n_mut=0, best_mutations='', best_grafted='',
            human_n_mut=0, human_mutations='', human_grafted='',
            custom_n_mut=0, custom_mutations='', custom_grafted='',
            error=regions.get('error', 'CDR extraction failed')
        )
    
    cdr1 = regions['CDR1']
    cdr2 = regions['CDR2']
    cdr3 = regions['CDR3']
    
    family = guess_family(cdr1, cdr2, cdr3)
    
    # Get predictions
    preds = get_predictions(cdr1, cdr2, cdr3, data, family, min_conf)
    
    # Score original
    orig_fw = {'FR1': regions['FR1'], 'FR2': regions['FR2'], 
               'FR3': regions['FR3'], 'FR4': regions['FR4']}
    orig_score, orig_mm = score_framework(orig_fw, preds)
    
    # Score all scaffolds
    results = []
    for name, info in SCAFFOLDS.items():
        scaffold_fw = {k: v for k, v in info.items() if k.startswith('FR')}
        score, mm = score_framework(scaffold_fw, preds)
        results.append((name, score, mm, info))
    
    results.sort(key=lambda x: -x[1])
    best_name, best_score, best_mm, best_info = results[0]
    
    # Graft into best
    best_fw = {k: v for k, v in best_info.items() if k.startswith('FR')}
    best_fw_mut = apply_mutations(best_fw, best_mm)
    best_grafted = assemble(best_fw_mut, cdr1, cdr2, cdr3)
    best_muts = '; '.join([f"{m['position']}:{m['current']}>{m['suggested']}" for m in best_mm])
    
    # Graft into humanized
    hum_info = SCAFFOLDS['humanized']
    hum_fw = {k: v for k, v in hum_info.items() if k.startswith('FR')}
    _, hum_mm = score_framework(hum_fw, preds)
    hum_fw_mut = apply_mutations(hum_fw, hum_mm)
    hum_grafted = assemble(hum_fw_mut, cdr1, cdr2, cdr3)
    hum_muts = '; '.join([f"{m['position']}:{m['current']}>{m['suggested']}" for m in hum_mm])
    
    # Graft into custom_hybrid
    cust_info = SCAFFOLDS['custom_hybrid']
    cust_fw = {k: v for k, v in cust_info.items() if k.startswith('FR')}
    _, cust_mm = score_framework(cust_fw, preds)
    cust_fw_mut = apply_mutations(cust_fw, cust_mm)
    cust_grafted = assemble(cust_fw_mut, cdr1, cdr2, cdr3)
    cust_muts = '; '.join([f"{m['position']}:{m['current']}>{m['suggested']}" for m in cust_mm])
    
    return Result(
        id=seq_id, original_seq=seq,
        cdr1=cdr1, cdr2=cdr2, cdr3=cdr3,
        fr1=regions['FR1'], fr2=regions['FR2'], fr3=regions['FR3'], fr4=regions['FR4'],
        valid=regions['valid'], family=family,
        orig_score=orig_score, orig_n_issues=len(orig_mm),
        best_scaffold=best_name, best_score=best_score, 
        best_n_mut=len(best_mm), best_mutations=best_muts, best_grafted=best_grafted,
        human_n_mut=len(hum_mm), human_mutations=hum_muts, human_grafted=hum_grafted,
        custom_n_mut=len(cust_mm), custom_mutations=cust_muts, custom_grafted=cust_grafted,
        error=''
    )

def read_seqs(path: str) -> List[Tuple[str, str]]:
    """Read sequences from file."""
    seqs = []
    with open(path) as f:
        content = f.read()
    
    if content.startswith('>'):
        for entry in content.split('>')[1:]:
            lines = entry.strip().split('\n')
            sid = lines[0].split()[0]
            seq = ''.join(lines[1:]).replace(' ', '').upper()
            if seq:
                seqs.append((sid, seq))
    elif ',' in content.split('\n')[0] or '\t' in content.split('\n')[0]:
        delim = ',' if ',' in content.split('\n')[0] else '\t'
        reader = csv.DictReader(content.strip().split('\n'), delimiter=delim)
        seq_col = None
        for col in reader.fieldnames or []:
            if col.lower() in ['sequence', 'seq', 'vhh', 'nanobody', 'aa_sequence']:
                seq_col = col
                break
        if seq_col:
            for i, row in enumerate(reader):
                seqs.append((f'seq_{i+1}', row[seq_col].upper()))
    else:
        for i, line in enumerate(content.strip().split('\n')):
            seq = line.strip().upper()
            if len(seq) > 80 and seq[0] in 'ACDEFGHIKLMNPQRSTVWY':
                seqs.append((f'seq_{i+1}', seq))
    
    return seqs

def write_csv(results: List[Result], path: str):
    """Write results."""
    fields = [
        'id', 'CDR1', 'CDR2', 'CDR3', 'CDR1_len', 'CDR2_len', 'CDR3_len',
        'predicted_family', 'original_score', 'original_n_issues',
        'best_scaffold', 'best_score', 'best_n_mutations', 'best_mutations', 'best_grafted',
        'humanized_n_mutations', 'humanized_mutations', 'humanized_grafted',
        'custom_n_mutations', 'custom_mutations', 'custom_grafted',
        'error'
    ]
    
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in results:
            w.writerow({
                'id': r.id,
                'CDR1': r.cdr1, 'CDR2': r.cdr2, 'CDR3': r.cdr3,
                'CDR1_len': len(r.cdr1), 'CDR2_len': len(r.cdr2), 'CDR3_len': len(r.cdr3),
                'predicted_family': r.family,
                'original_score': f'{r.orig_score:.1f}',
                'original_n_issues': r.orig_n_issues,
                'best_scaffold': r.best_scaffold,
                'best_score': f'{r.best_score:.1f}',
                'best_n_mutations': r.best_n_mut,
                'best_mutations': r.best_mutations,
                'best_grafted': r.best_grafted,
                'humanized_n_mutations': r.human_n_mut,
                'humanized_mutations': r.human_mutations,
                'humanized_grafted': r.human_grafted,
                'custom_n_mutations': r.custom_n_mut,
                'custom_mutations': r.custom_mutations,
                'custom_grafted': r.custom_grafted,
                'error': r.error
            })

def main():
    parser = argparse.ArgumentParser(description='VHH Framework Optimizer v2')
    parser.add_argument('--input', '-i', required=True)
    parser.add_argument('--output', '-o', required=True)
    parser.add_argument('--correlations', '-c', default='correlation_results_v3.pkl')
    parser.add_argument('--min-confidence', type=float, default=90.0,
                       help='Min confidence for mutations (default: 90)')
    
    args = parser.parse_args()
    
    print(f"Loading {args.correlations}...")
    with open(args.correlations, 'rb') as f:
        data = pickle.load(f)
    
    print(f"Reading {args.input}...")
    seqs = read_seqs(args.input)
    print(f"  Found {len(seqs)} sequences")
    
    print(f"Processing (min conf = {args.min_confidence}%)...")
    results = [process(sid, seq, data, args.min_confidence) for sid, seq in seqs]
    
    write_csv(results, args.output)
    
    # Summary
    valid = [r for r in results if r.valid]
    print(f"\n{'='*80}")
    print("SUMMARY")
    print(f"{'='*80}")
    print(f"Valid: {len(valid)}/{len(results)}")
    
    if valid:
        avg_score = sum(r.orig_score for r in valid) / len(valid)
        avg_issues = sum(r.orig_n_issues for r in valid) / len(valid)
        zero_issues = sum(1 for r in valid if r.orig_n_issues == 0)
        
        print(f"Average original score: {avg_score:.1f}%")
        print(f"Average issues (>90% conf): {avg_issues:.2f}")
        print(f"Sequences with 0 issues: {zero_issues}/{len(valid)} ({100*zero_issues/len(valid):.0f}%)")
        
        print(f"\nBest scaffold distribution:")
        for s, c in Counter(r.best_scaffold for r in valid).most_common():
            print(f"  {s}: {c}")
    
    print(f"\nOutput: {args.output}")

if __name__ == '__main__':
    main()