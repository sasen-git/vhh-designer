#!/usr/bin/env python3
"""
VHH Designer v2 - Mouse HC to VHH with Direct Mutations
========================================================

Converts a mouse heavy chain to VHH candidates using TWO approaches:
  1. UNIVERSAL: Graft CDRs into Universal VHH scaffold + apply mutations
  2. ORIGINAL: Keep original mouse FRs + apply VHH hallmark mutations

Data sources:
- correlation_results_v3_compensation.pkl: CDR->FR mutation rules (254 rules)
- epistasis_v2_full.pkl: Vernier cluster statistics for naturalness scoring

Usage:
  # Both approaches (default - 46 Universal + 46 Original)
  python vhh_designer_v2.py -s "EVQLVESGG..." -c correlations.pkl -e epistasis.pkl

  # Universal only
  python vhh_designer_v2.py -s "EVQLVESGG..." -c correlations.pkl -e epistasis.pkl --mode universal

  # Original only  
  python vhh_designer_v2.py -s "EVQLVESGG..." -c correlations.pkl -e epistasis.pkl --mode original
"""

import os
import re
import pickle
import argparse
import numpy as np
import pandas as pd
from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Optional, Set
from collections import defaultdict
from datetime import datetime

try:
    from antpack import SingleChainAnnotator
    ANTPACK_AVAILABLE = True
except ImportError:
    ANTPACK_AVAILABLE = False
    print("Warning: AntPack not available. Install: pip install antpack")

UNIVERSAL_SCAFFOLD = {
    'name': 'Universal Humanized VHH',
    'FR1': 'EVQLVESGGGLVQPGGSLRLSCAAS',
    'FR2': 'WFRQAPGKGLEWVS',
    'FR3': 'RFTISRDNSKNTLYLQMNSLRAEDTAVYYC',
    'FR4': 'WGQGTLVTVSS',
}

VHH_HALLMARKS = {
    'FERG': {'FR2[1]': 'F', 'FR2[8]': 'E', 'FR2[9]': 'R', 'FR2[11]': 'G'},
    'YERL': {'FR2[1]': 'Y', 'FR2[8]': 'E', 'FR2[9]': 'R', 'FR2[11]': 'L'},
    'FERF': {'FR2[1]': 'F', 'FR2[8]': 'E', 'FR2[9]': 'R', 'FR2[11]': 'F'},
    'FERA': {'FR2[1]': 'F', 'FR2[8]': 'E', 'FR2[9]': 'R', 'FR2[11]': 'A'},
}

AA_CATEGORIES = {
    'hydrophobic': set('AILMFVW'), 'polar': set('STNQ'),
    'basic': set('KRH'), 'acidic': set('DE'),
    'aromatic': set('FWY'), 'special': set('CGP'), 'reactive': set('CY'),
}

def categorize_aa(aa):
    for cat, aas in AA_CATEGORIES.items():
        if aa.upper() in aas:
            return cat
    return 'other'

@dataclass
class CDRSet:
    cdr1: str; cdr2: str; cdr3: str
    fr1: str = ''; fr2: str = ''; fr3: str = ''; fr4: str = ''

@dataclass
class Mutation:
    position: str; from_aa: str; to_aa: str
    rule: str; confidence: float; source: str
    def __str__(self): return f"{self.position}:{self.from_aa}->{self.to_aa}"

@dataclass 
class VHHCandidate:
    id: str; rank: int; sequence: str; framework_source: str
    cdr1: str; cdr2: str; cdr3: str
    fr1: str; fr2: str; fr3: str; fr4: str
    mutations: List[Mutation]
    naturalness_score: float = 0.0; confidence: float = 0.0
    cdr3_preserved: bool = True; strategy: str = ''
    warnings: List[str] = field(default_factory=list)

_annotator = None
def get_annotator():
    global _annotator
    if _annotator is None and ANTPACK_AVAILABLE:
        _annotator = SingleChainAnnotator(chains=['H'], scheme='imgt')
    return _annotator

def extract_cdrs(sequence):
    if not ANTPACK_AVAILABLE:
        print("ERROR: AntPack required"); return None
    annotator = get_annotator()
    seq = sequence.upper().replace('-', '').replace('.', '').replace(' ', '')
    try:
        numbering, pid, chain, err = annotator.analyze_seq(seq)
        # Only fail on actual errors, not warnings - check if numbering succeeded
        if numbering is None or len(numbering) == 0:
            print(f"  ANARCI failed: {err}")
            return None
        trimmed_seq, trimmed_nums, _, _ = annotator.trim_alignment(seq, (numbering, pid, chain, err))
        labels = annotator.assign_cdr_labels(trimmed_nums, chain)
        
        # Map AntPack labels to our expected format
        # AntPack returns: 'fmwk1', 'cdr1', 'fmwk2', etc.
        label_map = {
            # AntPack's actual labels (lowercase)
            'fmwk1': 'FR1', 'fmwk2': 'FR2', 'fmwk3': 'FR3', 'fmwk4': 'FR4',
            'cdr1': 'CDR1', 'cdr2': 'CDR2', 'cdr3': 'CDR3',
            # Other possible variations
            'fwh1': 'FR1', 'fwh2': 'FR2', 'fwh3': 'FR3', 'fwh4': 'FR4',
            'cdrh1': 'CDR1', 'cdrh2': 'CDR2', 'cdrh3': 'CDR3',
            'fwl1': 'FR1', 'fwl2': 'FR2', 'fwl3': 'FR3', 'fwl4': 'FR4',
            'cdrl1': 'CDR1', 'cdrl2': 'CDR2', 'cdrl3': 'CDR3',
            'fw1': 'FR1', 'fw2': 'FR2', 'fw3': 'FR3', 'fw4': 'FR4',
            'fr1': 'FR1', 'fr2': 'FR2', 'fr3': 'FR3', 'fr4': 'FR4',
            'FR1': 'FR1', 'FR2': 'FR2', 'FR3': 'FR3', 'FR4': 'FR4',
            'CDR1': 'CDR1', 'CDR2': 'CDR2', 'CDR3': 'CDR3',
        }
        
        # Initialize regions
        regions = {'FR1': '', 'CDR1': '', 'FR2': '', 'CDR2': '', 'FR3': '', 'CDR3': '', 'FR4': ''}
        for aa, label in zip(trimmed_seq, labels):
            # Normalize label
            normalized = label_map.get(label, label_map.get(label.lower(), None))
            if normalized and normalized in regions:
                regions[normalized] += aa
        
        # Debug: show what we got
        unique_labels = sorted(set(labels))
        if not regions['FR1'] or not regions['FR2']:
            print(f"  Warning: FRs not extracted. AntPack labels: {unique_labels}")
        
        return CDRSet(cdr1=regions['CDR1'], cdr2=regions['CDR2'], cdr3=regions['CDR3'],
                      fr1=regions['FR1'], fr2=regions['FR2'], fr3=regions['FR3'], fr4=regions['FR4'])
    except Exception as e:
        print(f"  CDR extraction error: {e}")
        return None

class CorrelationRulesEngine:
    def __init__(self, correlation_file):
        with open(correlation_file, 'rb') as f: self.data = pickle.load(f)
        self.multi_position_rules = self.data.get('multi_position_rules', [])
        self.vernier_archetypes = self.data.get('vernier_archetypes', {})
        print(f"Loaded {len(self.multi_position_rules)} rules")
    
    def get_cdr_features(self, cdrs):
        features = {}
        for i, aa in enumerate(cdrs.cdr1): features[f'cdr1[{i}]'] = aa
        for i in range(1, min(4, len(cdrs.cdr1)+1)): features[f'cdr1[{-i}]'] = cdrs.cdr1[-i] if i <= len(cdrs.cdr1) else ''
        for i, aa in enumerate(cdrs.cdr2): features[f'cdr2[{i}]'] = aa
        for i in range(1, min(4, len(cdrs.cdr2)+1)): features[f'cdr2[{-i}]'] = cdrs.cdr2[-i] if i <= len(cdrs.cdr2) else ''
        for i, aa in enumerate(cdrs.cdr3): features[f'cdr3[{i}]'] = aa
        for i in range(1, min(4, len(cdrs.cdr3)+1)): features[f'cdr3[{-i}]'] = cdrs.cdr3[-i] if i <= len(cdrs.cdr3) else ''
        return features
    
    def find_applicable_rules(self, cdrs, min_confidence=75.0):
        features = self.get_cdr_features(cdrs)
        applicable = []
        for rule in self.multi_position_rules:
            for pred in rule.get('top_predictors', []):
                if pred['conf'] < min_confidence: continue
                match = re.match(r'(cdr\d+)\[(-?\d+)\]=([A-Z])', pred.get('rule', ''))
                if match:
                    cdr, pos, expected = match.groups()
                    if features.get(f'{cdr}[{int(pos)}]') == expected:
                        applicable.append({'fw_position': rule['fw_position'], 'suggested_aa': rule['fw_residue'],
                                          'confidence': pred['conf'], 'rule': pred['rule'],
                                          'family': rule.get('family', ''), 'source': 'correlation'})
        by_pos = {}
        for r in applicable:
            if r['fw_position'] not in by_pos or r['confidence'] > by_pos[r['fw_position']]['confidence']:
                by_pos[r['fw_position']] = r
        return sorted(by_pos.values(), key=lambda x: -x['confidence'])
    
    def get_vernier_mutations(self, family='F_C2'):
        return [{'fw_position': pos, 'suggested_aa': info.get('dominant_residue', ''),
                 'confidence': 80.0, 'rule': f'Vernier {family}', 'source': 'vernier'}
                for pos, info in self.vernier_archetypes.get(family, {}).items()]

class NaturalnessScorer:
    def __init__(self, epistasis_file):
        with open(epistasis_file, 'rb') as f: self.data = pickle.load(f)
        self.clusters = self.data.get('analysis_2_vernier_clusters', {})
        lengths, charges = [], []
        for cluster in self.clusters.values():
            if isinstance(cluster, dict):
                cdr3_stats = cluster.get('cdr3_length', {})
                n = cdr3_stats.get('n', 0) if isinstance(cdr3_stats, dict) else 0
                if n >= 100:
                    lengths.extend([cdr3_stats.get('mean', 15)] * min(n, 500))
                    charge_stats = cluster.get('cdr3_charge', {})
                    if isinstance(charge_stats, dict):
                        charges.extend([charge_stats.get('mean', 0)] * min(n, 500))
        self.len_mean = np.mean(lengths) if lengths else 15
        self.len_std = np.std(lengths) if lengths else 5
        self.chg_mean = np.mean(charges) if charges else 0
        self.chg_std = np.std(charges) if charges else 3
    
    def score(self, cdr3):
        z_len = abs(len(cdr3) - self.len_mean) / max(self.len_std, 1)
        chg = sum(1 for a in cdr3 if a in 'KRH') - sum(1 for a in cdr3 if a in 'DE')
        z_chg = abs(chg - self.chg_mean) / max(self.chg_std, 0.5)
        return round(max(0, 100 - (z_len + z_chg) / 2 * 25), 1)

class VHHDesigner:
    def __init__(self, correlation_file, epistasis_file):
        self.rules = CorrelationRulesEngine(correlation_file)
        self.scorer = NaturalnessScorer(epistasis_file)
        self.scaffold = UNIVERSAL_SCAFFOLD.copy()
    
    def design(self, cdrs, n_candidates=92, mode='both', input_sequence=None):
        all_rules = self.rules.find_applicable_rules(cdrs, 70.0)
        print(f"\nFound {len(all_rules)} applicable rules")
        vernier_F = self.rules.get_vernier_mutations('F_C2')
        vernier_Y = self.rules.get_vernier_mutations('Y_C2')
        
        if mode == 'both': n_univ, n_orig = n_candidates // 2, n_candidates - n_candidates // 2
        elif mode == 'universal': n_univ, n_orig = n_candidates, 0
        else: n_univ, n_orig = 0, n_candidates
        
        candidates = []
        
        # Add original input sequence as first lead (for all modes)
        if input_sequence:
            orig_lead = VHHCandidate(
                id="Input_Original",
                rank=0,
                sequence=input_sequence,
                framework_source='input',
                cdr1=cdrs.cdr1, cdr2=cdrs.cdr2, cdr3=cdrs.cdr3,
                fr1=cdrs.fr1, fr2=cdrs.fr2, fr3=cdrs.fr3, fr4=cdrs.fr4,
                mutations=[],
                strategy="original_input"
            )
            candidates.append(orig_lead)
        
        if n_univ > 0:
            print(f"Generating {n_univ} Universal candidates...")
            candidates.extend(self._design_universal(cdrs, all_rules, vernier_F, vernier_Y, n_univ))
        if n_orig > 0 and cdrs.fr1 and cdrs.fr2:
            print(f"Generating {n_orig} Original candidates...")
            candidates.extend(self._design_original(cdrs, all_rules, vernier_F, vernier_Y, n_orig))
        elif n_orig > 0:
            print("No original FRs, using Universal")
            candidates.extend(self._design_universal(cdrs, all_rules, vernier_F, vernier_Y, n_orig))
        
        for c in candidates:
            c.naturalness_score = self.scorer.score(c.cdr3)
            if c.id == "Input_Original":
                c.confidence = 100.0  # Original always ranked first
            elif c.id == "Univ_Graft":
                c.confidence = 99.0  # Pure graft always ranked second
            elif c.mutations:
                c.confidence = 0.4 * c.naturalness_score + 0.6 * np.mean([m.confidence for m in c.mutations])
            else:
                c.confidence = c.naturalness_score * 0.85
        candidates.sort(key=lambda x: -x.confidence)
        for i, c in enumerate(candidates): c.rank = i + 1
        return candidates[:n_candidates]
    
    def _design_universal(self, cdrs, all_rules, vernier_F, vernier_Y, n):
        candidates, used = [], set()
        # Pure graft: CDRs into Universal scaffold with no additional mutations
        candidates.append(self._make("Univ_Graft", cdrs, [], "universal_graft", 'universal'))
        
        for i, r in enumerate(all_rules[:min(10, n//4)]):
            m = self._rule2mut(r, self.scaffold)
            if m: candidates.append(self._make(f"Univ_Rule_{i+1}", cdrs, [m], "single_rule", 'universal'))
        
        f_muts = [self._rule2mut(r, self.scaffold) for r in vernier_F]; f_muts = [m for m in f_muts if m]
        if f_muts: candidates.append(self._make("Univ_Vernier_F_C2", cdrs, f_muts, "vernier_F_C2", 'universal'))
        y_muts = [self._rule2mut(r, self.scaffold) for r in vernier_Y]; y_muts = [m for m in y_muts if m]
        if y_muts: candidates.append(self._make("Univ_Vernier_Y_C2", cdrs, y_muts, "vernier_Y_C2", 'universal'))
        
        high = [r for r in all_rules if r['confidence'] >= 85]
        for i in range(min(len(high), 8)):
            for j in range(i+1, min(len(high), 8)):
                if high[i]['fw_position'] != high[j]['fw_position']:
                    m1, m2 = self._rule2mut(high[i], self.scaffold), self._rule2mut(high[j], self.scaffold)
                    if m1 and m2:
                        key = tuple(sorted([str(m1), str(m2)]))
                        if key not in used:
                            used.add(key); candidates.append(self._make(f"Univ_Pair_{len(candidates)}", cdrs, [m1, m2], "pair", 'universal'))
                if len(candidates) >= n: break
            if len(candidates) >= n: break
        
        import random
        while len(candidates) < n and len(all_rules) >= 2:
            picked = random.sample(all_rules, min(3, len(all_rules)))
            if len({r['fw_position'] for r in picked}) == len(picked):
                muts = [self._rule2mut(r, self.scaffold) for r in picked]; muts = [m for m in muts if m]
                if muts:
                    key = tuple(sorted([str(m) for m in muts]))
                    if key not in used: used.add(key); candidates.append(self._make(f"Univ_Div_{len(candidates)}", cdrs, muts, "diverse", 'universal'))
        return candidates[:n]
    
    def _design_original(self, cdrs, all_rules, vernier_F, vernier_Y, n):
        candidates, used = [], set()
        orig = {'FR1': cdrs.fr1, 'FR2': cdrs.fr2, 'FR3': cdrs.fr3, 'FR4': cdrs.fr4}
        
        # Note: Input_Original is added in design() as the lead
        # Start with hallmark variants (VHH-izing mutations)
        for hname, hpos in VHH_HALLMARKS.items():
            muts = [self._pos2mut(p, a, orig, f"hallmark_{hname}", 95.0) for p, a in hpos.items()]
            muts = [m for m in muts if m]
            if muts: candidates.append(self._make(f"Orig_{hname}", cdrs, muts, f"hallmark_{hname}", 'original', orig))
        
        base_muts = [self._pos2mut(p, a, orig, "FERG", 95.0) for p, a in VHH_HALLMARKS['FERG'].items()]
        base_muts = [m for m in base_muts if m]
        
        for r in vernier_F:
            if r['fw_position'] not in [m.position for m in base_muts]:
                m = self._rule2mut(r, orig)
                if m: candidates.append(self._make(f"Orig_FERG_V_{r['fw_position']}", cdrs, base_muts + [m], "hallmark+vernier", 'original', orig))
        
        for i, r in enumerate(all_rules[:min(10, n//3)]):
            if r['fw_position'] in ['FR2[1]', 'FR2[8]', 'FR2[9]', 'FR2[11]']: continue
            m = self._rule2mut(r, orig)
            if m: candidates.append(self._make(f"Orig_FERG_R{i+1}", cdrs, base_muts + [m], "hallmark+rule", 'original', orig))
        
        high = [r for r in all_rules if r['confidence'] >= 85 and r['fw_position'] not in ['FR2[1]', 'FR2[8]', 'FR2[9]', 'FR2[11]']]
        for i in range(min(len(high), 6)):
            for j in range(i+1, min(len(high), 6)):
                if high[i]['fw_position'] != high[j]['fw_position']:
                    m1, m2 = self._rule2mut(high[i], orig), self._rule2mut(high[j], orig)
                    if m1 and m2:
                        key = tuple(sorted([str(m) for m in base_muts + [m1, m2]]))
                        if key not in used:
                            used.add(key); candidates.append(self._make(f"Orig_FERG_P{len(candidates)}", cdrs, base_muts + [m1, m2], "hallmark+pair", 'original', orig))
                if len(candidates) >= n: break
            if len(candidates) >= n: break
        
        import random
        while len(candidates) < n:
            hname = random.choice(['FERG', 'YERL', 'FERF'])
            hmuts = [self._pos2mut(p, a, orig, hname, 95.0) for p, a in VHH_HALLMARKS[hname].items()]
            hmuts = [m for m in hmuts if m]
            avail = [r for r in all_rules if r['fw_position'] not in ['FR2[1]', 'FR2[8]', 'FR2[9]', 'FR2[11]']]
            if avail:
                picked = random.sample(avail, min(2, len(avail)))
                extra = [self._rule2mut(r, orig) for r in picked]; extra = [m for m in extra if m]
                if extra:
                    all_m = hmuts + extra
                    key = tuple(sorted([str(m) for m in all_m]))
                    if key not in used: used.add(key); candidates.append(self._make(f"Orig_Div_{len(candidates)}", cdrs, all_m, f"diverse_{hname}", 'original', orig))
            if len(used) > n * 2: break
        return candidates[:n]
    
    def _rule2mut(self, rule, frs):
        match = re.match(r'(FR\d)\[(\d+)\]', rule['fw_position'])
        if not match: return None
        region, idx = match.groups(); idx = int(idx)
        fr_seq = frs.get(region, '')
        if idx >= len(fr_seq): return None
        cur = fr_seq[idx]
        if cur == rule['suggested_aa']: return None
        return Mutation(rule['fw_position'], cur, rule['suggested_aa'], rule.get('rule', ''), rule.get('confidence', 75), rule.get('source', 'correlation'))
    
    def _pos2mut(self, pos, to_aa, frs, rule, conf):
        match = re.match(r'(FR\d)\[(\d+)\]', pos)
        if not match: return None
        region, idx = match.groups(); idx = int(idx)
        fr_seq = frs.get(region, '')
        if idx >= len(fr_seq): return None
        cur = fr_seq[idx]
        if cur == to_aa: return None
        return Mutation(pos, cur, to_aa, rule, conf, 'hallmark')
    
    def _make(self, id, cdrs, mutations, strategy, fw_source, custom_frs=None):
        if fw_source == 'universal' or not custom_frs:
            fr1, fr2, fr3, fr4 = self.scaffold['FR1'], self.scaffold['FR2'], self.scaffold['FR3'], self.scaffold['FR4']
        else:
            fr1, fr2, fr3, fr4 = custom_frs['FR1'], custom_frs['FR2'], custom_frs['FR3'], custom_frs['FR4']
        
        frs = {'FR1': list(fr1), 'FR2': list(fr2), 'FR3': list(fr3), 'FR4': list(fr4)}
        for m in mutations:
            match = re.match(r'(FR\d)\[(\d+)\]', m.position)
            if match:
                region, idx = match.groups(); idx = int(idx)
                if region in frs and idx < len(frs[region]): frs[region][idx] = m.to_aa
        
        fr1, fr2, fr3, fr4 = ''.join(frs['FR1']), ''.join(frs['FR2']), ''.join(frs['FR3']), ''.join(frs['FR4'])
        seq = fr1 + cdrs.cdr1 + fr2 + cdrs.cdr2 + fr3 + cdrs.cdr3 + fr4
        return VHHCandidate(id, 0, seq, fw_source, cdrs.cdr1, cdrs.cdr2, cdrs.cdr3, fr1, fr2, fr3, fr4, mutations, strategy=strategy)

def format_mutations(mutations):
    if not mutations: return "None"
    return "; ".join(f"{m.position}:{m.from_aa}->{m.to_aa} [{m.source}, {m.confidence:.0f}%]" for m in mutations)

def to_dataframe(candidates, input_seq, cdrs):
    return pd.DataFrame([{
        'rank': c.rank, 'id': c.id, 'sequence': c.sequence, 'length': len(c.sequence),
        'framework_source': c.framework_source,
        'CDR1': c.cdr1, 'CDR2': c.cdr2, 'CDR3': c.cdr3, 'CDR3_preserved': c.cdr3 == cdrs.cdr3,
        'FR1': c.fr1, 'FR2': c.fr2, 'FR3': c.fr3, 'FR4': c.fr4,
        'confidence': round(c.confidence, 1), 'naturalness_score': round(c.naturalness_score, 1),
        'n_mutations': len(c.mutations), 'mutations': format_mutations(c.mutations),
        'rules_applied': "; ".join(m.rule for m in c.mutations) if c.mutations else "None",
        'strategy': c.strategy, 'input_sequence': input_seq,
    } for c in candidates])

def save_results(df, output_dir, name, mode, datetime_str):
    """Save results to folder with xlsx and MSA-ready csv."""
    import os
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Generate filenames with datetime prefix
    base_name = f"{datetime_str}_{name}_{mode}"
    
    xlsx_path = os.path.join(output_dir, f"{base_name}.xlsx")
    csv_path = os.path.join(output_dir, f"{base_name}_sequences.csv")
    
    # Save Excel with all sheets
    with pd.ExcelWriter(xlsx_path, engine='openpyxl') as w:
        df.to_excel(w, sheet_name='Candidates', index=False)
        pd.DataFrame({'Metric': ['Total', 'Input (original)', 'Universal', 'Original', 'CDR3 preserved', 'Avg conf'],
            'Value': [len(df), (df['framework_source']=='input').sum(), 
                      (df['framework_source']=='universal').sum(), (df['framework_source']=='original').sum(),
                      f"{df['CDR3_preserved'].sum()} ({100*df['CDR3_preserved'].mean():.0f}%)", f"{df['confidence'].mean():.1f}"]
        }).to_excel(w, sheet_name='Summary', index=False)
    
    # Save simple CSV for MSA (lead, ID, Seq)
    msa_df = pd.DataFrame()
    msa_df['lead'] = df['id'].apply(lambda x: 'true' if x == 'Input_Original' else 'false')
    msa_df['ID'] = df['id']
    msa_df['Seq'] = df['sequence']
    msa_df.to_csv(csv_path, index=False)
    
    print(f"\nOutput folder: {output_dir}/")
    print(f"  {base_name}.xlsx (full results)")
    print(f"  {base_name}_sequences.csv (for MSA)")
    
    return output_dir

def main():
    p = argparse.ArgumentParser(description='VHH Designer v2')
    p.add_argument('--sequence', '-s', help='Input sequence directly')
    p.add_argument('--input', '-i', help='Input FASTA file')
    p.add_argument('--name', help='Sequence name (auto-detected from FASTA header if not provided)')
    p.add_argument('--correlations', '-c', required=True, help='Correlation rules pickle')
    p.add_argument('--epistasis', '-e', required=True, help='Epistasis pickle')
    p.add_argument('--output-dir', '-o', help='Output directory (default: auto-generated)')
    p.add_argument('--n-candidates', '-n', type=int, default=92, help='Number of candidates')
    p.add_argument('--mode', '-m', choices=['both', 'universal', 'original'], default='both',
                   help='Design mode: both, universal, or original')
    args = p.parse_args()
    
    # Get input sequence and name
    seq_name = args.name
    if args.sequence:
        input_seq = args.sequence.upper().replace(' ', '').replace('\n', '')
        if not seq_name:
            seq_name = "input"
    elif args.input:
        with open(args.input) as f:
            content = f.read()
        if args.input.endswith(('.fasta', '.fa')):
            lines = content.strip().split('\n')
            # Extract name from FASTA header
            header = [l for l in lines if l.startswith('>')]
            if header and not seq_name:
                # Get name from header, clean it up
                seq_name = header[0][1:].split()[0].replace('/', '_').replace('\\', '_')
            input_seq = ''.join(l for l in lines if not l.startswith('>')).upper()
        else:
            input_seq = content.strip().split('\n')[0].upper()
        if not seq_name:
            # Use filename without extension
            seq_name = args.input.rsplit('/', 1)[-1].rsplit('.', 1)[0]
    else:
        p.error("Need --sequence or --input")
        return
    
    # Clean up sequence name for filesystem
    seq_name = re.sub(r'[^\w\-]', '_', seq_name)[:50]  # Limit length, safe chars only
    
    print("=" * 70)
    print("VHH DESIGNER v2")
    print("=" * 70)
    print(f"Sequence: {seq_name}")
    print(f"Mode: {args.mode.upper()}")
    print(f"Input ({len(input_seq)} aa): {input_seq[:50]}...")
    
    cdrs = extract_cdrs(input_seq)
    if not cdrs:
        print("ERROR: CDR extraction failed")
        return
    
    print(f"CDR1: {cdrs.cdr1}")
    print(f"CDR2: {cdrs.cdr2}")
    print(f"CDR3: {cdrs.cdr3}")
    print(f"FR1:  {cdrs.fr1[:20]}... ({len(cdrs.fr1)} aa)" if cdrs.fr1 else "FR1:  (empty)")
    print(f"FR2:  {cdrs.fr2} ({len(cdrs.fr2)} aa)" if cdrs.fr2 else "FR2:  (empty)")
    print(f"FR3:  {cdrs.fr3[:20]}... ({len(cdrs.fr3)} aa)" if cdrs.fr3 else "FR3:  (empty)")
    print(f"FR4:  {cdrs.fr4} ({len(cdrs.fr4)} aa)" if cdrs.fr4 else "FR4:  (empty)")
    
    if args.mode in ['original', 'both'] and not (cdrs.fr1 and cdrs.fr2 and cdrs.fr3 and cdrs.fr4):
        print("\nWARNING: Original mode requested but FRs not extracted!")
        print("         This usually means CDR extraction partially failed.")
        print("         Falling back to Universal mode for 'original' candidates.")
    
    designer = VHHDesigner(args.correlations, args.epistasis)
    candidates = designer.design(cdrs, args.n_candidates, args.mode, input_sequence=input_seq)
    df = to_dataframe(candidates, input_seq, cdrs)
    
    print(f"\nTotal: {len(df)} | Input: {(df['framework_source']=='input').sum()} | Univ: {(df['framework_source']=='universal').sum()} | Orig: {(df['framework_source']=='original').sum()}")
    print("Top 10:")
    for _, r in df.head(10).iterrows():
        print(f"  {r['rank']:2d}. {r['id'][:40]:40s} [{r['framework_source'][:4]}] conf={r['confidence']:.1f}")
    
    # Create output directory in results/analysis_runs/
    datetime_str = datetime.now().strftime('%Y%m%d_%H%M%S')
    script_name = "vhh_designer_v2"
    folder_name = f"{datetime_str}_{script_name}_{seq_name}_{args.mode}"
    
    if args.output_dir:
        output_dir = args.output_dir
    else:
        # Default to results/analysis_runs/ relative to script location or current dir
        results_base = "results/analysis_runs"
        if os.path.exists(results_base):
            output_dir = os.path.join(results_base, folder_name)
        else:
            output_dir = folder_name
    
    save_results(df, output_dir, f"{script_name}_{seq_name}", args.mode, datetime_str)
    print("\nDone!")

if __name__ == '__main__': main()