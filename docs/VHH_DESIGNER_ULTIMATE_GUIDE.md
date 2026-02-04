> **VHH_DESIGNER_ULTIMATE_GUIDE — Lightweight pointer**
>
> This document describes the technical architecture of VHH Designer.** For design philosophy, safeguards, and scope, see the project README**.

---

# VHH Designer v9.0 - Complete Technical Documentation

**Design & Architecture:** Sasen Efrem  
**Implementation & Documentation:** Claude (Anthropic)

## Table of Contents
1. [Executive Summary](#1-executive-summary)
2. [Background: The Biology](#2-background-the-biology)
3. [The Upstream Data Pipeline](#3-the-upstream-data-pipeline)
4. [Analysis: Learning Patterns from Nature](#4-analysis-learning-patterns-from-nature)
5. [VHH Designer Architecture](#5-vhh-designer-architecture)
6. [The M69 Example: A Real Case Study](#6-the-m69-example-a-real-case-study)
7. [The Track System](#7-the-track-system)
8. [Selection Pipeline](#8-selection-pipeline)
9. [Scoring System](#9-scoring-system)
10. [Output Files](#10-output-files)
11. [Command Line Reference](#11-command-line-reference)
12. [Caveats & Limitations](#12-caveats--limitations)
13. [Appendices](#appendices)

---

## 1. Executive Summary

**VHH Designer v9.0** converts conventional antibodies into nanobodies (single-domain antibodies) by intelligently mutating framework residues based on statistical patterns learned from **~12 million natural camelid sequences**.

This pipeline integrates large-scale biological data analysis, statistical pattern extraction, and machine learning validation (ESM2) into a unified protein engineering workflow—bridging computational biology with practical therapeutic design considerations.

### What Goes In
A conventional antibody sequence (VH domain) - mouse, human, or existing VHH.

### What Comes Out
**N nanobody variants** (user-specified via `--n-select`) organized as:
- **1 Lead**: Your original sequence (unchanged reference)
- **Controls**: Track 0 reference mutations (ranking-exempt, typically ~30-35)
- **Ranked**: Optimized candidates from Tracks 1-4 (scored and ranked, fills remaining slots)

### The Core Problem

You have a mouse or human antibody with valuable CDRs (binding loops). You want to convert its framework to a VHH-like scaffold while:
- Preserving the CDR sequences (your binding specificity)
- Making it stable without a light chain
- Making it look "natural" to a language model (ESM2)
- Following patterns observed in real VHHs

### What Makes This Different

Most VHH conversion tools apply generic rules. VHH Designer uses **hallmark-specific vernier consensus** - the framework mutations are tailored to the specific hallmark tetrad you're targeting, based on what's actually observed in nature for that hallmark.

Modern ML tools like AlphaFold, Boltz, and ESMFold have revolutionized protein science by accurately predicting 3D structures, while inverse folding methods like ProteinMPNN design sequences for given structures. These are powerful general-purpose tools. VHH Designer takes a different approach: rather than learning general protein patterns, it learns *VHH-specific* mutation patterns from 12 million real camelid sequences. It knows which framework changes actually accompany successful VH→VHH conversion in nature—knowledge that general-purpose models don't capture. The result is a domain-specific design tool that proposes biologically-grounded mutations, which can then be validated using structure predictors and scored using language models like ESM2.

### Design Philosophy

The pipeline emphasizes **interpretability and traceability**: every mutation can be traced back to its statistical support in the natural sequence database, enabling researchers to make informed decisions about which variants to pursue experimentally. Multiple validation layers—structural constraints, cysteine pairing checks, and ML-based naturalness scoring—help ensure that generated sequences are not just statistically plausible but biologically viable.

---

## 2. Background: The Biology

### 2.1 The Problem We're Solving

Human antibodies have a **heavy chain** (VH) that pairs with a **light chain** (VL). They need each other to be stable. Camelid antibodies evolved a special type that works alone - called **VHH** or **nanobodies**.

The key difference is in **four amino acid positions** in the framework. These positions normally form the VH-VL interface. In nanobodies, they've evolved different amino acids that make the antibody soluble and stable without a partner.

### 2.2 The Four Hallmark Positions

The "hallmark" positions at IMGT **42, 49, 50, and 52** define the VH→VHH transformation:

| Position | Human VH | Camelid VHH | Why It Matters |
|----------|----------|-------------|----------------|
| IMGT 42 | V (small) | F or Y (bulky aromatic) | Fills the gap where VL would sit |
| IMGT 49 | G (tiny) | E, Q, or K (charged/polar) | Adds solubility, fills cavity |
| IMGT 50 | L (hydrophobic) | R (positive charge) | **Critical** - provides solubility |
| IMGT 52 | W (aromatic) | G, A, L, or F (various) | VL contact removed |

**Position 50 is the most important.** The change from L (leucine, greasy) to R (arginine, charged) is what allows the nanobody to be soluble without a light chain partner.

### 2.3 One Scaffold, Many Families

VHH families are a bit like dog breeds. All dogs share the same basic body plan—four legs, a tail, ears—just as every VHH shares conserved structural elements (framework regions, CDRs, disulfide bonds). But breeds have distinct characteristics that make a Husky fundamentally different from a Labrador, and VHH families have signature patterns that distinguish them just as clearly:

| Dog Breed | Distinguishing Features | VHH Family Equivalent |
|-----------|------------------------|----------------------|
| **Labrador** | Short coat, floppy ears | **F_C2** (Phe-42, 2 cysteines) |
| **Poodle** | Curly coat, pointed muzzle | **Y_C2** (Tyr-42, 2 cysteines) |
| **Husky** | Blue eyes, thick coat | **F_C4** (Phe-42, 4 cysteines) |

You wouldn't breed a Chihuahua using Great Dane genetics—and we don't apply F_C2 mutations to a Y_C4 scaffold. The statistical patterns that make one family stable and functional are specific to that family's structural context.

### 2.4 The Hallmark Tetrad Notation

We represent the four hallmark positions as a **4-letter code**:
- **FERG** = F at 42, E at 49, R at 50, G at 52 (classic nanobody)
- **VGLW** = V at 42, G at 49, L at 50, W at 52 (human VH-like)
- **YQRL** = Y at 42, Q at 49, R at 50, L at 52 (specialized nanobody)

### 2.5 Literature Foundation

| Year | Authors | Key Finding |
|------|---------|-------------|
| 1993 | Hamers-Casterman et al. | Discovery of heavy-chain-only antibodies in camelids |
| 1994 | Muyldermans et al. | Identified hallmark positions 42, 49, 50, 52 |
| 1992 | Foote & Winter | Defined the "vernier zone" - framework positions supporting CDRs |
| 2009 | Vincke et al. | Comprehensive VHH sequence analysis, subfamily variation |
| 2012 | Saerens et al. | CDR3 disulfide patterns (C4 architecture) |

---

## 3. The Upstream Data Pipeline

### 3.1 Data Sources Overview

The consolidated database (`VHH_db_final.npz`) contains **12,049,840 unique sequences** from multiple sources:

| Source | Sequences | Description | Has Targets? |
|--------|-----------|-------------|--------------|
| **INDI_NGS** | 10,675,894 | NGS repertoires from llama, alpaca, camel | ❌ |
| **OAS_Camel** | 1,354,537 | OAS bulk download (camel heavy chains) | ❌ |
| **INDI_patent** | 14,356 | Patent-derived sequences | ✅ Biomolecule targets |
| **SAbDab_PDB** | 1,559 | Crystal structures | ✅ Antigen info |
| **INDI_manual** | 1,180 | Manually curated entries | ✅ Some |
| **INDI_genbank** | 1,128 | GenBank submissions | ❌ |
| **GenBank** | 433 | From earlier merge | Some |
| **INDI_structure** | 429 | PDB structures | ✅ PDB info |
| **SAbDab** | 309 | Structural database | ✅ |
| **TheraSAbDab** | 15 | Therapeutic nanobodies | ✅ |

**Key insight:** ~99% of sequences (INDI_NGS + OAS_Camel) are repertoire data without antigen info. The ~19,000 annotated sequences (patents, structures, curated) provide target/binding information.

### 3.2 OAS Camel Data: Nucleotide to Amino Acid Translation

The OAS_Camel subset required special processing because the pre-aligned amino acid column was **truncated**.

**The Problem:**
- OAS provides `sequence_alignment_aa` column pre-aligned to 196 positions
- This truncated sequences to fit a fixed-width canonical alignment
- V-domains appeared as 14-75 aa instead of the correct 100-130 aa
- CDR3 regions were incomplete or missing

**Raw Source Files (SRR = NCBI Sequence Read Archive):**
```
SRR3544217_Heavy_Bulk.csv  (307,632 sequences)
SRR3544218_Heavy_Bulk.csv  (280,073 sequences)
SRR3544219_Heavy_Bulk.csv  (235,315 sequences)
SRR3544220_Heavy_Bulk.csv  (235,620 sequences)
SRR3544221_Heavy_Bulk.csv  (273,732 sequences)
SRR3544222_Heavy_Bulk.csv  (268,598 sequences)
SRR3544217_Heavy_IGHM.csv  (130 sequences)
Total: ~1.6M raw nucleotide sequences
```

**The Solution Pipeline:**
```
Raw nucleotide sequences (300-450+ nt)
    ↓
Biopython translation (NT → AA)
    ↓  
ANARCI numbering (IMGT scheme)
    ↓
CDR extraction (IMGT boundaries)
    ↓
camel_vhh_clean_db.npz (~1.36M sequences)
```

### 3.3 INDI Database

INDI (Integrated Nanobody Database for Immunoinformatics) from Natural Antibody provided:

| File | Sequences | Metadata |
|------|-----------|----------|
| `patent_sequence.tsv` | 14,376 | ✅ `patent_meta.tsv` has biomolecule targets, patent IDs |
| `structure_sequence.tsv` | 535 | ✅ `structure_meta.tsv` has PDB info, organism |
| `manual_sequence.tsv` | 1,268 | ✅ `manual_meta.tsv` has URLs, titles |
| `abgenbank_sequence.tsv` | 1,858 | ❌ No metadata file |
| `ngs_sequence.tsv` | 10.7M | ❌ NGS repertoire (PRJNA/SRR IDs only) |

**INDI NGS species coverage:** Includes llama, alpaca, and camel sequences from multiple NGS studies.

### 3.4 Database Consolidation Pipeline

```
                OAS Bulk Download               INDI Database
                (Nucleotides)                   (Multiple sources)
                      │                               │
                      ▼                               ▼
                Translation                     Already AA
                (Biopython)                          │
                      │                               │
                      ▼                               ▼
                    ANARCI                      ANARCI
                    (IMGT numbering)           (IMGT numbering)
                      │                               │
                      ▼                               ▼
                camel_vhh_clean_db.npz         INDI sources
                (1.36M sequences)              (10.7M + 19K)
                      │                               │
                      └───────────┬───────────────────┘
                                  ▼
                          Deduplication
                          (keep best metadata)
                                  │
                                  ▼
                          VHH_db_final.npz
                          (12,049,840 sequences)
                                  │
                                  ▼
                          Sharding by source
                                  │
                    ┌─────────────┼─────────────────────┐
                    ▼             ▼                     ▼
            vhh_annotated.npz   vhh_oas_camel.npz   vhh_indi_ngs_*.npz
            (19,409 seqs)       (1.35M seqs)        (10.7M in 6 shards)
            Patents, PDB,       OAS Camel           INDI NGS
            curated             repertoire          repertoire
```

### 3.5 IMGT Numbering

The **IMGT numbering scheme** (Lefranc, 1997) provides a universal coordinate system:

```
Position ranges:
  FR1:  1-26
  CDR1: 27-38
  FR2:  39-55
  CDR2: 56-65
  FR3:  66-104
  CDR3: 105-117 (variable length, uses insertions like 111A, 111B)
  FR4:  118-128
```

**Critical VHH positions:**
- IMGT 42, 49, 50, 52: Hallmark positions
- IMGT 23, 104: Canonical cysteines (core disulfide bond)
- IMGT 55, 100: Extra cysteines (C4 families only)
- IMGT 66-94: Vernier zone (17 key positions)

### 3.6 CDR Extraction

**All sources were processed through ANARCI with IMGT numbering:**

Regardless of whether CDRs were pre-extracted in source files, every sequence was re-processed through ANARCI to ensure consistent IMGT numbering across the entire 12M database.

**IMGT CDR boundaries:**
```
CDR1: IMGT positions 27-38
CDR2: IMGT positions 56-65  
CDR3: IMGT positions 105-117 (with insertions like 111A, 111B for longer loops)
```

---

## 4. Analysis: Learning Patterns from Nature

The analysis phase represents the core data science work underlying VHH Designer: transforming 12 million raw sequences into actionable design rules. This required building custom pipelines in Python for large-scale sequence processing, statistical pattern extraction, and cross-validation of discovered rules against held-out data.

### 4.1 Family Classification

VHH sequences are classified based on hallmark positions and cysteine count:

```python
def classify_family(positions, cdr3, full_seq):
    pos42 = positions['IMGT42']
    pos49 = positions['IMGT49']
    pos50 = positions['IMGT50']
    pos52 = positions['IMGT52']
    n_cys = full_seq.count('C')
    
    # VH-like: has L at position 50 (human-like)
    if pos50 == 'L':
        return 'VH_like'
    
    # Classical VHH criteria
    is_classical = (pos42 in ['F','Y'] and 
                    pos49 in ['E','Q'] and 
                    pos50 == 'R')
    
    if is_classical:
        if pos42 == 'Y' and n_cys == 2: return 'Y_C2'
        if pos42 == 'F' and n_cys == 2: return 'F_C2'
        if pos42 == 'F' and n_cys == 4: return 'F_C4'
        if pos42 == 'Y' and n_cys == 4: return 'Y_C4'
        return 'Classical_other'
    
    return 'Non_classical'
```

**Family meanings:**
- **Y_C2**: Tyrosine at 42, 2 cysteines (standard disulfide only)
- **F_C2**: Phenylalanine at 42, 2 cysteines (most common)
- **Y_C4**: Tyrosine at 42, 4 cysteines (extra CDR3 disulfide)
- **F_C4**: Phenylalanine at 42, 4 cysteines (extra disulfide)
- **VH_like**: Has leucine at position 50 (humanized VHHs)

### 4.2 Hallmark Discovery: 128 Patterns Found

Analysis of 12M sequences revealed **128 distinct hallmark patterns** at positions 42/49/50/52:

**Breakdown by Position 50 (the critical VHH marker):**

| Position 50 | Count | Interpretation |
|-------------|-------|----------------|
| R (Arginine) | 73 | True VHH - can work without light chain |
| L (Leucine) | 38 | VH-like - designed to pair with VL |
| Other | 17 | Rare variants |

**Top 10 True VHH Hallmarks (R at position 50):**

| Rank | Hallmark | Sequences | Family | Notes |
|------|----------|-----------|--------|-------|
| 1 | FERG | 3,664,062 | F_C2 | Most common |
| 2 | YQRL | 1,105,608 | Y_C2 | Most constrained |
| 3 | FERF | 1,085,360 | F_C2 | F at pos52 |
| 4 | FERA | 362,407 | F_C2 | A at pos52 |
| 5 | YERL | 279,551 | Y_C2 | Y-type with L |
| 6 | YERG | 90,132 | Y_C2 | Y-type with G |
| 7 | FQRL | 84,041 | F_C2 | Q at pos49 |
| 8 | YERW | 69,577 | Y_C2 | W at pos52 (borderline) |
| 9 | FKRG | 65,818 | F_C2 | K at pos49 |
| 10 | YECL | 57,104 | Y_C4 | C4 cysteine pattern |

**Selected for VHH Designer:** FERG, YQRL, FERF, FERA, YERL, FKRG, FQRA

Selection criteria: statistical support (>50K sequences), family diversity (F and Y types), plus FQRA for vernier compatibility despite lower sequence count.

### 4.3 Hallmark Distribution in 12M Sequences

| Hallmark | Count | Family | Description |
|----------|-------|--------|-------------|
| FERG | 2.1M | F_C2 | F42, E49, R50, G52 - most common |
| FERF | 1.8M | F_C2 | F42, E49, R50, F52 |
| YQRL | 1.5M | Y_C2 | Y42, Q49, R50, L52 |
| YERL | 1.2M | Y_C2 | Y42, E49, R50, L52 |
| FKRG | 0.9M | F_C2 | F42, K49, R50, G52 |
| FERA | 0.7M | F_C2 | F42, E49, R50, A52 |

### 4.3 Vernier Zone Analysis

The **vernier zone** is arguably the most important discovery from our 12M sequence analysis. These 17 framework positions don't just "support" CDRs—they actively shape how CDRs are presented to antigens. Getting verniers right is often the difference between a functional nanobody and one that misfolds or loses binding.

**Key Finding: Remarkable Cross-Family Consistency**

When we analyzed vernier preferences across different VHH families (F_C2, Y_C2, F_C4, etc.), we found surprising consistency at core positions:

| IMGT Position | F_C2 | Y_C2 | F_C4 | Cross-Family? |
|---------------|------|------|------|---------------|
| 67 | Y (93%) | Y (94%) | Y (91%) | ✓ Highly conserved |
| 71 | V (91%) | V (94%) | V (89%) | ✓ Highly conserved |
| 76 | F (97%) | F (97%) | F (96%) | ✓ Near-universal |
| 89 | L (97%) | L (92%) | L (95%) | ✓ Highly conserved |
| 91 | M (93%) | M (93%) | M (90%) | ✓ Highly conserved |
| 94 | L (96%) | R (95%) | L (94%) | ⚠ Y_C2 differs |
| 68 | A (72%) | A (86%) | A (68%) | ~ Moderate |
| 69 | D (84%) | D (91%) | D (80%) | ~ Moderate |
| 66 | Y (41%) | N (71%) | Y (38%) | ✗ Family-specific |

This means positions 67, 71, 76, 89, 91 are "safe bets" across all families, while position 66 requires family-specific knowledge. This is why VHH Designer uses **hallmark-specific vernier consensus** rather than one-size-fits-all rules.

**Complete Vernier Zone (17 positions):**

| IMGT Position | What It Affects | F_C2 Consensus |
|---------------|-----------------|----------------|
| 66 | CDR2-FR3 junction | Y (41%) |
| 67 | FR3 stability | Y (93%) |
| 68 | CDR loop support | A (72%) |
| 69 | CDR loop support | D (84%) |
| 71 | Core packing | V (91%) |
| 76 | FR3-CDR3 interface | F (97%) |
| 78 | Core hydrophobics | I (82%) |
| 80 | Structural | - |
| 82 | Loop flexibility | N (78%) |
| 83 | Structural | - |
| 84 | Structural | - |
| 85 | Structural | - |
| 87 | CDR3 presentation | V (70%) |
| 89 | Core stability | L (97%) |
| 91 | CDR3 base | M (93%) |
| 94 | CDR3 support | L (96%) |
| 108 | C-terminal | - |

### 4.4 Compensation Rules (CDR↔Framework Correlations)

The analysis discovers **compensation rules** - when a CDR has certain features, specific framework residues are preferred:

```
"If CDR3 is long (>18 aa) in F_C2, prefer A at IMGT68"
"If CDR3 has negative charge in Y_C2, prefer D at IMGT69"
"If CDR3 has cysteines (extra disulfide), prefer specific FR3 residues"
```

These are learned by:
1. Grouping sequences by family
2. For each framework position, computing CDR feature statistics
3. Comparing residue distributions using Welford's online algorithm
4. Identifying positions where different residues associate with different CDR properties

**Statistical thresholds:**
- Minimum support: 500+ sequences
- Minimum effect size: 1.5 aa CDR3 length difference, or 0.3 charge difference

### 4.5 Epistatic Pairs (Higher-Order Rules)

Some framework mutations work best **in combination**:

```
"IMGT67=S AND IMGT71=V → prefer IMGT94=R (confidence: 87%)"
"IMGT68=A AND IMGT69=D → prefer IMGT76=F (confidence: 91%)"
```

**The 68-69-71 Module Example:**
```
In FERG: 68=A, 69=D, 71=V appears in 92% of sequences
In YERL: 68=A, 69=D, 71=V appears in 87% of sequences
In YQRL: 68=A, 69=D, 71=V appears in 91% of sequences
In FERF: 68=A, 69=D, 71=V appears in 89% of sequences

Conclusion: These three positions form a structural unit.
Test them TOGETHER, not separately.
```

### 4.6 Output Files from Analysis

| File | Contents |
|------|----------|
| `analysis_rules_v7_all_positions.json` | 6,898 CDR↔Framework correlation rules |
| `analysis_vernier_archetypes_v7.json` | Per-hallmark vernier consensus profiles |
| `comprehensive_subfamily_analysis_imgt.xlsx` | Complete statistical database (11 sheets) |

---

## 5. VHH Designer Architecture

### 5.1 Core Concepts

**Scaffolds:**

| Type | Description | Use Case |
|------|-------------|----------|
| **Original** | Your input sequence's framework | Preserve more of your sequence |
| **Universal** | Synthetic "ideal" VHH framework | Maximize VHH-like properties |

**The Universal Framework Sequences:**
```
FR1: QVQLVESGGGLVQPGGSLRLSCAASG  (26 aa)
FR2: WFRQAPGQGLEAVA              (14 aa, hallmark placeholder)
FR3: YYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYC (38 aa)
FR4: WGQGTLVTVSS                 (11 aa)
```

**Why Two Approaches?**

| Aspect | Original Mode | Universal Mode |
|--------|--------------|----------------|
| Mutations needed | Fewer | More |
| Preserves input features | Yes | No |
| Standardized starting point | No | Yes |
| Best when | Input has some VHH character | Comparing across inputs |

### 5.2 CDRs: Never Mutated

```
Antibody Structure:
═══════════════════════════════════════════════════════════
  FR1   │  CDR1  │  FR2   │ CDR2 │    FR3    │  CDR3   │ FR4
═══════════════════════════════════════════════════════════
        │ NEVER  │        │NEVER │           │ NEVER   │
        │MUTATED │        │MUTATE│           │ MUTATED │
```

**Critical Rule:** CDRs define your antibody's binding specificity. VHH Designer never modifies them.

### 5.3 Mutation Layers

Mutations are applied progressively from most essential to most exploratory:

| Layer | What Changes | Purpose |
|-------|--------------|---------|
| 1. Hallmarks | IMGT 42, 49, 50, 52 | Convert VH→VHH signature |
| 2. Verniers | 17 CDR-support positions | Optimize CDR presentation |
| 3. Framework Consensus | Other FR positions | Improve overall VHH-likeness |
| 4. Compensation Rules | CDR-conditional mutations | Handle specific CDR features |

### 5.4 Protected Positions

These positions are **never mutated**:
- CDR1, CDR2, CDR3 (all residues)
- IMGT 22, 23 (critical for folding)
- IMGT 55, 100 (canonical/extra cysteines for disulfide bonds)
- IMGT 103, 104 (CDR3 anchors)

### 5.5 Cysteine Validation

All candidates must pass cysteine validation:

**Odd Cysteine Filter:**
```python
def is_valid_cysteine_count(sequence):
    """Disulfides require pairs - odd count is always invalid."""
    return sequence.count('C') % 2 == 0
```

**C4 Family Enforcement:**
For C4 families (Y_C4, F_C4), both extra cysteines must be present:
```
Canonical cysteines (all VHH):
  - IMGT 23 (FR1) ────┐
  - IMGT 104 (CDR3)   │ Core disulfide
                      └──────────────────┘

Extra C4 cysteines:
  - IMGT 55 (FR2) ────┐
  - IMGT 100 (FR3)    │ Extra disulfide (CDR1-CDR3 bridge)
                      └──────────────────┘
```

---

## 6. The M69 Example: A Real Case Study

### 6.1 What is M69?

M69 is a **mouse antibody** that we want to convert into a nanobody. Let's trace through exactly what VHH Designer does with it.

### 6.2 Step 1: Identify the Input Hallmark

When we analyze M69's sequence, we find:
- Position 42: **I** (isoleucine)
- Position 49: **G** (glycine)
- Position 50: **L** (leucine)
- Position 52: **W** (tryptophan)

**Hallmark: IGLW** - This is a **VH-like** sequence (L at position 50, W at 52). It's designed to pair with a light chain.

### 6.3 Step 2: Identify the Input Verniers

We check all 17 vernier positions to see how M69 compares to typical VHH:

| IMGT | M69 Has | VHH Consensus | Match? | Notes |
|------|---------|---------------|--------|-------|
| 66 | K | Y/N | No | Family-specific position |
| 67 | A | Y | No | Highly conserved in VHH |
| 68 | N | A | No | CDR loop support |
| 69 | E | D | No | CDR loop support |
| 71 | F | V | No | Core packing |
| 76 | A | F | No | Near-universal F in VHH |
| 78 | F | I | No | Core hydrophobics |
| 80 | S | - | - | Variable |
| 82 | T | N | No | Loop flexibility |
| 83 | A | - | - | Variable |
| 84 | S | - | - | Variable |
| 85 | V | - | - | Variable |
| 87 | A | V | No | CDR3 presentation |
| 89 | M | L | No | Core stability |
| 91 | L | M | No | CDR3 base |
| 94 | R | L/R | Maybe | Family-dependent |
| 108 | - | - | - | Variable |

**Key Observation:** M69's verniers are typical human/mouse VH at almost every position. This is expected—M69 was designed to pair with a light chain, not work independently. Converting it to a nanobody requires changing most of these positions.

### 6.4 Step 3: Select Target Hallmarks

VHH Designer targets multiple hallmarks to maximize the chance of finding a compatible scaffold for your CDRs.

#### The Hallmark Universe

Our 12M sequence database contains **128 distinct hallmark patterns**:

```
                    ┌─────────────────────────────────────────┐
                    │     128 HALLMARKS IN DATABASE           │
                    │         (~12M sequences)                │
                    └──────────────────┬──────────────────────┘
                                       │
              ┌────────────────────────┴────────────────────────┐
              ▼                                                 ▼
    ┌─────────────────────┐                        ┌─────────────────────┐
    │   87 TRUE VHH       │                        │   41 VH-LIKE        │
    │   R@50, not W@52    │                        │   W@52 or not R@50  │
    │   (8.8M sequences)  │                        │   (3.2M sequences)  │
    └─────────┬───────────┘                        └─────────────────────┘
              │                                              ⬆
              │                                    (excluded - need light chain)
              ▼
    ┌─────────────────────────────────────────────────────────────────────┐
    │                    TRUE VHH BY FAMILY TYPE                          │
    │                                                                     │
    │   F-type ████████████████████████████████████████  43 hm, 6.7M seq │
    │   Y-type ████████████████                          23 hm, 1.8M seq │
    │   V-type ██                                         7 hm, 0.2M seq │
    │   L-type █                                          4 hm, 0.1M seq │
    │   Other  █                                         10 hm, 0.1M seq │
    └─────────────────────────────────────────────────────────────────────┘
```

#### Family Tree: From Position 42 to Hallmarks

Position 42 defines the major family type. Within each family, position 49 creates subfamilies:

```
F-type (6.7M sequences) ─── Most common, aromatic F fills VL cavity
│
├── FERx (5.8M) ─── FERG ★(3.66M)  FERF ★(1.09M)  FERA ★(362K)  FERE(159K)
├── FGRx (233K) ─── FGRG(174K)    FGRF(25K)      FGRA(18K)
├── FARx (174K) ─── FARG(153K)    FARF(14K)      FARA(7K)
├── FQRx (150K) ─── FQRL(84K)     FQRG(51K)      FQRF(9K)       FQRA ★(6K)
├── FDRx (145K) ─── FDRG(85K)     FDRF(29K)      FDRA(13K)
├── FKRx (83K)  ─── FKRG ★(66K)   FKRA(11K)      FKRF(6K)
└── Other (46K) ─── FRRG, FSRG, ...

Y-type (1.8M sequences) ─── Second most common, aromatic Y fills VL cavity
│
├── YQRx (1.2M) ─── YQRL ★(1.11M) YQRF(43K)      YQRV(22K)
├── YERx (462K) ─── YERL ★(280K)  YERG(90K)      YERF(53K)
└── Other (94K) ─── YARG, YKRL, YGRL, ...

V-type (179K sequences) ─── VH-like at position 42
│
└── VERx (134K) ─── VERG(96K)     VERF(21K)      VERA(10K)

L-type (80K sequences) ─── VH-like at position 42
│
└── LERx (74K)  ─── LERG(54K)     LERF(12K)      LERA(9K)
```

**Key:** ★ = in default pool

#### Sequence Distribution: The Long Tail

```
FERG  ████████████████████████████████████████████████████  3,664,062
YQRL  ███████████████                                       1,105,608
FERF  ██████████████                                        1,085,360
FERA  █████                                                   362,407
YERL  ████                                                    279,551
FGRG  ██                                                      173,979
FERE  ██                                                      159,262
FARG  ██                                                      152,596
FERR  ██                                                      141,830
FERL  █                                                       113,345
───────────────────────────────────────────────────────────────────────
FKRG  █                                                        65,818  ★
...   (47 more with 10K-65K)
...   (36 more with <10K)
FQRA                                                            5,941  ★
```

#### Selection Criteria: Why These 7 Hallmarks?

We balanced three factors:

| Factor | Why It Matters |
|--------|----------------|
| **Statistical support** | More sequences = better rules for vernier mutations |
| **Hallmark diversity** | Different amino acids at 42/49/52 = broader coverage |
| **Vernier compatibility** | Higher similarity to VH = fewer mutations needed |

**Top 15 True VHH Hallmarks by Vernier Similarity to M69:**

| Rank | Hallmark | Sequences | Vernier Match | Position 42 | Selected? |
|------|----------|-----------|---------------|-------------|-----------|
| 1 | **FQRA** | 5,941 | **100%** | F (VHH) | ★ Yes - best vernier match |
| 2 | LERA | 8,709 | 94% | L (VH-like) | No - L@42 won't work alone |
| 3 | VQRL | 24,677 | 82% | V (VH-like) | No - V@42 won't work alone |
| 4 | YQRG | 5,998 | 88% | Y (VHH) | No - low count |
| 5 | YERA | 9,459 | 88% | Y (VHH) | No - low count |
| 6 | VERA | 9,758 | 94% | V (VH-like) | No - V@42 won't work alone |
| 7 | FERI | 33,953 | 88% | F (VHH) | No - I@52 rare |
| 8 | **YQRL** | 1,105,608 | 82% | Y (VHH) | ★ Yes - most constrained |
| 9 | **FERF** | 1,085,360 | 88% | F (VHH) | ★ Yes - 2nd most common |
| 10 | FQRL | 84,041 | 82% | F (VHH) | No - YQRL covers Q@49+L@52 |
| 11 | YERF | 53,446 | 88% | Y (VHH) | No - FERF covers F@52 |
| 12 | YERG | 90,132 | 88% | Y (VHH) | No - FERG covers G@52 |
| 13 | **FERA** | 362,407 | 88% | F (VHH) | ★ Yes - A@52 diversity |
| 14 | **YERL** | 279,551 | 82% | Y (VHH) | ★ Yes - Y-type coverage |
| 15 | FERL | 113,345 | 88% | F (VHH) | No - YERL covers L@52 |
| ... | ... | ... | ... | ... | ... |
| 22 | **FERG** | 3,664,062 | 88% | F (VHH) | ★ Yes - most sequences |
| 28 | **FKRG** | 65,818 | 88% | F (VHH) | ★ Yes - K@49 diversity |

**Why some high-vernier hallmarks were excluded:**

- **LERA, VERA, VQRL** (ranks 2, 3, 6): Position 42 = L or V. These are VH-like—the aromatic (F/Y) at position 42 is critical for VHH stability without a light chain.
- **YQRG, YERA, FERI** (ranks 4, 5, 7): Low sequence counts (<35K) mean less reliable vernier rules.
- **FQRL, YERF, YERG, FERL** (ranks 10-15): Already covered by other hallmarks with same position 52.

**Why FERG despite rank #22 by vernier?**
FERG has **3.66 million sequences**—6x more than #2. This gives us the most reliable statistical patterns for framework mutations.

**Why FKRG despite rank #28?**
It's the only high-count hallmark with **K at position 49**, providing amino acid diversity.

**Why FQRA despite only 6K sequences?**
It's the only true VHH with **100% vernier match** to typical VH inputs. Fewer mutations needed = potentially easier conversion.

#### The Final Default Pool

```
DEFAULT_HALLMARK_POOL = [
    'FERG',   # 3.66M seqs - statistical champion
    'YQRL',   # 1.11M seqs - most constrained scaffold  
    'FERF',   # 1.09M seqs - F@52 diversity
    'FERA',   # 362K seqs  - A@52 diversity
    'YERL',   # 280K seqs  - Y-type coverage
    'FKRG',   # 66K seqs   - K@49 diversity
    'FQRA',   # 6K seqs    - vernier compatibility champion
]
```

**Coverage provided:**

| Position | Amino Acids Covered |
|----------|---------------------|
| 42 | F (5 hallmarks), Y (2 hallmarks) |
| 49 | E (5 hallmarks), Q (2 hallmarks), K (1 hallmark) |
| 52 | G (2), F (2), A (2), L (2) |

### 6.5 Step 4: Generate Candidates

VHH Designer generates candidates across multiple hallmarks and tracks:

**Track 0 Controls:** Minimal hallmark-only mutations to test basic compatibility
**Track 1:** Single vernier probes to identify which positions matter
**Track 2-3:** Vernier combinations (short and long motifs)
**Track 4:** Full optimization with all layers

---

## 7. The Track System

The track system organizes candidates by mutation strategy, from minimal changes (Track 0) to full optimization (Track 4). Each track answers a different experimental question.

### 7.1 Overview: What Each Track Tests

```
                        MUTATION COMPLEXITY
    ─────────────────────────────────────────────────►
    
    Track 0         Track 1         Track 2         Track 3         Track 4
    ────────        ────────        ────────        ────────        ────────
    CONTROLS        PROBES          SHORT MOTIFS    LONG MOTIFS     EXPANSIVE
    
    "Does the       "Which single   "Do pairs of    "Do larger      "What's the
    basic CDR       vernier         verniers work   combinations    best overall
    graft work?"    matters most?"  together?"      help?"          design?"
    
    0-15 mut        6-7 mut         7-9 mut         9-12 mut        12-25 mut
    IMMUTABLE       RANKED          RANKED          RANKED          RANKED
```

| Track | Purpose | Mutations | Status |
|-------|---------|-----------|--------|
| **0** | Reference baselines | 0-15 | Controls (never ranked) |
| **1** | Single vernier impact | 6-7 | Ranked |
| **2** | Vernier pairs | 7-9 | Ranked |
| **3** | Vernier triplets+ | 9-12 | Ranked |
| **4** | Full optimization | 12-25+ | Ranked |

---

### 7.2 Track 0: Controls (Ranking-Exempt)

Track 0 candidates are **never ranked** - they serve as experimental controls and reference points. They are marked as **IMMUTABLE** and cannot be modified after generation.

#### 7.2.1 Control Types

**Track 0a: GRAFT (0 framework mutations)**
The absolute baseline - CDRs grafted onto a VHH scaffold with NO framework changes.
This is the **only** track that doesn't include hallmark mutations.

```
Original M69:    EIQLQQSGAELM...GYTFSSYW...ILPGSGST...ARGDDYDEGFPS...
                 ├── FR1 ────┤  ├─CDR1──┤  ├─CDR2──┤  ├───CDR3────┤
                 
Track 0a GRAFT:  EIQLQQSGAELM...GYTFSSYW...ILPGSGST...ARGDDYDEGFPS...
                 └─ unchanged ─┘  └─ kept ─┘  └─ kept ─┘  └── kept ──┘
                 
Mutations: 0 (pure CDR graft, framework unchanged)
Purpose: "Do the CDRs work AT ALL in this scaffold without any VHH conversion?"
```

**Track 0b: Hallmarks + IMGT2 (5 mutations)**
The minimum VHH conversion - just the 4 hallmark positions plus IMGT position 2.

```
Position:        42    49    50    52    2
                 ▼     ▼     ▼     ▼     ▼
Original M69:    I     G     L     W     I
                 │     │     │     │     │
Target FERG:     F     E     R     G     V
                 
Mutations shown:
  IMGT42: I → F  (🔴 Hallmark - fills VL cavity)
  IMGT49: G → E  (🔴 Hallmark - polar contact)  
  IMGT50: L → R  (🔴 Hallmark - KEY VHH marker!)
  IMGT52: W → A  (🔴 Hallmark - removes bulky W)
  IMGT2:  I → V  (🟡 Conserved framework position)

Purpose: "Do minimal hallmark changes enable VHH function?"
```

**Track 0c: Full Vernier (14-17 mutations)**
Hallmarks + IMGT2 + ALL vernier positions changed to hallmark-specific consensus values.

```
Track 0c builds on Track 0b, adding all vernier positions:

  Track 0b base:  Hallmarks (4) + IMGT2 (1) = 5 mutations
  Track 0c adds:  + ALL verniers (~10-12)  = 14-17 total mutations

Vernier positions changed (example for FERG):
  
  Position:   2   42   49   50   52   66   67   68   69   71   76   78   82   87   89   91   94
             ─── ──── ──── ──── ──── ──── ──── ──── ──── ──── ──── ──── ──── ──── ──── ──── ────
  M69:        I    I    G    L    W    K    Y    N    E    F    A    F    T    A    M    L    R
              │    │    │    │    │    │    │    │    │    │    │    │    │    │    │    │    │
  FERG:       V    F    E    R    G    Y    Y    A    D    V    F    I    N    V    L    M    L
              
  Changes:   🟡   🔴   🔴   🔴   🔴   🟢   ·    🟢   🟢   🟢   🟢   🟢   🟢   🟢   🟢   🟢   🟢
  
Legend: 🔴 = Hallmark (4)  🟡 = IMGT2 (1)  🟢 = Vernier (~10-12)  · = No change

Purpose: "Does full vernier optimization help or hurt?"
```

**Track 0d-0p: Progressive Framework (YQRL only)**

YQRL is the **most constrained** VHH scaffold - its vernier positions show 93% average conservation, meaning nature has found ONE optimal solution and rarely deviates. This creates a problem: we can't just apply random combinations of mutations because almost any deviation breaks the scaffold.

**Why YQRL Is Different**

YQRL verniers are 10% MORE constrained on average than FERG, meaning YQRL tolerates almost no variation at these positions:

| Position | YQRL Consensus | FERG Consensus | Difference |
|----------|----------------|----------------|------------|
| IMGT67 | Y (94%) | Y (89%) | YQRL +5% |
| IMGT68 | A (86%) | A (68%) | YQRL +18% ← Much tighter! |
| IMGT69 | D (91%) | D (85%) | YQRL +6% |
| IMGT71 | V (94%) | V (92%) | YQRL +2% |
| IMGT76 | F (97%) | F (97%) | Same (both locked) |
| IMGT78 | I (94%) | I (81%) | YQRL +13% ← Much tighter! |
| IMGT89 | L (92%) | L (78%) | YQRL +14% ← Much tighter! |
| IMGT91 | M (93%) | M (72%) | YQRL +21% ← Much tighter! |
| IMGT94 | R (95%) | L (88%) | Different AA entirely! |
| **Average** | **93%** | **83%** | **YQRL +10%** |

**The Progressive Testing Strategy**

Because YQRL is so rigid, we can't use Tracks 1-4 (which test random combinations). Instead, we test framework positions **progressively** - adding mutations in small groups to identify exactly where tolerance runs out.

Each step ADDS to the previous (cumulative):

| Track | New Positions Added | Cumulative Mutations | What We Learn |
|-------|---------------------|----------------------|---------------|
| 0a | (none) | 0 | Pure CDR graft baseline |
| 0b | 42, 49, 50, 52, 2 | 5 | Minimal VHH conversion |
| 0c | + all verniers | ~15 | Full vernier effect |
| 0d | + 91, 94 | ~17 | FW3 core positions |
| 0e | + 24, 25 | ~19 | FR1 positions |
| 0f | + 83, 84 | ~21 | FW3 loop |
| 0g | + 74, 82 | ~23 | VL interface region |
| 0h | + 4, 6 | ~25 | N-terminal |
| 0i | + 19, 20 | ~27 | FR1 core |
| 0j | + 36, 39 | ~29 | CDR1 boundary |
| 0k | + 45, 47 | ~31 | FR2 positions |
| 0l | + 77, 79 | ~33 | CDR2 boundary |
| 0m | + 85, 86 | ~35 | FW3 extended |
| 0n | + 87, 88 | ~37 | FW3 C-terminal |
| 0o | + 40, 41 | ~39 | FR2 extended |
| 0p | + remaining | ~41 | Complete framework |

**Why Pairs Instead of Singles?**

Framework positions often work together structurally. Testing them in pairs reduces the number of controls needed (13 steps vs 27 individual), captures cooperative effects, and provides interpretable results since each step tests a structural region.

**Example Interpretation**

```
Experimental Results:
  Track 0b (5 mut):  pLDDT = 85  ✓ Good
  Track 0c (15 mut): pLDDT = 82  ✓ Good  
  Track 0d (17 mut): pLDDT = 80  ✓ Good
  Track 0e (19 mut): pLDDT = 78  ✓ Acceptable
  Track 0f (21 mut): pLDDT = 65  ⚠ Degraded
  Track 0g (23 mut): pLDDT = 45  ✗ Failed

Conclusion: Positions 83/84 (added in 0f) are tolerated but borderline.
            Positions 74/82 (added in 0g) break the structure.
            → For this CDR set, stop at Track 0e mutations for YQRL scaffold.
```
```

Purpose: "Which framework changes does YQRL tolerate for these specific CDRs?"
```

#### 7.2.2 Why Controls Are Immutable

Controls must remain unchanged so they serve as reliable reference points:
- Track 0a (GRAFT) always has exactly 0 framework mutations
- Track 0b (5-mut) always has exactly 5 mutations
- Track 0c (full vernier) always has all verniers

If these were modified, we couldn't compare experimental candidates against them.

---

### 7.3 Track 1: Single Vernier Probes

Track 1 tests each vernier position **individually** to identify which ones matter most for your specific CDRs.

#### 7.3.1 Load-Bearing Verniers

Not all 17 vernier positions are equally important. Through analysis of 12 million VHH sequences, we identified 12 positions that have the highest impact on CDR loop stability - the "load-bearing" verniers:

```python
LOAD_BEARING_VERNIERS = [15, 66, 68, 69, 71, 76, 78, 82, 87, 89, 91, 94]
```

**Why These 12 Positions?**

| Position | Region | Role | Why It's Critical |
|----------|--------|------|-------------------|
| IMGT15 | FR1 | Core | Packs against CDR3 base, anchors loop orientation |
| IMGT66 | FR3 | Hinge | CDR2-FR3 transition, controls CDR2 exit angle |
| IMGT67 | FR3 | Hinge | Works with 66, often Y (aromatic stacking) |
| IMGT68 | FR3 | Core | Hydrophobic core, packs under CDR2 |
| IMGT69 | FR3 | Core | Salt bridge network, often D (charged) |
| IMGT71 | FR3 | Core | Deep framework, packs against CDR1/CDR2 base |
| IMGT76 | FR3 | Anchor | Nearly invariant F - locks CDR3 N-terminus |
| IMGT78 | FR3 | Core | Hydrophobic packing, I preferred |
| IMGT82 | FR3 | Loop | CDR2 support, variable (N, S, T common) |
| IMGT87 | FR3 | Core | β-sheet hydrogen bonding |
| IMGT89 | FR3 | Core | Hydrophobic core, L highly conserved |
| IMGT91 | FR3 | Hinge | CDR3 N-terminal support |
| IMGT94 | FR3 | Anchor | Directly contacts CDR3, position-critical |

The other 5 vernier positions (2, 4, 35, 47, 93) are less critical: IMGT2 is always mutated anyway as part of Track 0b, while positions 4, 35, and 47 are more peripheral with lower conservation.

**CDR Support Architecture**

```
                    CDR1          CDR2           CDR3
                   ╭─────╮       ╭─────╮        ╭──────╮
                   │     │       │     │        │      │
    FR1────────────┴─────┴───FR2─┴─────┴───FR3──┴──────┴──FR4
         ▲                        ▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲
         │                        ││││││││││││││││
        15                       66│68│71│78│87│91│
                                  67 69 76 82 89 94
                                  
    Position 15 supports CDR3 from below (FR1)
    Positions 66-94 form the FR3 "platform" supporting CDR2 and CDR3
```

**How Track 1 Tests Load-Bearing Verniers**

Track 1 tests these 12 positions one at a time to identify which ones matter most for your specific CDRs:

```
Track 1 Probe for IMGT Position 68:

  Base (Track 0b):     ...F....E....R....G....V...     (5 mutations)
                          42   49   50   52    2
                          
  + Position 68:       ...F....E....R....G....V...A...  (6 mutations)
                          42   49   50   52    2   68
                                                    ▲
                                              Single vernier probe
                                              
  This tests: "Does changing position 68 from N→A help?"
```

#### 7.3.2 Example Track 1 Candidates

```
ID                          | Mutations                           | Purpose
──────────────────────────────────────────────────────────────────────────────
M69_Orig_FERG_T1_v68_r001   | Hallmarks + IMGT2 + IMGT68         | Test pos 68
M69_Orig_FERG_T1_v69_r002   | Hallmarks + IMGT2 + IMGT69         | Test pos 69  
M69_Orig_FERG_T1_v71_r003   | Hallmarks + IMGT2 + IMGT71         | Test pos 71
M69_Orig_FERG_T1_v76_r004   | Hallmarks + IMGT2 + IMGT76         | Test pos 76
...                         | ...                                 | ...
```

#### 7.3.3 Sprinkling: Probabilistic Framework Mutations

Sprinkling is a technique where we **probabilistically add framework consensus mutations** beyond the core vernier changes. Instead of deterministically applying all framework mutations, each position has an independent probability of being mutated based on its conservation frequency.

**The Problem Sprinkling Solves**

Without sprinkling, Track 1 tests ONLY the single vernier position. But that vernier might need framework support to work properly. We might falsely conclude "position 68 doesn't help" when the reality is "position 68 helps BUT needs positions 24+25 to stabilize it."

With sprinkling, some Track 1 candidates get position 68 alone, while others get position 68 plus a few framework mutations. This lets us detect synergies between verniers and framework.

**How Sprinkling Works**

For each candidate, sprinkling considers FR1/FR4 consensus positions. Each position is an independent Bernoulli trial - we roll a random number and compare it to the position's conservation frequency scaled by temperature.

| Position | Consensus AA | Conservation | Sprinkle Probability |
|----------|--------------|--------------|----------------------|
| IMGT1 | Q | 92% | 92% × temperature |
| IMGT24 | G | 88% | 88% × temperature |
| IMGT25 | A | 85% | 85% × temperature |
| IMGT128 | S | 95% | 95% × temperature |

**Example calculation with temperature = 0.8:**
- Position 128 (95% conserved): 95% × 0.8 = 76% chance of mutation
- Position 24 (88% conserved): 88% × 0.8 = 70% chance of mutation
- Position 25 (85% conserved): 85% × 0.8 = 68% chance of mutation

**The Temperature Parameter (--sprinkle-temp)**

Temperature controls how aggressively sprinkling applies mutations:

| Temperature | Effect | Use Case |
|-------------|--------|----------|
| 1.0 | Full | Apply mutations at natural frequency |
| 0.8 | Medium | Slightly more conservative (default) |
| 0.5 | Low | Only apply highly-conserved positions |
| 0.3 | Minimal | Very few sprinkled mutations |

**Why Probabilistic Instead of Deterministic?**

With a deterministic approach, every candidate gets the same framework mutations, making it impossible to separate "vernier effect" from "framework effect." The probabilistic approach creates diversity - candidates vary in which framework positions are mutated, allowing us to identify synergies like "vernier X works best with framework Y."

**Sprinkling by Track**

Track 1 probes have a **30% chance** of including sprinkling:

```
Without sprinkling:  Hallmarks + IMGT2 + one vernier  = 6 mutations
With sprinkling:     + 1-3 FR1/FR4 consensus positions = 7-9 mutations
```

This helps identify interactions between verniers and framework.

---

### 7.4 Track 2: Short Motifs (2-3 Verniers)

Track 2 tests **pairs of verniers** that are known to co-occur in natural VHH sequences.

#### 7.4.1 Why Pairs Matter

Some vernier positions work together - changing one without the other can be destabilizing:

```
Natural VHH database shows:
  Position 68 (A) and Position 69 (D) co-occur in 89% of FERG sequences
  Position 69 (D) and Position 71 (V) co-occur in 92% of FERG sequences
  
So we test these PAIRS together, not just individually.
```

#### 7.4.2 Example Track 2 Candidates

```
Motif: IMGT68 + IMGT69 (the "68-69 pair")

  Base:           ...F....E....R....G....V...         (5 mutations)
                     42   49   50   52    2
                     
  + Motif 68-69:  ...F....E....R....G....V...A...D... (7 mutations)
                     42   49   50   52    2   68  69
                                               ▲───▲
                                              Vernier pair

Mutation details:
  🔴 IMGT42: I→F  (hallmark)
  🔴 IMGT49: G→E  (hallmark)
  🔴 IMGT50: L→R  (hallmark)
  🔴 IMGT52: W→G  (hallmark)
  🟡 IMGT2:  I→V  (conserved)
  🟢 IMGT68: N→A  (vernier - part of motif)
  🟢 IMGT69: E→D  (vernier - part of motif)
```

#### 7.4.3 Common Track 2 Motifs

```
Motif ID        | Positions      | Biological Rationale
────────────────────────────────────────────────────────────────
68+69           | IMGT68, IMGT69 | Adjacent in structure, form H-bond network
69+71           | IMGT69, IMGT71 | Part of FR3 core
71+76           | IMGT71, IMGT76 | VL interface positions
76+78           | IMGT76, IMGT78 | Hydrophobic core pair
66+67           | IMGT66, IMGT67 | CDR2 boundary positions
```

#### 7.4.4 Sprinkling in Track 2

Track 2 candidates have a **50% chance** of receiving sprinkling (see Section 7.3.3 for full explanation).

The higher rate than Track 1 (50% vs 30%) reflects that vernier pairs create more structural perturbation than singles, so framework support helps stabilize the combined effect. We want roughly half the candidates to test "pair alone" vs "pair + framework."

```
Without sprinkling:  Hallmarks + IMGT2 + 2 verniers  = 7 mutations
With sprinkling:     + 1-3 FR1/FR4 positions         = 8-10 mutations
```

---

### 7.5 Track 3: Long Motifs (3-5 Verniers)

Track 3 tests **larger combinations** of verniers that form functional units.

#### 7.5.1 Epistatic Combinations

These motifs come from epistatic pair analysis - positions that show statistical co-evolution:

```
Motif: "Core vernier cluster" (68-69-71-76)

  Base:           ...F....E....R....G....V...              (5 mutations)
                     42   49   50   52    2
                     
  + Core cluster: ...F....E....R....G....V...A...D...V...F (9 mutations)
                     42   49   50   52    2   68  69  71  76
                                               ▲───▲───▲───▲
                                              4-position motif

These 4 positions form a structural unit in the VHH framework.
```

#### 7.5.2 Example Track 3 Candidates

```
ID                              | Mutations              | Motif Size
───────────────────────────────────────────────────────────────────────
M69_Orig_FERG_T3_motif3_r001    | HM + 2 + 68,69,71     | 3 verniers
M69_Orig_FERG_T3_motif4_r002    | HM + 2 + 68,69,71,76  | 4 verniers
M69_Orig_FERG_T3_motif5_r003    | HM + 2 + 66,68,69,71,78| 5 verniers
```

#### 7.5.3 Sprinkling in Track 3

Track 3 candidates have a **70% chance** of receiving sprinkling (see Section 7.3.3 for full explanation).

This is the highest rate among Tracks 1-3 because larger vernier motifs (3-5 positions) create significant structural changes that benefit most from framework support. The progression across tracks reflects increasing perturbation:

| Track | Vernier Count | Sprinkle Rate | Rationale |
|-------|---------------|---------------|-----------|
| Track 1 | 1 | 30% | Minimal perturbation |
| Track 2 | 2 | 50% | Moderate perturbation |
| Track 3 | 3-5 | 70% | Significant perturbation |
| Track 4 | 8+ | 100% | Full optimization |

```
Without sprinkling:  Hallmarks + IMGT2 + 3-5 verniers = 8-10 mutations
With sprinkling:     + 1-4 FR1/FR4 positions          = 9-14 mutations
```

---

### 7.6 Track 4: Expansive Exploration

Track 4 is the **full optimization** track - testing comprehensive vernier and framework combinations.

#### 7.6.1 Two Lanes of Exploration

Track 4 has two complementary strategies:

```
┌─────────────────────────────────────────────────────────────────────┐
│                         TRACK 4                                      │
├─────────────────────────────┬───────────────────────────────────────┤
│      LANE A: EXPANSIVE      │      LANE B: PRODUCTION              │
├─────────────────────────────┼───────────────────────────────────────┤
│ • Many verniers (8-15)      │ • Fewer verniers (5-10)              │
│ • Framework consensus tiers │ • CDR-conditional rules              │
│ • Broad exploration         │ • Targeted optimization              │
│ • Higher mutation count     │ • Balanced mutation count            │
└─────────────────────────────┴───────────────────────────────────────┘
```

#### 7.6.2 Example Track 4 Candidate (Expansive)

```
Full Track 4 candidate showing all mutation layers:

Original M69 sequence (positions shown):
  
  Position:    2   42   49   50   52   66   67   68   69   71   76   78   82   89   91
              ─── ──── ──── ──── ──── ──── ──── ──── ──── ──── ──── ──── ──── ──── ────
  M69:         I    I    G    L    W    K    Y    N    E    F    A    F    T    M    L
               │    │    │    │    │    │    │    │    │    │    │    │    │    │    │
  Track 4:     V    F    E    R    G    Y    Y    A    D    V    F    I    N    L    M
              ─── ──── ──── ──── ──── ──── ──── ──── ──── ──── ──── ──── ──── ──── ────
              🟡   🔴   🔴   🔴   🔴   🟢   ·    🟢   🟢   🟢   🟢   🟢   🟢   🟢   🟢
              
Legend:
  🔴 = Hallmark mutations (REQUIRED for VHH identity)
  🟢 = Vernier mutations (optimize CDR support)  
  🟡 = Conserved framework (IMGT2)
  ·  = No change (already consensus)

Total mutations: 14
```

#### 7.6.3 Framework Consensus Tiers

Track 4 can include framework consensus mutations beyond verniers. These are applied based on conservation tiers:

| Tier | Conservation | Behavior | Example |
|------|--------------|----------|---------|
| **S** | >95% | Almost always apply | Position 76 (97% F) |
| **A** | 90-95% | Usually apply | Position 24 (92% G) |
| **B** | 80-90% | Sometimes apply | Position 83 (85% A) |
| **C** | 70-80% | Rarely apply | Position 108 (72% L) |

#### 7.6.4 Universal vs Original Scaffold in Track 4

Track 4 generates candidates on **both** scaffolds:

| Scaffold | Starting Point | Best For |
|----------|----------------|----------|
| **Original** | Your input sequence (M69) | When your input has proven expression/stability |
| **Universal** | Standardized humanized VHH framework | Maximum consistency, known good expression |

**Original scaffold** applies mutations on top of your framework, preserving any beneficial features of your original antibody.

**Universal scaffold** applies mutations on top of a well-characterized baseline, providing consistent and predictable behavior.

---

### 7.7 YQRL: The Most Constrained Scaffold

YQRL is treated specially because it's the **most constrained** VHH scaffold. While other hallmarks like FERG show 80-85% vernier conservation, YQRL averages **93%** - meaning nature has converged on a single optimal solution.

#### 7.7.1 Why YQRL Is Special

| Position | Consensus | Frequency | Interpretation |
|----------|-----------|-----------|----------------|
| IMGT67 | Y | 94% | Almost always Y |
| IMGT68 | A | 86% | Strong preference for A |
| IMGT69 | D | 91% | Almost always D |
| IMGT71 | V | 94% | Almost always V |
| IMGT76 | F | 97% | Extremely constrained (locked) |
| IMGT78 | I | 94% | Almost always I |
| IMGT89 | L | 92% | Strong preference for L |
| IMGT91 | M | 93% | Almost always M |
| IMGT94 | R | 95% | Extremely constrained |
| **Average** | | **93%** | **Very constrained overall** |

What 93% conservation means: only 7% of natural YQRL sequences deviate from consensus. Mutations at these positions are rarely tolerated - the scaffold has ONE optimal configuration. Random exploration (Tracks 1-4) would mostly produce failures.

**Comparison: YQRL vs FERG Flexibility**

```
YQRL (rigid):     ████████████████████████████████████░░░ 93% locked
FERG (flexible):  ████████████████████████████░░░░░░░░░░░ 83% locked
```

This 10% difference is HUGE in practice: FERG can explore ~1000 viable combinations, while YQRL can explore only ~50.

#### 7.7.2 YQRL Gets Controls Only (No Tracks 1-4)

Because YQRL is so constrained, random Track 1-4 variations would mostly fail. Instead, YQRL gets **progressive framework controls** that systematically test tolerance.

**What YQRL gets:**
- ✓ Track 0a: GRAFT (0 mutations)
- ✓ Track 0b: Hallmarks + IMGT2 (5 mutations)
- ✓ Track 0c: Full vernier (15 mutations)
- ✓ Track 0d-0p: Progressive framework pairs (13 steps, ~40 mutations cumulative)

**What YQRL doesn't get:**
- ✗ Track 1: No single vernier probes (verniers are already at consensus!)
- ✗ Track 2: No short motifs (can't deviate from consensus)
- ✗ Track 3: No long motifs (can't deviate from consensus)
- ✗ Track 4: No expansive exploration (would produce failures)

Why no vernier probes? YQRL verniers are ALREADY at the consensus values. There's nothing to "probe" - we know position 68 should be A (86%), etc. The question isn't WHICH verniers to use, but WHETHER the CDRs fit at all.

#### 7.7.3 The YQRL "Acid Test"

YQRL serves as a **stress test** for your CDRs. If they work in the most constrained scaffold, they'll likely work anywhere.

**Decision Tree:**

```
                    Run YQRL Track 0b (5 mutations)
                              │
              ┌───────────────┴───────────────┐
              ▼                               ▼
         pLDDT ≥ 75                      pLDDT < 75
              │                               │
              ▼                               ▼
    "CDRs fit YQRL!"                 "CDRs don't fit YQRL"
              │                               │
              ▼                               ▼
    High confidence for              Focus on flexible scaffolds:
    ALL other scaffolds              FERG, FERF, FKRG
    (FERG, FERF, etc.)               (may still work fine!)
```

**Interpreting results:**
- If Track 0b (YQRL) passes → Your CDRs are highly compatible, any hallmark should work, YQRL gives maximum stability
- If Track 0b (YQRL) fails → Your CDRs have specific requirements, but FERG/FERF are more forgiving and will likely still work

#### 7.7.4 Reading YQRL Progressive Results

The Track 0d-0p progressive controls tell you exactly where YQRL's tolerance runs out.

**Example: Your CDRs with YQRL progressive framework**

| Track | Mutations | pLDDT | Status | Interpretation |
|-------|-----------|-------|--------|----------------|
| 0a | 0 | 88 | ✓ | CDRs graft cleanly |
| 0b | 5 | 85 | ✓ | VHH conversion works |
| 0c | 15 | 82 | ✓ | Full verniers OK |
| 0d | 17 | 80 | ✓ | FW3 core OK (91, 94) |
| 0e | 19 | 78 | ✓ | FR1 positions OK (24, 25) |
| 0f | 21 | 72 | ⚠ | FW3 loop marginal (83, 84) |
| 0g | 23 | 58 | ✗ | VL interface fails (74, 82) |
| 0h | 25 | 45 | ✗ | Continued degradation |

**Reading this example:** Your CDRs are compatible with YQRL up to Track 0e (19 mutations). Positions 83/84 are borderline, and positions 74/82 break it.

**Recommendation:** Use YQRL with Track 0e-level mutations for this CDR set, or use FERG which tolerates positions 74/82 better.

---

### 7.8 Visual Summary: Mutation Accumulation by Track

```
Track 0a (GRAFT):
  ■■■■■■■■■■■■■■■■■■■■ CDRs only, no FW changes
  Mutations: 0
  
Track 0b (5-mut):  
  ■■■■■■■■■■■■■■■■■■■■ + 🔴🔴🔴🔴🟡 (hallmarks + IMGT2)
  Mutations: 5
  
Track 0c (Full vernier):
  ■■■■■■■■■■■■■■■■■■■■ + 🔴🔴🔴🔴🟡 + 🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢 (hallmarks + IMGT2 + all verniers)
  Mutations: 14-17

Track 1 (Probe):
  ■■■■■■■■■■■■■■■■■■■■ + 🔴🔴🔴🔴🟡 + 🟢 (hallmarks + IMGT2 + one vernier)
  Mutations: 6-7

Track 2 (Short motif):
  ■■■■■■■■■■■■■■■■■■■■ + 🔴🔴🔴🔴🟡 + 🟢🟢 (hallmarks + IMGT2 + vernier pair)
  Mutations: 7-9

Track 3 (Long motif):
  ■■■■■■■■■■■■■■■■■■■■ + 🔴🔴🔴🔴🟡 + 🟢🟢🟢🟢 (hallmarks + IMGT2 + vernier group)
  Mutations: 9-12

Track 4 (Expansive):
  ■■■■■■■■■■■■■■■■■■■■ + 🔴🔴🔴🔴🟡 + 🟢🟢🟢🟢🟢🟢🟢🟢 + 🟣🟣🟣 (+ framework consensus)
  Mutations: 12-25

Legend:
  ■ = Original sequence (CDRs preserved)
  🔴 = Hallmark mutation (positions 42, 49, 50, 52)
  🟡 = IMGT2 (conserved framework position)
  🟢 = Vernier mutation (CDR-supporting positions)
  🟣 = Framework consensus mutation (beyond verniers)

NOTE: All tracks except 0a include the hallmark mutations - these are REQUIRED for VHH identity.
```

---

## 8. Selection Pipeline

The selection pipeline ensures that the final output is diverse, balanced, and informative. It's not just "pick the top N by score" - that would give you all FERG Universal candidates (which have the best scores but limited diversity).

### 8.1 The 6-Phase Process

```
┌─────────────────────────────────────────────────────────────────┐
│  PHASE 1: Track Minimums (Hallmark-Aware)                       │
│  Each track's minimum allocation distributed across hallmarks   │
│  Prevents any hallmark from dominating a track                  │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│  PHASE 2: Scaffold Balance                                      │
│  Ensure minimum 20% from Original scaffold                      │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│  PHASE 3: Bucket Fill                                           │
│  Fill remaining slots by (hallmark × scaffold × family) buckets │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│  PHASE 4: Hallmark Cap (35%)                                    │
│  No single hallmark exceeds 35% of ranked candidates            │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│  PHASE 5: Track Quota Repair                                    │
│  Restore track balance after hallmark cap enforcement           │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│  PHASE 6: Final Verification                                    │
│  Double-check all constraints satisfied                         │
└─────────────────────────────────────────────────────────────────┘
```

### 8.2 Hallmark Diversity Enforcement

**Why enforce diversity?**

Without enforcement, FERG would dominate because:
1. FERG has 3.66M sequences → best statistical models → best scores
2. Universal scaffold uses FERG → more FERG candidates generated
3. Top-by-score selection would give ~90% FERG

But we WANT diversity because:
- Different hallmarks may work better for different CDR sets
- Experimental diversity provides more information
- Hedges against single-hallmark failures

**How it works:**

**Phase 1: Track Minimums** distribute slots across hallmarks within each track:

| Track | Total Minimum | Per-Hallmark Minimum (7 hallmarks) |
|-------|---------------|-----------------------------------|
| Track 1 | 15 | ~2 per hallmark |
| Track 2 | 19 | ~2-3 per hallmark |
| Track 3 | 19 | ~2-3 per hallmark |
| Track 4 | 31 | ~4 per hallmark |

This ensures every non-rigid hallmark (FERG, FERF, FKRG, FERA, YERL, FERL, FQRA) gets representation in each track.

**Phase 4: Hallmark Cap** prevents any single hallmark from exceeding 35% of ranked candidates:

```
Example with 156 ranked slots:
  Maximum per hallmark: 156 × 35% = 54 candidates
  
  Before cap:
    FERG: 82 candidates  ← Exceeds cap!
    FERF: 35 candidates
    Others: 39 candidates
    
  After cap:
    FERG: 54 candidates  ← Capped at 35%
    FERF: 42 candidates  ← Filled with FERF
    FKRG: 33 candidates  ← Filled with FKRG
    Others: 27 candidates
```

### 8.3 Cross-Track Deduplication (v9.0)

When identical sequences appear in multiple tracks, v9.0 keeps only one:

**Priority:** Track 4 > Track 3 > Track 2 > Track 1

**Tie-breaker:** Better score (lower)

This reduces redundancy in wet lab panels while preserving the most informative version.

### 8.4 Final Output Distribution

| Category | Description |
|----------|-------------|
| **Lead** | Your original input (unchanged) |
| **Controls** | Track 0 references (ranking-exempt, ~30-35) |
| **Ranked** | Optimized candidates from Tracks 1-4 (fills to `--n-select`) |

### 8.5 Understanding the Rank Column

The `rank` column in the output reflects **score order within the selected candidates**, not the raw score order of all generated candidates.

**How rank is assigned:**
1. Controls get ranks [1], [2], [3]... (brackets indicate ranking-exempt)
2. Ranked candidates continue: 34, 35, 36... (sorted by `combined_score`, best first)

**What rank DOES reflect:**
- Score order among the candidates that passed all selection constraints
- Rank 34 has a better `combined_score` than rank 35

**What rank does NOT reflect:**
- Raw score position among ALL generated candidates
- A rank 34 candidate might have been the 500th best by raw score globally, but is the best among those that fit the diversity constraints

**To see raw score ranking**, sort the `_all.csv` file by `combined_score` ascending.

---

## 9. Scoring System

This section describes the **pure sequence ranking** - how we determine which sequences are "better" based solely on their properties. This is independent of hallmark diversity, track quotas, or selection caps (which are applied later in Section 8).

### 9.1 The Goal of Scoring

We generate thousands of candidates. We need to rank them by "how likely is this to work?"

**What "Work" Means:**
1. **Folds properly** - Adopts correct 3D structure
2. **Stays soluble** - Doesn't aggregate
3. **Retains binding** - CDRs still recognize the target
4. **Expresses well** - Can be produced in cells

### 9.2 The Combined Score Formula

The `combined_score` ranks sequences based on **two components** (three if ESMFold is enabled):

```
combined_score = (ESM2_weight × ESM2_normalized) + (Rule_weight × Rule_inverted)

Where:
  ESM2_weight = 0.29 (29%)
  Rule_weight = 0.71 (71%)
  
  ESM2_normalized = (esm2_loss - min_loss) / (max_loss - min_loss)
  Rule_inverted   = 1.0 - weighted_naturalness
```

**Lower combined_score = better candidate.**

| Component | Weight | Raw Value | Better Is | Contribution to Score |
|-----------|--------|-----------|-----------|----------------------|
| ESM2 Loss | 29% | 0.3 - 0.5 typical | Lower | Lower loss → lower score |
| Rule Compliance | 71% | 0.0 - 1.0 | Higher | Higher compliance → lower score |
| pLDDT (optional) | varies | 0 - 100 | Higher | Higher pLDDT → lower score |

**What's NOT in the combined_score:**
- ❌ Hallmark identity (FERG vs FERF etc)
- ❌ Track number
- ❌ Scaffold type (Original vs Universal)
- ❌ Any diversity bonuses

The combined_score is **purely about sequence quality**. A FERG candidate and a FERF candidate with identical ESM2 loss and rule compliance will have the **same combined_score**.

### 9.3 ESM2: The Language Model Score

ESM2 is a protein language model trained on millions of sequences. It predicts how "natural" a sequence looks.

**How It Works:**
- Feed in a sequence
- ESM2 predicts how "surprised" it is by each amino acid
- Lower surprise = more natural-looking = probably more stable

```python
def esm2_score(sequence):
    """
    Returns: (loss, perplexity)
    - loss: Lower is better (0.35 is good, 0.50 is concerning)
    """
    inputs = tokenizer(sequence, return_tensors="pt")
    with torch.no_grad():
        outputs = model(**inputs, labels=inputs["input_ids"])
    return outputs.loss.item()
```

**Typical Values:**

| Score Range | Interpretation |
|-------------|----------------|
| 0.30 - 0.38 | Excellent - looks like a natural nanobody |
| 0.38 - 0.42 | Good - reasonable conversion |
| 0.42 - 0.48 | Marginal - might work, might not |
| > 0.48 | Poor - sequence looks unnatural |

We use the `facebook/esm2_t6_8M_UR50D` model (smallest, fastest).

### 9.4 Rule Compliance Score (weighted_naturalness)

Checks each candidate against rules learned from 12M sequences. The rules are **hallmark-specific** - a FERG candidate is checked against FERG rules, a FERF candidate against FERF rules.

```python
def rule_compliance_score(candidate, rules):
    """
    Returns: 0.0 (violated all rules) to 1.0 (followed all rules)
    """
    passed = 0
    applicable = 0
    
    for rule in rules:
        if rule_applies(candidate, rule):
            applicable += 1
            if candidate_satisfies(candidate, rule):
                passed += 1
    
    return passed / applicable if applicable > 0 else 0.5
```

**Important:** Different hallmarks have different numbers of rules with different confidences. This can cause systematic differences:

| Hallmark | Database Size | # Rules | Typical Compliance |
|----------|---------------|---------|-------------------|
| FERG | 3.66M | 22 | 0.30 - 0.50 |
| FERF | varies | 18-22 | 0.40 - 0.60 |
| FQRA | 6K | 18 | 0.50 - 0.65 |

Because FERG has more sequences, its rules are more stringent (higher confidence thresholds), which can lead to lower raw compliance scores. This is normalized within each batch.

### 9.5 Putting It Together

```
Example ranking (all from same batch):

Rank │ ID                      │ ESM2  │ Rules │ Combined │ Hallmark
─────┼─────────────────────────┼───────┼───────┼──────────┼─────────
  1  │ M69_Univ_FERG_T4_001    │ 0.33  │ 0.85  │ 0.138    │ FERG
  2  │ M69_Univ_FERG_T4_002    │ 0.34  │ 0.83  │ 0.147    │ FERG
  3  │ M69_Orig_FERF_T4_001    │ 0.36  │ 0.80  │ 0.175    │ FERF
  4  │ M69_Univ_FERG_T3_001    │ 0.35  │ 0.78  │ 0.186    │ FERG
  ...
 50  │ M69_Orig_FQRA_T2_001    │ 0.42  │ 0.52  │ 0.463    │ FQRA

Note: FERG dominates top ranks because:
  1. Universal scaffold is FERG-derived → lower ESM2 loss
  2. FERG has most data → best-trained rules
  
This is why Section 8 enforces hallmark diversity AFTER scoring.
```

---

## 10. Output Files

### 10.1 Main Output Files

| File | Contents |
|------|----------|
| `{base}.csv` | Selected candidates with all metadata |
| `{base}.fasta` | Sequences in FASTA format |
| `{base}_summary.json` | Run statistics and parameters |
| `{base}_detailed.jsonl` | Full provenance for each candidate |
| `{base}_all.csv` | All generated candidates (before selection) |

### 10.2 MSA-Compatible ID Format (v9.0)

**Lead:**
```
M69_Orig_IGLW_LEAD
```

**Control:**
```
M69_Orig_YQRL_T0Ctrl_VnFw_s848
│   │    │    │      │    │
│   │    │    │      │    └── Score (0.848 → 848)
│   │    │    │      └─────── Mutations: Vn=verniers, Fw=framework
│   │    │    └────────────── T0Ctrl = Track 0, Control
│   │    └─────────────────── Hallmark pattern
│   └──────────────────────── Orig=Original, Univ=Universal
└──────────────────────────── Sequence prefix
```

**Ranked:**
```
M69_Univ_FERG_T4R_HVn_r001_s152
│   │    │    │   │   │    │
│   │    │    │   │   │    └── Score
│   │    │    │   │   └─────── Rank (001 = best)
│   │    │    │   └─────────── Mutations: H=hallmark, Vn=verniers
│   │    │    └─────────────── T4R = Track 4, Ranked
│   │    └──────────────────── Hallmark
│   └───────────────────────── Scaffold
└───────────────────────────── Prefix
```

### 10.3 Mutation Codes

| Code | Meaning |
|------|---------|
| `H` | Hallmark mutations applied |
| `V1` | Single vernier probe (Track 1) |
| `Vn` | Multiple verniers (Track 4) |
| `M23` | Motif length 2-3 (Track 2) |
| `M35` | Motif length 3-5 (Track 3) |
| `Fw` | Framework consensus applied |
| `X{n}` | Generic (n mutations) |

### 10.4 CSV Column Definitions

| Column | Description |
|--------|-------------|
| `rank` | Selection rank (1 = best scored) |
| `status` | LEAD, CONTROL, or RANKED |
| `id` | MSA-compatible ID for alignments |
| `seq_name` | Human-readable candidate name |
| `track` | Which track generated this candidate |
| `scaffold_type` | original or universal |
| `hallmarks` | The 4-letter hallmark code |
| `sequence` | Full amino acid sequence |
| `combined_score` | Final ranking score |
| `n_cysteines` | Number of cysteines (should be 2 or 4) |

---

## 11. Command Line Reference

### 11.1 Basic Usage

```bash
python vhh_designer_v9.0.py \
    --fasta input.fasta \
    --rules analysis_rules_v7_all_positions.json \
    --archetypes analysis_vernier_archetypes_v7.json \
    --n-generate 50000 \
    --n-select 192
```

### 11.2 Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--fasta` | Required | Input FASTA file |
| `--rules` | Required | Rules JSON from analysis |
| `--archetypes` | Required | Vernier archetypes JSON |
| `--n-generate` | 50000 | Candidates to generate |
| `--n-select` | 192 | Final candidates to select |
| `--hallmark-diversity-cap` | 0.35 | Max fraction per hallmark |
| `--sprinkle-temp` | 0.8 | Framework mutation aggressiveness |

### 11.3 Sprinkle Temperature

Controls how aggressively framework mutations are applied:

| Temperature | Effect |
|-------------|--------|
| **0.5** | More aggressive - mutations applied more readily |
| **1.0** | Neutral - uses natural probabilities as-is |
| **2.0** | Conservative - fewer mutations applied |

**Mathematical Formula:**
```
p' = sigmoid(logit(p) / T)

Where:
- p = natural probability (e.g., 0.91 for IMGT71:V)
- T = temperature
- p' = adjusted probability

Example with T=0.8:
- p=0.91 → p'=0.96 (more likely to mutate)
- p=0.60 → p'=0.72 (moderately more likely)
```

---

## 12. Caveats & Limitations

### 12.1 Data Biases

**Species Bias:**
- ~70% of data is from llama
- ~25% is from alpaca
- ~5% is from camel

Different camelid species may have different optimal patterns.

**Functional Bias:**
- 99% of data lacks antigen/target information
- We're learning "what looks natural," not "what binds well"

### 12.2 What This Tool CANNOT Do

1. **Predict binding affinity** - We preserve your CDRs, but can't guarantee binding is maintained
2. **Guarantee stability** - Statistical patterns don't guarantee biophysical success
3. **Replace experimental validation** - Computational design requires wet lab testing

### 12.3 Recommended Workflow

1. Run VHH Designer to generate candidates
2. Score/filter by ESM2 and rule compliance
3. Select diverse representatives for experimental testing
4. Validate binding and stability in the lab
5. Iterate based on experimental results

---

## Appendices

### A. Protected Positions Summary

```
NEVER MUTATE:
- CDR1 (IMGT 27-38)
- CDR2 (IMGT 56-65)
- CDR3 (IMGT 105-117)
- IMGT 22, 23 (folding critical)
- IMGT 55, 100 (disulfide cysteines)
- IMGT 103, 104 (CDR3 anchors)
```

### B. Hallmark Position Reference

| IMGT Position | Human VH | Classic VHH | Position 52 Variants |
|---------------|----------|-------------|---------------------|
| 42 | V | F or Y | - |
| 49 | G | E, Q, or K | - |
| 50 | L | R | - |
| 52 | W | G (FERG) | F (FERF), L (YERL), A (FERA) |

### C. Vernier Zone Positions

```
IMGT: 66, 67, 68, 69, 71, 76, 78, 80, 82, 83, 84, 85, 87, 89, 91, 94, 108
```

### D. Framework Consensus Tiers

| Tier | Confidence | Example Positions | Description |
|------|------------|-------------------|-------------|
| **S** (Super) | >95% | 8, 9, 41, 75, 98 | Almost universal in VHHs |
| **A** (High) | 90-95% | 7, 14, 22, 51 | Very common |
| **B** (Medium) | 85-90% | 17, 45, 48, 90 | Common but variable |
| **C** (Lower) | 80-85% | 86, 88, 120 | More variation allowed |

### E. Version History

**v9.0** (Current) - MAJOR RELEASE
- **CRITICAL: Fixed scoring direction bug** - combined_score was completely inverted!
  - Before: Higher rule compliance → worse rank (WRONG - ranked worst candidates first!)
  - After: Higher rule compliance → better rank (CORRECT)
  - Now correctly: lower ESM2 loss + higher pLDDT + higher rules = better rank
- **Fixed rank ordering bug** - rank now correctly reflects score order within selected candidates
  - Before: Candidates added during repair/fill phases were appended at end, breaking score order
  - After: Final sort ensures rank reflects score order among selected
- Added FQRA to default hallmark pool (100% vernier match to typical VH inputs)
- New `DEFAULT_HALLMARK_POOL`: FERG, FERF, YQRL, FERA, YERL, FKRG, FQRA (7 hallmarks)
- New `--target-hallmarks DEFAULT` mode (now the default instead of single FERG)
- Fixed `_all.csv` to save ALL scored candidates (was only saving selected ~192)
- **Fixed O(n²) performance bug in `to_dataframe()`** - was hanging on large candidate pools
- File naming now includes n-generate and temperature: `v9_0_n100000_t08_M69`
- Cross-track sequence deduplication
- MSA-compatible ID format (`M69_Orig_FERG_T4R_HVn_r001_s152`)
- Hallmark-aware track minimums
- Hallmark cap enforcement (35% max per hallmark)
- 6-phase selection pipeline for balanced output

| Version | Key Changes |
|---------|-------------|
| **v8.9.3** | Cross-track dedup, MSA IDs, hallmark-aware selection |
| **v8.9.2** | Track quota repair |
| **v8.9.1** | Cysteine protection, immutable validation, hallmark diversity cap |
| **v8.9** | YQRL framework pairs, simplified controls |
| **v8.0** | Universal scaffold, two-lane Track 4 |
| **v7.16** | Strict cysteine validation, C4 enforcement |

---

*Generated for VHH Designer v9.0*
*Based on analysis of 12 million natural VHH sequences*

---

### Acknowledgments

This work would not have been possible without the publicly available sequence databases (OAS, INDI, SAbDab) that have made large-scale antibody analysis accessible to the research community. The development of VHH Designer involved extensive computational biology work—from processing terabytes of raw sequence data through standardized pipelines, to developing statistical frameworks for extracting biologically meaningful patterns, to integrating modern ML validation (ESM2) into the design workflow.

### A Note on Responsible Use

Nanobody engineering, like all protein engineering, carries dual-use considerations. The sequences and methods described here are intended to support legitimate therapeutic development—converting existing antibodies to formats with improved stability, tissue penetration, and manufacturability. The statistical patterns derived from natural camelid sequences represent nature's solutions to the problem of stable, soluble single-domain antibodies; this tool simply makes those patterns accessible for rational design.

Researchers using this tool should follow institutional biosafety guidelines and consider the downstream applications of their engineered sequences. The transparency and traceability built into this pipeline—every mutation traceable to its statistical support in the natural database—is intended not just to help researchers make better design decisions, but also to maintain clear documentation of the rationale behind each modification.
