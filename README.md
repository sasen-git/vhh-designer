# VHH-Designer

**A constrained, evaluation-driven framework for antibody-to-nanobody conversion**
Strategy & system design by Sasen Efrem; implementation developed with Claude

![Version](https://img.shields.io/badge/version-9.1-blue)
![Python](https://img.shields.io/badge/python-3.8+-green)
![Sequences](https://img.shields.io/badge/training%20data-12M%20sequences-orange)

## Overview

VHH-Designer is a research framework for converting conventional antibodies (VH) into camelid-like single-domain antibodies (VHH/nanobodies) using **statistical patterns** derived from large-scale natural repertoires (approximately 12 million sequences that were manually compiled).

Rather than unconstrained sequence generation, the system emphasizes **disciplined exploration**: selectively mutating framework residues while preserving CDR binding loops, enforcing biological invariants, and evaluating outcomes through staged control tracks.

Design decisions throughout the system prioritize interpretability, traceability, and bounded optimization, ensuring that generated candidates remain evolutionarily plausible and biologically coherent.

## Design Philosophy: Safety, Constraints, and Discernment

VHH-Designer was intentionally built with explicit constraints, staged exploration, and built-in evaluation checkpoints to prevent uncontrolled or abberant optimization.

Instead of maximizing generative breadth, the framework prioritizes:

- Bounded exploration over brute-force generation
- Progressive capability evaluation rather than end-to-end optimization
- Early surfacing of failure modes through explicit controls and preservation of biological invariants
  - I.e. Controls in Track 0 1. Controls define the minimum viable perturbation. Track 0 asks:
    - If I change almost nothing, does the system still behave coherently?
  - This question is important because:
      - If small changes already cause instability, the system is fragile
      - If small changes produce large, unexplained effects, assumptions are wrong
      - If controls fail, optimization is meaningless
  - As such, controls and slowly-evolving tracks allow for early surfacing of what might have break the system
- Key architectural choices—such as immutable control tracks, scaffold-specific exclusion rules, and hallmark-aware mutation logic—were chosen to avoid overgeneralization and to distinguish genuine biological signal from statistical coincidence.

The goal is not autonomous design, but legible, stress-tested scientific exploration that supports downstream human judgment and validation.

## Scope and Non-Goals

This framework intentionally avoids:

- Fully autonomous or end-to-end optimization without intermediate evaluation
- Direct synthesis or deployment recommendations
- Inference of biological activity beyond structural and statistical plausibility

Generated sequences are not intended for experimental use without *independent review and validation*. Design constraints are chosen conservatively to limit extrapolation beyond observed biological distributions

## Key Features

- **Hallmark-based classification**: 128 distinct VHH subfamilies identified by tetrad positions (42, 49, 50, 52)
- **Track system**: Systematic experimental design with controls (T0) and optimization tiers (T1-T4)
- **ESM2 scoring**: Protein language model integration for naturalness assessment
- **CDR preservation**: Maintains binding loops while converting framework
- **Multi-scaffold support**: Original framework retention + Universal scaffold grafting

## Project Structure

```
vhh-designer/
├── vhh_designer/           # Main design tool (v9.1 - 7,664 lines)
│   └── designer.py
│
├── database/               # NPZ scanning and shard annotation
│   ├── npz_fullscan.py
│   └── annotate_shards.py
│
├── analysis/               # Statistical analysis pipeline
│   ├── epistasis_analyzer.py
│   ├── compensation_rules.py
│   └── naturalness_analyzer.py
│
├── visualization/          # MSA and alignment visualization
│   ├── msa_visualizer.py
│   └── csv2msa.py
│
├── utilities/              # Helper scripts
│   ├── dna_translator.py
│   ├── pull_cdrs.py
│   └── paths.py
│
├── docs/                   # Documentation
│   └── VHH_DESIGNER_ULTIMATE_GUIDE.md
│
└── archive/                # Key version milestones
    ├── designer/           # v2, v5, v7.1, v7.12, v8.0, v8.9, v9.0
    ├── scanner/            # v2, v5, v6
    └── epistasis/          # v1, v2
```

## Quick Start

```bash
# Install dependencies
pip install -r requirements.txt

# Run designer
python vhh_designer/designer.py \
    -i input.fasta \
    --rules analysis/rules_v7.json \
    --archetypes analysis/archetypes_v7.json \
    --n-generate 100000 \
    --n-select 192
```

## The Science

### What Makes VHH Different from VH?

| Position (IMGT) | Human VH | Camelid VHH | Function |
|-----------------|----------|-------------|----------|
| 42 | Conserved | F or Y | Hydrophobic core |
| 49 | G | E, Q, or K | Solubility |
| 50 | L | R | **Critical** - removes VL dependency |
| 52 | W | G or L | CDR3 packing |

### Track System

| Track | Purpose | Mutations | Use Case |
|-------|---------|-----------|----------|
| **T0** | Controls | 0-5 | Test fundamental compatibility |
| **T1** | Vernier probes | Single | Identify critical positions |
| **T2** | Motif pairs | 2-3 | Test co-occurring patterns |
| **T3** | Motif triplets | 3-5 | Larger structural motifs |
| **T4** | Optimized | 8-20 | Production candidates |

## Data Foundation

| Source | Sequences | Description |
|--------|-----------|-------------|
| INDI_NGS | 10.7M | High-throughput camelid repertoires |
| OAS_Camel | 1.36M | Observed Antibody Space camelid subset |
| **Total** | ~12M | Deduplicated, IMGT-numbered |

## Key Discoveries

1. **128 hallmark families** enumerated with sequence counts
2. **Hallmark-specific vernier preferences** (not just family-level)
3. **Position 50 L→R predicts CDR3 charge characteristics**
4. **YQRL paradox** solved with track-based experimental design
5. **Cysteine classification** requires position-based (not count-based) logic

## Code Evolution

| Phase | Dates | Versions | Focus |
|-------|-------|----------|-------|
| Database Building | Oct-Nov 2025 | scanner v1-6 | 12M sequence curation |
| Statistical Discovery | Nov-Dec 2025 | epistasis v1-2 | Hallmarks, compensation rules |
| Designer Development | Jan 2026 | v2-v8 | Track system, scoring |
| Production Release | Feb 2026 | v9.0-9.1 | Scoring fix, documentation |

**Total: 45 versions, 540 → 7,664 lines of code**

See [DEVELOPMENT_LOG.md](DEVELOPMENT_LOG.md) for the full journey.

## Documentation

- [DEVELOPMENT_LOG.md](DEVELOPMENT_LOG.md) - Project evolution and key discoveries
- [docs/VHH_DESIGNER_ULTIMATE_GUIDE.md](docs/VHH_DESIGNER_ULTIMATE_GUIDE.md) - Complete technical reference



Authorship & Attribution

Scientific strategy, system design, and evaluation framework: Sasen Efrem

Implementation and iterative development: Claude (LLM-assisted coding)

This project reflects a collaborative workflow in which human scientific judgment defined constraints, safeguards, and evaluation logic, with AI used as an implementation accelerator rather than an autonomous designer.

License

MIT
