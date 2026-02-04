# VHH-Designer

**Antibody-to-nanobody conversion using statistical patterns from 12M camelid sequences**

![Version](https://img.shields.io/badge/version-9.1-blue)
![Python](https://img.shields.io/badge/python-3.8+-green)
![Sequences](https://img.shields.io/badge/training%20data-12M%20sequences-orange)

## Overview

VHH-Designer converts conventional antibodies (VH) into camelid-like single-domain antibodies (VHH/nanobodies) by intelligently mutating framework residues while preserving CDR binding loops.

The tool uses statistical patterns learned from ~12 million natural camelid VHH sequences to guide the conversion, ensuring the resulting nanobodies are structurally sound and evolutionarily plausible.

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

## License

MIT

## Author

Developed as part of antibody engineering research, 2025-2026.
