# KA-Search VHH Analysis

Nanobody/VHH sequence analysis using KA-Search and epistasis modeling.

## Directory Structure

```
KA-Search/
├── active/          # Current working scripts
├── archive/         # Historical versions (organized by project)
├── data/            # Input data (sequences, databases, structures)
├── models/          # Trained models and checkpoints
├── results/         # Analysis outputs
├── docs/            # Documentation
└── legacy/          # Backup of original Archive/
```

## Quick Start

### Run Naturalness Analysis
```bash
python active/analysis/vhh_naturalness_analyzer_v3.py \
    --input data/raw/sequences/your_sequences.csv \
    --epistasis models/epistasis/current/epistasis_v2_full.pkl
```

### Run NPZ Scan
```bash
python active/database/npz_fullscan_v6_integrated.py \
    --input your_sequences.csv \
    --database data/databases/production/VHH_db_unified_v2.npz
```

## Key Models

- **Epistasis Model**: `models/epistasis/current/epistasis_v2_full.pkl`
  - Vernier cluster patterns from 11.5M camelid sequences
  - CDR3 length/charge/aromatic statistics per cluster
  
- **Correlation Results**: `models/correlations/correlation_results.pkl`
  - Position-specific CDR-framework coupling rules
  - Triplet co-occurrence statistics

## Documentation

See `docs/` for:
- `epistasis_methodology.md` - How the epistasis model works
- `naturalness_scoring.md` - Naturalness score interpretation
- `database_schema.md` - NPZ database format

## Maintenance

Use `kasearch_maintain.py` to:
- Check naming conventions: `python kasearch_maintain.py check`
- Create new scripts: `python kasearch_maintain.py new my_script`
- Archive old versions: `python kasearch_maintain.py archive script_name`

Reorganized: 2025-12-31