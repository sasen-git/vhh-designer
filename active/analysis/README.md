# active/analysis/

Current production scripts for the VHH design pipeline.

## Core Pipeline

| Script | Purpose |
|--------|---------|
| `vhh_designer_v9_1.py` | Main VHH designer — converts VH→VHH using learned rules (JSON input) |
| `vhh_analysis_unified_v7.7.py` | Unified epistasis + compensation + correlation analysis (IMGT, JSON output) |
| `vhh_naturalness_analyzer_v4.py` | Scores sequences for naturalness/evolutionary fitness |
| `annotate_all_shards_v8.py` | Annotates database shards with IMGT positions |

## Supporting Scripts

| Script | Purpose |
|--------|---------|
| `augment_imgt_positions_from_csvs_v2.py` | Augments CSVs with full IMGT position columns |
| `charge_balance_vhh_windows.py` | Charge balance analysis across VHH windows |
| `extract_script_headers.py` | Extracts docstrings/headers from all scripts |
| `patches.py` | Hotfix patches for data issues |
| `run_with_positions.py` | Runs analysis with specific position sets |

## Data Flow

```
NPZ shards (12M sequences)
    ↓
annotate_all_shards_v8.py → IMGT-augmented CSVs
    ↓
vhh_analysis_unified_v7.7.py → JSON rules (epistasis + compensation + correlation)
    ↓
vhh_designer_v9_1.py → designed VHH sequences (CSV + FASTA)
    ↓
../utilities/csv2msa_antpack_v6.py → MSA visualization
vhh_naturalness_analyzer_v4.py → fitness scoring
```

## Version History

All older iterations are in `archive/` organized by component:
- `archive/designer/` — v2 through v9_0
- `archive/analysis_unified/` — v7 through v7.6
- `archive/compensation/` — v1 through v5
- `archive/standalone_imgt/` — pre-unified IMGT scripts (epistasis_v3, compensation_v6, correlations_v3)
