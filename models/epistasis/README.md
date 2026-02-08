# models/epistasis/

Legacy epistasis analysis outputs from the standalone `vhh_epistasis_overnight_*.py`
scripts (pre-unified pipeline). These used substring-indexed positions (not IMGT)
and PKL output format.

**Superseded by** `models/analysis/` which uses IMGT numbering and JSON output.

## Subfolders

| Folder | Contents |
|--------|----------|
| `current/` | epistasis_v2_full.pkl (36MB, substring positions), imgt_rules.json, correlation_results_v4_epistasis_only.pkl |
| `checkpoints/` | v2 intermediate checkpoints (compensation, clusters, models, rules) |
| `legacy/` | epistasis_overnight_full.pkl (12MB, v1 from Dec 7) |
| `legacy_checkpoints/` | v1 intermediate checkpoints |
| `imgt_v1/` | First attempt at IMGT epistasis — 185-byte checkpoint PKLs (effectively empty/failed runs) |
| `imgt_light_v1/` | Light IMGT attempt — 318-byte PKL (also failed) |

## Note

The `imgt_v1/` and `imgt_light_v1/` folders contain PKLs that are only ~185 bytes,
meaning the runs likely failed or produced no meaningful data. These were the early
attempts that led to building the unified `vhh_analysis_unified_v7.7.py` pipeline,
which successfully completed the IMGT analysis.
