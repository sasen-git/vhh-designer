# models/correlations/

Legacy correlation and compensation analysis outputs from standalone scripts.

**Superseded by** `models/analysis/` unified pipeline.

## Files

| File | Size | What it is |
|------|------|-----------|
| `compensation_results.pkl` | 22KB | v1 compensation output (Dec 6) |
| `correlation_results_v3_compensation.pkl` | 194KB | v3 merged compensation + correlation (Jan 2) |
| `correlation_results_v4_merged.pkl` | 533KB | v4 merged with epistasis (Jan 3) |
| `imgt_csv_v3_sanity/correlation_results_v3_summary.json` | tracked | Sanity check summary from IMGT CSV correlation run |
| `imgt_csv_v3_sanity/correlation_rules_v3.json` | tracked | Rules from IMGT CSV correlation sanity check |

## Note

These PKLs use the old schema (`multi_position_rules`, `top_predictors`). The
current designer (v9.1) expects the new JSON schema from `models/analysis/`.
These files are kept for reference but are not used in the current pipeline.
