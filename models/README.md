# models/

Trained model outputs and statistical analysis results. All binary files (PKL, NPZ)
are gitignored. Only JSON config/rule files are tracked.

## Subfolders

| Folder | Contents | Status |
|--------|----------|--------|
| `analysis/` | Unified v7.x analysis run outputs (rules, archetypes, correlations) | **Production** — see analysis/README.md |
| `epistasis/` | Legacy epistasis pipeline outputs (pre-unified) | **Superseded** by analysis/ |
| `correlations/` | Legacy correlation/compensation pipeline outputs (pre-unified) | **Superseded** by analysis/ |
| `Testing-Walking-charge/` | Experimental charge-walking analysis — did not produce useful results | **Abandoned test** |

## Loose Files

| File | What it is |
|------|-----------|
| `family_fw_consensus.json` | Per-family framework consensus sequences (tracked) |
| `hallmark_fw_consensus.json` | Hallmark-grouped framework consensus sequences (tracked) |
| `claude.zip` | Archived Claude-generated scripts from early iterations |

## Current Production Model

The designer (`vhh_designer_v9_1.py`) reads from `analysis/`:

```bash
--rules   models/analysis/v7.7_2026-01-13_221057/analysis_rules_v7_all_positions.json
--archetypes models/analysis/v7.7_2026-01-13_221057/analysis_vernier_archetypes_v7.json
```

See `analysis/README.md` for full details on each JSON file.
