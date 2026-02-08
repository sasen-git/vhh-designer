# active/database/

## Scripts

| Script | Purpose |
|--------|---------|
| `npz_fullscan_v6_integrated.py` | Scans NPZ databases for sequence matches with integrated filtering |

## Database Location

Production databases are in `data/databases/`:
- `shards/` — 8 NPZ files (~12M sequences total)
- `production/` — consolidated VHH_db_unified_v2.npz
- `annotated/` — IMGT-augmented CSV versions
- `legacy/` — older database versions

See `data/VHH_DATABASE_CONSOLIDATION_README.md` for full details.
