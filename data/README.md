# data/ (51GB)

All sequence data for the VHH design pipeline. Entirely gitignored.

See also: `VHH_DATABASE_CONSOLIDATION_README.md` in this folder for detailed
provenance and consolidation history.

## Subfolders

| Folder | Size | Contents |
|--------|------|----------|
| `databases/` | ~47GB | Production NPZ shards, annotated CSVs, legacy databases |
| `archive/` | ~3GB | Older annotation versions (v3, pre-dedup) |
| `raw/` | ~1GB | Original source files before processing |

## Data Sources (~12M sequences total)

| Source | Sequences | Description |
|--------|-----------|-------------|
| OAS Camel | 1.35M | Observed Antibody Space — camelid repertoires |
| INDI NGS | 10.7M | INDI Molecular nanobody NGS data (6 shards) |
| Patents/Structures | 19K | Patent sequences, SAbDab/PDB structures, curated VHHs |

## Pipeline

```
raw/ sources → build scripts → databases/shards/*.npz
                             → annotate_all_shards_v8.py
                             → databases/annotated/*.csv (with IMGT positions)
                             → vhh_analysis_unified_v7.7.py → models/
```
