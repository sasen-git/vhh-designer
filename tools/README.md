# tools/

Maintenance and QC utilities.

## Scripts

| Script | Purpose |
|--------|---------|
| `dedup_annotated_imgt.py` | Deduplicates IMGT-annotated CSV shards by sequence identity |
| `dedupe_csv_headers.py` | Fixes duplicate column headers in CSVs |
| `qc_imgt_coverage_fast.py` | Quick QC check: what % of IMGT positions are filled per shard |
| `diagnose_unknown_files.py` | Identifies and categorizes unknown files in the repo |
| `kasearch_paths.py` | Centralized path management for all KA-Search scripts |
| `kasearch_reorganize.py` | Repo reorganization helper |
| `maintain.py` | General maintenance (cleanup temp files, check consistency) |
| `summarize_sweep.py` | Summarizes parameter sweep results across multiple runs |
| `pathlib_list_v2.py` | File listing with filtering |

## Other

| File | What it is |
|------|-----------|
| `.cache.json` | Cached metadata for maintenance scripts |
