# results/

All pipeline output — designer runs, NPZ scans, analysis results. Entirely gitignored.

## Subfolders

| Folder | Contents |
|--------|----------|
| `design_runs/` | 60 timestamped designer output directories (v7.5 through v9.1) |
| `analysis_runs/` | Early designer runs (v2-v7) + naturalness analysis outputs |
| `npz_scans/` | NPZ database scan results organized by date |
| `analysis_outputs/` | Loose analysis CSVs/XLSXs moved from root |

## Design Runs

Each `design_runs/YYYYMMDD_HHMMSS_vX_Y_nNNNNN_M69/` folder contains:
- `*.csv` — selected candidates
- `*_all.csv` — all generated candidates
- `*.fasta` — sequences in FASTA format
- `*_summary.json` — run configuration and statistics
- `*_csv2msa_antpack_*/` — MSA alignment subfolder (if post-processed)

## Notable Runs

| Run | Designer | Candidates |
|-----|----------|------------|
| `20260201_234207_v9_1_n100000_t08_M69/` | v9.1 (latest) | 100K generated |
| `20260201_221918_v9_0_n100000_t08_M69/` | v9.0 | 100K generated |
| `20260201_204445_v8_9_4_n100000_t08_M69/` | v8.9.4 | 100K generated |
