# data/raw/

Original source files before processing into the production database.

## Subfolders

| Folder | Contents |
|--------|----------|
| `oas_paper/compressed/` | OAS camelid bulk download (SRR3544217-SRR3544222 CSVs, gzipped) |
| `sequences/` | Input FASTA/XLSX files |

## Key Files

| File | What it is |
|------|-----------|
| `sequences/M69.fasta` | Primary input antibody (lead VH sequence for VHH conversion) |
| `sequences/M69_HC_LC.fasta` | Heavy + light chain version |
| `sequences/INDI_patent_sequences.xlsx` | Patent-derived VHH sequences |
| `sequences/INDI_structure_sequences.xlsx` | Structure-derived VHH sequences |
| `sequences/VHH_annotated_full.xlsx` | Curated annotated VHH collection |
| `oas_paper/compressed/vhh_sequences.csv.gz` | Processed VHH sequences from OAS |
| `oas_paper/compressed/all_nano_structures.zip` | Nanobody structural data |
