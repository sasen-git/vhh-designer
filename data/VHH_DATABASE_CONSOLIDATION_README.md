# VHH Database Consolidation History

## Overview

This document explains the evolution of VHH/nanobody databases in the KA-Search project, from raw data sources to the final consolidated database.

---

## Database Files (Chronological Order)

### 1. `camel_vhh_test_db.npz` (~10K sequences)
**Purpose:** Initial test run to validate the pipeline

**Source:** Subset of OAS camel heavy-chain data

**How it was built:**
```bash
python process_camel_vhh_pipeline.py \
    --input-dir /home/sasenefrem/KA-Search/extracted/oas-paper/ \
    --output-prefix camel_vhh_test \
    --batch-size 500 \
    --max-sequences 10000
```

**Key validation:**
- Confirmed ANARCI success rate: 100%
- V-domain lengths: 86-136 aa (correct range)
- CDR3 lengths: 3-30 aa (proper extraction)

---

### 2. `camel_vhh_clean_db.npz` (~1.36M sequences)
**Purpose:** Full OAS camel dataset with proper nucleotide→AA translation

**Source:** OAS Bulk Download - Camel Heavy Chain data
- `SRR3544217_Heavy_Bulk.csv` (307,632 rows)
- `SRR3544218_Heavy_Bulk.csv` (280,073 rows)
- `SRR3544219_Heavy_Bulk.csv` (235,315 rows)
- `SRR3544220_Heavy_Bulk.csv` (235,620 rows)
- `SRR3544221_Heavy_Bulk.csv` (273,732 rows)
- `SRR3544222_Heavy_Bulk.csv` (268,598 rows)
- `SRR3544217_Heavy_IGHM.csv` (130 rows)

**The Problem This Solved:**
The original OAS `sequence_alignment_aa` column was **truncated** to fit a fixed-width alignment (positions 16-190 only). This caused:
- V-domains appearing as 14-75 aa instead of 100-130 aa
- Missing or broken CDR3 regions
- Incomplete sequences

**The Solution:**
1. Extract raw nucleotide sequences from `sequence` column (300-450+ nt)
2. Translate NT → AA using Biopython (fresh translation)
3. Run ANARCI with IMGT numbering scheme
4. Extract complete V-domain and CDR1/2/3

**How it was built:**
```bash
python process_camel_vhh_pipeline.py \
    --input-dir /home/sasenefrem/KA-Search/extracted/oas-paper/ \
    --output-prefix camel_vhh_clean \
    --batch-size 5000
```

**Result:** 1,358,722 sequences with proper full-length V-domains

**Columns:**
| Column | Description |
|--------|-------------|
| seq_id | Unique identifier (camel_vhh_0000000, ...) |
| nt_sequence | Original nucleotide sequence |
| aa_v_full | Full V-domain amino acid sequence |
| cdr1 | CDR1 (IMGT 27-38) |
| cdr2 | CDR2 (IMGT 56-65) |
| cdr3 | CDR3 (IMGT 105-117) |
| len_v | V-domain length |
| len_cdr3 | CDR3 length |
| anarci_status | ok/fail |
| v_call, d_call, j_call | Germline assignments |
| _source_file | Original CSV file |

---

### 3. `VHH_db_unified.npz` (~1.36M sequences)
**Purpose:** First merge of OAS camel + annotated sources

**Sources merged:**
1. `camel_vhh_clean_db.npz` (1,354,538 - OAS Camel)
2. SAbDab nanobody structures (~1,500)
3. GenBank entries (~1,200)
4. TheraSAbDab therapeutics (~27)

**How it was built:**
```bash
python process_vhh_with_antigens.py \
    --vhh-csv vhh_sequences.csv \
    --oas-csv camel_vhh_clean_anarci_renumbered.csv \
    --output-prefix VHH_db_unified
```

**New columns added:**
| Column | Description |
|--------|-------------|
| source | OAS_Camel, SAbDab_PDB, GenBank, TheraSAbDab |
| targets | Antigen/binding partner (where available) |

---

### 4. `VHH_db_unified_v2.npz` (~1.36M sequences)
**Purpose:** Added SAbDab PDB structural data with deduplication

**What was added:**
- Extracted sequences directly from IMGT-numbered PDB files
- 1,681 unique VHH sequences from crystal structures
- Rich antigen annotations from SAbDab metadata

**How it was built:**
```bash
python process_sabdab_pdb.py \
    --tsv sabdab_nano_summary_all.tsv \
    --pdb-dir extracted/oas-paper/all_nano_structures \
    --existing-csv VHH_db_unified.csv \
    --output-prefix VHH_db_unified_v2
```

**Deduplication logic:**
When duplicate sequences found, keep the entry with:
1. Target/antigen info (highest priority)
2. More metadata (patent ID, PDB ID, etc.)
3. First occurrence (if equal)

**Result after deduplication:** 1,358,080 unique sequences

---

### 5. `VHH_db_final.npz` (~12M sequences) ⭐ PRODUCTION
**Purpose:** Complete consolidated database with all sources

**Sources merged:**
| Source | Sequences | Description |
|--------|-----------|-------------|
| INDI_NGS | 10,675,894 | Camelid repertoires (llama, alpaca, camel) |
| OAS_Camel | 1,354,537 | Original camel data |
| INDI_patent | 14,356 | Patent sequences with targets |
| SAbDab_PDB | 1,559 | Crystal structures |
| INDI_manual | 1,180 | Manually curated |
| INDI_genbank | 1,128 | GenBank submissions |
| GenBank | 433 | From earlier merge |
| INDI_structure | 429 | PDB structures |
| SAbDab | 309 | Structural database |
| TheraSAbDab | 15 | Therapeutic nanobodies |
| **TOTAL** | **12,049,840** | After deduplication |

**How it was built:**
```bash
python process_indi_merge_final.py \
    --indi-dir /home/sasenefrem/KA-Search/extracted/oas-paper/INDI \
    --existing-csv VHH_db_unified_v2.csv \
    --output-prefix VHH_db_final
```

**INDI data processing:**
- `INDI_patent`: Extracted biomolecule targets + patent IDs/titles
- `INDI_structure`: Extracted PDB info, organism, chain titles
- `INDI_manual`: Extracted curated titles/URLs
- `INDI_genbank`: Sequence only (no metadata)
- `ngs_sequence.tsv`: 10.7M NGS repertoire sequences

**Full column set:**
| Column | Description |
|--------|-------------|
| id | Unique identifier |
| source | Data source (INDI_patent, OAS_Camel, etc.) |
| aa_v_full | Full V-domain sequence |
| cdr1, cdr2, cdr3 | CDR sequences |
| len_v, len_cdr3 | Lengths |
| targets | Antigen/binding partner |
| patent_id | Patent number (INDI_patent only) |
| patent_title | Patent title |
| organism | Source organism |
| pdb_id | PDB ID (structures only) |
| reference | Literature reference |

---

## Database Lineage Diagram

```
RAW DATA SOURCES
================

OAS Bulk Download (Nucleotides)          INDI Database              SAbDab PDB
├── SRR3544217_Heavy_Bulk.csv            ├── INDI_patent (14K)      ├── ~2000 structures
├── SRR3544218_Heavy_Bulk.csv            ├── INDI_structure (500)   └── Antigen metadata
├── SRR3544219_Heavy_Bulk.csv            ├── INDI_manual (1.2K)
├── SRR3544220_Heavy_Bulk.csv            ├── INDI_genbank (1.8K)
├── SRR3544221_Heavy_Bulk.csv            └── INDI_NGS (10.7M)
├── SRR3544222_Heavy_Bulk.csv
└── Total: ~1.6M nucleotide sequences
         │
         │ PROBLEM: OAS amino acid column was truncated!
         │ SOLUTION: Translate from nucleotides + ANARCI
         ▼
┌─────────────────────────────────────┐
│  process_camel_vhh_pipeline.py      │
│  - Translate NT → AA (Biopython)    │
│  - ANARCI IMGT numbering            │
│  - Extract full V-domain + CDRs     │
└─────────────────────────────────────┘
         │
         ▼
┌─────────────────────────────────────┐
│  camel_vhh_clean_db.npz             │  ← First clean database
│  1,358,722 sequences                │
│  Full-length V-domains (100-130 aa) │
└─────────────────────────────────────┘
         │
         │ + SAbDab, GenBank, TheraSAbDab
         ▼
┌─────────────────────────────────────┐
│  VHH_db_unified.npz                 │
│  ~1.36M sequences                   │
│  + Source tracking                  │
│  + Target annotations               │
└─────────────────────────────────────┘
         │
         │ + SAbDab PDB extraction + dedup
         ▼
┌─────────────────────────────────────┐
│  VHH_db_unified_v2.npz              │
│  1,358,080 unique sequences         │
│  + Rich structural annotations      │
└─────────────────────────────────────┘
         │
         │ + INDI (patents, NGS, manual, genbank)
         ▼
┌─────────────────────────────────────┐
│  VHH_db_final.npz                   │  ← PRODUCTION DATABASE
│  12,049,840 unique sequences        │
│  - OAS_Camel: 1.35M                 │
│  - INDI_NGS: 10.7M                  │
│  - Annotated: ~19K (with targets)   │
└─────────────────────────────────────┘
         │
         │ Sharded for KA-Search
         ▼
┌─────────────────────────────────────┐
│  VHH_shards/                        │
│  ├── vhh_annotated.npz (19K)        │  ← Patents, structures, curated
│  ├── vhh_oas_camel.npz (1.35M)      │
│  └── vhh_indi_ngs_*.npz (10.7M)     │  ← Split into 2M chunks
└─────────────────────────────────────┘
```

---

## Key Scripts Used

| Script | Purpose |
|--------|---------|
| `process_camel_vhh_pipeline.py` | NT→AA translation + ANARCI numbering |
| `process_vhh_with_antigens.py` | Merge annotated sources |
| `process_sabdab_pdb.py` | Extract sequences from PDB files |
| `process_indi_merge_final.py` | Process INDI + final merge |
| `shard_database.py` | Split into shards for searching |

---

## Recommended Database for Different Use Cases

| Use Case | Database | Why |
|----------|----------|-----|
| **Full epitope search** | `VHH_db_final.npz` | Maximum coverage (12M) |
| **Quick annotated search** | `vhh_annotated.npz` | Only sequences with targets (19K) |
| **Camel-only analysis** | `vhh_oas_camel.npz` | Single species |
| **KA-Search production** | Sharded files | Memory efficient |

---

## Data Quality Notes

### What the NT→AA translation fixed:
- **Before:** V-domains 14-75 aa (truncated)
- **After:** V-domains 100-140 aa (complete)
- **Before:** CDR3 often missing or broken
- **After:** CDR3 properly extracted (5-30 aa typical for VHH)

### Deduplication strategy:
When same sequence appears in multiple sources:
1. Keep entry with **target/antigen info** (most valuable)
2. Keep entry with **patent ID** (traceable)
3. Keep entry with **PDB ID** (structural data)
4. Keep entry with **most metadata**

### Sequences with annotations:
- **Total database:** 12,049,840
- **With target info:** ~16,096 (0.13%)
- **With patent IDs:** ~14,356
- Most sequences are repertoire data (OAS_Camel, INDI_NGS) without functional annotation

---

## File Locations (Proposed Structure)

```
KA-Search/data/databases/
├── production/
│   └── VHH_db_final.npz        # Current production (12M)
├── shards/
│   ├── vhh_annotated.npz       # Annotated only (19K)
│   ├── vhh_oas_camel.npz       # OAS camel (1.35M)
│   └── vhh_indi_ngs_*.npz      # INDI NGS shards (10.7M)
└── legacy/
    ├── camel_vhh_test_db.npz   # Initial test (10K)
    ├── camel_vhh_clean_db.npz  # First full OAS (1.36M)
    ├── VHH_db_unified.npz      # First merge (1.36M)
    └── VHH_db_unified_v2.npz   # + SAbDab PDB (1.36M)
```

---

## Summary Table

| Database | Sequences | Sources | Key Addition |
|----------|-----------|---------|--------------|
| `camel_vhh_test_db.npz` | ~10K | OAS subset | Pipeline validation |
| `camel_vhh_clean_db.npz` | 1.36M | OAS Camel | NT→AA translation fix |
| `VHH_db_unified.npz` | 1.36M | + SAbDab, GenBank | Source tracking |
| `VHH_db_unified_v2.npz` | 1.36M | + PDB extraction | Deduplication |
| `VHH_db_final.npz` | **12M** | + INDI (NGS, patents) | **Full coverage** |

---

*Document created: January 2026*
*Based on work from November 2024 - January 2026*
