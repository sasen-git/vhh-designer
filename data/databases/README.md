# data/databases/

Production databases used by the pipeline.

## Subfolders

| Folder | Contents |
|--------|----------|
| `shards/` | 8 NPZ files (~12M sequences) — primary input for analysis pipeline |
| `production/` | VHH_db_unified_v2.npz (67MB) — consolidated single-file database |
| `annotated/` | IMGT-augmented CSV versions of shards |
| `legacy/` | Older database versions (VHH_db_final.npz, camel_vhh_clean_db.npz, etc.) |

## Shards (primary)

| File | Sequences | Source |
|------|-----------|--------|
| vhh_annotated.npz | 19,409 | Patents, structures, curated |
| vhh_oas_camel.npz | 1,354,537 | OAS Camel repertoire |
| vhh_indi_ngs_000.npz | 2,000,000 | INDI NGS batch 0 |
| vhh_indi_ngs_001.npz | 2,000,000 | INDI NGS batch 1 |
| vhh_indi_ngs_002.npz | 2,000,000 | INDI NGS batch 2 |
| vhh_indi_ngs_003.npz | 2,000,000 | INDI NGS batch 3 |
| vhh_indi_ngs_004.npz | 2,000,000 | INDI NGS batch 4 |
| vhh_indi_ngs_005.npz | 675,894 | INDI NGS batch 5 (remainder) |
| **Total** | **~12,049,840** | |

## Annotated CSV Versions

| Folder | Description |
|--------|-------------|
| `vhh_full_annotated_v4/` | Base annotation (AntPack numbering + family classification) |
| `vhh_full_annotated_v4_dedup/` | Deduplicated version + dedup temp buckets |
| `vhh_full_annotated_v4_dedup_indiv_imgtfull/` | Full IMGT position columns per shard — **primary input for unified analysis** |
