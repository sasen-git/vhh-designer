# Script Documentation Log

Generated: 2026-02-03 10:56:32
Directory: `/home/sasenefrem/KA-Search/`
Total scripts: 144

Folders: 28

- `(root)/` (7 scripts)
- `active/analysis/` (62 scripts)
- `active/analysis/annotate_all_shards/` (6 scripts)
- `active/analysis/vhh_compensation/` (3 scripts)
- `active/analysis/vhh_designer/` (6 scripts)
- `active/analysis/vhh_epistasis/` (3 scripts)
- `active/analysis/vhh_epistasis_to_imgt/` (3 scripts)
- `active/analysis/vhh_naturalness_analyzer/` (2 scripts)
- `active/database/` (1 scripts)
- `active/utilities/` (4 scripts)
- `archive/correlation_analysis/` (3 scripts)
- `archive/database_builders/` (4 scripts)
- `archive/epistasis_pipeline/v1_20251201/` (1 scripts)
- `archive/epistasis_pipeline/v2_20251205/` (1 scripts)
- `archive/full_analysis/` (5 scripts)
- `archive/naturalness_analyzer/` (2 scripts)
- `archive/naturalness_analyzer/v1_20251210/` (1 scripts)
- `archive/naturalness_analyzer/v2_20251212/` (1 scripts)
- `archive/npz_scanner/v2_20251103/` (1 scripts)
- `archive/npz_scanner/v3_20251103/` (1 scripts)
- `archive/npz_scanner/v5_20251104/` (1 scripts)
- `archive/npz_scanner/v6_20251105/` (1 scripts)
- `archive/npz_scanner/v7_20251106/` (1 scripts)
- `archive/npz_scanner/v8_20251107/` (1 scripts)
- `archive/npz_scanner/v9_20251108/` (1 scripts)
- `archive/one_off/` (11 scripts)
- `archive/visualizations/` (2 scripts)
- `tools/` (9 scripts)

---

## Summary by Folder

### 📁 `(root)/`

| Script | Lines | Size | Functions | Classes | Modified |
|--------|-------|------|-----------|---------|----------|
| align_vs_lead_clear3_antpack_legend_v7.py | 1032 | 38.2KB | 29 | 0 | 2026-01-31 23:21 |
| align_vs_lead_clear3_antpack_legend_v7_imgt_v2.py | 1176 | 43.5KB | 32 | 0 | 2026-02-01 00:08 |
| align_vs_lead_clear3_antpack_legend_v8.py | 945 | 33.7KB | 31 | 0 | 2026-01-31 23:11 |
| align_vs_lead_clear3_antpack_legend_v8_imgt.py | 1165 | 43.0KB | 32 | 0 | 2026-01-31 23:26 |
| csv2msa_antpack_v5.py | 536 | 19.7KB | 16 | 0 | 2026-01-31 23:55 |
| csv2msa_antpack_v5_imgt_v2.py | 536 | 19.7KB | 16 | 0 | 2026-01-31 23:54 |
| csv2msa_antpack_v6.py | 225 | 8.6KB | 11 | 0 | 2026-01-31 22:41 |

### 📁 `active/analysis/`

| Script | Lines | Size | Functions | Classes | Modified |
|--------|-------|------|-----------|---------|----------|
| analyze_cdr_framework_correlations_imgt_csv_v3.py | 664 | 24.7KB | 14 | 0 | 2026-01-09 20:42 |
| annotate_all_shards_v8.py | 651 | 22.4KB | 14 | 0 | 2026-01-07 15:36 |
| augment_imgt_positions_from_csvs_v1.py | 345 | 12.5KB | 11 | 0 | 2026-01-10 00:53 |
| augment_imgt_positions_from_csvs_v2.py | 208 | 6.8KB | 7 | 0 | 2026-01-10 02:50 |
| charge_balance_vhh_windows.py | 398 | 13.9KB | 17 | 1 | 2026-01-09 20:48 |
| convert_epistasis_to_designer.py | 603 | 20.8KB | 9 | 0 | 2026-01-03 19:15 |
| extract_script_headers.py | 379 | 12.4KB | 10 | 0 | 2026-02-03 10:52 |
| patches.py | 266 | 9.2KB | 16 | 2 | 2026-01-16 21:13 |
| run_with_positions.py | 0 | 0.0KB | 0 | 0 | 2026-01-13 13:05 |
| vhh_analysis_unified_v7.1.py | 1180 | 44.5KB | 33 | 6 | 2026-01-12 15:55 |
| vhh_analysis_unified_v7.2.py | 1277 | 49.4KB | 34 | 6 | 2026-01-12 16:17 |
| vhh_analysis_unified_v7.3.py | 1327 | 50.4KB | 38 | 7 | 2026-01-12 16:39 |
| vhh_analysis_unified_v7.4.py | 1414 | 54.4KB | 40 | 7 | 2026-01-13 13:12 |
| vhh_analysis_unified_v7.5.py | 1459 | 56.9KB | 40 | 7 | 2026-01-13 20:04 |
| vhh_analysis_unified_v7.6.py | 1517 | 59.2KB | 42 | 7 | 2026-01-13 21:05 |
| vhh_analysis_unified_v7.7.py | 1520 | 59.5KB | 42 | 7 | 2026-01-13 22:05 |
| vhh_analysis_unified_v7.py | 1058 | 38.7KB | 31 | 5 | 2026-01-13 19:19 |
| vhh_compensation_imgt_v5.py | 432 | 16.4KB | 13 | 1 | 2026-01-06 14:05 |
| vhh_compensation_imgt_v6_csvcols.py | 318 | 12.4KB | 14 | 1 | 2026-01-09 23:55 |
| vhh_designer_v7_1.py | 1721 | 68.6KB | 42 | 10 | 2026-01-05 15:31 |
| vhh_designer_v7_10.py | 3231 | 124.3KB | 67 | 11 | 2026-01-18 23:21 |
| vhh_designer_v7_11.py | 3240 | 124.8KB | 67 | 11 | 2026-01-19 00:27 |
| vhh_designer_v7_12.py | 3247 | 125.0KB | 67 | 11 | 2026-01-19 18:12 |
| vhh_designer_v7_13.py | 3539 | 140.3KB | 69 | 11 | 2026-01-20 13:58 |
| vhh_designer_v7_14.1.py | 3732 | 149.8KB | 71 | 11 | 2026-01-20 17:53 |
| vhh_designer_v7_14.py | 3654 | 146.4KB | 70 | 11 | 2026-01-20 15:46 |
| vhh_designer_v7_15.py | 3831 | 153.8KB | 71 | 11 | 2026-01-20 18:45 |
| vhh_designer_v7_16.1.py | 3978 | 158.8KB | 75 | 11 | 2026-01-20 22:33 |
| vhh_designer_v7_16.2.py | 4003 | 160.3KB | 75 | 11 | 2026-01-21 14:10 |
| vhh_designer_v7_16.py | 3978 | 158.8KB | 75 | 11 | 2026-01-20 21:41 |
| vhh_designer_v7_2.1.py | 1597 | 62.1KB | 27 | 7 | 2026-01-14 13:41 |
| vhh_designer_v7_2.py | 1122 | 40.4KB | 26 | 7 | 2026-01-14 11:11 |
| vhh_designer_v7_4_esmfold.py | 1937 | 72.2KB | 42 | 9 | 2026-01-14 16:24 |
| vhh_designer_v7_5.1.py | 2001 | 74.1KB | 41 | 10 | 2026-01-14 17:27 |
| vhh_designer_v7_5_esm.py | 1979 | 73.0KB | 41 | 10 | 2026-01-14 17:01 |
| vhh_designer_v7_6_esm.py | 2083 | 77.1KB | 45 | 10 | 2026-01-16 21:30 |
| vhh_designer_v7_7_vhhonly_debug.py | 2491 | 94.9KB | 59 | 11 | 2026-01-18 21:30 |
| vhh_designer_v7_8.py | 2891 | 109.4KB | 67 | 11 | 2026-01-18 21:41 |
| vhh_designer_v7_9.py | 2915 | 110.1KB | 67 | 11 | 2026-01-18 22:20 |
| vhh_designer_v7_91.py | 2961 | 111.9KB | 67 | 11 | 2026-01-18 22:36 |
| vhh_designer_v8_0.py | 4248 | 171.0KB | 80 | 11 | 2026-01-21 18:45 |
| vhh_designer_v8_1.py | 4309 | 173.8KB | 80 | 11 | 2026-01-21 19:16 |
| vhh_designer_v8_2.py | 4420 | 179.2KB | 80 | 11 | 2026-01-21 20:05 |
| vhh_designer_v8_3.py | 4669 | 190.2KB | 80 | 11 | 2026-01-21 22:09 |
| vhh_designer_v8_3_fixed.py | 4688 | 191.3KB | 80 | 11 | 2026-01-21 23:07 |
| vhh_designer_v8_4.py | 4900 | 200.4KB | 84 | 12 | 2026-01-26 17:44 |
| vhh_designer_v8_5.py | 5355 | 218.0KB | 86 | 12 | 2026-01-26 22:20 |
| vhh_designer_v8_6.py | 5863 | 236.7KB | 98 | 12 | 2026-01-27 19:54 |
| vhh_designer_v8_6_patched.py | 5606 | 224.8KB | 98 | 12 | 2026-01-27 21:11 |
| vhh_designer_v8_8.py | 5854 | 236.1KB | 98 | 12 | 2026-01-27 22:52 |
| vhh_designer_v8_8_1.py | 5966 | 241.6KB | 98 | 12 | 2026-01-28 12:22 |
| vhh_designer_v8_8_2.py | 6193 | 252.8KB | 98 | 12 | 2026-01-28 19:06 |
| vhh_designer_v8_9.py | 6715 | 274.9KB | 100 | 12 | 2026-01-30 00:15 |
| vhh_designer_v8_9_1.py | 7144 | 295.0KB | 101 | 12 | 2026-01-30 16:11 |
| vhh_designer_v8_9_2.py | 7270 | 300.5KB | 101 | 12 | 2026-01-31 00:00 |
| vhh_designer_v8_9_3.py | 7556 | 312.3KB | 103 | 12 | 2026-02-01 00:41 |
| vhh_designer_v8_9_4.py | 7563 | 312.7KB | 103 | 12 | 2026-02-01 15:51 |
| vhh_designer_v8_9_4_1.py | 7599 | 314.6KB | 103 | 12 | 2026-02-01 20:02 |
| vhh_designer_v9_0.py | 7652 | 317.1KB | 103 | 12 | 2026-02-01 21:42 |
| vhh_designer_v9_1.py | 7664 | 317.9KB | 103 | 12 | 2026-02-01 23:08 |
| vhh_epistasis_v3_imgt.py | 1608 | 61.3KB | 37 | 7 | 2026-01-02 23:07 |
| vhh_naturalness_analyzer_v4.py | 995 | 36.2KB | 17 | 5 | 2026-01-02 02:08 |

### 📁 `active/analysis/annotate_all_shards/`

| Script | Lines | Size | Functions | Classes | Modified |
|--------|-------|------|-----------|---------|----------|
| annotate_all_shards.py | 829 | 26.9KB | 14 | 0 | 2026-01-04 14:47 |
| annotate_all_shards_v2.py | 869 | 29.3KB | 14 | 0 | 2026-01-06 14:46 |
| annotate_all_shards_v4.py | 510 | 17.4KB | 12 | 0 | 2026-01-06 16:01 |
| annotate_all_shards_v5.py | 586 | 19.7KB | 14 | 0 | 2026-01-06 16:56 |
| annotate_all_shards_v6.py | 618 | 20.9KB | 14 | 0 | 2026-01-06 21:44 |
| annotate_all_shards_v7.py | 643 | 22.1KB | 14 | 0 | 2026-01-06 22:33 |

### 📁 `active/analysis/vhh_compensation/`

| Script | Lines | Size | Functions | Classes | Modified |
|--------|-------|------|-----------|---------|----------|
| vhh_compensation_imgt.py | 476 | 16.0KB | 10 | 1 | 2026-01-05 20:10 |
| vhh_compensation_imgt_v2.py | 391 | 13.9KB | 11 | 1 | 2026-01-05 20:29 |
| vhh_compensation_imgt_v4.py | 396 | 14.5KB | 12 | 1 | 2026-01-06 13:21 |

### 📁 `active/analysis/vhh_designer/`

| Script | Lines | Size | Functions | Classes | Modified |
|--------|-------|------|-----------|---------|----------|
| vhh_designer_v2.py | 540 | 25.9KB | 21 | 6 | 2026-01-02 21:09 |
| vhh_designer_v3.py | 1160 | 44.6KB | 33 | 8 | 2026-01-03 21:53 |
| vhh_designer_v4.py | 1318 | 51.5KB | 33 | 8 | 2026-01-03 22:09 |
| vhh_designer_v5.py | 1559 | 61.5KB | 36 | 9 | 2026-01-04 13:17 |
| vhh_designer_v5_1.py | 1561 | 61.9KB | 36 | 9 | 2026-01-04 14:25 |
| vhh_designer_v7.py | 1721 | 68.5KB | 42 | 10 | 2026-01-05 13:30 |

### 📁 `active/analysis/vhh_epistasis/`

| Script | Lines | Size | Functions | Classes | Modified |
|--------|-------|------|-----------|---------|----------|
| vhh_epistasis_overnight_final.py | 1193 | 45.2KB | 31 | 7 | 2025-12-06 23:02 |
| vhh_epistasis_overnight_v2.py | 1262 | 47.9KB | 31 | 7 | 2025-12-07 23:34 |
| visualize_epistasis.py | 345 | 11.9KB | 0 | 0 | 2025-12-07 13:16 |

### 📁 `active/analysis/vhh_epistasis_to_imgt/`

| Script | Lines | Size | Functions | Classes | Modified |
|--------|-------|------|-----------|---------|----------|
| vhh_epistasis_imgt_ligh_v3t.py | 867 | 30.1KB | 21 | 4 | 2026-01-05 00:03 |
| vhh_epistasis_imgt_light.py | 789 | 26.9KB | 21 | 4 | 2026-01-04 14:19 |
| vhh_epistasis_imgt_light_v2.py | 867 | 30.1KB | 21 | 4 | 2026-01-04 23:17 |

### 📁 `active/analysis/vhh_naturalness_analyzer/`

| Script | Lines | Size | Functions | Classes | Modified |
|--------|-------|------|-----------|---------|----------|
| vhh_foldability_analyzer.py | 915 | 33.1KB | 17 | 5 | 2025-12-10 21:58 |
| vhh_naturalness_analyzer_v3.py | 1036 | 37.9KB | 17 | 5 | 2026-01-02 01:03 |

### 📁 `active/database/`

| Script | Lines | Size | Functions | Classes | Modified |
|--------|-------|------|-----------|---------|----------|
| npz_fullscan_v6_integrated.py | 1060 | 39.5KB | 20 | 1 | 2026-01-02 00:57 |

### 📁 `active/utilities/`

| Script | Lines | Size | Functions | Classes | Modified |
|--------|-------|------|-----------|---------|----------|
| dna_translator.py | 312 | 11.3KB | 8 | 0 | 2025-12-08 15:56 |
| pathlib_list_v2.py | 60 | 1.8KB | 1 | 0 | 2026-01-02 19:21 |
| pull_cdrs.py | 530 | 15.8KB | 18 | 1 | 2026-01-02 01:53 |
| scfv_anarci.py | 100 | 3.9KB | 4 | 0 | 2026-01-18 14:53 |

### 📁 `archive/correlation_analysis/`

| Script | Lines | Size | Functions | Classes | Modified |
|--------|-------|------|-----------|---------|----------|
| analyze_cdr_framework_advanced.py | 738 | 27.9KB | 9 | 0 | 2025-12-01 21:14 |
| analyze_cdr_framework_correlations.py | 490 | 17.8KB | 6 | 0 | 2025-12-01 21:14 |
| vhh_global_compensation_analysis.py | 659 | 24.0KB | 19 | 4 | 2025-12-06 21:40 |

### 📁 `archive/database_builders/`

| Script | Lines | Size | Functions | Classes | Modified |
|--------|-------|------|-----------|---------|----------|
| build_camel_vhh_db.py | 715 | 23.2KB | 16 | 0 | 2025-11-25 20:49 |
| build_camel_vhh_db_one_step.py | 787 | 26.4KB | 15 | 0 | 2025-11-25 22:43 |
| build_models_only.py | 322 | 9.4KB | 9 | 1 | 2025-12-05 10:08 |
| build_pfr_cdr_models.py | 565 | 19.1KB | 9 | 0 | 2025-12-01 22:04 |

### 📁 `archive/epistasis_pipeline/v1_20251201/`

| Script | Lines | Size | Functions | Classes | Modified |
|--------|-------|------|-----------|---------|----------|
| vhh_full_epistasis_overnight_v1.py | 1172 | 44.3KB | 30 | 7 | 2025-12-06 22:17 |

### 📁 `archive/epistasis_pipeline/v2_20251205/`

| Script | Lines | Size | Functions | Classes | Modified |
|--------|-------|------|-----------|---------|----------|
| vhh_epistasis_overnight_v2.py | 1262 | 47.9KB | 31 | 7 | 2025-12-07 23:34 |

### 📁 `archive/full_analysis/`

| Script | Lines | Size | Functions | Classes | Modified |
|--------|-------|------|-----------|---------|----------|
| run_full_analysis_fast2_v3.py | 920 | 33.8KB | 14 | 0 | 2025-12-01 22:12 |
| run_full_analysis_fast_v2.py | 878 | 32.0KB | 14 | 0 | 2025-12-01 21:56 |
| run_full_analysis_lowmem_v4.py | 740 | 26.5KB | 18 | 3 | 2025-12-04 23:29 |
| run_full_analysis_lowmem_v5.py | 797 | 29.5KB | 15 | 2 | 2025-12-05 18:36 |
| run_full_analysis_v1.py | 974 | 37.1KB | 13 | 0 | 2025-12-01 21:13 |

### 📁 `archive/naturalness_analyzer/`

| Script | Lines | Size | Functions | Classes | Modified |
|--------|-------|------|-----------|---------|----------|
| score_framework.py | 460 | 15.5KB | 8 | 0 | 2025-12-01 21:13 |
| vhh_framework_optimizer.py | 503 | 17.5KB | 10 | 1 | 2025-12-06 20:56 |

### 📁 `archive/naturalness_analyzer/v1_20251210/`

| Script | Lines | Size | Functions | Classes | Modified |
|--------|-------|------|-----------|---------|----------|
| vhh_naturalness_analyzer.py | 844 | 30.3KB | 17 | 5 | 2026-01-02 01:18 |

### 📁 `archive/naturalness_analyzer/v2_20251212/`

| Script | Lines | Size | Functions | Classes | Modified |
|--------|-------|------|-----------|---------|----------|
| vhh_naturalness_analyzer_v2.py | 915 | 33.1KB | 17 | 5 | 2026-01-02 01:17 |

### 📁 `archive/npz_scanner/v2_20251103/`

| Script | Lines | Size | Functions | Classes | Modified |
|--------|-------|------|-----------|---------|----------|
| npz_fullscan_v2.py | 748 | 28.9KB | 16 | 1 | 2025-11-07 23:07 |

### 📁 `archive/npz_scanner/v3_20251103/`

| Script | Lines | Size | Functions | Classes | Modified |
|--------|-------|------|-----------|---------|----------|
| npz_fullscan_v3.py | 502 | 17.0KB | 14 | 1 | 2025-11-11 13:57 |

### 📁 `archive/npz_scanner/v5_20251104/`

| Script | Lines | Size | Functions | Classes | Modified |
|--------|-------|------|-----------|---------|----------|
| npz_fullscan_v5.py | 1075 | 38.7KB | 22 | 2 | 2025-11-11 22:41 |

### 📁 `archive/npz_scanner/v6_20251105/`

| Script | Lines | Size | Functions | Classes | Modified |
|--------|-------|------|-----------|---------|----------|
| npz_fullscan_v6_interactive_v6.py | 2365 | 94.7KB | 23 | 1 | 2025-11-23 15:35 |

### 📁 `archive/npz_scanner/v7_20251106/`

| Script | Lines | Size | Functions | Classes | Modified |
|--------|-------|------|-----------|---------|----------|
| npz_fullscan_v7_interactive_v7.py | 2574 | 102.9KB | 26 | 1 | 2025-11-24 22:01 |

### 📁 `archive/npz_scanner/v8_20251107/`

| Script | Lines | Size | Functions | Classes | Modified |
|--------|-------|------|-----------|---------|----------|
| npz_fullscan_v8_vhh_v8.py | 3117 | 120.0KB | 50 | 1 | 2025-11-29 23:18 |

### 📁 `archive/npz_scanner/v9_20251108/`

| Script | Lines | Size | Functions | Classes | Modified |
|--------|-------|------|-----------|---------|----------|
| npz_fullscan_v9_vhh_v9.py | 3160 | 121.9KB | 51 | 1 | 2025-11-30 14:07 |

### 📁 `archive/one_off/`

| Script | Lines | Size | Functions | Classes | Modified |
|--------|-------|------|-----------|---------|----------|
| debug_npz_format.py | 59 | 1.9KB | 1 | 0 | 2025-12-06 21:48 |
| diagnose_npz.py | 29 | 0.9KB | 0 | 0 | 2025-12-07 23:24 |
| diagnose_oas_aligned.py | 249 | 8.7KB | 4 | 0 | 2025-11-24 22:22 |
| diagnose_oas_data.py | 339 | 11.7KB | 4 | 0 | 2025-11-25 20:49 |
| download_oas_camel.py | 214 | 6.2KB | 3 | 0 | 2025-11-25 20:49 |
| interactive_comprehensive_cdr_analysis.py | 923 | 32.9KB | 18 | 1 | 2025-11-11 19:56 |
| process_camel_vhh_pipeline.py | 787 | 26.4KB | 15 | 0 | 2025-11-25 22:43 |
| process_indi_merge_final.py | 452 | 17.1KB | 8 | 0 | 2025-11-28 15:27 |
| process_sabdab_pdb.py | 435 | 15.1KB | 9 | 0 | 2025-11-27 00:36 |
| process_vhh_with_antigens.py | 486 | 15.4KB | 12 | 0 | 2025-11-26 23:45 |
| shard_database.py | 246 | 8.1KB | 6 | 0 | 2025-11-28 23:58 |

### 📁 `archive/visualizations/`

| Script | Lines | Size | Functions | Classes | Modified |
|--------|-------|------|-----------|---------|----------|
| visualize_correlations.py | 358 | 14.1KB | 7 | 0 | 2025-12-01 21:14 |
| visualize_pfr_models.py | 236 | 7.9KB | 7 | 0 | 2025-12-01 21:13 |

### 📁 `tools/`

| Script | Lines | Size | Functions | Classes | Modified |
|--------|-------|------|-----------|---------|----------|
| dedup_annotated_imgt.py | 131 | 5.1KB | 4 | 0 | 2026-01-09 18:17 |
| dedupe_csv_headers.py | 41 | 1.2KB | 2 | 0 | 2026-01-07 22:25 |
| diagnose_unknown_files.py | 275 | 9.6KB | 5 | 0 | 2025-12-31 15:37 |
| kasearch_paths.py | 434 | 11.9KB | 31 | 1 | 2025-12-31 21:38 |
| kasearch_reorganize.py | 894 | 30.4KB | 16 | 0 | 2025-12-31 16:45 |
| maintain.py | 749 | 25.9KB | 17 | 2 | 2026-01-01 18:31 |
| pathlib_list_v2.py | 60 | 1.8KB | 1 | 0 | 2026-01-02 01:23 |
| qc_imgt_coverage_fast.py | 105 | 2.9KB | 4 | 0 | 2026-01-09 18:53 |
| summarize_sweep.py | 37 | 1.1KB | 3 | 0 | 2026-01-13 16:48 |

---

## Script Details

## 📁 `(root)/`

### `align_vs_lead_clear3_antpack_legend_v7.py`
**Version:** 7.
**Path:** `align_vs_lead_clear3_antpack_legend_v7.py`
**Stats:** 1032 lines | 38.2KB | 29 functions | 0 classes
**Modified:** 2026-01-31 23:21

**Description:**
```
align_vs_lead_clear3_antpack_legend_v7.py

v7 Changes:
- Switched from Clustal Omega to MAFFT --localpair for better CDR alignment
- MAFFT's local pairwise alignment keeps motifs together (e.g., DYD aligns with DYD)
- Falls back to Clustal Omega if MAFFT is not installed

Install MAFFT: sudo apt install mafft

Updated to handle sorting even when no CDRs are detected.
Falls back to framework-only similarity scoring.
```

---

### `align_vs_lead_clear3_antpack_legend_v7_imgt_v2.py`
**Version:** 7.
**Path:** `align_vs_lead_clear3_antpack_legend_v7_imgt_v2.py`
**Stats:** 1176 lines | 43.5KB | 32 functions | 0 classes
**Modified:** 2026-02-01 00:08

**Description:**
```
align_vs_lead_clear3_antpack_legend_v7.py

v7 Changes:
- Switched from Clustal Omega to MAFFT --localpair for better CDR alignment
- MAFFT's local pairwise alignment keeps motifs together (e.g., DYD aligns with DYD)
- Falls back to Clustal Omega if MAFFT is not installed

Install MAFFT: sudo apt install mafft

Updated to handle sorting even when no CDRs are detected.
Falls back to framework-only similarity scoring.
```

---

### `align_vs_lead_clear3_antpack_legend_v8.py`
**Version:** 8.
**Path:** `align_vs_lead_clear3_antpack_legend_v8.py`
**Stats:** 945 lines | 33.7KB | 31 functions | 0 classes
**Modified:** 2026-01-31 23:11

**Description:**
```
align_vs_lead_clear3_antpack_legend_v8.py

v8 Changes:
- Added IMGT position track below consensus
- Shows H (hallmark) and V (vernier) markers at key VHH positions
- IMGT numbering derived from AntPack

v7 Changes:
- Switched from Clustal Omega to MAFFT --localpair for better CDR alignment
- MAFFT's local pairwise alignment keeps motifs together (e.g., DYD aligns with DYD)
- Falls back to Clustal Omega if MAFFT is not installed
```

---

### `align_vs_lead_clear3_antpack_legend_v8_imgt.py`
**Version:** 7.
**Path:** `align_vs_lead_clear3_antpack_legend_v8_imgt.py`
**Stats:** 1165 lines | 43.0KB | 32 functions | 0 classes
**Modified:** 2026-01-31 23:26

**Description:**
```
align_vs_lead_clear3_antpack_legend_v7.py

v7 Changes:
- Switched from Clustal Omega to MAFFT --localpair for better CDR alignment
- MAFFT's local pairwise alignment keeps motifs together (e.g., DYD aligns with DYD)
- Falls back to Clustal Omega if MAFFT is not installed

Install MAFFT: sudo apt install mafft

Updated to handle sorting even when no CDRs are detected.
Falls back to framework-only similarity scoring.
```

---

### `csv2msa_antpack_v5.py`
**Version:** 2
**Path:** `csv2msa_antpack_v5.py`
**Stats:** 536 lines | 19.7KB | 16 functions | 0 classes
**Modified:** 2026-01-31 23:55

**Description:**
```
csv2msa_antpack_v3.py

Updated wrapper that handles sorting even when no CDRs are detected.
Falls back to framework-only similarity scoring.
```

---

### `csv2msa_antpack_v5_imgt_v2.py`
**Version:** 2
**Path:** `csv2msa_antpack_v5_imgt_v2.py`
**Stats:** 536 lines | 19.7KB | 16 functions | 0 classes
**Modified:** 2026-01-31 23:54

**Description:**
```
csv2msa_antpack_v3.py

Updated wrapper that handles sorting even when no CDRs are detected.
Falls back to framework-only similarity scoring.
```

---

### `csv2msa_antpack_v6.py`
**Version:** 2
**Path:** `csv2msa_antpack_v6.py`
**Stats:** 225 lines | 8.6KB | 11 functions | 0 classes
**Modified:** 2026-01-31 22:41

**Description:**
```
csv2msa_antpack_v6.py

v6: Uses v8 renderer with IMGT H/V position markers
```

---

## 📁 `active/analysis/`

### `analyze_cdr_framework_correlations_imgt_csv_v3.py`
**Version:** 3
**Path:** `active/analysis/analyze_cdr_framework_correlations_imgt_csv_v3.py`
**Stats:** 664 lines | 24.7KB | 14 functions | 0 classes
**Modified:** 2026-01-09 20:42

**Description:**
```
Correlation Results v3 (CSV + IMGT) - Complete Breakdown
=======================================================

Goal
- Read the *dedup* annotated CSV with IMGT columns (imgt_##) and CDR strings.
- Produce:
  (1) correlation_results_v3_summary.json  (your "Complete Breakdown")
  (2) correlation_rules_v3.json           (compensation-style rules JSON)

Why this exists
- Your older correlations script is NPZ-based and extracts FRs via substring search,
  which is not IMGT-position-safe.
- This version treats IMGT columns as the coordinate system and emits the same rule
  schema style as the compensation pipeline:
    {condition, position, suggested_aa, confidence, support, lift, baseline_prob, source}

Outputs
- correlation_results_v3_summary.json:
  {
    "meta": {...},
    "germline_counts": {...},
    "hallmark_distributions": {...},
    "family_cdr_dependent": {...},
    "family_conservation": {...},
    "family_cdr_length_stats": {...},
    "indel_correlations": {...},
    "family_triplet_rules": {...}
  }

- correlation_rules_v3.json:
  list[ {condition, position, suggested_aa, confidence, support, lift, baseline_prob, source, rule_type} ]

Notes on performance
- Full-pass summaries are cheap (family counts, hallmark distributions, length stats).
- Rule mining can explode combinatorially; this script uses deterministic sampling for rules by default.
  Set --sample-per-million 1000000 to use all rows (slow).

USAGE
  python analyze_cdr_framework_correlations_imgt_csv_v3.py     --csv data/.../vhh_full_annotated_imgt_dedup.csv     -o models/correlations/imgt_csv_v3     --chunk-size 200000     --positions fr3fr4     --conditions minimal     --min-support 5000     --min-confidence 0.70     --min-lift 1.15     --sample-per-million 200000
```

---

### `annotate_all_shards_v8.py`
**Version:** 3
**Path:** `active/analysis/annotate_all_shards_v8.py`
**Stats:** 651 lines | 22.4KB | 14 functions | 0 classes
**Modified:** 2026-01-07 15:36

**Description:**
```
VHH Database ANARCI Annotation Pipeline (streaming, low-memory) - v3
===================================================================

Key differences vs v2:
- Writes ONE CSV per shard incrementally (no giant in-memory concat)
- Processes each NPZ shard in chunks, and ANARCI in batches inside each chunk
- Emits IMGT per-position columns for FR3/FR4 by default (optional FR2)
- Produces a summary JSON + manifest of outputs

Example:
  python annotate_all_shards_v3.py     -i data/databases/shards     -o data/databases/annotated/vhh_full_annotated_final_v2     --use-anarci     -b 2000     --chunk-size 50000
```

---

### `augment_imgt_positions_from_csvs_v1.py`
**Path:** `active/analysis/augment_imgt_positions_from_csvs_v1.py`
**Stats:** 345 lines | 12.5KB | 11 functions | 0 classes
**Modified:** 2026-01-10 00:53

**Description:**
```
Augment existing annotated VHH CSV shards with full IMGT position columns.

- Input: directory of per-shard CSVs (e.g., *.dedup.csv)
- Output: same shards with added imgt_1..imgt_128 columns (integer IMGT positions only)
- Uses ANARCI to compute IMGT numbering from aa_v_full (or fallback sequence column)

Notes:
- We only emit integer IMGT positions 1..128. Insertions like 27A are ignored.
- This is typically what you want for "framework-wide" and stable positional features.
```

---

### `augment_imgt_positions_from_csvs_v2.py`
**Path:** `active/analysis/augment_imgt_positions_from_csvs_v2.py`
**Stats:** 208 lines | 6.8KB | 7 functions | 0 classes
**Modified:** 2026-01-10 02:50

**Description:**
```
(No documentation found)
```

---

### `charge_balance_vhh_windows.py`
**Version:** 4
**Path:** `active/analysis/charge_balance_vhh_windows.py`
**Stats:** 398 lines | 13.9KB | 17 functions | 1 classes
**Modified:** 2026-01-09 20:48

**Description:**
```
charge_balance_vhh_windows.py

Streaming analysis of charge "balancing" across IMGT-position windows for VHH datasets.

Works with your current schema:
- IMGT columns present: typically 66-104, 118-128, plus 42,49,50,52 (may vary)
- Computes per-window net charge, then correlation matrices:
    1) overall raw correlation
    2) overall residualized on covariates (default: total_charge + len_cdr3 + len_v)
    3) within-family pooled raw correlation (controls for family)
    4) within-family pooled residual correlation (controls for family + covariates)

Outputs:
- windows.tsv
- corr_raw.overall.tsv
- corr_resid.overall.tsv
- corr_raw.within_family_pooled.tsv
- corr_resid.within_family_pooled.tsv
- family_counts.tsv (optional but useful)

Example:
  python3 charge_balance_vhh_windows.py     --inputs data/databases/annotated/vhh_full_annotated_v4/*dedup.csv     --outdir out_charge_bal     --window-size 8     --require-valid
```

---

### `convert_epistasis_to_designer.py`
**Version:** 2
**Path:** `active/analysis/convert_epistasis_to_designer.py`
**Stats:** 603 lines | 20.8KB | 9 functions | 0 classes
**Modified:** 2026-01-03 19:15

**Description:**
```
Convert Epistasis Analysis Output to VHH Designer Format

This script converts the epistasis_v2_full.pkl output into the format expected
by vhh_designer_v2.py, creating properly structured:
- multi_position_rules (from analysis_4_higher_order_rules)
- vernier_archetypes (from analysis_2_vernier_clusters)

Usage:
    python convert_epistasis_to_designer.py         --epistasis epistasis_v2_full.pkl         --output correlation_results_v3_compensation.pkl         [--existing correlation_results_v3.pkl]  # Optional: merge with existing

Author: Claude (Anthropic)
Date: January 2025
```

---

### `extract_script_headers.py`
**Path:** `active/analysis/extract_script_headers.py`
**Stats:** 379 lines | 12.4KB | 10 functions | 0 classes
**Modified:** 2026-02-03 10:52

**Description:**
```
Extract Script Headers/Docstrings
==================================
Scans a directory for Python files and extracts the module-level
docstrings and initial comments to create a documentation log.

Usage:
    python extract_script_headers.py /path/to/scripts
    python extract_script_headers.py /path/to/scripts --output my_log.md
    python extract_script_headers.py /path/to/scripts --format csv
```

---

### `patches.py`
**Path:** `active/analysis/patches.py`
**Stats:** 266 lines | 9.2KB | 16 functions | 2 classes
**Modified:** 2026-01-16 21:13

**Description:**
```
slight modification from https://github.com/ibivu/hydrophobic_patches/
@author: jeff-lafrence, gorantlal
```

---

### `run_with_positions.py`
**Path:** `active/analysis/run_with_positions.py`
**Stats:** 0 lines | 0.0KB | 0 functions | 0 classes
**Modified:** 2026-01-13 13:05

**Description:**
```
(No documentation found)
```

---

### `vhh_analysis_unified_v7.1.py`
**Version:** 7
**Path:** `active/analysis/vhh_analysis_unified_v7.1.py`
**Stats:** 1180 lines | 44.5KB | 33 functions | 6 classes
**Modified:** 2026-01-12 15:55

**Description:**
```
VHH Analysis / Compensation Unified v7 (CSV + IMGT)
====================================================

Single-script engine that:
  - Streams an IMGT-annotated VHH CSV
  - Applies truncation and position exclusions
  - Infers family and hallmarks (if missing)
  - Extracts CDR3 features and conditions
  - Learns CDR→FR position rules (support / confidence / lift)
  - Identifies Vernier archetypes per family
  - Computes mutual information (MI) between positions on a subsample
  - Computes family-aware conservation + correlations (CorrelationAnalyzer)

OUTPUTS:
  - analysis_summary_v7.json        : Metadata, family counts, stats
  - analysis_rules_v7.json          : All compensation rules (unified schema)
  - analysis_vernier_archetypes_v7.json : Per-family vernier consensus
  - analysis_mi_pairs_v7.json       : Top MI position pairs
  - analysis_correlations_v7.json   : Per-family conservation + CDR3↔FR correlations

USAGE:
  python vhh_analysis_unified_v7.py \
      --csv data/.../vhh_annotated_imgt.csv \
      -o models/analysis/v7 \
      --positions vernier \
      --conditions minimal \
      --filter-truncated \
      --min-support 5000 \
      --min-confidence 0.7 \
      --min-lift 1.15 \
      --mi-max-seqs 2000000
```

---

### `vhh_analysis_unified_v7.2.py`
**Version:** 7
**Path:** `active/analysis/vhh_analysis_unified_v7.2.py`
**Stats:** 1277 lines | 49.4KB | 34 functions | 6 classes
**Modified:** 2026-01-12 16:17

**Description:**
```
VHH Analysis / Compensation Unified v7 (CSV + IMGT)
====================================================

Single-script engine that:
  - Streams an IMGT-annotated VHH CSV
  - Applies truncation and position exclusions
  - Infers family and hallmarks (if missing)
  - Extracts CDR3 features and conditions
  - Learns CDR→FR position rules (support / confidence / lift)
  - Identifies Vernier archetypes per family
  - Computes mutual information (MI) between positions on a subsample
  - Computes family-aware conservation + correlations (CorrelationAnalyzer)

OUTPUTS:
  - analysis_summary_v7.json        : Metadata, family counts, stats
  - analysis_rules_v7.json          : All compensation rules (unified schema)
  - analysis_vernier_archetypes_v7.json : Per-family vernier consensus
  - analysis_mi_pairs_v7.json       : Top MI position pairs
  - analysis_correlations_v7.json   : Per-family conservation + CDR3↔FR correlations

USAGE:
  python vhh_analysis_unified_v7.py \
      --csv data/.../vhh_annotated_imgt.csv \
      -o models/analysis/v7 \
      --positions vernier \
      --conditions minimal \
      --filter-truncated \
      --min-support 5000 \
      --min-confidence 0.7 \
      --min-lift 1.15 \
      --mi-max-seqs 2000000
```

---

### `vhh_analysis_unified_v7.3.py`
**Version:** 7
**Path:** `active/analysis/vhh_analysis_unified_v7.3.py`
**Stats:** 1327 lines | 50.4KB | 38 functions | 7 classes
**Modified:** 2026-01-12 16:39

**Description:**
```
VHH Analysis / Compensation Unified v7 (CSV + IMGT)
====================================================

Single-script engine that:
  - Streams an IMGT-annotated VHH CSV
  - Applies truncation and position exclusions
  - Infers family and hallmarks (if missing)
  - Extracts CDR3 features and conditions
  - Learns CDR→FR position rules (support / confidence / lift)
  - Identifies Vernier archetypes per family
  - Computes mutual information (MI) between positions on a subsample
  - Computes family-aware conservation + correlations (CorrelationAnalyzer)

OUTPUTS:
  - analysis_summary_v7.json        : Metadata, family counts, stats
  - analysis_rules_v7.json          : All compensation rules (unified schema)
  - analysis_vernier_archetypes_v7.json : Per-family vernier consensus
  - analysis_mi_pairs_v7.json       : Top MI position pairs
  - analysis_correlations_v7.json   : Per-family conservation + CDR3↔FR correlations

USAGE:
  python vhh_analysis_unified_v7.py \
      --csv data/.../vhh_annotated_imgt.csv \
      -o models/analysis/v7 \
      --positions vernier \
      --conditions minimal \
      --filter-truncated \
      --min-support 5000 \
      --min-confidence 0.7 \
      --min-lift 1.15 \
      --mi-max-seqs 2000000
```

---

### `vhh_analysis_unified_v7.4.py`
**Version:** 7.4
**Path:** `active/analysis/vhh_analysis_unified_v7.4.py`
**Stats:** 1414 lines | 54.4KB | 40 functions | 7 classes
**Modified:** 2026-01-13 13:12

**Description:**
```
VHH Analysis / Compensation Unified v7.4 (CSV + IMGT)
====================================================

Single-script engine that:
  - Streams one or more IMGT-annotated VHH CSVs
  - Applies truncation and position exclusions
  - Infers family and hallmarks (if missing)
  - Extracts CDR3 features and conditions
  - Learns CDR→FR position rules (family-conditional lift, grouped output)
  - Identifies Vernier archetypes per family
  - Tracks true Vernier clusters with pattern-level CDR3 stats
  - Computes mutual information (MI) between positions on a subsample
  - Computes family-aware conservation + correlations

OUTPUTS (in dated subfolder: v7.4_YYYY-MM-DD/):
  - analysis_summary_v7.json             : Metadata, family counts, stats
  - analysis_rules_v7.json               : All rules (grouped with suggestions list)
  - analysis_vernier_archetypes_v7.json  : Per-family vernier consensus
  - analysis_vernier_clusters_v7.json    : True vernier pattern clusters with CDR3 stats
  - analysis_mi_pairs_v7.json            : [OPTIONAL] Top MI position pairs (--enable-mi)
  - analysis_correlations_v7.json        : Per-family conservation + CDR3↔FR correlations

LEGACY PICKLE OUTPUTS (--emit-legacy-pickles):
  - analysis_compensation_legacy.pkl     : For vhh_designer_v7 --compensation
  - analysis_epistasis_legacy.pkl        : For vhh_designer_v7 --epistasis / naturalness_analyzer

USAGE:
  # Multiple CSVs with legacy pickle output (creates models/analysis/v7.4_2026-01-12/)
  python vhh_analysis_unified_v7.py \
      --csv data/*.csv \
      -o models/analysis \
      --filter-truncated \
      --emit-legacy-pickles

v7.4 CHANGES:
  - Output folder now includes version + date (v7.4_YYYY-MM-DD)
  - Added --cluster-min-positions (default 12) to control cluster explosion
  - Added --no-date-folder to write directly to --output
  - Rules output grouped: 1 record per (family, condition, position) with suggestions list
  - Vernier archetypes optimized: dedicated streaming counter
  - True Vernier clusters: VernierClusterAnalyzer with pattern-level CDR3 stats
  - Picklable defaultdict factories (no lambdas)
```

---

### `vhh_analysis_unified_v7.5.py`
**Version:** 7.5
**Path:** `active/analysis/vhh_analysis_unified_v7.5.py`
**Stats:** 1459 lines | 56.9KB | 40 functions | 7 classes
**Modified:** 2026-01-13 20:04

**Description:**
```
VHH Analysis / Compensation Unified v7.5 (CSV + IMGT)
====================================================

Single-script engine that:
  - Streams one or more IMGT-annotated VHH CSVs
  - Applies truncation and position exclusions
  - Infers family and hallmarks (optional: --infer-family)
  - Extracts CDR3 features and conditions
  - Learns CDR→FR position rules (family-conditional lift, grouped output)
  - Identifies Vernier archetypes per family
  - Tracks true Vernier clusters with pattern-level CDR3 stats
  - Computes mutual information (MI) between positions on a subsample
  - Computes family-aware conservation + correlations

OUTPUTS (in dated subfolder: v7.5_YYYY-MM-DD/):
  - run_config.json                      : Parameters, system info (written at start)
  - analysis_summary_v7.json             : Metadata, family counts, stats
  - analysis_rules_v7.json               : All rules (grouped with suggestions list)
  - analysis_vernier_archetypes_v7.json  : Per-family vernier consensus
  - analysis_vernier_clusters_v7.json    : True vernier pattern clusters with CDR3 stats
  - analysis_mi_pairs_v7.json            : [OPTIONAL] Top MI position pairs (--enable-mi)
  - analysis_correlations_v7.json        : Per-family conservation + CDR3↔FR correlations
  - analysis_compensation_legacy.pkl     : For vhh_designer_v7 (default ON)
  - analysis_epistasis_legacy.pkl        : For vhh_naturalness_analyzer (default ON)

USAGE:
  # Standard run (sensible defaults)
  python vhh_analysis_unified_v7.py --csv data/*.csv -o models/analysis

  # With family re-inference (splits coarse families)
  python vhh_analysis_unified_v7.py --csv data/*.csv -o models/analysis --infer-family

v7.5 CHANGES:
  - New defaults: filter-truncated ON, emit-legacy-pickles ON, cluster-min-positions 16
  - New defaults: min-lift 1.1, min-aa-support 250
  - Added --infer-family to force family classification from sequence
  - Checkpoint now validates params match before resuming
  - run_config.json now includes vernier_positions, target_positions, excluded_positions
  - Removed VERNIER_OVERRIDE env var (use --positions or modify VERNIER_POSITIONS constant)
```

---

### `vhh_analysis_unified_v7.6.py`
**Version:** 7.5
**Path:** `active/analysis/vhh_analysis_unified_v7.6.py`
**Stats:** 1517 lines | 59.2KB | 42 functions | 7 classes
**Modified:** 2026-01-13 21:05

**Description:**
```
VHH Analysis / Compensation Unified v7.5 (CSV + IMGT)
====================================================

Single-script engine that:
  - Streams one or more IMGT-annotated VHH CSVs
  - Applies truncation and position exclusions
  - Infers family and hallmarks (optional: --infer-family)
  - Extracts CDR3 features and conditions
  - Learns CDR→FR position rules (family-conditional lift, grouped output)
  - Identifies Vernier archetypes per family
  - Tracks true Vernier clusters with pattern-level CDR3 stats
  - Computes mutual information (MI) between positions on a subsample
  - Computes family-aware conservation + correlations

OUTPUTS (in dated subfolder: v7.5_YYYY-MM-DD/):
  - run_config.json                      : Parameters, system info (written at start)
  - analysis_summary_v7.json             : Metadata, family counts, stats
  - analysis_rules_v7.json               : All rules (grouped with suggestions list)
  - analysis_vernier_archetypes_v7.json  : Per-family vernier consensus
  - analysis_vernier_clusters_v7.json    : True vernier pattern clusters with CDR3 stats
  - analysis_mi_pairs_v7.json            : [OPTIONAL] Top MI position pairs (--enable-mi)
  - analysis_correlations_v7.json        : Per-family conservation + CDR3↔FR correlations
  - analysis_compensation_legacy.pkl     : For vhh_designer_v7 (default ON)
  - analysis_epistasis_legacy.pkl        : For vhh_naturalness_analyzer (default ON)

USAGE:
  # Standard run (sensible defaults - infers family, all FR positions)
  python vhh_analysis_unified_v7.py --csv data/*.csv -o models/analysis

  # Use CSV family column instead of inferring
  python vhh_analysis_unified_v7.py --csv data/*.csv -o models/analysis --no-infer-family

v7.5 CHANGES:
  - New defaults: filter-truncated ON, emit-legacy-pickles ON, cluster-min-positions 16
  - New defaults: min-lift 1.1, min-aa-support 250, positions=all (85 FR positions)
  - New default: infer-family ON (computes F_C2/F_C4/Y_C2/Y_C4/VH_like/VHH_W52/Other_VHH)
  - Checkpoint now validates params match before resuming
  - run_config.json now includes vernier_positions, target_positions, excluded_positions
  - Rule normalization: all rules now have suggestions=[] list
```

---

### `vhh_analysis_unified_v7.7.py`
**Version:** 7.7
**Path:** `active/analysis/vhh_analysis_unified_v7.7.py`
**Stats:** 1520 lines | 59.5KB | 42 functions | 7 classes
**Modified:** 2026-01-13 22:05

**Description:**
```
VHH Analysis / Compensation Unified v7.7 (CSV + IMGT)
====================================================

Single-script engine that:
  - Streams one or more IMGT-annotated VHH CSVs
  - Applies truncation and position exclusions
  - Infers family and hallmarks (optional: --infer-family)
  - Extracts CDR3 features and conditions
  - Learns CDR→FR position rules (family-conditional lift, grouped output)
  - Identifies Vernier archetypes per family
  - Tracks true Vernier clusters with pattern-level CDR3 stats
  - Computes mutual information (MI) between positions on a subsample
  - Computes family-aware conservation + correlations

OUTPUTS (in dated subfolder: v7.7_YYYY-MM-DD/):
  - run_config.json                      : Parameters, system info (written at start)
  - analysis_summary_v7.json             : Metadata, family counts, stats
  - analysis_rules_v7.json               : All rules (grouped with suggestions list)
  - analysis_vernier_archetypes_v7.json  : Per-family vernier consensus
  - analysis_vernier_clusters_v7.json    : True vernier pattern clusters with CDR3 stats
  - analysis_mi_pairs_v7.json            : [OPTIONAL] Top MI position pairs (--enable-mi)
  - analysis_correlations_v7.json        : Per-family conservation + CDR3↔FR correlations
  - analysis_compensation_legacy.pkl     : For vhh_designer_v7 (default ON)
  - analysis_epistasis_legacy.pkl        : For vhh_naturalness_analyzer (default ON)

USAGE:
  # Standard run (sensible defaults - infers family, all FR positions)
  python vhh_analysis_unified_v7.py --csv data/*.csv -o models/analysis

  # Use CSV family column instead of inferring
  python vhh_analysis_unified_v7.py --csv data/*.csv -o models/analysis --no-infer-family

v7.7 CHANGES:
  - Output folder now includes time: v7.7_YYYY-MM-DD_HHMMSS (no overwrites)
  - New defaults: filter-truncated ON, emit-legacy-pickles ON, cluster-min-positions 16
  - New defaults: min-lift 1.1, min-aa-support 250, positions=all (85 FR positions)
  - New default: infer-family ON (computes F_C2/F_C4/Y_C2/Y_C4/VH_like/VHH_W52/Other_VHH)
  - Checkpoint validates params match before resuming
  - Rule normalization: all rules now have suggestions=[] list
```

---

### `vhh_analysis_unified_v7.py`
**Version:** 7
**Path:** `active/analysis/vhh_analysis_unified_v7.py`
**Stats:** 1058 lines | 38.7KB | 31 functions | 5 classes
**Modified:** 2026-01-13 19:19

**Description:**
```
VHH Analysis / Compensation Unified v7 (CSV + IMGT)
====================================================

Single-script engine that:
  - Streams an IMGT-annotated VHH CSV
  - Applies truncation and position exclusions
  - Infers family and hallmarks (if missing)
  - Extracts CDR3 features and conditions
  - Learns CDR→FR position rules (support / confidence / lift)
  - Identifies Vernier archetypes per family
  - Computes mutual information (MI) between positions on a subsample
  - Computes family-aware conservation + correlations (CorrelationAnalyzer)

OUTPUTS:
  - analysis_summary_v7.json        : Metadata, family counts, stats
  - analysis_rules_v7.json          : All compensation rules (unified schema)
  - analysis_vernier_archetypes_v7.json : Per-family vernier consensus
  - analysis_mi_pairs_v7.json       : Top MI position pairs
  - analysis_correlations_v7.json   : Per-family conservation + CDR3↔FR correlations

USAGE:
  python vhh_analysis_unified_v7.py \
      --csv data/.../vhh_annotated_imgt.csv \
      -o models/analysis/v7 \
      --positions vernier \
      --conditions minimal \
      --filter-truncated \
      --min-support 5000 \
      --min-confidence 0.7 \
      --min-lift 1.15 \
      --mi-max-seqs 2000000
```

---

### `vhh_compensation_imgt_v5.py`
**Version:** 2
**Path:** `active/analysis/vhh_compensation_imgt_v5.py`
**Stats:** 432 lines | 16.4KB | 13 functions | 1 classes
**Modified:** 2026-01-06 14:05

**Description:**
```
VHH Compensation Analysis - IMGT Positions (v2)

What this script does
- Streams an annotated VHH CSV (fr1/cdr1/fr2/cdr2/fr3/cdr3/fr4 + family + hallmarks)
- Learns simple conditional residue rules: condition -> (IMGT position -> preferred AA)
- Outputs JSON rules compatible with the vhh_designer_v7 IMGT-rule schema:
    {condition, position: "IMGT##", suggested_aa, confidence, support, lift, baseline_prob, source}

Why v2
- Adds sanity checks to ensure the FR strings in the CSV are aligned to true IMGT boundaries
- Greatly reduces runtime/memory by:
    * restricting analysis to a configurable set of FR positions
    * using a smaller, higher-signal condition set by default
    * avoiding pandas iterrows()

USAGE (recommended starting point)
python vhh_compensation_imgt_v2.py   --csv data/databases/annotated/vhh_full_annotated_v3.csv   --output models/compensation/imgt_v1   --chunk-size 200000   --positions fr3fr4   --conditions minimal   --min-support 5000   --min-confidence 0.70   --min-lift 1.15
```

---

### `vhh_compensation_imgt_v6_csvcols.py`
**Path:** `active/analysis/vhh_compensation_imgt_v6_csvcols.py`
**Stats:** 318 lines | 12.4KB | 14 functions | 1 classes
**Modified:** 2026-01-09 23:55

**Description:**
```
(No documentation found)
```

---

### `vhh_designer_v7_1.py`
**Version:** 7
**Path:** `active/analysis/vhh_designer_v7_1.py`
**Stats:** 1721 lines | 68.6KB | 42 functions | 10 classes
**Modified:** 2026-01-05 15:31

**Description:**
```
VHH Designer v7 - IMGT-based Rules with Full Integration
=========================================================

Combines:
- ChatGPT's FR2 IMGT mapping fix (W at IMGT41, hallmarks at 42,49,50,52)
- New IMGT rules from 12M sequence epistasis analysis (JSON format)
- v5's NaturalnessScorer, VernierEngine, diverse candidate strategies
- Compensation rules support

For UNIVERSAL mode (grafting CDRs onto humanized scaffold):
  1. Start with Universal scaffold (FR2 = WFRQAPGQGLEAVA, hallmarks = FGLA)
  2. Apply IMGT FR3 rules (from JSON) - FR2 core is PROTECTED
  3. Apply safe compensation rules to FR3/FR4
    
For ORIGINAL mode (keeping input FRs with VHH-izing mutations):
  PASS 1: HALLMARK mutations (IMGT 42,49,50,52 → target pattern like FERG, YERL)
  PASS 2: IMGT FR3 rules (family-gated by hallmark pattern)
  PASS 3: COMPENSATION rules + VERNIER archetype

Family Classification (based on IMGT positions):
  - F_C2: pos42=F, pos49=E, pos50=R, 2 cysteines
  - Y_C2: pos42=Y, pos49=E, pos50=R, 2 cysteines  
  - F_C4: pos42=F, pos49=E, pos50=R, 4 cysteines
  - Y_C4: pos42=Y, pos49=E, pos50=R, 4 cysteines
  - VH_like: pos50=L (human VH-like)
  - Non_classical: other combinations

CRITICAL FR2 IMGT MAPPING (corrected):
  FR2: W-F-R-Q-A-P-G-Q-G-L-E-A-V-A
  IMGT: 41-42-43-44-45-46-47-48-49-50-51-52-53-54
  So: IMGT42=FR2[1], IMGT49=FR2[8], IMGT50=FR2[9], IMGT52=FR2[11] (0-based)

Usage:
  python vhh_designer_v7.py -i M69.fasta \
      --imgt-rules models/epistasis/imgt_v1/imgt_rules.json \
      --compensation models/correlations/correlation_results_v3_compensation.pkl \
      --epistasis models/epistasis/current/epistasis_v2_full.pkl \
      --target-hallmarks FERG \
      --mode both

Author: Claude (Anthropic) + ChatGPT fixes
Date: January 2026
```

---

### `vhh_designer_v7_10.py`
**Version:** 7.5
**Path:** `active/analysis/vhh_designer_v7_10.py`
**Stats:** 3231 lines | 124.3KB | 67 functions | 11 classes
**Modified:** 2026-01-18 23:21

**Description:**
```
VHH Designer v7.5 - Fixed AntPack Parsing & IMGT Position Mapping
==================================================================

FIXES from v7.4:
1. Correct AntPack 4-tuple parsing + trim_alignment() for insertion safety
2. IMGT position-based lookup (no more hardcoded FR2 indices)
3. Returns {imgt_pos: aa} dict for safe position access throughout
4. Proper handling of insertions (52A, 52B, etc.)
5. Region extraction aligned to IMGT numbering, not raw sequence indices

This version includes:
- Candidate generation (N candidates, default 500)
- ESM2 language model scoring (pseudo-perplexity)  
- ESMFold structure prediction (pLDDT confidence)
- Multi-family probabilistic rule compliance scoring
- Comprehensive CSV output with all metrics
- Lead candidates preserved and kept at top

Usage:
  python vhh_designer_v7_5_fixed.py -i input.fasta \
      --rules analysis_rules_v7.json \
      --archetypes analysis_vernier_archetypes_v7.json \
      --n-generate 500 \
      --n-select 92

Author: Claude (Anthropic)
Date: January 2026
```

---

### `vhh_designer_v7_11.py`
**Version:** 7.5
**Path:** `active/analysis/vhh_designer_v7_11.py`
**Stats:** 3240 lines | 124.8KB | 67 functions | 11 classes
**Modified:** 2026-01-19 00:27

**Description:**
```
VHH Designer v7.5 - Fixed AntPack Parsing & IMGT Position Mapping
==================================================================

FIXES from v7.4:
1. Correct AntPack 4-tuple parsing + trim_alignment() for insertion safety
2. IMGT position-based lookup (no more hardcoded FR2 indices)
3. Returns {imgt_pos: aa} dict for safe position access throughout
4. Proper handling of insertions (52A, 52B, etc.)
5. Region extraction aligned to IMGT numbering, not raw sequence indices

This version includes:
- Candidate generation (N candidates, default 500)
- ESM2 language model scoring (pseudo-perplexity)  
- ESMFold structure prediction (pLDDT confidence)
- Multi-family probabilistic rule compliance scoring
- Comprehensive CSV output with all metrics
- Lead candidates preserved and kept at top

Usage:
  python vhh_designer_v7_5_fixed.py -i input.fasta \
      --rules analysis_rules_v7.json \
      --archetypes analysis_vernier_archetypes_v7.json \
      --n-generate 500 \
      --n-select 92

Author: Claude (Anthropic)
Date: January 2026
```

---

### `vhh_designer_v7_12.py`
**Version:** 7.5
**Path:** `active/analysis/vhh_designer_v7_12.py`
**Stats:** 3247 lines | 125.0KB | 67 functions | 11 classes
**Modified:** 2026-01-19 18:12

**Description:**
```
VHH Designer v7.5 - Fixed AntPack Parsing & IMGT Position Mapping
==================================================================

FIXES from v7.4:
1. Correct AntPack 4-tuple parsing + trim_alignment() for insertion safety
2. IMGT position-based lookup (no more hardcoded FR2 indices)
3. Returns {imgt_pos: aa} dict for safe position access throughout
4. Proper handling of insertions (52A, 52B, etc.)
5. Region extraction aligned to IMGT numbering, not raw sequence indices

This version includes:
- Candidate generation (N candidates, default 500)
- ESM2 language model scoring (pseudo-perplexity)  
- ESMFold structure prediction (pLDDT confidence)
- Multi-family probabilistic rule compliance scoring
- Comprehensive CSV output with all metrics
- Lead candidates preserved and kept at top

Usage:
  python vhh_designer_v7_5_fixed.py -i input.fasta \
      --rules analysis_rules_v7.json \
      --archetypes analysis_vernier_archetypes_v7.json \
      --n-generate 500 \
      --n-select 92

Author: Claude (Anthropic)
Date: January 2026
```

---

### `vhh_designer_v7_13.py`
**Version:** 7.5
**Path:** `active/analysis/vhh_designer_v7_13.py`
**Stats:** 3539 lines | 140.3KB | 69 functions | 11 classes
**Modified:** 2026-01-20 13:58

**Description:**
```
VHH Designer v7.5 - Fixed AntPack Parsing & IMGT Position Mapping
==================================================================

FIXES from v7.4:
1. Correct AntPack 4-tuple parsing + trim_alignment() for insertion safety
2. IMGT position-based lookup (no more hardcoded FR2 indices)
3. Returns {imgt_pos: aa} dict for safe position access throughout
4. Proper handling of insertions (52A, 52B, etc.)
5. Region extraction aligned to IMGT numbering, not raw sequence indices

This version includes:
- Candidate generation (N candidates, default 500)
- ESM2 language model scoring (pseudo-perplexity)  
- ESMFold structure prediction (pLDDT confidence)
- Multi-family probabilistic rule compliance scoring
- Comprehensive CSV output with all metrics
- Lead candidates preserved and kept at top

Usage:
  python vhh_designer_v7_5_fixed.py -i input.fasta \
      --rules analysis_rules_v7.json \
      --archetypes analysis_vernier_archetypes_v7.json \
      --n-generate 500 \
      --n-select 92

Author: Claude (Anthropic)
Date: January 2026
```

---

### `vhh_designer_v7_14.1.py`
**Version:** 7.14
**Path:** `active/analysis/vhh_designer_v7_14.1.py`
**Stats:** 3732 lines | 149.8KB | 71 functions | 11 classes
**Modified:** 2026-01-20 17:53

**Description:**
```
VHH Designer v7.14 - Within-Track Ranking & Explicit Immutability
==================================================================

UPDATES from v7.12/7.13:
1. Within-track ranking (track_rank field) for control/probe candidates
2. Explicit immutability flag (is_immutable) for Track 0 minimal hallmarks
3. Distance-first sorting for Track 0 (4mut > 5mut, then ESM)
4. ESM-only sorting for Track 1-3 probes
5. Capped probe complexity: 6 singles, 6 pairs, 3 triplets (was 10/8/4)

Track System (v7.12+):
  Track 0 (minimal_hallmark): IMMUTABLE controls - hallmarks only, no verniers
  Track 1 (single_vernier): Probe one vernier at a time
  Track 2 (paired_vernier): Test vernier pairs (epistatic)
  Track 3 (triplet_vernier): Test vernier triplets (rare)
  Track 4 (optimized): Full consensus + rules + diversity (competitive ranking)

Key Insight: This system separates three scientific questions:
  A) Does the hallmark itself work? → Track 0
  B) Which verniers matter for this hallmark? → Track 1-3
  C) What looks like the best "real" VHH? → Track 4

Within-Track Sorting Rules (v7.14):
  Track 0: n_mutations ASC (4>5), then ESM loss ASC
  Track 1: ESM loss ASC (distance constant ~5-6 mutations)
  Track 2: ESM loss ASC, tie-break by n_mutations ASC
  Track 3: ESM loss ASC

This version includes:
- Candidate generation (N candidates, default 500)
- ESM2 language model scoring (pseudo-perplexity)  
- ESMFold structure prediction (pLDDT confidence)
- Multi-family probabilistic rule compliance scoring
- Comprehensive CSV output with all metrics
- Lead candidates preserved and kept at top

Usage:
  python vhh_designer_v7_14.py -i input.fasta \
      --rules analysis_rules_v7.json \
      --archetypes analysis_vernier_archetypes_v7.json \
      --n-generate 20000 \
      --n-select 99

Author: Claude (Anthropic)
Date: January 2026
```

---

### `vhh_designer_v7_14.py`
**Version:** 7.14
**Path:** `active/analysis/vhh_designer_v7_14.py`
**Stats:** 3654 lines | 146.4KB | 70 functions | 11 classes
**Modified:** 2026-01-20 15:46

**Description:**
```
VHH Designer v7.14 - Within-Track Ranking & Explicit Immutability
==================================================================

UPDATES from v7.12/7.13:
1. Within-track ranking (track_rank field) for control/probe candidates
2. Explicit immutability flag (is_immutable) for Track 0 minimal hallmarks
3. Distance-first sorting for Track 0 (4mut > 5mut, then ESM)
4. ESM-only sorting for Track 1-3 probes
5. Capped probe complexity: 6 singles, 6 pairs, 3 triplets (was 10/8/4)

Track System (v7.12+):
  Track 0 (minimal_hallmark): IMMUTABLE controls - hallmarks only, no verniers
  Track 1 (single_vernier): Probe one vernier at a time
  Track 2 (paired_vernier): Test vernier pairs (epistatic)
  Track 3 (triplet_vernier): Test vernier triplets (rare)
  Track 4 (optimized): Full consensus + rules + diversity (competitive ranking)

Key Insight: This system separates three scientific questions:
  A) Does the hallmark itself work? → Track 0
  B) Which verniers matter for this hallmark? → Track 1-3
  C) What looks like the best "real" VHH? → Track 4

Within-Track Sorting Rules (v7.14):
  Track 0: n_mutations ASC (4>5), then ESM loss ASC
  Track 1: ESM loss ASC (distance constant ~5-6 mutations)
  Track 2: ESM loss ASC, tie-break by n_mutations ASC
  Track 3: ESM loss ASC

This version includes:
- Candidate generation (N candidates, default 500)
- ESM2 language model scoring (pseudo-perplexity)  
- ESMFold structure prediction (pLDDT confidence)
- Multi-family probabilistic rule compliance scoring
- Comprehensive CSV output with all metrics
- Lead candidates preserved and kept at top

Usage:
  python vhh_designer_v7_14.py -i input.fasta \
      --rules analysis_rules_v7.json \
      --archetypes analysis_vernier_archetypes_v7.json \
      --n-generate 20000 \
      --n-select 99

Author: Claude (Anthropic)
Date: January 2026
```

---

### `vhh_designer_v7_15.py`
**Version:** 7.14
**Path:** `active/analysis/vhh_designer_v7_15.py`
**Stats:** 3831 lines | 153.8KB | 71 functions | 11 classes
**Modified:** 2026-01-20 18:45

**Description:**
```
VHH Designer v7.14 - Global Deduplication & Positive Ranking
==================================================================

UPDATES from v7.12/7.13:
1. Within-track ranking (track_rank field) for control/probe candidates
2. Explicit immutability flag (is_immutable) for Track 0 minimal hallmarks
3. Distance-first sorting for Track 0 (4mut > 5mut, then ESM)
4. ESM-only sorting for Track 1-3 probes
5. Capped probe complexity: 6 singles, 6 pairs, 3 triplets (was 10/8/4)
6. GLOBAL DEDUPLICATION: Remove duplicate sequences across all tracks
7. POSITIVE RANKS: All candidates get positive ranks (no negatives)
   - Lead: rank 0
   - Exempt controls (Track 0-3): rank 1, 2, 3, ...
   - Ranked candidates (Track 4): continue from controls
   - Use 'ranking_exempt' column to distinguish controls from ranked
8. Deduplication stats captured in summary JSON

Track System (v7.12+):
  Track 0 (minimal_hallmark): IMMUTABLE controls - hallmarks only, no verniers
  Track 1 (single_vernier): Probe one vernier at a time
  Track 2 (paired_vernier): Test vernier pairs (epistatic)
  Track 3 (triplet_vernier): Test vernier triplets (rare)
  Track 4 (optimized): Full consensus + rules + diversity (competitive ranking)

Key Insight: This system separates three scientific questions:
  A) Does the hallmark itself work? → Track 0
  B) Which verniers matter for this hallmark? → Track 1-3
  C) What looks like the best "real" VHH? → Track 4

Within-Track Sorting Rules (v7.14):
  Track 0: n_mutations ASC (4>5), then ESM loss ASC
  Track 1: ESM loss ASC (distance constant ~5-6 mutations)
  Track 2: ESM loss ASC, tie-break by n_mutations ASC
  Track 3: ESM loss ASC

This version includes:
- Candidate generation (N candidates, default 500)
- ESM2 language model scoring (pseudo-perplexity)  
- ESMFold structure prediction (pLDDT confidence)
- Multi-family probabilistic rule compliance scoring
- Comprehensive CSV output with all metrics
- Lead candidates preserved and kept at top

Usage:
  python vhh_designer_v7_14.py -i input.fasta \
      --rules analysis_rules_v7.json \
      --archetypes analysis_vernier_archetypes_v7.json \
      --n-generate 20000 \
      --n-select 99

Author: Claude (Anthropic)
Date: January 2026
```

---

### `vhh_designer_v7_16.1.py`
**Version:** 7.16
**Path:** `active/analysis/vhh_designer_v7_16.1.py`
**Stats:** 3978 lines | 158.8KB | 75 functions | 11 classes
**Modified:** 2026-01-20 22:33

**Description:**
```
VHH Designer v7.16 - Cysteine Validation & C4 Enforcement
==================================================================

UPDATES from v7.14/7.15:
1. CYSTEINE VALIDATION (new in v7.16):
   - Reject candidates with odd number of cysteines (disulfides require pairs)
   - C4 family enforcement: require BOTH IMGT 55 AND IMGT 100 = C
   - Cysteine validation applied BEFORE deduplication
   - Accurate cysteine-class detection based on actual positions
2. All v7.14 features retained (global dedup, positive ranking, track system)

Filter Order (v7.16):
  Generate → Cysteine validation → Deduplication → Scoring → Selection

Cysteine Rules:
  - Odd Cys count = REJECT (unpaired cysteine = aggregation risk)
  - C4 families: require pos55='C' AND pos100='C' (extra disulfide pair)
  - C2 families: should NOT have extra C55/C100 pair

Track System (v7.12+):
  Track 0 (minimal_hallmark): IMMUTABLE controls - hallmarks only, no verniers
  Track 1 (single_vernier): Probe one vernier at a time
  Track 2 (paired_vernier): Test vernier pairs (epistatic)
  Track 3 (triplet_vernier): Test vernier triplets (rare)
  Track 4 (optimized): Full consensus + rules + diversity (competitive ranking)

Key Insight: This system separates three scientific questions:
  A) Does the hallmark itself work? → Track 0
  B) Which verniers matter for this hallmark? → Track 1-3
  C) What looks like the best "real" VHH? → Track 4

Within-Track Sorting Rules (v7.14):
  Track 0: n_mutations ASC (4>5), then ESM loss ASC
  Track 1: ESM loss ASC (distance constant ~5-6 mutations)
  Track 2: ESM loss ASC, tie-break by n_mutations ASC
  Track 3: ESM loss ASC

This version includes:
- Candidate generation (N candidates, default 500)
- ESM2 language model scoring (pseudo-perplexity)  
- ESMFold structure prediction (pLDDT confidence)
- Multi-family probabilistic rule compliance scoring
- Comprehensive CSV output with all metrics
- Lead candidates preserved and kept at top

Usage:
  python vhh_designer_v7_16.py -i input.fasta \
      --rules analysis_rules_v7.json \
      --archetypes analysis_vernier_archetypes_v7.json \
      --n-generate 20000 \
      --n-select 99

Author: Claude (Anthropic)
Date: January 2026
```

---

### `vhh_designer_v7_16.2.py`
**Version:** 7.16
**Path:** `active/analysis/vhh_designer_v7_16.2.py`
**Stats:** 4003 lines | 160.3KB | 75 functions | 11 classes
**Modified:** 2026-01-21 14:10

**Description:**
```
VHH Designer v7.16 - Cysteine Validation & C4 Enforcement
==================================================================

UPDATES from v7.14/7.15:
1. CYSTEINE VALIDATION (new in v7.16):
   - Reject candidates with odd number of cysteines (disulfides require pairs)
   - C4 family enforcement: require BOTH IMGT 55 AND IMGT 100 = C
   - Cysteine validation applied BEFORE deduplication
   - Accurate cysteine-class detection based on actual positions
2. CONSISTENT C2/C4 CLASSIFICATION (fixed in v7.16):
   - C2/C4 labeling now ALWAYS based on IMGT55+IMGT100 positions
   - NOT based on total cysteine count (which was inconsistent)
   - Output includes: detected_family, target_family, detected_cys_class, n_cysteines
   - detected_cys_class: 'C2' (no extra pair), 'C4' (has C55+C100), 'INVALID' (one but not both)
3. All v7.14 features retained (global dedup, positive ranking, track system)

Filter Order (v7.16):
  Generate → Cysteine validation → Deduplication → Scoring → Selection

Cysteine Rules:
  - Odd Cys count = REJECT (unpaired cysteine = aggregation risk)
  - C4 families: require pos55='C' AND pos100='C' (extra disulfide pair)
  - C2 families: should NOT have extra C55/C100 pair

Track System (v7.12+):
  Track 0 (minimal_hallmark): IMMUTABLE controls - hallmarks only, no verniers
  Track 1 (single_vernier): Probe one vernier at a time
  Track 2 (paired_vernier): Test vernier pairs (epistatic)
  Track 3 (triplet_vernier): Test vernier triplets (rare)
  Track 4 (optimized): Full consensus + rules + diversity (competitive ranking)

Key Insight: This system separates three scientific questions:
  A) Does the hallmark itself work? → Track 0
  B) Which verniers matter for this hallmark? → Track 1-3
  C) What looks like the best "real" VHH? → Track 4

Within-Track Sorting Rules (v7.14):
  Track 0: n_mutations ASC (4>5), then ESM loss ASC
  Track 1: ESM loss ASC (distance constant ~5-6 mutations)
  Track 2: ESM loss ASC, tie-break by n_mutations ASC
  Track 3: ESM loss ASC

This version includes:
- Candidate generation (N candidates, default 500)
- ESM2 language model scoring (pseudo-perplexity)  
- ESMFold structure prediction (pLDDT confidence)
- Multi-family probabilistic rule compliance scoring
- Comprehensive CSV output with all metrics
- Lead candidates preserved and kept at top

Usage:
  python vhh_designer_v7_16.py -i input.fasta \
      --rules analysis_rules_v7.json \
      --archetypes analysis_vernier_archetypes_v7.json \
      --n-generate 20000 \
      --n-select 99

Author: Claude (Anthropic)
Date: January 2026
```

---

### `vhh_designer_v7_16.py`
**Version:** 7.16
**Path:** `active/analysis/vhh_designer_v7_16.py`
**Stats:** 3978 lines | 158.8KB | 75 functions | 11 classes
**Modified:** 2026-01-20 21:41

**Description:**
```
VHH Designer v7.16 - Cysteine Validation & C4 Enforcement
==================================================================

UPDATES from v7.14/7.15:
1. CYSTEINE VALIDATION (new in v7.16):
   - Reject candidates with odd number of cysteines (disulfides require pairs)
   - C4 family enforcement: require BOTH IMGT 55 AND IMGT 100 = C
   - Cysteine validation applied BEFORE deduplication
   - Accurate cysteine-class detection based on actual positions
2. All v7.14 features retained (global dedup, positive ranking, track system)

Filter Order (v7.16):
  Generate → Cysteine validation → Deduplication → Scoring → Selection

Cysteine Rules:
  - Odd Cys count = REJECT (unpaired cysteine = aggregation risk)
  - C4 families: require pos55='C' AND pos100='C' (extra disulfide pair)
  - C2 families: should NOT have extra C55/C100 pair

Track System (v7.12+):
  Track 0 (minimal_hallmark): IMMUTABLE controls - hallmarks only, no verniers
  Track 1 (single_vernier): Probe one vernier at a time
  Track 2 (paired_vernier): Test vernier pairs (epistatic)
  Track 3 (triplet_vernier): Test vernier triplets (rare)
  Track 4 (optimized): Full consensus + rules + diversity (competitive ranking)

Key Insight: This system separates three scientific questions:
  A) Does the hallmark itself work? → Track 0
  B) Which verniers matter for this hallmark? → Track 1-3
  C) What looks like the best "real" VHH? → Track 4

Within-Track Sorting Rules (v7.14):
  Track 0: n_mutations ASC (4>5), then ESM loss ASC
  Track 1: ESM loss ASC (distance constant ~5-6 mutations)
  Track 2: ESM loss ASC, tie-break by n_mutations ASC
  Track 3: ESM loss ASC

This version includes:
- Candidate generation (N candidates, default 500)
- ESM2 language model scoring (pseudo-perplexity)  
- ESMFold structure prediction (pLDDT confidence)
- Multi-family probabilistic rule compliance scoring
- Comprehensive CSV output with all metrics
- Lead candidates preserved and kept at top

Usage:
  python vhh_designer_v7_16.py -i input.fasta \
      --rules analysis_rules_v7.json \
      --archetypes analysis_vernier_archetypes_v7.json \
      --n-generate 20000 \
      --n-select 99

Author: Claude (Anthropic)
Date: January 2026
```

---

### `vhh_designer_v7_2.1.py`
**Version:** 7.2
**Path:** `active/analysis/vhh_designer_v7_2.1.py`
**Stats:** 1597 lines | 62.1KB | 27 functions | 7 classes
**Modified:** 2026-01-14 13:41

**Description:**
```
VHH Designer v7.2 - Compatible with v7.7 Analysis Outputs
==========================================================

Updated to work directly with the JSON outputs from vhh_analysis_unified_v7.7.py:
  - analysis_rules_v7.json (all rules: compensation + triplet)
  - analysis_vernier_archetypes_v7.json (family patterns)
  - analysis_vernier_clusters_v7.json (for naturalness scoring)

For UNIVERSAL mode (grafting CDRs onto humanized scaffold):
  1. Start with Universal scaffold (FR2 = WFRQAPGQGLEAVA, hallmarks = FGLA)
  2. Apply rules from JSON gated by hallmarks/CDR features
  3. FR2 core is PROTECTED (only hallmark positions can change)
    
For ORIGINAL mode (keeping input FRs with VHH-izing mutations):
  PASS 1: HALLMARK mutations (IMGT 42,49,50,52 → target pattern like FERG, YERL)
  PASS 2: Apply rules based on conditions (hallmarks, cdr3_len, cdr3_charge, etc.)

Family Classification (based on IMGT positions):
  - F_C2: pos42=F, 2 cysteines
  - F_C4: pos42=F, 4+ cysteines
  - Y_C2: pos42=Y, 2 cysteines  
  - Y_C4: pos42=Y, 4+ cysteines
  - VH_like: pos50=L (human VH-like)
  - VHH_W52: pos52=W
  - Other_VHH: other combinations

Usage:
  python vhh_designer_v7_2.py -i mouse_hc.fasta \
      --rules models/analysis/v7.7_*/analysis_rules_v7.json \
      --archetypes models/analysis/v7.7_*/analysis_vernier_archetypes_v7.json \
      --target-hallmarks FERG \
      --mode both

Author: Claude (Anthropic)
Date: January 2026
```

---

### `vhh_designer_v7_2.py`
**Version:** 7.2
**Path:** `active/analysis/vhh_designer_v7_2.py`
**Stats:** 1122 lines | 40.4KB | 26 functions | 7 classes
**Modified:** 2026-01-14 11:11

**Description:**
```
VHH Designer v7.2 - Compatible with v7.7 Analysis Outputs
==========================================================

Updated to work directly with the JSON outputs from vhh_analysis_unified_v7.7.py:
  - analysis_rules_v7.json (all rules: compensation + triplet)
  - analysis_vernier_archetypes_v7.json (family patterns)
  - analysis_vernier_clusters_v7.json (for naturalness scoring)

For UNIVERSAL mode (grafting CDRs onto humanized scaffold):
  1. Start with Universal scaffold (FR2 = WFRQAPGQGLEAVA, hallmarks = FGLA)
  2. Apply rules from JSON gated by hallmarks/CDR features
  3. FR2 core is PROTECTED (only hallmark positions can change)
    
For ORIGINAL mode (keeping input FRs with VHH-izing mutations):
  PASS 1: HALLMARK mutations (IMGT 42,49,50,52 → target pattern like FERG, YERL)
  PASS 2: Apply rules based on conditions (hallmarks, cdr3_len, cdr3_charge, etc.)

Family Classification (based on IMGT positions):
  - F_C2: pos42=F, 2 cysteines
  - F_C4: pos42=F, 4+ cysteines
  - Y_C2: pos42=Y, 2 cysteines  
  - Y_C4: pos42=Y, 4+ cysteines
  - VH_like: pos50=L (human VH-like)
  - VHH_W52: pos52=W
  - Other_VHH: other combinations

Usage:
  python vhh_designer_v7_2.py -i mouse_hc.fasta \
      --rules models/analysis/v7.7_*/analysis_rules_v7.json \
      --archetypes models/analysis/v7.7_*/analysis_vernier_archetypes_v7.json \
      --target-hallmarks FERG \
      --mode both

Author: Claude (Anthropic)
Date: January 2026
```

---

### `vhh_designer_v7_4_esmfold.py`
**Version:** 7.4
**Path:** `active/analysis/vhh_designer_v7_4_esmfold.py`
**Stats:** 1937 lines | 72.2KB | 42 functions | 9 classes
**Modified:** 2026-01-14 16:24

**Description:**
```
VHH Designer v7.4 - ESMFold Structure Prediction & Comprehensive Scoring
=========================================================================

This version includes:
1. Candidate generation (N candidates, default 500)
2. ESM2 language model scoring (pseudo-perplexity)
3. ESMFold structure prediction (pLDDT confidence)
4. Multi-family probabilistic rule compliance scoring
5. Comprehensive CSV output with all metrics
6. Lead candidates preserved and kept at top

Output CSV includes:
  - generation_order: Original sequence number (1, 2, 3...)
  - detected_family: Family classification
  - vernier_matches: Number of vernier position agreements
  - rules_passed / rules_total: Rule compliance counts
  - framework_identity_pct: % of original framework retained
  - construction_method: How the candidate was built
  - esm2_loss, plddt_mean, plddt_cdr3, combined_score, etc.

Usage:
  python vhh_designer_v7_4_esmfold.py -i input.fasta \
      --rules analysis_rules_v7.json \
      --archetypes analysis_vernier_archetypes_v7.json \
      --n-generate 500 \
      --n-select 92

Author: Claude (Anthropic)
Date: January 2026
```

---

### `vhh_designer_v7_5.1.py`
**Version:** 7.5
**Path:** `active/analysis/vhh_designer_v7_5.1.py`
**Stats:** 2001 lines | 74.1KB | 41 functions | 10 classes
**Modified:** 2026-01-14 17:27

**Description:**
```
VHH Designer v7.5 - Fixed AntPack Parsing & IMGT Position Mapping
==================================================================

FIXES from v7.4:
1. Correct AntPack 4-tuple parsing + trim_alignment() for insertion safety
2. IMGT position-based lookup (no more hardcoded FR2 indices)
3. Returns {imgt_pos: aa} dict for safe position access throughout
4. Proper handling of insertions (52A, 52B, etc.)
5. Region extraction aligned to IMGT numbering, not raw sequence indices

This version includes:
- Candidate generation (N candidates, default 500)
- ESM2 language model scoring (pseudo-perplexity)  
- ESMFold structure prediction (pLDDT confidence)
- Multi-family probabilistic rule compliance scoring
- Comprehensive CSV output with all metrics
- Lead candidates preserved and kept at top

Usage:
  python vhh_designer_v7_5_fixed.py -i input.fasta \
      --rules analysis_rules_v7.json \
      --archetypes analysis_vernier_archetypes_v7.json \
      --n-generate 500 \
      --n-select 92

Author: Claude (Anthropic)
Date: January 2026
```

---

### `vhh_designer_v7_5_esm.py`
**Version:** 7.5
**Path:** `active/analysis/vhh_designer_v7_5_esm.py`
**Stats:** 1979 lines | 73.0KB | 41 functions | 10 classes
**Modified:** 2026-01-14 17:01

**Description:**
```
VHH Designer v7.5 - Fixed AntPack Parsing & IMGT Position Mapping
==================================================================

FIXES from v7.4:
1. Correct AntPack 4-tuple parsing + trim_alignment() for insertion safety
2. IMGT position-based lookup (no more hardcoded FR2 indices)
3. Returns {imgt_pos: aa} dict for safe position access throughout
4. Proper handling of insertions (52A, 52B, etc.)
5. Region extraction aligned to IMGT numbering, not raw sequence indices

This version includes:
- Candidate generation (N candidates, default 500)
- ESM2 language model scoring (pseudo-perplexity)  
- ESMFold structure prediction (pLDDT confidence)
- Multi-family probabilistic rule compliance scoring
- Comprehensive CSV output with all metrics
- Lead candidates preserved and kept at top

Usage:
  python vhh_designer_v7_5_fixed.py -i input.fasta \
      --rules analysis_rules_v7.json \
      --archetypes analysis_vernier_archetypes_v7.json \
      --n-generate 500 \
      --n-select 92

Author: Claude (Anthropic)
Date: January 2026
```

---

### `vhh_designer_v7_6_esm.py`
**Version:** 7.5
**Path:** `active/analysis/vhh_designer_v7_6_esm.py`
**Stats:** 2083 lines | 77.1KB | 45 functions | 10 classes
**Modified:** 2026-01-16 21:30

**Description:**
```
VHH Designer v7.5 - Fixed AntPack Parsing & IMGT Position Mapping
==================================================================

FIXES from v7.4:
1. Correct AntPack 4-tuple parsing + trim_alignment() for insertion safety
2. IMGT position-based lookup (no more hardcoded FR2 indices)
3. Returns {imgt_pos: aa} dict for safe position access throughout
4. Proper handling of insertions (52A, 52B, etc.)
5. Region extraction aligned to IMGT numbering, not raw sequence indices

This version includes:
- Candidate generation (N candidates, default 500)
- ESM2 language model scoring (pseudo-perplexity)  
- ESMFold structure prediction (pLDDT confidence)
- Multi-family probabilistic rule compliance scoring
- Comprehensive CSV output with all metrics
- Lead candidates preserved and kept at top

Usage:
  python vhh_designer_v7_5_fixed.py -i input.fasta \
      --rules analysis_rules_v7.json \
      --archetypes analysis_vernier_archetypes_v7.json \
      --n-generate 500 \
      --n-select 92

Author: Claude (Anthropic)
Date: January 2026
```

---

### `vhh_designer_v7_7_vhhonly_debug.py`
**Version:** 7.5
**Path:** `active/analysis/vhh_designer_v7_7_vhhonly_debug.py`
**Stats:** 2491 lines | 94.9KB | 59 functions | 11 classes
**Modified:** 2026-01-18 21:30

**Description:**
```
VHH Designer v7.5 - Fixed AntPack Parsing & IMGT Position Mapping
==================================================================

FIXES from v7.4:
1. Correct AntPack 4-tuple parsing + trim_alignment() for insertion safety
2. IMGT position-based lookup (no more hardcoded FR2 indices)
3. Returns {imgt_pos: aa} dict for safe position access throughout
4. Proper handling of insertions (52A, 52B, etc.)
5. Region extraction aligned to IMGT numbering, not raw sequence indices

This version includes:
- Candidate generation (N candidates, default 500)
- ESM2 language model scoring (pseudo-perplexity)  
- ESMFold structure prediction (pLDDT confidence)
- Multi-family probabilistic rule compliance scoring
- Comprehensive CSV output with all metrics
- Lead candidates preserved and kept at top

Usage:
  python vhh_designer_v7_5_fixed.py -i input.fasta \
      --rules analysis_rules_v7.json \
      --archetypes analysis_vernier_archetypes_v7.json \
      --n-generate 500 \
      --n-select 92

Author: Claude (Anthropic)
Date: January 2026
```

---

### `vhh_designer_v7_8.py`
**Version:** 7.5
**Path:** `active/analysis/vhh_designer_v7_8.py`
**Stats:** 2891 lines | 109.4KB | 67 functions | 11 classes
**Modified:** 2026-01-18 21:41

**Description:**
```
VHH Designer v7.5 - Fixed AntPack Parsing & IMGT Position Mapping
==================================================================

FIXES from v7.4:
1. Correct AntPack 4-tuple parsing + trim_alignment() for insertion safety
2. IMGT position-based lookup (no more hardcoded FR2 indices)
3. Returns {imgt_pos: aa} dict for safe position access throughout
4. Proper handling of insertions (52A, 52B, etc.)
5. Region extraction aligned to IMGT numbering, not raw sequence indices

This version includes:
- Candidate generation (N candidates, default 500)
- ESM2 language model scoring (pseudo-perplexity)  
- ESMFold structure prediction (pLDDT confidence)
- Multi-family probabilistic rule compliance scoring
- Comprehensive CSV output with all metrics
- Lead candidates preserved and kept at top

Usage:
  python vhh_designer_v7_5_fixed.py -i input.fasta \
      --rules analysis_rules_v7.json \
      --archetypes analysis_vernier_archetypes_v7.json \
      --n-generate 500 \
      --n-select 92

Author: Claude (Anthropic)
Date: January 2026
```

---

### `vhh_designer_v7_9.py`
**Version:** 7.5
**Path:** `active/analysis/vhh_designer_v7_9.py`
**Stats:** 2915 lines | 110.1KB | 67 functions | 11 classes
**Modified:** 2026-01-18 22:20

**Description:**
```
VHH Designer v7.5 - Fixed AntPack Parsing & IMGT Position Mapping
==================================================================

FIXES from v7.4:
1. Correct AntPack 4-tuple parsing + trim_alignment() for insertion safety
2. IMGT position-based lookup (no more hardcoded FR2 indices)
3. Returns {imgt_pos: aa} dict for safe position access throughout
4. Proper handling of insertions (52A, 52B, etc.)
5. Region extraction aligned to IMGT numbering, not raw sequence indices

This version includes:
- Candidate generation (N candidates, default 500)
- ESM2 language model scoring (pseudo-perplexity)  
- ESMFold structure prediction (pLDDT confidence)
- Multi-family probabilistic rule compliance scoring
- Comprehensive CSV output with all metrics
- Lead candidates preserved and kept at top

Usage:
  python vhh_designer_v7_5_fixed.py -i input.fasta \
      --rules analysis_rules_v7.json \
      --archetypes analysis_vernier_archetypes_v7.json \
      --n-generate 500 \
      --n-select 92

Author: Claude (Anthropic)
Date: January 2026
```

---

### `vhh_designer_v7_91.py`
**Version:** 7.5
**Path:** `active/analysis/vhh_designer_v7_91.py`
**Stats:** 2961 lines | 111.9KB | 67 functions | 11 classes
**Modified:** 2026-01-18 22:36

**Description:**
```
VHH Designer v7.5 - Fixed AntPack Parsing & IMGT Position Mapping
==================================================================

FIXES from v7.4:
1. Correct AntPack 4-tuple parsing + trim_alignment() for insertion safety
2. IMGT position-based lookup (no more hardcoded FR2 indices)
3. Returns {imgt_pos: aa} dict for safe position access throughout
4. Proper handling of insertions (52A, 52B, etc.)
5. Region extraction aligned to IMGT numbering, not raw sequence indices

This version includes:
- Candidate generation (N candidates, default 500)
- ESM2 language model scoring (pseudo-perplexity)  
- ESMFold structure prediction (pLDDT confidence)
- Multi-family probabilistic rule compliance scoring
- Comprehensive CSV output with all metrics
- Lead candidates preserved and kept at top

Usage:
  python vhh_designer_v7_5_fixed.py -i input.fasta \
      --rules analysis_rules_v7.json \
      --archetypes analysis_vernier_archetypes_v7.json \
      --n-generate 500 \
      --n-select 92

Author: Claude (Anthropic)
Date: January 2026
```

---

### `vhh_designer_v8_0.py`
**Version:** 8.0
**Path:** `active/analysis/vhh_designer_v8_0.py`
**Stats:** 4248 lines | 171.0KB | 80 functions | 11 classes
**Modified:** 2026-01-21 18:45

**Description:**
```
VHH Designer v8.0 - Universal Framework & Aggressive Vernier Testing
======================================================================

UPDATES from v7.16:

1. UNIVERSAL FRAMEWORK:
   - New "Universal" scaffold for CDR grafting
   - FR1: QVQLVESGGGLVQPGGSLRLSCAASG
   - FR2: WFRQAPGQGLEAVA  
   - FR3: YYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYC
   - FR4: WGQGTLVTVSS

2. EXPANDED TRACK 0 CONTROLS:
   For ALL tested hallmarks:
     a) Hallmarks only (4 mut)
     b) Hallmarks + IMGT2 (5 mut)
     c) Hallmarks + ALL verniers (full conversion)
   For Universal scaffold:
     d) Pure CDR graft into Universal
     e) Universal + hallmarks
     f) Universal + hallmarks + verniers
   ALL controls include IMGT2→V

3. AGGRESSIVE TRACK 1-4:
   Track 1: Load-bearing verniers only (12 positions)
   Track 2: Cross-family motif pairs (top 40)
   Track 3: Motif triplets + factorial design (2^6 = 64 combos)
   Track 4: Two-lane (Lane A: full consensus, Lane B: minimal 2-6 verniers)

4. RETAINED FROM v7.16:
   - Cysteine validation (odd-Cys rejection, C4 enforcement)
   - C2/C4 classification based on IMGT55+IMGT100
   - True VHH filtering, ESM2/ESMFold scoring
   - Immutable Track 0 controls

PHILOSOPHY: The goal is to find BINDERS. Each track explores a different
region of mutation space. We don't expect Track 0 to be "worse" than Track 4 -
we're casting a wide net to maximize chances of finding functional variants.

Track System:
  Track 0: Controls - test fundamental approaches (hallmarks only, full conversion, Universal)
  Track 1: Single vernier probes - which individual positions matter?
  Track 2: Paired vernier probes - which combinations co-occur naturally?
  Track 3: Triplet/factorial - systematic combination testing
  Track 4: Optimized candidates - full exploration with scoring

Usage:
  python vhh_designer_v8_0.py -i input.fasta \
      --rules analysis_rules_v7.json \
      --archetypes analysis_vernier_archetypes_v7.json \
      --n-generate 20000 \
      --n-select 99

Author: Claude (Anthropic)
Date: January 2026
```

---

### `vhh_designer_v8_1.py`
**Version:** 8.1
**Path:** `active/analysis/vhh_designer_v8_1.py`
**Stats:** 4309 lines | 173.8KB | 80 functions | 11 classes
**Modified:** 2026-01-21 19:16

**Description:**
```
VHH Designer v8.1 - FR1/FR4 Consensus & Track Budget Control
=============================================================

UPDATES from v8.0:

1. FR1/FR4 MANDATORY CONSENSUS (v8.1):
   - IMGT1: E→Q (~92% in VHH) - the "Q" in QVQL
   - IMGT128: A→S (~95% in VHH) - makes TVSS ending
   - Progressive sprinkling: Track 1 (30%) → Track 2 (50%) → Track 3 (70%) → Track 4 (100%)

2. TRACK BUDGET CONTROL (v8.1):
   - 92 controls (CONTROL status) across all tracks/families
   - 94 ranked candidates (RANKED status)
   - Reduced Track 1-3 counts to fit budget:
     * Track 1: 6 singles per hallmark (was 12)
     * Track 2: 8 pairs per hallmark (was 40)
     * Track 3: 5 triplets + 4 factorial samples (was 15+64)

3. CLEAR CSV LABELING (v8.1):
   - New 'status' column: LEAD, CONTROL, or RANKED
   - New 'track_subtype' column: lane_A, lane_B, factorial, motif_triplet
   - New 'scaffold_type' column: original, universal

4. UNIVERSAL FRAMEWORK (from v8.0):
   - FR1: QVQLVESGGGLVQPGGSLRLSCAASG (already starts with Q, ends with TVSS)
   - Generates controls: graft-only, +hallmarks, +full verniers

PHILOSOPHY: Find BINDERS by exploring different mutation strategies.
Track 0 controls test fundamentals, Tracks 1-3 are ranked controls with
progressive FR1/FR4 consensus, Track 4 is fully optimized (always has FR1/FR4).

Usage:
  python vhh_designer_v8_1.py -i input.fasta \
      --rules analysis_rules_v7.json \
      --archetypes analysis_vernier_archetypes_v7.json \
      --hallmark-db comprehensive_subfamily_analysis_imgt.xlsx \
      --target-hallmarks AUTO \
      --n-generate 20000 \
      --n-select 186

Author: Claude (Anthropic)
Date: January 2026
```

---

### `vhh_designer_v8_2.py`
**Version:** 8.2
**Path:** `active/analysis/vhh_designer_v8_2.py`
**Stats:** 4420 lines | 179.2KB | 80 functions | 11 classes
**Modified:** 2026-01-21 20:05

**Description:**
```
VHH Designer v8.2 - Universal as Subfamily & Equal Distribution
================================================================

UPDATES from v8.1:

1. UNIVERSAL AS SUBFAMILY (v8.2):
   - Universal scaffold now treated as a subfamily in Track 4
   - Gets equal allocation alongside YQRL, FERF, FERG, etc.
   - Vernier positions extrapolated from F_C2 family (most similar)

2. EQUAL TRACK 4 DISTRIBUTION (v8.2):
   - Track 4 candidates distributed equally across ALL subfamilies
   - Includes: Universal, YQRL, FERF, FERG, YERL, etc.
   - Ensures comprehensive coverage of all scaffold types

3. CLEAR RANKING FORMAT (v8.2):
   - Controls show "(Control)" in rank column
   - Ranked candidates show 1, 2, 3... (1 = best)
   - New 'overall_row' column for absolute position

4. OUTPUT BUDGET (v8.2):
   - Fixed at 198 total sequences
   - ~90-100 controls + ~90-100 ranked
   - Balanced across all subfamilies

5. FR1/FR4 from v8.1 (retained):
   - IMGT1 E→Q and IMGT128 A→S progressively applied
   - Track 1 (30%) → Track 2 (50%) → Track 3 (70%) → Track 4 (100%)

Usage:
  python vhh_designer_v8_2.py -i input.fasta \
      --rules analysis_rules_v7.json \
      --archetypes analysis_vernier_archetypes_v7.json \
      --hallmark-db comprehensive_subfamily_analysis_imgt.xlsx \
      --target-hallmarks AUTO \
      --n-generate 20000 \
      --n-select 198

Author: Claude (Anthropic)
Date: January 2026
```

---

### `vhh_designer_v8_3.py`
**Version:** 8.3
**Path:** `active/analysis/vhh_designer_v8_3.py`
**Stats:** 4669 lines | 190.2KB | 80 functions | 11 classes
**Modified:** 2026-01-21 22:09

**Description:**
```
VHH Designer v8.3 - Probabilistic Consensus Sprinkling & Motif Tracks
======================================================================

UPDATES from v8.2:

1. PROBABILISTIC CONSENSUS SPRINKLING (v8.3):
   - Build universal_prior[pos] = (top_aa, freq) aggregated across families
   - Tier A (≥90%): Apply with ~60% probability
   - Tier B (75-90%): Apply with ~25% probability
   - No more "all-or-nothing" - intermediate states now generated

2. MOTIF-BASED TRACKS (v8.3):
   - Track 2: Motif PAIRS (1-3 co-occurring pairs per candidate)
   - Track 3A: Motif TRIPLETS (1-2 triplets per candidate)
   - Track 3B: Partial consensus BUNDLES (3-6 Tier-A positions, sampled)
   - Not random verniers - uses cross-family enriched motifs

3. WITHIN-TRACK RANKING (v8.3):
   - Tracks 1-3 ranked by ESM+rules WITHIN their track
   - Track 4 doesn't dominate - diagnostic value preserved
   - Controls are still controls, not competing with optimized

4. VERSION IN FILENAMES (v8.3):
   - Output files now use VERSION constant (v8_3)
   - Example: 20260121_200000_v8_3_n20000_M69

5. RETAINED FROM v8.2:
   - Universal as subfamily with equal Track 4 allocation
   - Clear ranking format: "(Control)" vs 1/2/3...
   - 198 total output (99 controls + 99 ranked)

Usage:
  python vhh_designer_v8_3.py -i input.fasta \
      --rules analysis_rules_v7.json \
      --archetypes analysis_vernier_archetypes_v7.json \
      --hallmark-db comprehensive_subfamily_analysis_imgt.xlsx \
      --target-hallmarks AUTO \
      --n-generate 20000 \
      --n-select 198

Author: Claude (Anthropic)
Date: January 2026
```

---

### `vhh_designer_v8_3_fixed.py`
**Version:** 8.3
**Path:** `active/analysis/vhh_designer_v8_3_fixed.py`
**Stats:** 4688 lines | 191.3KB | 80 functions | 11 classes
**Modified:** 2026-01-21 23:07

**Description:**
```
VHH Designer v8.3 - Probabilistic Consensus Sprinkling & Motif Tracks
======================================================================

UPDATES from v8.2:

1. PROBABILISTIC CONSENSUS SPRINKLING (v8.3):
   - Build universal_prior[pos] = (top_aa, freq) aggregated across families
   - Tier A (≥90%): Apply with ~60% probability
   - Tier B (75-90%): Apply with ~25% probability
   - No more "all-or-nothing" - intermediate states now generated

2. MOTIF-BASED TRACKS (v8.3):
   - Track 2: Motif PAIRS (1-3 co-occurring pairs per candidate)
   - Track 3A: Motif TRIPLETS (1-2 triplets per candidate)
   - Track 3B: Partial consensus BUNDLES (3-6 Tier-A positions, sampled)
   - Not random verniers - uses cross-family enriched motifs

3. WITHIN-TRACK RANKING (v8.3):
   - Tracks 1-3 ranked by ESM+rules WITHIN their track
   - Track 4 doesn't dominate - diagnostic value preserved
   - Controls are still controls, not competing with optimized

4. VERSION IN FILENAMES (v8.3):
   - Output files now use VERSION constant (v8_3)
   - Example: 20260121_200000_v8_3_n20000_M69

5. RETAINED FROM v8.2:
   - Universal as subfamily with equal Track 4 allocation
   - Clear ranking format: "(Control)" vs 1/2/3...
   - 198 total output (99 controls + 99 ranked)

Usage:
  python vhh_designer_v8_3.py -i input.fasta \
      --rules analysis_rules_v7.json \
      --archetypes analysis_vernier_archetypes_v7.json \
      --hallmark-db comprehensive_subfamily_analysis_imgt.xlsx \
      --target-hallmarks AUTO \
      --n-generate 20000 \
      --n-select 198

Author: Claude (Anthropic)
Date: January 2026
```

---

### `vhh_designer_v8_4.py`
**Version:** 8.4
**Path:** `active/analysis/vhh_designer_v8_4.py`
**Stats:** 4900 lines | 200.4KB | 84 functions | 12 classes
**Modified:** 2026-01-26 17:44

**Description:**
```
VHH Designer v8.4 - CDR-Conditional Framework Tracks
======================================================================

UPDATES from v8.2:

1. PROBABILISTIC CONSENSUS SPRINKLING (v8.3):
   - Build universal_prior[pos] = (top_aa, freq) aggregated across families
   - Tier A (≥90%): Apply with ~60% probability
   - Tier B (75-90%): Apply with ~25% probability
   - No more "all-or-nothing" - intermediate states now generated

2. MOTIF-BASED TRACKS (v8.3):
   - Track 2: Motif PAIRS (1-3 co-occurring pairs per candidate)
   - Track 3A: Motif TRIPLETS (1-2 triplets per candidate)
   - Track 3B: Partial consensus BUNDLES (3-6 Tier-A positions, sampled)
   - Not random verniers - uses cross-family enriched motifs

3. WITHIN-TRACK RANKING (v8.3):
   - Tracks 1-3 ranked by ESM+rules WITHIN their track
   - Track 4 doesn't dominate - diagnostic value preserved
   - Controls are still controls, not competing with optimized

4. VERSION IN FILENAMES (v8.3):
   - Output files now use VERSION constant (v8_3)
   - Example: 20260121_200000_v8_3_n20000_M69

5. RETAINED FROM v8.2:
   - Universal as subfamily with equal Track 4 allocation
   - Clear ranking format: "(Control)" vs 1/2/3...
   - 198 total output (99 controls + 99 ranked)

Usage:
  python vhh_designer_v8_3.py -i input.fasta \
      --rules analysis_rules_v7.json \
      --archetypes analysis_vernier_archetypes_v7.json \
      --hallmark-db comprehensive_subfamily_analysis_imgt.xlsx \
      --target-hallmarks AUTO \
      --n-generate 20000 \
      --n-select 198

Author: Claude (Anthropic)
Date: January 2026
```

---

### `vhh_designer_v8_5.py`
**Version:** 8.5
**Path:** `active/analysis/vhh_designer_v8_5.py`
**Stats:** 5355 lines | 218.0KB | 86 functions | 12 classes
**Modified:** 2026-01-26 22:20

**Description:**
```
VHH Designer v8.5 - Framework Consensus Tiers & Expansive Track 4
==================================================================

UPDATES from v8.4:

1. FRAMEWORK CONSENSUS TIERS (v8.5):
   - Tier S (>=95%): 29 positions - near-universal, applied in Track 2
   - Tier A (90-95%): 11 positions - very high, applied in Track 3
   - Tier B (85-90%): 16 positions - high, applied probabilistically in Track 4
   - Tier C (75-85%): 3 positions - medium, applied probabilistically in Track 4
   
   These are cross-hallmark consensus derived from 12M VHH analysis.

2. CDR-CONDITIONAL RULES IN TRACK 3 (v8.5):
   - Track 3 now includes CDR-conditional rules based on CDR3 features
   - Rules triggered by: cdr3[-1]=X, cdr3_charge, cdr3_len, cdr3_cys
   - Applied to BOTH Original and Universal frameworks

3. EXPANSIVE TRACK 4 (v8.5):
   - Track 4 now SCALES complexity when mutation space is exhausted
   - If n_unique < n_requested, probability of adding mutations increases
   - Layers: CDR rules > R50 > Tier S > Tier A > Tier B > Tier C
   - Never violates: cysteine pairing, compensation rules

4. TRACK STRUCTURE (v8.5):
   
   ORIGINAL FRAMEWORK:
   - Track 0a: Hallmarks only (4 mut)
   - Track 0b: + IMGT2 (5 mut)
   - Track 0c: + All verniers
   - Track 1: Single vernier probes
   - Track 2: + Tier S framework (>=95%)  [NEW]
   - Track 3: + Tier A + CDR-conditional   [NEW]
   - Track 4: EXPANSIVE (Tier B/C + scaling)
   
   UNIVERSAL FRAMEWORK:
   - Track 0: Control (hallmarks only)
   - Track 1: High confidence CDR-conditional
   - Track 2: + Tier S framework
   - Track 3: + Tier A + override humanization
   - Track 4: EXPANSIVE (all tiers + scaling)

5. YQRL HANDLING (unchanged from v8.4):
   - Only Track 0 controls for YQRL (rigid scaffold)

Usage:
  python vhh_designer_v8_5.py -i input.fasta \
      --rules analysis_rules_v7_all_positions.json \
      --archetypes analysis_vernier_archetypes_v7.json \
      --hallmark-db comprehensive_subfamily_analysis_imgt.xlsx \
      --target-hallmarks AUTO \
      --n-generate 200000 \
      --n-select 198

Author: Claude (Anthropic)
Date: January 2026
```

---

### `vhh_designer_v8_6.py`
**Version:** 8.6
**Path:** `active/analysis/vhh_designer_v8_6.py`
**Stats:** 5863 lines | 236.7KB | 98 functions | 12 classes
**Modified:** 2026-01-27 19:54

**Description:**
```
VHH Designer v8.6 - Track Quota Selection & Bernoulli Sprinkling
==================================================================

UPDATES from v8.5:

1. TRACKS 1-3 ARE NOW RANKED LANES (v8.6):
   - Previously Tracks 1-3 were exempt controls only
   - Now Tracks 1-3 compete for selection alongside Track 4
   - Only Track 0 remains exempt (pure controls)
   - Selection uses per-track per-bucket quotas

2. PER-BUCKET PER-TRACK QUOTAS (v8.6):
   - Selection enforces quotas within each bucket (hallmark × scaffold × family)
   - Configurable quotas: --track0-quota, --track1-quota, etc.
   - Ensures representation from all tracks in final output
   - Final output count is EXACTLY --n-select (no more exceeding)

3. BERNOULLI SPRINKLING WITH TEMPERATURE (v8.6):
   - Independent probability per position (not all-or-nothing)
   - Temperature scaling: --sprinkle-temp (default 1.0)
   - Lower temperature = sharper probabilities
   - Higher temperature = softer probabilities

4. WEIGHTED CHOICE FOR POSITION 66 (v8.6):
   - Position 66 uses actual frequencies for consensus vs alternative
   - No more 50/50 randomization
   - weighted_choice([consensus, alt], [freq_cons, freq_alt])

5. DEDUP WITHIN TRACK (v8.6):
   - Deduplication preserves track intent
   - Same sequence in different tracks stays separate
   - Configurable: --dedup-mode [within_track|global|none]

6. ARCHETYPE MOTIF POOL (v8.6):
   - Motifs extracted from archetypes (2-5 positions)
   - Future: will replace enumerated pairs/triplets in Tracks 2-3

TRACK STRUCTURE:

  Track 0: Controls (exempt) - Hallmarks only
  Track 1: Single vernier probes (RANKED)
  Track 2: Paired/motif + Tier S framework (RANKED)
  Track 3: Triplet/motif + Tier A + CDR rules (RANKED)
  Track 4: Expansive optimization (RANKED)

Usage:
  python vhh_designer_v8_6.py -i input.fasta \
      --rules analysis_rules_v7_all_positions.json \
      --archetypes analysis_vernier_archetypes_v7.json \
      --hallmark-db comprehensive_subfamily_analysis_imgt.xlsx \
      --target-hallmarks AUTO \
      --n-generate 200000 \
      --n-select 198 \
      --sprinkle-temp 1.0 \
      --dedup-mode within_track

Author: Claude (Anthropic)
Date: January 2026
```

---

### `vhh_designer_v8_6_patched.py`
**Version:** 8.6
**Path:** `active/analysis/vhh_designer_v8_6_patched.py`
**Stats:** 5606 lines | 224.8KB | 98 functions | 12 classes
**Modified:** 2026-01-27 21:11

**Description:**
```
VHH Designer v8.6 - Track Quota Selection & Bernoulli Sprinkling
==================================================================

UPDATES from v8.5:

1. TRACKS 1-3 ARE NOW RANKED LANES (v8.6):
   - Previously Tracks 1-3 were exempt controls only
   - Now Tracks 1-3 compete for selection alongside Track 4
   - Only Track 0 remains exempt (pure controls)
   - Selection uses per-track per-bucket quotas

2. PER-BUCKET PER-TRACK QUOTAS (v8.6):
   - Selection enforces quotas within each bucket (hallmark × scaffold × family)
   - Configurable quotas: --track0-quota, --track1-quota, etc.
   - Ensures representation from all tracks in final output
   - Final output count is EXACTLY --n-select (no more exceeding)

3. BERNOULLI SPRINKLING WITH TEMPERATURE (v8.6):
   - Independent probability per position (not all-or-nothing)
   - Temperature scaling: --sprinkle-temp (default 1.0)
   - Lower temperature = sharper probabilities
   - Higher temperature = softer probabilities

4. WEIGHTED CHOICE FOR POSITION 66 (v8.6):
   - Position 66 uses actual frequencies for consensus vs alternative
   - No more 50/50 randomization
   - weighted_choice([consensus, alt], [freq_cons, freq_alt])

5. DEDUP WITHIN TRACK (v8.6):
   - Deduplication preserves track intent
   - Same sequence in different tracks stays separate
   - Configurable: --dedup-mode [within_track|global|none]

6. ARCHETYPE MOTIF POOL (v8.6):
   - Motifs extracted from archetypes (2-5 positions)
   - Future: will replace enumerated pairs/triplets in Tracks 2-3

TRACK STRUCTURE:

  Track 0: Controls (exempt) - Hallmarks only
  Track 1: Single vernier probes (RANKED)
  Track 2: Paired/motif + Tier S framework (RANKED)
  Track 3: Triplet/motif + Tier A + CDR rules (RANKED)
  Track 4: Expansive optimization (RANKED)

Usage:
  python vhh_designer_v8_6.py -i input.fasta \
      --rules analysis_rules_v7_all_positions.json \
      --archetypes analysis_vernier_archetypes_v7.json \
      --hallmark-db comprehensive_subfamily_analysis_imgt.xlsx \
      --target-hallmarks AUTO \
      --n-generate 200000 \
      --n-select 198 \
      --sprinkle-temp 1.0 \
      --dedup-mode within_track

Author: Claude (Anthropic)
Date: January 2026
```

---

### `vhh_designer_v8_8.py`
**Version:** 8.8
**Path:** `active/analysis/vhh_designer_v8_8.py`
**Stats:** 5854 lines | 236.1KB | 98 functions | 12 classes
**Modified:** 2026-01-27 22:52

**Description:**
```
VHH Designer v8.8 - Archetype Motif Sampling & Per-Position Sprinkling
========================================================================

UPDATES from v8.6/v8.7:

1. TRACKS 2-3 USE ARCHETYPE MOTIFS (v8.8):
   - Track 2: samples 2-3 position motifs from archetype pool
   - Track 3: samples 3-5 position motifs from archetype pool
   - Replaces enumerated MOTIF_PAIRS/MOTIF_TRIPLETS with learned patterns
   - Motif pools built from analysis_vernier_archetypes_v7.json

2. PER-POSITION BERNOULLI FOR MANDATORY CONSENSUS (v8.8):
   - Removed global 50%/70% sprinkling gate in Tracks 2-3
   - Each mandatory_consensus position now sprinkled independently
   - Probability = temp_scale_prob(mc.confidence, temperature)
   - No "all-or-nothing" behavior

3. STRICT OUTPUT SIZE ENFORCEMENT (v8.8):
   - Selection now raises RuntimeError if it cannot reach n_select_total
   - Fail-fast instead of silent shortfall
   - Clear error message with suggestions

4. TRACK 2/3 SCALING FOR QUOTAS (v8.8):
   - Track 2/3 generation scales with --n-generate (5% each, min 30)
   - Multiple sprinkling variants per motif for diversity
   - Vernier consensus also sprinkled (30% for Track 2, 50% for Track 3)
   - Makes track quotas actually achievable

5. LEAD HANDLING FIX (v8.8):
   - Lead candidate kept completely separate from selection pool
   - No more duplicate lead rows in output
   - Exact --n-select output count

INHERITED FROM v8.6:
- Tracks 1-3 are ranked lanes (not exempt)
- Per-bucket per-track quotas for selection
- Temperature scaling for probabilities
- Weighted choice for position 66
- Dedup within track

TRACK STRUCTURE:

  Track 0: Controls (exempt) - Hallmarks only
  Track 1: Single vernier probes (RANKED)
  Track 2: Archetype motifs 2-3 pos (RANKED)
  Track 3: Archetype motifs 3-5 pos (RANKED)
  Track 4: Expansive optimization (RANKED)

Usage:
  python vhh_designer_v8_8.py -i input.fasta \
      --rules analysis_rules_v7_all_positions.json \
      --archetypes analysis_vernier_archetypes_v7.json \
      --hallmark-db comprehensive_subfamily_analysis_imgt.xlsx \
      --target-hallmarks AUTO \
      --n-generate 200000 \
      --n-select 198

Author: Claude (Anthropic)
Date: January 2026
```

---

### `vhh_designer_v8_8_1.py`
**Version:** 8.8.1
**Path:** `active/analysis/vhh_designer_v8_8_1.py`
**Stats:** 5966 lines | 241.6KB | 98 functions | 12 classes
**Modified:** 2026-01-28 12:22

**Description:**
```
VHH Designer v8.8.1 - Archetype Motif Sampling & Per-Position Sprinkling
========================================================================

UPDATES from v8.6/v8.7:

1. TRACKS 2-3 USE ARCHETYPE MOTIFS (v8.8):
   - Track 2: samples 2-3 position motifs from archetype pool
   - Track 3: samples 3-5 position motifs from archetype pool
   - Replaces enumerated MOTIF_PAIRS/MOTIF_TRIPLETS with learned patterns
   - Motif pools built from analysis_vernier_archetypes_v7.json

2. PER-POSITION BERNOULLI FOR MANDATORY CONSENSUS (v8.8):
   - Removed global 50%/70% sprinkling gate in Tracks 2-3
   - Each mandatory_consensus position now sprinkled independently
   - Probability = temp_scale_prob(mc.confidence, temperature)
   - No "all-or-nothing" behavior

3. STRICT OUTPUT SIZE ENFORCEMENT (v8.8):
   - Selection now raises RuntimeError if it cannot reach n_select_total
   - Fail-fast instead of silent shortfall
   - Clear error message with suggestions

4. TRACK 2/3 SCALING FOR QUOTAS (v8.8):
   - Track 2/3 generation scales with --n-generate (5% each, min 30)
   - Multiple sprinkling variants per motif for diversity
   - Vernier consensus also sprinkled (30% for Track 2, 50% for Track 3)
   - Makes track quotas actually achievable
   - Both Original AND Universal scaffolds now generate Track 2/3

5. LEAD HANDLING FIX (v8.8):
   - Lead candidate kept completely separate from selection pool
   - No more duplicate lead rows in output
   - Exact --n-select output count

6. CONTROLS ALWAYS INCLUDED (v8.8):
   - ALL controls (ranking_exempt=True) are included in final output
   - Controls listed first after lead, then ranked candidates
   - Controls have bracket notation: [1], [2], [3]...
   - Ranked candidates have normal ranks: 4, 5, 6...

7. ENHANCED SELECTION SUMMARY (v8.8.1):
   - Terminal output shows breakdown of ranked candidates by track
   - Shows breakdown by hallmark (subfamily)
   - Shows breakdown by scaffold type (original vs universal)
   - All stats included in JSON summary

INHERITED FROM v8.6:
- Tracks 1-3 are ranked lanes (not exempt)
- Per-bucket per-track quotas for selection
- Temperature scaling for probabilities
- Weighted choice for position 66
- Dedup within track

TRACK STRUCTURE:

  Track 0: Controls (exempt) - Hallmarks only
  Track 1: Single vernier probes (RANKED)
  Track 2: Archetype motifs 2-3 pos (RANKED)
  Track 3: Archetype motifs 3-5 pos (RANKED)
  Track 4: Expansive optimization (RANKED)

Usage:
  python vhh_designer_v8_8.py -i input.fasta \
      --rules analysis_rules_v7_all_positions.json \
      --archetypes analysis_vernier_archetypes_v7.json \
      --hallmark-db comprehensive_subfamily_analysis_imgt.xlsx \
      --target-hallmarks AUTO \
      --n-generate 200000 \
      --n-select 198

Author: Claude (Anthropic)
Date: January 2026
```

---

### `vhh_designer_v8_8_2.py`
**Version:** 8.8.2
**Path:** `active/analysis/vhh_designer_v8_8_2.py`
**Stats:** 6193 lines | 252.8KB | 98 functions | 12 classes
**Modified:** 2026-01-28 19:06

**Description:**
```
VHH Designer v8.8.2 - Archetype Motif Sampling & Per-Position Sprinkling
========================================================================

UPDATES from v8.6/v8.7:

1. TRACKS 2-3 USE ARCHETYPE MOTIFS (v8.8):
   - Track 2: samples 2-3 position motifs from archetype pool
   - Track 3: samples 3-5 position motifs from archetype pool
   - Replaces enumerated MOTIF_PAIRS/MOTIF_TRIPLETS with learned patterns
   - Motif pools built from analysis_vernier_archetypes_v7.json

2. PER-POSITION BERNOULLI FOR MANDATORY CONSENSUS (v8.8):
   - Removed global 50%/70% sprinkling gate in Tracks 2-3
   - Each mandatory_consensus position now sprinkled independently
   - Probability = temp_scale_prob(mc.confidence, temperature)
   - No "all-or-nothing" behavior

3. STRICT OUTPUT SIZE ENFORCEMENT (v8.8):
   - Selection now raises RuntimeError if it cannot reach n_select_total
   - Fail-fast instead of silent shortfall
   - Clear error message with suggestions

4. TRACK 2/3 SCALING FOR QUOTAS (v8.8):
   - Track 2/3 generation scales with --n-generate (5% each, min 30)
   - Multiple sprinkling variants per motif for diversity
   - Vernier consensus also sprinkled (30% for Track 2, 50% for Track 3)
   - Makes track quotas actually achievable
   - Both Original AND Universal scaffolds now generate Track 2/3

5. LEAD HANDLING FIX (v8.8):
   - Lead candidate kept completely separate from selection pool
   - No more duplicate lead rows in output
   - Exact --n-select output count

6. CONTROLS ALWAYS INCLUDED (v8.8):
   - ALL controls (ranking_exempt=True) are included in final output
   - Controls listed first after lead, then ranked candidates
   - Controls have bracket notation: [1], [2], [3]...
   - Ranked candidates have normal ranks: 4, 5, 6...

7. ENHANCED SELECTION SUMMARY (v8.8.1):
   - Terminal output shows breakdown of ranked candidates by track
   - Shows breakdown by hallmark (subfamily)
   - Shows breakdown by scaffold type (original vs universal)
   - All stats included in JSON summary
   - Global minimum per track enforced (10-20% each)
   - Minimum original scaffold enforced (20% of ranked)
   - Track 1 generation increased (12 probes vs 6)

8. HALLMARK PROTECTION (v8.8.2):
   - Motifs can no longer override hallmark positions (42, 49, 50, 52)
   - Position 49 filter added: must be E, Q, or K (not G which is VH-like)
   - Ensures all candidates have proper VHH hallmarks

9. UNIVERSAL SCAFFOLD SINGLE HALLMARK (v8.8.2):
   - Universal scaffold now uses ONE fixed hallmark (FERG by default)
   - Variation comes from FW mutations only, not hallmark changes
   - Original scaffold is distributed across all target hallmarks
   - This ensures clean comparison: Universal tests "standard VHH framework"
     while Original tests "which hallmark works best with my FRs"

INHERITED FROM v8.6:
- Tracks 1-3 are ranked lanes (not exempt)
- Per-bucket per-track quotas for selection
- Temperature scaling for probabilities
- Weighted choice for position 66
- Dedup within track

TRACK STRUCTURE:

  Track 0: Controls (exempt) - Hallmarks only
  Track 1: Single vernier probes (RANKED)
  Track 2: Archetype motifs 2-3 pos (RANKED)
  Track 3: Archetype motifs 3-5 pos (RANKED)
  Track 4: Expansive optimization (RANKED)

Usage:
  python vhh_designer_v8_8.py -i input.fasta \
      --rules analysis_rules_v7_all_positions.json \
      --archetypes analysis_vernier_archetypes_v7.json \
      --hallmark-db comprehensive_subfamily_analysis_imgt.xlsx \
      --target-hallmarks AUTO \
      --n-generate 200000 \
      --n-select 198

Author: Claude (Anthropic)
Date: January 2026
```

---

### `vhh_designer_v8_9.py`
**Version:** 8.9
**Path:** `active/analysis/vhh_designer_v8_9.py`
**Stats:** 6715 lines | 274.9KB | 100 functions | 12 classes
**Modified:** 2026-01-30 00:15

**Description:**
```
VHH Designer v8.9 - Restructured Track System with Capped Layers
================================================================

MAJOR UPDATES from v8.8.3:

1. SIMPLIFIED YQRL CONTROL STRUCTURE (v8.9):
   - 0a: GRAFT (0 mutations) - pure CDR graft
   - 0b: Hallmarks + IMGT2 (5 mutations) - ALWAYS include IMGT2 with hallmarks
   - 0c: Hallmarks + IMGT2 + ALL verniers (~20 mutations)
   - 0d+: Progressive framework in PAIRS (2 at a time), all positions >80% consensus
   - Removed separate "hallmarks only (4 mut)" - always include IMGT2

2. ALL HALLMARKS GET GRAFT FIRST (v8.9):
   - 0a: GRAFT (0 mutations) - always first for ALL hallmarks
   - 0b: Hallmarks + IMGT2 (5 mutations)
   - 0c: + ALL verniers (only if hallmark avg vernier consensus ≥60%)

3. TRACK 1 EXPANDED (v8.9):
   - Now includes ALL verniers, not just "load-bearing"
   - Expanded from max 12 to max 20 probes per hallmark
   - Tests each vernier position individually

4. TRACK 2 REDEFINED: "MINIMAL STRUCTURED PERTURBATIONS" (v8.9):
   - Base: Hallmarks + IMGT2 (IMGT2 now protected like hallmarks)
   - Add: Exactly ONE 2-3 position motif (weighted by family support)
   - Add: 0-K extra vernier mutations (K_EXTRA_VERNIERS_T2 = 2, capped)
   - NO framework tier sprinkling - keeps it interpretable

5. TRACK 3 REDEFINED: "MOTIF ATTRIBUTION" (v8.9):
   - Base: Hallmarks + IMGT2
   - Add: Exactly ONE 3-5 position motif (weighted)
   - Add: 0-K extra vernier mutations (K_EXTRA_VERNIERS_T3 = 2, max 3)
   - NO framework tier sprinkling
   - Goal: Interpretable - can analyze "motif X worked/didn't work"

6. TRACK 4 REDEFINED: "PRODUCTION OPTIMIZER" with CAPPED LAYERS (v8.9):
   - Full probabilistic layering for production candidates
   - Vernier layer: Bernoulli + cap (K_VERNIER_MAX = 8)
   - Framework tiers applied in order with individual caps:
     - Tier S (≥95%): K_S_MAX = 10
     - Tier A (90-95%): K_A_MAX = 6
     - Tier B (85-90%): K_B_MAX = 3
     - Tier C (75-85%): K_C_MAX = 2 (only when desperate)
   - CDR-conditional rules: capped at K_CDR_RULES_MAX = 2
   - Adaptive expansion when uniqueness drops

7. UNIVERSAL SCAFFOLD CONTROLS (v8.9):
   - 0a: Complete CDR graft into Universal (0 mutations)
   - 0b: Hallmarks + IMGT2 (5 mutations)
   - 0c: Hallmarks + IMGT2 + ALL verniers

8. IMGT2 NOW PROTECTED LIKE HALLMARKS (v8.9):
   - PROTECTED_POSITIONS = {2, 42, 49, 50, 52}
   - Motifs cannot override these positions

KEY DISTINCTION:
- Track 3 = "Motif attribution": one motif + minimal extras → interpretable
- Track 4 = "Production optimizer": full probabilistic layering → coverage

TRACK STRUCTURE (v8.9):

  ORIGINAL FRAMEWORK:
  Track 0: Controls (exempt)
    - 0a: GRAFT (0 mutations) - ALL hallmarks get this
    - 0b: Hallmarks + IMGT2 (5 mutations)
    - 0c: + ALL verniers (~20 mutations) - if avg vernier ≥60%
    - 0d+: YQRL only: + Progressive FW pairs (>80% consensus)
  Track 1: Single vernier probes - ALL verniers (RANKED)
  Track 2: Minimal structured perturbations - 1 motif + 0-2 verniers (RANKED)
  Track 3: Motif attribution - 1 motif (3-5 pos) + 0-3 verniers (RANKED)
  Track 4: Production optimizer - capped probabilistic layers (RANKED)

  UNIVERSAL FRAMEWORK:
  Track 0: Controls (exempt)
    - 0a: Complete CDR graft (0 mutations)
    - 0b: Hallmarks + IMGT2 (5 mutations)
    - 0c: + ALL verniers (~20 mutations)
  Track 1-4: Similar to Original but CDR-conditional focused

Usage:
  python vhh_designer_v8_9.py -i input.fasta \
      --rules analysis_rules_v7_all_positions.json \
      --archetypes analysis_vernier_archetypes_v7.json \
      --hallmark-db comprehensive_subfamily_analysis_imgt.xlsx \
      --target-hallmarks AUTO \
      --n-generate 200000 \
      --n-select 198

Author: Claude (Anthropic)
Date: January 2026
```

---

### `vhh_designer_v8_9_1.py`
**Version:** 8.9.1
**Path:** `active/analysis/vhh_designer_v8_9_1.py`
**Stats:** 7144 lines | 295.0KB | 101 functions | 12 classes
**Modified:** 2026-01-30 16:11

**Description:**
```
VHH Designer v8.9.1 - Restructured Track System with Capped Layers
===================================================================

v8.9.1 CRITICAL FIXES:
1. CYSTEINE PROTECTION (v8.9.1):
   - NEVER mutate FROM any cysteine position (disulfide bonds are essential)
   - Added CONSERVED_CYSTEINE_POSITIONS = {22, 23, 103, 104}
   - Three layers of protection:
     a) Position exclusion: Skip CONSERVED_CYSTEINE_POSITIONS in FRAMEWORK_CONSENSUS_TIERS loops
     b) Mutation filter: if scaffold_aa == 'C': continue (in all mutation generation)
     c) Builder safety: Filter in _build_*_candidate() catches any missed mutations
   - Added validate_canonical_disulfide() to reject candidates missing disulfide bond
   - Fixed bug where FRAMEWORK_CONSENSUS_TIERS['A'][22]='S' and ['B'][103]='Y' 
     would mutate Universal scaffold cysteines (C22→S, C103→Y)
   - NOTE: Cysteine positions NOT in EXCLUDE_FROM_FRAMEWORK because YQRL needs to
     mutate TO cysteine at 23/104 - protection is via individual checks instead

2. YQRL RIGID HALLMARK FIX (v8.9.1):
   - Added is_rigid check to _generate_original (was only in _generate_universal)
   - YQRL now correctly gets only Track 0 controls on Original scaffold

3. REDUCED YQRL FRAMEWORK STEPS (v8.9.1):
   - Reduced from 27 individual steps to 13 (5 pairs + 8 groups)
   - Maintains interpretability while reducing control explosion

4. HALLMARK DIVERSITY ENFORCEMENT (v8.9.1):
   - Added diversity enforcement in select_with_track_quotas
   - Prevents single hallmark (e.g., FERG) from dominating ranked candidates

MAJOR UPDATES from v8.8.3:

1. SIMPLIFIED YQRL CONTROL STRUCTURE (v8.9):
   - 0a: GRAFT (0 mutations) - pure CDR graft
   - 0b: Hallmarks + IMGT2 (5 mutations) - ALWAYS include IMGT2 with hallmarks
   - 0c: Hallmarks + IMGT2 + ALL verniers (~20 mutations)
   - 0d+: Progressive framework in PAIRS (2 at a time), all positions >80% consensus
   - Removed separate "hallmarks only (4 mut)" - always include IMGT2

2. ALL HALLMARKS GET GRAFT FIRST (v8.9):
   - 0a: GRAFT (0 mutations) - always first for ALL hallmarks
   - 0b: Hallmarks + IMGT2 (5 mutations)
   - 0c: + ALL verniers (only if hallmark avg vernier consensus ≥60%)

3. TRACK 1 EXPANDED (v8.9):
   - Now includes ALL verniers, not just "load-bearing"
   - Expanded from max 12 to max 20 probes per hallmark
   - Tests each vernier position individually

4. TRACK 2 REDEFINED: "MINIMAL STRUCTURED PERTURBATIONS" (v8.9):
   - Base: Hallmarks + IMGT2 (IMGT2 now protected like hallmarks)
   - Add: Exactly ONE 2-3 position motif (weighted by family support)
   - Add: 0-K extra vernier mutations (K_EXTRA_VERNIERS_T2 = 2, capped)
   - NO framework tier sprinkling - keeps it interpretable

5. TRACK 3 REDEFINED: "MOTIF ATTRIBUTION" (v8.9):
   - Base: Hallmarks + IMGT2
   - Add: Exactly ONE 3-5 position motif (weighted)
   - Add: 0-K extra vernier mutations (K_EXTRA_VERNIERS_T3 = 2, max 3)
   - NO framework tier sprinkling
   - Goal: Interpretable - can analyze "motif X worked/didn't work"

6. TRACK 4 REDEFINED: "PRODUCTION OPTIMIZER" with CAPPED LAYERS (v8.9):
   - Full probabilistic layering for production candidates
   - Vernier layer: Bernoulli + cap (K_VERNIER_MAX = 8)
   - Framework tiers applied in order with individual caps:
     - Tier S (≥95%): K_S_MAX = 10
     - Tier A (90-95%): K_A_MAX = 6
     - Tier B (85-90%): K_B_MAX = 3
     - Tier C (75-85%): K_C_MAX = 2 (only when desperate)
   - CDR-conditional rules: capped at K_CDR_RULES_MAX = 2
   - Adaptive expansion when uniqueness drops

7. UNIVERSAL SCAFFOLD CONTROLS (v8.9):
   - 0a: Complete CDR graft into Universal (0 mutations)
   - 0b: Hallmarks + IMGT2 (5 mutations)
   - 0c: Hallmarks + IMGT2 + ALL verniers

8. IMGT2 NOW PROTECTED LIKE HALLMARKS (v8.9):
   - PROTECTED_POSITIONS = {2, 42, 49, 50, 52}
   - Motifs cannot override these positions

KEY DISTINCTION:
- Track 3 = "Motif attribution": one motif + minimal extras → interpretable
- Track 4 = "Production optimizer": full probabilistic layering → coverage

TRACK STRUCTURE (v8.9):

  ORIGINAL FRAMEWORK:
  Track 0: Controls (exempt)
    - 0a: GRAFT (0 mutations) - ALL hallmarks get this
    - 0b: Hallmarks + IMGT2 (5 mutations)
    - 0c: + ALL verniers (~20 mutations) - if avg vernier ≥60%
    - 0d+: YQRL only: + Progressive FW pairs (>80% consensus)
  Track 1: Single vernier probes - ALL verniers (RANKED)
  Track 2: Minimal structured perturbations - 1 motif + 0-2 verniers (RANKED)
  Track 3: Motif attribution - 1 motif (3-5 pos) + 0-3 verniers (RANKED)
  Track 4: Production optimizer - capped probabilistic layers (RANKED)

  UNIVERSAL FRAMEWORK:
  Track 0: Controls (exempt)
    - 0a: Complete CDR graft (0 mutations)
    - 0b: Hallmarks + IMGT2 (5 mutations)
    - 0c: + ALL verniers (~20 mutations)
  Track 1-4: Similar to Original but CDR-conditional focused

Usage:
  python vhh_designer_v8_9.py -i input.fasta \
      --rules analysis_rules_v7_all_positions.json \
      --archetypes analysis_vernier_archetypes_v7.json \
      --hallmark-db comprehensive_subfamily_analysis_imgt.xlsx \
      --target-hallmarks AUTO \
      --n-generate 200000 \
      --n-select 198

Author: Claude (Anthropic)
Date: January 2026
```

---

### `vhh_designer_v8_9_2.py`
**Version:** 8.9.1
**Path:** `active/analysis/vhh_designer_v8_9_2.py`
**Stats:** 7270 lines | 300.5KB | 101 functions | 12 classes
**Modified:** 2026-01-31 00:00

**Description:**
```
VHH Designer v8.9.1 - Restructured Track System with Capped Layers
===================================================================

v8.9.1 CRITICAL FIXES:
1. CYSTEINE PROTECTION (v8.9.1):
   - NEVER mutate FROM any cysteine position (disulfide bonds are essential)
   - Added CONSERVED_CYSTEINE_POSITIONS = {22, 23, 103, 104}
   - Three layers of protection:
     a) Position exclusion: Skip CONSERVED_CYSTEINE_POSITIONS in FRAMEWORK_CONSENSUS_TIERS loops
     b) Mutation filter: if scaffold_aa == 'C': continue (in all mutation generation)
     c) Builder safety: Filter in _build_*_candidate() catches any missed mutations
   - Added validate_canonical_disulfide() to reject candidates missing disulfide bond
   - Fixed bug where FRAMEWORK_CONSENSUS_TIERS['A'][22]='S' and ['B'][103]='Y' 
     would mutate Universal scaffold cysteines (C22→S, C103→Y)
   - NOTE: Cysteine positions NOT in EXCLUDE_FROM_FRAMEWORK because YQRL needs to
     mutate TO cysteine at 23/104 - protection is via individual checks instead

2. YQRL RIGID HALLMARK FIX (v8.9.1):
   - Added is_rigid check to _generate_original (was only in _generate_universal)
   - YQRL now correctly gets only Track 0 controls on Original scaffold

3. REDUCED YQRL FRAMEWORK STEPS (v8.9.1):
   - Reduced from 27 individual steps to 13 (5 pairs + 8 groups)
   - Maintains interpretability while reducing control explosion

4. HALLMARK DIVERSITY ENFORCEMENT (v8.9.1):
   - Added diversity enforcement in select_with_track_quotas
   - Prevents single hallmark (e.g., FERG) from dominating ranked candidates

MAJOR UPDATES from v8.8.3:

1. SIMPLIFIED YQRL CONTROL STRUCTURE (v8.9):
   - 0a: GRAFT (0 mutations) - pure CDR graft
   - 0b: Hallmarks + IMGT2 (5 mutations) - ALWAYS include IMGT2 with hallmarks
   - 0c: Hallmarks + IMGT2 + ALL verniers (~20 mutations)
   - 0d+: Progressive framework in PAIRS (2 at a time), all positions >80% consensus
   - Removed separate "hallmarks only (4 mut)" - always include IMGT2

2. ALL HALLMARKS GET GRAFT FIRST (v8.9):
   - 0a: GRAFT (0 mutations) - always first for ALL hallmarks
   - 0b: Hallmarks + IMGT2 (5 mutations)
   - 0c: + ALL verniers (only if hallmark avg vernier consensus ≥60%)

3. TRACK 1 EXPANDED (v8.9):
   - Now includes ALL verniers, not just "load-bearing"
   - Expanded from max 12 to max 20 probes per hallmark
   - Tests each vernier position individually

4. TRACK 2 REDEFINED: "MINIMAL STRUCTURED PERTURBATIONS" (v8.9):
   - Base: Hallmarks + IMGT2 (IMGT2 now protected like hallmarks)
   - Add: Exactly ONE 2-3 position motif (weighted by family support)
   - Add: 0-K extra vernier mutations (K_EXTRA_VERNIERS_T2 = 2, capped)
   - NO framework tier sprinkling - keeps it interpretable

5. TRACK 3 REDEFINED: "MOTIF ATTRIBUTION" (v8.9):
   - Base: Hallmarks + IMGT2
   - Add: Exactly ONE 3-5 position motif (weighted)
   - Add: 0-K extra vernier mutations (K_EXTRA_VERNIERS_T3 = 2, max 3)
   - NO framework tier sprinkling
   - Goal: Interpretable - can analyze "motif X worked/didn't work"

6. TRACK 4 REDEFINED: "PRODUCTION OPTIMIZER" with CAPPED LAYERS (v8.9):
   - Full probabilistic layering for production candidates
   - Vernier layer: Bernoulli + cap (K_VERNIER_MAX = 8)
   - Framework tiers applied in order with individual caps:
     - Tier S (≥95%): K_S_MAX = 10
     - Tier A (90-95%): K_A_MAX = 6
     - Tier B (85-90%): K_B_MAX = 3
     - Tier C (75-85%): K_C_MAX = 2 (only when desperate)
   - CDR-conditional rules: capped at K_CDR_RULES_MAX = 2
   - Adaptive expansion when uniqueness drops

7. UNIVERSAL SCAFFOLD CONTROLS (v8.9):
   - 0a: Complete CDR graft into Universal (0 mutations)
   - 0b: Hallmarks + IMGT2 (5 mutations)
   - 0c: Hallmarks + IMGT2 + ALL verniers

8. IMGT2 NOW PROTECTED LIKE HALLMARKS (v8.9):
   - PROTECTED_POSITIONS = {2, 42, 49, 50, 52}
   - Motifs cannot override these positions

KEY DISTINCTION:
- Track 3 = "Motif attribution": one motif + minimal extras → interpretable
- Track 4 = "Production optimizer": full probabilistic layering → coverage

TRACK STRUCTURE (v8.9):

  ORIGINAL FRAMEWORK:
  Track 0: Controls (exempt)
    - 0a: GRAFT (0 mutations) - ALL hallmarks get this
    - 0b: Hallmarks + IMGT2 (5 mutations)
    - 0c: + ALL verniers (~20 mutations) - if avg vernier ≥60%
    - 0d+: YQRL only: + Progressive FW pairs (>80% consensus)
  Track 1: Single vernier probes - ALL verniers (RANKED)
  Track 2: Minimal structured perturbations - 1 motif + 0-2 verniers (RANKED)
  Track 3: Motif attribution - 1 motif (3-5 pos) + 0-3 verniers (RANKED)
  Track 4: Production optimizer - capped probabilistic layers (RANKED)

  UNIVERSAL FRAMEWORK:
  Track 0: Controls (exempt)
    - 0a: Complete CDR graft (0 mutations)
    - 0b: Hallmarks + IMGT2 (5 mutations)
    - 0c: + ALL verniers (~20 mutations)
  Track 1-4: Similar to Original but CDR-conditional focused

Usage:
  python vhh_designer_v8_9.py -i input.fasta \
      --rules analysis_rules_v7_all_positions.json \
      --archetypes analysis_vernier_archetypes_v7.json \
      --hallmark-db comprehensive_subfamily_analysis_imgt.xlsx \
      --target-hallmarks AUTO \
      --n-generate 200000 \
      --n-select 198

Author: Claude (Anthropic)
Date: January 2026
```

---

### `vhh_designer_v8_9_3.py`
**Version:** 8.9.1
**Path:** `active/analysis/vhh_designer_v8_9_3.py`
**Stats:** 7556 lines | 312.3KB | 103 functions | 12 classes
**Modified:** 2026-02-01 00:41

**Description:**
```
VHH Designer v8.9.1 - Restructured Track System with Capped Layers
===================================================================

v8.9.1 CRITICAL FIXES:
1. CYSTEINE PROTECTION (v8.9.1):
   - NEVER mutate FROM any cysteine position (disulfide bonds are essential)
   - Added CONSERVED_CYSTEINE_POSITIONS = {22, 23, 103, 104}
   - Three layers of protection:
     a) Position exclusion: Skip CONSERVED_CYSTEINE_POSITIONS in FRAMEWORK_CONSENSUS_TIERS loops
     b) Mutation filter: if scaffold_aa == 'C': continue (in all mutation generation)
     c) Builder safety: Filter in _build_*_candidate() catches any missed mutations
   - Added validate_canonical_disulfide() to reject candidates missing disulfide bond
   - Fixed bug where FRAMEWORK_CONSENSUS_TIERS['A'][22]='S' and ['B'][103]='Y' 
     would mutate Universal scaffold cysteines (C22→S, C103→Y)
   - NOTE: Cysteine positions NOT in EXCLUDE_FROM_FRAMEWORK because YQRL needs to
     mutate TO cysteine at 23/104 - protection is via individual checks instead

2. YQRL RIGID HALLMARK FIX (v8.9.1):
   - Added is_rigid check to _generate_original (was only in _generate_universal)
   - YQRL now correctly gets only Track 0 controls on Original scaffold

3. REDUCED YQRL FRAMEWORK STEPS (v8.9.1):
   - Reduced from 27 individual steps to 13 (5 pairs + 8 groups)
   - Maintains interpretability while reducing control explosion

4. HALLMARK DIVERSITY ENFORCEMENT (v8.9.1):
   - Added diversity enforcement in select_with_track_quotas
   - Prevents single hallmark (e.g., FERG) from dominating ranked candidates

MAJOR UPDATES from v8.8.3:

1. SIMPLIFIED YQRL CONTROL STRUCTURE (v8.9):
   - 0a: GRAFT (0 mutations) - pure CDR graft
   - 0b: Hallmarks + IMGT2 (5 mutations) - ALWAYS include IMGT2 with hallmarks
   - 0c: Hallmarks + IMGT2 + ALL verniers (~20 mutations)
   - 0d+: Progressive framework in PAIRS (2 at a time), all positions >80% consensus
   - Removed separate "hallmarks only (4 mut)" - always include IMGT2

2. ALL HALLMARKS GET GRAFT FIRST (v8.9):
   - 0a: GRAFT (0 mutations) - always first for ALL hallmarks
   - 0b: Hallmarks + IMGT2 (5 mutations)
   - 0c: + ALL verniers (only if hallmark avg vernier consensus ≥60%)

3. TRACK 1 EXPANDED (v8.9):
   - Now includes ALL verniers, not just "load-bearing"
   - Expanded from max 12 to max 20 probes per hallmark
   - Tests each vernier position individually

4. TRACK 2 REDEFINED: "MINIMAL STRUCTURED PERTURBATIONS" (v8.9):
   - Base: Hallmarks + IMGT2 (IMGT2 now protected like hallmarks)
   - Add: Exactly ONE 2-3 position motif (weighted by family support)
   - Add: 0-K extra vernier mutations (K_EXTRA_VERNIERS_T2 = 2, capped)
   - NO framework tier sprinkling - keeps it interpretable

5. TRACK 3 REDEFINED: "MOTIF ATTRIBUTION" (v8.9):
   - Base: Hallmarks + IMGT2
   - Add: Exactly ONE 3-5 position motif (weighted)
   - Add: 0-K extra vernier mutations (K_EXTRA_VERNIERS_T3 = 2, max 3)
   - NO framework tier sprinkling
   - Goal: Interpretable - can analyze "motif X worked/didn't work"

6. TRACK 4 REDEFINED: "PRODUCTION OPTIMIZER" with CAPPED LAYERS (v8.9):
   - Full probabilistic layering for production candidates
   - Vernier layer: Bernoulli + cap (K_VERNIER_MAX = 8)
   - Framework tiers applied in order with individual caps:
     - Tier S (≥95%): K_S_MAX = 10
     - Tier A (90-95%): K_A_MAX = 6
     - Tier B (85-90%): K_B_MAX = 3
     - Tier C (75-85%): K_C_MAX = 2 (only when desperate)
   - CDR-conditional rules: capped at K_CDR_RULES_MAX = 2
   - Adaptive expansion when uniqueness drops

7. UNIVERSAL SCAFFOLD CONTROLS (v8.9):
   - 0a: Complete CDR graft into Universal (0 mutations)
   - 0b: Hallmarks + IMGT2 (5 mutations)
   - 0c: Hallmarks + IMGT2 + ALL verniers

8. IMGT2 NOW PROTECTED LIKE HALLMARKS (v8.9):
   - PROTECTED_POSITIONS = {2, 42, 49, 50, 52}
   - Motifs cannot override these positions

KEY DISTINCTION:
- Track 3 = "Motif attribution": one motif + minimal extras → interpretable
- Track 4 = "Production optimizer": full probabilistic layering → coverage

TRACK STRUCTURE (v8.9):

  ORIGINAL FRAMEWORK:
  Track 0: Controls (exempt)
    - 0a: GRAFT (0 mutations) - ALL hallmarks get this
    - 0b: Hallmarks + IMGT2 (5 mutations)
    - 0c: + ALL verniers (~20 mutations) - if avg vernier ≥60%
    - 0d+: YQRL only: + Progressive FW pairs (>80% consensus)
  Track 1: Single vernier probes - ALL verniers (RANKED)
  Track 2: Minimal structured perturbations - 1 motif + 0-2 verniers (RANKED)
  Track 3: Motif attribution - 1 motif (3-5 pos) + 0-3 verniers (RANKED)
  Track 4: Production optimizer - capped probabilistic layers (RANKED)

  UNIVERSAL FRAMEWORK:
  Track 0: Controls (exempt)
    - 0a: Complete CDR graft (0 mutations)
    - 0b: Hallmarks + IMGT2 (5 mutations)
    - 0c: + ALL verniers (~20 mutations)
  Track 1-4: Similar to Original but CDR-conditional focused

Usage:
  python vhh_designer_v8_9.py -i input.fasta \
      --rules analysis_rules_v7_all_positions.json \
      --archetypes analysis_vernier_archetypes_v7.json \
      --hallmark-db comprehensive_subfamily_analysis_imgt.xlsx \
      --target-hallmarks AUTO \
      --n-generate 200000 \
      --n-select 198

Author: Claude (Anthropic)
Date: January 2026
```

---

### `vhh_designer_v8_9_4.py`
**Version:** 8.9.1
**Path:** `active/analysis/vhh_designer_v8_9_4.py`
**Stats:** 7563 lines | 312.7KB | 103 functions | 12 classes
**Modified:** 2026-02-01 15:51

**Description:**
```
VHH Designer v8.9.1 - Restructured Track System with Capped Layers
===================================================================

v8.9.1 CRITICAL FIXES:
1. CYSTEINE PROTECTION (v8.9.1):
   - NEVER mutate FROM any cysteine position (disulfide bonds are essential)
   - Added CONSERVED_CYSTEINE_POSITIONS = {22, 23, 103, 104}
   - Three layers of protection:
     a) Position exclusion: Skip CONSERVED_CYSTEINE_POSITIONS in FRAMEWORK_CONSENSUS_TIERS loops
     b) Mutation filter: if scaffold_aa == 'C': continue (in all mutation generation)
     c) Builder safety: Filter in _build_*_candidate() catches any missed mutations
   - Added validate_canonical_disulfide() to reject candidates missing disulfide bond
   - Fixed bug where FRAMEWORK_CONSENSUS_TIERS['A'][22]='S' and ['B'][103]='Y' 
     would mutate Universal scaffold cysteines (C22→S, C103→Y)
   - NOTE: Cysteine positions NOT in EXCLUDE_FROM_FRAMEWORK because YQRL needs to
     mutate TO cysteine at 23/104 - protection is via individual checks instead

2. YQRL RIGID HALLMARK FIX (v8.9.1):
   - Added is_rigid check to _generate_original (was only in _generate_universal)
   - YQRL now correctly gets only Track 0 controls on Original scaffold

3. REDUCED YQRL FRAMEWORK STEPS (v8.9.1):
   - Reduced from 27 individual steps to 13 (5 pairs + 8 groups)
   - Maintains interpretability while reducing control explosion

4. HALLMARK DIVERSITY ENFORCEMENT (v8.9.1):
   - Added diversity enforcement in select_with_track_quotas
   - Prevents single hallmark (e.g., FERG) from dominating ranked candidates

MAJOR UPDATES from v8.8.3:

1. SIMPLIFIED YQRL CONTROL STRUCTURE (v8.9):
   - 0a: GRAFT (0 mutations) - pure CDR graft
   - 0b: Hallmarks + IMGT2 (5 mutations) - ALWAYS include IMGT2 with hallmarks
   - 0c: Hallmarks + IMGT2 + ALL verniers (~20 mutations)
   - 0d+: Progressive framework in PAIRS (2 at a time), all positions >80% consensus
   - Removed separate "hallmarks only (4 mut)" - always include IMGT2

2. ALL HALLMARKS GET GRAFT FIRST (v8.9):
   - 0a: GRAFT (0 mutations) - always first for ALL hallmarks
   - 0b: Hallmarks + IMGT2 (5 mutations)
   - 0c: + ALL verniers (only if hallmark avg vernier consensus ≥60%)

3. TRACK 1 EXPANDED (v8.9):
   - Now includes ALL verniers, not just "load-bearing"
   - Expanded from max 12 to max 20 probes per hallmark
   - Tests each vernier position individually

4. TRACK 2 REDEFINED: "MINIMAL STRUCTURED PERTURBATIONS" (v8.9):
   - Base: Hallmarks + IMGT2 (IMGT2 now protected like hallmarks)
   - Add: Exactly ONE 2-3 position motif (weighted by family support)
   - Add: 0-K extra vernier mutations (K_EXTRA_VERNIERS_T2 = 2, capped)
   - NO framework tier sprinkling - keeps it interpretable

5. TRACK 3 REDEFINED: "MOTIF ATTRIBUTION" (v8.9):
   - Base: Hallmarks + IMGT2
   - Add: Exactly ONE 3-5 position motif (weighted)
   - Add: 0-K extra vernier mutations (K_EXTRA_VERNIERS_T3 = 2, max 3)
   - NO framework tier sprinkling
   - Goal: Interpretable - can analyze "motif X worked/didn't work"

6. TRACK 4 REDEFINED: "PRODUCTION OPTIMIZER" with CAPPED LAYERS (v8.9):
   - Full probabilistic layering for production candidates
   - Vernier layer: Bernoulli + cap (K_VERNIER_MAX = 8)
   - Framework tiers applied in order with individual caps:
     - Tier S (≥95%): K_S_MAX = 10
     - Tier A (90-95%): K_A_MAX = 6
     - Tier B (85-90%): K_B_MAX = 3
     - Tier C (75-85%): K_C_MAX = 2 (only when desperate)
   - CDR-conditional rules: capped at K_CDR_RULES_MAX = 2
   - Adaptive expansion when uniqueness drops

7. UNIVERSAL SCAFFOLD CONTROLS (v8.9):
   - 0a: Complete CDR graft into Universal (0 mutations)
   - 0b: Hallmarks + IMGT2 (5 mutations)
   - 0c: Hallmarks + IMGT2 + ALL verniers

8. IMGT2 NOW PROTECTED LIKE HALLMARKS (v8.9):
   - PROTECTED_POSITIONS = {2, 42, 49, 50, 52}
   - Motifs cannot override these positions

KEY DISTINCTION:
- Track 3 = "Motif attribution": one motif + minimal extras → interpretable
- Track 4 = "Production optimizer": full probabilistic layering → coverage

TRACK STRUCTURE (v8.9):

  ORIGINAL FRAMEWORK:
  Track 0: Controls (exempt)
    - 0a: GRAFT (0 mutations) - ALL hallmarks get this
    - 0b: Hallmarks + IMGT2 (5 mutations)
    - 0c: + ALL verniers (~20 mutations) - if avg vernier ≥60%
    - 0d+: YQRL only: + Progressive FW pairs (>80% consensus)
  Track 1: Single vernier probes - ALL verniers (RANKED)
  Track 2: Minimal structured perturbations - 1 motif + 0-2 verniers (RANKED)
  Track 3: Motif attribution - 1 motif (3-5 pos) + 0-3 verniers (RANKED)
  Track 4: Production optimizer - capped probabilistic layers (RANKED)

  UNIVERSAL FRAMEWORK:
  Track 0: Controls (exempt)
    - 0a: Complete CDR graft (0 mutations)
    - 0b: Hallmarks + IMGT2 (5 mutations)
    - 0c: + ALL verniers (~20 mutations)
  Track 1-4: Similar to Original but CDR-conditional focused

Usage:
  python vhh_designer_v8_9.py -i input.fasta \
      --rules analysis_rules_v7_all_positions.json \
      --archetypes analysis_vernier_archetypes_v7.json \
      --hallmark-db comprehensive_subfamily_analysis_imgt.xlsx \
      --target-hallmarks AUTO \
      --n-generate 200000 \
      --n-select 198

Author: Claude (Anthropic)
Date: January 2026
```

---

### `vhh_designer_v8_9_4_1.py`
**Version:** 8.9.1
**Path:** `active/analysis/vhh_designer_v8_9_4_1.py`
**Stats:** 7599 lines | 314.6KB | 103 functions | 12 classes
**Modified:** 2026-02-01 20:02

**Description:**
```
VHH Designer v8.9.1 - Restructured Track System with Capped Layers
===================================================================

v8.9.1 CRITICAL FIXES:
1. CYSTEINE PROTECTION (v8.9.1):
   - NEVER mutate FROM any cysteine position (disulfide bonds are essential)
   - Added CONSERVED_CYSTEINE_POSITIONS = {22, 23, 103, 104}
   - Three layers of protection:
     a) Position exclusion: Skip CONSERVED_CYSTEINE_POSITIONS in FRAMEWORK_CONSENSUS_TIERS loops
     b) Mutation filter: if scaffold_aa == 'C': continue (in all mutation generation)
     c) Builder safety: Filter in _build_*_candidate() catches any missed mutations
   - Added validate_canonical_disulfide() to reject candidates missing disulfide bond
   - Fixed bug where FRAMEWORK_CONSENSUS_TIERS['A'][22]='S' and ['B'][103]='Y' 
     would mutate Universal scaffold cysteines (C22→S, C103→Y)
   - NOTE: Cysteine positions NOT in EXCLUDE_FROM_FRAMEWORK because YQRL needs to
     mutate TO cysteine at 23/104 - protection is via individual checks instead

2. YQRL RIGID HALLMARK FIX (v8.9.1):
   - Added is_rigid check to _generate_original (was only in _generate_universal)
   - YQRL now correctly gets only Track 0 controls on Original scaffold

3. REDUCED YQRL FRAMEWORK STEPS (v8.9.1):
   - Reduced from 27 individual steps to 13 (5 pairs + 8 groups)
   - Maintains interpretability while reducing control explosion

4. HALLMARK DIVERSITY ENFORCEMENT (v8.9.1):
   - Added diversity enforcement in select_with_track_quotas
   - Prevents single hallmark (e.g., FERG) from dominating ranked candidates

MAJOR UPDATES from v8.8.3:

1. SIMPLIFIED YQRL CONTROL STRUCTURE (v8.9):
   - 0a: GRAFT (0 mutations) - pure CDR graft
   - 0b: Hallmarks + IMGT2 (5 mutations) - ALWAYS include IMGT2 with hallmarks
   - 0c: Hallmarks + IMGT2 + ALL verniers (~20 mutations)
   - 0d+: Progressive framework in PAIRS (2 at a time), all positions >80% consensus
   - Removed separate "hallmarks only (4 mut)" - always include IMGT2

2. ALL HALLMARKS GET GRAFT FIRST (v8.9):
   - 0a: GRAFT (0 mutations) - always first for ALL hallmarks
   - 0b: Hallmarks + IMGT2 (5 mutations)
   - 0c: + ALL verniers (only if hallmark avg vernier consensus ≥60%)

3. TRACK 1 EXPANDED (v8.9):
   - Now includes ALL verniers, not just "load-bearing"
   - Expanded from max 12 to max 20 probes per hallmark
   - Tests each vernier position individually

4. TRACK 2 REDEFINED: "MINIMAL STRUCTURED PERTURBATIONS" (v8.9):
   - Base: Hallmarks + IMGT2 (IMGT2 now protected like hallmarks)
   - Add: Exactly ONE 2-3 position motif (weighted by family support)
   - Add: 0-K extra vernier mutations (K_EXTRA_VERNIERS_T2 = 2, capped)
   - NO framework tier sprinkling - keeps it interpretable

5. TRACK 3 REDEFINED: "MOTIF ATTRIBUTION" (v8.9):
   - Base: Hallmarks + IMGT2
   - Add: Exactly ONE 3-5 position motif (weighted)
   - Add: 0-K extra vernier mutations (K_EXTRA_VERNIERS_T3 = 2, max 3)
   - NO framework tier sprinkling
   - Goal: Interpretable - can analyze "motif X worked/didn't work"

6. TRACK 4 REDEFINED: "PRODUCTION OPTIMIZER" with CAPPED LAYERS (v8.9):
   - Full probabilistic layering for production candidates
   - Vernier layer: Bernoulli + cap (K_VERNIER_MAX = 8)
   - Framework tiers applied in order with individual caps:
     - Tier S (≥95%): K_S_MAX = 10
     - Tier A (90-95%): K_A_MAX = 6
     - Tier B (85-90%): K_B_MAX = 3
     - Tier C (75-85%): K_C_MAX = 2 (only when desperate)
   - CDR-conditional rules: capped at K_CDR_RULES_MAX = 2
   - Adaptive expansion when uniqueness drops

7. UNIVERSAL SCAFFOLD CONTROLS (v8.9):
   - 0a: Complete CDR graft into Universal (0 mutations)
   - 0b: Hallmarks + IMGT2 (5 mutations)
   - 0c: Hallmarks + IMGT2 + ALL verniers

8. IMGT2 NOW PROTECTED LIKE HALLMARKS (v8.9):
   - PROTECTED_POSITIONS = {2, 42, 49, 50, 52}
   - Motifs cannot override these positions

KEY DISTINCTION:
- Track 3 = "Motif attribution": one motif + minimal extras → interpretable
- Track 4 = "Production optimizer": full probabilistic layering → coverage

TRACK STRUCTURE (v8.9):

  ORIGINAL FRAMEWORK:
  Track 0: Controls (exempt)
    - 0a: GRAFT (0 mutations) - ALL hallmarks get this
    - 0b: Hallmarks + IMGT2 (5 mutations)
    - 0c: + ALL verniers (~20 mutations) - if avg vernier ≥60%
    - 0d+: YQRL only: + Progressive FW pairs (>80% consensus)
  Track 1: Single vernier probes - ALL verniers (RANKED)
  Track 2: Minimal structured perturbations - 1 motif + 0-2 verniers (RANKED)
  Track 3: Motif attribution - 1 motif (3-5 pos) + 0-3 verniers (RANKED)
  Track 4: Production optimizer - capped probabilistic layers (RANKED)

  UNIVERSAL FRAMEWORK:
  Track 0: Controls (exempt)
    - 0a: Complete CDR graft (0 mutations)
    - 0b: Hallmarks + IMGT2 (5 mutations)
    - 0c: + ALL verniers (~20 mutations)
  Track 1-4: Similar to Original but CDR-conditional focused

Usage:
  python vhh_designer_v8_9.py -i input.fasta \
      --rules analysis_rules_v7_all_positions.json \
      --archetypes analysis_vernier_archetypes_v7.json \
      --hallmark-db comprehensive_subfamily_analysis_imgt.xlsx \
      --target-hallmarks AUTO \
      --n-generate 200000 \
      --n-select 198

Author: Claude (Anthropic)
Date: January 2026
```

---

### `vhh_designer_v9_0.py`
**Version:** 9.0
**Path:** `active/analysis/vhh_designer_v9_0.py`
**Stats:** 7652 lines | 317.1KB | 103 functions | 12 classes
**Modified:** 2026-02-01 21:42

**Description:**
```
VHH Designer v9.0 - Corrected Scoring & Multi-Hallmark Pool
============================================================

v9.0 MAJOR CHANGES:

1. **CRITICAL: SCORING DIRECTION FIX**
   - FIXED: Combined_score now correctly ranks better sequences higher
   - Was inverted: higher rule compliance → WORSE rank (WRONG!)
   - Now correct: higher rule compliance → BETTER rank
   - ESM2 loss: lower loss (more natural) → better rank
   - pLDDT: higher pLDDT (better structure) → better rank  
   - Rules: higher compliance → better rank
   - This was a fundamental bug - previous versions ranked the WORST candidates first!

2. DEFAULT HALLMARK POOL:
   - New DEFAULT_HALLMARK_POOL with 7 curated hallmarks
   - FERG, FERF, YQRL, FERA, YERL, FKRG, FQRA
   - --target-hallmarks now defaults to 'DEFAULT' (was 'FERG')
   - FQRA added for 100% vernier match to typical VH inputs

3. PERFORMANCE FIX - O(n²) BUG:
   - Fixed to_dataframe() hanging on large candidate pools (85K+)
   - Was: O(n) scan per candidate = O(n²) total = 7.2B comparisons
   - Now: O(1) dict lookup per candidate = O(n) total

4. _all.csv FIX:
   - Now saves ALL scored candidates before selection
   - Copies candidates list immediately after scoring

5. FILE NAMING:
   - Includes n-generate and temperature: v9_0_n100000_t08_M69
   - Format: {datetime}_{version}_n{n_generate}_t{temp}_{seq_name}

v8.9.3 CHANGES:
- Cross-track sequence deduplication (keeps higher track version)
- MSA-compatible ID format (M69_Orig_FERG_T4R_HVn_r001_s152)
- Hallmark-aware track minimums
- Hallmark cap enforcement (35% max per hallmark)
- 6-phase selection pipeline

v8.9.1 CRITICAL FIXES:
1. CYSTEINE PROTECTION (v8.9.1):
   - NEVER mutate FROM any cysteine position (disulfide bonds are essential)
   - Added CONSERVED_CYSTEINE_POSITIONS = {22, 23, 103, 104}
   - Three layers of protection:
     a) Position exclusion: Skip CONSERVED_CYSTEINE_POSITIONS in FRAMEWORK_CONSENSUS_TIERS loops
     b) Mutation filter: if scaffold_aa == 'C': continue (in all mutation generation)
     c) Builder safety: Filter in _build_*_candidate() catches any missed mutations
   - Added validate_canonical_disulfide() to reject candidates missing disulfide bond
   - Fixed bug where FRAMEWORK_CONSENSUS_TIERS['A'][22]='S' and ['B'][103]='Y' 
     would mutate Universal scaffold cysteines (C22→S, C103→Y)
   - NOTE: Cysteine positions NOT in EXCLUDE_FROM_FRAMEWORK because YQRL needs to
     mutate TO cysteine at 23/104 - protection is via individual checks instead

2. YQRL RIGID HALLMARK FIX (v8.9.1):
   - Added is_rigid check to _generate_original (was only in _generate_universal)
   - YQRL now correctly gets only Track 0 controls on Original scaffold

3. REDUCED YQRL FRAMEWORK STEPS (v8.9.1):
   - Reduced from 27 individual steps to 13 (5 pairs + 8 groups)
   - Maintains interpretability while reducing control explosion

4. HALLMARK DIVERSITY ENFORCEMENT (v8.9.1):
   - Added diversity enforcement in select_with_track_quotas
   - Prevents single hallmark (e.g., FERG) from dominating ranked candidates

MAJOR UPDATES from v8.8.3:

1. SIMPLIFIED YQRL CONTROL STRUCTURE (v8.9):
   - 0a: GRAFT (0 mutations) - pure CDR graft
   - 0b: Hallmarks + IMGT2 (5 mutations) - ALWAYS include IMGT2 with hallmarks
   - 0c: Hallmarks + IMGT2 + ALL verniers (~20 mutations)
   - 0d+: Progressive framework in PAIRS (2 at a time), all positions >80% consensus
   - Removed separate "hallmarks only (4 mut)" - always include IMGT2

2. ALL HALLMARKS GET GRAFT FIRST (v8.9):
   - 0a: GRAFT (0 mutations) - always first for ALL hallmarks
   - 0b: Hallmarks + IMGT2 (5 mutations)
   - 0c: + ALL verniers (only if hallmark avg vernier consensus ≥60%)

3. TRACK 1 EXPANDED (v8.9):
   - Now includes ALL verniers, not just "load-bearing"
   - Expanded from max 12 to max 20 probes per hallmark
   - Tests each vernier position individually

4. TRACK 2 REDEFINED: "MINIMAL STRUCTURED PERTURBATIONS" (v8.9):
   - Base: Hallmarks + IMGT2 (IMGT2 now protected like hallmarks)
   - Add: Exactly ONE 2-3 position motif (weighted by family support)
   - Add: 0-K extra vernier mutations (K_EXTRA_VERNIERS_T2 = 2, capped)
   - NO framework tier sprinkling - keeps it interpretable

5. TRACK 3 REDEFINED: "MOTIF ATTRIBUTION" (v8.9):
   - Base: Hallmarks + IMGT2
   - Add: Exactly ONE 3-5 position motif (weighted)
   - Add: 0-K extra vernier mutations (K_EXTRA_VERNIERS_T3 = 2, max 3)
   - NO framework tier sprinkling
   - Goal: Interpretable - can analyze "motif X worked/didn't work"

6. TRACK 4 REDEFINED: "PRODUCTION OPTIMIZER" with CAPPED LAYERS (v8.9):
   - Full probabilistic layering for production candidates
   - Vernier layer: Bernoulli + cap (K_VERNIER_MAX = 8)
   - Framework tiers applied in order with individual caps:
     - Tier S (≥95%): K_S_MAX = 10
     - Tier A (90-95%): K_A_MAX = 6
     - Tier B (85-90%): K_B_MAX = 3
     - Tier C (75-85%): K_C_MAX = 2 (only when desperate)
   - CDR-conditional rules: capped at K_CDR_RULES_MAX = 2
   - Adaptive expansion when uniqueness drops

7. UNIVERSAL SCAFFOLD CONTROLS (v8.9):
   - 0a: Complete CDR graft into Universal (0 mutations)
   - 0b: Hallmarks + IMGT2 (5 mutations)
   - 0c: Hallmarks + IMGT2 + ALL verniers

8. IMGT2 NOW PROTECTED LIKE HALLMARKS (v8.9):
   - PROTECTED_POSITIONS = {2, 42, 49, 50, 52}
   - Motifs cannot override these positions

KEY DISTINCTION:
- Track 3 = "Motif attribution": one motif + minimal extras → interpretable
- Track 4 = "Production optimizer": full probabilistic layering → coverage

TRACK STRUCTURE (v8.9):

  ORIGINAL FRAMEWORK:
  Track 0: Controls (exempt)
    - 0a: GRAFT (0 mutations) - ALL hallmarks get this
    - 0b: Hallmarks + IMGT2 (5 mutations)
    - 0c: + ALL verniers (~20 mutations) - if avg vernier ≥60%
    - 0d+: YQRL only: + Progressive FW pairs (>80% consensus)
  Track 1: Single vernier probes - ALL verniers (RANKED)
  Track 2: Minimal structured perturbations - 1 motif + 0-2 verniers (RANKED)
  Track 3: Motif attribution - 1 motif (3-5 pos) + 0-3 verniers (RANKED)
  Track 4: Production optimizer - capped probabilistic layers (RANKED)

  UNIVERSAL FRAMEWORK:
  Track 0: Controls (exempt)
    - 0a: Complete CDR graft (0 mutations)
    - 0b: Hallmarks + IMGT2 (5 mutations)
    - 0c: + ALL verniers (~20 mutations)
  Track 1-4: Similar to Original but CDR-conditional focused

Usage:
  python vhh_designer_v8_9.py -i input.fasta \
      --rules analysis_rules_v7_all_positions.json \
      --archetypes analysis_vernier_archetypes_v7.json \
      --hallmark-db comprehensive_subfamily_analysis_imgt.xlsx \
      --target-hallmarks AUTO \
      --n-generate 200000 \
      --n-select 198

Author: Claude (Anthropic)
Date: January 2026
```

---

### `vhh_designer_v9_1.py`
**Version:** 9.1
**Path:** `active/analysis/vhh_designer_v9_1.py`
**Stats:** 7664 lines | 317.9KB | 103 functions | 12 classes
**Modified:** 2026-02-01 23:08

**Description:**
```
VHH Designer v9.1 - Corrected Scoring & Multi-Hallmark Pool
============================================================

v9.1 MAJOR CHANGES:

1. **CRITICAL: SCORING DIRECTION FIX**
   - FIXED: Combined_score now correctly ranks better sequences higher
   - Was inverted: higher rule compliance → WORSE rank (WRONG!)
   - Now correct: higher rule compliance → BETTER rank
   - ESM2 loss: lower loss (more natural) → better rank
   - pLDDT: higher pLDDT (better structure) → better rank  
   - Rules: higher compliance → better rank
   - This was a fundamental bug - previous versions ranked the WORST candidates first!

2. DEFAULT HALLMARK POOL:
   - New DEFAULT_HALLMARK_POOL with 7 curated hallmarks
   - FERG, FERF, YQRL, FERA, YERL, FKRG, FQRA
   - --target-hallmarks now defaults to 'DEFAULT' (was 'FERG')
   - FQRA added for 100% vernier match to typical VH inputs

3. PERFORMANCE FIX - O(n²) BUG:
   - Fixed to_dataframe() hanging on large candidate pools (85K+)
   - Was: O(n) scan per candidate = O(n²) total = 7.2B comparisons
   - Now: O(1) dict lookup per candidate = O(n) total

4. _all.csv FIX:
   - Now saves ALL scored candidates before selection
   - Copies candidates list immediately after scoring

5. FILE NAMING:
   - Includes n-generate and temperature: v9_0_n100000_t08_M69
   - Format: {datetime}_{version}_n{n_generate}_t{temp}_{seq_name}

6. RANK ORDERING FIX:
   - Fixed: rank now correctly reflects score order within selected candidates
   - Was: candidates added during repair/fill phases were appended at end, breaking score order
   - Now: final sort ensures rank 1 = best score among selected, rank 2 = second best, etc.

v8.9.3 CHANGES:
- Cross-track sequence deduplication (keeps higher track version)
- MSA-compatible ID format (M69_Orig_FERG_T4R_HVn_r001_s152)
- Hallmark-aware track minimums
- Hallmark cap enforcement (35% max per hallmark)
- 6-phase selection pipeline

v8.9.1 CRITICAL FIXES:
1. CYSTEINE PROTECTION (v8.9.1):
   - NEVER mutate FROM any cysteine position (disulfide bonds are essential)
   - Added CONSERVED_CYSTEINE_POSITIONS = {22, 23, 103, 104}
   - Three layers of protection:
     a) Position exclusion: Skip CONSERVED_CYSTEINE_POSITIONS in FRAMEWORK_CONSENSUS_TIERS loops
     b) Mutation filter: if scaffold_aa == 'C': continue (in all mutation generation)
     c) Builder safety: Filter in _build_*_candidate() catches any missed mutations
   - Added validate_canonical_disulfide() to reject candidates missing disulfide bond
   - Fixed bug where FRAMEWORK_CONSENSUS_TIERS['A'][22]='S' and ['B'][103]='Y' 
     would mutate Universal scaffold cysteines (C22→S, C103→Y)
   - NOTE: Cysteine positions NOT in EXCLUDE_FROM_FRAMEWORK because YQRL needs to
     mutate TO cysteine at 23/104 - protection is via individual checks instead

2. YQRL RIGID HALLMARK FIX (v8.9.1):
   - Added is_rigid check to _generate_original (was only in _generate_universal)
   - YQRL now correctly gets only Track 0 controls on Original scaffold

3. REDUCED YQRL FRAMEWORK STEPS (v8.9.1):
   - Reduced from 27 individual steps to 13 (5 pairs + 8 groups)
   - Maintains interpretability while reducing control explosion

4. HALLMARK DIVERSITY ENFORCEMENT (v8.9.1):
   - Added diversity enforcement in select_with_track_quotas
   - Prevents single hallmark (e.g., FERG) from dominating ranked candidates

MAJOR UPDATES from v8.8.3:

1. SIMPLIFIED YQRL CONTROL STRUCTURE (v8.9):
   - 0a: GRAFT (0 mutations) - pure CDR graft
   - 0b: Hallmarks + IMGT2 (5 mutations) - ALWAYS include IMGT2 with hallmarks
   - 0c: Hallmarks + IMGT2 + ALL verniers (~20 mutations)
   - 0d+: Progressive framework in PAIRS (2 at a time), all positions >80% consensus
   - Removed separate "hallmarks only (4 mut)" - always include IMGT2

2. ALL HALLMARKS GET GRAFT FIRST (v8.9):
   - 0a: GRAFT (0 mutations) - always first for ALL hallmarks
   - 0b: Hallmarks + IMGT2 (5 mutations)
   - 0c: + ALL verniers (only if hallmark avg vernier consensus ≥60%)

3. TRACK 1 EXPANDED (v8.9):
   - Now includes ALL verniers, not just "load-bearing"
   - Expanded from max 12 to max 20 probes per hallmark
   - Tests each vernier position individually

4. TRACK 2 REDEFINED: "MINIMAL STRUCTURED PERTURBATIONS" (v8.9):
   - Base: Hallmarks + IMGT2 (IMGT2 now protected like hallmarks)
   - Add: Exactly ONE 2-3 position motif (weighted by family support)
   - Add: 0-K extra vernier mutations (K_EXTRA_VERNIERS_T2 = 2, capped)
   - NO framework tier sprinkling - keeps it interpretable

5. TRACK 3 REDEFINED: "MOTIF ATTRIBUTION" (v8.9):
   - Base: Hallmarks + IMGT2
   - Add: Exactly ONE 3-5 position motif (weighted)
   - Add: 0-K extra vernier mutations (K_EXTRA_VERNIERS_T3 = 2, max 3)
   - NO framework tier sprinkling
   - Goal: Interpretable - can analyze "motif X worked/didn't work"

6. TRACK 4 REDEFINED: "PRODUCTION OPTIMIZER" with CAPPED LAYERS (v8.9):
   - Full probabilistic layering for production candidates
   - Vernier layer: Bernoulli + cap (K_VERNIER_MAX = 8)
   - Framework tiers applied in order with individual caps:
     - Tier S (≥95%): K_S_MAX = 10
     - Tier A (90-95%): K_A_MAX = 6
     - Tier B (85-90%): K_B_MAX = 3
     - Tier C (75-85%): K_C_MAX = 2 (only when desperate)
   - CDR-conditional rules: capped at K_CDR_RULES_MAX = 2
   - Adaptive expansion when uniqueness drops

7. UNIVERSAL SCAFFOLD CONTROLS (v8.9):
   - 0a: Complete CDR graft into Universal (0 mutations)
   - 0b: Hallmarks + IMGT2 (5 mutations)
   - 0c: Hallmarks + IMGT2 + ALL verniers

8. IMGT2 NOW PROTECTED LIKE HALLMARKS (v8.9):
   - PROTECTED_POSITIONS = {2, 42, 49, 50, 52}
   - Motifs cannot override these positions

KEY DISTINCTION:
- Track 3 = "Motif attribution": one motif + minimal extras → interpretable
- Track 4 = "Production optimizer": full probabilistic layering → coverage

TRACK STRUCTURE (v8.9):

  ORIGINAL FRAMEWORK:
  Track 0: Controls (exempt)
    - 0a: GRAFT (0 mutations) - ALL hallmarks get this
    - 0b: Hallmarks + IMGT2 (5 mutations)
    - 0c: + ALL verniers (~20 mutations) - if avg vernier ≥60%
    - 0d+: YQRL only: + Progressive FW pairs (>80% consensus)
  Track 1: Single vernier probes - ALL verniers (RANKED)
  Track 2: Minimal structured perturbations - 1 motif + 0-2 verniers (RANKED)
  Track 3: Motif attribution - 1 motif (3-5 pos) + 0-3 verniers (RANKED)
  Track 4: Production optimizer - capped probabilistic layers (RANKED)

  UNIVERSAL FRAMEWORK:
  Track 0: Controls (exempt)
    - 0a: Complete CDR graft (0 mutations)
    - 0b: Hallmarks + IMGT2 (5 mutations)
    - 0c: + ALL verniers (~20 mutations)
  Track 1-4: Similar to Original but CDR-conditional focused

Usage:
  python vhh_designer_v8_9.py -i input.fasta \
      --rules analysis_rules_v7_all_positions.json \
      --archetypes analysis_vernier_archetypes_v7.json \
      --hallmark-db comprehensive_subfamily_analysis_imgt.xlsx \
      --target-hallmarks AUTO \
      --n-generate 200000 \
      --n-select 198

Author: Claude (Anthropic)
Date: January 2026
```

---

### `vhh_epistasis_v3_imgt.py`
**Version:** 3
**Path:** `active/analysis/vhh_epistasis_v3_imgt.py`
**Stats:** 1608 lines | 61.3KB | 37 functions | 7 classes
**Modified:** 2026-01-02 23:07

**Description:**
```
VHH Full Epistasis Analysis - OVERNIGHT VERSION v3 (IMGT NUMBERING)
====================================================================

VERSION 3 CHANGES:
1. PROPER IMGT NUMBERING using AntPack/ANARCI
   - Positions stored as IMGT numbers (e.g., IMGT42, IMGT49)
   - No longer depends on FR2 starting at exact position
   - Handles insertions/deletions correctly

2. DUAL EXTRACTION MODE:
   - --use-imgt: Use AntPack for accurate IMGT numbering (slower)
   - Default: Use regex-based extraction (faster, backward compatible)

3. IMGT POSITION NAMES:
   - IMGT42 = Hallmark 1 (F/Y)
   - IMGT49 = Hallmark 2 (E/G)
   - IMGT50 = Hallmark 3 (R/L)
   - IMGT52 = Hallmark 4 (G/L/W)
   - Plus Vernier zone positions

RETAINED FROM v2:
- All 5 analysis modules (Compensation, Vernier, Conditional, Rules, MI)
- Streaming stats with Welford's algorithm
- Memory-efficient integer encoding for MI
- CLI flags for light vs full mode
- Progress bars with tqdm

Position Mapping Reference (regex mode):
  FR2: W  F  R  Q  x  P  x  x  E  R  x  G  L  x
  Idx: 0  1  2  3  4  5  6  7  8  9  10 11 12 13
  IMGT:41 42 43 44 45 46 47 48 49 50 51 52 53 54
  
  Old name -> IMGT name:
  FR2_2  -> IMGT42 (hallmark 1)
  FR2_9  -> IMGT49 (hallmark 2)
  FR2_10 -> IMGT50 (hallmark 3)
  FR2_12 -> IMGT52 (hallmark 4)
```

---

### `vhh_naturalness_analyzer_v4.py`
**Version:** 4
**Path:** `active/analysis/vhh_naturalness_analyzer_v4.py`
**Stats:** 995 lines | 36.2KB | 17 functions | 5 classes
**Modified:** 2026-01-02 02:08

**Description:**
```
VHH CDR-Framework Naturalness Analyzer v4
=========================================

Uses camelid epistasis data as a PRIOR on CDR3-framework compatibility,
NOT as a binary fold/don't fold gate.

Key insight: For humanized scaffolds, foldability is known from experiments.
The question is: "Given this framework, how natural is this CDR3?"

Two modes:
1. Full sequence input: Extracts CDRs and frameworks automatically
2. Framework + CDR input: Uses provided framework regions with grafted CDRs

Output:
- Original sequence (if full sequence mode)
- Naturalness score (lower = more natural)
- Per-feature z-scores (length, charge, aromatics)
- Framework annotation (classical VHH, humanized, VH-like)
- Risk assessment based on CDR-framework compatibility
- Suggested optimizations

Usage:
  # Full sequences
  python vhh_naturalness_analyzer_v4.py --input sequences.csv --epistasis epistasis_v2_full.pkl
  
  # Framework library + CDRs
  python vhh_naturalness_analyzer_v4.py --input designs.csv --frameworks VHH_frameworks.xlsx --epistasis epistasis_v2_full.pkl

Changelog:
  v4 (2026-01-02):
    - FIXED: CDR3 extraction no longer includes conserved Cys (position 104)
      CDR3 now correctly starts at position 105 per IMGT
    - FIXED: FR3 extraction now correctly INCLUDES the conserved Cys at end
    - ADDED: original_sequence column in output (right after id)
    - FIXED: Removed duplicated fr4_patterns definitions
  v3 (2025-12-13):
    - Fixed CDR3 extraction edge cases
    - 100% extraction success on test set
  v2 (2025-12-12):
    - CDR3 extraction bug fixes
  v1 (2025-12-10):
    - Initial naturalness scoring implementation
```

---

## 📁 `active/analysis/annotate_all_shards/`

### `annotate_all_shards.py`
**Path:** `active/analysis/annotate_all_shards/annotate_all_shards.py`
**Stats:** 829 lines | 26.9KB | 14 functions | 0 classes
**Modified:** 2026-01-04 14:47

**Description:**
```
VHH Database ANARCI Annotation Pipeline
========================================

Processes all VHH database shards with ANARCI to create a unified, fully-annotated database.

Features:
- Runs ANARCI for IMGT numbering (if needed)
- Extracts CDR1, CDR2, CDR3 and FR1, FR2, FR3, FR4
- Classifies VHH family (F_C2, Y_C2, VH_like, etc.)
- Preserves metadata (targets, source, patent info, etc.)
- Outputs both CSV and NPZ formats

Estimated time:
- If shards have CDRs already: ~30-60 min for 12M sequences (no ANARCI needed)
- If ANARCI needed: ~8-12 hours for 12M sequences (batch processing)

Usage:
    python annotate_all_shards.py         --input-dir data/databases/shards/         --output-prefix data/databases/annotated/vhh_annotated_full         --batch-size 2000

Author: Claude (Anthropic)
Date: January 2026
```

---

### `annotate_all_shards_v2.py`
**Path:** `active/analysis/annotate_all_shards/annotate_all_shards_v2.py`
**Stats:** 869 lines | 29.3KB | 14 functions | 0 classes
**Modified:** 2026-01-06 14:46

**Description:**
```
VHH Database ANARCI Annotation Pipeline
========================================

Processes all VHH database shards with ANARCI to create a unified, fully-annotated database.

Features:
- Runs ANARCI for IMGT numbering (if needed)
- Extracts CDR1, CDR2, CDR3 and FR1, FR2, FR3, FR4
- Classifies VHH family (F_C2, Y_C2, VH_like, etc.)
- Preserves metadata (targets, source, patent info, etc.)
- Outputs both CSV and NPZ formats

Estimated time:
- If shards have CDRs already: ~30-60 min for 12M sequences (no ANARCI needed)
- If ANARCI needed: ~8-12 hours for 12M sequences (batch processing)

Usage:
    python annotate_all_shards.py         --input-dir data/databases/shards/         --output-prefix data/databases/annotated/vhh_annotated_full         --batch-size 2000

Author: Claude (Anthropic)
Date: January 2026
```

---

### `annotate_all_shards_v4.py`
**Version:** 3
**Path:** `active/analysis/annotate_all_shards/annotate_all_shards_v4.py`
**Stats:** 510 lines | 17.4KB | 12 functions | 0 classes
**Modified:** 2026-01-06 16:01

**Description:**
```
VHH Database ANARCI Annotation Pipeline (streaming, low-memory) - v3
===================================================================

Key differences vs v2:
- Writes ONE CSV per shard incrementally (no giant in-memory concat)
- Processes each NPZ shard in chunks, and ANARCI in batches inside each chunk
- Emits IMGT per-position columns for FR3/FR4 by default (optional FR2)
- Produces a summary JSON + manifest of outputs

Example:
  python annotate_all_shards_v3.py     -i data/databases/shards     -o data/databases/annotated/vhh_full_annotated_final_v2     --use-anarci     -b 2000     --chunk-size 50000
```

---

### `annotate_all_shards_v5.py`
**Version:** 3
**Path:** `active/analysis/annotate_all_shards/annotate_all_shards_v5.py`
**Stats:** 586 lines | 19.7KB | 14 functions | 0 classes
**Modified:** 2026-01-06 16:56

**Description:**
```
VHH Database ANARCI Annotation Pipeline (streaming, low-memory) - v3
===================================================================

Key differences vs v2:
- Writes ONE CSV per shard incrementally (no giant in-memory concat)
- Processes each NPZ shard in chunks, and ANARCI in batches inside each chunk
- Emits IMGT per-position columns for FR3/FR4 by default (optional FR2)
- Produces a summary JSON + manifest of outputs

Example:
  python annotate_all_shards_v3.py     -i data/databases/shards     -o data/databases/annotated/vhh_full_annotated_final_v2     --use-anarci     -b 2000     --chunk-size 50000
```

---

### `annotate_all_shards_v6.py`
**Version:** 3
**Path:** `active/analysis/annotate_all_shards/annotate_all_shards_v6.py`
**Stats:** 618 lines | 20.9KB | 14 functions | 0 classes
**Modified:** 2026-01-06 21:44

**Description:**
```
VHH Database ANARCI Annotation Pipeline (streaming, low-memory) - v3
===================================================================

Key differences vs v2:
- Writes ONE CSV per shard incrementally (no giant in-memory concat)
- Processes each NPZ shard in chunks, and ANARCI in batches inside each chunk
- Emits IMGT per-position columns for FR3/FR4 by default (optional FR2)
- Produces a summary JSON + manifest of outputs

Example:
  python annotate_all_shards_v3.py     -i data/databases/shards     -o data/databases/annotated/vhh_full_annotated_final_v2     --use-anarci     -b 2000     --chunk-size 50000
```

---

### `annotate_all_shards_v7.py`
**Version:** 3
**Path:** `active/analysis/annotate_all_shards/annotate_all_shards_v7.py`
**Stats:** 643 lines | 22.1KB | 14 functions | 0 classes
**Modified:** 2026-01-06 22:33

**Description:**
```
VHH Database ANARCI Annotation Pipeline (streaming, low-memory) - v3
===================================================================

Key differences vs v2:
- Writes ONE CSV per shard incrementally (no giant in-memory concat)
- Processes each NPZ shard in chunks, and ANARCI in batches inside each chunk
- Emits IMGT per-position columns for FR3/FR4 by default (optional FR2)
- Produces a summary JSON + manifest of outputs

Example:
  python annotate_all_shards_v3.py     -i data/databases/shards     -o data/databases/annotated/vhh_full_annotated_final_v2     --use-anarci     -b 2000     --chunk-size 50000
```

---

## 📁 `active/analysis/vhh_compensation/`

### `vhh_compensation_imgt.py`
**Version:** 7.
**Path:** `active/analysis/vhh_compensation/vhh_compensation_imgt.py`
**Stats:** 476 lines | 16.0KB | 10 functions | 1 classes
**Modified:** 2026-01-05 20:10

**Description:**
```
VHH Compensation Analysis - IMGT Positions
===========================================

Analyzes CDR→FR correlations using true IMGT positions from the annotated database.
Outputs rules in JSON format compatible with vhh_designer_v7.

Input: vhh_full_annotated_v3.csv (12M sequences with FR/CDR columns)
Output: compensation_imgt_rules.json

Key difference from old compensation analysis:
- OLD: FR2[3], FR3[6] (substring positions, dataset-specific)
- NEW: IMGT42, IMGT71 (universal IMGT positions)

Usage:
    python vhh_compensation_imgt.py \
        --csv data/databases/annotated/vhh_full_annotated_v3.csv \
        --output models/compensation/imgt_v1 \
        --chunk-size 500000

Author: Claude (Anthropic)
Date: January 2026
```

---

### `vhh_compensation_imgt_v2.py`
**Version:** 2
**Path:** `active/analysis/vhh_compensation/vhh_compensation_imgt_v2.py`
**Stats:** 391 lines | 13.9KB | 11 functions | 1 classes
**Modified:** 2026-01-05 20:29

**Description:**
```
VHH Compensation Analysis - IMGT Positions (v2)

What this script does
- Streams an annotated VHH CSV (fr1/cdr1/fr2/cdr2/fr3/cdr3/fr4 + family + hallmarks)
- Learns simple conditional residue rules: condition -> (IMGT position -> preferred AA)
- Outputs JSON rules compatible with the vhh_designer_v7 IMGT-rule schema:
    {condition, position: "IMGT##", suggested_aa, confidence, support, lift, baseline_prob, source}

Why v2
- Adds sanity checks to ensure the FR strings in the CSV are aligned to true IMGT boundaries
- Greatly reduces runtime/memory by:
    * restricting analysis to a configurable set of FR positions
    * using a smaller, higher-signal condition set by default
    * avoiding pandas iterrows()

USAGE (recommended starting point)
python vhh_compensation_imgt_v2.py   --csv data/databases/annotated/vhh_full_annotated_v3.csv   --output models/compensation/imgt_v1   --chunk-size 200000   --positions fr3fr4   --conditions minimal   --min-support 5000   --min-confidence 0.70   --min-lift 1.15
```

---

### `vhh_compensation_imgt_v4.py`
**Version:** 2
**Path:** `active/analysis/vhh_compensation/vhh_compensation_imgt_v4.py`
**Stats:** 396 lines | 14.5KB | 12 functions | 1 classes
**Modified:** 2026-01-06 13:21

**Description:**
```
VHH Compensation Analysis - IMGT Positions (v2)

What this script does
- Streams an annotated VHH CSV (fr1/cdr1/fr2/cdr2/fr3/cdr3/fr4 + family + hallmarks)
- Learns simple conditional residue rules: condition -> (IMGT position -> preferred AA)
- Outputs JSON rules compatible with the vhh_designer_v7 IMGT-rule schema:
    {condition, position: "IMGT##", suggested_aa, confidence, support, lift, baseline_prob, source}

Why v2
- Adds sanity checks to ensure the FR strings in the CSV are aligned to true IMGT boundaries
- Greatly reduces runtime/memory by:
    * restricting analysis to a configurable set of FR positions
    * using a smaller, higher-signal condition set by default
    * avoiding pandas iterrows()

USAGE (recommended starting point)
python vhh_compensation_imgt_v2.py   --csv data/databases/annotated/vhh_full_annotated_v3.csv   --output models/compensation/imgt_v1   --chunk-size 200000   --positions fr3fr4   --conditions minimal   --min-support 5000   --min-confidence 0.70   --min-lift 1.15
```

---

## 📁 `active/analysis/vhh_designer/`

### `vhh_designer_v2.py`
**Version:** 2
**Path:** `active/analysis/vhh_designer/vhh_designer_v2.py`
**Stats:** 540 lines | 25.9KB | 21 functions | 6 classes
**Modified:** 2026-01-02 21:09

**Description:**
```
VHH Designer v2 - Mouse HC to VHH with Direct Mutations
========================================================

Converts a mouse heavy chain to VHH candidates using TWO approaches:
  1. UNIVERSAL: Graft CDRs into Universal VHH scaffold + apply mutations
  2. ORIGINAL: Keep original mouse FRs + apply VHH hallmark mutations

Data sources:
- correlation_results_v3_compensation.pkl: CDR->FR mutation rules (254 rules)
- epistasis_v2_full.pkl: Vernier cluster statistics for naturalness scoring

Usage:
  # Both approaches (default - 46 Universal + 46 Original)
  python vhh_designer_v2.py -s "EVQLVESGG..." -c correlations.pkl -e epistasis.pkl

  # Universal only
  python vhh_designer_v2.py -s "EVQLVESGG..." -c correlations.pkl -e epistasis.pkl --mode universal

  # Original only  
  python vhh_designer_v2.py -s "EVQLVESGG..." -c correlations.pkl -e epistasis.pkl --mode original
```

---

### `vhh_designer_v3.py`
**Version:** 3
**Path:** `active/analysis/vhh_designer/vhh_designer_v3.py`
**Stats:** 1160 lines | 44.6KB | 33 functions | 8 classes
**Modified:** 2026-01-03 21:53

**Description:**
```
VHH Designer v3 - Three-Pass Rule Application
==============================================

Converts a mouse heavy chain to VHH candidates using proper rule ordering:

  PASS 1: COMPENSATION rules (CDR features → FR residues)
          Source: correlation_results_v3_compensation.pkl
          
  PASS 2: EPISTASIS rules (FR-FR interactions given CDR context)
          Source: epistasis_v2_full.pkl → analysis_4_higher_order_rules
          
  PASS 3: VERNIER archetype (family-specific FR pattern)
          Source: vernier_archetypes from either file

Family Classification (based on IMGT positions):
  - F_C2: pos42=F, pos49=E, pos50=R, 2 cysteines
  - Y_C2: pos42=Y, pos49=E, pos50=R, 2 cysteines  
  - F_C4: pos42=F, pos49=E, pos50=R, 4 cysteines
  - Y_C4: pos42=Y, pos49=E, pos50=R, 4 cysteines
  - VH_like: pos50=L (human VH-like)
  - Non_classical: other combinations

Usage:
  python vhh_designer_v3.py -i M69.fasta \
      --compensation correlation_results_v3_compensation.pkl \
      --epistasis epistasis_v2_full.pkl \
      --mode original
```

---

### `vhh_designer_v4.py`
**Version:** 3
**Path:** `active/analysis/vhh_designer/vhh_designer_v4.py`
**Stats:** 1318 lines | 51.5KB | 33 functions | 8 classes
**Modified:** 2026-01-03 22:09

**Description:**
```
VHH Designer v3 - Three-Pass Rule Application
==============================================

Converts a mouse heavy chain to VHH candidates using proper rule ordering:

  PASS 1: COMPENSATION rules (CDR features → FR residues)
          Source: correlation_results_v3_compensation.pkl
          
  PASS 2: EPISTASIS rules (FR-FR interactions given CDR context)
          Source: epistasis_v2_full.pkl → analysis_4_higher_order_rules
          
  PASS 3: VERNIER archetype (family-specific FR pattern)
          Source: vernier_archetypes from either file

Family Classification (based on IMGT positions):
  - F_C2: pos42=F, pos49=E, pos50=R, 2 cysteines
  - Y_C2: pos42=Y, pos49=E, pos50=R, 2 cysteines  
  - F_C4: pos42=F, pos49=E, pos50=R, 4 cysteines
  - Y_C4: pos42=Y, pos49=E, pos50=R, 4 cysteines
  - VH_like: pos50=L (human VH-like)
  - Non_classical: other combinations

Usage:
  python vhh_designer_v3.py -i M69.fasta \
      --compensation correlation_results_v3_compensation.pkl \
      --epistasis epistasis_v2_full.pkl \
      --mode original
```

---

### `vhh_designer_v5.py`
**Version:** 5
**Path:** `active/analysis/vhh_designer/vhh_designer_v5.py`
**Stats:** 1559 lines | 61.5KB | 36 functions | 9 classes
**Modified:** 2026-01-04 13:17

**Description:**
```
VHH Designer v5 - Protected Framework Rule Application
==============================================

Converts a mouse heavy chain to VHH candidates using proper rule ordering:

  For UNIVERSAL mode (grafting CDRs onto humanized scaffold):
    1. Start with Universal scaffold (FR2 = WFRQAPGQGLEAVA, hallmarks = FGLA)
    2. Apply safe rules (FR3, FR4 positions) - FR2 core is PROTECTED
    3. Hallmark positions (IMGT 42,49,50,52) only changed via explicit patterns
    
  For ORIGINAL mode (keeping input FRs with VHH-izing mutations):
    PASS 1: HALLMARK mutations (IMGT 42,49,50,52 → VHH patterns like FERG, YERL)
    PASS 2: VERNIER archetype (family-specific FR pattern for stability)
    PASS 3: COMPENSATION/EPISTASIS rules (CDR-conditioned fine-tuning)

Family Classification (based on IMGT positions):
  - F_C2: pos42=F, pos49=E, pos50=R, 2 cysteines
  - Y_C2: pos42=Y, pos49=E, pos50=R, 2 cysteines  
  - F_C4: pos42=F, pos49=E, pos50=R, 4 cysteines
  - Y_C4: pos42=Y, pos49=E, pos50=R, 4 cysteines
  - VH_like: pos50=L (human VH-like)
  - Non_classical: other combinations

CRITICAL: Substring indices (FR2_4, FR3_6) are 1-based and NOT equivalent to IMGT positions.
This script maps them correctly and protects FR2 core in universal mode.

Usage:
  python vhh_designer_v5.py -i M69.fasta \
      --compensation correlation_results_v3_compensation.pkl \
      --epistasis epistasis_v2_full.pkl \
      --mode original
```

---

### `vhh_designer_v5_1.py`
**Version:** 5
**Path:** `active/analysis/vhh_designer/vhh_designer_v5_1.py`
**Stats:** 1561 lines | 61.9KB | 36 functions | 9 classes
**Modified:** 2026-01-04 14:25

**Description:**
```
VHH Designer v5 - Protected Framework Rule Application
==============================================

Converts a mouse heavy chain to VHH candidates using proper rule ordering:

  For UNIVERSAL mode (grafting CDRs onto humanized scaffold):
    1. Start with Universal scaffold (FR2 = WFRQAPGQGLEAVA, hallmarks = FGLA)
    2. Apply safe rules (FR3, FR4 positions) - FR2 core is PROTECTED
    3. Hallmark positions (IMGT 42,49,50,52) only changed via explicit patterns
    
  For ORIGINAL mode (keeping input FRs with VHH-izing mutations):
    PASS 1: HALLMARK mutations (IMGT 42,49,50,52 → VHH patterns like FERG, YERL)
    PASS 2: VERNIER archetype (family-specific FR pattern for stability)
    PASS 3: COMPENSATION/EPISTASIS rules (CDR-conditioned fine-tuning)

Family Classification (based on IMGT positions):
  - F_C2: pos42=F, pos49=E, pos50=R, 2 cysteines
  - Y_C2: pos42=Y, pos49=E, pos50=R, 2 cysteines  
  - F_C4: pos42=F, pos49=E, pos50=R, 4 cysteines
  - Y_C4: pos42=Y, pos49=E, pos50=R, 4 cysteines
  - VH_like: pos50=L (human VH-like)
  - Non_classical: other combinations

CRITICAL: Substring indices (FR2_4, FR3_6) are 1-based and NOT equivalent to IMGT positions.
This script maps them correctly and protects FR2 core in universal mode.

Usage:
  python vhh_designer_v5.py -i M69.fasta \
      --compensation correlation_results_v3_compensation.pkl \
      --epistasis epistasis_v2_full.pkl \
      --mode original
```

---

### `vhh_designer_v7.py`
**Version:** 7
**Path:** `active/analysis/vhh_designer/vhh_designer_v7.py`
**Stats:** 1721 lines | 68.5KB | 42 functions | 10 classes
**Modified:** 2026-01-05 13:30

**Description:**
```
VHH Designer v7 - IMGT-based Rules with Full Integration
=========================================================

Combines:
- ChatGPT's FR2 IMGT mapping fix (W at IMGT41, hallmarks at 42,49,50,52)
- New IMGT rules from 12M sequence epistasis analysis (JSON format)
- v5's NaturalnessScorer, VernierEngine, diverse candidate strategies
- Compensation rules support

For UNIVERSAL mode (grafting CDRs onto humanized scaffold):
  1. Start with Universal scaffold (FR2 = WFRQAPGQGLEAVA, hallmarks = FGLA)
  2. Apply IMGT FR3 rules (from JSON) - FR2 core is PROTECTED
  3. Apply safe compensation rules to FR3/FR4
    
For ORIGINAL mode (keeping input FRs with VHH-izing mutations):
  PASS 1: HALLMARK mutations (IMGT 42,49,50,52 → target pattern like FERG, YERL)
  PASS 2: IMGT FR3 rules (family-gated by hallmark pattern)
  PASS 3: COMPENSATION rules + VERNIER archetype

Family Classification (based on IMGT positions):
  - F_C2: pos42=F, pos49=E, pos50=R, 2 cysteines
  - Y_C2: pos42=Y, pos49=E, pos50=R, 2 cysteines  
  - F_C4: pos42=F, pos49=E, pos50=R, 4 cysteines
  - Y_C4: pos42=Y, pos49=E, pos50=R, 4 cysteines
  - VH_like: pos50=L (human VH-like)
  - Non_classical: other combinations

CRITICAL FR2 IMGT MAPPING (corrected):
  FR2: W-F-R-Q-A-P-G-Q-G-L-E-A-V-A
  IMGT: 41-42-43-44-45-46-47-48-49-50-51-52-53-54
  So: IMGT42=FR2[1], IMGT49=FR2[8], IMGT50=FR2[9], IMGT52=FR2[11] (0-based)

Usage:
  python vhh_designer_v7.py -i M69.fasta \
      --imgt-rules models/epistasis/imgt_v1/imgt_rules.json \
      --compensation models/correlations/correlation_results_v3_compensation.pkl \
      --epistasis models/epistasis/current/epistasis_v2_full.pkl \
      --target-hallmarks FERG \
      --mode both

Author: Claude (Anthropic) + ChatGPT fixes
Date: January 2026
```

---

## 📁 `active/analysis/vhh_epistasis/`

### `vhh_epistasis_overnight_final.py`
**Version:** 2
**Path:** `active/analysis/vhh_epistasis/vhh_epistasis_overnight_final.py`
**Stats:** 1193 lines | 45.2KB | 31 functions | 7 classes
**Modified:** 2025-12-06 23:02

**Description:**
```
VHH Full Epistasis Analysis - OVERNIGHT VERSION v2
===================================================

Fixed based on ChatGPT's feedback:
1. Fixed cond_results NameError bug
2. Streaming stats for FeatureCompensationScanner (Welford's algorithm)
3. Integer encoding for MI analyzer (memory efficient)
4. Subsample cap for MI (default 2M sequences)
5. CLI flags for light vs full mode
6. Removed dead code in MI analyzer
7. Memory-efficient HigherOrderRuleFinder (streaming counts)

Implements ALL of ChatGPT's recommendations at maximum depth.
```

---

### `vhh_epistasis_overnight_v2.py`
**Version:** 2
**Path:** `active/analysis/vhh_epistasis/vhh_epistasis_overnight_v2.py`
**Stats:** 1262 lines | 47.9KB | 31 functions | 7 classes
**Modified:** 2025-12-07 23:34

**Description:**
```
VHH Full Epistasis Analysis - OVERNIGHT VERSION v2 (CORRECTED)
===============================================================

CORRECTIONS FROM v1:
1. FIXED POSITION MAPPING:
   - FR2_2 (index 1) = IMGT 42 = Y/F hallmark (was incorrectly FR2_4)
   - FR2_9 (index 8) = IMGT 49 = E/Q hallmark
   - FR2_10 (index 9) = IMGT 50 = R/L hallmark (was incorrectly FR2_12)
   - FR2_12 (index 11) = IMGT 52 = G/W hallmark (was incorrectly FR2_14)

2. CORRECTED FAMILY CLASSIFICATION:
   - Classical VHH: pos42=F/Y AND pos49=E/Q AND pos50=R
   - Y_C2, F_C2, F_C4: Based on pos42 + total cysteine count
   - VH_like: pos50=L (includes humanized VHHs)
   - VHH_W52: has W at pos52
   - Non_classical: everything else

3. Uses TOTAL cysteine count (not just CDR3) for C2/C4 classification

RETAINED FROM v1:
- Fixed cond_results NameError bug
- Streaming stats for FeatureCompensationScanner (Welford's algorithm)
- Integer encoding for MI analyzer (memory efficient)
- Subsample cap for MI (default 2M sequences)
- CLI flags for light vs full mode
- Progress bars with tqdm

Position Mapping Reference:
  FR2: W  F  R  Q  x  P  x  x  E  R  x  G  L  x
  Idx: 0  1  2  3  4  5  6  7  8  9  10 11 12 13
  IMGT:41 42 43 44 45 46 47 48 49 50 51 52 53 54
```

---

### `visualize_epistasis.py`
**Path:** `active/analysis/vhh_epistasis/visualize_epistasis.py`
**Stats:** 345 lines | 11.9KB | 0 functions | 0 classes
**Modified:** 2025-12-07 13:16

**Description:**
```
Visualize VHH Epistasis Analysis Results
```

---

## 📁 `active/analysis/vhh_epistasis_to_imgt/`

### `vhh_epistasis_imgt_ligh_v3t.py`
**Version:** 1
**Path:** `active/analysis/vhh_epistasis_to_imgt/vhh_epistasis_imgt_ligh_v3t.py`
**Stats:** 867 lines | 30.1KB | 21 functions | 4 classes
**Modified:** 2026-01-05 00:03

**Description:**
```
VHH IMGT Epistasis Analysis - LIGHTWEIGHT VERSION
==================================================

Skips heavy computations to run faster and use less memory:
- SKIPS: Random Forest models (analysis_3) - saves ~10GB RAM, hours of compute
- SKIPS: Mutual Information matrices (analysis_5) - saves ~5GB RAM
- KEEPS: Compensation stats (analysis_1) - lightweight streaming
- KEEPS: Vernier clusters (analysis_2) - essential for designer
- KEEPS: Higher-order rules (analysis_4) - essential for designer

Expected runtime: ~2-4 hours for 12M sequences (vs 2-5 days for full)
Expected memory: ~4GB peak (vs 20GB+ for full)

Usage:
    python vhh_epistasis_imgt_light.py         --npz data/processed/full_12m/full_12m_imgt.npz         --output models/epistasis/imgt_light_v1         --checkpoint-interval 1000000
```

---

### `vhh_epistasis_imgt_light.py`
**Version:** 1
**Path:** `active/analysis/vhh_epistasis_to_imgt/vhh_epistasis_imgt_light.py`
**Stats:** 789 lines | 26.9KB | 21 functions | 4 classes
**Modified:** 2026-01-04 14:19

**Description:**
```
VHH IMGT Epistasis Analysis - LIGHTWEIGHT VERSION
==================================================

Skips heavy computations to run faster and use less memory:
- SKIPS: Random Forest models (analysis_3) - saves ~10GB RAM, hours of compute
- SKIPS: Mutual Information matrices (analysis_5) - saves ~5GB RAM
- KEEPS: Compensation stats (analysis_1) - lightweight streaming
- KEEPS: Vernier clusters (analysis_2) - essential for designer
- KEEPS: Higher-order rules (analysis_4) - essential for designer

Expected runtime: ~2-4 hours for 12M sequences (vs 2-5 days for full)
Expected memory: ~4GB peak (vs 20GB+ for full)

Usage:
    python vhh_epistasis_imgt_light.py         --npz data/processed/full_12m/full_12m_imgt.npz         --output models/epistasis/imgt_light_v1         --checkpoint-interval 1000000
```

---

### `vhh_epistasis_imgt_light_v2.py`
**Version:** 1
**Path:** `active/analysis/vhh_epistasis_to_imgt/vhh_epistasis_imgt_light_v2.py`
**Stats:** 867 lines | 30.1KB | 21 functions | 4 classes
**Modified:** 2026-01-04 23:17

**Description:**
```
VHH IMGT Epistasis Analysis - LIGHTWEIGHT VERSION
==================================================

Skips heavy computations to run faster and use less memory:
- SKIPS: Random Forest models (analysis_3) - saves ~10GB RAM, hours of compute
- SKIPS: Mutual Information matrices (analysis_5) - saves ~5GB RAM
- KEEPS: Compensation stats (analysis_1) - lightweight streaming
- KEEPS: Vernier clusters (analysis_2) - essential for designer
- KEEPS: Higher-order rules (analysis_4) - essential for designer

Expected runtime: ~2-4 hours for 12M sequences (vs 2-5 days for full)
Expected memory: ~4GB peak (vs 20GB+ for full)

Usage:
    python vhh_epistasis_imgt_light.py         --npz data/processed/full_12m/full_12m_imgt.npz         --output models/epistasis/imgt_light_v1         --checkpoint-interval 1000000
```

---

## 📁 `active/analysis/vhh_naturalness_analyzer/`

### `vhh_foldability_analyzer.py`
**Version:** 2
**Path:** `active/analysis/vhh_naturalness_analyzer/vhh_foldability_analyzer.py`
**Stats:** 915 lines | 33.1KB | 17 functions | 5 classes
**Modified:** 2025-12-10 21:58

**Description:**
```
VHH CDR-Framework Naturalness Analyzer
======================================

Uses camelid epistasis data as a PRIOR on CDR3-framework compatibility,
NOT as a binary fold/don't fold gate.

Key insight: For humanized scaffolds, foldability is known from experiments.
The question is: "Given this framework, how natural is this CDR3?"

Two modes:
1. Full sequence input: Extracts CDRs and frameworks automatically
2. Framework + CDR input: Uses provided framework regions with grafted CDRs

Output:
- Naturalness score (lower = more natural)
- Per-feature z-scores (length, charge, aromatics)
- Framework annotation (classical VHH, humanized, VH-like)
- Risk assessment based on CDR-framework compatibility
- Suggested optimizations

Usage:
  # Full sequences
  python vhh_naturalness_analyzer.py --input sequences.csv --epistasis epistasis_v2_full.pkl
  
  # Framework library + CDRs
  python vhh_naturalness_analyzer.py --input designs.csv --frameworks VHH_frameworks.xlsx --epistasis epistasis_v2_full.pkl
```

---

### `vhh_naturalness_analyzer_v3.py`
**Version:** 2
**Path:** `active/analysis/vhh_naturalness_analyzer/vhh_naturalness_analyzer_v3.py`
**Stats:** 1036 lines | 37.9KB | 17 functions | 5 classes
**Modified:** 2026-01-02 01:03

**Description:**
```
VHH CDR-Framework Naturalness Analyzer
======================================

Uses camelid epistasis data as a PRIOR on CDR3-framework compatibility,
NOT as a binary fold/don't fold gate.

Key insight: For humanized scaffolds, foldability is known from experiments.
The question is: "Given this framework, how natural is this CDR3?"

Two modes:
1. Full sequence input: Extracts CDRs and frameworks automatically
2. Framework + CDR input: Uses provided framework regions with grafted CDRs

Output:
- Naturalness score (lower = more natural)
- Per-feature z-scores (length, charge, aromatics)
- Framework annotation (classical VHH, humanized, VH-like)
- Risk assessment based on CDR-framework compatibility
- Suggested optimizations

Usage:
  # Full sequences
  python vhh_naturalness_analyzer.py --input sequences.csv --epistasis epistasis_v2_full.pkl
  
  # Framework library + CDRs
  python vhh_naturalness_analyzer.py --input designs.csv --frameworks VHH_frameworks.xlsx --epistasis epistasis_v2_full.pkl
```

---

## 📁 `active/database/`

### `npz_fullscan_v6_integrated.py`
**Version:** 6
**Path:** `active/database/npz_fullscan_v6_integrated.py`
**Stats:** 1060 lines | 39.5KB | 20 functions | 1 classes
**Modified:** 2026-01-02 00:57

**Description:**
```
NPZ FULLSCAN v6 INTEGRATED - Combined Interactive Analysis & Full Scan

This version integrates:
1. Interactive threshold analysis (sampling mode)
2. Full exhaustive scan with all v6 improvements

NEW: Asks if you want to sample first to determine optimal thresholds!

Usage:
  python3 npz_fullscan_v6_integrated.py --db-root /path/to/db --query-seq "EIQLQQ..."
  
  The script will ask if you want to:
  1. Run a quick analysis to determine optimal thresholds (recommended)
  2. Proceed directly with full scan
```

---

## 📁 `active/utilities/`

### `dna_translator.py`
**Path:** `active/utilities/dna_translator.py`
**Stats:** 312 lines | 11.3KB | 8 functions | 0 classes
**Modified:** 2025-12-08 15:56

**Description:**
```
DNA to Amino Acid Translator
Translates DNA sequences in all 6 reading frames and identifies the best one.
```

---

### `pathlib_list_v2.py`
**Path:** `active/utilities/pathlib_list_v2.py`
**Stats:** 60 lines | 1.8KB | 1 functions | 0 classes
**Modified:** 2026-01-02 19:21

**Description:**
```
(No documentation found)
```

---

### `pull_cdrs.py`
**Path:** `active/utilities/pull_cdrs.py`
**Stats:** 530 lines | 15.8KB | 18 functions | 1 classes
**Modified:** 2026-01-02 01:53

**Description:**
```
(No documentation found)
```

---

### `scfv_anarci.py`
**Path:** `active/utilities/scfv_anarci.py`
**Stats:** 100 lines | 3.9KB | 4 functions | 0 classes
**Modified:** 2026-01-18 14:53

**Description:**
```
JD
```

---

## 📁 `archive/correlation_analysis/`

### `analyze_cdr_framework_advanced.py`
**Path:** `archive/correlation_analysis/analyze_cdr_framework_advanced.py`
**Stats:** 738 lines | 27.9KB | 9 functions | 0 classes
**Modified:** 2025-12-01 21:14

**Description:**
```
Advanced VHH CDR-Framework Correlation Analyzer
================================================
Accounts for germline preferences by:
1. Clustering sequences by framework similarity (proxy for V-gene family)
2. Computing associations WITHIN each cluster
3. Identifying UNIVERSAL associations (consistent across clusters)
4. Flagging GERMLINE-SPECIFIC associations (only in certain clusters)

Additional analyses:
- Vernier zone / contact position flagging
- Mutual information for position pairs
- Conditional probability p(FR | CDR) estimates

Usage:
    python analyze_cdr_framework_advanced.py /path/to/shards/*.npz --output ./results

Output:
    - correlation_results_advanced.pkl: Full analysis with germline control
    - correlation_summary_advanced.txt: Human-readable summary
    - germline_vs_universal.txt: Which rules are germline-specific vs universal
```

---

### `analyze_cdr_framework_correlations.py`
**Path:** `archive/correlation_analysis/analyze_cdr_framework_correlations.py`
**Stats:** 490 lines | 17.8KB | 6 functions | 0 classes
**Modified:** 2025-12-01 21:14

**Description:**
```
VHH CDR-Framework Correlation Analyzer
======================================
Processes NPZ shards to find correlations between:
- Single amino acid positions (CDR pos X ↔ FW pos Y)
- Triplet motifs (CDR triplet ↔ FW motif)
- Lengths (CDR length ↔ FW length)

Usage:
    python analyze_cdr_framework_correlations.py /path/to/shards/*.npz
    
    # Or specify output directory:
    python analyze_cdr_framework_correlations.py --output ./results /path/to/shards/*.npz

Output:
    - correlation_results.pkl: Raw counts for all correlations
    - correlation_summary.txt: Human-readable summary
    - (Optional) Visualizations if --visualize flag is used
```

---

### `vhh_global_compensation_analysis.py`
**Path:** `archive/correlation_analysis/vhh_global_compensation_analysis.py`
**Stats:** 659 lines | 24.0KB | 19 functions | 4 classes
**Modified:** 2025-12-06 21:40

**Description:**
```
VHH Global Compensation Analysis
================================

Following ChatGPT's recommendations for finding "counterbalance" patterns:
1. Feature-level compensation scans (FW position vs CDR global features)
2. Vernier zone clustering + CDR profiling per cluster
3. Mutual information between positions

Input: NPZ files with ANARCI-numbered VHH sequences
Output: Compensation rules, Vernier clusters, MI matrices
```

---

## 📁 `archive/database_builders/`

### `build_camel_vhh_db.py`
**Path:** `archive/database_builders/build_camel_vhh_db.py`
**Stats:** 715 lines | 23.2KB | 16 functions | 0 classes
**Modified:** 2025-11-25 20:49

**Description:**
```
Build Clean Camelid VHH Database from OAS Nucleotide Data

This script:
1. Loads camel heavy-chain nucleotide sequences from an OAS TSV/CSV export
2. Translates them to amino acids (properly, from raw nucleotides)
3. Renumbers them with ANARCI (scheme="imgt")
4. Extracts full V-domain and CDR1/2/3
5. Saves a clean output (CSV and optionally NPZ) for downstream similarity search

This bypasses the truncated sequence_alignment_aa column in OAS and gives you
complete V-domains with proper IMGT numbering.

Usage:
    python build_camel_vhh_db.py         --input-tsv camel_heavy_oas.tsv         --output-prefix camel_vhh_anarci         --batch-size 2000         --min-nt-len 150

Author: Claude (Anthropic)
Date: 2025
```

---

### `build_camel_vhh_db_one_step.py`
**Path:** `archive/database_builders/build_camel_vhh_db_one_step.py`
**Stats:** 787 lines | 26.4KB | 15 functions | 0 classes
**Modified:** 2025-11-25 22:43

**Description:**
```
OAS Camel VHH Complete Pipeline

This script processes OAS CSV/TSV files containing camel heavy-chain nucleotide 
sequences and creates a clean VHH database with proper ANARCI numbering.

Pipeline steps:
1. Load all CSV/TSV files from a directory (or single file)
2. Extract nucleotide sequences from first column or 'sequence' column
3. Translate nucleotides → amino acids
4. Run ANARCI for IMGT numbering
5. Extract CDR1, CDR2, CDR3 and full V-domain
6. Save clean CSV and NPZ outputs

Usage:
    # Process all CSV files in a directory
    python process_camel_vhh_pipeline.py         --input-dir /path/to/csv/files/         --output-prefix camel_vhh_clean         --batch-size 2000

    # Process a single file
    python process_camel_vhh_pipeline.py         --input-file /path/to/file.csv         --output-prefix camel_vhh_clean

Author: Claude (Anthropic)
Date: 2025
```

---

### `build_models_only.py`
**Path:** `archive/database_builders/build_models_only.py`
**Stats:** 322 lines | 9.4KB | 9 functions | 1 classes
**Modified:** 2025-12-05 10:08

**Description:**
```
Build models from existing shards WITHOUT recomputing correlations.
Much faster - just samples 100k sequences for model training.

Usage:
    pip install scikit-learn
    python build_models_only.py VHH_shards/*.npz --output ./results/2_models
```

---

### `build_pfr_cdr_models.py`
**Path:** `archive/database_builders/build_pfr_cdr_models.py`
**Stats:** 565 lines | 19.1KB | 9 functions | 0 classes
**Modified:** 2025-12-01 22:04

**Description:**
```
Conditional p(FR | CDR) Model Builder
=====================================
Trains predictive models for each framework position given CDR features.

For each FR position, answers:
- "Given these CDRs, what amino acid does nature choose here?"
- "How 'weird' is residue X at this position for these CDRs?"

Features used:
- CDR lengths (CDR1, CDR2, CDR3)
- CDR physicochemical properties (charge, hydrophobicity, aromaticity)
- CDR boundary amino acids (first/last 3 of each CDR)

Usage:
    # Build models from NPZ files
    python build_pfr_cdr_models.py /path/to/shards/*.npz --output ./models
    
    # Then use score_framework.py to evaluate your framework

Output:
    - pfr_cdr_models.pkl: Trained models for each FR position
    - model_summary.txt: Feature importances and model quality metrics
```

---

## 📁 `archive/epistasis_pipeline/v1_20251201/`

### `vhh_full_epistasis_overnight_v1.py`
**Version:** 2
**Path:** `archive/epistasis_pipeline/v1_20251201/vhh_full_epistasis_overnight_v1.py`
**Stats:** 1172 lines | 44.3KB | 30 functions | 7 classes
**Modified:** 2025-12-06 22:17

**Description:**
```
VHH Full Epistasis Analysis - OVERNIGHT VERSION v2
===================================================

Fixed based on ChatGPT's feedback:
1. Fixed cond_results NameError bug
2. Streaming stats for FeatureCompensationScanner (Welford's algorithm)
3. Integer encoding for MI analyzer (memory efficient)
4. Subsample cap for MI (default 2M sequences)
5. CLI flags for light vs full mode
6. Removed dead code in MI analyzer
7. Memory-efficient HigherOrderRuleFinder (streaming counts)

Implements ALL of ChatGPT's recommendations at maximum depth.
```

---

## 📁 `archive/epistasis_pipeline/v2_20251205/`

### `vhh_epistasis_overnight_v2.py`
**Version:** 2
**Path:** `archive/epistasis_pipeline/v2_20251205/vhh_epistasis_overnight_v2.py`
**Stats:** 1262 lines | 47.9KB | 31 functions | 7 classes
**Modified:** 2025-12-07 23:34

**Description:**
```
VHH Full Epistasis Analysis - OVERNIGHT VERSION v2 (CORRECTED)
===============================================================

CORRECTIONS FROM v1:
1. FIXED POSITION MAPPING:
   - FR2_2 (index 1) = IMGT 42 = Y/F hallmark (was incorrectly FR2_4)
   - FR2_9 (index 8) = IMGT 49 = E/Q hallmark
   - FR2_10 (index 9) = IMGT 50 = R/L hallmark (was incorrectly FR2_12)
   - FR2_12 (index 11) = IMGT 52 = G/W hallmark (was incorrectly FR2_14)

2. CORRECTED FAMILY CLASSIFICATION:
   - Classical VHH: pos42=F/Y AND pos49=E/Q AND pos50=R
   - Y_C2, F_C2, F_C4: Based on pos42 + total cysteine count
   - VH_like: pos50=L (includes humanized VHHs)
   - VHH_W52: has W at pos52
   - Non_classical: everything else

3. Uses TOTAL cysteine count (not just CDR3) for C2/C4 classification

RETAINED FROM v1:
- Fixed cond_results NameError bug
- Streaming stats for FeatureCompensationScanner (Welford's algorithm)
- Integer encoding for MI analyzer (memory efficient)
- Subsample cap for MI (default 2M sequences)
- CLI flags for light vs full mode
- Progress bars with tqdm

Position Mapping Reference:
  FR2: W  F  R  Q  x  P  x  x  E  R  x  G  L  x
  Idx: 0  1  2  3  4  5  6  7  8  9  10 11 12 13
  IMGT:41 42 43 44 45 46 47 48 49 50 51 52 53 54
```

---

## 📁 `archive/full_analysis/`

### `run_full_analysis_fast2_v3.py`
**Path:** `archive/full_analysis/run_full_analysis_fast2_v3.py`
**Stats:** 920 lines | 33.8KB | 14 functions | 0 classes
**Modified:** 2025-12-01 22:12

**Description:**
```
VHH CDR-Framework Complete Analysis Pipeline (OPTIMIZED)
=========================================================
Vectorized version - 10-50x faster than original

Usage:
    python run_full_analysis_fast.py /path/to/shards/*.npz --output ./results
```

---

### `run_full_analysis_fast_v2.py`
**Path:** `archive/full_analysis/run_full_analysis_fast_v2.py`
**Stats:** 878 lines | 32.0KB | 14 functions | 0 classes
**Modified:** 2025-12-01 21:56

**Description:**
```
VHH CDR-Framework Complete Analysis Pipeline (OPTIMIZED)
=========================================================
Vectorized version - 10-50x faster than original

Usage:
    python run_full_analysis_fast.py /path/to/shards/*.npz --output ./results
```

---

### `run_full_analysis_lowmem_v4.py`
**Path:** `archive/full_analysis/run_full_analysis_lowmem_v4.py`
**Stats:** 740 lines | 26.5KB | 18 functions | 3 classes
**Modified:** 2025-12-04 23:29

**Description:**
```
LOW MEMORY VERSION - CDR-Framework Analysis Pipeline

Processes shards ONE AT A TIME, accumulating statistics without
holding all sequences in memory. Suitable for large datasets on
memory-constrained systems (e.g., WSL with limited RAM).

Usage:
    python run_full_analysis_lowmem.py /path/to/shards/*.npz --output ./results
    python run_full_analysis_lowmem.py /path/to/shards/*.npz --output ./results --max-shards 10
```

---

### `run_full_analysis_lowmem_v5.py`
**Version:** 3
**Path:** `archive/full_analysis/run_full_analysis_lowmem_v5.py`
**Stats:** 797 lines | 29.5KB | 15 functions | 2 classes
**Modified:** 2025-12-05 18:36

**Description:**
```
LOW MEMORY VERSION v3 - Germline-Aware CDR-Framework Analysis

Features:
- Germline classification using FR2 hallmark positions (42/49/50/52)
- Cysteine counting for C2 vs C4 classification
- Species/germline family grouping (Alpaca Y_C2/F_C2/F_C4, Llama VHH1-5, Camel, VH-like)
- TRUE CDR-dependence analysis WITHIN germline families
- Separates germline effects from CDR-framework compatibility

Usage:
    python run_full_analysis_lowmem_v3.py /path/to/shards/*.npz --output ./results_v3
```

---

### `run_full_analysis_v1.py`
**Path:** `archive/full_analysis/run_full_analysis_v1.py`
**Stats:** 974 lines | 37.1KB | 13 functions | 0 classes
**Modified:** 2025-12-01 21:13

**Description:**
```
VHH CDR-Framework Complete Analysis Pipeline
=============================================
Master script that runs all analyses sequentially.

Usage:
    python run_full_analysis.py /path/to/shards/*.npz --output ./results
    
    # With custom framework to score:
    python run_full_analysis.py /path/to/shards/*.npz --output ./results         --fr1 "EVQLVESGGGLVQPGGSLRLSCAAS"         --fr2 "WFRQAPGKGREFVA"         --fr3 "YYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYC"         --fr4 "WGQGTQVTVSS"
    
    # With specific CDRs to score against:
    python run_full_analysis.py /path/to/shards/*.npz --output ./results         --cdr1 "GFTFSSYA" --cdr2 "ISYDGSNK" --cdr3 "ARDLLVRY"

Output structure:
    ./results/
    ├── 1_correlations/
    │   ├── correlation_results.pkl
    │   ├── correlation_summary.txt
    │   └── figures/
    ├── 2_advanced_correlations/
    │   ├── correlation_results_advanced.pkl
    │   ├── correlation_summary_advanced.txt
    │   └── germline_vs_universal.txt
    ├── 3_pfr_cdr_models/
    │   ├── pfr_cdr_models.pkl
    │   ├── model_summary.txt
    │   └── figures/
    ├── 4_framework_scoring/
    │   └── framework_report.txt
    └── FINAL_REPORT.md
```

---

## 📁 `archive/naturalness_analyzer/`

### `score_framework.py`
**Path:** `archive/naturalness_analyzer/score_framework.py`
**Stats:** 460 lines | 15.5KB | 8 functions | 0 classes
**Modified:** 2025-12-01 21:13

**Description:**
```
Framework Scorer
================
Evaluates a universal framework against p(FR|CDR) models.

For each FR position, shows:
- How well your framework residue matches the repertoire preference
- Alternative residues that might work better for specific CDR panels
- Positions that are "out of distribution" for your target CDRs

Usage:
    # Score against all training data
    python score_framework.py --models pfr_cdr_models.pkl

    # Score for specific CDR sequences
    python score_framework.py --models pfr_cdr_models.pkl         --cdr1 "GFTFSSYA" --cdr2 "ISYDGSNK" --cdr3 "ARDLLVRY"

    # Score for a CDR panel (file with one CDR set per line)
    python score_framework.py --models pfr_cdr_models.pkl         --cdr-file my_cdrs.txt

    # Custom framework (default is Vincke universal)
    python score_framework.py --models pfr_cdr_models.pkl         --fr1 "EVQLVESGGGLVQPGGSLRLSCAAS"         --fr2 "WFRQAPGKGREFVA"         --fr3 "YYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYC"         --fr4 "WGQGTQVTVSS"
```

---

### `vhh_framework_optimizer.py`
**Version:** 2
**Path:** `archive/naturalness_analyzer/vhh_framework_optimizer.py`
**Stats:** 503 lines | 17.5KB | 10 functions | 1 classes
**Modified:** 2025-12-06 20:56

**Description:**
```
VHH Framework Optimizer v2 - Fixed CDR extraction and scoring

Key fixes from v1:
1. Better CDR extraction (excludes FR1 anchor from CDR1)
2. Higher default confidence threshold (90%)
3. Ignores cdr1[0:3] rules (often contaminated by FR anchor)
4. Added user's custom scaffold
```

---

## 📁 `archive/naturalness_analyzer/v1_20251210/`

### `vhh_naturalness_analyzer.py`
**Version:** 2
**Path:** `archive/naturalness_analyzer/v1_20251210/vhh_naturalness_analyzer.py`
**Stats:** 844 lines | 30.3KB | 17 functions | 5 classes
**Modified:** 2026-01-02 01:18

**Description:**
```
VHH CDR-Framework Naturalness Analyzer
======================================

Uses camelid epistasis data as a PRIOR on CDR3-framework compatibility,
NOT as a binary fold/don't fold gate.

Key insight: For humanized scaffolds, foldability is known from experiments.
The question is: "Given this framework, how natural is this CDR3?"

Two modes:
1. Full sequence input: Extracts CDRs and frameworks automatically
2. Framework + CDR input: Uses provided framework regions with grafted CDRs

Output:
- Naturalness score (lower = more natural)
- Per-feature z-scores (length, charge, aromatics)
- Framework annotation (classical VHH, humanized, VH-like)
- Risk assessment based on CDR-framework compatibility
- Suggested optimizations

Usage:
  # Full sequences
  python vhh_naturalness_analyzer.py --input sequences.csv --epistasis epistasis_v2_full.pkl
  
  # Framework library + CDRs
  python vhh_naturalness_analyzer.py --input designs.csv --frameworks VHH_frameworks.xlsx --epistasis epistasis_v2_full.pkl
```

---

## 📁 `archive/naturalness_analyzer/v2_20251212/`

### `vhh_naturalness_analyzer_v2.py`
**Version:** 2
**Path:** `archive/naturalness_analyzer/v2_20251212/vhh_naturalness_analyzer_v2.py`
**Stats:** 915 lines | 33.1KB | 17 functions | 5 classes
**Modified:** 2026-01-02 01:17

**Description:**
```
VHH CDR-Framework Naturalness Analyzer
======================================

Uses camelid epistasis data as a PRIOR on CDR3-framework compatibility,
NOT as a binary fold/don't fold gate.

Key insight: For humanized scaffolds, foldability is known from experiments.
The question is: "Given this framework, how natural is this CDR3?"

Two modes:
1. Full sequence input: Extracts CDRs and frameworks automatically
2. Framework + CDR input: Uses provided framework regions with grafted CDRs

Output:
- Naturalness score (lower = more natural)
- Per-feature z-scores (length, charge, aromatics)
- Framework annotation (classical VHH, humanized, VH-like)
- Risk assessment based on CDR-framework compatibility
- Suggested optimizations

Usage:
  # Full sequences
  python vhh_naturalness_analyzer.py --input sequences.csv --epistasis epistasis_v2_full.pkl
  
  # Framework library + CDRs
  python vhh_naturalness_analyzer.py --input designs.csv --frameworks VHH_frameworks.xlsx --epistasis epistasis_v2_full.pkl
```

---

## 📁 `archive/npz_scanner/v2_20251103/`

### `npz_fullscan_v2.py`
**Path:** `archive/npz_scanner/v2_20251103/npz_fullscan_v2.py`
**Stats:** 748 lines | 28.9KB | 16 functions | 1 classes
**Modified:** 2025-11-07 23:07

**Description:**
```
NPZ Fullscan with Corrected IMGT CDR Extraction
Fixed to properly extract CDRs according to strict IMGT boundaries
CDR-H2 should be ILPGSGST not QILPGSGSTK
```

---

## 📁 `archive/npz_scanner/v3_20251103/`

### `npz_fullscan_v3.py`
**Version:** 3
**Path:** `archive/npz_scanner/v3_20251103/npz_fullscan_v3.py`
**Stats:** 502 lines | 17.0KB | 14 functions | 1 classes
**Modified:** 2025-11-11 13:57

**Description:**
```
NPZ FULLSCAN v3 - FIXED with ASCII decoding
- Live dual progress bars (global shards + per-shard mini-bar)
- Direct CDR extraction from numbering arrays (no ANARCI on DB entries)
- Smart FULL-region optimization
- Persistent per-shard stats lines (hits, max ID, top 25% cutoff, skipped, time)
- Per-run CSV summary (timestamped) saved to output dir

Usage example:
python -u npz_fullscan_v3_fixed.py   --db-root "/home/user/KA-Search/extracted/oas-paper/oasdb_20230111/Heavy/Human"   --query-seq "EIQLQQ...TVSA"   --regions full   --min-id-full 0.30   --outdir "/home/user/KA-Search/runs_fullscan"   --tag "full_30_human"
```

---

## 📁 `archive/npz_scanner/v5_20251104/`

### `npz_fullscan_v5.py`
**Version:** 5
**Path:** `archive/npz_scanner/v5_20251104/npz_fullscan_v5.py`
**Stats:** 1075 lines | 38.7KB | 22 functions | 2 classes
**Modified:** 2025-11-11 22:41

**Description:**
```
NPZ FULLSCAN v5 - Complete Feature Set

NEW in v5:
- Sequence length pre-filtering for full scans (automatic speedup)
- CDR length pre-filtering for full scans (automatic when using full mode)
- Framework extraction and searching
- Position-specific amino acid constraints with similarity groups

Features:
- Live dual progress bars
- Direct CDR/FR extraction from numbering arrays
- Multiple filtering strategies for speed
- Position-specific matching constraints
- Persistent per-shard stats
- Per-run CSV summary

Usage examples:

# Full sequence with automatic optimizations:
python3 npz_fullscan_v5.py   --db-root "/path/to/Heavy/Human"   --query-seq "EIQLQQ..."   --regions full   --min-id-full 0.30   --outdir "./results"   --tag "full_optimized"

# CDRs with position constraints:
python3 npz_fullscan_v5.py   --db-root "/path/to/Heavy/Camel"   --query-seq "EIQLQQ..."   --regions cdr1,cdr2,cdr3   --min-id-cdr1 0.35 --min-id-cdr2 0.35 --min-id-cdr3 0.40   --position-constraint "cdr3:4:H:similar"   --position-constraint "cdr3:7:C:exact"   --outdir "./results"   --tag "cdr_with_constraints"

# Framework search:
python3 npz_fullscan_v5.py   --db-root "/path/to/Heavy/Human"   --query-seq "EIQLQQ..."   --regions frameworks   --min-id-frameworks 0.35   --outdir "./results"   --tag "framework_search"
```

---

## 📁 `archive/npz_scanner/v6_20251105/`

### `npz_fullscan_v6_interactive_v6.py`
**Version:** 6
**Path:** `archive/npz_scanner/v6_20251105/npz_fullscan_v6_interactive_v6.py`
**Stats:** 2365 lines | 94.7KB | 23 functions | 1 classes
**Modified:** 2025-11-23 15:35

**Description:**
```
NPZ FULLSCAN v6 INTERACTIVE - Enhanced Interactive Mode
Full-featured antibody search with interactive setup

FEATURES:
- Run without arguments for fully interactive mode
- Interactive threshold analysis (sampling mode)
- Full exhaustive scan with all v6 improvements
- Automatic recommendation of optimal thresholds

USAGE:
1. Interactive mode (just hit enter):
   python3 npz_fullscan_v6_interactive.py
   
2. Command line mode (as before):
   python3 npz_fullscan_v6_interactive.py --db-root /path/to/db --query-seq "EIQLQQ..."
```

---

## 📁 `archive/npz_scanner/v7_20251106/`

### `npz_fullscan_v7_interactive_v7.py`
**Version:** 6
**Path:** `archive/npz_scanner/v7_20251106/npz_fullscan_v7_interactive_v7.py`
**Stats:** 2574 lines | 102.9KB | 26 functions | 1 classes
**Modified:** 2025-11-24 22:01

**Description:**
```
NPZ FULLSCAN v6 INTERACTIVE - Enhanced Interactive Mode
Full-featured antibody search with interactive setup

FEATURES:
- Run without arguments for fully interactive mode
- Interactive threshold analysis (sampling mode)
- Full exhaustive scan with all v6 improvements
- Automatic recommendation of optimal thresholds

USAGE:
1. Interactive mode (just hit enter):
   python3 npz_fullscan_v6_interactive.py
   
2. Command line mode (as before):
   python3 npz_fullscan_v6_interactive.py --db-root /path/to/db --query-seq "EIQLQQ..."

FIXES IN V7:
- Fixed ANARCI unpacking error
- Added query file memory to avoid retyping
- Fixed database paths for extracted location
```

---

## 📁 `archive/npz_scanner/v8_20251107/`

### `npz_fullscan_v8_vhh_v8.py`
**Version:** 8
**Path:** `archive/npz_scanner/v8_20251107/npz_fullscan_v8_vhh_v8.py`
**Stats:** 3117 lines | 120.0KB | 50 functions | 1 classes
**Modified:** 2025-11-29 23:18

**Description:**
```
NPZ FULLSCAN v8 - VHH Database with Target/Patent Metadata
Enhanced antibody search supporting both OAS and VHH unified databases

NEW IN V8:
- Support for VHH unified database shards (12M+ sequences)
- Displays target/antigen binding information  
- Shows patent IDs and titles for patent-derived sequences
- Source tracking (OAS_Camel, INDI_patent, SAbDab, TheraSAbDab, etc.)
- Shard-aware searching with shard_index.json support
- Annotated-first search option (search high-value sequences first)

DATABASE FORMATS SUPPORTED:
1. OAS format: NPZ with "numberings" array (IMGT-numbered integer arrays)
2. VHH format: NPZ with string arrays (aa_v_full, cdr1, cdr2, cdr3, targets, etc.)

USAGE:
1. Interactive mode:
   python3 npz_fullscan_v8_vhh.py
   
2. VHH database search:
   python3 npz_fullscan_v8_vhh.py --vhh-shards /path/to/VHH_shards --query-seq "EIQLQQ..."
   
3. Search annotated sequences only (fast target lookup):
   python3 npz_fullscan_v8_vhh.py --vhh-shards /path/to/VHH_shards --annotated-only --query-seq "..."
```

---

## 📁 `archive/npz_scanner/v9_20251108/`

### `npz_fullscan_v9_vhh_v9.py`
**Version:** 8
**Path:** `archive/npz_scanner/v9_20251108/npz_fullscan_v9_vhh_v9.py`
**Stats:** 3160 lines | 121.9KB | 51 functions | 1 classes
**Modified:** 2025-11-30 14:07

**Description:**
```
NPZ FULLSCAN v8 - VHH Database with Target/Patent Metadata
Enhanced antibody search supporting both OAS and VHH unified databases

NEW IN V8:
- Support for VHH unified database shards (12M+ sequences)
- Displays target/antigen binding information  
- Shows patent IDs and titles for patent-derived sequences
- Source tracking (OAS_Camel, INDI_patent, SAbDab, TheraSAbDab, etc.)
- Shard-aware searching with shard_index.json support
- Annotated-first search option (search high-value sequences first)

DATABASE FORMATS SUPPORTED:
1. OAS format: NPZ with "numberings" array (IMGT-numbered integer arrays)
2. VHH format: NPZ with string arrays (aa_v_full, cdr1, cdr2, cdr3, targets, etc.)

USAGE:
1. Interactive mode:
   python3 npz_fullscan_v8_vhh.py
   
2. VHH database search:
   python3 npz_fullscan_v8_vhh.py --vhh-shards /path/to/VHH_shards --query-seq "EIQLQQ..."
   
3. Search annotated sequences only (fast target lookup):
   python3 npz_fullscan_v8_vhh.py --vhh-shards /path/to/VHH_shards --annotated-only --query-seq "..."
```

---

## 📁 `archive/one_off/`

### `debug_npz_format.py`
**Path:** `archive/one_off/debug_npz_format.py`
**Stats:** 59 lines | 1.9KB | 1 functions | 0 classes
**Modified:** 2025-12-06 21:48

**Description:**
```
Debug script to check NPZ file format.
```

---

### `diagnose_npz.py`
**Path:** `archive/one_off/diagnose_npz.py`
**Stats:** 29 lines | 0.9KB | 0 functions | 0 classes
**Modified:** 2025-12-07 23:24

**Description:**
```
Quick diagnostic for NPZ structure
```

---

### `diagnose_oas_aligned.py`
**Path:** `archive/one_off/diagnose_oas_aligned.py`
**Stats:** 249 lines | 8.7KB | 4 functions | 0 classes
**Modified:** 2025-11-24 22:22

**Description:**
```
Advanced diagnostic for OAS-aligned NPZ files.
Based on the paper's description of 200 unique positions in canonical alignment.
```

---

### `diagnose_oas_data.py`
**Path:** `archive/one_off/diagnose_oas_data.py`
**Stats:** 339 lines | 11.7KB | 4 functions | 0 classes
**Modified:** 2025-11-25 20:49

**Description:**
```
OAS/KA-Search Data Diagnostic Script

This script examines your data directory and determines:
1. What format the data is in (NPZ vs TSV/CSV)
2. Whether it contains raw nucleotide sequences
3. What the sequence lengths look like
4. Whether the data is truncated or complete

Run this on your system:
    python diagnose_oas_data.py --data-dir /path/to/your/data
```

---

### `download_oas_camel.py`
**Path:** `archive/one_off/download_oas_camel.py`
**Stats:** 214 lines | 6.2KB | 3 functions | 0 classes
**Modified:** 2025-11-25 20:49

**Description:**
```
Download OAS Camel Heavy Chain Data

This script downloads raw camel heavy-chain sequences from the OAS bulk API.
It downloads the sequences with nucleotide data, which can then be processed
by build_camel_vhh_db.py to create a clean VHH database.

Usage:
    python download_oas_camel.py --output-dir ./oas_camel_data

The script will download all available camel heavy-chain data units from OAS.
This may take considerable time depending on the size of the data.

Note: OAS data is organized by study. Each study is downloaded as a separate file.
```

---

### `interactive_comprehensive_cdr_analysis.py`
**Path:** `archive/one_off/interactive_comprehensive_cdr_analysis.py`
**Stats:** 923 lines | 32.9KB | 18 functions | 1 classes
**Modified:** 2025-11-11 19:56

**Description:**
```
Interactive Comprehensive CDR Analysis with Threshold Recommendations

This script:
1. Asks you questions about your search preferences
2. Runs a comprehensive analysis sampling sequences from the database
3. Calculates optimal thresholds to achieve your target hit count
4. Provides ready-to-run KA-Search commands
```

---

### `process_camel_vhh_pipeline.py`
**Path:** `archive/one_off/process_camel_vhh_pipeline.py`
**Stats:** 787 lines | 26.4KB | 15 functions | 0 classes
**Modified:** 2025-11-25 22:43

**Description:**
```
OAS Camel VHH Complete Pipeline

This script processes OAS CSV/TSV files containing camel heavy-chain nucleotide 
sequences and creates a clean VHH database with proper ANARCI numbering.

Pipeline steps:
1. Load all CSV/TSV files from a directory (or single file)
2. Extract nucleotide sequences from first column or 'sequence' column
3. Translate nucleotides → amino acids
4. Run ANARCI for IMGT numbering
5. Extract CDR1, CDR2, CDR3 and full V-domain
6. Save clean CSV and NPZ outputs

Usage:
    # Process all CSV files in a directory
    python process_camel_vhh_pipeline.py         --input-dir /path/to/csv/files/         --output-prefix camel_vhh_clean         --batch-size 2000

    # Process a single file
    python process_camel_vhh_pipeline.py         --input-file /path/to/file.csv         --output-prefix camel_vhh_clean

Author: Claude (Anthropic)
Date: 2025
```

---

### `process_indi_merge_final.py`
**Path:** `archive/one_off/process_indi_merge_final.py`
**Stats:** 452 lines | 17.1KB | 8 functions | 0 classes
**Modified:** 2025-11-28 15:27

**Description:**
```
Process INDI database and merge with existing VHH unified database.

Creates a comprehensive VHH database with:
- Detailed source tracking (OAS_Camel, INDI_patent, INDI_structure, SAbDab_PDB, etc.)
- Target/binding partner information where available
- Patent numbers for patent-derived sequences
- Deduplication keeping entries with most metadata

INDI files:
- patent_sequence.tsv + patent_meta.tsv (targets from biomolecules)
- structure_sequence.tsv + structure_meta.tsv (PDB info)
- abgenbank_sequence.tsv (no metadata file)
- manual_sequence.tsv + manual_meta.tsv (curated entries)
- ngs_sequence.tsv (SKIP - NGS repertoire data, no functional annotation)
```

---

### `process_sabdab_pdb.py`
**Path:** `archive/one_off/process_sabdab_pdb.py`
**Stats:** 435 lines | 15.1KB | 9 functions | 0 classes
**Modified:** 2025-11-27 00:36

**Description:**
```
Extract VHH sequences from SAbDab IMGT-numbered PDB files and merge with unified database.

This script:
1. Reads sabdab_nano_summary_all.tsv for metadata
2. Extracts sequences from IMGT-numbered PDB files
3. Extracts CDR regions using IMGT positions
4. Optionally merges with existing unified VHH database
```

---

### `process_vhh_with_antigens.py`
**Path:** `archive/one_off/process_vhh_with_antigens.py`
**Stats:** 486 lines | 15.4KB | 12 functions | 0 classes
**Modified:** 2025-11-26 23:45

**Description:**
```
Process VHH Sequences from GenBank/PLAbDab and Merge with OAS Camel Data

This script:
1. Loads VHH sequences from CSV (GenBank/PLAbDab format)
2. Extracts CDRs (uses provided CDRs or re-runs ANARCI)
3. Preserves antigen/target information
4. Optionally merges with OAS camel data
5. Creates unified CSV and NPZ outputs

Key feature: Preserves antigen info for display in search results!

Usage:
    # Process VHH CSV only
    python process_vhh_with_antigens.py         --vhh-csv /path/to/vhh_sequences.csv         --output-prefix vhh_annotated

    # Merge with OAS camel data
    python process_vhh_with_antigens.py         --vhh-csv /path/to/vhh_sequences.csv         --oas-csv /path/to/camel_vhh_clean_anarci_renumbered.csv         --output-prefix VHH_db_unified

Author: Claude (Anthropic)
Date: 2025
```

---

### `shard_database.py`
**Path:** `archive/one_off/shard_database.py`
**Stats:** 246 lines | 8.1KB | 6 functions | 0 classes
**Modified:** 2025-11-28 23:58

**Description:**
```
Shard the VHH database for efficient KA-Search queries.

Creates multiple NPZ shards from a single large database file,
plus an index file for coordinating searches across shards.
```

---

## 📁 `archive/visualizations/`

### `visualize_correlations.py`
**Path:** `archive/visualizations/visualize_correlations.py`
**Stats:** 358 lines | 14.1KB | 7 functions | 0 classes
**Modified:** 2025-12-01 21:14

**Description:**
```
Visualization Script for CDR-Framework Correlations
===================================================
Creates comprehensive visualizations from correlation_results.pkl

Usage:
    python visualize_correlations.py correlation_results.pkl
    
    # Or specify output directory:
    python visualize_correlations.py correlation_results.pkl --output ./figures
```

---

### `visualize_pfr_models.py`
**Path:** `archive/visualizations/visualize_pfr_models.py`
**Stats:** 236 lines | 7.9KB | 7 functions | 0 classes
**Modified:** 2025-12-01 21:13

**Description:**
```
Visualize p(FR|CDR) Model Results
=================================
Creates visualizations from the trained models showing:
1. Feature importance by FR position
2. Most CDR-dependent positions
3. CDR property influences on FR

Usage:
    python visualize_pfr_models.py pfr_cdr_models.pkl --output ./figures
```

---

## 📁 `tools/`

### `dedup_annotated_imgt.py`
**Path:** `tools/dedup_annotated_imgt.py`
**Stats:** 131 lines | 5.1KB | 4 functions | 0 classes
**Modified:** 2026-01-09 18:17

**Description:**
```
(No documentation found)
```

---

### `dedupe_csv_headers.py`
**Path:** `tools/dedupe_csv_headers.py`
**Stats:** 41 lines | 1.2KB | 2 functions | 0 classes
**Modified:** 2026-01-07 22:25

**Description:**
```
(No documentation found)
```

---

### `diagnose_unknown_files.py`
**Path:** `tools/diagnose_unknown_files.py`
**Stats:** 275 lines | 9.6KB | 5 functions | 0 classes
**Modified:** 2025-12-31 15:37

**Description:**
```
Diagnose Unknown Files in KA-Search Directory
==============================================
Extracts docstrings, function names, and key patterns to help categorize files.

Run this on your local machine:
    python diagnose_unknown_files.py /path/to/KA-Search

Output: unknown_files_report.txt with summaries of each file
```

---

### `kasearch_paths.py`
**Version:** 2
**Path:** `tools/kasearch_paths.py`
**Stats:** 434 lines | 11.9KB | 31 functions | 1 classes
**Modified:** 2025-12-31 21:38

**Description:**
```
KA-Search Paths Configuration
==============================
Central path resolution for all KA-Search scripts.

IMPORT THIS MODULE instead of hardcoding paths!

Usage in your scripts:
    from kasearch_paths import PATHS, get_path, find_file
    
    # Get database path (works before AND after reorganization)
    db_path = PATHS.database_production
    
    # Or use flexible finder
    epistasis_model = find_file('epistasis_v2_full.pkl')
    
    # Get path with fallback
    shards_dir = get_path('shards', fallback='Archive/VHH_shards')

Benefits:
    - Works before AND after reorganization
    - Auto-detects which structure is in use
    - Single place to update if paths change
    - Graceful fallbacks
```

---

### `kasearch_reorganize.py`
**Path:** `tools/kasearch_reorganize.py`
**Stats:** 894 lines | 30.4KB | 16 functions | 0 classes
**Modified:** 2025-12-31 16:45

**Description:**
```
KA-Search Reorganization Script
================================
Moves files from current structure to proposed new organization.

IMPORTANT: Run diagnose_unknown_files.py FIRST to understand unknown files!

Usage:
    python kasearch_reorganize.py --dry-run     # Preview changes
    python kasearch_reorganize.py --execute     # Actually move files
    python kasearch_reorganize.py --backup      # Create backup first, then move

Steps:
1. Creates new directory structure
2. Copies current Archive/ to legacy/ (safety backup)
3. Moves files to new locations
4. Generates CHANGELOGs for script families
5. Creates README files
```

---

### `maintain.py`
**Path:** `tools/maintain.py`
**Stats:** 749 lines | 25.9KB | 17 functions | 2 classes
**Modified:** 2026-01-01 18:31

**Description:**
```
KA-Search Maintain
===========================
One script to rule them all.

Just run:
    python maintain.py

It will automatically:
    1. Clean Zone.Identifier files
    2. Detect new/unorganized files and offer to move them
    3. Delete duplicate files (source files that already exist in destination)
    4. Detect modified scripts and auto-version them
    5. Generate changelog entries from detected changes
    6. Clean up old folders (Archive/, PKL Files/, etc.)
    7. Show project status

Usage:
    python maintain.py              # Full maintenance pass (interactive)
    python maintain.py --auto       # Non-interactive (auto-approve safe changes)
    python maintain.py --status     # Just show status, don't change anything
    python maintain.py --help       # Show this help
```

---

### `pathlib_list_v2.py`
**Path:** `tools/pathlib_list_v2.py`
**Stats:** 60 lines | 1.8KB | 1 functions | 0 classes
**Modified:** 2026-01-02 01:23

**Description:**
```
(No documentation found)
```

---

### `qc_imgt_coverage_fast.py`
**Path:** `tools/qc_imgt_coverage_fast.py`
**Stats:** 105 lines | 2.9KB | 4 functions | 0 classes
**Modified:** 2026-01-09 18:53

**Description:**
```
(No documentation found)
```

---

### `summarize_sweep.py`
**Path:** `tools/summarize_sweep.py`
**Stats:** 37 lines | 1.1KB | 3 functions | 0 classes
**Modified:** 2026-01-13 16:48

**Description:**
```
(No documentation found)
```

---
