# VHH Designer v7.14 - Change Summary

## Key Changes from v7.12/v7.13

### 1. Within-Track Ranking (`track_rank` field)

Added a new field `track_rank` to enable sorting within the control bucket without mixing objectives. This is a local rank (1-indexed) within each track+hallmark combination.

**Output columns now include:**
- `track_rank`: Position within track (e.g., "best ESM within minimal YQRL")
- `is_immutable`: Flag indicating no post-generation modifications allowed

### 2. Explicit Immutability for Track 0

Minimal hallmark controls are now explicitly flagged as immutable:
- `is_immutable = True` for all Track 0 candidates
- Clear comments: "If this breaks, the hallmark itself is incompatible"
- No compensation rules applied
- No consensus sprinkling

### 3. Distance-First Sorting for Tracks 0-3

Different tracks now have different sorting objectives:

| Track | Sort Order | Rationale |
|-------|-----------|-----------|
| **Track 0** (minimal_hallmark) | `n_mutations` ASC, then `esm_loss` ASC | 4mut beats 5mut unless ESM dramatically worse |
| **Track 1** (single_vernier) | `esm_loss` ASC | Distance is constant (~5-6 mutations) |
| **Track 2** (paired_vernier) | `esm_loss` ASC, tie-break by `n_mutations` ASC | Prefer fewer mutations at same ESM |
| **Track 3** (triplet_vernier) | `esm_loss` ASC | ESM is primary signal |
| **Track 4** (optimized) | `combined_score` DESC | Full optimization |

### 4. Capped Probe Complexity

Reduced probe counts to focus on information quality over volume:

| Track | Old Max | New Max | Rationale |
|-------|---------|---------|-----------|
| Track 1 (single) | 10 | **6** | Focus on top-confidence verniers |
| Track 2 (paired) | 8 | **6** | Only high-confidence pairs |
| Track 3 (triplet) | 4 | **3** | Testing causality, not coverage |

### Expected Output Structure

For a run with 6 hallmarks and 20K generation:
- **~12 Track 0 controls** (2 per hallmark: 4mut + 5mut)
- **~36 Track 1 probes** (6 singles × 6 hallmarks)
- **~36 Track 2 probes** (6 pairs × 6 hallmarks)
- **~18 Track 3 probes** (3 triplets × 6 hallmarks, if available)
- **~99 Track 4 ranked** (per --n-select)

Total: ~150-200 sequences including controls

### New Output Columns

```
rank              # Global rank (-N for exempt, +N for ranked)
track_rank        # Within-track rank (1, 2, 3... within track+hallmark)
design_track      # lead, minimal_hallmark, single_vernier, paired_vernier, triplet_vernier, optimized
ranking_exempt    # True for Track 0-3
track_info        # Details (e.g., "YQRL_4mut", "IMGT69_E>D")
is_immutable      # True for Track 0
```

### Sorting Example (Track 0 - Minimal Hallmark)

Before scoring (after ESM):
```
| track_info  | n_mutations | esm_loss | track_rank |
|-------------|-------------|----------|------------|
| YQRL_4mut   | 4           | 2.31     | 1          |
| YQRL_5mut   | 5           | 2.29     | 2          |
| FERV_4mut   | 4           | 2.15     | 1          |
| FERV_5mut   | 5           | 2.18     | 2          |
```

4mut always ranks higher than 5mut within the same hallmark (distance-first), then sorted by ESM within each mutation count.

## Usage

```bash
python vhh_designer_v7_14.py -i M69.fasta \
  --rules analysis_rules_v7.json \
  --archetypes analysis_vernier_archetypes_v7.json \
  --hallmark-db comprehensive_subfamily_analysis_imgt.xlsx \
  --target-hallmarks AUTO \
  --target-families F_C2 F_C4 Y_C2 Y_C4 \
  --mode original \
  --n-generate 20000 \
  --n-select 99 \
  --no-esmfold
```

## What This Fixes

The "YQRL paradox" is now resolved by separating:
1. **Scientific probe controls** (Tracks 0-3) - hypothesis tests, ranked by distance/ESM
2. **Optimized natural-looking candidates** (Track 4) - compete on combined score

This prevents a single scoring objective from having to serve both purposes.
