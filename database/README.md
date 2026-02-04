# Database Tools

## Add these files:

```
npz_fullscan.py      ← NPZ database scanner (latest version)
annotate_shards.py   ← IMGT annotation pipeline
```

## From your KA-Search directory:

```bash
cp ~/KA-Search/npz_fullscan_v6_interactive.py ./npz_fullscan.py
cp ~/KA-Search/annotate_shards_imgt.py ./annotate_shards.py
```

These tools:
- Scan KA-Search NPZ files for sequence similarity
- Annotate sequences with IMGT numbering
- Extract CDRs and frameworks
