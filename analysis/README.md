# Analysis Pipeline

## Add these files:

```
epistasis_analyzer.py     ← Overnight epistasis analysis
compensation_rules.py     ← Compensation/conditional rule extraction  
naturalness_analyzer.py   ← ESM2 naturalness scoring

# Also include the JSON rule files if desired:
rules_v7.json             ← 1,730 conditional rules
archetypes_v7.json        ← Family-specific vernier patterns
```

## From your KA-Search directory:

```bash
cp ~/KA-Search/vhh_epistasis_overnight_final.py ./epistasis_analyzer.py
cp ~/KA-Search/vhh_compensation_rules_v7.py ./compensation_rules.py
cp ~/KA-Search/vhh_naturalness_analyzer.py ./naturalness_analyzer.py

# Optional: include the derived rule files
cp ~/KA-Search/analysis_rules_v7.json ./rules_v7.json
cp ~/KA-Search/analysis_vernier_archetypes_v7.json ./archetypes_v7.json
```

These tools derive statistical patterns from the 12M sequence database.
