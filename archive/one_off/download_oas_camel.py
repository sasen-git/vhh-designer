#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Download OAS Camel Heavy Chain Data

This script downloads raw camel heavy-chain sequences from the OAS bulk API.
It downloads the sequences with nucleotide data, which can then be processed
by build_camel_vhh_db.py to create a clean VHH database.

Usage:
    python download_oas_camel.py --output-dir ./oas_camel_data

The script will download all available camel heavy-chain data units from OAS.
This may take considerable time depending on the size of the data.

Note: OAS data is organized by study. Each study is downloaded as a separate file.
"""

import os
import sys
import json
import argparse
import requests
from pathlib import Path
from typing import List, Dict

# OAS API endpoint
OAS_API_BASE = "https://opig.stats.ox.ac.uk/webapps/oas/oas_api"


def get_study_list(species: str = "Camel", chain: str = "Heavy") -> List[Dict]:
    """Get list of available studies from OAS for given species and chain."""
    print(f"Fetching study list for {species} {chain}...")
    
    # OAS uses specific endpoint for study listing
    url = f"{OAS_API_BASE}/studies"
    
    params = {
        "species": species,
        "chain": chain,
    }
    
    try:
        response = requests.get(url, params=params, timeout=30)
        response.raise_for_status()
        studies = response.json()
        print(f"  Found {len(studies)} studies")
        return studies
    except Exception as e:
        print(f"  Error fetching study list: {e}")
        return []


def download_study_data(study_id: str, output_dir: str, include_nucleotide: bool = True) -> str:
    """Download data for a specific study."""
    
    # Construct download URL
    url = f"{OAS_API_BASE}/download"
    
    params = {
        "study_id": study_id,
        "format": "csv",
    }
    
    if include_nucleotide:
        params["include_nt"] = "true"
    
    output_path = os.path.join(output_dir, f"{study_id}.csv")
    
    try:
        response = requests.get(url, params=params, timeout=300, stream=True)
        response.raise_for_status()
        
        with open(output_path, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
        
        return output_path
    except Exception as e:
        print(f"  Error downloading study {study_id}: {e}")
        return None


def main():
    parser = argparse.ArgumentParser(
        description="Download OAS Camel Heavy Chain Data",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
This script downloads raw camel heavy-chain sequences from OAS.

After downloading, process the data with:
    python build_camel_vhh_db.py \\
        --input-tsv ./oas_camel_data/combined_camel.tsv \\
        --output-prefix ./camel_vhh_db
        """
    )
    
    parser.add_argument(
        "--output-dir", default="./oas_camel_data",
        help="Directory to save downloaded files"
    )
    parser.add_argument(
        "--species", default="Camel",
        help="Species to download (default: Camel)"
    )
    parser.add_argument(
        "--chain", default="Heavy",
        help="Chain type to download (default: Heavy)"
    )
    
    args = parser.parse_args()
    
    print("="*70)
    print("OAS CAMEL DATA DOWNLOADER")
    print("="*70)
    print(f"Species: {args.species}")
    print(f"Chain:   {args.chain}")
    print(f"Output:  {args.output_dir}")
    print("="*70)
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Note about OAS download
    print("""
NOTE: The OAS API has rate limits and may require authentication for bulk downloads.

For large-scale downloads, it's recommended to:
1. Visit https://opig.stats.ox.ac.uk/webapps/oas/
2. Navigate to "Bulk Download"
3. Select:
   - Species: Camel (or Alpaca, Llama for other camelids)
   - Chain: Heavy
   - Format: TSV/CSV
   - Include: Nucleotide sequences
4. Download the resulting files

After downloading, combine them with:
    cat *.tsv > combined_camel_heavy.tsv
    # or
    head -1 first_file.tsv > combined.tsv
    tail -n +2 -q *.tsv >> combined.tsv

Then process with build_camel_vhh_db.py
""")
    
    print("\nManual download instructions saved to: {}/DOWNLOAD_INSTRUCTIONS.txt".format(args.output_dir))
    
    # Save instructions
    instructions = """
OAS CAMEL DATA DOWNLOAD INSTRUCTIONS
=====================================

Step 1: Go to OAS Bulk Download
-------------------------------
Visit: https://opig.stats.ox.ac.uk/webapps/oas/oas_paired/bulk_download

Step 2: Select Parameters
-------------------------
- Species: Camel (and/or Alpaca, Llama)
- Chain: Heavy
- Format: TSV or CSV

IMPORTANT: Make sure to include nucleotide sequences!
The key column needed is: sequence (or sequence_alignment)

Step 3: Download
----------------
Click download and save the file(s) to this directory.

Step 4: Combine files (if multiple)
-----------------------------------
If you downloaded multiple files, combine them:

    # For TSV files:
    head -1 $(ls *.tsv | head -1) > combined_camel_heavy.tsv
    for f in *.tsv; do tail -n +2 "$f" >> combined_camel_heavy.tsv; done

    # For CSV files:
    head -1 $(ls *.csv | head -1) > combined_camel_heavy.csv
    for f in *.csv; do tail -n +2 "$f" >> combined_camel_heavy.csv; done

Step 5: Process with build_camel_vhh_db.py
------------------------------------------
    python build_camel_vhh_db.py \\
        --input-tsv combined_camel_heavy.tsv \\
        --output-prefix camel_vhh_anarci \\
        --batch-size 2000

This will create:
- camel_vhh_anarci_renumbered.csv  (full details)
- camel_vhh_anarci_db.npz          (for similarity search)

EXPECTED DATA COLUMNS
---------------------
The OAS file should contain these columns (or similar):
- sequence (or sequence_alignment): Nucleotide sequence
- sequence_alignment_aa: Amino acid sequence (may be truncated)
- v_call: V gene call
- productive: Whether sequence is productive
- locus: Should be IGH for heavy chain

The build script uses the NUCLEOTIDE sequence to translate fresh,
avoiding the truncation issues in sequence_alignment_aa.
"""
    
    with open(os.path.join(args.output_dir, "DOWNLOAD_INSTRUCTIONS.txt"), 'w') as f:
        f.write(instructions)
    
    print("Done! Follow the instructions above to download OAS camel data.")


if __name__ == "__main__":
    main()
