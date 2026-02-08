#!/usr/bin/env python3
"""Quick diagnostic for NPZ structure"""
import os
import glob
import numpy as np

db_root = "VHH_shards"

# Find NPZ files
npz_files = sorted(glob.glob(os.path.join(db_root, "**", "*.npz"), recursive=True))
print(f"Found {len(npz_files)} NPZ files")

for npz_file in npz_files[:3]:
    print(f"\n=== {npz_file} ===")
    arr = np.load(npz_file, allow_pickle=True)
    print(f"Keys: {list(arr.keys())}")
    
    for key in arr.keys():
        data = arr[key]
        print(f"  {key}: shape={getattr(data, 'shape', 'N/A')}, dtype={getattr(data, 'dtype', type(data))}")
        
        if hasattr(data, '__len__'):
            print(f"    len={len(data)}")
            if len(data) > 0:
                print(f"    first item type: {type(data[0])}")
                if hasattr(data[0], 'shape'):
                    print(f"    first item shape: {data[0].shape}")
    
    arr.close()
