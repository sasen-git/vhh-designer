#!/usr/bin/env python3
"""Debug script to check NPZ file format."""

import numpy as np
import sys
import os

def inspect_npz(path):
    print(f"\n{'='*60}")
    print(f"FILE: {path}")
    print('='*60)
    
    try:
        data = np.load(path, allow_pickle=True)
        
        print(f"\nKeys: {list(data.keys())}")
        
        for key in data.keys():
            arr = data[key]
            print(f"\n  {key}:")
            print(f"    type: {type(arr)}")
            print(f"    dtype: {arr.dtype}")
            print(f"    shape: {arr.shape}")
            
            # Show sample
            if arr.size > 0:
                if arr.dtype == object:
                    print(f"    first item type: {type(arr.flat[0])}")
                    sample = arr.flat[0]
                    if isinstance(sample, (bytes, str)):
                        print(f"    first item: {str(sample)[:100]}...")
                    elif isinstance(sample, np.ndarray):
                        print(f"    first item shape: {sample.shape}")
                        print(f"    first item dtype: {sample.dtype}")
                        print(f"    first item sample: {sample[:10]}...")
                    elif isinstance(sample, dict):
                        print(f"    first item keys: {list(sample.keys())[:10]}")
                    else:
                        print(f"    first item: {sample}")
                else:
                    print(f"    sample: {arr.flat[:5]}...")
                    
    except Exception as e:
        print(f"  ERROR: {e}")

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python debug_npz_format.py <npz_file_or_directory>")
        sys.exit(1)
    
    path = sys.argv[1]
    
    if os.path.isdir(path):
        # Process first 2 files
        files = [f for f in os.listdir(path) if f.endswith('.npz')][:2]
        for f in files:
            inspect_npz(os.path.join(path, f))
    else:
        inspect_npz(path)
