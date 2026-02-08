#!/usr/bin/env python
# JD

import sys
try:
    import anarci

except:
    print("I can't find anarci! Make sure you're in the right python environment - for this script, that is `biopy` or `pmpnn`. Try `conda env list` to make sure you're in the right environment, `conda list | grep anarci` to make ensure your current environment has anarci installed, and `whereis python` to ensure you only have 1 active python environment as initial trouble shooting steps.")
    sys.exit()
import random
import argparse
import os


# I tried to be really fancy and point the script to the appropriate environment without requiring being in the correct environment, but idk how to call an executable from that environment with dependencies in that environment.
#def anarci_out(scfv):
#    ap = '/home/tagteam/anaconda3/envs/biopy/bin/ANARCI'
#    command = f'{ap} -h'
#    #os.system('ANARCI -h')
#    result = subprocess.run(command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#    return result.stdout, result.stderr

def anarci_out(antibody_fasta):
    carrot_counter = 0
    with open(antibody_fasta) as file:
        for line in file:
            if line[0] == '>':
                carrot_counter += 1

    if carrot_counter > 2:
        print(f'ERROR: {antibody_fasta} has more than 2 sequences. Please supply 1 sequence that is an scfv to see if we can optimize it, or provide 2 sequences to combine them into an scfv. I detected {carrot_counter} sequences.')
        sys.exit()
    N = random.randint(0,10000)
    anarci_file = f'tmp{N}.tmp'
    os.system(f'ANARCI -i {antibody_fasta} -s kabat > {anarci_file}')
    return anarci_file


def parse_sequences_from_file(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    heavy_chain_sequence = ""
    light_chain_sequence = ""

    # Extracting the heavy and light chain sequences based on the first letter of each line
    for line in lines:
        if line.startswith("H "):
            parts = line.strip().split()
            heavy_chain_sequence += parts[-1]
        elif line.startswith("L "):
            parts = line.strip().split()
            light_chain_sequence += parts[-1]

    return heavy_chain_sequence, light_chain_sequence


def concatenate_sequences(file_path, linker="GGGGSGGGGSGGGGS", polyG=False):
    heavy_chain_sequence, light_chain_sequence = parse_sequences_from_file(file_path)
    out_seq = f"{heavy_chain_sequence}{linker}{light_chain_sequence}"
    out_seq = out_seq.replace('-','')
    if out_seq[0] != 'M':
        out_seq = 'M' + out_seq
    
    if polyG:
        out_seq += ':GGGGGGGGGG'
    return out_seq

# Example usage
def main():
    parser = argparse.ArgumentParser(description="turn an antibody sequence into an scFv via ANARCI.")
    parser.add_argument("fasta_file", type=str,  help="Path to the fasta file.")
    parser.add_argument("--linker-seq", type=str,  default="GGGGSGGGGSGGGGS", help="linker sequence used to combine the VH and VL chains.")
    parser.add_argument("--output-file", type=str, default=None, help="if specified, create fasta file of new scFv")
    parser.add_argument("--polyG", action='store_true', help="True/False: if included, add (polyG)x10 to fasta sequence. Intended for colabfold structure prediction.")

    args = parser.parse_args()
   
    anarci_file = anarci_out(args.fasta_file) 
    concatenated_sequence = concatenate_sequences(anarci_file, args.linker_seq, args.polyG)
    os.system(f'rm {anarci_file}')

    if not args.output_file:
        print(concatenated_sequence)

    else:
        if args.output_file[-6:] != '.fasta':
            file_name  = args.output_file + '.fasta'
            fasta_name = args.output_file
        else:
            file_name  = args.output_file
            fasta_name = args.output_file[:-6]
        
        with open(file_name, 'w') as file:
            file.write(f'>{fasta_name}\n')
            file.write(concatenated_sequence)
            
if __name__ == "__main__":
    main()
