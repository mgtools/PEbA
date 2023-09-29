"""This script gets the pdb files for each sequence in the balibase references.

__author__ = "Ben Iovino"
__date__ = "09/28/23"

"""

import os


def main():


    seq_dir = 'data/sequences/RV11'
    os.mkdir('data/pdb')
    for msa in os.listdir(seq_dir):
        os.mkdir(f'data/pdb/{msa}')

        # For each seq, get pdb file (ignore chain)
        for seq in os.listdir(f'{seq_dir}/{msa}'):
            seq = seq.split('_')[0]
            os.system(
                f'wget https://files.rcsb.org/download/{seq}.pdb '
                    f'-O data/pdb/{msa}/{seq}.pdb')


if __name__ == '__main__':
    main()
