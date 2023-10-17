"""This script gets the pdb files for each sequence in the balibase references.

__author__ = "Ben Iovino"
__date__ = "09/28/23"

"""

import os
import re


def get_pdbs(seq_dir: str):
    """Writes all fasta files in a directory to pdb format

    :param seq_dir: directory of sequences
    """

    for msa in os.listdir(seq_dir):
        os.mkdir(f'data/pdb/{msa}')

        # For each seq, get pdb file (ignore chain)
        for seq in os.listdir(f'{seq_dir}/{msa}'):
            chain = re.search(r'_([A-Z])', seq)  # Chain is first letter right after '_'
            seq = seq.split('_')[0]
            os.system(
                f'wget https://files.rcsb.org/download/{seq}.pdb '
                    f'-O data/pdb/{msa}/{seq}_.pdb')

            # with pdb_tools select chain, if necessary
            if chain is not None:
                chain = chain.group(1)
                os.system(
                    f'pdb_selchain -{chain} data/pdb/{msa}/{seq}_.pdb'
                    f'> data/pdb/{msa}/{seq}_{chain}.pdb')
                os.system(f'rm data/pdb/{msa}/{seq}_.pdb')


def main():
    """Initialize sequence directory
    """

    seq_dir = 'data/sequences/RV11'
    os.mkdir('data/pdb')
    get_pdbs(seq_dir)


if __name__ == '__main__':
    main()
