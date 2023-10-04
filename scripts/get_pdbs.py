"""This script gets the pdb files for each sequence in the balibase references.

__author__ = "Ben Iovino"
__date__ = "09/28/23"

"""

import os
import re
import subprocess
from Bio import SeqIO


def get_chain(msa: str, seq: str):
    """Selects chain from pdb file if one exists

    :param msa: msa name
    :param seq: sequence name
    """

    chain = re.search(r'_([A-Z])', seq)  # Chain is first letter right after '_'
    if chain is not None:
        chain = chain.group(1)
        os.system(
            f'pdb_selchain -{chain} data/pdb/{msa}/{seq}_.pdb'
            f'> data/pdb/{msa}/{seq}_{chain}.pdb')
        os.system(f'rm data/pdb/{msa}/{seq}_.pdb')


def get_matches(pdb_seq: str, fa_seq: str) -> list:
    """Returns a list of matching positions between two sequences

    :param pdb_seq: sequence from pdb file
    :param fa_seq: sequence from fasta file
    :return list: list of matching positions
    """

    # Find fasta sequence in pdb sequence
    beg = pdb_seq.find(fa_seq)
    while beg == -1:  # Remove characters from end until it matches
        fa_seq = fa_seq[:-1]
        beg = pdb_seq.find(fa_seq)
    print(beg)


def get_res(msa: str, seq: str):
    """Selects residues shared between sequence fasta and pdb

    :param msa: msa name
    :param seq: sequence name
    """

    # Get sequence from pdb file
    pdb_id = seq.split('.')[0]
    pdb_seq = subprocess.getoutput(
        f'pdb_tofasta data/pdb/{msa}/{pdb_id}.pdb')
    pdb_seq = ''.join(pdb_seq.split()[1:])  # Remove id from sequence

    # Get sequence from fasta file
    for record in SeqIO.parse(f'data/sequences/RV11/{msa}/{seq}', 'fasta'):
        fa_seq = str(record.seq)

    # Get all matching residues between pdb and fasta
    get_matches(pdb_seq, fa_seq)


def get_pdbs(seq_dir: str):
    """Writes all fasta files in a directory to pdb format

    :param seq_dir: directory of sequences
    """

    for msa in os.listdir(seq_dir):
        if msa in ['BB11006', 'BB11037']:
            continue
        if not os.path.exists(f'data/pdb/{msa}'):
            os.mkdir(f'data/pdb/{msa}')

        # For each seq, get pdb file (ignore chain)
        for seq in os.listdir(f'{seq_dir}/{msa}'):
            seq_ac = seq.split('_')[0]  # accession code
            #if seq_ac != '1lln':
                #continue

            '''
            os.system(
                f'wget https://files.rcsb.org/download/{seq_ac}.pdb '
                    f'-O data/pdb/{msa}/{seq_ac}_.pdb')

            # with pdb_tools, select chain and residues
            get_chain(msa, seq)
            '''

            get_res(msa, seq)


def main():
    """Initialize sequence directory
    """

    seq_dir = 'data/sequences/RV11'
    if not os.path.exists('data/pdb'):
        os.mkdir('data/pdb')
    get_pdbs(seq_dir)


if __name__ == '__main__':
    main()
