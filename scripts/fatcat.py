"""Calls fatcat on two pdb files and writes the alignment in the desired format

__author__ = "Ben Iovino"
__date__ = 10/02/23
"""

import argparse
import subprocess
import utility as ut
from Bio import SeqIO


def get_pos(align: str, pdb: str) -> tuple:
    """Returns beginning and end positions of an aligned sequence

    :param align: aligned sequence
    :param pdb: pdb file
    :return (int, int): beginning and end positions
    """

    # Get fasta sequence from file
    pdb = pdb.split('/')
    name, direc = pdb[-1].split('.')[0], pdb[-2]
    file = f'data/sequences/RV11/{direc}/{name}.fa'
    seq = SeqIO.read(file, 'fasta').seq

    # Remove gaps and find beginning and end positions
    align = align.replace('.', '')
    beg = seq.find(align)
    tmp_align = align
    while beg == -1:  # Some aligns don't match exactly
        tmp_align = tmp_align[:-1]
        beg = seq.find(tmp_align)
    end = beg + len(align)-1

    return beg, end


def clean_align(align: str) -> tuple:
    """Returns aligned sequences from fatcat output

    :param align: fatcat stdout output
    :return (str, str): aligned sequences with gaps
    """

    align1, align2 = '', ''
    for line in align.split('\n'):
        if line.startswith('Chain 1:'):
            line1 = line.split()
            align1 += line1[3].replace('-', '.')
        elif line.startswith('Chain 2:'):
            line2 = line.split()
            align2 += line2[3].replace('-', '.')

    return align1, align2


def fatcat(pdb1: str, pdb2: str) -> tuple:
    """Returns FATCAT alignment of two pdb files

    :param pdb1: first pdb file
    :param pdb2: second pdb file
    :return (str, str, list, list): aligned sequences, beg/end positions of each seq
    """

    align = subprocess.getoutput(f'FATCAT -p1 {pdb1} -p2 {pdb2} -q')
    align1, align2 = clean_align(align)
    beg1, end1 = get_pos(align1, pdb1)
    beg2, end2 = get_pos(align2, pdb2)
    beg = [beg1, beg2]
    end = [end1, end2]

    return align1, align2, beg, end


def main():
    """Takes two pdb files and writes the alignment in the desired format
    """

    parser = argparse.ArgumentParser()
    parser.add_argument('-f1', type=str)
    parser.add_argument('-f2', type=str)
    parser.add_argument('-o', '--output', type=str, default='msf')
    parser.add_argument('-sf', '--savefile', type=str)
    args = parser.parse_args()

    align1, align2, beg, end = fatcat(args.f1, args.f2)
    fatcat(args.f1, args.f2)
    id1 = args.f1.split('/')[-1].split('.')[0]
    id2 = args.f2.split('/')[-1].split('.')[0]

    # Write align based on desired output format
    if args.output == 'msf':
        ut.write_msf(align1, align2, id1, id2, 'fatcat',
                      0, 0, args.savefile, beg, end)
    if args.output == 'fa':
        ut.write_fasta(align1, align2, id1, id2, args.savefile)


if __name__ == '__main__':
    main()
