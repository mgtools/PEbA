"""Calls fatcat on two pdb files and writes the alignment in the desired format

__author__ = "Ben Iovino"
__date__ = 10/02/23
"""

import argparse
import subprocess
import utility as ut


def clean_align(align: str) -> tuple:
    """Returns aligned sequences from fatcat output

    :param align: fatcat stdout output
    :return (str, str, list, list): aligned sequences, beg/end positions of each seq
    """

    align1, align2 = '', ''
    al1_pos, al2_pos = [], []
    for line in align.split('\n'):
        if line.startswith('Chain 1:'):
            line1 = line.split()
            align1 += line1[3].replace('-', '.')
            al1_pos.append(line1[2])
        elif line.startswith('Chain 2:'):
            line2 = line.split()
            align2 += line2[3].replace('-', '.')
            al2_pos.append(line2[2])

    # Beg/end positions of each seq
    beg = [al1_pos[0], al2_pos[0]]
    end1 = int(line1[2]) + len(line1[3].strip('-'))-1
    end2 = int(line2[2]) + len(line2[3].strip('-'))-1
    end = [str(end1), str(end2)]

    return align1, align2, beg, end


def fatcat(pdb1: str, pdb2: str) -> tuple:
    """Returns FATCAT alignment of two pdb files

    :param pdb1: first pdb file
    :param pdb2: second pdb file
    :return (str, str, list, list): aligned sequences, beg/end positions of each seq
    """

    align = subprocess.getoutput(f'FATCAT -p1 {pdb1} -p2 {pdb2} -q')

    return clean_align(align)


def main():
    """Takes two pdb files and writes the alignment in the desired format
    """

    parser = argparse.ArgumentParser()
    parser.add_argument('-f1', type=str, default='data/pdb/BB11001/1aab_.pdb')
    parser.add_argument('-f2', type=str, default='data/pdb/BB11001/1j46_A.pdb')
    parser.add_argument('-o', '--output', type=str, default='msf')
    parser.add_argument('-sf', '--savefile', type=str, default='/home/ben/Desktop')
    args = parser.parse_args()

    align1, align2, beg, end = fatcat(args.f1, args.f2)
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
