"""================================================================================================
This script compares two different alignments and calculates the 'percentage shared columns' (PSC)
between them. The PSC is the number of columns that are identical in both alignments divided by the
the total number of columns in the first alignment.

Ben Iovino  03/9/23  VecAligns
================================================================================================"""

import argparse
import logging

logger = logging.getLogger('compute_psc')
logger.setLevel(logging.INFO)
ch = logging.StreamHandler()
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
ch.setFormatter(formatter)
logger.addHandler(ch)


def parse_align(filename, id1, id2):
    """=============================================================================================
    This function accepts an alignment file and two sequence ids and returns the PCS between them.

    :param filename: name of file
    :param id1: first sequence id
    :param id2: second sequence id
    return string: seq1, seq2 - each sequence in the alignment
    ============================================================================================="""

    seq1 = []
    seq2 = []
    with open(filename, 'r', encoding='utf8') as file:
        for line in file:  # Append sequences to corresponding lists
            if line.startswith(id1):
                seq1.append(''.join(line.split()[1:]))
            elif line.startswith(id2):
                seq2.append(''.join(line.split()[1:]))

    # Join the two sequences and return them
    seq1 = ''.join(seq1)
    seq2 = ''.join(seq2)
    return seq1, seq2


def compute_score(aligns, beg, end):
    """=============================================================================================
    This function accepts two alignments and returns the PSC between them.

    :param aligns: dict containing two alignments
    :param beg: beginning of alignment
    :param end: end of alignment
    return float: PSC between the two alignments
    ============================================================================================="""

    align1 = aligns[list(aligns.keys())[0]]
    align2 = aligns[list(aligns.keys())[1]]

    # End of comparison is either specified or the length of the shortest sequence
    if end == 0:
        end = min(len(align1[0]), len(align2[0]))

    # Compare columns in specified range
    shared_cols = 0
    for i in range(beg, end):
        pair1 = (align1[0][i], align1[1][i])
        pair2 = (align2[0][i], align2[1][i])
        if pair1 == pair2:
            shared_cols += 1

    # PSC is shared_cols / length of comparison
    return shared_cols / (end-beg)


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-align1', default='peba0.msf', help='First alignment')
    parser.add_argument('-align2', default='ref0.msf', help='Second alignment')
    parser.add_argument('-id1', default='1aab_', help='Sequence ID for first sequence')
    parser.add_argument('-id2', default='1j46_A', help='Sequence ID for second sequence')
    parser.add_argument('-beg', default=0, help='Beginning of alignment')
    parser.add_argument('-end', default=0, help='End of alignment')
    args = parser.parse_args()

    # Parse aligns and store in a dict - {align1: (seq1, seq2), align2: (seq1, seq2)}
    aligns = {}
    aligns[args.align1] = parse_align(args.align1, args.id1, args.id2)
    aligns[args.align2] = parse_align(args.align2, args.id1, args.id2)

    # Compute PSC
    psc = compute_score(aligns, args.beg, args.end)
    print(psc)


if __name__ == '__main__':
    main()
