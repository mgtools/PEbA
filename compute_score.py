"""================================================================================================
This script compares two different alignments in msf format and calculates the 'total column score'
(TCS) between them. The TCS is the number of residues that are aligned to the same position
in both alignments divided by the total number of residues in the reference (first) alignment.

Ben Iovino  03/9/23  VecAligns
================================================================================================"""

import argparse
import logging
import sys

logger = logging.getLogger('compute_tcs')
logger.setLevel(logging.INFO)
ch = logging.StreamHandler(stream=sys.stdout)
formatter = logging.Formatter('%(message)s')
ch.setFormatter(formatter)
logger.addHandler(ch)


def parse_align(filename):
    """=============================================================================================
    This function accepts an alignment in msf format and returns each sequence with gaps as well
    as the character positions that the alignment begins and ends at.

    :param filename: name of file
    return seq1, seq2, id1, id2: sequences and names of sequences
    ============================================================================================="""

    seq1 = []
    seq2 = []
    in_align = False
    count = 0
    with open(filename, 'r', encoding='utf8') as file:
        lines = file.readlines()
        for i, line in enumerate(lines):

            # We are "in the alignment" after the '//', read each aligned sequence
            if in_align is True:
                if line.startswith('\n'):
                    continue
                if count == 0 or count % 2 == 0:  # Even count indicates first sequence
                    seq1.append(''.join(line.split()[1:]))
                else:  # Odd count indicates second sequence
                    seq2.append(''.join(line.split()[1:]))
                count += 1
                continue

            # Get names of sequences
            if i == 6:
                id1 = line.split()[1]
            if i == 7:
                id2 = line.split()[1]

            # Alignment starts after '//'
            if line.startswith('//'):
                in_align = True

    # Join the two sequences and return them
    seq1 = ''.join(seq1)
    seq2 = ''.join(seq2)

    return seq1, seq2, id1, id2


def get_pairs(align, beg1, beg2):
    """=============================================================================================
    This function accepts an alignment and returns a dictionary of the pairs that are matched.

    :param aligns: dict of an alignment
    :param beg1: beginning char of first sequence
    :param beg2: beginning char of second sequence
    return dict: pairs
    ============================================================================================="""

    seq1_count = beg1
    seq2_count = beg2
    pairs = {}

    # Iterate through each column in the alignment
    for i in range(len(align[0])):
        seq1_char = align[0][i]
        seq2_char = align[1][i]

        # If there is a gap, increment the count of the other sequence
        if seq1_char == '.':
            seq2_count += 1
            continue
        if seq2_char == '.':
            seq1_count += 1
            continue

        # If there is a match, add the pair to the dict
        pairs[f'{seq1_char}{seq1_count}'] = f'{seq2_char}{seq2_count}'
        seq1_count += 1
        seq2_count += 1

    return pairs


def compute_score(aligns, score):
    """=============================================================================================
    This function accepts two alignments and returns either the TCS or F1 score between them. The
    score and other information is logged to stdout.

    :param aligns: dict containing two alignments
    :param score: tcs or f1
    ============================================================================================="""

    # Get alignments from files and get their respective pairs
    align1 = aligns[list(aligns.keys())[0]]
    align2 = aligns[list(aligns.keys())[1]]

    # Find the beginning characters of ref align and find pairs
    beg1 = int(align2[0].replace('.', '').find(align1[0].replace('.', '')))
    beg2 = int(align2[1].replace('.', '').find(align1[1].replace('.', '')))
    if beg1 == -1:
        beg1 = 0
    if beg2 == -1:
        beg2 = 0
    align1_pairs = get_pairs(align1, beg1, beg2)

    # Find the beginning characters of test align and find pairs
    beg3 = int(align1[0].replace('.', '').find(align2[0].replace('.', '')))
    beg4 = int(align1[1].replace('.', '').find(align2[1].replace('.', '')))
    if beg3 == -1:
        beg3 = 0
    if beg4 == -1:
        beg4 = 0
    align2_pairs = get_pairs(align2, beg3, beg4)

    # Beginning of comparison is the max beg position of the two aligns
    beg = max(beg1, beg3)

    # End of comparison is the min end position of the two aligns
    end1 = list(align1_pairs.keys())[-1]
    end2 = list(align2_pairs.keys())[-1]
    end = min(int(end1[1:]), int(end2[1:]))

    # Compare pairs between the two alignments
    shared_pairs = 0
    seq_sim = 0
    for key, value in align1_pairs.items():
        if int(key[1:]) < beg or int(key[1:]) > end:  # Check if pair is in the comparison region
            continue
        if key[0] == value[0]:  # First index is the character
            seq_sim += 1
        if key in align2_pairs:  # Check if same pair is in the other alignment
            if align2_pairs[key] == value:
                shared_pairs += 1

    # TCS is shared_pairs / number of pairs in first align
    if score == 'tcs':
        sc = round(shared_pairs / len(align1_pairs), 3)

    # F1 is 2 * (precision*recall) / (precision + recall)
    elif score == 'f1':
        precision = shared_pairs / (shared_pairs + len(align2_pairs) - shared_pairs)
        recall = shared_pairs / (shared_pairs + len(align1_pairs) - shared_pairs)

        # Check if F1 is zero, can't divide by zero
        if precision == 0 and recall == 0:
            sc = 0
        else:
            sc = 2 * (precision * recall) / (precision + recall)

    # Report the score and other info
    sim = round(seq_sim / len(align1_pairs)*100, 2)
    logger.info('%s: %s   ref_length: %s   comparison_length: %s   similarity: %s',
            score.upper(), sc, len(align1[0]), len(align1_pairs), sim)


def main():
    """=============================================================================================
    This function takes two alignments, the first being the reference, and the second being the
    alignment to be compared to the reference. It calls parse_align() to extract the sequences
    and relevant info from the alignments and then calls compute_score() to get the TCS.
    ============================================================================================="""

    parser = argparse.ArgumentParser()
    parser.add_argument('-align1', type=str, help='First alignment')
    parser.add_argument('-align2', type=str, help='Second alignment')
    parser.add_argument('-score', type=str, default='tcs', help='Comparison score (tcs or f1)')
    args = parser.parse_args()

    # Parse aligns and store in a dict - {align: (seq1, seq2, beg, end)}
    aligns = {}
    aligns[args.align1] = parse_align(args.align1)
    aligns[args.align2] = parse_align(args.align2)

    # Compute score
    compute_score(aligns, args.score)


if __name__ == '__main__':
    main()
