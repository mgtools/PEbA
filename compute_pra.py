"""================================================================================================
This script compares two different alignments in msf format and calculates the 'percentage residues
aligned' (PRA) between them. The PRA is the number of residues that are aligned to the same position
in both alignments divided by the total number of residues in the first alignment.

Ben Iovino  03/9/23  VecAligns
================================================================================================"""

import argparse
import logging
import sys

logger = logging.getLogger('compute_pra')
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
    return: seq1, seq2 - aligned sequences; beg, end - beginning and end chars of aligns
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

            # Line 5 contains align info, get beginning and end chars of alig (if specified)
            if i == 4:
                try:
                    beg = line.split()[3].split('-')[0]
                    end = line.split()[3].split('-')[1]
                except IndexError:
                    beg = 0
                    end = 0
            if line.startswith('//'):  # Alignment starts after '//'
                in_align = True

    # Join the two sequences and return them
    seq1 = ''.join(seq1)
    seq2 = ''.join(seq2)

    return seq1, seq2, beg, end


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


def compute_score(aligns):
    """=============================================================================================
    This function accepts two alignments and returns the PRA between them.

    :param aligns: dict containing two alignments
    return float: PRA between the two alignments
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
    beg1 = int(align1[0].replace('.', '').find(align2[0].replace('.', '')))
    beg2 = int(align1[1].replace('.', '').find(align2[1].replace('.', '')))
    if beg1 == -1:
        beg1 = 0
    if beg2 == -1:
        beg2 = 0
    align2_pairs = get_pairs(align2, beg1, beg2)

    # Beginning character of comparison is the max beginning position
    beg = max(int(align1[2]), int(align2[2]))

    # If end character is 0 for both, set it to the length of ref alignment
    if (int(align1[3]) and int(align2[3])) == 0:
        end = len(align1[0])
    if int(align1[3]) == 0:
        if int(align2[3]) > 0:
            end = int(align2[3])
    if int(align2[3]) == 0:
        if int(align1[3]) > 0:
            end = int(align1[3])

    # Compare pairs between the two alignments
    shared_pairs = 0
    seq_sim = 0
    comparison_length = 0
    for key, value in align1_pairs.items():
        if int(key[1:]) < beg or int(key[1:]) > end:  # Check if pair is in the comparison region
            continue
        comparison_length += 1
        if key[0] == value[0]:  # First index is the character
            seq_sim += 1
        if key in align2_pairs:  # Check if same pair is in the other alignment
            if align2_pairs[key] == value:
                shared_pairs += 1

    # PRA is shared_pairs / length of comparison
    pra = round(shared_pairs / len(align1_pairs)*100, 2)
    sim = round(seq_sim / len(align1_pairs)*100, 2)
    logger.info('PRA: %s   ref_length: %s   comparison_length: %s   similarity: %s',
                pra, len(align1[0]), comparison_length, sim)

    return pra


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-align1', type=str, default='ref1.msf', help='First alignment')
    parser.add_argument('-align2', type=str, default='peba1.msf', help='Second alignment')
    args = parser.parse_args()

    # Parse aligns and store in a dict - {align: (seq1, seq2, beg, end)}
    aligns = {}
    aligns[args.align1] = parse_align(args.align1)
    aligns[args.align2] = parse_align(args.align2)

    # Compute PRA
    compute_score(aligns)


if __name__ == '__main__':
    main()
