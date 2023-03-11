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
formatter = logging.Formatter('%(message)s')
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


def get_pairs(align):
    """=============================================================================================
    This function accepts an alignment and returns a dictionary of the pairs that are matched.

    :param aligns: dict of an alignment
    return dict: pairs
    ============================================================================================="""

    seq1_count = 0
    seq2_count = 0
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


def compute_score(aligns, beg, end):
    """=============================================================================================
    This function accepts two alignments and returns the PSC between them.

    :param aligns: dict containing two alignments
    :param beg: beginning of alignment
    :param end: end of alignment
    return float: PSC between the two alignments
    ============================================================================================="""

    # Get alignments from files and get their respective pairs
    align1 = aligns[list(aligns.keys())[0]]
    align2 = aligns[list(aligns.keys())[1]]
    align1_pairs = get_pairs(align1)
    align2_pairs = get_pairs(align2)

    # If no end is specified, use the length of the first alignment
    if end == 0:
        end = len(align1[0])

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

    # PSC is shared_pairs / length of comparison
    psc = round(shared_pairs / len(align1_pairs), 4)
    sim = round(seq_sim / len(align1_pairs), 4)
    logger.info('PCS: %s   ref_length: %s   comparison_length: %s   comparison_region: %s-%s   similarity: %s',
                psc, len(align1[0]), comparison_length, beg, end, sim)

    return psc


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-align1', type=str, default='ref0.msf', help='First alignment')
    parser.add_argument('-align2', type=str, default='peba0.msf', help='Second alignment')
    parser.add_argument('-id1', type=str, default='Q00X80', help='Sequence ID for first sequence')
    parser.add_argument('-id2', type=str, default='Q8L936', help='Sequence ID for second sequence')
    parser.add_argument('-beg', type=int, default=0, help='Beginning of alignment')
    parser.add_argument('-end', type=int, default=0, help='End of alignment')
    args = parser.parse_args()

    # Parse aligns and store in a dict - {align1: (seq1, seq2), align2: (seq1, seq2)}
    aligns = {}
    aligns[args.align1] = parse_align(args.align1, args.id1, args.id2)
    aligns[args.align2] = parse_align(args.align2, args.id1, args.id2)

    # Compute PSC
    compute_score(aligns, args.beg, args.end)


if __name__ == '__main__':
    main()
