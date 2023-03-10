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


def match_pairs(align):
    """=============================================================================================
    This function accepts an alignment and returns a dictionary of the pairs that are matched.

    :param aligns: dict of an alignment
    return dict: pairs
    ============================================================================================="""

    seq1_count = 0
    seq2_count = 0
    pairs = {}
    for i in range(len(align[0])):
        seq1_char = align[0][i]
        seq2_char = align[1][i]
        if seq1_char == '.':
            seq2_count += 1
            continue
        if seq2_char == '.':
            seq1_count += 1
            continue
        pairs[f'{seq1_char}{seq1_count}'] = f'{seq2_char}{seq2_count}'
        seq1_count += 1
        seq2_count += 1
    print(pairs)
    return pairs


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

    count1 = 0
    count2 = 0
    for i in range(beg, end):
        ref_pair = (align1[0][i], align1[1][i])
        test_pair = (align2[0][i], align2[1][i])
        print(ref_pair[0], count1, test_pair[0], count2)
        if ref_pair[0] == '.' and test_pair[0] == '.':
            count1 += 1
            count2 += 1
            continue
        if '.' in ref_pair[0] or '.' in test_pair[0]:
            if '.' in ref_pair:
                count1 += 1
            if '.' in test_pair:
                count2 += 1
        elif ref_pair[0] != '.' and test_pair[0] != '.':
            count1 += 1
            count2 += 1
    print(count1)


    # PSC is shared_cols / length of comparison
    #psc = round(shared_cols / length, 4)
    #sim = round(similarity / length, 4)
    #logger.info('PCS: %s   ref_length: %s   comparison_length: %s   comparison_region: %s-%s   similarity: %s',
                #psc, len(align1[0]), length, beg, end, sim)
    #return psc


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-align1', default='ref0.msf', help='First alignment')
    parser.add_argument('-align2', default='peba0.msf', help='Second alignment')
    parser.add_argument('-id1', default='Q00X80', help='Sequence ID for first sequence')
    parser.add_argument('-id2', default='Q8L936', help='Sequence ID for second sequence')
    parser.add_argument('-beg', default=0, help='Beginning of alignment')
    parser.add_argument('-end', default=0, help='End of alignment')
    args = parser.parse_args()

    # Parse aligns and store in a dict - {align1: (seq1, seq2), align2: (seq1, seq2)}
    aligns = {}
    aligns[args.align1] = parse_align(args.align1, args.id1, args.id2)
    aligns[args.align2] = parse_align(args.align2, args.id1, args.id2)

    # Compute PSC
    compute_score(aligns, args.beg, args.end)


if __name__ == '__main__':
    main()


"""=============================================================================================
    # Compare columns in specified range
    shared_cols = 0
    length = 0
    similarity = 0
    for i in range(beg, end):

        # Define pairs of residues in each alignment
        pair1 = (align1[0][i], align1[1][i])
        pair2 = (align2[0][i], align2[1][i])
        if '.' in pair1:  # Skip columns with gaps in reference alignment
            continue
        if pair1[0] == pair1[1]:  # Increment similarity if residues are identical
            similarity += 1
        if pair1 == pair2:
            shared_cols += 1
        length += 1  # Increment length of comparison
    ============================================================================================="""
