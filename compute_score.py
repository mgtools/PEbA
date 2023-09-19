"""Compares two pairwise alignments and returns the Sum of Pairs (SP), Total Column Score (TCS),
or F1 score between them.

__author__ = "Ben Iovino"
__date__ = 09/18/23
"""

import argparse
import logging
import sys

logger = logging.getLogger('compute_tcs')
logger.setLevel(logging.INFO)
ch = logging.StreamHandler(stream=sys.stdout)
formatter = logging.Formatter('%(message)s')
ch.setFormatter(formatter)
logger.addHandler(ch)


def parse_align(filename: str) -> tuple:
    """Returns a tuple with the two aligned sequences and their names.

    :param filename: name of file (msf format)
    return (str, str): sequences with gaps
    """

    seq1, seq2 = [], []
    in_align = False
    count = 0
    with open(filename, 'r', encoding='utf8') as file:
        lines = file.readlines()
    for line in lines:

        # We are "in the alignment" after the '//', read each aligned sequence
        if in_align is True:
            if line.startswith('\n'):
                continue
            if count % 2 == 0:  # Even count indicates first sequence
                seq1.append(''.join(line.split()[1:]))
            else:  # Odd count indicates second sequence
                seq2.append(''.join(line.split()[1:]))
            count += 1
            continue

        # Alignment starts after '//'
        if line.startswith('//'):
            in_align = True

    # Join the two sequences and return them
    seq1 = ''.join(seq1)
    seq2 = ''.join(seq2)

    return seq1, seq2


def get_pairs(filename: str) -> dict:
    """Returns dictionary where keys are the positions of the aligned characters in the alignment
    and values are the matched positions (including gaps).

    :param filename: name of file (msf format)
    return dict: dict where keys are positions of aligned characters and values are matched positions
    """

    seq1, seq2 = parse_align(filename)
    pairs = {}
    seq1_count, seq2_count = 0, 0
    for i, char1 in enumerate(seq1):
        char2 = seq2[i]  # Character in second sequence at same position in alignment
        if char1 != '.':  # If not a gap, add the position to the count
            char1 = f'{char1}{seq1_count}'
            seq1_count += 1
        if seq2[i] != '.':
            char2 = f'{seq2[i]}{seq2_count}'
            seq2_count += 1
        pairs[i] = (char1, char2)

    return pairs


def sp_score(al1: dict, al2: dict) -> tuple:
    """Returns sum of pairs (sp) score between two alignments.

    :param al1: dict where keys are positions of aligned characters and values are matched positions
    :param al2: dict same as al1
    :return (float, float): similarity and sp scores
    """

    score, total, sim = 0, 0, 0
    for pair1 in al1.values():
        if '.' in pair1:  # if gap in first align, skip comparison
            continue
        total += 1  # if no gap, add to total number of columns compared
        if pair1[0][0] == pair1[1][0]:  # if same character, add to similarity
            sim += 1
        if pair1 in al2.values():
            score += 1

    # sp is (shared pairs between ref/test align) / (total number of pairs in ref align)
    score = round(score/total, 3)
    sim = round(sim/total, 3)
    logger.info('SP: %s   ref_length: %s   comparison_length: %s   similarity: %s',
            score, len(al1.values()), total, sim)


def tc_score(al1: dict, al2: dict) -> tuple:
    """Returns TC (total column) score between two alignments.

    :param al1: dict where keys are positions of aligned characters and values are matched positions
    :param al2: dict same as al1
    :return (float, float): similarity and TC scores
    """

    score, sim, total = 0, 0, 0
    for pair1 in al1.values():
        if '.' not in pair1:  # for calculating similarity
            total += 1
        if pair1 in al2.values():
            if pair1[0][0] == pair1[1][0]:  # if same character, add to similarity
                sim += 1
            score += 1

    # TC is (shared number of columns between ref/test align) / (total number of cols in ref align)
    score = round(score/len(al1), 3)
    sim = round(sim/total, 3)
    logger.info('TCS: %s   ref_length: %s   comparison_length: %s   similarity: %s',
            score, len(al1.values()), len(al1.values()), sim)


def f1_score(al1: dict, al2: dict) -> tuple:
    """Returns F1 score between two alignments.

    :param al1: dict where keys are positions of aligned characters and values are matched positions
    :param al2: dict same as al1
    :return (float, float): 
    """

    shared_pairs, sim, total = 0, 0, 0
    for pair1 in al1.values():
        if '.' not in pair1:  # for calculating similarity
            total += 1
        if pair1 in al2.values():
            if pair1[0][0] == pair1[1][0]:  # if same character, add to similarity
                sim += 1
            shared_pairs += 1

    # prec/rec
    precision = shared_pairs / (shared_pairs + len(al2.values()) - shared_pairs)
    recall = shared_pairs / (shared_pairs + len(al1.values()) - shared_pairs)

    # Check if F1 is zero, can't divide by zero
    if precision == 0 and recall == 0:
        score = 0
    else:
        score = 2 * (precision * recall) / (precision + recall)
        score = round(score, 3)

    sim = round(sim/total, 3)
    logger.info('F1: %s   ref_length: %s   comparison_length: %s   similarity: %s',
            score, len(al1.values()), len(al1.values()), sim)


def compare_aligns(args: argparse.Namespace) -> float:
    """Returns score between two alignments.

    :param args.align1: first alignment
    :param args.align2: second alignment
    :param args.score: score to return (sp/tcs/f1)
    :return float: score
    """

    al1 = get_pairs(args.align1)
    al2 = get_pairs(args.align2)

    if args.score == 'sp':
        sp_score(al1, al2)
    if args.score == 'tcs':
        tc_score(al1, al2)
    if args.score == 'f1':
        f1_score(al1, al2)


def main():
    """Main
    """

    parser = argparse.ArgumentParser()
    parser.add_argument('-align1', type=str, help='First alignment', default='/home/ben/Desktop/BOX192_0.msf')
    parser.add_argument('-align2', type=str, help='Second alignment', default='/home/ben/Desktop/alignment_0.msf')
    parser.add_argument('-score', type=str, default='tcs', help='Comparison score (sp/tcs/f1)')
    args = parser.parse_args()

    compare_aligns(args)


if __name__ == '__main__':
    main()
