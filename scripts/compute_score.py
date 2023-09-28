"""Compares two pairwise alignments and returns the Sum of Pairs (SP) or F1 score between them.

__author__ = "Ben Iovino"
__date__ = 09/18/23
"""

import argparse
import logging
import sys

logger = logging.getLogger('compute_score')
logger.setLevel(logging.INFO)
ch = logging.StreamHandler(stream=sys.stdout)
formatter = logging.Formatter('%(message)s')
ch.setFormatter(formatter)
logger.addHandler(ch)


def parse_align(filename: str) -> tuple:
    """Returns a tuple with the two aligned sequences and their names.

    :param filename: name of file (msf format)
    return (str, str, dict): sequences with gaps and dict with beg/end positions
    """

    seq1, seq2 = [], []
    in_align = False
    count = 0
    pos = {}
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

        # Get beginning/end positions of aligned sequences
        if line.startswith(' Name:'):
            line = line.split()
            pos[line[1]] = line[6]
            continue

        # Alignment starts after '//'
        if line.startswith('//'):
            in_align = True

    # Join the two sequences and return them
    seq1 = ''.join(seq1)
    seq2 = ''.join(seq2)

    return seq1, seq2, pos


def get_pos(positions: dict) -> tuple:
    """Returns beginning positions of aligned sequences

    :param positions: dict with beg/end positions
    return (int, int): beginning positions of aligned sequences
    """

    seq1_count = list(positions.values())[0]
    seq2_count = list(positions.values())[1]
    seq1_count = int(seq1_count.split(',')[0])
    seq2_count = int(seq2_count.split(',')[0])

    return seq1_count, seq2_count


def get_pairs(filename: str) -> dict:
    """Returns dictionary where keys are the positions of the aligned characters in the alignment
    and values are the matched positions (including gaps).

    :param filename: name of file (msf format)
    return dict: dict where keys are positions of aligned residues and values are matched positions
    """

    seq1, seq2, pos = parse_align(filename)
    pairs = {}
    seq1_count, seq2_count = get_pos(pos)
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


def get_shared(al1: dict, al2: dict) -> tuple:
    """Returns the number of shared pairs between two alignments, as well as the total number of
    pairs compared and the pairwise similarity in the first alignment.

    :param al1: dict where keys are positions of aligned residues and values are matched positions
    :param al2: dict same as al1
    :return (int, int, float): shared pairs, total number of pairs, and pairwise similarity
    """

    shared_pairs, total, sim = 0, 0, 0
    for pair1 in al1.values():
        if '.' in pair1:  # if gap in first align, skip comparison
            continue
        total += 1  # if no gap, add to total number of columns compared
        if pair1[0][0] == pair1[1][0]:  # if same character, add to similarity
            sim += 1
        if pair1 in al2.values():
            shared_pairs += 1

    return shared_pairs, total, sim


def sp_score(al1: dict, al2: dict) -> tuple:
    """Returns sum of pairs (sp) score between two alignments.

    :param al1: dict where keys are positions of aligned residues and values are matched positions
    :param al2: dict same as al1
    :return (float, float): similarity and sp scores
    """

    shared_pairs, total, sim = get_shared(al1, al2)

    # sp is (shared pairs between ref/test align) / (total number of pairs in ref align)
    score = round(shared_pairs/total, 3)
    sim = round(sim/total, 3)
    logger.info('SP: %s   ref_length: %s   comparison_length: %s   similarity: %s',
            score, len(al1.values()), total, sim)


def f1_score(al1: dict, al2: dict) -> tuple:
    """Returns F1 score between two alignments.

    :param al1: dict where keys are positions of aligned residues and values are matched positions
    :param al2: dict same as al1
    :return (float, float): 
    """

    shared_pairs, total, sim = get_shared(al1, al2)

    # Number of non-gapped pairs in each alignment
    al1_pairs = sum([1 for pair in al1.values() if '.' not in pair])
    al2_pairs = sum([1 for pair in al2.values() if '.' not in pair])

    # prec/rec
    precision = shared_pairs / (shared_pairs + al2_pairs - shared_pairs)
    recall = shared_pairs / (shared_pairs + al1_pairs - shared_pairs)

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
    :param args.score: score to return (sp/f1)
    :return float: score
    """

    al1 = get_pairs(args.align1)
    al2 = get_pairs(args.align2)

    if args.score == 'sp':
        sp_score(al1, al2)
    if args.score == 'f1':
        f1_score(al1, al2)


def main():
    """Initialize two alignments and compare them.
    """

    parser = argparse.ArgumentParser()
    parser.add_argument('-align1', type=str, help='First alignment')
    parser.add_argument('-align2', type=str, help='Second alignment')
    parser.add_argument('-score', type=str, default='sp', help='Comparison score (sp/f1)')
    args = parser.parse_args()

    compare_aligns(args)


if __name__ == '__main__':
    main()
