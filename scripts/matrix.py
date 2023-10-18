"""This script takes two protein sequences of varying length and finds the highest scoring
alignment between the two using a substitution matrix.

__author__ = "Ben Iovino"
__date__ = 09/19/23
"""

import argparse
import blosum as bl
import utility as ut


def score_align(seq1: str, seq2: str, subs_matrix, gopen: float, gext: float, align: str) -> tuple:
    """Returns scoring and traceback matrices of optimal scores for the alignment of sequences

    :param seq1: first sequence
    :param seq2: second sequence
    :param subs_matrix: substitution scoring matrix (i.e. BLOSUM62)
    :param gopen: gap penalty for opening a new gap
    :param gext: gap penalty for extending a gap
    :param align: alignment type (global or local)
    return (np.ndarray, np.ndarray):
      scoring and traceback matrices of optimal scores for the NW or SW alignment of sequences
    """

    # Initialize scoring and traceback matrix based on sequence lengths
    score_m, trace_m = ut.initialize_matrices(seq1, seq2, align, gopen, gext)

    # Score matrix by moving through each index
    gap = False
    for i, char in enumerate(seq1):
        seq1_char = char  # Character in 1st sequence
        for j, char in enumerate(seq2):
            seq2_char = char  # Character in 2nd sequence

            # Preceding scoring matrix values
            diagonal = score_m[i][j]
            horizontal = score_m[i+1][j]
            vertical = score_m[i][j+1]

            # Get score from substitution matrix and add to scoring matrix values
            matrix_score = subs_matrix[f'{seq1_char}{seq2_char}']
            diagonal += matrix_score
            horizontal, vertical = ut.gap_penalty(gap, horizontal, vertical, gopen, gext)

            # Assign value to traceback matrix and update gap status
            score = max(diagonal, horizontal, vertical)
            trace_m, gap = ut.assign_score(score, diagonal, horizontal, vertical, trace_m, gap, i, j)

            # Assign max value to scoring matrix
            if align == 'local':
                score_m[i+1][j+1] = max(score, 0)
            if align == 'global':
                score_m[i+1][j+1] = score

    return score_m, trace_m


def main():
    """Initializes two protein sequences and a scoring matrix, calls SW_align() to get
    the scoring and traceback matrix from SW alignment, calls traceback() to get the local
    alignment, and then writes the alignment to a file in the desired format.
    """

    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--align', type=str, default='local', help='Alignment algorithm to use (global or local)')
    parser.add_argument('-f1', '--file1', type=str, help='Name of first fasta file')
    parser.add_argument('-f2', '--file2', type=str, help='Name of second fasta file')
    parser.add_argument('-go', '--gopen', type=float, default=-11, help='Penalty for opening a gap')
    parser.add_argument('-ge', '--gext', type=float, default=-1, help='Penalty for extending a gap')
    parser.add_argument('-m', '--matrix', type=str, default='blosum', help='Substitution matrix to use (blosum or pfasum)')
    parser.add_argument('-s', '--score', type=int, default=62, help='Log odds score of subsitution matrix')
    parser.add_argument('-o', '--output', type=str, default='msf', help='Output format (msf or fa)')
    parser.add_argument('-sf', '--savefile', type=str, help='Filename to save alignment to')
    args = parser.parse_args()

    # Parse fasta files for sequences and ids
    seq1, id1 = ut.parse_fasta(args.file1)
    seq2, id2 = ut.parse_fasta(args.file2)

    # Intialize scoring matrix
    if args.matrix == 'blosum':
        matrix = bl.BLOSUM(args.score)
    if args.matrix == 'pfasum':
        matrix = ut.parse_matrix('data/PFASUM60.txt')

    # Align and traceback
    score_m, trace_m = score_align(seq1, seq2, matrix, args.gopen, args.gext, args.align)
    if args.align == 'global':
        align1, align2 = ut.global_traceback(trace_m, seq1, seq2)
        beg, end = [0, 0], [len(seq1), len(seq2)]
    if args.align == 'local':
        align1, align2, beg, end = ut.local_traceback(score_m, trace_m, seq1, seq2)

    # Write align based on desired output format
    if args.output == 'msf':
        ut.write_msf(align1, align2, id1, id2, args.matrix+str(args.score),
                args.gopen, args.gext, args.savefile, beg, end)
    if args.output == 'fa':
        ut.write_fasta(align1, align2, id1, id2, args.savefile)


if __name__ == '__main__':
    main()
