"""================================================================================================
This script takes two protein sequences of varying length and finds the highest scoring local
alignment between the two using a substitution matrix.

Ben Iovino  01/23/23   VecAligns
================================================================================================"""

import argparse
import numpy as np
import blosum as bl
from utility import parse_fasta, write_msf, parse_matrix, write_fasta


def initialize_matrices(seq1: str, seq2: str, align: str, gopen: float, gext: float) -> tuple:
    """Returns scoring and traceback matrices initialized based on alignment type

    :param seq1: first sequence
    :param seq2: second sequence
    :param align: alignment type (global or local)
    :param gopen: gap open penalty for global matrix intialization
    :param gext: gap extension penalty for global matrix intialization
    """

    # Initialize scoring and traceback matrix based on sequence lengths
    row_length = len(seq1)+1
    col_length = len(seq2)+1
    score_m = np.full((row_length, col_length), 0)
    trace_m = np.full((row_length, col_length), 0)

    # Initialize scoring matrix based on alignment type
    if align == 'global':
        for i in range(1, len(score_m[0])):
            score_m[0][i] = gopen+gext*i+1  # +1 to offset i starting at 1
            trace_m[0][i] = -1
        for i in range(1, len(score_m.T[0])):
            score_m.T[0][i] = gopen+gext*i+1
            trace_m.T[0][i] = 1

    return score_m, trace_m


def gap_penalty(gap: bool, horizontal: int, vertical: int, gopen: float, gext: float) -> tuple:
    """Returns score with gap penalty applied

    :param gap: boolean indicating if there is a gap
    :param horizontal: score if gap is horizontal
    :param vertical: score if gap is vertical
    :param gopen: gap open penalty
    :param gext: gap extension penalty
    """

    if gap is False:  # Apply gap open penalty if there is no gap
        horizontal += gopen
        vertical += gopen
    if gap is True:  # Apply gap extension penalty if there is a gap
        horizontal += gext
        vertical += gext

    return horizontal, vertical


def assign_score(score: int, diagonal: int, horizontal: int,
                  vertical: int, trace_m: np.ndarray, gap: bool, i: int, j: int) -> tuple:
    """Returns traceback matrix with value based on highest score and gap status.

    :param score: score of current residue comparison
    :param diagonal: score of diagonal cell
    :param horizontal: score of horizontal cell
    :param vertical: score of vertical cell
    :param trace_m: traceback matrix
    :param gap: boolean indicating if there is a gap
    :param i: row index of current cell
    :param j: column index of current cell
    return (np.ndarray, bool): traceback matrix with value based on highest score and gap status
    """

    if score == diagonal:
        trace_m[i+1][j+1] = 0
        gap = False
    if score == horizontal:
        trace_m[i+1][j+1] = -1
        gap = True
    if score == vertical:
        trace_m[i+1][j+1] = 1
        gap = True

    return trace_m, gap


def score_align(seq1: str, seq2: str, subs_matrix, gopen: float, gext: float, align: str) -> tuple:
    """Returns scoring and traceback matrices of optimal scores for the SW-alignment of sequences

    :param seq1: first sequence
    :param seq2: second sequence
    :param subs_matrix: substitution scoring matrix (i.e. BLOSUM62)
    :param gopen: gap penalty for opening a new gap
    :param gext: gap penalty for extending a gap
    :param align: alignment type (global or local)
    return: scoring and traceback matrices of optimal scores for the SW-alignment of sequences
    """

    # Initialize scoring and traceback matrix based on sequence lengths
    score_m, trace_m = initialize_matrices(seq1, seq2, align, gopen, gext)

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
            horizontal, vertical = gap_penalty(gap, horizontal, vertical, gopen, gext)

            # Assign value to traceback matrix and update gap status
            score = max(diagonal, horizontal, vertical)
            trace_m, gap = assign_score(score, diagonal, horizontal, vertical, trace_m, gap, i, j)

            # Assign max value to scoring matrix
            if align == 'local':
                score_m[i+1][j+1] = max(score, 0)
            if align == 'global':
                score_m[i+1][j+1] = score

    return score_m, trace_m


def get_aligns(rev_seq1: list, rev_seq2: list, index1: int, index0: int) -> tuple:
    """Returns aligned sequences with gaps inserted based on traceback matrix

    :param rev_seq1: first sequence reversed
    :param rev_seq2: second sequence reversed
    :param index1: final traceback matrix row index
    :param index0: final traceback matrix column index
    """

    seq1 = ''.join(rev_seq1)
    seq2 = ''.join(rev_seq2)
    seq1 = seq1[::-1]
    seq2 = seq2[::-1]

    # Introduce gaps at beginning of either sequence based off final index positions
    seq1 = "."*index1+seq1
    seq2 = "."*index0+seq2

    # Introduce gaps at end of either sequence based off length of other sequence
    align1 = seq1+"."*max(0, len(seq2)-len(seq1))
    align2 = seq2+"."*max(0, len(seq1)-len(seq2))

    return align1, align2


def global_traceback(trace_m: np.ndarray, seq1: str, seq2: str) -> tuple:
    """Returns global alignment of two sequences with gaps inserted based on traceback matrix

    :param score_m: scoring matrix
    :param trace_m: traceback matrix
    :param seq1: first sequence
    :param seq2: second sequence
    return (str, str): seq1 with gaps, seq2 with gaps
    """

    # Reverse strings and convert to lists so gaps can be inserted
    rev_seq1 = list(seq1[::-1])
    rev_seq2 = list(seq2[::-1])

    # Move through matrix starting at bottom right
    rows, cols = trace_m.shape
    index = [rows-1, cols-1]
    count = 0
    while index != [0, 0]:
        val = trace_m[index[0], index[1]]
        if val == 1:  # If cell is equal to 1, insert a gap into the second sequence
            index[0] = max(index[0] - 1, 0)  # Taking max of new index and 0 so index never below 0
            rev_seq2.insert(count, '.')
        if val == -1:  # If cell is equal to -1, insert a gap into the first sequence
            index[1] = max(index[1] - 1, 0)
            rev_seq1.insert(count, '.')
        if val == 0:  # If cell is equal to 0, there is no gap
            index[0] = max(index[0] - 1, 0)
            index[1] = max(index[1] - 1, 0)
        count += 1

    # Join lists and reverse strings again
    align1, align2 = get_aligns(rev_seq1, rev_seq2, index[1], index[0])

    return align1, align2


def local_traceback(score_m: np.ndarray, trace_m: np.ndarray, seq1: str, seq2: str) -> tuple:
    """Returns local alignment of two sequences with gaps inserted based on traceback matrix

    :param score_m: scoring matrix
    :param trace_m: traceback matrix
    :param seq1: first sequence
    :param seq2: second sequence
    return (str, str): seq1 with gaps, seq2 with gaps
    """

    # Find index of highest score in scoring matrix, start traceback at this matrix
    high_score_ind = np.unravel_index(np.argmax(score_m, axis=None), score_m.shape)

    # Reverse strings and convert to lists so gaps can be inserted
    rev_seq1 = list(seq1[::-1])
    rev_seq2 = list(seq2[::-1])

    # Move through matrix starting at highest scoring cell
    index = [high_score_ind[0], high_score_ind[1]]

    # Gaps are inserted based on count increasing as we move through sequences
    # If we don't start at bottom right, need to adjust position at which gaps inserted
    count_adjust1 = len(seq1) - high_score_ind[0]
    count_adjust2 = len(seq2) - high_score_ind[1]
    count = 0
    while (index[0] and index[1]) != 0:
        val = trace_m[index[0], index[1]]

        if val == 1:  # If cell is equal to 1, insert a gap into the second sequence
            index[0] = index[0] - 1
            rev_seq2.insert(count+count_adjust2, '.')
        if val == -1:  # If cell is equal to -1, insert a gap into the first sequence
            index[1] = index[1] - 1
            rev_seq1.insert(count+count_adjust1, '.')
        if val == 0:  # If cell is equal to 0, there is no gap
            index[0] = index[0] - 1
            index[1] = index[1] - 1
        count += 1

    # Get aligned sequences with gaps inserted
    align1, align2 = get_aligns(rev_seq1, rev_seq2, index[1], index[0])

    return align1, align2


def main():
    """Initializes two protein sequences and a scoring matrix, calls SW_align() to get
    the scoring and traceback matrix from SW alignment, calls traceback() to get the local
    alignment, and then writes the alignment to a file in the desired format.
    """

    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--align', type=str, help='Alignment algorithm to use i.e. global or local', default='local')
    parser.add_argument('-f1', '--file1', type=str, help='Name of first fasta file', default='Data/Example/1j46_A.fa')
    parser.add_argument('-f2', '--file2', type=str, help='Name of second fasta file', default='Data/Example/1k99_A.fa')
    parser.add_argument('-go', '--gopen', type=float, default=-11, help='Penalty for opening a gap')
    parser.add_argument('-ge', '--gext', type=float, default=-1, help='Penalty for extending a gap')
    parser.add_argument('-m', '--matrix', type=str, default='blosum', help='Substitution matrix to use i.e. blosum or pfasum')
    parser.add_argument('-s', '--score', type=int, default=62, help='Log odds score of subsitution matrix')
    parser.add_argument('-o', '--output', type=str, default='msf', help='Output format i.e. msf or fa')
    parser.add_argument('-sf', '--savefile', type=str, default='n', help='Filename to save alignment to')
    args = parser.parse_args()

    # Parse fasta files for sequences and ids
    seq1, id1 = parse_fasta(args.file1)
    seq2, id2 = parse_fasta(args.file2)

    # Intialize scoring matrix
    if args.matrix == 'blosum':
        matrix = bl.BLOSUM(args.score)
    if args.matrix == 'pfasum':
        matrix = parse_matrix('Data/PFASUM60.txt')

    # Align and traceback
    score_m, trace_m = score_align(seq1, seq2, matrix, args.gopen, args.gext, args.align)
    if args.align == 'global':
        align1, align2 = global_traceback(trace_m, seq1, seq2)
    if args.align == 'local':
        align1, align2 = local_traceback(score_m, trace_m, seq1, seq2)

    # Write align based on desired output format
    if args.output == 'msf':
        write_msf(align1, align2, id1, id2, args.matrix+str(args.score),
                args.gopen, args.gext, args.savefile)
    if args.output == 'fa':
        write_fasta(align1, align2, id1, id2, args.savefile)


if __name__ == '__main__':
    main()
