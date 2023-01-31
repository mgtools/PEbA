"""================================================================================================
This script takes two protein sequences of varying length and finds the highest scoring global
alignment between the two.

Ben Iovino  01/25/23   VecAligns
================================================================================================"""

import argparse
import numpy as np
import blosum as bl
from utility import parse_fasta, write_align  # pylint: disable=E0611


def global_align(seq1, seq2, subs_matrix, gopen, gext):
    """=============================================================================================
    This function accepts two sequences, creates a matrix corresponding to their lengths, and
    calculates the score of the alignments for each index. A second matrix is scored so that the
    best alignment can be tracebacked.

    :param seq1: first sequence
    :param seq2: second sequence
    :param subs_matrix: substitution scoring matrix (i.e. BLOSUM62)
    :param gopen: gap penalty for opening a new gap
    :param gext: gap penalty for extending a gap
    return: traceback matrix
    ============================================================================================="""

    # Initialize scoring and traceback matrix based on sequence lengths
    row_length = len(seq1)+1
    col_length = len(seq2)+1
    score_m = np.full((row_length, col_length), 0)
    trace_m = np.full((row_length, col_length), 0)

    # Initialize first row and column with gap values for S matrix, traceback values for T matrix
    for i in range(1, len(score_m[0])):
        score_m[0][i] = gopen+gext*i+1  # +1 to offset i starting at 1
        trace_m[0][i] = -1
    for i in range(1, len(score_m.T[0])):
        score_m.T[0][i] = gopen+gext*i+1
        trace_m.T[0][i] = 1

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

            # Score pair of residues based off BLOSUM matrix
            matrix_score = subs_matrix[f'{seq1_char}{seq2_char}']

            # Add to matrix values via scoring method
            diagonal += matrix_score
            if gap is False:  # Apply gap_open penalty if there is no gap
                horizontal += gopen
                vertical += gopen
            if gap is True:  # Apply gap_extension penalty if there is a gap
                horizontal += gext
                vertical += gext

            # Update gap status
            score = max(diagonal, horizontal, vertical)
            if score == horizontal:
                gap = True
            if score == vertical:
                gap = True
            if score == diagonal:
                gap = False

            # Assign value to traceback matrix
            if score == diagonal:
                trace_m[i+1][j+1] = 0
            if score == horizontal:
                trace_m[i+1][j+1] = -1
            if score == vertical:
                trace_m[i+1][j+1] = 1

            # Assign value to scoring matrix
            score_m[i+1][j+1] = score

    return trace_m


def traceback(trace_m, seq1, seq2):
    """=============================================================================================
    This function accepts a scoring and a traceback matrix and two sequences and returns global
    alignment between the two sequences

    :param trace_m: traceback matrix
    :param seq1: first sequence
    :param seq2: second sequence
    return: seq1 with gaps, seq2 with gaps
    ============================================================================================="""

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
    seq1 = ''.join(rev_seq1)
    seq2 = ''.join(rev_seq2)
    seq1 = seq1[::-1]
    seq2 = seq2[::-1]

    # Introduce gaps at end of either sequence based off length of other sequence
    align1 = seq1+"."*max(0, len(seq2)-len(seq1))
    align2 = seq2+"."*max(0, len(seq1)-len(seq2))
    return align1, align2


def main():
    """=============================================================================================
    This function initializes the BLOSUM62 matrix and two protein sequences, calls global_align() to
    get the scoring and traceback matrix from NW alignment, calls traceback() to get the global
    alignment, and then write_align() to write the alignment to a file in MSF format.
    ============================================================================================="""

    parser = argparse.ArgumentParser()
    parser.add_argument('-file1', type=str, default='test1.fa', help='Name of first fasta file')
    parser.add_argument('-file2', type=str, default='test2.fa', help='Name of second fasta file')
    parser.add_argument('-gopen', type=int, default=-3, help='Penalty for opening a gap')
    parser.add_argument('-gext', type=int, default=-1, help='Penalty for extending a gap')
    parser.add_argument('-blosum', type=int, default=62, help='BLOSUM matrix to use')
    args = parser.parse_args()

    # Parse fasta files for sequences and ids
    seq1, id1 = parse_fasta(args.file1)
    seq2, id2 = parse_fasta(args.file2)

    # Intialize BLOSUM matrix
    blosum = bl.BLOSUM(args.blosum)

    # Call global_align() to get traceback matrix
    trace_m = global_align(seq1, seq2, blosum, args.gopen, args.gext)

    # Get global alignment between seq1 and seq2 and write to file
    align1, align2 = traceback(trace_m, seq1, seq2)
    write_align(align1, align2, id1, id2)  # pylint: disable=E1121


if __name__ == '__main__':
    main()
