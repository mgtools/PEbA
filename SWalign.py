"""================================================================================================
This script takes two protein sequences of varying length and finds the highest scoring local
alignment between the two.

Ben Iovino  01/23/23   VecAligns
================================================================================================"""

import argparse
import re
import numpy as np


def local_align(seq1, seq2, subs_matrix):
    """=============================================================================================
    This function accepts two sequences, creates a matrix corresponding to their lengths, and  
    calculates the score of the alignments for each index. A second matrix is scored so that the
    best alignment can be tracebacked.

    :param seq1: first sequence
    :param seq2: second sequence
    :param subs_matrix: substitution scoring matrix (i.e. BLOSUM62)
    return: scoring and traceback matrices of optimal scores for the SW-alignment of sequences
    ============================================================================================="""

    # NCBI default gap costs
    gap_open = -11
    gap_ext = -1
    gap = False

    # Protein alphabet
    chars = 'ACDEFGHIKLMNPQRSTVWY'

    # Initialize scoring and traceback matrix based on sequence lengths
    row_length = len(seq1)+1
    col_length = len(seq2)+1
    score_m = np.full((row_length, col_length), 0)
    trace_m = np.full((row_length, col_length), 0)

    # Score matrix by moving through each index
    for i, char in enumerate(seq1):
        seq1_char = char  # Character in 1st sequence
        seq1_index = chars.index(seq1_char)  # Corresponding row in BLOSUM matrix
        for j, char in enumerate(seq2):
            seq2_char = char

            # Preceding scoring matrix values
            diagonal = score_m[i][j]
            horizontal = score_m[i+1][j]
            vertical = score_m[i][j+1]

            # Score residues based off BLOSUM matrix
            # print(f'Char1: {seq1_char}, Char2: {seq2_char}, BLOSUM score: {score}')
            seq2_index = chars.index(seq2_char)  # Corresponding column in BLOSUM matrix
            matrix_score = subs_matrix[seq2_index][seq1_index]

            # Add to matrix values via scoring method
            diagonal += matrix_score
            if gap is False:  # Apply gap_open penalty if there is no gap
                horizontal += gap_open
                vertical += gap_open
            if gap is True:  # Apply gap_extension penalty if there is a gap
                horizontal += gap_ext
                vertical += gap_ext

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
                trace_m[i+1][j+1] = -1

            # Assign max value to scoring matrix
            score_m[i+1][j+1] = max(score, 0)

    return score_m, trace_m


def traceback(score_m, trace_m, seq1, seq2):
    """=============================================================================================
    This function accepts a scoring and a traceback matrix and two sequences and returns the highest 
    scoring local alignment between the two sequences

    :param score_m: scoring matrix
    :param trace_m: traceback matrix
    :param seq1: first sequence
    :param seq2: second sequence
    return: highest scoring local alignment of the two sequences
    ============================================================================================="""

    # Find index of highest score in scoring matrix, start traceback at this matrix
    high_score_ind = np.unravel_index(np.argmax(score_m, axis=None), score_m.shape)

    # Reverse strings and convert to lists so gaps can be inserted
    rev_seq1 = list(seq1[::-1])
    rev_seq2 = list(seq2[::-1])

    # Move through matrix starting at bottom right
    index = [high_score_ind[0], high_score_ind[1]]
    count = 0
    while (index[0] and index[1]) != 0:
        val = trace_m[index[0], index[1]]

        if val == 1:  # If cell is equal to 1, insert a gap into the second sequence
            index[0] = index[0] - 1
            rev_seq2.insert(count, '_')
        if val == -1:  # If cell is equal to -1, insert a gap into the first sequence
            index[1] = index[1] - 1
            rev_seq1.insert(count, '_')
        if val == 0:  # If cell is equal to 0, there is no gap
            index[0] = index[0] - 1
            index[1] = index[1] - 1
        count += 1

    # Join lists and reverse strings again
    seq1 = ''.join(rev_seq1)
    seq2 = ''.join(rev_seq2)

    # Store alignment
    space = ' '
    print(index[1]*space+seq1[::-1])
    print(index[0]*space+seq2[::-1])


def main():
    """=============================================================================================
    This function initializes the BLOSUM62 matrix and two protein sequences, calls SW_align() to get
    the scoring and traceback matrix from SW alignment, and then calls traceback() to print the
    local alignment.
    ============================================================================================="""

    # Take fasta sequences for arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-seq1', type=str, default='AGFISVISKKQGEYLEDEWY')
    parser.add_argument('-seq2', type=str, default='QVLDKFGS')
    args = parser.parse_args()

    # Intialize BLOSUM62 matrix
    blosum = [[4,0,-2,-1,-2,0,-2,-1,-1,-1,-1,-2,-1,-1,-1,1,0,0,-3,-2],
    [0,9,-3,-4,-2,-3,-3,-1,-3,-1,-1,-3,-3,-3,-3,-1,-1,-1,-2,-2],
    [-2,-3,6,2,-3,-1,-1,-3,-1,-4,-3,1,-1,0,-2,0,-1,-3,-4,-3],
    [-1,-4,2,5,-3,-2,0,-3,1,-3,-2,0,-1,2,0,0,-1,-2,-3,-2],
    [-2,-2,-3,-3,6,-3,-1,0,-3,0,0,-3,-4,-3,-3,-2,-2,-1,1,3],
    [0,-3,-1,-2,-3,6,-2,-4,-2,-4,-3,0,-2,-2,-2,0,-2,-3,-2,-3],
    [-2,-3,-1,0,-1,-2,8,-3,-1,-3,-2,1,-2,0,0,-1,-2,-3,-2,2],
    [-1,-1,-3,-3,0,-4,-3,4,-3,2,1,-3,-3,-3,-3,-2,-1,3,-3,-1],
    [-1,-3,-1,1,-3,-2,-1,-3,5,-2,-1,0,-1,1,2,0,-1,-2,-3,-2],
    [-1,-1,-4,-3,0,-4,-3,2,-2,4,2,-3,-3,-2,-2,-2,-1,1,-2,-1],
    [-1,-1,-3,-2,0,-3,-2,1,-1,2,5,-2,-2,0,-1,-1,-1,1,-1,-1],
    [-2,-3,1,0,-3,0,1,-3,0,-3,-2,6,-2,0,0,1,0,-3,-4,-2],
    [-1,-3,-1,-1,-4,-2,-2,-3,-1,-3,-2,-2,7,-1,-2,-1,-1,-2,-4,-3],
    [-1,-3,0,2,-3,-2,0,-3,1,-2,0,0,-1,5,1,0,-1,-2,-2,-1],
    [-1,-3,-2,0,-3,-2,0,-3,2,-2,-1,0,-2,1,5,-1,-1,-3,-3,-2],
    [1,-1,0,0,-2,0,-1,-2,0,-2,-1,1,-1,0,-1,4,1,-2,-3,-2],
    [0,-1,-1,-1,-2,-2,-2,-1,-1,-1,-1,0,-1,-1,-1,1,5,0,-2,-2],
    [0,-1,-3,-2,-1,-3,-3,3,-2,1,1,-3,-2,-2,-3,-2,0,4,-3,-1],
    [-3,-2,-4,-3,1,-2,-2,-3,-3,-2,-1,-4,-4,-2,-3,-3,-2,-3,11,2],
    [-2,-2,-3,-2,3,-3,2,-1,-2,-1,-1,-2,-3,-1,-2,-2,-2,-1,2,7]]

    # Call SW_align() to get scoring and traceback matrix
    score_m, trace_m = local_align(args.seq1, args.seq2, blosum)

    # Call traceback() to get highest scoring local alignment between seq1 and seq2
    traceback(score_m, trace_m, args.seq1, args.seq2)


if __name__ == '__main__':
    main()
