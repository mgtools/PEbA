"""================================================================================================
This script takes two protein sequences of varying length and finds the highest scoring local
alignment between the two.

Ben Iovino  01/23/23   VecAligns
================================================================================================"""

import argparse
import numpy as np
from Bio import SeqIO


def parse_fasta(filename):
    """=============================================================================================
    This function accepts a fasta file name and returns the sequence.

    :param filename: name of file
    return: sequence
    ============================================================================================="""

    # Parse fasta file
    seq = ''
    with open(filename, 'r', encoding='utf8') as file:
        for seq in SeqIO.parse(file, 'fasta'):
            seq = str(seq.seq)
    return seq


def local_align(seq1, seq2, subs_matrix, gopen, gext):
    """=============================================================================================
    This function accepts two sequences, creates a matrix corresponding to their lengths, and
    calculates the score of the alignments for each index. A second matrix is scored so that the
    best alignment can be tracebacked.

    :param seq1: first sequence
    :param seq2: second sequence
    :param subs_matrix: substitution scoring matrix (i.e. BLOSUM62)
    :param gopen: gap penalty for opening a new gap
    :param gext: gap penalty for extending a gap
    return: scoring and traceback matrices of optimal scores for the SW-alignment of sequences
    ============================================================================================="""

    # Protein alphabet
    chars = 'ACDEFGHIKLMNPQRSTVWY'

    # Initialize scoring and traceback matrix based on sequence lengths
    row_length = len(seq1)+1
    col_length = len(seq2)+1
    score_m = np.full((row_length, col_length), 0)
    trace_m = np.full((row_length, col_length), 0)

    # Score matrix by moving through each index
    gap = False
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
            if gap is False:  # Apply gap open penalty if there is no gap
                horizontal += gopen
                vertical += gopen
            if gap is True:  # Apply gap extension penalty if there is a gap
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

            # Assign max value to scoring matrix
            score_m[i+1][j+1] = max(score, 0)

    return score_m, trace_m


def write_align(seq1, seq2):
    """=============================================================================================
    This function accepts two sequences after gaps have been introduced and writes them to a file
    in no particular format (yet).

    :param seq1: first aligned sequence
    :param seq2: second aligned sequence
    ============================================================================================="""

    # Add space every 10 characters
    seq1 = [seq1[i:i+10] for i in range(0, len(seq1), 10)]
    seq1 = ' '.join(seq1)
    seq2 = [seq2[i:i+10] for i in range(0, len(seq2), 10)]
    seq2 = ' '.join(seq2)

    # Split sequences every 50 characters
    seq1_split = [seq1[i:i+55] for i in range(0, len(seq1), 55)]
    seq2_split = [seq2[i:i+55] for i in range(0, len(seq2), 55)]

    # Find max length sequence and write to file based on its length
    name1 = 'seque1'
    name2 = 'seque2'
    with open('local_alignment.txt', 'w', encoding='utf8') as file:
        file.write('PileUp\n\n\n')
        file.write(f'   MSF:  {len(seq1)}  Type:  P\n\n')
        file.write(f'Name: {name1} oo  Len:  {len(seq1)}\n')
        file.write(f'Name: {name2} oo  Len:  {len(seq2)}\n\n//\n\n\n\n')
        for i in range(len(seq1_split)):
            file.write(f'{name1}      {seq1_split[i]}\n')
            file.write(f'{name2}      {seq2_split[i]}\n\n')


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

    # Join lists and reverse strings again
    seq1 = ''.join(rev_seq1)
    seq2 = ''.join(rev_seq2)
    seq1 = seq1[::-1]
    seq2 = seq2[::-1]

    # Introduce gaps at beginning of either sequence based off final index positions
    seq1 = "."*index[1]+seq1
    seq2 = "."*index[0]+seq2

    # Introduce gaps at end of either sequence based off length of other sequence
    seq1 = seq1+"."*max(0, len(seq2)-len(seq1))
    seq2 = seq2+"."*max(0, len(seq1)-len(seq2))
    write_align(seq1, seq2)


def main():
    """=============================================================================================
    This function initializes the BLOSUM62 matrix and two protein sequences, calls SW_align() to get
    the scoring and traceback matrix from SW alignment, and then calls traceback() to print the
    local alignment.
    ============================================================================================="""

    # Take fasta sequences for arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-file1', type=str, default='test1.fa', help='Name of first fasta file')
    parser.add_argument('-file2', type=str, default='test2.fa', help='Name of second fasta file')
    parser.add_argument('-gopen', type=int, default=-3, help='Penalty for opening a gap')
    parser.add_argument('-gext', type=int, default=-1, help='Penalty for extending a gap')
    args = parser.parse_args()

    # Parse fasta files
    seq1 = parse_fasta(args.file1)
    seq2 = parse_fasta(args.file2)

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

    # Call local_align() to get scoring and traceback matrix
    score_m, trace_m = local_align(seq1, seq2, blosum, args.gopen, args.gext)

    # Call traceback() to get highest scoring local alignment between seq1 and seq2
    traceback(score_m, trace_m, seq1, seq2)


if __name__ == '__main__':
    main()