"""================================================================================================
This script aligns two protein sequences using either PEbA or BLOSUM and plots the distribution
of scores for the alignments. Place two fasta files and their embeddings to output a graph.

Ben Iovino  04/18/23   VecAligns
================================================================================================"""

import sys
import os
import argparse
import blosum as bl
import matplotlib.pyplot as plt
import numpy as np
sys.path.append(os.path.abspath(os.path.join(os.path.pardir, 'VecAligns')))
from utility import parse_fasta


def peba_align(seq1, seq2, vecs1, vecs2, gopen, gext):
    """=============================================================================================
    This function accepts two sequences, creates a matrix corresponding to their lengths, and
    calculates the score of the alignments for each index. A second matrix is scored so that the
    best alignment can be tracebacked.

    :param seq1: first sequence
    :param seq2: second sequence
    :param vecs1: first sequence's amino acid vectors
    :param vecs2: second sequence's amino acid vectors
    :param gopen: gap penalty for opening a new gap
    :param gext: gap penalty for extending a gap
    return: list of scores from matching pairs
    ============================================================================================="""

    # Initialize scoring and traceback matrix based on sequence lengths
    row_length = len(seq1)+1
    col_length = len(seq2)+1
    score_m = np.full((row_length, col_length), 0)
    trace_m = np.full((row_length, col_length), 0)

    # List of scores from matching pairs
    scores = []

    # Score matrix by moving through each index
    gap = False
    for i in range(len(seq1)):
        seq1_vec = vecs1[i]  # Corresponding amino acid vector in 1st sequence
        for j in range(len(seq2)):

            # Preceding scoring matrix values
            diagonal = score_m[i][j]
            horizontal = score_m[i+1][j]
            vertical = score_m[i][j+1]

            # Score pair of residues based off cosine similarity
            seq2_vec = vecs2[j]  # Corresponding amino acid vector in 2nd sequence
            cos_sim = np.dot(seq1_vec,seq2_vec)/(np.linalg.norm(seq1_vec)*np.linalg.norm(seq2_vec))
            cos_sim *= 10

            # Add to scoring matrix values via scoring method
            diagonal += cos_sim
            if gap is False:  # Apply gap open penalty if there is no gap
                horizontal += gopen
                vertical += gopen
            if gap is True:  # Apply gap extension penalty if there is a gap
                horizontal += gext
                vertical += gext

            # Assign value to traceback matrix and update gap status
            score = max(diagonal, horizontal, vertical)
            if score == diagonal:
                scores.append(cos_sim)
                trace_m[i+1][j+1] = 0
                gap = False
            if score == horizontal:
                trace_m[i+1][j+1] = -1
                gap = True
            if score == vertical:
                trace_m[i+1][j+1] = 1
                gap = True

            # Assign max value to scoring matrix
            score_m[i+1][j+1] = max(score, 0)

    return scores


def matrix_align(seq1, seq2, subs_matrix, gopen, gext):
    """=============================================================================================
    This function accepts two sequences, creates a matrix corresponding to their lengths, and
    calculates the score of the alignments for each index. A second matrix is scored so that the
    best alignment can be tracebacked.

    :param seq1: first sequence
    :param seq2: second sequence
    :param subs_matrix: substitution scoring matrix (i.e. BLOSUM62)
    :param gopen: gap penalty for opening a new gap
    :param gext: gap penalty for extending a gap
    return: list of scores from matching pairs
    ============================================================================================="""

    # Initialize scoring and traceback matrix based on sequence lengths
    row_length = len(seq1)+1
    col_length = len(seq2)+1
    score_m = np.full((row_length, col_length), 0)
    trace_m = np.full((row_length, col_length), 0)

    # List of scores from matching pairs
    scores = []

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
            if gap is False:  # Apply gap open penalty if there is no gap
                horizontal += gopen
                vertical += gopen
            if gap is True:  # Apply gap extension penalty if there is a gap
                horizontal += gext
                vertical += gext

            # Assign value to traceback matrix and update gap status
            score = max(diagonal, horizontal, vertical)
            if score == diagonal:
                scores.append(matrix_score)
                trace_m[i+1][j+1] = 0
                gap = False
            if score == horizontal:
                trace_m[i+1][j+1] = -1
                gap = True
            if score == vertical:
                trace_m[i+1][j+1] = 1
                gap = True

            # Assign max value to scoring matrix
            score_m[i+1][j+1] = max(score, 0)

    return scores

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-file1', type=str, default='./test1.fa', help='Name of first fasta file')
    parser.add_argument('-file2', type=str, default='./test2.fa', help='Name of second fasta file')
    parser.add_argument('-embed1', type=str, default='./test1.txt', help='Name of first embedding')
    parser.add_argument('-embed2', type=str, default='./test2.txt', help='Name of second embedding')
    parser.add_argument('-gopen', type=float, default=-11, help='Penalty for opening a gap')
    parser.add_argument('-gext', type=float, default=-1, help='Penalty for extending a gap')
    parser.add_argument('-encoder', type=str, default='ProtT5', help='Encoder to use')
    parser.add_argument('-score', type=int, default=62, help='log odds score of subsitution matrix')
    args = parser.parse_args()

    # Load fasta files and ids
    seq1, id1 = parse_fasta(args.file1)  #pylint: disable=W0612
    seq2, id2 = parse_fasta(args.file2)  #pylint: disable=W0612

    # Load embedding files
    vecs1 = np.loadtxt(args.embed1)
    vecs2 = np.loadtxt(args.embed2)

    # Load substitution matrix
    matrix = bl.BLOSUM(args.score)

    # Call align functions to get scores
    peba_scores = peba_align(seq1, seq2, vecs1, vecs2, args.gopen, args.gext)
    blosum_scores = matrix_align(seq1, seq2, matrix, args.gopen, args.gext)

    # Plot histogram of scores
    plt.hist(peba_scores, bins=10, alpha=0.5, label='Cosine Similarity')
    plt.hist(blosum_scores, bins=10, alpha=0.5, label=f'BLOSUM{args.score}')
    plt.legend(loc='upper right')
    plt.title("Scores of Matching Pairs")
    plt.xlabel("Score")
    plt.ylabel("Frequency")
    plt.show()


if __name__ == '__main__':
    main()
