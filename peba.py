"""================================================================================================
This script takes two protein sequences of any length and finds the highest scoring local
alignment between the two using the cosine similarity between embedded amino acids.

Ben Iovino  01/23/23   VecAligns
================================================================================================"""

import argparse
import torch
import numpy as np
from transformers import T5EncoderModel, T5Tokenizer
from utility import parse_fasta, write_msf, write_fasta
from embed_seqs import prot_t5xl_embed


def local_align(seq1: str, seq2: str, vecs1: list, vecs2:list, gopen: float, gext: float):
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
    :return list: scoring and traceback matrices of optimal scores for the SW-alignment of sequences
    ============================================================================================="""

    # Initialize scoring and traceback matrix based on sequence lengths
    row_length = len(seq1)+1
    col_length = len(seq2)+1
    score_m = np.full((row_length, col_length), 0)
    trace_m = np.full((row_length, col_length), 0)

    # Score matrix by moving through each index
    gap = False
    for i in range(len(seq1)):
        seq1_vec = vecs1[i]  # Corresponding amino acid vector in 1st sequence
        seq1_norm = np.linalg.norm(seq1_vec)
        for j in range(len(seq2)):

            # Preceding scoring matrix values
            diagonal = score_m[i][j]
            horizontal = score_m[i+1][j]
            vertical = score_m[i][j+1]

            # Score pair of residues based off cosine similarity
            seq2_vec = vecs2[j]  # Corresponding amino acid vector in 2nd sequence
            cos_sim = np.dot(seq1_vec,seq2_vec)/(seq1_norm*np.linalg.norm(seq2_vec))
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

    return score_m, trace_m


def traceback(score_m: list, trace_m: list, seq1: str, seq2: str):
    """=============================================================================================
    This function accepts a scoring and a traceback matrix and two sequences and returns the highest
    scoring local alignment between the two sequences.

    :param score_m: scoring matrix
    :param trace_m: traceback matrix
    :param seq1: first sequence
    :param seq2: second sequence
    return: seq1 with gaps, seq2 with gaps
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
        val = trace_m[index[0], index[1]]  # \\NOSONAR

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
    align1 = seq1+"."*max(0, len(seq2)-len(seq1))
    align2 = seq2+"."*max(0, len(seq1)-len(seq2))
    return align1, align2


def main():
    """=============================================================================================
    Main initializes two protein sequences, embeds them if embeddings are not provided, calls
    local_align() to obtain the scoring and traceback matrix from local alignment (with cosine
    similarity as the scoring method), calls traceback() to get the highest scoring local alignment,
    and then writes the alignment in the desired format to either the console or to a specified file.
    ============================================================================================="""

    parser = argparse.ArgumentParser()
    parser.add_argument('-f1', '--file1', type=str, help='First fasta file')
    parser.add_argument('-f2', '--file2', type=str, help='Second fasta file')
    parser.add_argument('-e1', '--embed1', type=str, default='n', help='First embedding file')
    parser.add_argument('-e2', '--embed2', type=str, default='n', help='Second embedding file')
    parser.add_argument('-go', '--gopen', type=float, default=-11, help='Gap open penalty')
    parser.add_argument('-ge', '--gext', type=float, default=-1, help='Gap extension penalty')
    parser.add_argument('-e', '--encoder', type=str, default='ProtT5', help='Encoder used')
    parser.add_argument('-o', '--output', type=str, default='msf', help='Output format')
    parser.add_argument('-s', '--savefile', type=str, default='n', help='File to save alignment')
    args = parser.parse_args()

    # Load fasta files and ids
    seq1, id1 = parse_fasta(args.file1)
    seq2, id2 = parse_fasta(args.file2)

    # Load models, embed sequences
    if args.embed1=='n':
        device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')  #pylint: disable=E1101
        tokenizer = T5Tokenizer.from_pretrained('Rostlab/prot_t5_xl_uniref50', do_lower_case=False)
        model = T5EncoderModel.from_pretrained("Rostlab/prot_t5_xl_uniref50")
        model.to(device)  # Loads to GPU if available
        vecs1 = prot_t5xl_embed(seq1, tokenizer, model)
        vecs2 = prot_t5xl_embed(seq2, tokenizer, model)

    # Load numpy arrays
    else:
        vecs1 = np.loadtxt(args.embed1)
        vecs2 = np.loadtxt(args.embed2)

    # Align and traceback
    score_m, trace_m = local_align(seq1, seq2, vecs1, vecs2, args.gopen, args.gext)
    align1, align2 = traceback(score_m, trace_m, seq1, seq2)

    # Write align based on desired output format
    if args.output == 'msf':
        write_msf(align1, align2, id1, id2, f'{args.encoder}_Sim',
                   args.gopen, args.gext, args.savefile)
    if args.output == 'fa':
        write_fasta(align1, align2, id1, id2, args.savefile)


if __name__ == '__main__':
    main()
