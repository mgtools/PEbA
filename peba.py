"""This script takes two protein sequences of any length and finds the highest scoring alignment
between the two using the cosine similarity between embedded amino acids.

__author__ = "Ben Iovino"
__date__ = 09/19/23
"""

import argparse
import torch
import numpy as np
from transformers import T5EncoderModel, T5Tokenizer
from embed_seqs import prot_t5xl_embed
import utility as ut


def score_align(seq1: str, seq2: str, vecs1: list, vecs2:list,
                 gopen: float, gext: float, align: str) -> tuple:
    """Returns scoring and traceback matrices of optimal scores for the alignment of sequences

    :param seq1: first sequence
    :param seq2: second sequence
    :param vecs1: first sequence's amino acid vectors
    :param vecs2: second sequence's amino acid vectors
    :param gopen: gap penalty for opening a new gap
    :param gext: gap penalty for extending a gap
    :param align: alignment type (global or local)
    :return (np.ndarray, np.ndarray):
      scoring and traceback matrices of optimal scores for the NW or SW alignment of sequences
    """

    # Initialize scoring and traceback matrix based on sequence lengths
    score_m, trace_m = ut.initialize_matrices(seq1, seq2, align, gopen, gext)

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
    """Initializes two protein sequences, embeds them if embeddings are not provided, calls
    local_align() to obtain the scoring and traceback matrix from local alignment (with cosine
    similarity as the scoring method), calls traceback() to get the highest scoring local alignment,
    and then writes the alignment in the desired format to either the console or to a specified file.
    """

    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--align', type=str, default='local', help='Alignment algorithm to use (global or local)')
    parser.add_argument('-f1', '--file1', type=str, help='First fasta file', default='Data/Example/1j46_A.fa')
    parser.add_argument('-f2', '--file2', type=str, help='Second fasta file', default='Data/Example/1k99_A.fa')
    parser.add_argument('-e1', '--embed1', type=str, default='Data/Example/1j46_A.txt', help='First embedding file')
    parser.add_argument('-e2', '--embed2', type=str, default='Data/Example/1k99_A.txt', help='Second embedding file')
    parser.add_argument('-go', '--gopen', type=float, default=-11, help='Gap open penalty')
    parser.add_argument('-ge', '--gext', type=float, default=-1, help='Gap extension penalty')
    parser.add_argument('-e', '--encoder', type=str, default='ProtT5', help='Encoder used')
    parser.add_argument('-o', '--output', type=str, default='msf', help='Output format')
    parser.add_argument('-s', '--savefile', type=str, default='n', help='File to save alignment')
    args = parser.parse_args()

    # Load fasta files and ids
    seq1, id1 = ut.parse_fasta(args.file1)
    seq2, id2 = ut.parse_fasta(args.file2)

    # Load models, embed sequences
    if args.embed1=='n':
        device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')  #pylint: disable=E1101
        tokenizer = T5Tokenizer.from_pretrained('Rostlab/prot_t5_xl_uniref50', do_lower_case=False)
        model = T5EncoderModel.from_pretrained("Rostlab/prot_t5_xl_uniref50")
        model.to(device)  # Loads to GPU if available
        vecs1 = prot_t5xl_embed(seq1, tokenizer, model, device)
        vecs2 = prot_t5xl_embed(seq2, tokenizer, model, device)

    # Load numpy arrays
    else:
        vecs1 = np.loadtxt(args.embed1)
        vecs2 = np.loadtxt(args.embed2)

    # Align and traceback
    score_m, trace_m = score_align(seq1, seq2, vecs1, vecs2, args.gopen, args.gext, args.align)
    if args.align == 'global':
        align1, align2 = ut.global_traceback(trace_m, seq1, seq2)
    if args.align == 'local':
        align1, align2 = ut.local_traceback(score_m, trace_m, seq1, seq2)

    # Write align based on desired output format
    if args.output == 'msf':
        ut.write_msf(align1, align2, id1, id2, f'{args.encoder}_Sim',
                   args.gopen, args.gext, args.savefile)
    if args.output == 'fa':
        ut.write_fasta(align1, align2, id1, id2, args.savefile)


if __name__ == '__main__':
    main()
