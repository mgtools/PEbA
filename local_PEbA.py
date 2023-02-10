"""================================================================================================
This script takes two protein sequences of varying length and finds the highest scoring local
alignment between the two.

Ben Iovino  01/23/23   VecAligns
================================================================================================"""

import os
import argparse
import torch
import numpy as np
from transformers import T5EncoderModel, T5Tokenizer
from utility import parse_fasta, write_align


def embed_seq(seq, tokenizer, encoder):
    """=============================================================================================
    This function accepts a protein sequence and returns a list of vectors, each vector representing
    a single amino acid.

    :param seq: protein sequence
    :param: tokenizer: tokenizer model
    :param encoder: encoder model
    return: list of vectors
    ============================================================================================="""

    # Add space after each amino acid so each residue is vectorized
    seq = [' '.join([*seq])]

    # Tokenize, encode, and load sequence
    ids = tokenizer.batch_encode_plus(seq, add_special_tokens=True, padding=True)
    input_ids = torch.tensor(ids['input_ids'])  # pylint: disable=E1101
    attention_mask = torch.tensor(ids['attention_mask'])  # pylint: disable=E1101

    # Extract sequence features
    with torch.no_grad():
        embedding = encoder(input_ids=input_ids,attention_mask=attention_mask)
    embedding = embedding.last_hidden_state.cpu().numpy()

    # Remove padding and special tokens
    features = []
    for seq_num in range(len(embedding)):  # pylint: disable=C0200
        seq_len = (attention_mask[seq_num] == 1).sum()
        seq_emd = embedding[seq_num][:seq_len-1]
        features.append(seq_emd)
    return features[0]


def local_align(seq1, seq2, vecs1, vecs2, gopen, gext):
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

    # Initialize scoring and traceback matrix based on sequence lengths
    row_length = len(seq1)+1
    col_length = len(seq2)+1
    score_m = np.full((row_length, col_length), 0)
    trace_m = np.full((row_length, col_length), 0)

    # Score matrix by moving through each index
    gap = False
    for i in range(len(seq1)):
        seq1_vec = vecs1[i]  # Corresponding amino acid vector in 1st sequence
        for j in range(len(seq2)):

            # Preceding scoring matrix values
            diagonal = score_m[i][j]
            horizontal = score_m[i+1][j]
            vertical = score_m[i][j+1]

            # Score pair of residues based off BLOSUM matrix
            seq2_vec = vecs2[j]  # Corresponding amino acid vector in 2nd sequence
            cos_sim = np.dot(seq1_vec,seq2_vec)/(np.linalg.norm(seq1_vec)*np.linalg.norm(seq2_vec))
            cos_sim = (cos_sim*10)

            # Add to matrix values via scoring method
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


def traceback(score_m, trace_m, seq1, seq2):
    """=============================================================================================
    This function accepts a scoring and a traceback matrix and two sequences and returns the highest
    scoring local alignment between the two sequences

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
    align1 = seq1+"."*max(0, len(seq2)-len(seq1))
    align2 = seq2+"."*max(0, len(seq1)-len(seq2))
    return align1, align2


def main():
    """=============================================================================================
    This function initializes two protein sequences, calls embed_seq() to vectorize each sequence,
    calls local_align() to obtain the scoring and traceback matrix from SW alignment (with vector
    similarity in the scoring system), calls traceback() to get the local alignment, and
    then write_align() to write the alignment to a file in MSF format.
    ============================================================================================="""

    parser = argparse.ArgumentParser()
    parser.add_argument('-file1', type=str, default='./test1.fa', help='Name of first fasta file')
    parser.add_argument('-file2', type=str, default='./test2.fa', help='Name of second fasta file')
    parser.add_argument('-gopen', type=int, default=-11, help='Penalty for opening a gap')
    parser.add_argument('-gext', type=int, default=-1, help='Penalty for extending a gap')
    args = parser.parse_args()

    # Parse fasta files for sequences and ids
    seq1, id1 = parse_fasta(args.file1)
    seq2, id2 = parse_fasta(args.file2)

    # Load model tokenizer and encoder models
    if os.path.exists('tok.pt'):
        tokenizer = torch.load('tok.pt')
    else:
        tokenizer = T5Tokenizer.from_pretrained("Rostlab/prot_t5_xl_uniref50", do_lower_case=False)
        torch.save(tokenizer, 'tok.pt')
    if os.path.exists('prottrans.pt'):
        model = torch.load('prottrans.pt')
    else:
        model = T5EncoderModel.from_pretrained("Rostlab/prot_t5_xl_uniref50")
        torch.save(model, 'prottrans.pt')

    # Vectorize sequences
    vecs1 = embed_seq(seq1, tokenizer, model)
    vecs2 = embed_seq(seq2, tokenizer, model)

    # Call local_align() to get scoring and traceback matrix
    score_m, trace_m = local_align(seq1, seq2, vecs1, vecs2, args.gopen, args.gext)

    # Get highest scoring local alignment between seq1 and seq2 and write to file
    align1, align2 = traceback(score_m, trace_m, seq1, seq2)
    write_align(align1, align2, id1, id2, 'local_PEbA', 'None', args.gopen, args.gext, args.file1)  #pylint: disable=E1121


if __name__ == '__main__':
    main()
