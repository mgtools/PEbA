"""================================================================================================
This script takes two protein sequences of varying length and finds the highest scoring global
alignment between the two. Instead of using typical scoring matrices, such as BLOSUM, each amino
acid is embedded into a vector and the cosine similarity between vectors is used to determine if
two residues should be aligned.

Ben Iovino  01/26/23   VecAligns
================================================================================================"""

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


def global_align(seq1, seq2, vecs1, vecs2, gopen, gext):
    """=============================================================================================
    This function accepts two sequences, creates a matrix corresponding to their lengths, and
    calculates the score of the alignments for each index. A second matrix is scored so that the
    best alignment can be tracebacked.

    :param seq1: first sequence
    :param seq2: second sequence
    :param vecs1: list of vectors for first sequence
    :param vecs2: list of vectors for second sequence
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
    for i in range(len(seq1)):
        seq1_vec = vecs1[i]  # Corresponding amino acid vector
        for j in range(len(seq2)):

            # Preceding scoring matrix values
            diagonal = score_m[i][j]
            horizontal = score_m[i+1][j]
            vertical = score_m[i][j+1]

            # Score residues based off cosine similarity between vectors
            seq2_vec = vecs2[j]  # Corresponding amino acid vector
            cos_sim = np.dot(seq1_vec,seq2_vec)/(np.linalg.norm(seq1_vec)*np.linalg.norm(seq2_vec))
            cos_sim = (cos_sim*10)

            # Add to matrix values via scoring method
            diagonal += cos_sim
            if gap is False:  # Apply gap_open penalty if there is no gap
                horizontal += gopen
                vertical += gopen
            if gap is True:  # Apply gap_extension penalty if there is a gap
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
    This function initializes two protein sequences, calls embed_seq() to vectorize each sequence,
    calls global_align() to obtain the scoring and traceback matrix from NW alignment (with vector
    similarity in the scoring system), calls traceback() to get the global alignment, and
    then write_align() to write the alignment to a file in MSF format.
    ============================================================================================="""

    parser = argparse.ArgumentParser()
    parser.add_argument('-file1', type=str, default='sequences/BB30004_1.fa', help='Name of first fasta file')
    parser.add_argument('-file2', type=str, default='sequences/BB30004_2.fa', help='Name of second fasta file')
    parser.add_argument('-gopen', type=int, default=-11, help='Penalty for opening a gap')
    parser.add_argument('-gext', type=int, default=-1, help='Penalty for extending a gap')
    args = parser.parse_args()

    # Parse fasta files for sequences and ids
    seq1, id1 = parse_fasta(args.file1)
    seq2, id2 = parse_fasta(args.file2)

    # Load model tokenizer and encoder models
    tokenizer = T5Tokenizer.from_pretrained("Rostlab/prot_t5_xl_uniref50", do_lower_case=False)
    model = T5EncoderModel.from_pretrained("Rostlab/prot_t5_xl_uniref50")

    # Vectorize sequences
    vecs1 = embed_seq(seq1, tokenizer, model)
    vecs2 = embed_seq(seq2, tokenizer, model)

    # Call global_align() to get traceback matrix
    trace_m = global_align(seq1, seq2, vecs1, vecs2, args.gopen, args.gext)

    # Get global alignment between seq1 and seq2 and write to file
    align1, align2 = traceback(trace_m, seq1, seq2)
    write_align(align1, align2, id1, id2, 'PEbA')


if __name__ == '__main__':
    main()
