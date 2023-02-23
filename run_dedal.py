"""================================================================================================
Test script for DEDAL model. Place dedal folder in the same directory as this script.

Ben Iovino  02/21/23   VecAligns
================================================================================================"""

import tensorflow as tf
import argparse
#from dedal import infer  #pylint: disable=E0401
from utility import parse_fasta


def dedal(model, seq1, seq2):
    """=============================================================================================
    Runs the DEDAL model to get a pairwise alignment between two proteins.

    :param model: DEDAL model
    :param seq1: First protein sequence
    :param seq2: Second protein sequence
    :return: Alignment object
    ============================================================================================="""

    inputs = infer.preprocess(seq1, seq2)
    align_out = model(inputs)
    output = infer.expand(
        [align_out['sw_scores'], align_out['paths'], align_out['sw_params']])
    output = infer.postprocess(output, len(seq1), len(seq2))
    alignment = infer.Alignment(seq1, seq2, *output)
    return alignment


def parse_align(file):
    """=============================================================================================
    This function gathers the truncated sequences and their beggining and ending indices from the
    alignment file.

    :param file: alignment file
    return: truncated sequences and their positions
    ============================================================================================="""

    # Gather beginning position, truncated sequence, and ending position
    tseq1 = [0, '', 0]
    tseq2 = [0, '', 0]
    with open(file, 'r', encoding='utf8') as f:
        count = 0
        for line in f:
            split_line = line.split()
            if count == 0:  # First line contains first sequence
                tseq1[0] = int(split_line[0])
                tseq1[1] = split_line[1].replace('-', '.')
                tseq1[2] = int(split_line[2])
            if count == 2:  # Third line contains second sequence
                tseq2[0] = int(split_line[0])
                tseq2[1] = split_line[1].replace('-', '.')
                tseq2[2] = int(split_line[2])
            count += 1

    return tseq1, tseq2


def match_seq(seq, tseq):
    """=============================================================================================
    This function adds any missing characters from an original protein sequence to the truncated 
    sequence so that they can be accurately compared.

    :param seq: original protein sequence
    :param trunc_seq: list with beginning position, truncated sequence, and ending position
    return: full sequence
    ============================================================================================="""

    full_seq = ''
    beg, end = tseq[0], tseq[2]
    beg_chars = seq[:beg]  # Characters before truncated sequence
    end_chars = seq[end:]  # Characters after truncated sequence
    full_seq = beg_chars + tseq[1] + end_chars

    return full_seq


def pad_seq(seq1, seq2):
    """=============================================================================================
    This function adds gaps to the beginning and end of either sequence so their lengths are equal.

    :param seq1: list with beginning position, full first sequence, and ending position
    :param seq2: list with beginning position, full second sequence, and ending position
    return: padded sequences
    ============================================================================================="""

    # Pad the beginning of the shorter sequence
    fseq1, fseq2 = seq1[1], seq2[1]
    beg1, beg2 = seq1[0], seq2[0]
    if beg1 < beg2:
        pad = beg2 - beg1  # Number of gaps to add
        fseq1 = '.' * pad + fseq1
    elif beg2 < beg1:
        pad = beg1 - beg2
        fseq2 = '.' * pad + fseq2

    # Pad the end of the shorter sequence
    end1, end2 = seq1[2], seq2[2]
    if end1 > end2:
        pad = end1 - end2
        fseq1 += '.' * pad
    elif end2 > end1:
        pad = end2 - end1
        fseq2 += '.' * pad

    print(fseq1)
    print()
    print(fseq2)


def main():
    """=============================================================================================
    Run the DEDAL model to get a pairwise alignment between two proteins.
    ============================================================================================="""

    parser = argparse.ArgumentParser()
    parser.add_argument('-file1', type=str, default='./test1.fa', help='Name of first fasta file')
    parser.add_argument('-file2', type=str, default='./test2.fa', help='Name of second fasta file')
    args = parser.parse_args()

    # Load fasta files and ids
    seq1, id1 = parse_fasta(args.file1)
    seq2, id2 = parse_fasta(args.file2)

    '''
    # Load model and preprocess proteins
    dedal_model = tf.saved_model.load('dedal_3')
    alignment = dedal(dedal_model, seq1, seq2)
    with open('dedal_output.txt', 'w', encoding='utf8') as f:
        f.write(str(alignment))
    '''

    # Parse alignment file, match to original sequences, and write to msf file
    tseq1, tseq2 = parse_align('dedal_output.txt')
    tseq1[1] = match_seq(seq1, tseq1)
    tseq2[1] = match_seq(seq2, tseq2)
    pad_seq(tseq1, tseq2)


if __name__ == '__main__':
    main()
