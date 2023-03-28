"""================================================================================================
Test script for DEDAL model. Place dedal folder in the same directory as this script.

Ben Iovino  02/21/23   VecAligns
================================================================================================"""

import os
import argparse
import logging
import tensorflow as tf
from dedal import infer  #pylint: disable=E0401
from utility import parse_fasta, write_align

logging.basicConfig(filename='dedal.log', filemode='w', format='%(asctime)s - %(levelname)s - %(message)s', level=logging.INFO)


def dedal(model, seq1, seq2):
    """=============================================================================================
    Runs the DEDAL model to get a pairwise alignment between two proteins.

    :param model: DEDAL model
    :param seq1: First protein sequence
    :param seq2: Second protein sequence
    :return: Alignment object
    ============================================================================================="""

    inputs = infer.preprocess(seq1, seq2)
    logging.info('DEDAL tokenization complete')
    align_out = model(inputs)
    logging.info('DEDAL alignment complete')
    output = infer.expand(
        [align_out['sw_scores'], align_out['paths'], align_out['sw_params']])
    output = infer.postprocess(output, len(seq1), len(seq2))
    alignment = infer.Alignment(seq1, seq2, *output)
    logging.info('DEDAL postprocessing complete')
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

    logging.info('DEDAL starting to run')

    # Load model and preprocess proteins
    dedal_model = tf.saved_model.load('dedal_3')
    logging.info('DEDAL model loaded')
    alignment = dedal(dedal_model, seq1, seq2)
    with open('dedal_output.txt', 'w', encoding='utf8') as f:
        f.write(str(alignment))

    # Parse alignment file, match to original sequences, and write to msf file
    tseq1, tseq2 = parse_align('dedal_output.txt')
    write_align(tseq1[1], tseq2[1], id1, id2, 'DEDAL', 'None', 'None', args.file1)
    os.remove('dedal_output.txt')


if __name__ == '__main__':
    main()
