"""================================================================================================
Test script for DEDAL model. Place dedal folder in the same directory as this script.

Ben Iovino  02/21/23   VecAligns
================================================================================================"""

import tensorflow as tf
import tensorflow_hub as hub
import argparse
from dedal import infer  #pylint: disable=E0401
from utility import parse_fasta


def dedal(model, seq1, seq2):
    """=============================================================================================
    Runs the DEDAL model to get a pairwise alignment between two proteins.
    ============================================================================================="""

    inputs = infer.preprocess(seq1, seq2)
    align_out = model(inputs)
    output = infer.expand(
        [align_out['sw_scores'], align_out['paths'], align_out['sw_params']])
    output = infer.postprocess(output, len(seq1), len(seq2))
    alignment = infer.Alignment(seq1, seq2, *output)
    return alignment


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

    # Load model and preprocess proteins
    dedal_model = tf.saved_model.load('dedal_3')
    alignment = dedal(dedal_model, seq1, seq2)
    with open('dedal_output2.txt', 'w', encoding='utf8') as f:
        f.write(str(alignment))


if __name__ == '__main__':
    main()
