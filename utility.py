"""================================================================================================
This script contains utility functions to be imported into other scripts.

Ben Iovino  01/27/23   VecAligns
================================================================================================"""

import os
from datetime import datetime
from Bio import SeqIO


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

    # Create directory for file if it does not exist
    if not os.path.isdir('alignments'):
        os.makedirs('alignments')

    # Find max length sequence and write to file based on its length
    name1 = 'seque1'
    name2 = 'seque2'
    with open(f'alignments/{datetime.now()}.txt', 'w', encoding='utf8') as file:
        file.write('PileUp\n\n\n')
        file.write(f'   MSF:  {len(seq1)}  Type:  P\n\n')
        file.write(f' Name: {name1} oo  Len:  {len(seq1)}\n')
        file.write(f' Name: {name2} oo  Len:  {len(seq2)}\n\n//\n\n\n\n')
        for i in range(len(seq1_split)):
            file.write(f'{name1}      {seq1_split[i]}\n')
            file.write(f'{name2}      {seq2_split[i]}\n\n')


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
