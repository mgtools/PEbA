"""================================================================================================
This script contains utility functions to be imported into other scripts.

Ben Iovino  01/27/23   VecAligns
================================================================================================"""

import os
from datetime import datetime
from Bio import SeqIO


def write_align(seq1, seq2, id1, id2, script, matrix, gopen, gext):
    """=============================================================================================
    This function accepts two sequences after gaps have been introduced and writes them to a file
    in MSF format, with some extra information about the alignment parameters.

    :param seq1: first aligned sequence
    :param seq2: second aligned sequence
    :param id1: first sequence id
    :param id2: second sequence id
    :param script: type of alignment performed
    :param matrix: scoring matrix used
    :param gopen: gap penalty for opening a new gap
    :param gext: gap penalty for extending a gap
    ============================================================================================="""

    # Add space every 10 characters
    seq1 = [seq1[i:i+10] for i in range(0, len(seq1), 10)]
    seq1 = ' '.join(seq1)
    seq2 = [seq2[i:i+10] for i in range(0, len(seq2), 10)]
    seq2 = ' '.join(seq2)

    # Split sequences every 50 characters
    seq1_split = [seq1[i:i+55] for i in range(0, len(seq1), 55)]
    seq2_split = [seq2[i:i+55] for i in range(0, len(seq2), 55)]

    # Add extra spaces to either id if they are not the same length
    if len(id1) != len(id2):
        if len(id1) > len(id2):
            id2 = id2 + ' ' * (len(id1) - len(id2))
        else:
            id1 = id1 + ' ' * (len(id2) - len(id1))

    # Write to a new line for every index in the split list i.e. every 55 characters
    if not os.path.isdir('alignments'):
        os.makedirs('alignments')
    with open(f'alignments/{script}-{datetime.now()}.msf', 'w', encoding='utf8') as file:
        file.write('PileUp\n\n\n\n')
        file.write(f'   MSF:  {len(seq1)}  Type:  P  Matrix:  {matrix}  Gap Open:  {gopen}  Gap Ext:  {gext}\n\n')
        file.write(f' Name: {id1} oo  Len:  {len(seq1)}\n')
        file.write(f' Name: {id2} oo  Len:  {len(seq2)}\n\n//\n\n\n\n')
        for i in range(len(seq1_split)):  # pylint: disable=C0200
            file.write(f'{id1}      {seq1_split[i]}\n')
            file.write(f'{id2}      {seq2_split[i]}\n\n')


def parse_fasta(filename):
    """=============================================================================================
    This function accepts a fasta file name and returns the sequence and its ids.

    :param filename: name of file
    return: sequence and id
    ============================================================================================="""

    seq = ''
    with open(filename, 'r', encoding='utf8') as file:
        for seq in SeqIO.parse(file, 'fasta'):
            name = seq.id
            seq = str(seq.seq)
    return seq, name
