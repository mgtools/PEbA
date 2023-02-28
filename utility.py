"""================================================================================================
This script contains utility functions to be imported into other scripts.

Ben Iovino  01/27/23   VecAligns
================================================================================================"""

import os
from Bio import SeqIO


def write_align(seq1, seq2, id1, id2, script, method, gopen, gext, path):
    """=============================================================================================
    This function accepts two sequences after gaps have been introduced and writes them to a file
    in MSF format, with some extra information about the alignment parameters.

    :param seq1: first aligned sequence
    :param seq2: second aligned sequence
    :param id1: first sequence id
    :param id2: second sequence id
    :param script: type of alignment performed
    :param method: scoring method used
    :param gopen: gap penalty for opening a new gap
    :param gext: gap penalty for extending a gap
    :param path: directory to write alignment to
    ============================================================================================="""

    # Length of alignment
    length = len(seq1)

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

    # Keep track of number of alignments in directory to name them
    path = '/'.join(path.split('/')[:-1])  # Remove fa file from path
    count = 0
    script = script.split('_')
    for file in os.listdir(path):
        if file.startswith('MATRIX') and script[1] == 'MATRIX':
            count += 1
        if file.startswith('PEbA_T5') and method == 'ProtT5_Sim':
            count += 1
        if file.startswith('PEbA_ESM2') and method == 'ESM2_Sim':
            count += 1
        if file.startswith('DEDAL') and script[1] == 'DEDAL':
            count += 1


    # Name file
    if script[1] == 'PEbA':
        if method == 'ProtT5_Sim':
            filename = f'{path}/{script[1]}_T5_{count}.msf'
        if method == 'ESM2_Sim':
            filename = f'{path}/{script[1]}_ESM2_{count}.msf'
    else:
        filename = f'{path}/{script[1]}_{count}.msf'

    # Write to a new line for every index in the split list i.e. every 55 characters
    with open(filename, 'w', encoding='utf8') as file:
        file.write('PileUp\n\n\n\n')
        file.write(f'   MSF: {length}  Type: P  Algorithm: {script[0]}  Method: {method}  Gopen: {gopen}  Gext: {gext}\n\n')
        file.write(f' Name: {id1} oo  Len:  {length}\n')
        file.write(f' Name: {id2} oo  Len:  {length}\n\n//\n\n\n\n')
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


def parse_matrix(file):
    """=============================================================================================
    This function takes a text file containing a substitution matrix and parses it into a dict.

    :param file: text file containing substitution matrix
    :return subs_matrix: dictionary of substitution matrix values
    ============================================================================================="""

    # Open file and read lines
    with open(file, 'r', encoding='utf8') as f:
        lines = f.readlines()

    # Parse first line to get amino acids
    amino_acids = lines[0].split()

    # Parse remaining lines to get substitution scores
    subs_matrix = {}
    for line in lines[1:]:
        line = line.split()
        for j, score in enumerate(line[1:]):
            subs_matrix[f'{line[0]}{amino_acids[j]}'] = int(score)
            subs_matrix[f'{amino_acids[j]}{line[0]}'] = int(score)

    return subs_matrix
