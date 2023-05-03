"""================================================================================================
This script contains utility functions to be imported into other scripts.

Ben Iovino  01/27/23   VecAligns
================================================================================================"""

from Bio import SeqIO


def write_msf(seq1, seq2, id1, id2, method, gopen, gext, path):
    """=============================================================================================
    This function accepts two sequences after gaps have been introduced and writes them to a file
    in MSF format, with some extra information about the alignment parameters.

    :param seq1: first aligned sequence
    :param seq2: second aligned sequence
    :param id1: first sequence id
    :param id2: second sequence id
    :param method: scoring method used
    :param gopen: gap penalty for opening a new gap
    :param gext: gap penalty for extending a gap
    :param path: directory to write alignment to
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

    # Put alignment in string format
    length = len(seq1)
    alignment = 'PileUp\n\n\n\n'
    alignment += f'   MSF: {length}  Type: P  Method: {method}  Gopen: {gopen}  Gext: {gext}\n\n'
    alignment += f' Name: {id1} oo  Len:  {length}\n'
    alignment += f' Name: {id2} oo  Len:  {length}\n\n//\n\n\n\n'
    for i in range(len(seq1_split)): # pylint: disable=C0200
        alignment += f'{id1}      {seq1_split[i]}\n'
        alignment += f'{id2}      {seq2_split[i]}\n\n'

    # If no path is determined then print to console, otherwise write to file
    if path == 'n':
        print(alignment)
    else:
        with open(path, 'w', encoding='utf8') as file:
            file.write(alignment)


def write_fasta(seq1, seq2, id1, id2, path):
    """=============================================================================================
    This function accepts two sequences after gaps have been introduced and writes them to a file
    in fasta format.

    :param seq1: first aligned sequence
    :param seq2: second aligned sequence
    :param id1: first sequence id
    :param id2: second sequence id
    :param path: directory to write alignment to
    ============================================================================================="""

    # Replace '.' characters with '-' characters
    seq1 = seq1.replace('.', '-')
    seq2 = seq2.replace('.', '-')

    # Split sequences every 50 characters
    seq1_split = [seq1[i:i+50] for i in range(0, len(seq1), 50)]
    seq2_split = [seq2[i:i+50] for i in range(0, len(seq2), 50)]

    # Put alignment in string format
    alignment = f'>{id1}\n'
    for i in range(len(seq1_split)):  # pylint: disable=C0200
        alignment += f'{seq1_split[i]}\n'
    alignment += f'>{id2}\n'
    for i in range(len(seq2_split)):  # pylint: disable=C0200
        alignment += f'{seq2_split[i]}\n'

    # If no path is determined then print to console, otherwise write to file
    if path == 'n':
        print(alignment)
    else:
        with open(path, 'w', encoding='utf8') as file:
            file.write(alignment)


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
