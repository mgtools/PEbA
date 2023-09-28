"""Contains functions used in the main alignment scripts.

__author__ = "Ben Iovino"
__date__ = 09/19/23
"""

import numpy as np
from Bio import SeqIO


def write_msf(seq1: str, seq2: str, id1: str, id2: str, method: str,
               gopen: float, gext: float, path: str, beg: list, end: list):
    """Writes alignment to file in msf format

    :param seq1: first aligned sequence
    :param seq2: second aligned sequence
    :param id1: first sequence id
    :param id2: second sequence id
    :param method: scoring method used
    :param gopen: gap penalty for opening a new gap
    :param gext: gap penalty for extending a gap
    :param path: directory to write alignment to
    :param list: beginning positions of seqs in alignment
    :param list: end positions of seqs in alignment
    """

    fpath = f'{path}/{id1}-{id2}.msf'

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
    alignment += f' Name: {id1} oo  Len:  {length}  Start/End:  {beg[0]},{end[0]}\n'
    alignment += f' Name: {id2} oo  Len:  {length}  Start/End:  {beg[1]},{end[1]}\n\n//\n\n\n\n'
    for i in range(len(seq1_split)): # pylint: disable=C0200
        alignment += f'{id1}      {seq1_split[i]}\n'
        alignment += f'{id2}      {seq2_split[i]}\n\n'

    # If no path is determined then print to console, otherwise write to file
    if path == 'n':
        print(alignment)
    else:
        with open(fpath, 'w', encoding='utf8') as file:
            file.write(alignment)


def write_fasta(seq1: str, seq2: str, id1: str, id2: str, path: str):
    """Writes alignment to file in fasta format

    :param seq1: first aligned sequence
    :param seq2: second aligned sequence
    :param id1: first sequence id
    :param id2: second sequence id
    :param path: directory to write alignment to
    """

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


def parse_fasta(filename: str) -> tuple:
    """Returns sequence and id from fasta file

    :param filename: name of file
    return (str, str): sequence and id
    """

    seq = ''
    with open(filename, 'r', encoding='utf8') as file:
        for seq in SeqIO.parse(file, 'fasta'):
            name = seq.id
            seq = str(seq.seq)
    return seq, name


def parse_matrix(file: str) -> dict:
    """Returns substitution matrix as dictionary

    :param file: text file containing substitution matrix
    :return dict: dict where keys are pairs of amino acids and values are substitution scores
    """

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


def initialize_matrices(seq1: str, seq2: str, align: str, gopen: float, gext: float) -> tuple:
    """Returns scoring and traceback matrices initialized based on alignment type

    :param seq1: first sequence
    :param seq2: second sequence
    :param align: alignment type (global or local)
    :param gopen: gap open penalty for global matrix intialization
    :param gext: gap extension penalty for global matrix intialization
    return (np.ndarray, np.ndarray): scoring and traceback matrices initialized based on align
    """

    # Initialize scoring and traceback matrix based on sequence lengths
    row_length = len(seq1)+1
    col_length = len(seq2)+1
    score_m = np.full((row_length, col_length), 0)
    trace_m = np.full((row_length, col_length), 0)

    # Initialize scoring matrix based on alignment type
    if align == 'global':
        for i in range(1, len(score_m[0])):
            score_m[0][i] = gopen+gext*i+1  # +1 to offset i starting at 1
            trace_m[0][i] = -1
        for i in range(1, len(score_m.T[0])):
            score_m.T[0][i] = gopen+gext*i+1
            trace_m.T[0][i] = 1

    return score_m, trace_m


def gap_penalty(gap: bool, horizontal: int, vertical: int, gopen: float, gext: float) -> tuple:
    """Returns score with gap penalty applied

    :param gap: boolean indicating if there is a gap
    :param horizontal: score if gap is horizontal
    :param vertical: score if gap is vertical
    :param gopen: gap open penalty
    :param gext: gap extension penalty
    return (float, float): score with gap penalty applied
    """

    if gap is False:  # Apply gap open penalty if there is no gap
        horizontal += gopen
        vertical += gopen
    if gap is True:  # Apply gap extension penalty if there is a gap
        horizontal += gext
        vertical += gext

    return horizontal, vertical


def assign_score(score: int, diagonal: int, horizontal: int,
                  vertical: int, trace_m: np.ndarray, gap: bool, i: int, j: int) -> tuple:
    """Returns traceback matrix with value based on highest score and gap status.

    :param score: score of current residue comparison
    :param diagonal: score of diagonal cell
    :param horizontal: score of horizontal cell
    :param vertical: score of vertical cell
    :param trace_m: traceback matrix
    :param gap: boolean indicating if there is a gap
    :param i: row index of current cell
    :param j: column index of current cell
    return (np.ndarray, bool): traceback matrix with value based on highest score and gap status
    """

    if score == diagonal:
        trace_m[i+1][j+1] = 0
        gap = False
    if score == horizontal:
        trace_m[i+1][j+1] = -1
        gap = True
    if score == vertical:
        trace_m[i+1][j+1] = 1
        gap = True

    return trace_m, gap


def global_traceback(trace_m: np.ndarray, seq1: str, seq2: str) -> tuple:
    """Returns global alignment of two sequences with gaps inserted based on traceback matrix

    :param score_m: scoring matrix
    :param trace_m: traceback matrix
    :param seq1: first sequence
    :param seq2: second sequence
    return (str, str): seq1 with gaps, seq2 with gaps
    """

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
    align1 = ''.join(rev_seq1)[::-1]
    align2 = ''.join(rev_seq2)[::-1]

    return align1, align2


def get_align(seq1: list, seq2: list, index: list) -> tuple:
    """Returns aligned sequences based on index

    :param seq1: first sequence
    :param seq2: second sequence
    :param index: final index of traceback matrix
    return (str, str): aligned sequences with gaps inserted
    """

    # Join and reverse lists to get sequences with gaps
    seq1 = ''.join(seq1)[::-1]
    seq2 = ''.join(seq2)[::-1]

    # Index sequences so they begin at final cell of traceback
    align1 = seq1[index[0]:]  # pylint: disable=E1127
    align2 = seq2[index[1]:]

    # Sometimes the initial cell of traceback is incorrect if there are two equally scoring cells
    # In this case, we will truncate the sequences so they begin at the first aligned residue
    if len(align1) != len(align2):
        if len(align1) > len(align2):
            align1 = align1[0:len(align2)]
        else:
            align2 = align2[0:len(align1)]

    return align1, align2


def local_traceback(score_m: np.ndarray, trace_m: np.ndarray, seq1: str, seq2: str) -> tuple:
    """Returns local alignment of two sequences with gaps inserted based on traceback matrix

    :param score_m: scoring matrix
    :param trace_m: traceback matrix
    :param seq1: first sequence
    :param seq2: second sequence
    return (str, str, list, list): seq1 with gaps, seq2 with gaps, index of highest score, and
        index of final cell of traceback matrix
    """

    # Find index of highest score in scoring matrix, start traceback at this matrix
    high_score_ind = np.unravel_index(np.argmax(score_m, axis=None), score_m.shape)
    seq1 = seq1[0:high_score_ind[0]+1]
    seq2 = seq2[0:high_score_ind[1]+1]

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
    score = score_m[index[0], index[1]]
    while score != 0:  # Stop traceback when we reach cell in score_m with 0
        val = trace_m[index[0], index[1]]  # Get value of cell in traceback matrix

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
        score = score_m[index[0], index[1]]  # Get score of next cell
    align1, align2 = get_align(rev_seq1, rev_seq2, index)

    return align1, align2, index, high_score_ind
