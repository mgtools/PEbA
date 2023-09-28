"""This script parses individual sequences from the BALIBASE database and writes them to a single
file for the corresponding reference. It also parses the pairwise alignments from each MSA and
writes each one to its own file.

__author__ = "Ben Iovino"
__date__ = 09/21/23
"""

import os
from Bio import SeqIO


def parse_ref_folder(path: str) -> tuple:
    """Returns two lists, one containing the paths to all MSF files and the other containing the
    paths to all FASTA files.

    :param path: path to reference folder
    :return (list, list): list of msf files, list of fasta files
    """

    # Get MSF and FASTA files
    msf_files = []
    fasta_files = []
    for file in os.listdir(path):
        if file.endswith('.msf'):  # Append to msf list
            msf_files.append(f'{path}/{file}')
        if file.endswith('.tfa'):  # Append to fasta list
            fasta_files.append(f'{path}/{file}')
        if file.endswith('.in_tfa'):
            new_file = file.split('.in_tfa')[0] + '.tfa'
            os.rename(f'{path}/{file}',f'{path}/{new_file}')
            fasta_files.append(f'{path}/{new_file}')

    return msf_files, fasta_files


def parse_fasta(filename: str, bb_dir: str) -> list:
    """Writes each sequence in a fasta file to its own file

    :param filename: name of file
    :param bb_dir: directory to place parsed files
    return list: list of sequences
    """

    # Get reference folder name and folder for the correpsonding fasta files
    famname = filename.split('/')[-1]  # Get fasta file name
    famname = famname.split('.tfa')[0]  # Remove file extension
    if not os.path.isdir(f'{bb_dir}/{famname}'):
        os.makedirs(f'{bb_dir}/{famname}')

    # Parse fasta file and write each sequence to its own file in the corresponding folder
    seq = ''
    seqs = []
    with open(filename, 'r', encoding='utf8') as file:
        for seq in SeqIO.parse(file, 'fasta'):
            SeqIO.write(seq, f'{bb_dir}/{famname}/{seq.id}.fa', 'fasta')
            seqs.append(f'{seq.id}.fa')

    return seqs


def parse_msf(filename: str, id1: str, id2: str) -> tuple:
    """Returns the pairwise alignment of two sequences from an MSF file

    :param filename: name of file
    :param id1: first sequence id
    :param id2: second sequence id
    return (str, str): align1, align2 - corresponding pairwise alignments
    """

    seq1 = []
    seq2 = []
    with open(filename, 'r', encoding='utf8') as file:
        for line in file:  # Get each sequence
            if line.startswith(id1):
                seq1.append(''.join(line.split()[1:]))
            elif line.startswith(id2):
                seq2.append(''.join(line.split()[1:]))

    # Go through both sequences and remove positions with gaps in both
    seq1 = list(''.join(seq1))
    seq2 = list(''.join(seq2))
    for i in range(len(seq1)):  # pylint: disable=C0200
        if seq1[i] == '.' and seq2[i] == '.':
            seq1[i] = ''
            seq2[i] = ''
    align1 = ''.join(seq1)
    align2 = ''.join(seq2)

    return align1, align2


def write_align(seq1: str, seq2: str, id1: str, id2: str, path: str):
    """Writes two alignments to a file in the msf format

    :param seq1: first aligned sequence
    :param seq2: second aligned sequence
    :param id1: first sequence id
    :param id2: second sequence id
    :param path: directory
    """

    path += f'{id1}-{id2}.msf'
    alength = len(seq1)  # Length of alignment
    length1, length2 = len(seq1.strip('.')), len(seq2.strip('.'))  # Length of seqs

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
    with open(f'{path}', 'w', encoding='utf8') as file:
        file.write('PileUp\n\n\n\n')
        file.write(f'   MSF:  {alength}  Type:  P\n\n')
        file.write(f' Name: {id1} oo  Len:  {alength}  Start/End:  {0},{length1}\n')
        file.write(f' Name: {id2} oo  Len:  {alength}  Start/End:  {0},{length2}\n\n//\n\n\n\n')
        for i in range(len(seq1_split)):  # pylint: disable=C0200
            file.write(f'{id1}      {seq1_split[i]}\n')
            file.write(f'{id2}      {seq2_split[i]}\n\n')


def get_aligns(seqs):
    """Returns a list of all pairs in a list of sequences

    :param pairs: list of sequences 
    """

    # Get all pairwise alignments from the fasta files correpsonding to the MSF file
    pw_aligns = []
    for i, seq in enumerate(seqs):
        loop_count = i  # Keep track of number of loops so no repeats occur
        while loop_count != len(seqs):
            if seq != seqs[loop_count]:  # Don't want to align a sequence to itself
                pw_aligns.append(sorted([seq, seqs[loop_count]]))
            loop_count+=1

    return pw_aligns


def parse_align_files(msf_files: list, fasta_files: list, seq_dir: str, align_dir: str):
    """Writes fasta files to a single file for each reference alignment and parses the pairwise
    alignments from each MSA and writes each one to its own file.

    :param msf_files: list of msf files
    :param fasta_files: list of fasta files
    :param seq_dir: directory to place fasta files
    :param align_dir: directory to place pairwise alignments
    """

    # Parse each fasta file, store names of each for subsequent msf parsing
    seqs = []
    for file in fasta_files:
        new_seqs = parse_fasta(file, seq_dir)
        seqs.append(new_seqs)  # Store in nested list to access only relevant fa files for each msf

    # Each MSF files corresponds to a set of fasta files
    for i, file in enumerate(msf_files):
        ref_align = file.rsplit('/', maxsplit=1)[-1].strip('.msf')  # Get name of ref alignment
        pw_aligns = get_aligns(seqs[i])

        # For the selected pairs, get PW alignment from ref MSA
        for pair in pw_aligns:
            seq1 = pair[0]
            seq2 = pair[1]

            # Grab pairwise alignment from reference MSA
            seq1, seq2 = seq1.split('.')[0], seq2.split('.')[0]  # Remove fa
            align1, align2 = parse_msf(file, seq1, seq2)  # Gather pairwise alignment
            dir_path = f'{align_dir}/{ref_align}'
            if not os.path.isdir(dir_path):
                os.makedirs(dir_path)
            file_path = f'{align_dir}/{ref_align}/'
            write_align(align1, align2, seq1, seq2, file_path)  # Write pairwise alignment


def main():
    """Main
    """

    refs = os.listdir('data/BAliBASE_R1-5')
    for ref in refs:

        # Create directory for sequences and alignments
        ref_dir = f'data/BAliBASE_R1-5/{ref}'
        seq_dir = f'data/sequences/{ref}'
        align_dir = f'data/alignments/refs/{ref}'

        # Get all msf and fasta files in ref, sort, and parse/align
        msf_files, fasta_files = parse_ref_folder(ref_dir)
        msf_files.sort()
        fasta_files.sort()
        parse_align_files(msf_files, fasta_files, seq_dir, align_dir)


if __name__ == '__main__':
    main()
