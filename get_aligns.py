"""This script aligns all pairwise combinations of sequences in a reference folder
using the specified method.

__author__ = "Ben Iovino"
__date__ = "09/22/23"
"""

import logging
import os

log_filename = 'data/logs/get_aligns.log'  #pylint: disable=C0103
os.makedirs(os.path.dirname(log_filename), exist_ok=True)
logging.basicConfig(filename=log_filename, filemode='w',
                     level=logging.INFO, format='%(message)s')


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


def local_blosum(pw_aligns: list, ref: str, direc: str):
    """Writes all pairwise alignments to a file in the msf format

    :param pw_aligns: list of all pairwise combinations of sequences
    :param ref: reference folder
    :param direc: subfolder
    """

    method_direc = 'data/alignments/local_blosum'
    if not os.path.isdir(method_direc):
        os.makedirs(method_direc)
    if not os.path.isdir(f'{method_direc}/{ref}'):
        os.makedirs(f'{method_direc}/{ref}')
    if not os.path.isdir(f'{method_direc}/{ref}/{direc}'):
        os.makedirs(f'{method_direc}/{ref}/{direc}')

    # Align each pair of sequences
    for pair in pw_aligns:
        seq1 = pair[0]
        seq2 = pair[1]

        args = (f'--align local '
                f'--file1 data/sequences/{ref}/{direc}/{seq1} '
                f'--file2 data/sequences/{ref}/{direc}/{seq2} '
                f'--gopen -11 '
                f'--gext -1 '
                f'--matrix blosum '
                f'--score 62 '
                f'--savefile {method_direc}/{ref}/{direc}')
        os.system(f'python matrix.py {args}')


def global_blosum(pw_aligns: list, ref: str, direc: str):
    """Writes all pairwise alignments to a file in the msf format

    :param pw_aligns: list of all pairwise combinations of sequences
    :param ref: reference folder
    :param direc: subfolder
    """

    method_direc = 'data/alignments/global_blosum'
    if not os.path.isdir(method_direc):
        os.makedirs(method_direc)
    if not os.path.isdir(f'{method_direc}/{ref}'):
        os.makedirs(f'{method_direc}/{ref}')
    if not os.path.isdir(f'{method_direc}/{ref}/{direc}'):
        os.makedirs(f'{method_direc}/{ref}/{direc}')

    # Align each pair of sequences
    for pair in pw_aligns:
        seq1 = pair[0]
        seq2 = pair[1]

        args = (f'--align global '
                f'--file1 data/sequences/{ref}/{direc}/{seq1} '
                f'--file2 data/sequences/{ref}/{direc}/{seq2} '
                f'--gopen -11 '
                f'--gext -1 '
                f'--matrix blosum '
                f'--score 62 '
                f'--savefile {method_direc}/{ref}/{direc}')
        os.system(f'python matrix.py {args}')


def main():
    """Main
    """

    refs = os.listdir('data/BAliBASE_R1-5')
    for ref in refs:
        logging.info('Working on %s', ref)

        # Get all fasta files in reference folder
        seq_dir = os.listdir(f'data/sequences/{ref}')
        for direc in seq_dir:
            logging.info('Working on %s', direc)
            files = os.listdir(f'data/sequences/{ref}/{direc}')

            # Get pairwise alignments for each pair of sequences
            pw_aligns = get_aligns(files)
            local_blosum(pw_aligns, ref, direc)

        break



if __name__ == '__main__':
    main()
