"""This script aligns all pairwise combinations of sequences in a reference folder
using the specified method.

__author__ = "Ben Iovino"
__date__ = "09/22/23"
"""

import argparse
import datetime
import os
import logging
#import tensorflow as tf
import utility as ut

log_filename = 'data/logs/get_aligns.log'  #pylint: disable=C0103
os.makedirs(os.path.dirname(log_filename), exist_ok=True)
logging.basicConfig(filename=log_filename, filemode='w',
                     level=logging.INFO, format='%(message)s')


def log_align(ref: str, direc: str, pair: str, method: str) -> bool:
    """Returns True if alignment needs to be done, False if not

    :param ref: reference folder
    :param direc: subfolder
    :param pair: pair of sequences
    :param method: alignment method
    :return bool: True or False
    """

    # Check if alignment already exists
    seq1, seq2 = pair[0].split('.')[0], pair[1].split('.')[0]
    if os.path.isfile(f'data/alignments/{method}/{ref}/{direc}/{seq1}-{seq2}.msf'):
        logging.info('%s\t%s/%s/%s-%s %s',
                      datetime.datetime.now(), ref, direc, seq1, seq2, 'already aligned')

        return False
    
    # If not, need to align
    logging.info('%s\t%s %s/%s/%s-%s',
                  datetime.datetime.now(), 'aligning', ref, direc, seq1, seq2)

    return True
    

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


def blosum(pw_aligns: list, ref: str, direc: str, method: str):
    """Writes all pairwise SW blosum alignments to a file in the msf format

    :param pw_aligns: list of all pairwise combinations of sequences
    :param ref: reference folder
    :param direc: subfolder
    :param method: alignment method
    """

    method_direc = f'data/alignments/{method}_blosum'
    if not os.path.isdir(method_direc):
        os.makedirs(method_direc)
    if not os.path.isdir(f'{method_direc}/{ref}'):
        os.makedirs(f'{method_direc}/{ref}')
    if not os.path.isdir(f'{method_direc}/{ref}/{direc}'):
        os.makedirs(f'{method_direc}/{ref}/{direc}')

    # Align each pair of sequences
    for pair in pw_aligns:
        if not log_align(ref, direc, pair, f'{method}_blosum'):  # Ignore if already aligned
            continue
        seq1, seq2 = pair[0], pair[1]

        args = (f'--align {method} '
                f'--file1 data/sequences/{ref}/{direc}/{seq1} '
                f'--file2 data/sequences/{ref}/{direc}/{seq2} '
                f'--gopen -11 '
                f'--gext -1 '
                f'--matrix blosum '
                f'--score 62 '
                f'--savefile {method_direc}/{ref}/{direc}')
        os.system(f'python scripts/matrix.py {args}')


def peba(pw_aligns: list, ref: str, direc: str, method: str):
    """Writes all pairwise SW peba alignments to a file in the msf format

    :param pw_aligns: list of all pairwise combinations of sequences
    :param ref: reference folder
    :param direc: subfolder
    :param method: alignment method
    """

    method_direc = f'data/alignments/{method}_peba'
    if not os.path.isdir(method_direc):
        os.makedirs(method_direc)
    if not os.path.isdir(f'{method_direc}/{ref}'):
        os.makedirs(f'{method_direc}/{ref}')
    if not os.path.isdir(f'{method_direc}/{ref}/{direc}'):
        os.makedirs(f'{method_direc}/{ref}/{direc}')

    # Align each pair of sequences
    for pair in pw_aligns:
        if not log_align(ref, direc, pair, f'{method}_peba'):  # Ignore if already aligned
            continue
        seq1, seq2 = pair[0], pair[1]
        embed1 = seq1.split('.')[0] + '.txt'
        embed2 = seq2.split('.')[0] + '.txt'

        args = (f'--align {method} '
                f'--file1 data/sequences/{ref}/{direc}/{seq1} '
                f'--file2 data/sequences/{ref}/{direc}/{seq2} '
                f'--embed1 data/embeddings/{ref}/{direc}/{embed1} '
                f'--embed2 data/embeddings/{ref}/{direc}/{embed2} '
                f'--gopen -11 '
                f'--gext -1 '
                f'--savefile {method_direc}/{ref}/{direc}')
        os.system(f'python scripts/peba.py {args}')


def dedal_run(pw_aligns: list, ref: str, direc: str, model):
    """Writes all pairwise dedal alignments to a file in the msf format

    :param pw_aligns: list of all pairwise combinations of sequences
    :param ref: reference folder
    :param direc: subfolder
    :param model: dedal model
    """

    from dedal import infer  #pylint: disable-all

    method_direc = 'data/alignments/dedal'
    if not os.path.isdir(method_direc):
        os.makedirs(method_direc)
    if not os.path.isdir(f'{method_direc}/{ref}'):
        os.makedirs(f'{method_direc}/{ref}')
    if not os.path.isdir(f'{method_direc}/{ref}/{direc}'):
        os.makedirs(f'{method_direc}/{ref}/{direc}')

    # Align each pair of sequences
    for pair in pw_aligns:
        if not log_align(ref, direc, pair, 'dedal'):  # Ignore if already aligned
            continue
        seq1, id1 = ut.parse_fasta(f'data/sequences/{ref}/{direc}/{pair[0]}')
        seq2, id2 = ut.parse_fasta(f'data/sequences/{ref}/{direc}/{pair[1]}')

        # Get alignment from dedal
        inputs = infer.preprocess(seq1, seq2)
        align_out = model(inputs)
        output = infer.expand(
            [align_out['sw_scores'], align_out['paths'], align_out['sw_params']])
        output = infer.postprocess(output, len(seq1), len(seq2))
        alignment = infer.Alignment(seq1, seq2, *output)

        # Want first and third lines of alignment
        align1 = str(alignment).split('\n')[0]
        align2 = str(alignment).split('\n')[2]
        beg = [align1.split()[0], align2.split()[0]]
        end = [align1.split()[2], align2.split()[2]]
        align1 = align1.split()[1].replace('-', '.')
        align2 = align2.split()[1].replace('-', '.')
        ut.write_msf(align1, align2, id1, id2, 'dedal',
                     0, 0, f'{method_direc}/{ref}/{direc}', beg, end)
        

def fatcat(pw_aligns: list, ref: str, direc: str):
    """Writes all pairwise fatcat alignments to a file in the desired format

    :param pw_aligns: list of all pairwise combinations of sequences
    :param ref: reference folder
    :param direc: subfolder
    :param method: alignment method
    """

    if ref != 'RV11':  # only works for RV11, rest don't have pdb files
        return
    if direc in ('BB11006', 'BB11037'):  # unable to get correct pdb files
        return

    method_direc = 'data/alignments/fatcat'
    if not os.path.isdir(method_direc):
        os.makedirs(method_direc)
    if not os.path.isdir(f'{method_direc}/{ref}'):
        os.makedirs(f'{method_direc}/{ref}')
    if not os.path.isdir(f'{method_direc}/{ref}/{direc}'):
        os.makedirs(f'{method_direc}/{ref}/{direc}')

    # Align each pair of sequences
    for pair in pw_aligns:
        if not log_align(ref, direc, pair, 'fatcat'):  # Ignore if already aligned
            continue
        seq1, seq2 = pair[0], pair[1]
        pdb1, pdb2 = seq1.split('.')[0], seq2.split('.')[0]
        args = (f'-f1 data/pdb/{direc}/{pdb1}.pdb '
                f'-f2 data/pdb/{direc}/{pdb2}.pdb '
                f'-sf {method_direc}/{ref}/{direc}')
        os.system(f'python scripts/fatcat.py {args}')


def vcmsa(pw_aligns: list, ref: str, direc: str):
    """Writes all pairwise vcmsa alignments to a file in the desired format
    """

    method_direc = 'data/alignments/vcmsa'
    if not os.path.isdir(method_direc):
        os.makedirs(method_direc)
    if not os.path.isdir(f'{method_direc}/{ref}'):
        os.makedirs(f'{method_direc}/{ref}')
    if not os.path.isdir(f'{method_direc}/{ref}/{direc}'):
        os.makedirs(f'{method_direc}/{ref}/{direc}')

    # Align each pair of sequences
    for pair in pw_aligns:
        if not log_align(ref, direc, pair, 'vcmsa'):  # Ignore if already aligned
            continue
        seq1, seq2 = pair[0], pair[1]
        args = (f'-f1 data/sequences/{ref}/{direc}/{seq1} '
                f'-f2 data/sequences/{ref}/{direc}/{seq2} '
                f'-sf {method_direc}/{ref}/{direc}')
        os.system(f'python scripts/vcmsa.py {args}')
               

def main():
    """Main
    """

    parser = argparse.ArgumentParser()
    parser.add_argument('-a', type=str, default='blosum')
    parser.add_argument('-m', type=str, default='global')
    args = parser.parse_args()

    # Don't want to load if not using
    if args.a == 'dedal':
        dedal_model = tf.saved_model.load('dedal_3')

    refs = os.listdir('data/BAliBASE_R1-5')
    for ref in refs:

        # Get all fasta files in reference folder
        seq_dir = os.listdir(f'data/sequences/{ref}')
        for direc in seq_dir:

            # Get pairwise alignments for each pair of sequences
            files = os.listdir(f'data/sequences/{ref}/{direc}')
            pw_aligns = get_aligns(files)
            if args.a == 'blosum':
                blosum(pw_aligns, ref, direc, args.m)
            elif args.a == 'peba':
                peba(pw_aligns, ref, direc, args.m)
            elif args.a == 'dedal':
                dedal_run(pw_aligns, ref, direc, dedal_model)
            elif args.a == 'fatcat':
                fatcat(pw_aligns, ref, direc)
            elif args.a == 'vcmsa':
                vcmsa(pw_aligns, ref, direc)


if __name__ == '__main__':
    main()
