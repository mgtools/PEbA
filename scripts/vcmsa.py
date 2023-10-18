"""Calls vcmsa on two fasta files and writes the alignment in the desired format

__author__ = "Ben Iovino"
__date__ = 10/03/23
"""

import argparse
import os
import utility as ut
from Bio import SeqIO
from transformers import T5Tokenizer, T5Model


def get_model():
    """Saves prott5 model and tokenizer if not already saved
    """

    model_direc = 'data/prot_t5_xl_uniref50'
    if not os.path.exists(model_direc):
        os.mkdir(model_direc)
        tokenizer = T5Tokenizer.from_pretrained('Rostlab/prot_t5_xl_uniref50')
        model = T5Model.from_pretrained('Rostlab/prot_t5_xl_uniref50')
        model.save_pretrained(model_direc)
        tokenizer.save_pretrained(model_direc)


def write_fasta(fa1: str, fa2: str) -> tuple:
    """Writes input fasta sequences to one file and returns sequences and ids

    :param f1: fasta file 1
    :param f2: fasta file 2
    :return (str, str, str, str): sequences and ids
    """

    if not os.path.exists('data/tmp'):
        os.mkdir('data/tmp')

    # Get fasta sequences and ids from file
    for record in SeqIO.parse(fa1, 'fasta'):
        seq1 = str(record.seq)
        id1 = record.id
    for record in SeqIO.parse(fa2, 'fasta'):
        seq2 = str(record.seq)
        id2 = record.id

    # Write to file
    with open('data/tmp/seqs.fa', 'w', encoding='utf8') as file:
        file.write(f'>{id1}\n{seq1}\n>{id2}\n{seq2}')

    return seq1, seq2, id1, id2


def clean_align(align: str, id1: str, id2: str) -> tuple:
    """Returns aligned sequences from vcmsa output

    :param align: vcmsa output from file
    :param id1: fasta id 1
    :param id2: fasta id 2
    """

    align1, align2 = '', ''
    for line in align.split('\n'):
        if line.startswith(id1):
            align1 += line.split()[1].replace('-', '.')
        elif line.startswith(id2):
            align2 += line.split()[1].replace('-', '.')

    return align1, align2


def get_pos(align: str, seq: str) -> tuple:
    """Returns beginning and end positions of an aligned sequence

    :param align: aligned sequence
    :param seq: full fasta sequence
    :return (int, int): beginning and end positions
    """

    # Get beginning and end indices with find()
    align = align.replace('.', '')
    beg = seq.find(align)
    end = beg + len(align)-1

    return beg, end


def align_seqs(seq1: str, seq2: str, id1: str, id2: str):
    """Calls vcmsa on a fasta file and writes the alignment

    :param seq1: fasta sequence 1
    :param seq2: fasta sequence 2
    :param id1: fasta id 1
    :param id2: fasta id 2
    :return (str, str, list, list): aligned sequences, beg/end positions of each seq
    """

    os.system('vcmsa -i data/tmp/seqs.fa -o data/tmp/align -m data/prot_t5_xl_uniref50')

    # Read output file and get aligned sequences
    with open('data/tmp/align', 'r', encoding='utf8') as file:
        align = file.read()
    align1, align2 = clean_align(align, id1, id2)
    beg1, end1 = get_pos(align1, seq1)
    beg2, end2 = get_pos(align2, seq2)
    beg = [beg1, beg2]
    end = [end1, end2]

    return align1, align2, beg, end


def main():
    """Takes two pdb files and writes the alignment in the desired format
    """

    parser = argparse.ArgumentParser()
    parser.add_argument('-f1', type=str, default='data/sequences/RV11/BB11002/1abo_A.fa')
    parser.add_argument('-f2', type=str, default='data/sequences/RV11/BB11002/1bb9_.fa')
    parser.add_argument('-o', '--output', type=str, default='msf')
    parser.add_argument('-sf', '--savefile', type=str)
    args = parser.parse_args()

    get_model()
    seq1, seq2, id1, id2 = write_fasta(args.f1, args.f2)
    align1, align2, beg, end = align_seqs(seq1, seq2, id1, id2)

    # Write align based on desired output format
    if args.output == 'msf':
        ut.write_msf(align1, align2, id1, id2, 'vcmsa',
                      0, 0, args.savefile, beg, end)
    if args.output == 'fa':
        ut.write_fasta(align1, align2, id1, id2, args.savefile)

    os.system('rm -r data/tmp')


if __name__ == '__main__':
    main()
