"""================================================================================================
This script embes protein sequences using a pretrained T5 model and outputs a numpy array of the
embedded sequences to a file.

Ben Iovino  02/16/23  VecAligns
================================================================================================"""

import os
import re
import torch
import numpy as np
from transformers import T5EncoderModel, T5Tokenizer
from Bio import SeqIO


def parse_ref_folder(path):
    """=============================================================================================
    This function accepts a folder name parses every MSF and FASTA file to write each pairwise
    alignment and fasta sequence to their own files.

    :param path: directory path to folder
    :return: list of fasta files
    ============================================================================================="""

    # Get FASTA files
    fasta_files = []
    for file in os.listdir(path):
        if file.endswith('.tfa'):  # Append to fasta list
            fasta_files.append(f'{path}/{file}')
    return fasta_files


def embed_seq(seq, tokenizer, encoder):
    """=============================================================================================
    This function accepts a protein sequence and returns a list of vectors, each vector representing
    a single amino acid.

    :param seq: protein sequence
    :param: tokenizer: tokenizer model
    :param encoder: encoder model
    return: list of vectors
    ============================================================================================="""

    # Remove special chars, add space after each amino acid so each residue is vectorized
    seq = re.sub(r"[UZOB]", "X", seq)
    seq = [' '.join([*seq])]

    # Tokenize, encode, and load sequence
    ids = tokenizer.batch_encode_plus(seq, add_special_tokens=True, padding=True)
    input_ids = torch.tensor(ids['input_ids'])  # pylint: disable=E1101
    attention_mask = torch.tensor(ids['attention_mask'])  # pylint: disable=E1101

    # Extract sequence features
    with torch.no_grad():
        embedding = encoder(input_ids=input_ids,attention_mask=attention_mask)
    embedding = embedding.last_hidden_state.cpu().numpy()

    # Remove padding and special tokens
    features = []
    for seq_num in range(len(embedding)):  # pylint: disable=C0200
        seq_len = (attention_mask[seq_num] == 1).sum()
        seq_emd = embedding[seq_num][:seq_len-1]
        features.append(seq_emd)
    return features[0]


def parse_fasta(filename, tokenizer, model):
    """=============================================================================================
    This function accepts a fasta file with multiple sequences in each one and writes each sequence
    to its own file in the corresponding folder.

    :param filename: name of file
    return: sequence and id
    ============================================================================================="""

    # Get reference folder name and folder for the correpsonding fasta files
    refname = filename.split('/')[-2:]  # First index is ref folder, second is fa file
    refname[1] = refname[1].split('.tfa')[0]  # Remove file extension
    if not os.path.isdir(f'bb_embed/{refname[0]}/{refname[1]}'):
        os.makedirs(f'bb_embed/{refname[0]}/{refname[1]}')

    # Parse fasta file and write each sequence to its own file in the corresponding folder
    with open(filename, 'r', encoding='utf8') as file:
        for seq in SeqIO.parse(file, 'fasta'):
            vec = embed_seq(str(seq.seq), tokenizer, model)
            seqname = seq.id
            with open(f'bb_embed/{refname[0]}/{refname[1]}/{seqname}.txt', 'w', encoding='utf8') as seqfile:
                np.savetxt(seqfile, vec, fmt='%4.6f', delimiter=' ')


def main():
    """=============================================================================================
    ============================================================================================="""

    # Parse reference folder of interest
    path = 'BAliBASE_R1-5/bb3_release/RV11'
    ref_dir = path.rsplit('/', maxsplit=1)[-1]  # Get last directory in path
    fasta_files = parse_ref_folder(path)
    if not os.path.isdir(f'bb_embed/{ref_dir}'):
        os.makedirs(f'bb_embed/{ref_dir}')

    # Load model tokenizer and encoder models
    if os.path.exists('tok.pt'):
        tokenizer = torch.load('tok.pt')
    else:
        tokenizer = T5Tokenizer.from_pretrained("Rostlab/prot_t5_xl_uniref50", do_lower_case=False)
        torch.save(tokenizer, 'tok.pt')
    if os.path.exists('prottrans.pt'):
        model = torch.load('prottrans.pt')
    else:
        model = T5EncoderModel.from_pretrained("Rostlab/prot_t5_xl_uniref50")
        torch.save(model, 'prottrans.pt')

    # Parse each fasta file and write each embedding to its own file
    for file in fasta_files:
        parse_fasta(file, tokenizer, model)


if __name__ == '__main__':
    main()
