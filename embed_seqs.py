"""================================================================================================
This script embes protein sequences using a pretrained T5 model and outputs a numpy array of the
embedded sequences to a file.

Ben Iovino  02/16/23  VecAligns
================================================================================================"""

import os
import re
import torch
import numpy as np
from transformers import T5EncoderModel, T5Tokenizer, AutoTokenizer, EsmModel
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
        if file.endswith('.in_tfa'):
            new_file = file.split('.in_tfa')[0] + '.tfa'
            os.rename(f'{path}/{file}',f'{path}/{new_file}')
            fasta_files.append(f'{path}/{new_file}')
    return fasta_files


def prot_t5xl_embed(seq, tokenizer, encoder):
    """=============================================================================================
    This function accepts a protein sequence and returns a list of vectors, each vector representing
    a single amino acid using RostLab's ProtT5_XL_UniRef50 model.

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


def esm2_embed(seq, tokenizer, encoder):
    """=============================================================================================
    This function accepts a protein sequence and returns a list of vectors, each vector representing
    a single amino acid using Facebook's ESM-2 model.

    :param seq: protein sequence
    :param: tokenizer: tokenizer model
    :param encoder: encoder model
    return: list of vectors
    ============================================================================================="""

    inputs = tokenizer(seq, return_tensors="pt")
    outputs = encoder(**inputs)
    last_hidden_states = outputs.last_hidden_state
    return last_hidden_states[0][1:-1]  # First and last tokens are BOS and EOS tokens


def parse_fasta(filename, encoder, tokenizer, model):
    """=============================================================================================
    This function accepts a fasta file with multiple sequences in each one and writes each sequence
    to its own file in the corresponding folder.

    :param filename: name of file
    :param encoder: encoder used
    :param tokenizer: tokenizer model
    :param model: model
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

            # Embed with ProtT5_XL_UniRef50
            if encoder == 'prott5':
                vec = prot_t5xl_embed(str(seq.seq), tokenizer, model)

            # Embed with ESM-2
            if encoder == 'esm2':
                vec = esm2_embed(str(seq.seq), tokenizer, model)
                vec = vec.detach().numpy()

            # Write embeddings to file
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


    # Set an encoder and process each fasta file
    encoder = 'prott5'

    # ProtT5_XL_UniRef50
    if encoder == 'prott5':
        if os.path.exists('t5_tok.pt'):
            tokenizer = torch.load('t5_tok.pt')
        else:
            tokenizer = T5Tokenizer.from_pretrained("Rostlab/prot_t5_xl_uniref50", do_lower_case=False)
            torch.save(tokenizer, 't5_tok.pt')
        if os.path.exists('prot_t5_xl.pt'):
            model = torch.load('prot_t5_xl.pt')
        else:
            model = T5EncoderModel.from_pretrained("Rostlab/prot_t5_xl_uniref50")
            torch.save(model, 'prot_t5_xl.pt')

        for file in fasta_files:
            parse_fasta(file, encoder, tokenizer, model)

    # ESM-2_t36_3B
    if encoder == 'esm2':
        if os.path.exists('auto_tok.pt'):
            tokenizer = torch.load('auto_tok.pt')
        else:
            tokenizer = AutoTokenizer.from_pretrained("facebook/esm2_t36_3B_UR50D")
            torch.save(tokenizer, 'auto_tok.pt')
        if os.path.exists('esm2_t36_3B.pt'):
            model = torch.load('esm2_t36_3B.pt')
        else:
            model = EsmModel.from_pretrained("facebook/esm2_t36_3B_UR50D")
            torch.save(model, 'esm2_t36_3B.pt')

        for file in fasta_files:
            parse_fasta(file, encoder, tokenizer, model)


if __name__ == '__main__':
    main()
