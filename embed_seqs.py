"""This script embes protein sequences using a pretrained language model and writes them
each to a file.

__author__ = "Ben Iovino"
__date__ = 09/22/23
"""

import esm
import os
import re
import torch
import numpy as np
from transformers import T5EncoderModel, T5Tokenizer
from Bio import SeqIO


def parse_ref_folder(path: str):
    """Returns list of fasta files in a directory

    :param path: directory path to folder
    :return: list of fasta files
    """

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


def prot_t5xl_embed(seq: str, tokenizer, encoder, device: str) -> np.ndarray:
    """Returns an embedding of a protein sequence using ProtT5_XL_UniRef50

    :param seq: protein sequence
    :param: tokenizer: tokenizer model
    :param encoder: encoder model
    :param device: cpu/gpu
    return np.ndarray: list of vectors
    """

    # Remove special chars, add space after each amino acid so each residue is vectorized
    seq = re.sub(r"[UZOB]", "X", seq)
    seq = [' '.join([*seq])]

    # Tokenize, encode, and load sequence
    ids = tokenizer.batch_encode_plus(seq, add_special_tokens=True, padding=True)
    input_ids = torch.tensor(ids['input_ids']).to(device)  # pylint: disable=E1101
    attention_mask = torch.tensor(ids['attention_mask']).to(device)  # pylint: disable=E1101

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


def esm2_embed(seq: str, tokenizer, encoder, device: str) -> np.ndarray:
    """Returns an embedding of a protein sequence using ESM-2

    :param seq: protein sequence
    :param: tokenizer: tokenizer model
    :param encoder: encoder model
    :param device: cpu/gpu
    return np.ndarray: list of vectors
    """

    # Embed sequence
    seq_str = str(seq.seq).upper()  # tok does not convert to uppercase
    seq = np.array([seq.id, seq_str], dtype=object)
    _, _, batch_tokens = tokenizer([seq])  # id and seq are batched together
    batch_tokens = batch_tokens.to(device)  # send tokens to gpu

    # Encode
    with torch.no_grad():
        results = encoder(batch_tokens)
    embed = results['logits'].cpu().numpy()

    return embed[0]


def parse_fasta(filename: str, encoder, tokenizer, model, device: str):
    """Writes each sequence in a fasta file as an embedding to a file

    :param filename: name of file
    :param encoder: encoder used
    :param tokenizer: tokenizer model
    :param model: model
    :param device: cpu/gpu
    """

    # Get reference folder name and folder for the correpsonding fasta files
    refname = filename.split('/')[-2:]  # First index is ref folder, second is fa file
    refname[1] = refname[1].split('.tfa')[0]  # Remove file extension
    if not os.path.isdir(f'data/embeddings/{refname[0]}/{refname[1]}'):
        os.makedirs(f'data/embeddings/{refname[0]}/{refname[1]}')

    # Parse fasta file and write each sequence to its own file in the corresponding folder
    with open(filename, 'r', encoding='utf8') as file:
        for seq in SeqIO.parse(file, 'fasta'):

            # Embed with ProtT5_XL_UniRef50
            if encoder == 'prott5':
                vec = prot_t5xl_embed(str(seq.seq), tokenizer, model, device)

            # Embed with ESM-2
            if encoder == 'esm2':
                vec = esm2_embed(seq, tokenizer, model, device)

            # Write embeddings to file
            seqname = seq.id
            with open(f'data/embeddings/{refname[0]}/{refname[1]}/{seqname}.txt',
                       'w', encoding='utf8') as seqfile:
                np.savetxt(seqfile, vec, fmt='%4.6f', delimiter=' ')


def main():
    """Main
    """

    # Parse reference folder of interest
    path = 'data/BAliBASE_R1-5/RV11'
    ref_dir = path.rsplit('/', maxsplit=1)[-1]  # Get last directory in path
    fasta_files = parse_ref_folder(path)
    if not os.path.isdir(f'data/embeddings/{ref_dir}'):
        os.makedirs(f'data/embeddings/{ref_dir}')

    # Set an encoder and process each fasta file
    encoder = 'prott5'
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')  #pylint: disable=E1101

    # ProtT5_XL_UniRef50
    if encoder == 'prott5':
        tokenizer = T5Tokenizer.from_pretrained('Rostlab/prot_t5_xl_uniref50', do_lower_case=False)
        model = T5EncoderModel.from_pretrained("Rostlab/prot_t5_xl_uniref50")
        model.to(device)  # Loads to GPU if available
        for file in fasta_files:
            parse_fasta(file, encoder, tokenizer, model, device)

    # ESM-2_t36_3B
    if encoder == 'esm2':
        model, alphabet = esm.pretrained.esm2_t36_3B_UR50D()
        tokenizer = alphabet.get_batch_converter()
        model.eval()  # disables dropout for deterministic results
        model.to(device)
        for file in fasta_files:
            parse_fasta(file, encoder, tokenizer, model, device)


if __name__ == '__main__':
    main()
