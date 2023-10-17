"""This script aligns takes a random sample of residues from alignments and plots their distribution
of cosine similarity and substitution scores to other random residues in the same alignment.

__author__ = "Ben Iovino"
__date__ = 10/17/23
"""

import os
import blosum as bl
import matplotlib.pyplot as plt
import numpy as np
from random import sample
from Bio import SeqIO


def sample_seqs() -> list:
    """Returns a list of (at most) 5 random sequences from each reference MSA

    :return: list of 5 random sequences
    """

    direc = 'data/sequences'
    refs = os.listdir(direc)
    aligns = []

    # For each MSA directory, take (at most) 5 random sequences
    for ref in refs:
        msas = os.listdir(f'{direc}/{ref}')
        for msa in msas:
            msa_dir = f'{direc}/{ref}/{msa}'
            for seq in sample(os.listdir(msa_dir), min(5, len(os.listdir(msa_dir)))):
                aligns.append(f'{direc}/{ref}/{msa}/{seq}')

    return aligns


def sample_pos(seqs: list) -> dict:
    """Returns a dictionary where each key is a sequence and each value is a list of 4 random
    residues

    :param seqs: list of sequences
    :return dict: dictionary of sequences and random residues from each
    """

    # Read fasta file
    fastas = {}
    for seq in seqs:
        sequence = SeqIO.read(seq, 'fasta')
        sample_indices = sample(range(len(sequence)), 4)
        sample_chars = [sequence[i]+str(i) for i in sample_indices]
        fastas[seq] = sample_chars

    return fastas


def get_embeds(seq: str, positions: list) -> list:
    """Returns a list of embeddings for a given sequence and positions

    :param seq: path to fasta sequence
    :param positions: list of residue positions
    :return list: list of embeddings
    """

    # Get positions and embedding file
    positions = [int(pos[1:]) for pos in positions]
    embed = f'data/prot_t5_embed/{"/".join(seq.split("/")[-3:]).split(".")[0]}.txt'

    # Get embedding for each residue
    embeds = []
    with open(embed, 'r', encoding='utf8') as file:
        lines = file.readlines()
    for i in positions:
        embeds.append(lines[i].split()[1:])

    return embeds


def cos_sim(embeds: list) -> list:
    """Returns a list of cosine similarity scores for each pair of embeddings

    :param embeds: list of embeddings
    :return list: list of cosine similarity scores
    """

    sims = []
    for i in range(len(embeds)):  #pylint: disable=C0200
        for j in range(i+1, len(embeds)):
            if i == j:
                continue
            vec1 = [float(x) for x in embeds[i]]
            vec2 = [float(x) for x in embeds[j]]
            sim = np.dot(vec1, vec2) / (np.linalg.norm(vec1) * np.linalg.norm(vec2))
            sims.append(sim*10)

    return sims


def sub_score(residues: list) -> list:
    """Returns a list of substitution scores for each pair of residues

    :param residues: list of residues
    :return list: list of substitution scores
    """

    matrix = bl.BLOSUM(62)
    scores = []
    for i in range(len(residues)):  #pylint: disable=C0200
        for j in range(i+1, len(residues)):
            if i == j:
                continue
            scores.append(matrix[f'{residues[i]}{residues[j]}'])

    return scores


def get_scores(random_pos: dict):
    """Returns two lists of scores for each random residue in each sequence
    """

    # Get cosine similarity and substitution scores for each pair
    sims, subs = [], []
    for seq, positions in random_pos.items():

        # Get residues and embeddings for each sequence
        embeds = get_embeds(seq, positions)
        residues = [pos[0] for pos in positions]

        # Get cosine similarity and substitution scores for each pair
        sims.append(cos_sim(embeds))
        subs.append(sub_score(residues))

    # Join lists of lists
    sims = [item for sublist in sims for item in sublist]
    subs = [item for sublist in subs for item in sublist]

    return sims, subs


def graph_scores(sims: list, subs: list):
    """Plots the distribution of cosine similarity and substitution scores

    :param sims: list of cosine similarity scores
    :param subs: list of substitution scores
    """

    plt.hist(sims, bins=10, alpha=0.5, label='Cosine Similarity')
    plt.hist(subs, bins=10, alpha=0.5, label='BLOSUM62 Score')
    plt.legend(loc='upper right')
    plt.xlabel('Score')
    plt.ylabel('Frequency')
    plt.title('Cosine Similarity vs. BLOSUM62 Scores for Random Residues')
    plt.show()


def main():
    """Main
    """

    seq_sample = sample_seqs()
    random_pos = sample_pos(seq_sample)
    sims, subs = get_scores(random_pos)
    graph_scores(sims, subs)


if __name__ == '__main__':
    main()
