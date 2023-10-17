"""This script aligns takes residue pairs from alignments, determines their cosine similarity and 
substitution score, and plots the distribution of scores against each other.

__author__ = "Ben Iovino"
__date__ = 10/17/23
"""

import os
import blosum as bl
import matplotlib.pyplot as plt
import numpy as np
from random import sample


def parse_align(filename: str) -> tuple:
    """Returns a tuple with the two aligned sequences and their names.

    :param filename: name of file (msf format)
    return (str, str, dict): sequences with gaps and dict with beg/end positions
    """

    seq1, seq2 = [], []
    in_align = False
    count = 0
    pos = {}
    with open(filename, 'r', encoding='utf8') as file:
        lines = file.readlines()
    for line in lines:

        # We are "in the alignment" after the '//', read each aligned sequence
        if in_align is True:
            if line.startswith('\n'):
                continue
            if count % 2 == 0:  # Even count indicates first sequence
                seq1.append(''.join(line.split()[1:]))
            else:  # Odd count indicates second sequence
                seq2.append(''.join(line.split()[1:]))
            count += 1
            continue

        # Get beginning/end positions of aligned sequences
        if line.startswith(' Name:'):
            line = line.split()
            pos[line[1]] = line[6]
            continue

        # Alignment starts after '//'
        if line.startswith('//'):
            in_align = True

    # Join the two sequences and return them
    seq1 = ''.join(seq1)
    seq2 = ''.join(seq2)

    return seq1, seq2, pos


def sample_aligns() -> list:
    """Returns a list of 5 random pairwise alignments from each reference MSA

    :return list: list of pairwise alignments
    """

    direc = 'data/alignments/refs'
    refs = os.listdir(direc)
    aligns = []

    # For each MSA directory, take 5 random pairwise alignments
    for ref in refs:
        msas = os.listdir(f'{direc}/{ref}')
        for msa in msas:
            for align in sample(os.listdir(f'{direc}/{ref}/{msa}'), 5):
                aligns.append(f'{direc}/{ref}/{msa}/{align}')

    return aligns


def get_pairs(align1: str, align2: str) -> list:
    """Returns a list of residue pairs from two alignments

    :param align1: first aligned sequence
    :param align2: second aligned sequence
    :return list: list of residue pairs
    """

    pairs = []
    seq1_count, seq2_count = 0, 0
    for i, pos in enumerate(align1):
        if pos != '.' and align2[i] != '.':
            pairs.append((pos+str(seq1_count), align2[i]+str(seq2_count)))
        if pos != '.':
            seq1_count += 1
        if align2[i] != '.':
            seq2_count += 1

    return pairs


def sample_pairs(align_sample: list) -> dict:
    """Returns a dictionary of (at most) 100 random residue pairs from each alignment

    :param align_sample: list of pairwise alignments
    :return dict: dict where key is alignment and value is a list of residue pairs
    """

    # Randomly sample 100 residue pairs from each alignment
    random_pairs = {}
    for align in align_sample:
        seq1, seq2, _ = parse_align(align)
        pairs = get_pairs(seq1, seq2)

        # Remove first five pairs from list due to abnormally high cosine sim scores
        pairs = pairs[4:]
        if len(pairs) > 100:
            pairs = sample(pairs, 100)
        random_pairs[align] = pairs

    return random_pairs


def get_embeds(embed_path: str, seqid: str, pos: list) -> list:
    """Returns a list of embeddings for a given sequence and positions

    :param embed_path: path to embedding file
    :param seqid: sequence id
    :param pos: list of positions
    :return list: list of embeddings
    """

    embeds = []
    with open(f'{embed_path}/{seqid}.txt', 'r', encoding='utf8') as file:
        lines = file.readlines()
    for i in pos:
        embeds.append(lines[i].split()[1:])

    return embeds


def cos_sim(embed1: list, embed2: list) -> list:
    """Returns the cosine similarity for each position between two embeddings

    :param embed1: first embedding
    :param embed2: second embedding
    :return float: list of cosine similarity scores
    """

    cos_sims = []
    for i, vec in enumerate(embed1):

        # Convert strings to lists of floats
        vec1 = [float(x) for x in vec]
        vec2 = [float(x) for x in embed2[i]]
        sim = np.dot(vec1, vec2) / (np.linalg.norm(vec1) * np.linalg.norm(vec2))
        cos_sims.append(sim*10)

    return cos_sims


def sub_score(seq1: str, seq2: str) -> list:
    """Returns the substitution score for each residue pair between two sequences

    :param seq1: first sequence
    :param seq2: second sequence
    :return list: list of substitution scores
    """

    blosum = bl.BLOSUM(62)
    scores = []
    for i, char1 in enumerate(seq1):
        char2 = seq2[i]
        scores.append(blosum[f'{char1}{char2}'])

    return scores


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
    plt.title('Cosine Similarity vs. BLOSUM62 Scores for Aligned Residues')
    plt.show()


def get_scores(random_pairs: dict) -> tuple:
    """Returns two lists of scores for each residue pair

    :param random_pairs: dict where key is alignment and value is a list of residue pairs
    :return (list, list): lists of cosine similarity and substitution scores
    """

    # Get cosine similarity and substitution scores for each pair
    sims, subs = [], []
    for align, pairs in random_pairs.items():

        # Get positions of aligned residue pairs
        seq1, seq2, pos1, pos2 = [], [], [], []
        for pair in pairs:
            seq1.append(pair[0][0])
            seq2.append(pair[1][0])
            pos1.append(int(pair[0][1:]))
            pos2.append(int(pair[1][1:]))

        # Get embeddings for each residue
        id1 = align.split('/')[-1].split('-')[0]
        id2 = align.split('/')[-1].split('-')[1].strip('.msf')
        embed_path = f'data/prot_t5_embed/{"/".join(align.split("/")[-3:-1])}'
        embed1 = get_embeds(embed_path, id1, pos1)
        embed2 = get_embeds(embed_path, id2, pos2)

        # Add scores to lists
        sims.append(cos_sim(embed1, embed2))
        subs.append(sub_score(seq1, seq2))

    # Join lists of lists
    sims = [item for sublist in sims for item in sublist]
    subs = [item for sublist in subs for item in sublist]

    return sims, subs


def main():
    """Main
    """

    align_sample = sample_aligns()
    random_pairs = sample_pairs(align_sample)
    sims, subs = get_scores(random_pairs)
    graph_scores(sims, subs)


if __name__ == '__main__':
    main()
