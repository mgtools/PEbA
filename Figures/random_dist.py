"""================================================================================================
This script aligns takes a random sample of residues from alignments and plots their distribution
of cosine similarity and substitution scores to other random residues in the same alignment.

Ben Iovino  04/20/23   VecAligns
================================================================================================"""

import blosum as bl
import numpy as np
import os
import shutil
import sys
from Bio import SeqIO
from random import sample
sys.path.append(os.path.abspath(os.path.join(os.path.pardir, 'VecAligns')))
from compare_aligns import parse_ref_folder
from aligned_dist import parse_aligns, graph_scores


def get_fastas(path):
    """=============================================================================================
    This function takes a directory path with more directories and returns a dict of fasta seqs and
    their ids for each directory.

    :param path: directory path
    :return: dict of fasta files, keys are dirs and values are list of fastas in dir
    ============================================================================================="""

    fastas = {}
    for direc in os.listdir(path):

        # Count number of fasta files in directory
        count = 0
        fastas[direc] = []
        for file in os.listdir(f'{path}/{direc}'):
            if file.endswith('.fa'):
                if count == 5:  # Stop early so we don't have overrepresentation from one MSF
                    break
                count += 1

                # Read in fasta file
                for seq in SeqIO.parse(f'{path}/{direc}/{file}', 'fasta'):
                    sequence = seq.seq
                    seq_id = seq.id

                # Take random sample of residues from sequence
                sample_indices = sample(range(len(sequence)), 10)
                sample_chars = [sequence[i]+str(i) for i in sample_indices]
                fastas[direc].append((sample_chars, seq_id))

    return fastas


def get_embeds(path, fastas):
    """=============================================================================================
    This function takes a directory path and a dictionary of fasta files and returns a dict of
    embeddings for each residue in each fasta file.

    :param path: directory path
    :param fastas: dict of fasta files, keys are dirs and values are list of fastas in dir
    :return: dict of embeddings, keys are dirs and values are list of embeddings
    ============================================================================================="""

    path = f'prot_t5_embed/{path}'

    # Need to open each fasta file in each directory and get embeddings for each residue
    embeds = {}
    for msf, fastas in fastas.items():
        embeds[msf] = []
        for fasta in fastas:

            # Get list of lines that we need to extract from embedding files
            lines = [int(i[1:]) for i in fasta[0]]

            # Open fasta file
            embeddings = []
            with open(f'{path}/{msf}/{fasta[1]}.txt', 'r', encoding='utf8') as file:

                # Check if line number is in list of lines we want to extract
                for i, line in enumerate(file):
                    if i in lines:
                        line = [float(x) for x in line.strip('\n').split()]
                        embeddings.append(line)
            embeds[msf].append(embeddings)

    return embeds


def blosum_scores(fastas, blosum):
    """=============================================================================================
    This function takes a dictionary of fasta sequences and a blosum matrix and returns a list
    of substitution scores.

    :param fastas: dict of fasta files, keys are dirs and values are list of fastas in dir
    :param blosum: blosum matrix
    :return: list of substitution scores
    ============================================================================================="""

    # For each set of fastas from each msf
    scores = []
    for value in fastas.values():

        # Pairwise comparison of each fasta in set
        for fasta in value:
            for fasta2 in value:

                # Ignore if same fasta
                if fasta == fasta2:
                    continue

                # Get substitution scores for each possible pair of residues
                for i in fasta[0]:
                    for j in fasta2[0]:
                        i, j = i[0], j[0]  # Removing index from residue
                        score = blosum[f'{i}{j}']
                        scores.append(score)

    return scores


def cos_scores(embeds):
    """=============================================================================================
    This function takes a dictionary of sequence embeddings and returns a list of cosine similarity
    scores.

    :param fastas: dict of fasta files, keys are dirs and values are list of fastas in dir
    :return: dict of embeddings, keys are dirs and values are list of embeddings
    ============================================================================================="""

    scores = []
    for value in embeds.values():

        # Pairwise comparison of each embedding in set
        for embed in value:
            for embed2 in value:

                # Ignore if same vector
                if embed == embed2:
                    continue

                # Get cosine sim for each possible pair of residues
                for i in embed:
                    for j in embed2:
                        cos_sim = np.dot(i,j)/(np.linalg.norm(i)*np.linalg.norm(j))
                        cos_sim *= 10
                        scores.append(cos_sim)

    return scores


def main():
    """=============================================================================================
    ============================================================================================="""

    # Read reference alignments from file
    path = 'BAliBASE_R1-5/bb3_release/RV913'
    ref_dir = path.rsplit('/', maxsplit=1)[-1]

    # Create directory for storing parsed alignments
    bb_dir = f'bb_data/{ref_dir}'
    os.makedirs(bb_dir)

    # Parse ref folder
    msf_files, fasta_files = parse_ref_folder(path)

    # Sort each list of files to ensure they match up for msf parsing
    msf_files.sort()
    fasta_files.sort()
    parse_aligns(msf_files, fasta_files, bb_dir)

    # Read in sequences from parsed alignments
    fastas = get_fastas(bb_dir)

    # Get corresponding embeddings for each set of residues
    ref_dir = path.rsplit('/', maxsplit=1)[-1]
    embeds = get_embeds(ref_dir, fastas)

    # For each set of sequences, get substitution scores and cosine similarities
    blosum = bl.BLOSUM(62)
    bl_scores = blosum_scores(fastas, blosum)
    sim_scores = cos_scores(embeds)

    # Graph scores
    graph_scores(sim_scores, bl_scores, ref_dir)
    shutil.rmtree(bb_dir)
    os.rmdir('bb_data')


if __name__ == "__main__":
    main()
