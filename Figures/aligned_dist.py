"""================================================================================================
This script aligns takes residue pairs from alignments, determines their cosine similarity and 
substitution score, and plots the distribution of scores against each other.

Ben Iovino  04/20/23   VecAligns
================================================================================================"""

import os
import shutil
import sys
import regex as re
import blosum as bl
import matplotlib.pyplot as plt
import numpy as np
from random import sample
sys.path.append(os.path.abspath(os.path.join(os.path.pardir, 'VecAligns')))
from compare_aligns import parse_ref_folder, parse_fasta_ca, parse_msf, write_align_ca
from compute_score import parse_align


def parse_aligns(msf_files, fasta_files, bb_dir):
    """=============================================================================================
    This function accepts lists of two sets of files and a directory to place them in where they
    are parsed correspondingly.

    :param msf_files: list of msf files
    :param fasta_files: list of fasta files
    :param bb_dir: directory to place files in
    ============================================================================================="""

    # Parse each fasta file, store names of each for subsequent msf parsing
    seqs = []
    for file in fasta_files:
        new_seqs = parse_fasta_ca(file, bb_dir)
        seqs.append(new_seqs)  # Store in nested list to access only relevant fa files for each msf

    # Each MSF files corresponds to a set of fasta files
    for i, file in enumerate(msf_files):
        ref_align = file.rsplit('/', maxsplit=1)[-1].strip('.msf')  # Get name of ref alignment

        # Get all pairwise alignments from the fasta files correpsonding to the MSF file
        sequences = seqs[i]
        pairwise_aligns = []
        for i, seq in enumerate(sequences):
            loop_count = i  # Keep track of number of loops so no repeats occur
            while loop_count != len(sequences):
                if seq != sequences[loop_count]:  # Don't want to align a sequence to itself
                    pairwise_aligns.append([seq, sequences[loop_count]])
                loop_count+=1

        # For the selected pairs, take the pairwise alignment from the reference MSA
        for j, pair in enumerate(pairwise_aligns):
            seq1, seq2 = pair[0], pair[1]
            seq1, seq2 = seq1.split('.')[0], seq2.split('.')[0]  # Remove fa
            align1, align2 = parse_msf(file, seq1, seq2)  # Gather pairwise alignment
            file_path = f'{bb_dir}/{ref_align}/{ref_align}_{j}'
            write_align_ca(align1, align2, seq1, seq2, file_path)  # Write pairwise alignment


def return_pairs(bb_dir):
    """=============================================================================================
    This function takes a directory of pairwise alignments and returns a list of residue pairs from
    the first five pairwise alignments from each MSA.

    :param bb_dir: directory of pairwise alignments
    :return: list of residue pairs
    ============================================================================================="""

    # Get to first pairwise alignment from each MSA
    pairs_dict = {}
    for direc in os.listdir(bb_dir):

        # Not all files are msf's, need a counter instead of loop iterator
        count = 0
        for file in os.listdir(f'{bb_dir}/{direc}'):

            # Break at 5 so some MSA's are not over represented
            if count == 5:
                break

            if file.endswith('msf'):
                count += 1
                seq1, seq2, id1, id2 = parse_align(f'{bb_dir}/{direc}/{file}')

                # Get residue pairs - residues aligned to other residues
                pairs = []
                seq1_count, seq2_count = 0, 0
                for i, res in enumerate(seq1):

                    # Check for aligned pairs
                    if res != '.' and seq2[i] != '.':
                        pairs.append([res+str(seq1_count), seq2[i]+str(seq2_count)])

                    # Add to counts after checking for alignment so we start at position 0
                    if res != '.':
                        seq1_count += 1
                    if seq2[i] != '.':
                        seq2_count += 1

                # Remove first five pairs from list due to abnormally high similarity scores
                pairs = pairs[4:]

                # Take a random sample of 100 pairs from list (if eligible)
                if len(pairs) > 100:
                    pairs = sample(pairs, 100)

                # Add pairs and ids and dict
                pairs_dict[file] = [pairs, id1, id2]

    return pairs_dict


def get_embeddings(path, pos):
    """=============================================================================================
    This function takes directory to a set of protein embeddings and a list of positions that are to
    be extracted. It returns a list of embeddings corresponding to the positions.

    :param path: directory to embeddings
    :param pos: list of positions
    :return: list of embeddings
    ============================================================================================="""

    # Remove msf file number and extension from path
    path = re.sub(r'_+\d+\.msf', '', path)

    embed = []
    with open(path, 'r', encoding='utf8') as file:
        for i, line in enumerate(file):
            if i in pos:
                embed.append(line)
                continue
            if i > max(pos):  # Break after reading last position
                break

    return embed


def cos_sim(vecs1, vecs2):
    """=============================================================================================
    This function takes two equal length lists of vectors and returns a list of their cosine
    similarities.

    :param vecs1: first list of vectors
    :param vecs2: second list of vectors
    :return: list of cosine similarities
    ============================================================================================="""

    cos_sims = []
    for i, vec in enumerate(vecs1):

        # Convert strings to lists of floats
        vec1 = [float(x) for x in vec.strip('\n').split()]
        vec2 = [float(x) for x in vecs2[i].strip('\n').split()]
        sim = np.dot(vec1, vec2) / (np.linalg.norm(vec1) * np.linalg.norm(vec2))
        cos_sims.append(sim*10)

    return cos_sims


def sub_score(matrix, seq1, seq2):
    """=============================================================================================
    This function takes two equal length lists of residues and returns a list of their BLOSUM
    substitution scores.

    :param matrix: substitution matrix
    :param seq1: first list of residues
    :param seq2: second list of residues
    :return: list of substitution scores
    ============================================================================================="""

    scores = []
    for i, res in enumerate(seq1):
        scores.append(matrix[f'{res}{seq2[i]}'])

    return scores


def get_scores(pairs, matrix, path):
    """=============================================================================================
    This function takes a set of residue pairs, a score for the BLOSUM substitution matrix and a
    directory to a set of protein embeddings corresponding to the residue pairs. It returns a list
    of substitution scores and cosine similarities from all of the pairs.

    :param pairs: dict with directory names as keys and residue pairs/sequence ids as values
    :param matrix: score for BLOSUM substitution matrix
    :param path: directory to embeddings
    ============================================================================================="""

    blosum = bl.BLOSUM(matrix)
    sims, subs = [], []

    # For each pair, get the cosine similarity and substitution score
    for key, item in pairs.items():

        # Get the positions of the aligned residue pairs from the sequences
        seq1, seq2, pos1, pos2 = [], [], [], []
        for pair in item[0]:
            seq1.append(pair[0][0])
            seq2.append(pair[1][0])
            pos1.append(int(pair[0][1:]))
            pos2.append(int(pair[1][1:]))

        # Get embeddings from files
        embed_path = f'{path}/{key}/{item[1]}.txt'
        embed1 = get_embeddings(embed_path, pos1)
        embed_path = f'{path}/{key}/{item[2]}.txt'
        embed2 = get_embeddings(embed_path, pos2)

        # Calculate cosine similarity for each pair of vectors
        sim = cos_sim(embed1, embed2)
        sims.append(sim)

        # Calculate substitution score for each pair of residues
        sub = sub_score(blosum, seq1, seq2)
        subs.append(sub)

    # Join lists of lists
    sims = [item for sublist in sims for item in sublist]
    subs = [item for sublist in subs for item in sublist]

    return sims, subs


def graph_scores(cos, sub):
    """=============================================================================================
    This function takes two lists of values and plots their distributions on a histogram.

    :param cos: list of cosine similarity values
    :param scores2: list of substitution scores
    ============================================================================================="""

    plt.hist(cos, bins=10, alpha=0.5, label='Cosine Similarity')
    plt.hist(sub, bins=10, alpha=0.5, label='BLOSUM62')
    plt.legend(loc='upper right')
    plt.title('Cosine Similarity vs. BLOSUM62 Scores for Aligned Residues')
    plt.xlabel('Score')
    plt.ylabel('Frequency')
    plt.show()


def main():
    """=============================================================================================
    This function takes a directory of multiple sequence alignments, parses them into pairwise
    alignments, and calculates cosine similarities and substitution scores for the first 100 residue
    pairs for the first five pairwise alignments from each MSA. It then plots the distributions of
    the scores.
    ============================================================================================="""

    paths = ['Data/BAliBASE_R1-5/bb3_release/RV11', 'Data/BAliBASE_R1-5/bb3_release/RV12',
             'Data/BAliBASE_R1-5/bb3_release/RV911', 'Data/BAliBASE_R1-5/bb3_release/RV912',
             'Data/BAliBASE_R1-5/bb3_release/RV913']

    bl_scores = []
    sim_scores = []

    for path in paths:
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

        # Get first 100 residue pairs from first PW align in each MSA
        pairs = return_pairs(bb_dir)

        # For each set of pairs, get cosine similarities and substitution scores
        embed_path = f'prot_t5_embed/{ref_dir}'
        sims, subs = get_scores(pairs, 62, embed_path)
        bl_scores.append(subs)
        sim_scores.append(sims)
        shutil.rmtree(bb_dir)

    # Combine lists of scores
    bl_scores = [item for sublist in bl_scores for item in sublist]
    sim_scores = [item for sublist in sim_scores for item in sublist]

    # Graph cosine similarities and substitution scores
    graph_scores(sim_scores, bl_scores)
    os.rmdir('bb_data')


if __name__ == "__main__":
    main()
