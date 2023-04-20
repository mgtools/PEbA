"""================================================================================================
This script aligns takes residue pairs from alignments, determines their cosine similarity and 
substitution score, and plots the distribution of scores against each other.

Ben Iovino  04/20/23   VecAligns
================================================================================================"""

import os
import sys
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
        file_count = 0
        for pair in pairwise_aligns:
            seq1 = pair[0]
            seq2 = pair[1]
            seq1, seq2 = seq1.split('.')[0], seq2.split('.')[0]  # Remove fa
            align1, align2 = parse_msf(file, seq1, seq2)  # Gather pairwise alignment
            file_path = f'{bb_dir}/{ref_align}/{ref_align}_{file_count}'
            write_align_ca(align1, align2, seq1, seq2, file_path)  # Write pairwise alignment
            file_count += 1


def return_pairs(bb_dir):
    """=============================================================================================
    This function takes a directory of pairwise alignments and returns a list of residue pairs from
    the first pairwise alignment in each MSA.

    :param bb_dir: directory of pairwise alignments
    :return: list of residue pairs
    ============================================================================================="""

    # Get to first pairwise alignment from each MSA
    for direc in os.listdir(bb_dir):
        for file in os.listdir(f'{bb_dir}/{direc}'):
            if file.endswith('msf'):
                seq1, seq2, id1, id2 = parse_align(f'{bb_dir}/{direc}/{file}')

                # Get residue pairs - residues aligned to other residues
                pairs = []
                for i, res in enumerate(seq1):
                    if res != '.' and seq2[i] != '.':  # If aligned
                            pairs.append([res, seq2[i]])
                print(pairs)
                sys.exit()

def main():

    # Read reference alignments from file
    path = 'BAliBASE_R1-5/bb3_release/RV11'

    # Get directory of reference alignments i.e. 'RV11'
    ref_dir = path.rsplit('/', maxsplit=1)[-1]

    # Create unique directory for results, this allows for parallel runs of the script i.e. 'bb_data0'
    bb_ct = 0
    for direc in os.listdir():
        if direc.startswith('bb_data'):
            bb_ct += 1
    bb_dir = f'bb_data{bb_ct}/{ref_dir}'
    os.makedirs(bb_dir)

    # Parse ref folder
    msf_files, fasta_files = parse_ref_folder(path)

    # Sort each list of files to ensure they match up for msf parsing
    msf_files.sort()
    fasta_files.sort()
    parse_aligns(msf_files, fasta_files, bb_dir)

    # Get first 100 residue pairs from first PW align in each MSA
    return_pairs(bb_dir)


if __name__ == "__main__":
    main()
