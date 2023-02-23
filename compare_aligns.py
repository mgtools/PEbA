"""================================================================================================
This script parses BAliBASE reference folders. MSF files are written pairwise alignments between
each pair of two sequences. It also parses FASTA files to then write each sequence to an individual
FASTA file. It then creates global BLOSUM and PEbA alignments between each pair of sequences for the
purpose of comparison to the reference BAliBASE pairwise alignments. Each PW align is compared to
the ref, the results are stored in a csv, and then a scatterplot is created to show the difference
between the similarity score of PEbA vs. BLOSUM.

Ben Iovino  02/10/23   VecAligns
================================================================================================"""

import os
import sys
from time import strftime
from random import sample
from Bio import SeqIO
import matplotlib.pyplot as plt


def parse_ref_folder(path):
    """=============================================================================================
    This function accepts a folder name parses every MSF and FASTA file to write each pairwise
    alignment and fasta sequence to their own files.

    :param path: directory path to folder
    :return: list of msf files, list of fasta files
    ============================================================================================="""

    # Get MSF and FASTA files
    msf_files = []
    fasta_files = []
    for file in os.listdir(path):
        if file.endswith('.msf'):  # Append to msf list
            msf_files.append(f'{path}/{file}')
        if file.endswith('.tfa'):  # Append to fasta list
            fasta_files.append(f'{path}/{file}')
        if file.endswith('.in_tfa'):
            new_file = file.split('.in_tfa')[0] + '.tfa'
            os.rename(f'{path}/{file}',f'{path}/{new_file}')
            fasta_files.append(f'{path}/{new_file}')
    return msf_files, fasta_files


def parse_fasta(filename, bb_dir):
    """=============================================================================================
    This function accepts a fasta file with multiple sequences in each one and writes each sequence
    to its own file in the corresponding folder.

    :param filename: name of file
    :param bb_dir: directory to place parsed files
    return: sequence and id
    ============================================================================================="""

    # Get reference folder name and folder for the correpsonding fasta files
    famname = filename.split('/')[-1]  # Get fasta file name
    famname = famname.split('.tfa')[0]  # Remove file extension
    if not os.path.isdir(f'{bb_dir}/{famname}'):
        os.makedirs(f'{bb_dir}/{famname}')

    # Parse fasta file and write each sequence to its own file in the corresponding folder
    seq = ''
    seqs = []
    with open(filename, 'r', encoding='utf8') as file:
        for seq in SeqIO.parse(file, 'fasta'):
            SeqIO.write(seq, f'{bb_dir}/{famname}/{seq.id}.fa', 'fasta')
            seqs.append(f'{seq.id}.fa')
    return seqs


def parse_msf(filename, id1, id2):
    """=============================================================================================
    This function accepts an MSF file and returns the pairwise alignment between two sequences.

    :param filename: name of file
    :param id1: first sequence id
    :param id2: second sequence id
    return: align1, align2 - corresponding pairwise alignments
    ============================================================================================="""

    seq1 = []
    seq2 = []
    with open(filename, 'r', encoding='utf8') as file:
        for line in file:
            if line.startswith(id1):
                seq1.append(''.join(line.split()[1:]))
            elif line.startswith(id2):
                seq2.append(''.join(line.split()[1:]))

    # Go through both sequences and remove positions with gaps in both
    seq1 = list(''.join(seq1))
    seq2 = list(''.join(seq2))
    for i in range(len(seq1)):  # pylint: disable=C0200
        if seq1[i] == '.' and seq2[i] == '.':
            seq1[i] = ''
            seq2[i] = ''
    align1 = ''.join(seq1)
    align2 = ''.join(seq2)
    return align1, align2


def write_align(seq1, seq2, id1, id2, path):
    """=============================================================================================
    This function accepts two sequences after gaps have been introduced and writes them to a file
    in MSF format, with some extra information about the alignment parameters.

    :param seq1: first aligned sequence
    :param seq2: second aligned sequence
    :param id1: first sequence id
    :param id2: second sequence id
    :param path: directory
    ============================================================================================="""

    # Length of alignment
    length = len(seq1)

    # Add space every 10 characters
    seq1 = [seq1[i:i+10] for i in range(0, len(seq1), 10)]
    seq1 = ' '.join(seq1)
    seq2 = [seq2[i:i+10] for i in range(0, len(seq2), 10)]
    seq2 = ' '.join(seq2)

    # Split sequences every 50 characters
    seq1_split = [seq1[i:i+55] for i in range(0, len(seq1), 55)]
    seq2_split = [seq2[i:i+55] for i in range(0, len(seq2), 55)]

    # Add extra spaces to either id if they are not the same length
    if len(id1) != len(id2):
        if len(id1) > len(id2):
            id2 = id2 + ' ' * (len(id1) - len(id2))
        else:
            id1 = id1 + ' ' * (len(id2) - len(id1))

    # Write to a new line for every index in the split list i.e. every 55 characters
    with open(f'{path}.msf', 'w', encoding='utf8') as file:
        file.write('PileUp\n\n\n\n')
        file.write(f'   MSF:  {length}  Type:  P\n\n')
        file.write(f' Name: {id1} oo  Len:  {length}\n')
        file.write(f' Name: {id2} oo  Len:  {length}\n\n//\n\n\n\n')
        for i in range(len(seq1_split)):  # pylint: disable=C0200
            file.write(f'{id1}      {seq1_split[i]}\n')
            file.write(f'{id2}      {seq2_split[i]}\n\n')


def parse_align_files(msf_files, fasta_files, bb_dir):
    """=============================================================================================
    This function accepts lists of two sets of files and a directory to place them in where they
    are parsed correspondingly. As they are parsed, they are also aligned using global_align.py
    and PEbA_align.py.

    :param msf_files: list of msf files
    :param fasta_files: list of fasta files
    :param bb_dir: directory to place files in
    :param tokenizer: loaded tokenizer
    :param model: loaded encoder
    :return: arguments used to call matrix and PEbA alignments
    ============================================================================================="""

    # Parse each fasta file, store names of each for subsequent msf parsing
    seqs = []
    for file in fasta_files:
        new_seqs = parse_fasta(file, bb_dir)
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

        # Set a sample size for the PW aligns - sometimes there are 1000+ pairs
        sample_size = 1
        if len(pairwise_aligns) > sample_size:
            pairwise_aligns = sample(pairwise_aligns, sample_size)

        # For the selected pairs, align them using local programs
        file_count = 0
        for pair in pairwise_aligns:
            seq1 = pair[0]  #pylint: disable=E1136
            seq2 = pair[1]  #pylint: disable=E1136
            ref_dir = bb_dir.split('/')[1]  # Embedding dirs don't change like results dirs
            args = (f'-file1 {bb_dir}/{ref_align}/{seq1} '
                    f'-file2 {bb_dir}/{ref_align}/{seq2} '
                    f'-embed1 bb_embed/{ref_dir}/{ref_align}/{seq1.split(".")[0]}.txt '
                    f'-embed2 bb_embed/{ref_dir}/{ref_align}/{seq2.split(".")[0]}.txt '
                    f'-gopen {-11} '
                    f'-gext {-1} ')
            print(f'{strftime("%H:%M:%S")} PEbA: {ref_align}/{seq1} and {ref_align}/{seq2}\n',
                           file=sys.stdout)
            os.system(f"python local_PEbA.py {args}")

            args = (f'-file1 {bb_dir}/{ref_align}/{seq1} '
                    f'-file2 {bb_dir}/{ref_align}/{seq2}')
                    # f'-gopen {-11} '
                    # f'-gext {-1} '
                    # f'-matrix blosum '
                    # f'-score {45}')
            print(f'{strftime("%H:%M:%S")} MATRIX: {ref_align}/{seq1} and {ref_align}/{seq2}\n',
                           file=sys.stdout)
            os.system(f"python run_DEDAL.py {args}")

            # Grab alignment from reference MSA
            seq1, seq2 = seq1.split('.')[0], seq2.split('.')[0]  # Remove fa
            align1, align2 = parse_msf(file, seq1, seq2)  # Gather pairwise alignment
            file_path = f'{bb_dir}/{ref_align}/{ref_align}_{file_count}'
            write_align(align1, align2, seq1, seq2, file_path)  # Write pairwise alignment
            file_count += 1
    return args


def compare_aligns(path):
    """=============================================================================================
    This function takes a directory and compares matrix and PEbA alignments to the reference
    alignment using t_coffee's aln_compare function. The output is stored in a text file.

    :param path: directory where alignments exist
    ============================================================================================="""

    folders = os.listdir(path)
    for folder in folders:
        files = os.listdir(f'{path}/{folder}')
        matrix_aligns = []
        peba_aligns = []
        ref_aligns = []
        for file in files: # Add each alignment to list of alignments
            if file.startswith('MATRIX'):
                matrix_aligns.append(f'{path}/{folder}/{file}')
            if file.startswith('PEbA'):
                peba_aligns.append(f'{path}/{folder}/{file}')
            if file.startswith('BB'):
                ref_aligns.append(f'{path}/{folder}/{file}')

        # Sort so that correct alignments are compared
        matrix_aligns.sort()
        peba_aligns.sort()
        ref_aligns.sort()

        # Call t_coffee to compare global and peba aligns to refs
        for i, ref_align in enumerate(ref_aligns):
            matrix_align = matrix_aligns[i]
            peba_align = peba_aligns[i]

            # Get names of alignments for output file
            matrix_name = matrix_align.split('/')[-1].strip('.msf')
            peba_name = peba_align.split('/')[-1].strip('.msf')
            print(f'Comparing {matrix_name} and {peba_name} to {ref_align}')
            os.system(f't_coffee -other_pg aln_compare -al1 {matrix_align} -al2 {ref_align} > {path}/{folder}/matrix_{matrix_name}_compare.txt')
            os.system(f't_coffee -other_pg aln_compare -al1 {peba_align} -al2 {ref_align} > {path}/{folder}/peba_{peba_name}_compare.txt')


def parse_compare(path):
    """=============================================================================================
    This function takes a directory and compares the global and PEbA alignment comparisons to the
    reference alignment.

    :param path: directory where alignments exist
    ============================================================================================="""

    folders = os.listdir(path)
    for folder in folders:
        files = os.listdir(f'{path}/{folder}')
        matrix_compares = []
        peba_compares = []
        for file in files:
            if file.startswith('matrix'):
                matrix_compares.append(f'{path}/{folder}/{file}')
            if file.startswith('peba'):
                peba_compares.append(f'{path}/{folder}/{file}')

        # Sort so that correct comparisons are compared
        matrix_compares.sort()
        peba_compares.sort()

        # For both compare files, store their third line in a csv
        matrix_lines = []
        peba_lines = []
        for i, compare in enumerate(matrix_compares):
            with open(compare, 'r', encoding='utf8') as file1:
                lines = file1.readlines()
                third_line = lines[2].split()
                third_line = f'{", ".join(third_line[0:4])}, {third_line[-1].strip("]")}\n'  # Columns used not important
                matrix_lines.append(third_line)
                os.remove(compare)
            with open(peba_compares[i], 'r', encoding='utf8') as file2:
                lines = file2.readlines()
                third_line = lines[2].split()
                third_line = f'{", ".join(third_line[0:4])}, {third_line[-1].strip("]")}\n'  # Columns used not important
                peba_lines.append(third_line)
                os.remove(peba_compares[i])
        with open(f'{path}/{folder}/compare.csv', 'w', encoding='utf8') as file3:
            for i, line in enumerate(matrix_lines):
                file3.write(line)
                file3.write(peba_lines[i])


def graph_compare(path, matrix):
    """=============================================================================================
    This function takes a directory and makes two graphs. The first graphs the differences between 
    the global and PEbA alignments compared the reference alignment and the second graphs the sim
    scores against each other with the better BLOSUM scores on the left side of the diagonal line
    and the better PEbA scores on the right side of the diagonal line.

    :param path: directory where compare.csv files exist
    :param matrix: matrix used for substitution matrix alignments
    ============================================================================================="""

    # Get the similarity scores from the compare files
    matrix_sim = []
    peba_sim = []
    folders = os.listdir(path)
    for folder in folders:
        with open(f'{path}/{folder}/compare.csv', 'r', encoding='utf8') as file:
            for line in file:
                line = line.split(',')
                if 'MATRIX' in line[0]:
                    matrix_sim.append(float(line[3]))
                if 'PEbA' in line[0]:
                    peba_sim.append(float(line[3]))

    # Average the similarity scores for graphing
    blosum_avg = round(sum(matrix_sim)/len(matrix_sim), 1)
    peba_avg = round(sum(peba_sim)/len(peba_sim), 1)

    # Graph the difference between similarity scores for each alignment
    fig = plt.figure()
    ax = fig.add_subplot()
    sim_diff = []
    for i, mat_sim in enumerate(matrix_sim):
        sim_diff.append(peba_sim[i]-mat_sim)
    ax.scatter(list(range(1, len(sim_diff) + 1)), sim_diff)
    ax.set_title(f'Difference in TCS PEbA (Avg={peba_avg}) vs. {matrix} (Avg={blosum_avg})')
    ax.set_xlabel('Alignment Number')
    ax.set_ylabel('Similarity Difference')
    ax.set_ylim(-20, 80)
    ax.axhline(0, color='black')
    plt.savefig(f'{path}/differences.png')

    # Graph the similarity scores against each other, better score on either side of diag line
    fig = plt.figure()
    ax = fig.add_subplot()
    peba_scores = []
    matrix_scores = []
    for i, mat_sim in enumerate(matrix_sim):
        if mat_sim < peba_sim[i]:
            peba_scores.append([peba_sim[i], mat_sim])
        else:
            matrix_scores.append([mat_sim, peba_sim[i]])
    ax.scatter([i[0] for i in matrix_scores], [i[1] for i in matrix_scores], color='blue')
    ax.scatter([i[1] for i in peba_scores], [i[0] for i in peba_scores], color='red')
    ax.set_title(f'PEbA Alignment (Avg={peba_avg}) vs. {matrix} Alignment (Avg={blosum_avg})')
    ax.set_xlabel('TCS MATRIX')
    ax.set_ylabel('TCS PEbA')
    plt.plot([0, 100], [0, 100], color='black')
    plt.savefig(f'{path}/comparison.png')


def main():
    """=============================================================================================
    This function calls parse_ref_folder to get lists of all the msf and tfa fasta files in the
    reference directory of interest. It then calls parse_align_files to parse each tfa file and msf
    file, while also aligning each pariwise comparison of fasta sequences. These alignments are
    compared using t_coffee's 'aln_compare' function and the results are parsed and graphed.
    ============================================================================================="""

    # Get directory of reference alignments
    path = 'BAliBASE_R1-5/bb3_release/RV11'
    ref_dir = path.rsplit('/', maxsplit=1)[-1]

    # Create unique directory for results, this allows for parallel runs of the script
    bb_ct = 0
    for direc in os.listdir():
        if direc.startswith('bb_data'):
            bb_ct += 1
    bb_dir = f'bb_data{bb_ct}/{ref_dir}'
    os.makedirs(bb_dir)

    # Parse reference folder of interest
    print(f'{strftime("%H:%M:%S")} Parsing and computing alignments...\n', file=sys.stdout)
    msf_files, fasta_files = parse_ref_folder(path)

    # Sort each list of files to ensure they match up for msf parsing
    msf_files.sort()
    fasta_files.sort()
    args = parse_align_files(msf_files, fasta_files, bb_dir)

    # Get type of matrix used from args
    split_args = args.split('-')
    matrix = split_args[7].split(' ')[1]
    score = split_args[8].split(' ')[1]
    matrix = matrix+score

    # Compare alignments using t_coffee
    print(f'{strftime("%H:%M:%S")} Comparing alignments...\n', file=sys.stdout)
    compare_aligns(bb_dir)
    parse_compare(bb_dir)
    graph_compare(bb_dir, matrix)
    print(f'{strftime("%H:%M:%S")} Program Complete!\n', file=sys.stdout)


if __name__ == '__main__':
    main()
