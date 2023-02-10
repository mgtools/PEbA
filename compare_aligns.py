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
from random import sample
from Bio import SeqIO
import matplotlib.pyplot as plt


def parse_ref_folder(path):
    """=============================================================================================
    This function accepts a folder name parses every MSF and FASTA file to write each pairwise
    alignment and fasta sequence to their own files.

    :param path: directory path to folder
    ============================================================================================="""

    # Get MSF and FASTA files
    msf_files = []
    fasta_files = []
    for file in os.listdir(path):
        if file.endswith('.msf'):  # Append to msf list
            msf_files.append(f'{path}/{file}')
        if file.endswith('.tfa'):  # Append to fasta list
            fasta_files.append(f'{path}/{file}')
    return msf_files, fasta_files


def parse_fasta(filename):
    """=============================================================================================
    This function accepts a fasta file with multiple sequences in each one and writes each sequence
    to its own file in the corresponding folder.

    :param filename: name of file
    return: sequence and id
    ============================================================================================="""

    # Get reference folder name and folder for the correpsonding fasta files
    refname = filename.split('/')[-2:]  # First index is ref folder, second is fa file
    refname[1] = refname[1].split('.tfa')[0]  # Remove file extension
    if not os.path.isdir(f'bb_data/{refname[0]}/{refname[1]}'):
        os.makedirs(f'bb_data/{refname[0]}/{refname[1]}')

    # Parse fasta file and write each sequence to its own file in the corresponding folder
    seq = ''
    seqs = []
    with open(filename, 'r', encoding='utf8') as file:
        for seq in SeqIO.parse(file, 'fasta'):
            SeqIO.write(seq, f'bb_data/{refname[0]}/{refname[1]}/{seq.id}.fa', 'fasta')
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

            
def parse_align_files(msf_files, fasta_files, ref_dir):
    """=============================================================================================
    This function accepts lists of two sets of files and a directory to place them in where they
    are parsed correspondingly. As they are parsed, they are also aligned using global_align.py
    and PEbA_align.py.

    :param msf_files: list of msf files
    :param fasta_files: list of fasta files
    :param ref_dir: directory to place files in
    :param tokenizer: loaded tokenizer
    :param model: loaded encoder
    ============================================================================================="""

    # Parse each fasta file, store names of each for subsequent msf parsing
    seqs = []
    for file in fasta_files:
        new_seqs = parse_fasta(file)
        seqs.append(new_seqs)  # Store in nested list to access only relevant fa files for each msf

    # To speed up the testing process, two random sequences from each group of sequences will be
    # selected for pairwise alignment
    random_seqs = []
    for seq in seqs:
        random_seqs.append(sample(seq, 2))

    # Parse each msf file
    for i, file in enumerate(msf_files):
        ref_align = file.rsplit('/', maxsplit=1)[-1].strip('.msf')  # Get name of ref alignment

        # Get corresponding fasta files for this msf file
        # CHANGE 'random_seqs[i]' TO 'seqs[i]' TO USE ALL SEQUENCES
        sequences = random_seqs[i]

        # Only want to align each sequence to every other sequence once
        file_count = 0  # Keep track of number of files for naming purposes
        for i, seq in enumerate(sequences):
            loop_count = i  # Keep track of number of loops so no repeats occur
            while loop_count != len(sequences):
                if seq != sequences[loop_count]:  # Don't want to align a sequence to itself

                    # Align sequences with local programs
                    args = (f'-file1 bb_data/{ref_dir}/{ref_align}/{seq} '
                            f'-file2 bb_data/{ref_dir}/{ref_align}/{sequences[loop_count]} '
                            f'-gopen {-11} '
                            f'-gext {-1} '
                            f'-blosum {45}')
                    os.system(f"python global_align.py {args}")

                    # Embed sequences in this script to save time on loading models
                    args = (f'-file1 bb_data/{ref_dir}/{ref_align}/{seq} '
                            f'-file2 bb_data/{ref_dir}/{ref_align}/{sequences[loop_count]} '
                            f'-gopen {-11} '
                            f'-gext {-1} ')
                    os.system(f"python PEbA_align.py {args}")

                    # Grab alignment from reference MSA
                    seq1, seq2 = seq.split('.')[0], sequences[loop_count].split('.')[0]  # Remove fa
                    align1, align2 = parse_msf(file, seq1, seq2)  # Gather pairwise alignment
                    file_path = f'bb_data/{ref_dir}/{ref_align}/{ref_align}_{file_count}'
                    write_align(align1, align2, seq1, seq2, file_path)  # Write pairwise alignment
                    file_count += 1
                loop_count+=1


def compare_aligns(path):
    """=============================================================================================
    This function takes a directory and compares global and PEbA alignments to the reference
    alignment using t_coffee's aln_compare function. The output is stored in a text file.

    :param path: directory where alignments exist
    ============================================================================================="""

    folders = os.listdir(path)
    for folder in folders:
        files = os.listdir(f'{path}/{folder}')
        if 'global_0.msf' in files:  # Only parse folders with alignments
            global_aligns = []
            peba_aligns = []
            ref_aligns = []
            for file in files: # Add each alignment to list of alignments
                if file.startswith('global'):
                    global_aligns.append(f'{path}/{folder}/{file}')
                if file.startswith('PEbA'):
                    peba_aligns.append(f'{path}/{folder}/{file}')
                if file.startswith('BB'):
                    ref_aligns.append(f'{path}/{folder}/{file}')

            # Sort so that correct alignments are compared
            global_aligns.sort()
            peba_aligns.sort()
            ref_aligns.sort()

            # Call t_coffee to compare global and peba aligns to refs
            for i, ref_align in enumerate(ref_aligns):
                global_align = global_aligns[i]
                peba_align = peba_aligns[i]

                # Get names of alignments for output file
                global_name = global_align.split('/')[-1].split('.')[0]
                peba_name = peba_align.split('/')[-1].split('.')[0]
                os.system(f't_coffee -other_pg aln_compare -al1 {global_align} -al2 {ref_align} > {path}/{folder}/{global_name}_compare.txt')
                os.system(f't_coffee -other_pg aln_compare -al1 {peba_align} -al2 {ref_align} > {path}/{folder}/{peba_name}_compare.txt')


def parse_compare(path):
    """=============================================================================================
    This function takes a directory and compares the global and PEbA alignment comparisons to the
    reference alignment.

    :param path: directory where alignments exist
    ============================================================================================="""

    folders = os.listdir(path)
    for folder in folders:
        files = os.listdir(f'{path}/{folder}')
        if 'global_0.msf' in files:  # Only parse folders with alignments
            global_compares = []
            peba_compares = []
            for file in files:
                if 'compare' in file:
                    if file.startswith('global'):
                        global_compares.append(f'{path}/{folder}/{file}')
                    if file.startswith('peba'):
                        peba_compares.append(f'{path}/{folder}/{file}')

            # Sort so that correct comparisons are compared
            global_compares.sort()
            peba_compares.sort()

            # For both compare files, store their third line in a csv
            global_lines = []
            peba_lines = []
            for i, compare in enumerate(global_compares):
                with open(compare, 'r', encoding='utf8') as file:
                    lines = file.readlines()
                    third_line = lines[2].split()
                    third_line = f'{", ".join(third_line[0:4])}, {third_line[-1].strip("]")}\n'  # Columns used not important
                    global_lines.append(third_line)
                    os.remove(compare)
                with open(peba_compares[i], 'r', encoding='utf8') as file:
                    lines = file.readlines()
                    third_line = lines[2].split()
                    third_line = f'{", ".join(third_line[0:4])}, {third_line[-1].strip("]")}\n'  # Columns used not important
                    peba_lines.append(third_line)
                    os.remove(peba_compares[i])
            with open(f'{path}/{folder}/compare.csv', 'w', encoding='utf8') as file:
                for i, line in enumerate(global_lines):
                    file.write(line)
                    file.write(peba_lines[i])


def graph_compare(path):
    """=============================================================================================
    This function takes a directory and graphs the differences between the global and PEbA
    alignments compared the reference alignment.
    ============================================================================================="""

    global_sim = []
    peba_sim = []
    folders = os.listdir(path)
    for folder in folders:
        files = os.listdir(f'{path}/{folder}')
        if 'compare.csv' in files:  # Only parse folders with comparison csvs
            with open(f'{path}/{folder}/compare.csv', 'r', encoding='utf8') as file:
                for line in file:
                    line = line.split(',')
                    if line[0].startswith('global'):
                        global_sim.append(float(line[3]))
                    if line[0].startswith('PEbA'):
                        peba_sim.append(float(line[3]))

    # Graph the the similarity scores on the same plot to compare
    fig = plt.figure()
    ax = fig.add_subplot()
    sim_diff = []
    for i, gsim in enumerate(global_sim):
        sim_diff.append(peba_sim[i]-gsim)
    ax.scatter(list(range(1, len(sim_diff) + 1)), sim_diff)
    ax.set_title('Difference in Total Column Score PEbA vs BLOSUM')
    ax.set_xlabel('Alignment Number')
    ax.set_ylabel('Similarity Difference')
    ax.axhline(0, color='black')
    plt.savefig(f'{path}/compare.png')


def main():
    """=============================================================================================
    This function calls parse_ref_folder to get lists of all the msf and tfa fasta files in the
    reference directory of interest. It then calls parse_align_files to parse each tfa file and msf
    file, while also aligning each pariwise comparison of fasta sequences. These alignments are
    compared using t_coffee's 'aln_compare' function and the results are parsed and graphed.
    ============================================================================================="""

    # Parse reference folder of interest
    path = 'BAliBASE_R1-5/bb3_release/RV11'
    ref_dir = path.rsplit('/', maxsplit=1)[-1]  # Get last directory in path
    msf_files, fasta_files = parse_ref_folder(path)
    if not os.path.isdir(f'bb_data/{ref_dir}'):
        os.makedirs(f'bb_data/{ref_dir}')

    # Sort each list of files to ensure they match up for msf parsing
    msf_files.sort()
    fasta_files.sort()
    parse_align_files(msf_files, fasta_files, ref_dir)

    # Compare alignments using t_coffee
    path = 'bb_data/RV11'
    compare_aligns(path)
    parse_compare(path)
    graph_compare(path)


if __name__ == '__main__':
    main()
