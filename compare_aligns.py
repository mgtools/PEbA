"""================================================================================================
This script takes multiple sequence alignments (MSAs) from BAliBASE and parses each one to get each
pairwise (PW) alignment. It takes the FASTA sequences from these MSAs and aligns them using
specified methods. The results from these methods are compared to the reference alignment using
compute_pra.py to get a similarity score. A scatterplot is then generated to show the difference
between the similarity score of the two specified methods.

Ben Iovino  03/07/23   VecAligns
================================================================================================"""

import os
import sys
import argparse
from time import strftime
from random import sample
from Bio import SeqIO
import matplotlib.pyplot as plt
import tensorflow as tf
#from dedal import infer  #pylint: disable=E0401
from utility import parse_fasta, write_align


def dedal(model, seq1, seq2):
    """=============================================================================================
    Runs the DEDAL model to get a pairwise alignment between two proteins.

    :param model: DEDAL model
    :param seq1: First protein sequence
    :param seq2: Second protein sequence
    :return: Alignment object
    ============================================================================================="""

    inputs = infer.preprocess(seq1, seq2)
    align_out = model(inputs)
    output = infer.expand(
        [align_out['sw_scores'], align_out['paths'], align_out['sw_params']])
    output = infer.postprocess(output, len(seq1), len(seq2))
    alignment = infer.Alignment(seq1, seq2, *output)
    return alignment


def parse_align(file):
    """=============================================================================================
    This function gathers the truncated sequences and their beggining and ending indices from the
    alignment file.

    :param file: alignment file
    return: truncated sequences and their positions
    ============================================================================================="""

    # Gather beginning position, truncated sequence, and ending position
    tseq1 = [0, '', 0]
    tseq2 = [0, '', 0]
    with open(file, 'r', encoding='utf8') as f:
        count = 0
        for line in f:
            split_line = line.split()
            if count == 0:  # First line contains first sequence
                tseq1[0] = int(split_line[0])
                tseq1[1] = split_line[1].replace('-', '.')
                tseq1[2] = int(split_line[2])
            if count == 2:  # Third line contains second sequence
                tseq2[0] = int(split_line[0])
                tseq2[1] = split_line[1].replace('-', '.')
                tseq2[2] = int(split_line[2])
            count += 1

    return tseq1, tseq2


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


def parse_fasta_ca(filename, bb_dir):
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


def write_align_ca(seq1, seq2, id1, id2, path):
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


def run_PEbA(bb_dir, ref_align, seq1, ref_dir, seq2, matrix, gopen, gext, encoder):
    """=============================================================================================
    This function accepts a list of args and runs PEbA on two sequences. Args explained in main().
    ============================================================================================="""

    if encoder == 'ProtT5':
        embed = 'prot_t5_embed'
    if encoder == 'ESM2':
        embed = 'esm2_t36_embed'

    args = (f'-file1 {bb_dir}/{ref_align}/{seq1} '
            f'-file2 {bb_dir}/{ref_align}/{seq2} '
            f'-embed1 {embed}/{ref_dir}/{ref_align}/{seq1.split(".")[0]}.txt '
            f'-embed2 {embed}/{ref_dir}/{ref_align}/{seq2.split(".")[0]}.txt '
            f'-matrix {matrix} '
            f'-gopen {gopen} '
            f'-gext {gext} '
            f'-encoder {encoder}')

    print(f'{strftime("%H:%M:%S")} PEbA: {ref_align}/{seq1} and {ref_align}/{seq2}\n',
                           file=sys.stdout)
    os.system(f"python local_PEbA.py {args}")


def run_matrix(bb_dir, ref_align, seq1, seq2, gopen, gext, matrix, value):
    """=============================================================================================
    This function accepts a list of args and runs substitution matrix scoring based alignment on
    two sequences. Args explained in main().
    ============================================================================================="""

    args = (f'-file1 {bb_dir}/{ref_align}/{seq1} '
            f'-file2 {bb_dir}/{ref_align}/{seq2} '
            f'-gopen {gopen} '
            f'-gext {gext} '
            f'-matrix {matrix} '
            f'-score {value}')

    print(f'{strftime("%H:%M:%S")} MATRIX: {ref_align}/{seq1} and {ref_align}/{seq2}\n',
                           file=sys.stdout)
    os.system(f"python local_MATRIX.py {args}")


def dedal_run(bb_dir, ref_align, seq1, seq2, dedal_model):
    """=============================================================================================
    This function accepts a list of args and runs dedal on two sequences. Args explained in main().
    ============================================================================================="""

    print(f'{strftime("%H:%M:%S")} DEDAL: {ref_align}/{seq1} and {ref_align}/{seq2}\n',
                           file=sys.stdout)
    
    seq1, id1 = parse_fasta(f'{bb_dir}/{ref_align}/{seq1}')
    seq2, id2 = parse_fasta(f'{bb_dir}/{ref_align}/{seq2}')

    alignment = dedal(dedal_model, seq1, seq2)
    with open('dedal_output.txt', 'w', encoding='utf8') as f:
        f.write(str(alignment))

    # Parse alignment file, match to original sequences, and write to msf file
    tseq1, tseq2 = parse_align('dedal_output.txt')
    write_align(tseq1[1], tseq2[1], id1, id2, 'DEDAL', 'None', 'None', f'{bb_dir}/{ref_align}/{seq1}')


def parse_align_files(msf_files, fasta_files, bb_dir, methods, samp, dedal_model):
    """=============================================================================================
    This function accepts lists of two sets of files and a directory to place them in where they
    are parsed correspondingly. As they are parsed, they are also aligned using two different
    alignment methods.

    :param msf_files: list of msf files
    :param fasta_files: list of fasta files
    :param bb_dir: directory to place files in
    :param methods: dictionary containing methods and their parameters
    :param samp: number of PW alignments to sample from each MSA
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
        
        # Set a sample size for the PW aligns for testing purposes
        if len(pairwise_aligns) > samp:
            pairwise_aligns = sample(pairwise_aligns, samp)

        # For the selected pairs, align them using local programs
        file_count = 0
        for pair in pairwise_aligns:
            seq1 = pair[0]
            seq2 = pair[1]
            ref_dir = bb_dir.split('/')[1]  # Embedding dirs don't change like results dirs

            for method, pars in methods.items():  #pylint: disable=W0612
                if pars[0] == 'PEbA':
                    run_PEbA(bb_dir, ref_align, seq1, ref_dir, seq2, pars[2], pars[3], pars[4], pars[5])
                if pars[0] == 'matrix':
                    run_matrix(bb_dir, ref_align, seq1, seq2, pars[3], pars[4], pars[1], pars[2])
                if pars[0] == 'dedal':
                    dedal_run(bb_dir, ref_align, seq1, seq2, dedal_model)

            # Grab pairwise alignment from reference MSA
            seq1, seq2 = seq1.split('.')[0], seq2.split('.')[0]  # Remove fa
            align1, align2 = parse_msf(file, seq1, seq2)  # Gather pairwise alignment
            file_path = f'{bb_dir}/{ref_align}/{ref_align}_{file_count}'
            write_align_ca(align1, align2, seq1, seq2, file_path)  # Write pairwise alignment
            file_count += 1


def compare_aligns(path):
    """=============================================================================================
    This function takes a directory and compares the two alignment methods to the reference
    alignment using compute_pra.py The output is stored in a text file.

    :param path: directory where alignments exist
    ============================================================================================="""

    folders = os.listdir(path)
    for folder in folders:
        files = os.listdir(f'{path}/{folder}')
        method1_aligns = []
        method2_aligns = []
        ref_aligns = []
        for file in files: # Add each alignment to list of alignments

            # Even numbered alignments are from method1, odd numbered from method2
            if file.startswith('alignment'):
                align_number = int(file.split('_')[1].strip('.msf'))
                if align_number == 0 or align_number % 2 == 0:
                    method1_aligns.append(f'{path}/{folder}/{file}')
                else:
                    method2_aligns.append(f'{path}/{folder}/{file}')

            # Reference alignments
            if file.startswith('BB'):
                ref_aligns.append(f'{path}/{folder}/{file}')
            if file.startswith('BOX'):
                ref_aligns.append(f'{path}/{folder}/{file}')

        # Sort so that correct alignments are compared
        method1_aligns = sorted(method1_aligns, key=lambda x: int(x.split('_')[2].strip('.msf')))
        method2_aligns = sorted(method2_aligns, key=lambda x: int(x.split('_')[2].strip('.msf')))
        ref_aligns = sorted(ref_aligns, key=lambda x: int(x.split('_')[2].strip('.msf')))

        # Call compute_pra.py to compare global and peba aligns to refs
        for i, ref_align in enumerate(ref_aligns):
            method1_align = method1_aligns[i]
            method2_align = method2_aligns[i]

            # Get names of alignments for output file and call t-coffee's aln_compare
            method1_name = method1_align.split('/')[-1].strip('.msf')
            method2_name = method2_align.split('/')[-1].strip('.msf')
            print(f'Comparing {method1_name} and {method2_name} to {ref_align}')
            os.system(f'python compute_pra.py -align1 {ref_align} -align2 {method1_align} > '
                      f'{path}/{folder}/method1_{method1_name}_compare.txt')
            os.system(f'python compute_pra.py -align1 {ref_align} -align2 {method2_align} > '
                      f'{path}/{folder}/method2_{method2_name}_compare.txt')


def parse_compare(path):
    """=============================================================================================
    This function takes a directory and compares the two alignment methods to the reference.

    :param path: directory where alignments exist
    ============================================================================================="""

    folders = os.listdir(path)
    for folder in folders:
        files = os.listdir(f'{path}/{folder}')
        method1_compares = []
        method2_compares = []
        for file in files:
            if file.startswith('method1'):
                method1_compares.append(f'{path}/{folder}/{file}')
            if file.startswith('method2'):
                method2_compares.append(f'{path}/{folder}/{file}')

        # Sort so that correct comparisons are compared
        method1_compares = sorted(method1_compares, key=lambda x: int(x.split('_')[3]))
        method2_compares = sorted(method2_compares, key=lambda x: int(x.split('_')[3]))

        # For both compare files, join only numerical values to lists
        method1_vals = []
        method2_vals = []
        for i, compare in enumerate(method1_compares):
            with open(compare, 'r', encoding='utf8') as file1:
                vals = file1.readline().split()
                method1_vals.append(f"{', '.join(vals[1::2])}\n")
                os.remove(compare)
            with open(method2_compares[i], 'r', encoding='utf8') as file2:
                vals = file2.readline().split()
                method2_vals.append(f"{', '.join(vals[1::2])}\n")
                os.remove(method2_compares[i])

        # Write the comma separated values to a csv
        with open(f'{path}/{folder}/compare.csv', 'w', encoding='utf8') as file3:
            for i, val in enumerate(method1_vals):
                file3.write(val)
                file3.write(method2_vals[i])


def graph_compare(path, methods):
    """=============================================================================================
    This function takes a directory and makes two graphs. The first graphs the differences between 
    the two methods compared the reference alignment and the second graphs the sim scores against
    each other with the better method2 scores on the left side of the diagonal line and the better 
    method1 scores on the right side of the diagonal line.

    :param path: directory where compare.csv files exist
    :param methods: dict of methods used and their parameters
    ============================================================================================="""

    # Get the similarity scores from the compare files
    method1_sim = []
    method2_sim = []
    folders = os.listdir(path)
    for folder in folders:
        with open(f'{path}/{folder}/compare.csv', 'r', encoding='utf8') as file:
            for i, line in enumerate(file):
                line = line.split(',')
                if i == 0 or i % 2 == 0:
                    method1_sim.append(float(line[0]))
                if i % 2 != 0:
                    method2_sim.append(float(line[0]))

    # Average the similarity scores to put on the graph
    m1_avg = round(sum(method1_sim)/len(method1_sim), 1)
    m2_avg = round(sum(method2_sim)/len(method2_sim), 1)

    # Get names and parameters of methods used for titles of graph
    titles = []
    for method, pars in methods.items():  #pylint: disable=W0612
        if pars[0] == 'PEbA':
            titles.append(f'PEbA_{pars[5]}')
        if pars[0] == 'matrix':
            titles.append(f'{pars[1]}{pars[2]}')
        if pars[0] == 'dedal':
            titles.append('DEDAL')

    # Graph the difference between similarity scores for each alignment
    fig = plt.figure()
    ax = fig.add_subplot()
    sim_diff = []
    for i, mat_sim in enumerate(method1_sim):
        sim_diff.append(mat_sim-method2_sim[i])
    ax.scatter(list(range(1, len(sim_diff) + 1)), sim_diff)
    ax.set_title(f'Difference in {titles[0]} (Avg={m1_avg}) vs. {titles[1]} (Avg={m2_avg})')
    ax.set_xlabel('Alignment Number')
    ax.set_ylabel('Similarity Difference')
    ax.set_ylim(-20, 80)
    ax.axhline(0, color='black')
    plt.savefig(f'{path}/differences.png')

    # Graph the similarity scores against each other, better score on either side of diag line
    fig = plt.figure()
    ax = fig.add_subplot()
    method1_scores = []
    method2_scores = []
    for i, pra in enumerate(method2_sim):  # pra = percent residues aligned
        if pra < method1_sim[i]:  # Make method1 prca the x value if greater than method2 pra
            method1_scores.append([method1_sim[i], pra])
        else:  # Make method2 pra the x value if greater than method1 pra
            method2_scores.append([pra, method1_sim[i]])
    ax.scatter([i[0] for i in method2_scores], [i[1] for i in method2_scores], color='blue')
    ax.scatter([i[1] for i in method1_scores], [i[0] for i in method1_scores], color='red')
    ax.set_title(f'{titles[0]} Alignment (Avg={m1_avg}) vs. {titles[1]} Alignment (Avg={m2_avg})')
    ax.set_xlabel(f'PRA {titles[1]} (%)')
    ax.set_ylabel(f'PRA {titles[0]} (%)')
    plt.plot([0, 100], [0, 100], color='black')
    plt.savefig(f'{path}/comparison.png')


def main():
    """=============================================================================================
    This function calls parse_ref_folder to get lists of all the msf and tfa fasta files in the
    reference directory of interest. It then calls parse_align_files to parse each tfa file and msf
    file, while also aligning each pariwise comparison of fasta sequences. These alignments are
    compared using compute_pra.py and the results are parsed and graphed.

    Listed below are the current supported methods for comparison. You can put either method for
    -method1 and -method2, the only major difference is how the results are presented in the
    final scatterplot (method1 goes on the y-axis and method2 goes on the x-axis).

    If using pre-embedded sequences, you can use the -encoder to specify which encoder was used.
    Place them into the VecAligns/ folder and name the folder 'x_embed'.

    PEbA with ProtT5 encodings: -compare1 PEbA -encoder ProtT5
    PEbA with ESM2 encodings: -compare1 PEbA -encoder ESM2
    Substitution matrix with BLOSUM: -compare1 matrix -matrix1 blosum - value1 <integer>
    Substitution matrix with PFASUM: -compare1 matrix -matrix1 pfasum - value1 60  (60 only)
    dedal: -compare1 dedal
    ============================================================================================="""

    parser = argparse.ArgumentParser()
    parser.add_argument('-path', type=str, default='BAliBASE_R1-5/bb3_release/RV11', help='Ref direc')
    parser.add_argument('-sample', type=int, default=1, help='MSA sample size')
    parser.add_argument('-method1', type=str, default='PEbA', help='First method for comparison')
    parser.add_argument('-matrix1', type=str, default='blosum', help='Substution matrix')
    parser.add_argument('-value1', type=int, default=45, help='Sub matrix value')
    parser.add_argument('-gopen1', type=float, default=-11, help='Gap open score')
    parser.add_argument('-gext1', type=float, default=-1, help='Gap ext score')
    parser.add_argument('-encoder1', type=str, default='ProtT5', help='Model used for embeddings')
    parser.add_argument('-method2', type=str, default='matrix', help='Second method for comparison')
    parser.add_argument('-matrix2', type=str, default='blosum', help='Substution matrix')
    parser.add_argument('-value2', type=int, default=45, help='Sub matrix value')
    parser.add_argument('-gopen2', type=float, default=-11, help='Gap open score')
    parser.add_argument('-gext2', type=float, default=-1, help='Gap ext score')
    parser.add_argument('-encoder2', type=str, default='ESM2', help='Model used for embeddings')
    args = parser.parse_args()

    # Places methods and their arguments into dict for easy access and readability
    methods = {'method1': [args.method1, args.matrix1, args.value1, args.gopen1, args.gext1, args.encoder1],
               'method2': [args.method2, args.matrix2, args.value2, args.gopen2, args.gext2, args.encoder2]}

    # Get directory of reference alignments i.e. 'RV11'
    ref_dir = args.path.rsplit('/', maxsplit=1)[-1]

    # Create unique directory for results, this allows for parallel runs of the script i.e. 'bb_data0'
    bb_ct = 0
    for direc in os.listdir():
        if direc.startswith('bb_data'):
            bb_ct += 1
    bb_dir = f'bb_data{bb_ct}/{ref_dir}'
    os.makedirs(bb_dir)

    # Load DEDAL if necessary
    if args.method1 == 'dedal' or args.method2 == 'dedal':
        dedal_model = tf.saved_model.load('dedal_3')
    else:
        dedal_model = 'n'

    # Parse reference folder of interest
    print(f'{strftime("%H:%M:%S")} Parsing and computing alignments...\n', file=sys.stdout)
    msf_files, fasta_files = parse_ref_folder(args.path)

    # Sort each list of files to ensure they match up for msf parsing
    msf_files.sort()
    fasta_files.sort()
    parse_align_files(msf_files, fasta_files, bb_dir, methods, args.sample, dedal_model)

    # Compare alignments to get PRA and graph results
    print(f'{strftime("%H:%M:%S")} Comparing alignments...\n', file=sys.stdout)
    compare_aligns(bb_dir)
    parse_compare(bb_dir)
    graph_compare(bb_dir, methods)
    print(f'{strftime("%H:%M:%S")} Program Complete!\n', file=sys.stdout)


if __name__ == '__main__':
    main()
