"""================================================================================================
This script uses t_coffee to compare the pairwise alignments between pairs of sequences created
from global NW and PEbA alignments to the reference pairwise alignments. The output is stored in
a text file for later analysis.

Ben Iovino  02/06/23   VecAligns
================================================================================================"""

import os
import matplotlib.pyplot as plt


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
                os.system(f't_coffee -other_pg aln_compare -al1 {global_align} -al2 {ref_align} > {path}/{folder}/global_{global_name}_compare.txt')
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
    ax.set_title('Difference in Similarity Scores Between PEbA and NW')
    ax.set_xlabel('Alignment Number')
    ax.set_ylabel('Similarity Difference')
    ax.legend('PEbA vs. Global')
    ax.axhline(0, color='black')
    plt.savefig(f'{path}/compare.png')


def main():
    """=============================================================================================
    Sets path and calls functions.
    ============================================================================================="""

    path = 'bb_data/RV11'
    compare_aligns(path)
    parse_compare(path)
    graph_compare(path)

if __name__ == '__main__':
    main()
