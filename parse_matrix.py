"""================================================================================================
This script parses a substituation matrix stored in a text file and returns a dictionary of the
matrix values. The matrix must be in the following format: the first line contains the amino acids
in the matrix, the first column contains the amino acids in the matrix, and the remaining values are
the substitution scores.

Ben Iovino  02/15/23   VecAligns
================================================================================================"""


def parse_matrix(file):
    """=============================================================================================
    This function takes a text file containing a substitution matrix and parses it into a dict.

    :param file: text file containing substitution matrix
    ============================================================================================="""

    # Open file and read lines
    with open(file, 'r', encoding='utf8') as f:
        lines = f.readlines()

    # Parse first line to get amino acids
    amino_acids = lines[0].split()

    # Parse remaining lines to get substitution scores
    subs_matrix = {}
    for line in lines[1:]:
        line = line.split()
        for j, score in enumerate(line[1:]):
            subs_matrix[f'{line[0]}{amino_acids[j]}'] = int(score)

    print(subs_matrix)


def main():
    """=============================================================================================
    ============================================================================================="""

    file = '/home/ben/Documents/PFASUM60.txt'
    parse_matrix(file)


if __name__ == '__main__':
    main()
