"""================================================================================================
This script parses MSF files to return alignments between particular sequences.

Ben Iovino  02/01/23   VecAligns
================================================================================================"""

from utility import write_align


def parse_msf(filename, id1, id2):
    """=============================================================================================
    This function accepts an MSF file and returns the pairwise alignment between two sequences.

    :param filename: name of file
    :param id1: first sequence id
    :param id2: second sequence id
    return: msf file with particular alignments
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
    seq1 = ''.join(seq1)
    seq2 = ''.join(seq2)
    return seq1, seq2


def main():
    """=============================================================================================
    Parse MSF and write to new MSF.
    ============================================================================================="""

    filename = '/home/ben/Desktop/BB30004.msf'
    id1 = '1buc_A'
    id2 = '1ivh_A'
    seq1, seq2 = parse_msf(filename, id1, id2)
    write_align(seq1, seq2, id1, id2, f'ref_alignment{filename.rsplit("/", maxsplit=-1)[-1]}')

if __name__ == '__main__':
    main()
