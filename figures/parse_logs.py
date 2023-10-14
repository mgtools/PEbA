"""Various log parsing.

__author__ = "Ben Iovino"
__date__ = "10/12/23"
"""

import argparse
from datetime import datetime
from Bio import SeqIO


def parse_score_log(logfile: str) -> dict:
    """Returns a dictionary of SP scores from a log file

    :param logfile: path to log file
    :return dict: dict where key is alignment and value is SP score
    """

    # Get SP scores from log file
    scores = {}
    with open(logfile, 'r', encoding='utf8') as lfile:
        for line in lfile:
            line = line.split()
            scores[f'{line[0].strip(",")}/{line[1]}'] = line[3]

    return scores


def compare_zeros(scores1: dict, scores2: dict):
    """Prints some information about the zero values in the two dictionaries

    :param scores1: dictionary of SP scores from first method
    :param scores2: dictionary of SP scores from second method
    """

    zeros1, msas1 = [], {}
    for align, score in scores1.items():
        if float(score) == 0:
            zeros1.append(align)
            msas1[align.split('/')[0]] = msas1.get(align.split('/')[0], []) + [align]
    print(f'There are {len(zeros1)} zero values in the first method')

    zeros2, msas2 = [], {}
    for align, score in scores2.items():
        if float(score) == 0:
            zeros2.append(align)
            msas2[align.split('/')[0]] = msas2.get(align.split('/')[0], []) + [align]
    print(f'There are {len(zeros2)} zero values in the second method')

    count = 0
    for align in zeros1:
        if align not in zeros2:
            count += 1
    print(f'There are {count} zero values in the first method that are not in the second method')


def get_length(align: str, line: str) -> float:
    """Returns the average length of the two sequences in an alignment

    :param align: alignment
    :param line: line from log file
    :return int: average length of the two sequences in the alignment
    """

    seq1, seq2 = align.split('-')
    direc = '/'.join(line.split()[3].split('/')[:2])
    seq1 = SeqIO.read(f'data/sequences/{direc}/{seq1}.fa', 'fasta')
    seq2 = SeqIO.read(f'data/sequences/{direc}/{seq2}.fa', 'fasta')

    return (len(seq1) + len(seq2)) / 2


def parse_align_log(logfile: str) -> dict:
    """Returns a dict of time to align each pair of sequences

    :param logfile: path to log file
    :return dict: dict where key is alignment and value is time to align
    """

    # Get times from log file
    times = {}
    with open(logfile, 'r', encoding='utf8') as lfile:
        prev_pair, prev_time, lengths = '', 0, (0, 0)
        for line in lfile:
            if prev_pair == '':
                prev_pair = line.split()[3].split('/')[-1]
                lengths = get_length(prev_pair, line)
                times[prev_pair] = (lengths, 0)
                prev_time = datetime.strptime(line.split()[1], '%H:%M:%S.%f').time()
                continue

            # Get new time, subtract from previous time, and add to dict
            time = datetime.strptime(line.split()[1], '%H:%M:%S.%f').time()
            delta = datetime.combine(datetime.min, time) - datetime.combine(datetime.min, prev_time)
            times[prev_pair] = (lengths, delta.total_seconds())

            # Update previous pair and time
            prev_pair = line.split()[3].split('/')[-1]
            prev_time = time
            lengths = get_length(prev_pair, line)

    return times


def main():
    """Main
    """

    # Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-sl1', type=str, default='data/alignments/local_peba_t5/RV911_sp.log')
    parser.add_argument('-sl2', type=str, default='data/alignments/local_blosum/RV911_sp.log')
    parser.add_argument('-al1', type=str, default='/home/ben/Desktop/peba_vcmsa_time/RV911/get_aligns_peba.log')
    parser.add_argument('-al2', type=str, default='/home/ben/Desktop/blosum_vcmsa_time/RV911/get_aligns_vcmsa.log')
    args = parser.parse_args()

    # Get SP scores from log files
    #scores1 = parse_score_log(args.sl1)
    #scores2 = parse_score_log(args.sl2)

    # Compare zero values
    #compare_zeros(scores1, scores2)

    # Get times from log files
    times1 = parse_align_log(args.al1)
    print(times1)
    
    # Plot times vs length
    import matplotlib.pyplot as plt
    x, y = [], []
    for (length, time) in times1.values():
        if time < 0:
            continue
        x.append(length)
        y.append(time)

    plt.scatter(x, y)
    plt.xlabel('Length')
    plt.ylabel('Time')
    plt.show()


if __name__ == '__main__':
    main()
