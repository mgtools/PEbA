"""Returns a table of the scores for the desired method. Requires log files from
compare_refs.py to already exist and be in the correct directory.

__author__ = "Ben Iovino"
__date__ = 10/10/23
"""

import os
import argparse
import pandas as pd


def parse_ref(method: str, bucket: str) -> dict:
    """Returns a dict of ref alignments and each alignment's sim score OR length and SP score.

    :param method: directory containing alignments
    :param bucket: value of interest for bucketing (id or len)
    :return dict: dict where key is reference name and value is a list of tuples
    """

    # For each log file in directory, get sim score OR ref align length and SP score
    refs = {}
    for file in os.listdir(method):
        if file.endswith('.log'):
            refs[file.split('_')[0]] = []
            with open(f'{method}/{file}', 'r', encoding='utf8') as lfile:
                for line in lfile:
                    line = line.split()

                    # SP on line[3], length on line[5], sim on line[9]
                    if bucket == 'len':
                        refs[file.split('_')[0]].append((line[5], line[3]))
                    if bucket == 'id':
                        refs[file.split('_')[0]].append((line[9], line[3]))

    return refs


def parse_scores_id(scores: dict):  #\\NOSONAR
    """Prints a pandas table of the average SP score for each bucket of pairwise id

    :param scores: dict where key is reference name and value is a list of tuples
    """

    # Create pandas table with columns for each bucket and rows for each reference
    table = pd.DataFrame(columns=['0-9', '10-19', '20-29', '30-39', '40-49',
                                    '50-59', '60-69', '70-79', '80-89', '90-99'],
                            index=scores.keys())

    # For each reference in dict, add SP score to appropriate bucket
    for ref in scores:

        # Create dict of lists for each bucket
        id_dict = {9: [], 19: [], 29: [], 39: [], 49: [],
                59: [], 69: [], 79: [], 89: [], 99: []}

        for score in scores[ref]:

            # Add SP score to appropriate bucket
            for bucket in id_dict:  #pylint: disable=C0206
                if float(score[0])*100 <= bucket:
                    id_dict[bucket].append(float(score[1]))
                    break

        # Get average SP score for each bucket
        for bucket, score_list in id_dict.items():
            if len(score_list) <= 10:
                id_dict[bucket] = 0
                continue
            id_dict[bucket] = round(sum(score_list)/len(score_list), 2)

        # Add row to table, each value in dict is a column
        for i, bucket in enumerate(id_dict):
            table.iloc[table.index.get_loc(ref), i] = id_dict[bucket]

    # Print table
    print(table)


def parse_scores_len(scores: dict):  #\\NOSONAR
    """Prints a pandas table of the average SP score for each bucket of length

    :param scores: dict where key is reference name and value is a list of tuples
    """

    # Create pandas table with columns for each bucket and rows for each reference
    table = pd.DataFrame(columns=['0-499', '500-999', '1000-1499', '1500-1999', '2000-2499'],
                            index=scores.keys())

    # For each reference in dict, add SP score to appropriate bucket
    for ref in scores:

        # Create dict of lists for each bucket
        len_dict = {499: [], 999: [], 1499: [], 1999: [], 2499: []}

        for score in scores[ref]:

            # Add SP score to appropriate bucket
            for bucket in len_dict:  #pylint: disable=C0206
                if int(score[0]) <= bucket:
                    len_dict[bucket].append(float(score[1]))
                    break

        # Get average SP score for each bucket
        for bucket, score_list in len_dict.items():
            if len(score_list) <= 10:
                len_dict[bucket] = 0
                continue
            len_dict[bucket] = round(sum(score_list)/len(score_list), 2)

        # Add row to table, each value in dict is a column
        for i, bucket in enumerate(len_dict):
            table.iloc[table.index.get_loc(ref), i] = len_dict[bucket]

    # Print table
    print(table)


def main():
    """Main
    """

    parser = argparse.ArgumentParser()
    parser.add_argument('-d', type=str, default='data/alignments/local_blosum')
    parser.add_argument('-t', type=str, default='len', help='pairwise id (id) or length (len)')
    args = parser.parse_args()

    scores = parse_ref(args.d, args.t)
    scores = dict(sorted(scores.items()))
    if args.t == 'id':
        parse_scores_id(scores)
    if args.t == 'len':
        parse_scores_len(scores)


if __name__ == "__main__":
    main()
