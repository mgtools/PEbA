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
    refs = dict(sorted(refs.items()))

    return refs


def get_table(value: str, scores: dict) -> pd.DataFrame:
    """Returns a pandas dataframe with appropriate column names for the desired value

    :param value: value of interest for bucketing (id or len)
    :param scores: dict where key is reference name and value is a list of tuples
    :return pd.DataFrame: pandas dataframe
    """

    if value == 'len':
        table = pd.DataFrame(columns=['0-499', '500-999', '1000-1499', '1500-1999', '2000-2499'],
                            index=scores.keys())
    if value == 'id':
        table = pd.DataFrame(columns=['0-9', '10-19', '20-29', '30-39', '40-49',
                                    '50-59', '60-69', '70-79', '80-89', '90-99'],
                            index=scores.keys())

    return table


def put_scores(id_dict: dict, scores: dict, ref: str, value: str) -> dict:
    """Adds SP scores to appropriate bucket in dict

    :param id_dict: dict where key is bucket and value is list of SP scores
    :param scores: dict where key is reference name and value is a list of tuples
    :param ref: reference name
    :param value: value of interest for bucketing (id or len)
    :return dict: dict where key is bucket and value is list of SP scores
    """

    for score in scores[ref]:

        # Change value to int if bucketing by length
        if value == 'len':
            score = (int(score[0]), float(score[1]))
        if value == 'id':  # Otherwise change value to a percentage
            score = (float(score[0])*100, float(score[1]))

        # Add SP score to appropriate bucket
        for bucket in id_dict:  #pylint: disable=C0206
            if score[0] <= bucket:
                id_dict[bucket].append(float(score[1]))
                break  # Only add to one bucket

    return id_dict


def parse_scores(scores: dict, value: str) -> pd.DataFrame:  #\\NOSONAR
    """Prints a pandas table of the average SP score for each bucket of pairwise id

    :param scores: dict where key is reference name and value is a list of tuples
    :param value: value of interest for bucketing (id or len)
    :return pd.DataFrame: pandas table of average SP score for each bucket
    """

    # For each reference in dict, add SP score to appropriate bucket
    table = get_table(value, scores)
    for ref in scores:

        # Create dict of lists for each bucket
        if value == 'len':
            value_dict = {499: [], 999: [], 1499: [], 1999: [], 2499: []}
        if value == 'id':
            value_dict = {9: [], 19: [], 29: [], 39: [], 49: [],
                    59: [], 69: [], 79: [], 89: [], 99: []}

        # Add SP score to appropriate bucket
        value_dict = put_scores(value_dict, scores, ref, value)

        # Get average SP score for each bucket
        for bucket, score_list in value_dict.items():
            if len(score_list) <= 10:
                value_dict[bucket] = 0
                continue
            value_dict[bucket] = round(sum(score_list)/len(score_list), 2)

        # Add row to table, each value in dict is a column
        for i, bucket in enumerate(value_dict):
            table.iloc[table.index.get_loc(ref), i] = value_dict[bucket]

    return table


def main():
    """Prints a table of average scores for the desired method
    """

    parser = argparse.ArgumentParser()
    parser.add_argument('-d', type=str, default='data/alignments/local_blosum')
    parser.add_argument('-t', type=str, default='len', help='pairwise id (id) or length (len)')
    args = parser.parse_args()

    scores = parse_ref(args.d, args.t)
    table = parse_scores(scores, args.t)
    print(table)


if __name__ == "__main__":
    main()
