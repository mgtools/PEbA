"""Plots a line graph of each method's average similarity score to the reference

__author__ = "Ben Iovino"
__date__ = 10/10/23
"""

import argparse
import os
import matplotlib.pyplot as plt
import pandas as pd
from make_table import parse_ref, parse_scores


def get_table(method: str, bucket: str) -> pd.DataFrame:
    """Returns a pandas table of the average SP score for each bucket of length or pairwise id

    :param method: directory containing alignments
    :param bucket: value of interest for bucketing (id or len)
    :return pd.DataFrame: pandas table of average SP score for each bucket
    """

    scores = parse_ref(f'data/alignments/{method}', bucket)

    # Combine scores for RV11/RV12 and RV911/RV912/RV913
    new_scores = {}
    for ref, score_list in scores.items():
        if ref in ['RV11', 'RV12']:
            new_scores['RV11/RV12'] = new_scores.get('RV11/RV12', []) + score_list  #\\NOSONAR
        elif ref in ['RV911', 'RV912', 'RV913']:
            new_scores['RV911/RV912/RV913'] = new_scores.get('RV911/RV912/RV913', []) + score_list  #\\NOSONAR
    table = parse_scores(new_scores, bucket)

    return table


def graph_id(ref_values: list, refs: str):
    """Plots performance of each method on each reference

    :param ref_values: list of lists of SP scores for each method
    :param refs: string of reference names
    """

    # X axis labels
    if refs == 'RV11/RV12':
        savefile = 'refs1'
        x_labels = [i*10 for i in range(1, 6)]
    if refs == 'RV911/RV912/RV913':
        savefile = 'refs2'
        x_labels = [i*10 for i in range(1, 9)]


    # Make single line graph
    methods = ['PEbA', 'vcMSA', 'BLOSUM62', 'DEDAL']
    fig = plt.figure()
    ax = fig.add_subplot()
    for i, method in enumerate(methods):
        ax.plot(x_labels, ref_values[i], label=method, marker='o')
    ax.set_xlabel('Pairwise Identity')
    ax.set_ylabel('Average SP Score')
    ax.set_title('Average SP Score for Each Reference')
    ax.legend()
    ax.grid(color='grey', linestyle='-', linewidth=0.25)
    plt.savefig(f'figures/graphs/{savefile}_id.png')


def graph_len(ref_values: list, refs: str):
    """Plots performance of each method on each reference

    :param ref_values: list of lists of SP scores for each method
    :param refs: string of reference names
    """

    # X axis labels
    if refs == 'RV11/RV12':
        savefile = 'refs1'
        x_labels = [(500+i*500) for i in range(0, 3)]
    if refs == 'RV911/RV912/RV913':
        savefile = 'refs2'
        x_labels = [(500+i*500) for i in range(0, 5)]

    # Make single line graph
    methods = ['PEbA', 'vcMSA', 'BLOSUM62', 'DEDAL']
    fig = plt.figure()
    ax = fig.add_subplot()
    for i, method in enumerate(methods):
        ax.plot(x_labels, ref_values[i], label=method, marker='o')
    ax.set_xlabel('Length')
    ax.set_ylabel('Average SP Score')
    ax.set_title('Average SP Score for Each Reference')
    ax.legend()
    ax.grid(color='grey', linestyle='-', linewidth=0.25)
    plt.savefig(f'figures/graphs/{savefile}_len.png')


def main():
    """Main
    """

    parser = argparse.ArgumentParser()
    parser.add_argument('-t', type=str, default='len', help='pairwise id (id) or length (len)')
    args = parser.parse_args()

    if not os.path.isdir('figures/graphs'):
        os.makedirs('figures/graphs')

    methods = ['local_peba_t5', 'vcmsa', 'local_blosum', 'dedal']
    ref1_values = []
    ref2_values = []

    # Get scores for each reference
    for method in methods:
        table = get_table(method, args.t)

        # Table has two rows, one for RV11/RV12 and one for RV911/RV912/RV913
        if args.t == 'id':
            ref1_values.append(table.iloc[0, :].to_list()[:5])
            ref2_values.append(table.iloc[1, :].to_list()[:8])

        # Ignore last column for RV11/12 because no sequences are 1500+ long
        if args.t == 'len':
            ref1_values.append(table.iloc[0, :-2].to_list()[:4])
            ref2_values.append(table.iloc[1, :].to_list()[:5])

    # Graphing id and length differently
    if args.t == 'id':
        graph_id(ref1_values, 'RV11/RV12')
        graph_id(ref2_values, 'RV911/RV912/RV913')
    if args.t == 'len':
        graph_len(ref1_values, 'RV11/RV12')
        graph_len(ref2_values, 'RV911/RV912/RV913')


if __name__ == '__main__':
    main()
