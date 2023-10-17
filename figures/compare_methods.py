"""Compares the performance of different methods on the balibase reference alignments.

__author__ = "Ben Iovino"
__date__ = 09/25/23
"""

import argparse
import os
import matplotlib.pyplot as plt


def get_scores(method: str, ref: str, score: str) -> dict:
    """Returns a list of scores for each alignment in the given directory

    :param method: path to directory containing alignments
    :param ref: particular reference
    :param score: score to compare (sp/f1)
    :return dict: dict where key is alignment name and value is sim score to ref align
    """

    scores = {}
    with open(f'{method}/{ref}_{score}.log', 'r', encoding='utf8') as file:
        for line in file:
            line = line.split()
            scores[line[1]] = float(line[3])

    return scores


def graph_compare(m1: str, m2: str, ref: str, metric: str, m1_scores: list, m2_scores: list):
    """Plots a graph comparing the scores of the two methods to the reference

    :param m1: path to directory containing alignments for method 1
    :param m2: path to directory containing alignments for method 2
    :param ref: particular reference
    :param metric: metric used to compare alignments
    :param m1_scores: list of scores for method 1
    :param m2_scores: list of scores for method 2
    """

    # Average the similarity scores to put on the graph, rounding to same decimals as in paper
    m1_avg = round(sum(m1_scores)/len(m1_scores), 3)
    m2_avg = round(sum(m2_scores)/len(m2_scores), 3)
    metric = metric.upper()

    # Get method names
    m1_title = m1.split('/')[-1]
    m2_title = m2.split('/')[-1]

    fig = plt.figure()
    ax = fig.add_subplot()
    method1_scores = []
    method2_scores = []
    for i, score in enumerate(m2_scores):
        if score < m1_scores[i]:  # Make method1 score the x value if greater than method2 score
            method1_scores.append([m1_scores[i], score])
        else:  # Make method2 score the x value if greater than method1 score
            method2_scores.append([score, m1_scores[i]])
    ax.scatter([i[0] for i in method2_scores], [i[1] for i in method2_scores],
                color='blue', alpha=0.5)
    ax.scatter([i[1] for i in method1_scores], [i[0] for i in method1_scores],
                color='red', alpha=0.5)
    ax.set_title(f'{m1_title} Alignments vs. {m2_title} Alignments in {ref}')
    ax.set_xlabel(f'{m2_title}')
    ax.set_ylabel(f'{m1_title}')
    ax.legend([f'{m2_title} Avg {metric}: {m2_avg}', f'{m1_title} Avg {metric}: {m1_avg}'])
    plt.plot([0, 1], [0, 1], color='black')
    plt.savefig(f'figures/graphs/{m1_title}-{m2_title}-{ref}-comparison-{metric}.png')


def main():
    """Main
    """

    parser = argparse.ArgumentParser()
    parser.add_argument('-m1', type=str, default='data/alignments/local_peba_t5')
    parser.add_argument('-m2', type=str, default='data/alignments/local_blosum')
    parser.add_argument('-r', type=str, default='RV913')
    parser.add_argument('-s', type=str, default='sp')
    args = parser.parse_args()

    if not os.path.isdir('figures/graphs'):
        os.makedirs('figures/graphs')

    # Get comparison scores for each set of alignments
    m1_dict = get_scores(args.m1, args.r, args.s)
    m2_dict = get_scores(args.m2, args.r, args.s)

    # Only take scores from alignments that are in both methods
    m1_scores = [m1_dict[align] for align in m1_dict if align in m2_dict]
    m2_scores = [m2_dict[align] for align in m2_dict if align in m1_dict]
    graph_compare(args.m1, args.m2, args.r, args.s, m1_scores, m2_scores)


if __name__ == '__main__':
    main()
