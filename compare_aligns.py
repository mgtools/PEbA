"""Compares pairwise alignments produced by different methods to balibase
reference alignments.

__author__ = "Ben Iovino"
__date__ = 09/22/23
"""

import argparse
import logging
import os

log_filename = 'data/logs/compare_aligns.log'  #pylint: disable=C0103
os.makedirs(os.path.dirname(log_filename), exist_ok=True)
logging.basicConfig(filename=log_filename, filemode='w',
                     level=logging.INFO, format='%(message)s')


def main():
    """Main
    """

    parser = argparse.ArgumentParser()
    parser.add_argument('-m', type=str, default='data/alignments/global_blosum', help='path to directory containing alignments')
    parser.add_argument('-r', type=str, default='RV913', help='reference to compare to')
    args = parser.parse_args()

    # Compare alignments in directory to corresponding reference alignments
    for direc in os.listdir(f'{args.m}/{args.r}'):
        for align in os.listdir(f'{args.m}/{args.r}/{direc}'):

            # Call compute_score
            os.system(f'python compute_score.py -align1 data/alignments/refs/{args.r}/{direc}/{align} '
                      f' -align2 {args.m}/{args.r}/{direc}/{align}')


if __name__ == '__main__':
    main()
