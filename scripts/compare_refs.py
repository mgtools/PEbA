"""Compares pairwise alignments produced by different methods to balibase
reference alignments.

__author__ = "Ben Iovino"
__date__ = 09/22/23
"""

import argparse
import logging
import os
import subprocess


def main():
    """Main
    """

    parser = argparse.ArgumentParser()
    parser.add_argument('-m', type=str, default='data/alignments/vcmsa')
    parser.add_argument('-r', type=str, default='RV913', help='reference to compare to')
    parser.add_argument('-s', type=str, default='sp', help='score to compare (sp/f1)')
    args = parser.parse_args()

    log_filename = f'{args.m}/{args.r}_{args.s}.log'  #pylint: disable=C0103
    logging.basicConfig(filename=log_filename, filemode='w',
                     level=logging.INFO, format='%(message)s')

    # Compare alignments in directory to corresponding reference alignments
    direcs = sorted(os.listdir(f'{args.m}/{args.r}'))
    for direc in direcs:
        for align in os.listdir(f'{args.m}/{args.r}/{direc}'):

            # Call compute_score
            score = subprocess.getoutput(
                f'python scripts/compute_score.py '
                f'-align1 data/alignments/refs/{args.r}/{direc}/{align} '
                f'-align2 {args.m}/{args.r}/{direc}/{align} '
                f'-score {args.s}')

            # Write results to log
            logging.info('%s, %s, %s', direc, align, score)


if __name__ == '__main__':
    main()
