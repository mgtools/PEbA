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
    parser.add_argument('-m', type=str, default='data/alignments/local_blosum', help='path to directory containing alignments')
    parser.add_argument('-r', type=str, default='RV913', help='reference to compare to')
    args = parser.parse_args()

    log_filename = f'{args.m}//compare_aligns_{args.r}.log'  #pylint: disable=C0103
    logging.basicConfig(filename=log_filename, filemode='w',
                     level=logging.INFO, format='%(message)s')

    # Compare alignments in directory to corresponding reference alignments
    for direc in os.listdir(f'{args.m}/{args.r}'):
        for align in os.listdir(f'{args.m}/{args.r}/{direc}'):

            # Call compute_score
            score = subprocess.getoutput(
                f'python compute_score.py '
                f'-align1 data/alignments/refs/{args.r}/{direc}/{align} '
                f'-align2 {args.m}/{args.r}/{direc}/{align}')

            # Write results to log
            logging.info('%s, %s, %s', direc, align, score)


if __name__ == '__main__':
    main()
