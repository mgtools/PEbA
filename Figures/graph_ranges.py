"""================================================================================================
This script reads csv files and graphs values of interest. The directories that are read in this
script are a result from compare_aligns.py. A set of runs comes from using two methods, i.e. PEbA
and BLOSUM, on each BAliBASE reference used in this project, of which there are five. The five refs
are split into two separate graphs due to differences in sequence length and pairwise identity.

Ben Iovino  04/24/23   VecAligns
================================================================================================"""

import csv
import os
import argparse
import matplotlib.pyplot as plt


def read_csv(filename):
    """=============================================================================================
    This function accepts a csv filename and returns a list of lists with each line being a list.

    :param filename: csv file to be read
    :return: list of lists
    ============================================================================================="""

    with open(filename, 'r', encoding='utf8') as file:
        reader = csv.reader(file)
        data = list(reader)
    return data


def parse_data(data, parse):
    """=============================================================================================
    This function accepts a list of lists and parses for relevant info.

    :param data: list of lists
    :param parse: output to be parsed
    :return: list of lists
    ============================================================================================="""

    # Initialize dict based on type of info to be parsed
    if parse == 'id':
        dict_m1 = {9: [0, 0], 19: [0, 0], 29: [0, 0], 39: [0, 0], 49: [0, 0],
             59: [0, 0], 69: [0, 0], 79: [0, 0], 89: [0, 0], 99: [0, 0]}
        dict_m2 = {9: [0, 0], 19: [0, 0], 29: [0, 0], 39: [0, 0], 49: [0, 0],
             59: [0, 0], 69: [0, 0], 79: [0, 0], 89: [0, 0], 99: [0, 0]}
    elif parse == 'len':
        dict_m1 = {499: [0, 0], 999: [0, 0], 1499: [0, 0], 1999: [0, 0], 2499: [0, 0]}
        dict_m1 = {499: [0, 0], 999: [0, 0], 1499: [0, 0], 1999: [0, 0], 2499: [0, 0]}

    count = 0
    for line in data:

        # Even lines are from method 1
        if count == 0 or count % 2 == 0:
            for key, value in dict_m1.items():
                if parse == 'id':
                    if float(line[3]) <= key:
                        value[0] += float(line[0])/100
                        value[1] += 1
                        break
                elif parse == 'len':
                    if float(line[1]) <= key:
                        value[0] += float(line[0])/100
                        value[1] += 1
                        break

        # Odd lines are from method 2
        else:
            for key, value in dict_m2.items():
                if parse == 'id':
                    if float(line[3]) <= key:
                        value[0] += float(line[0])/100
                        value[1] += 1
                        break
                elif parse == 'len':
                    if float(line[1]) <= key:
                        value[0] += float(line[0])/100
                        value[1] += 1
                        break
        count += 1

    return dict_m1, dict_m2


def parse_run(run, parse):
    """=============================================================================================
    This function accepts a run directory and parses for relevant info.

    :param run: directory to be parsed
    :param parse: information of interest i.e. pairwise identity or alignment length
    :return: dictionary of lists
    ============================================================================================="""

    # Initialize dict based on type of info to be parsed
    if parse == 'id':
        dict_m1 = {9: [0, 0], 19: [0, 0], 29: [0, 0], 39: [0, 0], 49: [0, 0],
             59: [0, 0], 69: [0, 0], 79: [0, 0], 89: [0, 0], 99: [0, 0]}
        dict_m2 = {9: [0, 0], 19: [0, 0], 29: [0, 0], 39: [0, 0], 49: [0, 0],
             59: [0, 0], 69: [0, 0], 79: [0, 0], 89: [0, 0], 99: [0, 0]}
    elif parse == 'len':
        dict_m1 = {499: [0, 0], 999: [0, 0], 1499: [0, 0], 1999: [0, 0], 2499: [0, 0]}
        dict_m1 = {499: [0, 0], 999: [0, 0], 1499: [0, 0], 1999: [0, 0], 2499: [0, 0]}

    # Initialize dictionary
    data = {}
    for ref in os.listdir(f'{run}'):
        for msa in os.listdir(f'{run}/{ref}'):
            if msa.startswith('B'):
                for file in os.listdir(f'{run}/{ref}/{msa}'):
                    if file.endswith('csv'):

                        # Read csv and parse
                        data = read_csv(f'{run}/{ref}/{msa}/{file}')
                        dict2_m1, dict2_m2 = parse_data(data, parse)

                        # Add to running total
                        for key, value in dict_m1.items():
                            value[0] += dict2_m1[key][0]
                            value[1] += dict2_m1[key][1]

                        for key, value in dict_m2.items():
                            value[0] += dict2_m2[key][0]
                            value[1] += dict2_m2[key][1]

    return dict_m1, dict_m2


def main():
    """=============================================================================================
    ============================================================================================="""

    parser = argparse.ArgumentParser()
    parser.add_argument('-p', type=str, default='/home/ben/Desktop/PEbA_Data/Runs/PEBA-BLOSUM')
    parser.add_argument('-r', type=str, default='1')
    parser.add_argument('-t', type=str, default='id')
    args = parser.parse_args()

    # Initialize dict based on type of info to be parsed
    if args.t == 'id':
        dict_m1 = {9: [0, 0], 19: [0, 0], 29: [0, 0], 39: [0, 0], 49: [0, 0],
             59: [0, 0], 69: [0, 0], 79: [0, 0], 89: [0, 0], 99: [0, 0]}
        dict_m2 = {9: [0, 0], 19: [0, 0], 29: [0, 0], 39: [0, 0], 49: [0, 0],
             59: [0, 0], 69: [0, 0], 79: [0, 0], 89: [0, 0], 99: [0, 0]}
    elif args.t == 'len':
        dict_m1 = {499: [0, 0], 999: [0, 0], 1499: [0, 0], 1999: [0, 0], 2499: [0, 0]}
        dict_m1 = {499: [0, 0], 999: [0, 0], 1499: [0, 0], 1999: [0, 0], 2499: [0, 0]}

    # Run directories
    runs = sorted(os.listdir(args.p))
    if args.r == '1':  # '1' stands for references 11 and 12
        runs = runs[:2]
    elif args.r == '9':  # '9' stands for references 911, 912, and 913
        runs = runs[2:]

    # Parse csvs in each run
    for run in runs:
        dict2_m1, dict2_m2 = parse_run(f'{args.p}/{run}', args.t)

        # Add to running total
        for key, value in dict_m1.items():
            value[0] += dict2_m1[key][0]
            value[1] += dict2_m1[key][1]

        for key, value in dict_m2.items():
            value[0] += dict2_m2[key][0]
            value[1] += dict2_m2[key][1]

    # Print results
    print(dict_m1)
    print(dict_m2)


if __name__ == '__main__':
    main()

