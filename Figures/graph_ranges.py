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


def initialize_dict(parse):
    """=============================================================================================
    This function accepts a type of information to be parsed and returns a dictionary with
    corresponding keys.

    :param parse: information of interest i.e. pairwise identity or alignment length
    :return: two dicts
    ============================================================================================="""

    if parse == 'id':
        dct = {9: [0, 0], 19: [0, 0], 29: [0, 0], 39: [0, 0], 49: [0, 0],
             59: [0, 0], 69: [0, 0], 79: [0, 0], 89: [0, 0], 99: [0, 0]}
    elif parse == 'len':
        dct = {499: [0, 0], 999: [0, 0], 1499: [0, 0], 1999: [0, 0], 2499: [0, 0]}

    return dct


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
    dict_m1 = initialize_dict(parse)
    dict_m2 = initialize_dict(parse)

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
    dict_m1 = initialize_dict(parse)
    dict_m2 = initialize_dict(parse)

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


def avg_dict(dct):
    """=============================================================================================
    This function takes a dictionary with a list of values and returns a dictionary with the average
    of the values.

    :return: dict
    ============================================================================================="""

    # Dicts are global
    for key, value in dct.items():
        if value[1] > 10: # Want at least 10 alignments in this range before we average
            dct[key] = value[0]/value[1]*100
        else:
            dct[key] = 0


def build_graph(dict_m1, dict_m2, dict_m3, ref, parse):
    """=============================================================================================
    This function takes a set of dictionaries and graphs them.

    :param dict_m1: dictionary of method 1
    :param dict_m2: dictionary of method 2
    :param dict_m3: dictionary of method 3
    :param ref: reference alignments
    :param parse: information of interest
    ============================================================================================="""

    # Vars for graph title
    if ref == '1':
        ref = 'RV11/12'
    elif ref == '9':
        ref = 'RV911/912/913'
    if parse == 'id':
        parse = 'Pairwise Identity'
    elif parse == 'len':
        parse = 'Alignment Length'

    # Ignore keys with 0 values
    x = [key for key in dict_m1.keys() if dict_m1[key] != 0]
    y1 = [val for val in dict_m1.values() if val != 0]
    y2 = [val for val in dict_m2.values() if val != 0]

    # When making alignment length plot for references 911/912/913, have to include all values
    # because the TCS values are 0 for lengths in last bucket
    y3 = [val for val in dict_m3.values()]

    # Plot
    plt.plot(x, y1, label='PEbA')
    plt.plot(x, y2, label='BLOSUM62')
    plt.plot(x, y3, label='DEDAL')
    plt.xlabel(f'{parse}')
    plt.ylabel('Average TCS')
    plt.title(f'Average TCS vs {parse} in {ref}')
    plt.legend()
    plt.grid()
    plt.show()


def main():
    """=============================================================================================
    ============================================================================================="""

    parser = argparse.ArgumentParser()
    parser.add_argument('-p1', type=str, default='/home/ben/Desktop/PEbA_Data/Runs/PEBA-BLOSUM')
    parser.add_argument('-p2', type=str, default='/home/ben/Desktop/PEbA_Data/Runs/PEBA-DEDAL')
    parser.add_argument('-r', type=str, default='9')
    parser.add_argument('-t', type=str, default='len')
    args = parser.parse_args()

    # Initialize dict based on type of info to be parsed
    dict_m1 = initialize_dict(args.t)
    dict_m2 = initialize_dict(args.t)
    dict_m3 = initialize_dict(args.t)

    # Run directories
    runs = sorted(os.listdir(args.p1))
    if args.r == '1':  # '1' stands for references 11 and 12
        runs = runs[:2]
    elif args.r == '9':  # '9' stands for references 911, 912, and 913
        runs = runs[2:]

    # Parse csv's in each run
    for run in runs:
        dict2_m1, dict2_m2 = parse_run(f'{args.p1}/{run}', args.t)

        # Add to running total
        for key, value in dict_m1.items():
            value[0] += dict2_m1[key][0]
            value[1] += dict2_m1[key][1]
        for key, value in dict_m2.items():
            value[0] += dict2_m2[key][0]
            value[1] += dict2_m2[key][1]

    # Same as above but for DEDAL
    runs = sorted(os.listdir(args.p2))
    if args.r == '1':
        runs = runs[:2]
    elif args.r == '9':
        runs = runs[2:]

    # Parse csv's in each run
    for run in runs:
        dict2_m1, dict2_m2 = parse_run(f'{args.p2}/{run}', args.t)
        for key, value in dict_m3.items():
            value[0] += dict2_m2[key][0]
            value[1] += dict2_m2[key][1]

    # Average values and build graph
    for dct in [dict_m1, dict_m2, dict_m3]:
        avg_dict(dct)
        print(dct)
    build_graph(dict_m1, dict_m2, dict_m3, args.r, args.t)


if __name__ == '__main__':
    main()
