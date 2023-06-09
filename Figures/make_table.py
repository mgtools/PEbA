"""================================================================================================
This script reads csv files and builds tables with the results. The directories that are read
in this script are a result from compare_aligns.py. A set of runs comes from using two methods,
i.e. PEbA and BLOSUM, on each BAliBASE reference used in this project, of which there are five.
The resulting tables have a column for reach run and a row for each range of values, for either
pairwise identity or alignment length.

Ben Iovino  04/05/23   VecAligns
================================================================================================"""

import csv
import os
import pickle
import argparse
import pandas as pd


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

    avg_align = 0
    count = 0
    zeros = 0
    zeros2 = 0
    zero_avg = 0
    in_zero = False
    for line in data:

        # Even lines are from method 1
        if count == 0 or count % 2 == 0:
            for key, value in COMPARE_DICT_M1.items():
                if parse == 'id':
                    if float(line[3]) <= key:
                        if float(line[0]) == 0:
                            in_zero = True
                            zeros += 1
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
            for key, value in COMPARE_DICT_M2.items():
                if parse == 'id':
                    if float(line[3]) <= key:
                        if in_zero is True:
                            if float(line[0])>0:
                                zeros2 += 1
                                zero_avg += float(line[0])/100
                        in_zero = False
                        value[0] += float(line[0])/100
                        value[1] += 1
                        break
                elif parse == 'len':
                    if float(line[1]) <= key:
                        value[0] += float(line[0])/100
                        value[1] += 1
                        break

        avg_align += float(line[1])
        count += 1

    # Get average BLOSUM score for TC scores where PEbA is 0 and BLOSUM is > 0
    if zeros2 > 0:
        print(zero_avg/zeros2)

    print(zeros, zeros2)
    return avg_align, count


def avg_dict():
    """=============================================================================================
    This function takes a dictionary with a list of values and returns a dictionary with the average
    of the values.

    :return: dict
    ============================================================================================="""

    # Dicts are global
    for key, value in COMPARE_DICT_M1.items():
        if value[1] > 10: # Want at least 10 alignments in this range before we average
            COMPARE_DICT_M1[key] = value[0]/value[1]*100
        else:
            COMPARE_DICT_M1[key] = 0
    for key, value in COMPARE_DICT_M2.items():
        if value[1] > 10:
            COMPARE_DICT_M2[key] = value[0]/value[1]*100
        else:
            COMPARE_DICT_M2[key] = 0

    # Count number of pkl files
    count = 0
    for file in os.listdir('Figures'):
        if file.endswith('.pkl'):
            count += 1

    # Save to pickle
    with open(f'Figures/compare_dict_m1_{count}.pkl', 'wb') as file:
        pickle.dump(COMPARE_DICT_M1, file)
    with open(f'Figures/compare_dict_m2_{count}.pkl', 'wb') as file:
        pickle.dump(COMPARE_DICT_M2, file)


def build_table():
    """=============================================================================================
    This function reads pickle files and builds a table with the results.

    :return: Pandas dataframe
    ============================================================================================="""

    # Read pickle files
    pfiles = []
    for file in os.listdir('Figures'):
        if file.endswith('.pkl'):
            pfiles.append(file)

    # Sort files so that they are added to dict sequentially
    pfiles.sort()

    # Build dict
    dicts = {}
    for file in pfiles:
        with open(f'Figures/{file}', 'rb') as f:
            dicts[file] = pickle.load(f)
        os.remove(f'Figures/{file}')

    # Build table
    table = pd.DataFrame(dicts)

    # Split table into two
    table1 = table.iloc[:, :len(table.columns)//2]
    table2 = table.iloc[:, len(table.columns)//2:]

    # Rename columns
    refs = ['RV11', 'RV12', 'RV911', 'RV912', 'RV913']
    table1.columns = refs
    table2.columns = refs

    print(table1)
    print()
    print(table2)

    # Save to csv
    #table1.to_csv('Figures/table1.csv')
    #table2.to_csv('Figures/table2.csv')


def main():
    """=============================================================================================
    This function initializes global dictionaries where values from csv files are stored in buckets.
    It then reads every csv file in the directory structure and parses for info of interest. It then
    finds the average of each bucket and builds a dataframe with the results.
    ============================================================================================="""

    parser = argparse.ArgumentParser()
    parser.add_argument('-p', type=str, default='Data/Alignments/PEBA-BLOSUM')
    parser.add_argument('-t', type=str, default='id')
    args = parser.parse_args()

    # Values for each run
    global COMPARE_DICT_M1  #pylint: disable=W0601
    global COMPARE_DICT_M2  #pylint: disable=W0601

    # Sort runs so dictionaries are built sequentially
    runs = []
    for run in os.listdir(args.p):
        runs.append(args.p+'/'+run)
    runs.sort()

    # Directory structure -> set/run/ref/msa/compare.csv
    # Want to read every single csv
    avg_align, count = 0, 0
    for run in runs:  #pylint: disable=R1702
        for ref in os.listdir(f'{run}'):
            print(ref)

            # Initialize dicts for each run
            if args.t == 'id':
                COMPARE_DICT_M1 = {9: [0, 0], 19: [0, 0], 29: [0, 0], 39: [0, 0], 49: [0, 0],
                     59: [0, 0], 69: [0, 0], 79: [0, 0], 89: [0, 0], 99: [0, 0]}
                COMPARE_DICT_M2 = {9: [0, 0], 19: [0, 0], 29: [0, 0], 39: [0, 0], 49: [0, 0],
                     59: [0, 0], 69: [0, 0], 79: [0, 0], 89: [0, 0], 99: [0, 0]}
            elif args.t == 'len':
                COMPARE_DICT_M1 = {499: [0, 0], 999: [0, 0], 1499: [0, 0], 1999: [0, 0], 2499: [0, 0]}
                COMPARE_DICT_M2 = {499: [0, 0], 999: [0, 0], 1499: [0, 0], 1999: [0, 0], 2499: [0, 0]}

            for msa in os.listdir(f'{run}/{ref}'):
                print(msa)
                if msa.startswith('B'):
                    for file in os.listdir(f'{run}/{ref}/{msa}'):
                        if file.endswith('csv'):

                            # Read csv and parse
                            data = read_csv(f'{run}/{ref}/{msa}/{file}')
                            a, b = parse_data(data, args.t)
                            avg_align += a
                            count += b

        # Find average for each key
        avg_dict()

    # Build table with averaged dictionaries
    build_table()


if __name__ == "__main__":
    main()
