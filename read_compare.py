"""================================================================================================
This script reads csv files containing results from compute_pra and parses for relevant info.

Ben Iovino  04/05/23   VecAligns
================================================================================================"""

import csv
import os


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


def parse_data(data):
    """=============================================================================================
    This function accepts a list of lists and parses for relevant info.

    :param data: list of lists
    :return: list of lists
    ============================================================================================="""

    count = 0
    for line in data:
        if count == 0 or count % 2 == 0:
            for key, value in COMPARE_DICT_M1.items():
                if float(line[3]) <= key:
                    value[0] += float(line[0])/100
                    value[1] += 1
                    break
        else:
            for key, value in COMPARE_DICT_M2.items():
                if float(line[3]) <= key:
                    value[0] += float(line[0])/100
                    value[1] += 1
                    break
        count += 1


def avg_dict():
    """=============================================================================================
    This function takes a dictionary with a list of values and returns a dictionary with the average
    of the values.

    :return: dict
    ============================================================================================="""

    # Dicts are global
    for key, value in COMPARE_DICT_M1.items():
        COMPARE_DICT_M1[key] = value[0]/value[1]*100
    for key, value in COMPARE_DICT_M2.items():
        COMPARE_DICT_M2[key] = value[0]/value[1]*100

    print(f'M1: {COMPARE_DICT_M1}')
    print()
    print(f'M2: {COMPARE_DICT_M2}')


def main():
    """=============================================================================================
    This function initializes global dictionaries where values from csv files are stored in buckets.
    It then reads every csv file in the directory structure and parses for info of interest. It then
    finds the average of each bucket and prints the results.
    ============================================================================================="""

    # Store values from csv
    global COMPARE_DICT_M1  #pylint: disable=W0601
    global COMPARE_DICT_M2  #pylint: disable=W0601

    # Identity buckets
    COMPARE_DICT_M1 = {9: [0, 0], 19: [0, 0], 29: [0, 0], 39: [0, 0], 49: [0, 0],
                     59: [0, 0], 69: [0, 0], 79: [0, 0], 89: [0, 0], 99: [0, 0]}
    COMPARE_DICT_M2 = {9: [0, 0], 19: [0, 0], 29: [0, 0], 39: [0, 0], 49: [0, 0],
                     59: [0, 0], 69: [0, 0], 79: [0, 0], 89: [0, 0], 99: [0, 0]}
    
    # Length buckets
    # COMPARE_DICT_M1 = {499: [0, 0], 999: [0, 0], 1499: [0, 0], 1999: [0, 0], 2499: [0, 0]}
    # COMPARE_DICT_M2 = {499: [0, 0], 999: [0, 0], 1499: [0, 0], 1999: [0, 0], 2499: [0, 0]}

    # Directory structure -> set/run/ref/msa/compare.csv
    # Want to read every single csv
    path = '/home/ben/Desktop/PEbA_Data/Runs/gen3/PROTT5-ESM2'
    for run in os.listdir(path):  #pylint: disable=R1702
        for ref in os.listdir(f'{path}/{run}'):
            for msa in os.listdir(f'{path}/{run}/{ref}'):
                if msa.startswith('B'):
                    for file in os.listdir(f'{path}/{run}/{ref}/{msa}'):
                        if file.endswith('csv'):

                            # Read csv and parse
                            data = read_csv(f'{path}/{run}/{ref}/{msa}/{file}')
                            parse_data(data)

    # Find average for each key
    avg_dict()


if __name__ == "__main__":
    main()
