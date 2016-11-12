import config
import re
import os,sys
import numpy as np
import pandas as pd

def count_liangd_num(input_file):
    '''
    count the number of ligands for given ligands in mol2 foramt.

    :param input_file: path of input file eg. '/n/scratch2/xl198/data/source/ligands/1a2b/1a2b_1234_liagnd.mol2'
    :return:  [filename,atom_num] eg [1a2b_1234_liagnd,10]
    '''

    # count the number of the single ligand in one file
    ligand_count = 0
    with open(input_file) as input:
        for line in input:
            if line == '@<TRIPOS>MOLECULE\n':
                ligand_count += 1

    return [os.path.basename(input_file),ligand_count]


def count_atom_num(input_file):
    '''
    count the number of atoms for given ligand in mol2 foramt.

    :param input_file: path of input file eg. '/n/scratch2/xl198/data/source/ligands/1a2b/1a2b_1234_liagnd.mol2'
    :return: [filename,atom_num] eg [1a2b_1234_liagnd,10]
    '''



    # count the number of the atoms in ligand
    flag = False
    atom_num = 0
    with open(input_file) as input:
        for line in input:
            if not flag and line == '@<TRIPOS>ATOM\n':
                flag = True
            elif line == '@<TRIPOS>BOND\n':
                break
            elif flag:
                atom_num += 1


    return [os.path.basename(input_file),atom_num]



def read_file_path(input_path):
    '''
    Get the path of file used as input of count function
    :param input_path:
    :return:
    '''

    for dirname, dirnames, filenames in os.walk(input_path):
        for filename in filenames:
            file_path = os.path.join(dirname, filename)
            yield  file_path



def count_and_report(input_path,report,counter):

    pathes = read_file_path(input_path)
    result = [ counter(path) for path in pathes ]
    df = pd.DataFrame(data = result, columns=['file','atom_num'])
    df.to_csv(report,index=False)

def main():
    count_and_report(config.BASE_YI,
                     os.path.join(config.REPORT,'liangd_count.csv'),
                     count_liangd_num)

if __name__ == '__main__':
    main()
