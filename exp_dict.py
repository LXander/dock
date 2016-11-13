import numpy as np
import pandas as pd
import re
import os,sys



def read_single_path(input_path):
    for dirname, dirnames, filenames in os.walk(input_path):
        for filename in filenames:
            file_path = os.path.join(dirname, filename)
            return file_path


def read_file_path(input_path):
    '''
    Get the path of file used as input of count function
    :param input_path: path of data
    :return: generator every time return one path of file
    '''

    for dirname, dirnames, filenames in os.walk(input_path):
        for filename in filenames:
            file_path = os.path.join(dirname, filename)
            yield  file_path

def get_remark_columns(file_path):
    '''
    extract the remark from experimental


    :param file_path: file path
    :return: columns, list
    '''
    name = os.path.basename(file_path)
    brand = '_'.join(name.split('_')[0:2])

    with open(file_path) as fr:
        remark = fr.readline()

    if not re.search('^# Remark', remark):
        columns = None

    else:
        remark = remark.split('*@*')
        remark = [r.split(':') for r in remark[1:]]
        columns = [r[0] for r in remark]
        columns.insert(0, 'ID')


    return columns

def get_remark_data(file_path):
    '''
    extract the remark from experimental

    :param file_path: file path
    :return:
    '''
    name = os.path.basename(file_path)
    brand = '_'.join(name.split('_')[0:2])

    with open(file_path) as fr:
        remark = fr.readline()

    if not re.search('^# Remark.',remark):

        data = None
    else:
        remark = remark.split('*@*')
        remark = [ r.split(':') for r in remark[1:] ]
        columns = [ r[0] for r in remark ]
        columns.insert(0,'ID')
        data = [ r[1] for r in remark ]
        data.insert(0,brand)

    return data


def get_remarks():
    '''
    extract remark from exp data

    :return:
    '''
    input_path = '/n/scratch2/yw174/result/experiment'
    output_file_path = '/n/scratch2/xl198/data/remark/exp.csv'
    columns = get_remark_columns(read_single_path(input_path))

    walk = read_file_path(input_path)
    data = [ get_remark_data(file_path) for file_path in walk ]

    df = pd.DataFrame( data = data, columns = columns )
    df.to_csv(output_file_path,index=False)

def main():
    get_remarks()


if __name__ == '__main__':
    main()
