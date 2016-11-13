import numpy as np
import pandas as pd
import re
import os,sys


def get_id(filename):
    '''
    extract the ID of a file
    :param filename: eg '1a2b_1234.pdb'
    :return: eg '1a2b_1234'
    '''
    ID = '_'.join(filename.split('_')[0:2])
    return ID

def standarlize_exp( data):
    '''
    only first two columns is strng, and the rest are float,
    and some of the rest would have multiple value, should split them
    it's hard to convert with operator like gt or lt so we discard them at this time
    :param data: list of the remark
    :return: standarlized data
    '''

    result = []
    for i in range(len(data)):
        if i < 2:
            # they're string so keep it as original
            result.append(data[i])
        elif data[i] == '':
            # Nothing here
            result.append(None)
        else:

            m = re.match('\[(.*)\]', data[i])
            if m:
                # it was a list
                nums = m.group(1)
                nums = nums.split(',')
                nums = [float(n) for n in nums]
                result.append(nums)
            else:
                nums = data[i].split('|')
                nums = [float(num) for num in nums if not re.search('>|<', num)]
                result.append(nums)
    return result

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
        data = standarlize_exp(data)
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
