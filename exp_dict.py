import numpy as np
import pandas as pd
import re
import os, sys

'''
This is to used formalize the remard in experience data.

'''
def get_id(filename):
    '''
    extract the ID of a file
    :param filename: eg '1a2b_1234.pdb'
    :return: eg '1a2b_1234'
    '''
    ID = '_'.join(filename.split('_')[0:2])
    return ID

def standarlize_exp(data):
    '''
    only first two columns is strng, and the rest are float,
    and some of the rest would have multiple value, should split them
    it's hard to convert with operator like gt or lt so we discard them at this time
    :param data: list of the remark
    :return: standarlized data
    '''

    result = []
    for i in range(len(data)):
        if i < 2 or i == 9:
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
                print i
                print data[i]
                nums = data[i].split('|')
                nums = [float(num) for num in nums if not re.search('>|<|NA', num)]
                result.append(nums)
    return result

def standarlize_dock(data):
    result = []
    for i in range(len(data)):
        if i == 0:
            result.append(data[i])
        elif i == 1:
            result.append(data[i])
            result.append(int(data[i].split('_')[1]))
	elif data[i] == 'NA':
	    result.append('')
        else:

            m = re.match('\[(.*)\]', data[i])
            if m:
                # it was a list
                nums = m.group(1)
                nums = nums.split(',')
                nums = [float(n) for n in nums]
                result.append(nums)
            else:
		
                result.append(float(data[i]))
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
            yield file_path

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

    if not re.search('^# Remark.', remark):

        data = None
    else:
        remark = remark.split('*@*')
        remark = [r.split(':') for r in remark[1:]]
        columns = [r[0] for r in remark]
        columns.insert(0, 'ID')
        data = [r[1] for r in remark]
        data = standarlize_exp(data)
        data.insert(0, brand)

    return data

def get_dock_remark_columns(file_path):
    name = os.path.basename(file_path)
    brand = '_'.join(name.split('_')[0:2])

    remark = None
    with open(file_path) as fr:
        for line in fr:
            if re.search('^# Remark', line):
                remark = line.strip('\n')



    if remark:
        remark = remark.split('_{')[1:]
        remark = [ r.strip('}') for r in remark]
        remark = [ r.split(':') for r in remark]
        columns = [ r[0] for r in remark]
        columns.insert(0, 'ID')
        columns.insert(3,'inner_index')

    return columns

def get_remarks_exp():
    '''
    extract remark from exp data

    :return:
    '''
    input_path = '/n/scratch2/yw174/result/experiment'
    output_file_path = '/n/scratch2/xl198/data/remark/exp.csv'
    columns = get_remark_columns(read_single_path(input_path))

    walk = read_file_path(input_path)
    data = [get_remark_data(file_path) for file_path in walk]

    df = pd.DataFrame(data=data, columns=columns)
    df.to_csv(output_file_path, index=False)

def get_dock_remark_data(file_path):
    '''
        extract the remark from experimental

        :param file_path: file path
        :return:
        '''
    name = os.path.basename(file_path)
    brand = '_'.join(name.split('_')[0:2])

    remarks = []
    with open(file_path) as fr:
        for line in fr:
            if re.search('^# Remark.', line):
                remarks.append(line.strip('\n'))

    data = []

    for remark in remarks:
        remark = remark.split('{')[1:]
        remark = [r.strip('}_') for r in remark]
        remark = [r.split(':') for r in remark]
        datum = [r[1] for r in remark]
        result = standarlize_dock(datum)
	result.insert(0,brand)
        data.append(result)

    return data

def get_remarks_dock():
    '''
    extract remark from dock result
    :return:
    '''
    input_path = '/n/scratch2/yw174/result/rigor_so'
    output_file_path = '/n/scratch2/xl198/data/remark/rigor_so.csv'
    columns = get_dock_remark_columns(read_single_path(input_path))

    walk = read_file_path(input_path)
    raw_data = [get_dock_remark_data(file_path) for file_path in walk]
    data = []
    for raw in raw_data:
        for r in raw:
            data.append(r)

    df = pd.DataFrame(data = data,columns = columns)
    df.to_csv(output_file_path,index=False)



def main():
    get_remarks_dock()


if __name__ == '__main__':
    main()
