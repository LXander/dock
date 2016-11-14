import os,sys

def path_check(path):
    # to check if a path is exists
    if not os.path.exists(path):
        print "%s not exists!"%path
        exit(1)
    else:
        return path

def debug(debug_flag,output,content):
    # print debug information into stdout
    if debug_flag:
        print output%content


def get_id(self, filename):
    '''
    extract the ID of a file
    :param filename: eg '1a2b_1234.pdb'
    :return: eg '1a2b_1234'
    '''
    ID = '_'.join(filename.split('_')[0:2])
    return ID

