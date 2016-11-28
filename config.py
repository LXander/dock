import os,sys

'''
Ususally I put all the path in this file
but sometime is't not so convenitent...
'''

#smina's path
SMINA = '/home/xl198/dock/Orchestra/dataJob/smina.static'

#base path for the data on Orchestra
_BASE = '/n/scratch2/xl198/data/H'

# commands have been generated and store here
_SCRIPT = '/home/xl198/dock/Orchestra/dataJob'

# every line store a command
COMMANDS_FILE = [
    os.path.join(_SCRIPT, 'filter_new_wp_fast.txt'),
    os.path.join(_SCRIPT, 'filter_new_wp_rigorous.txt'),
    os.path.join(_SCRIPT, 'filter_new_so_rigorous.txt')
]

# FILTER will fill the path
ROW_FILE = [
    'new_wp_fast.txt',
    'new_wp_rigorous.txt',
    'new_so_rigorous.txt'
]

JOBARRAY = '/home/xl198/code/bucket'

FILTER = os.path.join(_SCRIPT,'run.py')
RUN = os.path.join(os.getcwd(),'run.py')
CONVERT = os.path.join(os.getcwd(),'get_split_ligands.py')
INSERT = os.path.join(os.getcwd(),'new_line_inserter.py')
_BASE_SCRIPT = os.getcwd()

# _BASE_SELECT and ROW_FOLDER are used in select.py
_BASE_SELECT = '/n/scratch2/xl198/data/select'
COMMANDS_SELECT = {
    'temp' :  os.path.join(_BASE_SELECT,'temp'),
    'final' : os.path.join(_BASE_SELECT,'final')
}

ROW_FOLDER = {
    'wp_fast': os.path.join(_BASE,'wp_fast'),
    'wp_rigorous': os.path.join(_BASE,'wp_rigorous'),
    'so_rigorous' : os.path.join(_BASE,'so_rigorous')
}

# BASE_CONVERT is used in new_line_inserter.py
BASE_YI = '/n/scratch2/yw174/result/fast'

BASE_CONVERT = '/n/scratch2/xl198/data/result'
BASE_CONVERT2PDB = '/n/scratch2/xl198/data/pdbs'

# BASE_DATA is the root directory of almost all the operation
BASE_DATA = '/n/scratch2/xl198/data'

# REPORT is the folder to save the result
REPORT = '/home/xl198/report'

AVAIL = '/n/scratch2/xl198/data/remark/avail.csv'
