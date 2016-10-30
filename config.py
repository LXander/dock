import os,sys

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
