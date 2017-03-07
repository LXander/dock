import os,sys
import numpy as np
import pandas as pd
import re
from glob import glob


class FLAGS:
    sourceBase = '/n/scratch2/xl198/dude/data/all'
    destBase = '/n/scratch2/xl198/dude/data/dude/pdbs'

def get_repair_list():
    '''
    Some receptor in dude doesn't has mol2 for its ligands
    find out those file
    :return:
    '''

    actives_files = glob(os.path.join(FLAGS.sourceBase,'*','actives_final.mol2.gz'))
    actives_folder = map(lambda x:os.path.dirname(x),actives_files)

    decoys_files = glob(os.patt.join(FLAGS.sourceBase,'*','decoys_final.mol2.gz'))

