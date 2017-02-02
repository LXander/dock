import os,sys
import re

import numpy as np
import pandas as pd
sys.path.append(os.path.dirname(sys.path[0]))
from util.createfolder import try_create_chain_parent_folder

dude_index_csv = ''

def get_crystal():
    df = pd.read_csv(dude_index_csv)
    crystal_ligand = list(set(df['ligand_box']))

    crystal = pd.DataFrame(crystal_ligand,columns=['sourcePath'])
    crystal['destPath'] = crystal['sourcePath'].apply(re.sub('/n/scratch2/xl198/dude/data/dude/pdbs/','/n/scratch2/xl198/dude/data/dude/pdbs/'))
