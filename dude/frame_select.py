import os,sys
import numpy as np
from glob import glob
sys.path.append(os.path.dirname(sys.path[0]))
from util.createfolder import try_create_chain_parent_folder

def select():
    ligand_list = glob(os.path.join(FLAGS.sourcePath,'*','*.pdb'))
    for ligand in ligand_list:
        dest_ligand = ligand.replace(FLAGS.sourcePath,FLAGS.destPath)
        dest_ligand = dest_ligand.replace('.','_1.')
        try_create_chain_parent_folder(dest_ligand)
        cmd = 'obabel -ipdb {} -f 1 -l 1 -opdb -O {}'.format(ligand,dest_ligand)
        print cmd




class FLAGS:
    sourcePath = '/n/scratch2/xl198/dude/data/dude/pdbs/ligands'
    destPath = '/n/scratch2/xl198/dude/data/dude_1'


if __name__ == '__main__':
    select()