import os,sys
import numpy as np
from glob import glob
sys.path.append(os.path.dirname(sys.path[0]))
from util.createfolder import try_create_chain_parent_folder

def select(index=None):
    count = 0
    ligand_list = glob(os.path.join(FLAGS.sourcePath,'*','*.pdb'))
    if index:
        segment = len(ligand_list)/1000
        indexes = range((index-1)*segment,index*segment)
    else:
        indexes = range(len(ligand_list))
    for ind in indexes:
        ligand = ligand_list[ind]
        count +=1
        print count
        dest_ligand = ligand.replace(FLAGS.sourcePath,FLAGS.destPath)
        dest_ligand = dest_ligand.replace('.','_1.')
        try_create_chain_parent_folder(dest_ligand)
        cmd = 'obabel -ipdb {} -f 1 -l 1 -opdb -O {}'.format(ligand,dest_ligand)
        if not os.path.exists(dest_ligand):
            os.system(cmd)




class FLAGS:
    sourcePath = '/n/scratch2/xl198/dude/data/dude/docked_dude/docked_ligands'
    destPath = '/n/scratch2/xl198/dude/data/dude_1'


if __name__ == '__main__':
    args = sys.argv
    if len(args)>1:
        select(int(args[1]))
    else:
        select()
