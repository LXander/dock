import os,sys,re
from glob import glob
import pandas as pd
sys.path.append(os.path.dirname(sys.path[0]))
from util.createfolder import try_create_chain_parent_folder


def convert():
    source = '/home/ubuntu/xiao/data/overlap/crystal_ligands'
    dest = '/home/ubuntu/xiao/data/overlap/crystal'
    for sourcePath in glob(os.path.join(source,'*','*.pdb')):
        destPath = re.sub(source,dest,sourcePath)
        cmd = 'obabel -ipdb {} -opdb {}'.format(sourcePath,destPath)
        try_create_chain_parent_folder(destPath)
        os.system(cmd)

def index():
    decoys_source = '/home/ubuntu/xiao/data/newkaggle/dude/train_set/decoys_av4'
    crystal_source = '/home/ubuntu/xiao/data/overlap_av4'

    records = []
    for decoys in glob(os.path.join(decoys_source,'*','*[_]*.av4')):
        receptor = decoys.split('/')[-2]
        crystal = os.path.join(crystal_source,receptor,receptor+'_crystal.av4')
        records.append(decoys,crystal)
        assert os.path.exists(decoys) and os.path.exists(crystal)

    df = pd.DataFrame(records,columns=['decoys','crystal'])


    df.to_csv("decoy_index.csv",index=False)

index()