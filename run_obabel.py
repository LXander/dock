import numpy as np
import pandas as pd
import config
import os,sys
import shutil

def get_pdb_and_crystal(input_file):
    # source to place crystal ligand
    crystal_source = os.path.join(config.BASE_DATA, 'H', 'addH')
    # dest to place converted ligand
    crystal_dest = os.path.join(config.BASE_DATA, 'filter_rmsd', 'crystal')

    # source to place pdb
    pdb_source = os.path.join(config.BASE_DATA, 'H', 'data')
    # dest to place pdb
    pdb_dest = os.path.join(config.BASE_DATA, 'filter_rmsd', 'receptor')

    #dirname = os.path.dirname(input_file)
    #receptor = os.path.basename(dirname)
    receptor = input_file.split('_')[0]
    filename = input_file+'_ligand'

    receptor_in = os.path.join(pdb_source,receptor,receptor+'.pdb')
    shutil.copy(receptor_in,pdb_dest)

    crystal_in = os.path.join(crystal_source,receptor,filename+'.pdb')
    crystal_path = os.path.join(crystal_dest,receptor)
    if not os.path.exists(crystal_path):
        os.mkdir(crystal_path)

    crystal_out = os.path.join(crystal_path,filename+'.pdb')

    cmd = 'obabel -ipdb %s -opdb -O %s -d'%(crystal_in,crystal_out)
    os.system(cmd)

def convert(item):
    PDB =item['PDBname']
    RS = item['PDBResId']
    RES,Id = RS.split('_')

    source_base = '/n/scratch2/xl198/data/result'
    source_file_path= os.path.join(source_base,PDB,'_'.join(PDB,RES,'ligand','fast.mol'))

    dest_base = '/n/scratch2/xl198/data/filter_rmsd/ligands'
    dest_path = os.path.join(dest_base,PDB)
    if not os.path.exists(dest_path):
        os.mkdir(dest_path)
    dest_file_path = os.path.join(dest_path,'_'.join(PDB,RES,'ligand','fast',Id+'.pdb'))

    cmd = 'obabel -imol2 %s -f %s -l %s -opdb -O %s '%(source_file_path,Id,Id,dest_file_path)

    os.system(cmd)

def run(base, offset):
    df = pd.read_csv('/n/scratch2/xl198/data/remark/filter.csv')
    convert(df.ix[base*1000+offset-1])


def test():
    base=0
    offset = 1
    df = pd.read_csv('/n/scratch2/xl198/data/remark/filter.csv')
    convert(df.ix[base * 1000 + offset - 1])

def get():
    df = pd.read_csv('/n/scratch2/xl198/data/remark/avail.csv')
    for d in df:
        get_pdb_and_crystal(d['ID'])

def main():
    args = sys.argv
    if len(args) >= 3:
        base = int(args[1])
        offset = int(args[2])
        print 'base %d offset %d' % (base, offset)
        run(base, offset)

if __name__ == '__main__':
    get()