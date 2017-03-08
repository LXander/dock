import numpy as np
import pandas as pd
import config
import os,sys
import shutil

'''
In our current database docking result stored in one file
This code used obabel to convert specific ligand from file
and store them orderly in given path
The souce is a csv file contain columns ['PDBname','PDBResId']
'''

dest_pdb_folder = 'mar_06'
dest_csv_index = '/n/scratch2/xl198/YI/rigor_so/small_rmsd_rigor_so.csv'

def get_pdb_and_crystal(input_file):
    # source to place crystal ligand
    crystal_source = os.path.join(config.BASE_DATA, 'H', 'addH')
    # dest to place converted ligand
    crystal_dest = os.path.join(config.BASE_DATA, dest_pdb_folder, 'crystal_ligands')
    if not os.path.exists(crystal_dest):
        os.mkdir(crystal_dest)

    # source to place pdb
    pdb_source = os.path.join(config.BASE_DATA, 'H', 'data')
    # dest to place pdb
    pdb_dest = os.path.join(config.BASE_DATA, dest_pdb_folder, 'receptors')
    if not os.path.exists(pdb_dest):
        os.mkdir(pdb_dest)

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
    if not os.path.exists(crystal_out):
        cmd = 'obabel -ipdb %s -opdb -O %s'%(crystal_in,crystal_out)
        os.system(cmd)

def convert(item):
    PDB =item['PDBname']
    RS = item['PDBResId']
    RES,Id = RS.split('_')

    source_base = '/n/scratch2/xl198/YI/rigor_so/final'
    source_file_path= os.path.join(source_base,PDB,'_'.join([PDB,RES,'ligand','rigor_so.pdb']))

    dest_base = os.path.join(config.BASE_DATA, dest_pdb_folder, 'docked_ligands')

    if not os.path.exists(dest_base):
        os.mkdir(dest_base)
    dest_path = os.path.join(dest_base,PDB)
    if not os.path.exists(dest_path):
        os.mkdir(dest_path)
    dest_file_path = os.path.join(dest_path,'_'.join([PDB,RES,'ligand','rigor_so',Id+'.pdb']))

    cmd = 'obabel -ipdb %s -f %s -l %s -opdb -O %s '%(source_file_path,Id,Id,dest_file_path)

    os.system(cmd)

def run(base, offset):
    df = pd.read_csv(dest_csv_index)
    print 'num',len(df)
    total = int(len(df)*1.0/1000)+1
    print 'total',total
    for i in range(total):
    	convert(df.ix[base*1000+(offset-1)*total+i])
        get_pdb_and_crystal(df.ix[base*1000+(offset-1)*total+i]['ID'])

def main():
    args = sys.argv
    if len(args) >= 3:
        base = int(args[1])
        offset = int(args[2])
        print 'base %d offset %d' % (base, offset)
        run(base, offset)



if __name__ == '__main__':
    main()
