import pandas as pd
import config
import os,sys
import shutil
import numpy as np
import os,re, random,prody
from itertools import chain
from av2_atomdict import atom_dictionary

'''
Convert pdb data into npy

'''

database_path = '/n/scratch2/xl198/data/jan_01'
destination_path = '/n/scratch2/xl198/data/jan_01_npy'

def get_pdb_and_crystal(input_file):
    # source to place crystal ligand
    crystal_source = os.path.join(config.BASE_DATA, 'H', 'addH')
    # dest to place converted ligand
    crystal_dest = os.path.join(config.BASE_DATA, 'dec_17_large_npy', 'crystal_ligands')
    if not os.path.exists(crystal_dest):
        os.mkdir(crystal_dest)

    # source to place pdb
    pdb_source = os.path.join(config.BASE_DATA, 'H', 'data')
    # dest to place pdb
    pdb_dest = os.path.join(config.BASE_DATA, 'dec_17_large_npy', 'receptors')
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
        cmd = 'obabel -ipdb %s -opdb -O %s -d'%(crystal_in,crystal_out)
        os.system(cmd)


def ligand_atom_to_number(atomname):
    atomic_tag_number = atom_dictionary.LIG[atomname.lower()]
    return atomic_tag_number


def receptor_atom_to_number(atomname):
    atomic_tag_number = atom_dictionary.REC[atomname.lower()]
    return atomic_tag_number

def convert_ligand(item):
    print "convert",item
    try:
        prody_ligand = prody.parsePDB(item)
    except Exception:
        pass

    try:
        atom_numbers = map(ligand_atom_to_number, prody_ligand.getElements())
        coordinates_and_atoms = np.hstack((prody_ligand.getCoords(), np.reshape(atom_numbers, (-1, 1))))
        destination = re.sub(".pdb$",'',item)
        destination = re.sub(database_path,destination_path,destination)
        if not os.path.exists(os.path.dirname(destination)):
            os.mkdir(os.path.dirname(destination))
        np.save(destination, coordinates_and_atoms)
    except Exception:
        pass

def convert_receptor(item):
    try:
        prody_ligand = prody.parsePDB(item)
    except Exception:
        pass

    try:
        atom_numbers = map(receptor_atom_to_number, prody_ligand.getElements())
        coordinates_and_atoms = np.hstack((prody_ligand.getCoords(), np.reshape(atom_numbers, (-1, 1))))
        destination = re.sub(".pdb$",'',item)
        destination = re.sub(database_path,destination_path,destination)
        if not os.path.exists(os.path.dirname(destination)):
            os.mkdir(os.path.dirname(destination))
        np.save(destination, coordinates_and_atoms)
    except Exception:
        pass

def create_folder(folder_path):
    if not os.path.exists(folder_path):
        os.mkdir(folder_path)

def run(base, offset):
    dest_ligands_folder = os.path.join(destination_path,'docked_ligands')
    dest_crystal_folder = os.path.join(destination_path,'crystal_ligands')
    dest_receptor_folder = os.path.join(destination_path,'receptors')

    source_ligands_folder = os.path.join(database_path,'docked_ligands')
    source_crystal_folder = os.path.join(database_path,'crystal_ligands')
    source_receptor_folder = os.path.join(database_path,'receptors')


    map(create_folder,[dest_crystal_folder,dest_ligands_folder,dest_receptor_folder])

    cur = base*1000+offset-1
    con = 0
     

    for dirpath,dirname,filenames in chain(os.walk(source_crystal_folder),os.walk(source_ligands_folder)):
        for f in filenames:
             
            convert_ligand(os.path.join(dirpath,f))
            con += 1
            sys.stderr.write("convert %d\n"%(con))

    for dirpath,dirname,filenames in chain(os.walk(source_receptor_folder)):
        for f in filenames:
            convert_receptor(os.path.join(dirpath,f))
            con += 1
            sys.stderr.write("convert %d\n"%(con))


def main():
    args = sys.argv
    if len(args) >= 3:
        base = int(args[1])
        offset = int(args[2])
        print 'base %d offset %d' % (base, offset)
        run(base, offset)

if __name__ == '__main__':
    main()
