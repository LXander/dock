import os,sys,re
import prody
import numpy as np


def create_chain_parent_folder(filePath):
    dirPath = os.path.dirname(filePath)
    create_chain_folder(dirPath)

def create_chain_folder(folderPath):
    dirPath = os.path.dirname(folderPath)
    if not os.path.exists(dirPath):
        create_chain_folder(dirPath)
    if not os.path.exists(folderPath):
        os.mkdir(folderPath)

def create_folder(folderPath):
    # If folder doesn't exists create it, if exists ignore
    if not os.path.exists(folderPath):
        os.mkdir(folderPath)

def create_parent_folder(filePath):
    dirPath = os.path.dirname(filePath)
    create_folder(dirPath)


def path_iterator(chain):
    for dirpath,dirname,filenames in chain:
        for filename in filenames:
            yield os.path.join(dirpath,filename)

def docked_ligand_overlaps_with_crystal(docked_lig_path,cryst_lig_path,tanimoto_cutoff=0.75,clash_cutoff_A=4,clash_size_cutoff=0.3):
    """This script takes docked ligand, and searches crystal ligands to find if doked position has a chance to be correct
    (overlaps). And in that case returns one."""

    def tanimoto_similarity(cryst_lig_path, docked_lig_path):
        command = os.popen('babel -d {} {} -ofpt -xfFP4'.format(cryst_lig_path, docked_lig_path))
        ls = command.read()
        tanimoto_similarity = re.split('=|\n', ls)[2]
        return tanimoto_similarity

    def calculate_clash(cryst_lig_path, docked_lig_path,clash_cutoff_A=clash_cutoff_A):
        """calculates the ratio of atoms between two molecules that appear within van Der Waals cutoff of 4A"""
        # calculate atom to atom distances
        prody_docked_ligand = prody.parsePDB(docked_lig_path).getCoords()
        prody_cryst_ligand = prody.parsePDB(cryst_lig_path).getCoords()

        def atom_clashing(index):
            # calculate distance to all atoms in the other molecule and check if any of them falls below cutoff
            if any(np.sum((prody_cryst_ligand[index,:] - prody_docked_ligand)**2,axis=1)**0.5 < clash_cutoff_A):
                return True
            else:
                return False

        clash_mask = map(atom_clashing,np.arange(len(prody_cryst_ligand[:, 0])))
        return np.average(clash_mask)

    clash_found = False

    if tanimoto_similarity(cryst_lig_path, docked_lig_path) > tanimoto_cutoff:
        if calculate_clash(cryst_lig_path,docked_lig_path)>clash_size_cutoff:
            clash_found = True

    '''
    clash_found = False
    # walk through /crystal/ligands and see if any of them are similar
    # there should be at least one
    for dirpath, dirnames, filenames in os.walk(database_path + "/crystal_ligands/" + docked_lig_path.split("/")[-2]):
        for filename in filenames:
            if re.search('.pdb$', filename):
                cryst_lig_path = str(os.path.abspath(dirpath) + "/" + filename)
                if tanimoto_similarity(cryst_lig_path,docked_lig_path) > tanimoto_cutoff:
                    # if ligand in crystal is similar to the docked ligand, calculate clash
                    if (calculate_clash(cryst_lig_path,docked_lig_path) > clash_size_cutoff):
                        clash_found = True
    '''

    return clash_found