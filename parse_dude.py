import os,sys
import numpy as np
import pandas as pd
import re

# convert dude's data in the way we use


sourceBase = '/n/scratch2/xl198/dude/data/all'
destBase = '/n/scratch2/xl198/dude/data/dude/pdbs'

def createFolder(folderPath):
    if not os.path.exists(folderPath):
        os.mkdir(folderPath)

def createBaseFolder():
    BaseFolders =['crystal_ligands' , 'docked_ligands' , 'receptors']
    map(lambda folder:createFolder(os.path.join(destBase,folder)),BaseFolders)

def convert(receptor):

    # for every folder in dude database
    receptorSource = os.path.join(sourceBase,receptor,'receptor.pdb')
    receptorDest = os.path.join(destBase,'receptors',receptor+'.pdb')
    receptorCmd = "obabel -ipdb {} -opdb -O {} ".format(receptorSource,receptorDest)
    os.popen(receptorCmd)

    crystalSource = os.path.join(sourceBase,receptor,'crystal_ligand.mol2')
    crystalDestFolder = os.path.join(destBase,'crystal_ligands',receptor)
    createFolder(crystalDestFolder)

    crystalDest = os.path.join(crystalDestFolder,receptor+'_crystal.pdb')
    crystalCmd = 'obabel -imol2 {} -opdb -O {}'.format(crystalSource,crystalDest)
    os.popen(crystalCmd)

    #activesSource = os.path.join(sourceBase,receptor,'')

def run():
    map(lambda receptor:convert(receptor),os.listdir(sourceBase))

createBaseFolder()
run()