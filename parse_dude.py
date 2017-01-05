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

def gunzip(gzfile):
    if os.path.exists(gzfile):
        os.popen("gunzip -d {}".format(gzfile))
    tmp = re.search('(.*)\.gz',gzfile)
    if not tmp:
        raise Exception("can't parse mol2 file from gz file")
    mol2file = tmp.groups()[0]
    if not os.path.exists(mol2file):
        mol2file = None
        #raise Exception("can't fild active/decoy mol2 file from",mol2file)
    return mol2file

def convert(receptor):

    # for every folder in dude database
    #

    # Receptor
    receptorSource = os.path.join(sourceBase,receptor,'receptor.pdb')
    createFolder(os.path.join(destBase,'original_receptors'))
    receptorDest = os.path.join(destBase,'original_receptors',receptor+'.pdb')
    #receptorCmd = "obabel -ipdb {} -opdb -O {} ".format(receptorSource,receptorDest)
    # snima can't parse receptor converted by obabel
    # prody can't parse original receptor from dude
    receptorCmd = "cp {} {}".format(receptorSource,receptorDest)
    os.popen(receptorCmd)

    # Crystal
    crystalSource = os.path.join(sourceBase,receptor,'crystal_ligand.mol2')
    crystalDestFolder = os.path.join(destBase,'crystal_ligands',receptor)
    createFolder(crystalDestFolder)

    crystalDest = os.path.join(crystalDestFolder,receptor+'_crystal.pdb')
    crystalCmd = 'obabel -imol2 {} -opdb -O {}'.format(crystalSource,crystalDest)
    os.popen(crystalCmd)


    # Actives
    activesGz = os.path.join(sourceBase,receptor,'actives_final.mol2.gz')
    activesMol2 = gunzip(activesGz)
    if not activesMol2:
        return None
    activesDestFolder = os.path.join(destBase,'docked_ligands',receptor)
    createFolder(activesDestFolder)
    activesDestFile = os.path.join(activesDestFolder,receptor+'_actives_.pdb')
    activesCmd = 'obabel -imol2 {} -opdb -O {} -m'.format(activesMol2,activesDestFile)
    os.popen(activesCmd)

    # Decoys
    decoysGz = os.path.join(sourceBase,receptor,'decoys_final.mol2.gz')
    decoysMol2 = gunzip(decoysGz)
    if not decoysMol2:
        return None
    decoysDestFolder = os.path.join(destBase,'docked_ligands',receptor)
    createFolder(decoysDestFolder)
    decoysDestFile = os.path.join(decoysDestFolder,receptor+'_decoys_.pdb')
    decoysCmd = 'obabel -imol2 {} -opdb -O {} -m'.format(decoysMol2,decoysDestFile)
    os.popen(decoysCmd)

def runs():
    map(lambda receptor:convert(receptor),os.listdir(sourceBase))



#run()
#convert('comt')

if __name__=='__main__':
    createBaseFolder()
    receptors = os.listdir(sourceBase)
    args = sys.argv
    if len(args)<2:
        print "not enouth args, need 2 input {}".format(len(args))
    else:
       offset = int(args[1])-1
       if offset<len(receptors):
            convert(receptors[offset])
       else:
            print "offset overflow, total size {} input {}".format(len(receptors),offset)


