import os,sys
import numpy as np
import pandas as pd
import re

# convert dude's data in the way we use


sourceBase = '/n/scratch2/xl198/dude/data/all'
destBase = '/n/scratch2/xl198/dude/data/dude/pdbs'
dockingIndex = '/n/scratch2/xl198/dude/code/dude.csv'
smina = '/home/xl198/program/smina/smina.static'

def createFolder(folderPath):
    if not os.path.exists(folderPath):
        os.mkdir(folderPath)

def createBaseFolder():
    BaseFolders =['crystal_ligands' , 'ligands' , 'receptors']
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
    createFolder(os.path.join(destBase,'addh_receptors'))
    receptorDest = os.path.join(destBase,'addh_receptors',receptor+'.pdb')
    receptorCmd = "obabel -ipdb {} -opdb -O {} -h ".format(receptorSource,receptorDest)
    # snima can't parse receptor converted by obabel
    # prody can't parse original receptor from dude
    #receptorCmd = "cp {} {}".format(receptorSource,receptorDest)
    os.popen(receptorCmd)

    # Crystal
    crystalSource = os.path.join(sourceBase,receptor,'crystal_ligand.mol2')
    crystalDestFolder = os.path.join(destBase,'crystal_ligands',receptor)
    createFolder(crystalDestFolder)

    crystalDest = os.path.join(crystalDestFolder,receptor+'_crystal.pdb')
    crystalCmd = 'obabel -imol2 {} -opdb -O {}'.format(crystalSource,crystalDest)
    #os.popen(crystalCmd)


    # Actives
    activesGz = os.path.join(sourceBase,receptor,'actives_final.mol2.gz')
    activesMol2 = gunzip(activesGz)
    if not activesMol2:
        return None
    activesDestFolder = os.path.join(destBase,'ligands',receptor)
    createFolder(activesDestFolder)
    activesDestFile = os.path.join(activesDestFolder,receptor+'_actives_.pdb')
    activesCmd = 'obabel -imol2 {} -opdb -O {} -m'.format(activesMol2,activesDestFile)
    os.popen(activesCmd)

    # Decoys
    decoysGz = os.path.join(sourceBase,receptor,'decoys_final.mol2.gz')
    decoysMol2 = gunzip(decoysGz)
    if not decoysMol2:
        return None
    decoysDestFolder = os.path.join(destBase,'ligands',receptor)
    createFolder(decoysDestFolder)
    decoysDestFile = os.path.join(decoysDestFolder,receptor+'_decoys_.pdb')
    decoysCmd = 'obabel -imol2 {} -opdb -O {} -m'.format(decoysMol2,decoysDestFile)
    print decoysCmd
    #os.popen(decoysCmd)

def runs():
    map(lambda receptor:convert(receptor),os.listdir(sourceBase))

def run_single(receptor):
    convert(receptor)

def docking():
    pass

def docking_offset(offset):
    print "docking ",offset
    df = pd.read_csv(dockingIndex)
    remain = []
    if offset:
        index = range(offset,len(df)+1,1000)
    else:
        index = range(1,len(df)+1)
    print "index size ",len(index)
    
    for i in index:
        receptor, ligand, ligand_box, output = df.ix[i-1]
        command = "{} -r {} -l {} --autobox_ligand {} -o {} --num_modes=1000 --energy_range=100 --cpu=1 " \
            .format(smina, receptor, ligand, ligand_box, output)
        createFolder(os.path.dirname(output))
        print "dock..."
        if not os.path.exists(output):
            if offset:
                os.popen(command)
                #print command
            else:
                remain.append([receptor,ligand,ligand_box,output])

    # when no offset, store the index for all unconverted ligands' information
    if not offset:
        remain_df = pd.DataFrame(data=remain, columns=['receptor', 'ligand', 'ligand_box', 'output'])
        remain_df.to_csv('/n/scratch2/xl198/dude/frame/remain.csv', index=False)



if __name__=='__main__':
    receptors = os.listdir(sourceBase)
    index = int(sys.argv[1])
    run_single(receptors[index])
    sys.exit()
    createBaseFolder()
    receptors = os.listdir(sourceBase)
    args = sys.argv
    if len(args)<2:
        print "not enouth args, need 2 input {}".format(len(args))
    else:
       offset = int(args[1])
       docking_offset(offset)


