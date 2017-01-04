import pandas as pd
import numpy as np
import sys,os
import re


    # get folder path for docked file and receptor's name
    # run obabel to convert the first result into single file
    # rename the converted one
    # move crystal ligand and receptor into folder

sourceBase = '/n/scratch2/xl198/dude/data/all'
destSource = '/n/scratch2/xl198/dude/data/jan_04/pdbs'

destFolders = ['crystal_ligands',  'docked_ligands' , 'receptors']

def createFolder(folderPath):
    if not os.path.exists(folderPath):
        os.mkdir(folderPath)


def prepareTarget(filePath,receptor):
    baseName = os.path.basename(filePath)
    baseTarget = re.sub('\.','_.',baseName)
    targetName = receptor+'_'+baseTarget
    targetFolder = os.path.join(destSource,'docked_ligands',receptor)
    createFolder(targetFolder)
    targetFile = os.path.join(targetFolder,targetName)


    return targetFile

def getSourceFolder():
    receptor_list = ['cxcr4', 'comt', 'inha', 'fabp4', 'ampc']
    return map(lambda x:[x,os.path.join(sourceBase,x,'docked')],receptor_list)

def getReceptor():
    receptorFolder = os.path.join(destSource,'receptors')
    createFolder(receptorFolder)
    receptor_list = ['cxcr4', 'comt', 'inha', 'fabp4', 'ampc']
    source = map(lambda x:os.path.join(sourceBase,x,'receptor.pdb'),receptor_list)
    dest = map(lambda x:os.path.join(destSource,'receptors',x+'.pdb'),receptor_list)
    pairs = zip(source,dest)

    cmds = map(lambda x:"obabel -ipdb {} -opdb -O  {}".format(x[0],x[1]),pairs)
    map(lambda x:os.popen(x),cmds)




def convert():
    folders = getSourceFolder()
    for receptor,folder in folders:
        for dirpath,dirname,filenames in os.walk(folder):
            for filename in filenames:
                filePath = os.path.join(dirpath,filename)
                
                targetPath = prepareTarget(filePath,receptor)
                for i in range(1,11):
                    numberTarget = re.sub('\.',str(i)+'.',targetPath)
                    command = 'obabel -ipdb {} -f {} -l {} -opdb -O {} '.format(filePath,i,i,numberTarget)
                    if not os.path.exists(numberTarget):
                        os.popen(command)

convert()
getReceptor()
