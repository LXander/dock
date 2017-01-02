import os,sys
from functools import partial
import re
import pandas as pd

'''
This code is used to dock dude dataset

Create three folder to store the file

First extract activate and decoy mol2 file into folder
and then index them.



'''


folders = ['actives','decoys','docked']
receptor_list = ['cxcr4','comt','inha','fabp4','ampc']

def unzip(filename,path):
    filePath = os.path.join(path,filename)
    if os.path.exists(filePath):
        os.popen("gunzip -d {}".format(filePath))

def folderCreate(path):
    for foldername in folders:
        folderPath = os.path.join(path,foldername)
        if not os.path.exists(folderPath):
            os.mkdir(folderPath)

def split_mol2(folderPath):

    source_actives = os.path.join(folderPath,'actives_final.mol2')
    dest_actives = os.path.join(folderPath,'actives','actives_.pdb')

    source_decoys = os.path.join(folderPath,'decoys_final.mol2')
    dest_decoys = os.path.join(folderPath,'decoys','decoys_.pdb')



    os.popen("obabel -imol2 {} -opdb -O {} -m".format(source_actives,dest_actives))
    os.popen("obabel -imol2 {} -opdb -O {} -m".format(source_decoys,dest_decoys))


def get_decoy_list(folder,path):
    folderPath = os.path.join(folder,path)
    unzip('decoys_final.mol2.gz',folderPath)
    folderCreate()




def make_list_for_onefolder(folder,path='/n/scratch2/xl198/dude/data/all'):
    folderPath = os.path.join(path,folder)
    print "folder path",folderPath
    unzip('actives_final.mol2.gz',folderPath)
    unzip('decoys_final.mol2.gz',folderPath)
    folderCreate(folderPath)
    split_mol2(folderPath)
    actives = os.listdir(os.path.join(folderPath,'actives'))
    decoys = os.listdir(os.path.join(folderPath,'decoys'))
    actives = map(lambda x:os.path.join(folderPath,'actives',x),actives)
    decoys =map(lambda x:os.path.join(folderPath,'decoys',x),decoys)

    crystal = os.path.join(folderPath,'crystal_ligand.mol2')
    receptor = os.path.join(folderPath,'receptor.pdb')

    total_list = map(lambda x:[receptor,x,crystal],actives+decoys)

    return total_list



def create_docking_list():
    lists = map(lambda x:make_list_for_onefolder(x),receptor_list)
    total = reduce(lambda x,y:x+y,lists)
    df = pd.DataFrame(data = total,columns=['receptor','ligand','ligand_box'])
    df.to_csv('/n/scratch2/xl198/dude/frame/list.csv',index=False)



create_docking_list()



