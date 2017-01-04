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


def output_file(path):
    tmp = re.sub('/actives/','/docked/',path)
    return re.sub('/decoys/','/docked/',tmp)

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
    
    total_list = map(lambda x:[receptor,x,crystal,output_file(x)],actives+decoys)

    return total_list



def create_docking_list():
    lists = map(lambda x:make_list_for_onefolder(x),receptor_list)
    total = reduce(lambda x,y:x+y,lists)
    df = pd.DataFrame(data = total,columns=['receptor','ligand','ligand_box','output'])
    df.to_csv('/n/scratch2/xl198/dude/frame/list.csv',index=False)





smina = '/home/xl198/program/smina/smina.static'

def score(offset):
    df = pd.read_csv("/n/scratch2/xl198/dude/code/crystal.csv")
    indexes = range(len(df))
    result = []
    for i in indexes:
        receptor,ligand,ligand_box,output = df.ix[i]
        command = "{} -r {} -l {} --autobox_ligand {}  --num_modes=1000 --energy_range=100 --cpu=1 --score_only"\
            .format(smina,receptor,ligand,ligand_box)

        l = os.popen(command)
        ls = l.read()
        m = re.search('Affinity: (.*)\(',ls)
        if m:
	    affinity = float(m.groups()[0])
            result.append([os.path.basename(ligand),affinity])
    df = pd.DataFrame(data = result,columns = ['filename','energy'])
    df.to_csv("/n/scratch2/xl198/dude/code/crystal_score.csv",index=False)

def run(offset):
    df = pd.read_csv('/n/scratch2/xl198/dude/code/crystal.csv')
    remain = []
    if offset:
    	indexes = range(offset-1,len(df),1000)
    else:
        indexes = range(len(df))
    for i in indexes:
        receptor,ligand,ligand_box,output = df.ix[i]
        command = "{} -r {} -l {} --autobox_ligand {} -o {} --num_modes=1000 --energy_range=100 --cpu=1 "\
            .format(smina,receptor,ligand,ligand_box,output)

        #print command
        if not os.path.exists(output):
            if offset:
                os.popen(command)
            else:
                remain.append([receptor,ligand,ligand_box,output])
        
        
        

    remain_df = pd.DataFrame(data = remain,columns=['receptor','ligand','ligand_box','output'])
    remain_df.to_csv('/n/scratch2/xl198/dude/frame/remain.csv',index=False)

if __name__ == '__main__':

    argv = sys.argv
    if len(argv)>=2:
        score(int(argv[1]))
        #run(int(argv[1]))
    else:
        print "not enough args, need 2 receive {}".format(len(argv))
        print "e.g. python docking_dude.py 1"



