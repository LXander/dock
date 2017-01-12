import pandas as pd
import numpy as np
import random
import os,sys,re
import shutil
from functools import partial
from util import *
from obabel import obabelOperation
import multiprocessing
import threading

#############################
#
#   according to the input csv file,
#   divide data into trainset and testset.
#   if needed, coded testset data into random code
#   use obabel conver them directly
#
#
#############################

class kaggleDataset:
    source_base = '/n/scratch2/xl198/data/pdbs'
    crystal_base = '/n/scratch2/xl198/data/H/addH'
    receptor_base = '/n/scratch2/xl198/data/H/data'
    workplace = '/n/scratch2/xl198/data'
    ligandFolderName   = 'docked_ligands'
    crystalFolderName  = 'crystal_ligands'
    receptorFolderName = 'receptors'
    thread_num = 16
    process_num = 12

    def __init__(self,folerName):
        self.basePath = os.path.join(self.workplace,folerName)

        self.tempFolderPath = os.path.join(self.basePath,'temp')

        self.kaggleBasePath = os.path.join(self.basePath,'kaggle')

        create_chain_folder(self.tempFolderPath)
        create_chain_folder(self.kaggleBasePath)


    def write_dataframe(self,dataframe,filename):
        dataframe.to_csv(os.path.join(self.tempFolderPath,filename+'.csv'),index=False)

    def generate_database_index(self,dataBasePath):

        ligandsFolderPath  = os.path.join(dataBasePath,self.ligandFolderName)
        crystalFolderPath  = os.path.join(dataBasePath,self.crystalFolderName)
        receptorFolderPath = os.path.join(dataBasePath,self.receptorFolderName)

        # collect information for docked ligands
        raw = []
        for dirpath, dirname, filenames in os.walk(ligandsFolderPath):
            for filename in filenames:
                row = []
                # ID abcd_123
                row.append('_'.join(filename.split('_')[:2]))
                # Preceptor_folder_named
                row.append(filename.split('_')[0])
                # path
                row.append(os.path.join(self.ligandFolderName, filename.split('_')[0], filename.split('.')[0]))
                # innder_index
                row.append(int(filename.split('.')[0].split('_')[-1]))
                raw.append(row)
        dataBase = pd.DataFrame(data=raw, columns=['ID', 'PDBname', 'path', 'inner_index'])

        # filter dataBase by crystal ligands
        cry = []
        for dirpath, dirname, filenames in os.walk(crystalFolderPath):
            for filename in filenames:
                cry.append('_'.join(filename.split('_')[:2]))
        crystal = pd.DataFrame(data=cry, columns=['ID'])
        dataBase = dataBase.merge(crystal)

        # filter dataBase by receptor
        rec = []
        for dirpath, dirname, filenames in os.walk(receptorFolderPath):
            for filename in filenames:
                rec.append(filename.split('.')[0])
        receptor = pd.DataFrame(data=rec, columns=['PDBname'])
        dataBase = dataBase.merge(receptor)

        self.write_dataframe(dataBase,'database_index')

        return dataBase

    def database_from_csv(self,sourceCsv):

        df = pd.read_csv(sourceCsv)
        # select the port that we needed to reduce dataframe size.
        df = df[['PBDname','PDBResId']]
        df['RES'] = df['PDBResId'].apply(lambda x:x.split('_')[0])
        df['FrameId'] = df['PDBResId'].apply(lambda x:x.split('_')[1])

        trainset,testset = self.divide_by_PDBname(df,testset_ratio=0.1)

        trainset = trainset.merge(df)
        testset  = testset.merge(df)


        #self.write_dataframe(trainset,'train_set')
        #self.write_dataframe(testset,'test_set')

        train_set,train_receptor_crystal = self.train_set(trainset)

        test_set,test_receptor_crystal = self.code_test_set(testset,True)

        self.write_dataframe(train_set,'train_set')
        self.write_dataframe(test_set,'test_set')
        self.write_dataframe(train_receptor_crystal,'train_receptor_crystal')
        self.write_dataframe(test_receptor_crystal,'test_receptor_crystal')

        #self.write_dataframe(coded_testset,'coded_testset')

    def receptor_crystal_convert(self,item,coded):

        receptor_source_file_path   = item['receptor_sourcepath']
        crystal_source_file_path    = item['crystal_sourcepath']
        if coded:
            receptor_dest_file_path = item['code_crystal_destpath']
            crystal_dest_file_path  = item['code_crystal_destpath']
        else:
            receptor_dest_file_path = item['receptor_destpath']
            crystal_dest_file_path  = item['crystal_destpath']

        crystal_cmd = "cp {} {}".format(crystal_source_file_path, crystal_dest_file_path)
        receptor_cmd = "cp {} {}".format(receptor_source_file_path,receptor_dest_file_path)

        create_chain_parent_folder(crystal_dest_file_path)
        create_chain_parent_folder(receptor_dest_file_path)

        if not os.path.exists(crystal_dest_file_path):
            os.popen(crystal_cmd)
        if not os.path.exists(receptor_dest_file_path):
            os.popen(receptor_cmd)

    def entry_convert(self,item,coded):

        source_file_path = item["SourcePath"]
        if coded:
            dest_file_path = item['Coded_destpath']
        else:
            dest_file_path = item['DestPath']
        Id = item['FrameId']

        cmd = 'obabel -ipdb %s -f %s -opdb -O %s'.format(source_file_path,Id,Id,dest_file_path)

        create_chain_parent_folder(dest_file_path)
        if not os.path.exists(dest_file_path):
            os.popen(cmd)

    def thread_convert(self,func,dataframe,coded,index):
        for i in index:
            func(dataframe[i],coded)

    def process_convert(self,func,dataframe,coded,index):
        edge = np.linspace(0,len(index),self.thread_num+1)
        thread_list = [ threading.Thread(target=self.thread_convert,
                                         args=(func,
                                               dataframe,
                                               coded,
                                               range(index[edge[i]],
                                                     index[edge[i+1]])))
                        for i in range(self.thread_num) ]

        for t in thread_list:
            print "thread start: ",t
            t.start()

        for t in thread_list:
            t.join()

    def convert(self,dataframe,coded=False):
        edge = np.linspace(0,len(dataframe),self.process_num+1)
        process_list = [ multiprocessing.Process(target=self.process_convert,
                                                 args=(self.entry_convert,
                                                       dataframe,
                                                       coded,
                                                       range(edge[i],
                                                             edge[i+1])))
                         for i in range(self.process_num)]

        for p in process_list:
            print "process start: ",p
            p.start()

        for p in process_list:
            p.join()

    def train_set(self,trainset):


        trainset['SourcePath'] = trainset.apply(
            lambda item: os.path.join(self.source_base,
                                      item['PDBname'],
                                      '_'.join(item['PDBname'],
                                               item['RES'],
                                               'ligand',
                                               'fast.pdb')),
            axis=1
        )

        trainset['DestPath'] = trainset.apply(
            lambda item: os.path.join(self.kaggleBasePath,
                                      "labeled_pdb",
                                      self.ligandFolderName,
                                      item['PDBname'],
                                      '_'.join(item['PDBname'],
                                               item['RES'],
                                               'ligand',
                                               'fast',
                                               item['FrameId'] + '.pdb')),
            asix=1
        )

        crystalLigands = list(set(zip(trainset['PDBname'], trainset['RES'])))
        receptorAndCrystal = pd.dataFrame(data=crystalLigands, columns=['PDBname', 'RES'])

        receptorAndCrystal['receptor_sourcepath'] = receptorAndCrystal.apply(
            lambda item: os.path.join(self.receptor_base,
                                      receptorAndCrystal['PDBname'],
                                      receptorAndCrystal['PDBname'] + '.pdb'),
            axis=1
        )

        receptorAndCrystal['receptor_destpath'] = receptorAndCrystal.apply(
            lambda item: os.path.join(self.kaggleBasePath,
                                      "labeled_pdb",
                                      self.receptorFolderName,
                                      item['PDBname']+'.pdb'),
            axis = 1
        )

        receptorAndCrystal['crystal_sourcepath'] = receptorAndCrystal.apply(
            lambda item: os.path.join(self.crystal_base,
                                      receptorAndCrystal['PDBname'],
                                      '_'.join(receptorAndCrystal['PDBname'],
                                               receptorAndCrystal['RES'],
                                               'ligand.pdb')),
            axis=1
        )

        receptorAndCrystal['crystal_destpath'] = receptorAndCrystal.apply(
            lambda item: os.path.join(self.kaggleBasePath,
                                      "labeled_pdb",
                                      self.crystalFolderName,
                                      item['PDBname'],
                                      '_'.join(item['PDBname'],
                                               item['RES'],
                                               'ligand.pdb'))
        )

        return trainset,receptorAndCrystal

    def code_test_set(self,testset,coded=False):
        '''
        code the test set into random number label, and mix crystal ligands inside
        :param testset: pandas.DataFrame test_set
        :return:
        '''

        testset['SourcePath'] = testset.apply(
            lambda item: os.path.join(self.source_base,
                                      item['PDBname'],
                                      '_'.join(item['PDBname'],
                                               item['RES'],
                                               'ligand',
                                               'fast.pdb')),
            axis=1
        )

        testset['DestPath'] = testset.apply(
            lambda item: os.path.join(self.kaggleBasePath,
                                      "unlabeled_pdb",
                                      "ligands",
                                      item["PDBname"],
                                      '_'.join(item['PDBname'],
                                               item['RES'],
                                               'ligand',
                                               'fast',
                                               item['FrameId'] + '.pdb')),
            axis=1
        )

        crystalLigands = list(set(zip(testset['PDBname'],testset['RES'])))
        receptorAndCrystal = pd.dataFrame(data=crystalLigands, columns=['PDBname','RES'])

        receptorAndCrystal['receptor_sourcepath'] = receptorAndCrystal.apply(
            lambda item: os.path.join(self.receptor_base,
                                      receptorAndCrystal['PDBname'],
                                      receptorAndCrystal['PDBname'] + '.pdb'),
            axis=1
        )


        receptorAndCrystal['crystal_sourcepath'] = receptorAndCrystal.apply(
            lambda item: os.path.join(self.crystal_base,
                                      receptorAndCrystal['PDBname'],
                                      '_'.join(receptorAndCrystal['PDBname'],
                                               receptorAndCrystal['RES'],
                                               'ligand.pdb')),
            axis=1
        )




        if coded:



            code = np.arange(len(testset)+len(receptorAndCrystal)*2)
            random.shuffle(code)

            receptorAndCrystal['receptor_code'] = code[:len(receptorAndCrystal)]
            receptorAndCrystal['crystal_code']  = code[len(receptorAndCrystal):len(receptorAndCrystal)*2]
            testset['ligand_code']              = code[len(receptorAndCrystal)*2:]

            testset['Coded_destpath'] = testset.apply(
                lambda item:os.path.join(self.kaggleBasePath,
                                         "unlabeled_pdb",
                                         "ligands",
                                         'L'+item['ligand_code']+'.pdb'),
                axis=1
            )




            receptorAndCrystal['code_receptor_destpath'] = receptorAndCrystal.apply(
                lambda item:os.path.join(self.kaggleBasePath,
                                         "unlabeledd_pdb",
                                         "receptors",
                                         'P'+item['receptor_code']+'.pdb'),
                axis = 1
            )




            receptorAndCrystal['code_crystal_destpath'] = receptorAndCrystal.apply(
                lambda item:os.path.join(self.kaggleBasePath,
                                         "unlabeled_pdb",
                                         "ligands",
                                         'L'+item['crystal_code']+'.pdb'),
                axis=1
            )


        return testset,receptorAndCrystal

    def divide_by_PDBname(self,dataset,testset_ratio=0.1):
        '''
        Divide data into two part, each of them have differenc receptor pdb('PDBname')

        :param dataset: pandas.DataFrame
        :param testset_ratio: float [0,1)
        :return: list of pdbname
        '''

        pdbs     = list(set(dataset['PDBname']))
        trainSet = []
        testSet  = []

        for pdb in pdbs:
            if random.random() < testset_ratio:
                testSet.append(pdb)
            else:
                trainSet.append(pdb)

        train = pd.DataFrame(data = trainSet,columns=['PDBname'])
        test  = pd.DataFrame(data = testSet ,columns=['PDBname'])
        return train,test


if __name__ == '__main__':
    kaggle = kaggleDataset('jan_11')
    kaggle.database_from_csv('/home/xl198/remark/dec_17_small.csv')
