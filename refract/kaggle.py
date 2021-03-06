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
from itertools import chain
import prody
from av2_atomdict import atom_dictionary
import getopt

#############################
#
#   according to the input csv file,
#   divide data into trainset and testset.
#   if needed, coded testset data into random code
#   use obabel conver them directly
#
#
#############################

class FLAGS:
    orchestra_arrayjob = False

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
        '''
        Create folder for kaggle dataset
        tempFolderPath used to store intermediate csv file
        kaggleBasePath used to store converted and coded data
        :param folerName:
        '''
        self.basePath = os.path.join(self.workplace,folerName)

        self.tempFolderPath = os.path.join(self.basePath,'temp')

        self.kaggleBasePath = os.path.join(self.basePath,'kaggle')

        create_chain_folder(self.tempFolderPath)
        create_chain_folder(self.kaggleBasePath)


    def write_dataframe(self,dataframe,filename):
        '''
        Write dataframe to tempFolderPath
        :param dataframe: pandas DataFrame
        :param filename:  name of output file without suffix
        :return:
        '''
        dataframe.to_csv(os.path.join(self.tempFolderPath,filename+'.csv'),index=False)

    def generate_database_index(self,dataBasePath):
        '''
        build database index by scan a folder
        and save index file
        :param dataBasePath:
        :return: pandas.DataFrame
        '''

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
        '''
        parse database from csv file, each entry of database
        should have columns [PDBname] and [PDBResId]
        PDBname is the name of receptor e.g. 1ab3
        PDBResId is the combinate of Res number and Frame Id e.g. 1234_1
        file name can be derived by these value.



        :param sourceCsv: path of csv file
        :return:
        '''

        df = pd.read_csv(sourceCsv)

        # select the part that we need to reduce dataframe size.
        df = df[['PDBname','PDBResId']]
        df['RES'] = df['PDBResId'].apply(lambda x:x.split('_')[0])
        df['FrameId'] = df['PDBResId'].apply(lambda x:x.split('_')[1])

        # devided dataset into train set and test set
        trainset,testset = self.divide_by_PDBname(df,testset_ratio=0.1)
        trainset = trainset.merge(df)
        testset  = testset.merge(df)


        # generate source file path and dest filepath
        # for receptor and crystal ligands we dont have a individual entry
        # in dataframe so we generate another dataframe for them
        train_set,train_receptor_crystal = self.train_set(trainset)

        # don't code test set at this time
        test_set,test_receptor_crystal = self.test_set(testset)

        # write down all the dataframe into temp folder
        self.write_dataframe(train_set,'train_set')
        self.write_dataframe(test_set,'test_set')
        self.write_dataframe(train_receptor_crystal,'train_receptor_crystal')
        self.write_dataframe(test_receptor_crystal,'test_receptor_crystal')

        #self.write_dataframe(coded_testset,'coded_testset')

    def receptor_crystal_convert(self,item,coded):
        '''
        get receptor and crystal ligand to dest folder


        :param item: row of dataframe
        :param coded: bool
        :return:
        '''

        receptor_source_file_path   = item['receptor_sourcepath']
        crystal_source_file_path    = item['crystal_sourcepath']
        if coded:
            receptor_dest_file_path = item['code_receptor_destpath']
            crystal_dest_file_path  = item['code_crystal_destpath']
        else:
            receptor_dest_file_path = item['receptor_destpath']
            crystal_dest_file_path  = item['crystal_destpath']


        # copying them from source to dest
        crystal_cmd = "cp {} {}".format(crystal_source_file_path, crystal_dest_file_path)
        receptor_cmd = "cp {} {}".format(receptor_source_file_path,receptor_dest_file_path)

        create_chain_parent_folder(crystal_dest_file_path)
        create_chain_parent_folder(receptor_dest_file_path)

        if not os.path.exists(crystal_dest_file_path):
            os.popen(crystal_cmd)
        if not os.path.exists(receptor_dest_file_path):
            os.popen(receptor_cmd)

    def entry_convert(self,item,coded):
        '''
        original pdb file contains multiple frames
        use obabel for select one frame from it

        :param item: row of dataframe
        :param coded: bool
        :return:
        '''
        source_file_path = item["SourcePath"]
        if coded:
            dest_file_path = item['code_destpath']
        else:
            dest_file_path = item['DestPath']
        Id = item['FrameId']


        cmd = 'obabel -ipdb {} -f {} -l {} -opdb -O {}'.format(source_file_path,Id,Id,dest_file_path)

        try:
            create_chain_parent_folder(dest_file_path)
        except:
            pass
        if not os.path.exists(dest_file_path):
            os.popen(cmd)

    def thread_convert(self,func,dataframe,coded,index):
        '''
        job running in thread
        :param func: function running in single thread
        :param dataframe: pandas.DataFrame
        :param coded: bool
        :param index: list of int
        :return:
        '''
        for i in index:
            #print "dataframe size: {}, ix {}".format(len(dataframe),i)
            func(dataframe.iloc[i],coded)

    def process_convert(self,func,dataframe,coded,index):

        # linspace contain end value but range don't
        # so we use edge[i+1] to select value in index
        # end should be len(index)-1
        

        if len(index)<self.thread_num:
            for i in index:
                func(dataframe.iloc[i],coded)
            return 
        

        edge = np.linspace(0,len(index)-1,self.thread_num+1).astype(int)

        thread_list = [ threading.Thread(target=self.thread_convert,
                                         args=(func,
                                               dataframe,
                                               coded,
                                               range(index[edge[i]],
                                                     index[edge[i+1]])))
                        for i in range(self.thread_num) ]

        for t in thread_list:
            #print "thread start: ",t
            t.start()

        for t in thread_list:
            t.join()

    def convert(self,dataframe,coded=False,is_docked=True,rm_overlap=False):
        '''
        according the result of 'database_from_csv'
        running multiprocess to get result

        :param dataframe: pandas.DataFrame
                          string
        :param coded: bool
        :return:
        '''


        convert_func = self.entry_convert if is_docked else self.receptor_crystal_convert

        if rm_overlap:
            convert_func = self.remove_overlap

        if hasattr(FLAGS,'cores'):
            self.process_num = FLAGS.cores

        # if get a str, read csv
        if type(dataframe) == str:
            try:
                dataframe = pd.read_csv(os.path.join(self.tempFolderPath,dataframe))
            except Exception as e:
                print e

        if FLAGS.orchestra_arrayjob:
            jobsize = FLAGS.orchestra_jobsize
            jobid = FLAGS.orchestra_jobid

            # using iloc to slice dataframe result doesn't contains end [start,end)
            linspace = np.linspace(0,len(dataframe),jobsize+1).astype(int)
            dataframe = dataframe.iloc[linspace[jobid-1]:linspace[jobid]]

        print "dataframe size ",len(dataframe)
        
        # when there's not enough entry to comvert , decrease thread's num
        if len(dataframe)<self.process_num*self.thread_num:
            
            for i in range(len(dataframe)):
                convert_func(dataframe.iloc[i],coded)
            return
        edge = np.linspace(0,len(dataframe),self.process_num+1).astype(int)
        process_list = [ multiprocessing.Process(target=self.process_convert,
                                                 args=(convert_func,
                                                       dataframe,
                                                       coded,
                                                       range(edge[i],
                                                             edge[i+1])))
                         for i in range(self.process_num)]

        for p in process_list:
            print "process start: ",p
            p.start()

        for p in process_list:
            print "process end: ",p
            p.join()

    def train_set(self,trainset):


        trainset['SourcePath'] = trainset.apply(
            lambda item: os.path.join(self.source_base,
                                      item['PDBname'],
                                      '_'.join([item['PDBname'],
                                               item['RES'],
                                               'ligand',
                                               'fast.pdb'])),
            axis=1
        )

        trainset['DestPath'] = trainset.apply(
            lambda item: os.path.join(self.kaggleBasePath,
                                      "train_pdb",
                                      self.ligandFolderName,
                                      item['PDBname'],
                                      '_'.join([item['PDBname'],
                                               item['RES'],
                                               'ligand',
                                               'fast',
                                               item['FrameId'] + '.pdb'])),
            axis=1
        )

        crystalLigands = list(set(zip(trainset['PDBname'], trainset['RES'])))
        receptorAndCrystal = pd.DataFrame(data=crystalLigands, columns=['PDBname', 'RES'])

        receptorAndCrystal['receptor_sourcepath'] = receptorAndCrystal.apply(
            lambda item: os.path.join(self.receptor_base,
                                      item['PDBname'],
                                      item['PDBname'] + '.pdb'),
            axis=1
        )

        receptorAndCrystal['receptor_destpath'] = receptorAndCrystal.apply(
            lambda item: os.path.join(self.kaggleBasePath,
                                      "train_pdb",
                                      self.receptorFolderName,
                                      item['PDBname']+'.pdb'),
            axis = 1
        )

        receptorAndCrystal['crystal_sourcepath'] = receptorAndCrystal.apply(
            lambda item: os.path.join(self.crystal_base,
                                      item['PDBname'],
                                      '_'.join([item['PDBname'],
                                               item['RES'],
                                               'ligand.pdb'])),
            axis=1
        )

        receptorAndCrystal['crystal_destpath'] = receptorAndCrystal.apply(
            lambda item: os.path.join(self.kaggleBasePath,
                                      "train_pdb",
                                      self.crystalFolderName,
                                      item['PDBname'],
                                      '_'.join([item['PDBname'],
                                               item['RES'],
                                               'ligand.pdb'])),
            axis=1
        )

        docked_dest = trainset[['PDBname', 'RES', 'DestPath']]
        crystal_dest = receptorAndCrystal[['PDBname', 'RES', 'crystal_destpath']]
        docked_crystal_pair = docked_dest.merge(crystal_dest, on=['PDBname', 'RES'])
        self.write_dataframe(docked_crystal_pair, 'train_docked_crystal_pair')

        return trainset,receptorAndCrystal

    def test_set(self, testset):


        testset['SourcePath'] = testset.apply(
            lambda item: os.path.join(self.source_base,
                                      item['PDBname'],
                                      '_'.join([item['PDBname'],
                                               item['RES'],
                                               'ligand',
                                               'fast.pdb'])),
            axis=1
        )

        testset['DestPath'] = testset.apply(
            lambda item: os.path.join(self.kaggleBasePath,
                                      "test_pdb",
                                      self.ligandFolderName,
                                      item['PDBname'],
                                      '_'.join([item['PDBname'],
                                               item['RES'],
                                               'ligand',
                                               'fast',
                                               item['FrameId'] + '.pdb'])),
            axis=1
        )

        crystalLigands = list(set(zip(testset['PDBname'], testset['RES'])))
        receptorAndCrystal = pd.DataFrame(data=crystalLigands, columns=['PDBname', 'RES'])

        receptorAndCrystal['receptor_sourcepath'] = receptorAndCrystal.apply(
            lambda item: os.path.join(self.receptor_base,
                                      item['PDBname'],
                                      item['PDBname'] + '.pdb'),
            axis=1
        )

        receptorAndCrystal['receptor_destpath'] = receptorAndCrystal.apply(
            lambda item: os.path.join(self.kaggleBasePath,
                                      "test_pdb",
                                      self.receptorFolderName,
                                      item['PDBname']+'.pdb'),
            axis = 1
        )

        receptorAndCrystal['crystal_sourcepath'] = receptorAndCrystal.apply(
            lambda item: os.path.join(self.crystal_base,
                                      item['PDBname'],
                                      '_'.join([item['PDBname'],
                                               item['RES'],
                                               'ligand.pdb'])),
            axis=1
        )

        receptorAndCrystal['crystal_destpath'] = receptorAndCrystal.apply(
            lambda item: os.path.join(self.kaggleBasePath,
                                      "test_pdb",
                                      self.crystalFolderName,
                                      item['PDBname'],
                                      '_'.join([item['PDBname'],
                                               item['RES'],
                                               'ligand.pdb'])),
            axis=1
        )

        docked_dest = testset[['PDBname', 'RES', 'DestPath']]
        crystal_dest = receptorAndCrystal[['PDBname', 'RES', 'crystal_destpath']]
        docked_crystal_pair = docked_dest.merge(crystal_dest, on=['PDBname', 'RES'])
        self.write_dataframe(docked_crystal_pair, 'test_docked_crystal_pair')

        return testset, receptorAndCrystal

    def code_test_set(self,testset,coded=True):
        '''
        code the test set into random number label, and mix crystal ligands inside
        :param testset: pandas.DataFrame test_set
        :return:
        '''

        testset['SourcePath'] = testset.apply(
            lambda item: os.path.join(self.source_base,
                                      item['PDBname'],
                                      '_'.join([item['PDBname'],
                                               item['RES'],
                                               'ligand',
                                               'fast.pdb'])),
            axis=1
        )

        testset['DestPath'] = testset.apply(
            lambda item: os.path.join(self.kaggleBasePath,
                                      "unlabeled_pdb",
                                      "ligands",
                                      item["PDBname"],
                                      '_'.join([item['PDBname'],
                                               item['RES'],
                                               'ligand',
                                               'fast',
                                               item['FrameId'] + '.pdb'])),
            axis=1
        )

        crystalLigands = list(set(zip(testset['PDBname'],testset['RES'])))
        receptorAndCrystal = pd.DataFrame(data=crystalLigands, columns=['PDBname','RES'])



        if coded:
            # assign a code for each docked liagnd and crystal ligand
            # assign a code for receptor in each recptor-crystal ligand pair
            # usually a receptor pair with multiple crystal ligand,
            # it will have more than one code

            code = np.arange(len(testset)+len(receptorAndCrystal)*2)
            random.shuffle(code)

            receptorAndCrystal['receptor_code'] = code[:len(receptorAndCrystal)]
            receptorAndCrystal['crystal_code']  = code[len(receptorAndCrystal):len(receptorAndCrystal)*2]
            testset['ligand_code']              = code[len(receptorAndCrystal)*2:]

            testset = testset.merge(receptorAndCrystal,on=['PDBname','RES'])

            testset['code_destpath'] = testset.apply(
                lambda item:os.path.join(self.kaggleBasePath,
                                         "unlabeled_pdb",
                                         "ligands",
                                         'P'+str(item['receptor_code']),
                                         'L'+str(item['ligand_code'])+'.pdb'),
                axis=1
            )




            receptorAndCrystal['code_receptor_destpath'] = receptorAndCrystal.apply(
                lambda item:os.path.join(self.kaggleBasePath,
                                         "unlabeled_pdb",
                                         "receptors",
                                         'P'+str(item['receptor_code'])+'.pdb'),
                axis = 1
            )




            receptorAndCrystal['code_crystal_destpath'] = receptorAndCrystal.apply(
                lambda item:os.path.join(self.kaggleBasePath,
                                         "unlabeled_pdb",
                                         "ligands",
                                         'P'+str(item['receptor_code']),
                                         'L'+str(item['crystal_code'])+'.pdb'),
                axis=1
            )

        receptorAndCrystal['receptor_sourcepath'] = receptorAndCrystal.apply(
            lambda item: os.path.join(self.receptor_base,
                                      item['PDBname'],
                                      item['PDBname'] + '.pdb'),
            axis=1
        )

        receptorAndCrystal['crystal_sourcepath'] = receptorAndCrystal.apply(
            lambda item: os.path.join(self.crystal_base,
                                      item['PDBname'],
                                      '_'.join([item['PDBname'],
                                                item['RES'],
                                                'ligand.pdb'])),
            axis=1
        )


        print "test set -------------------------------------"
        print testset.columns

        docked_dest= testset[['PDBname','RES','code_destpath']]
        crystal_dest = receptorAndCrystal[['PDBname','RES','code_crystal_destpath']]
        docked_crystal_pair = docked_dest.merge(crystal_dest,on=['PDBname','RES'])
        self.write_dataframe(docked_crystal_pair,'test_docked_crystal_pair')




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

    def PDB_2_npy(self,sourceFilePath,coded,is_receptor):
        '''
        convert pdb into npy, only need atom coordinate and atom type
        :param sourceFilePath:
        :param coded: if file coded
        :param is_receptor:
        :return:
        '''
        atom_dict = atom_dictionary.REC if is_receptor else atom_dictionary.LIG

        if coded:
            destFilePath = re.sub('unlabeled_pdb', 'unlabeled_npy', sourceFilePath)
            destFilePath = re.sub('.pdb$', '', destFilePath)
        else:
            destFilePath = re.sub('labeled_pdb', 'labeled_npy', sourceFilePath)
            destFilePath = re.sub('.pdb$', '', destFilePath)

        try:
            prody_ligand = prody.parsePDB(sourceFilePath)
        except Exception as e:
            print e

        try:
            atom_numbers = map(lambda atomname: atom_dict[atomname.lower()],
                               prody_ligand.getElements())
            coordinates_and_atoms = np.hstack((prody_ligand.getCoords(),
                                               np.reshape(atom_numbers, (-1, 1))))
            create_chain_parent_folder(destFilePath)
            np.save(destFilePath,coordinates_and_atoms)
        except Exception as e:
            print e

    def clashdetect(self,item,coded):
        '''
        if docked ligand doesn't overlap with crystal ligand
        convert docked ligand into npy format
        :param item: row of dataframe
        :param coded: bool
        :return:
        '''
        if coded:
            docked_ligand_path = item['code_destpath']
            crystal_ligand_path = item['code_crystal_destpath']
        else:
            docked_ligand_path = item['DestPath']
            crystal_ligand_path =item['crystal_destpath']
        if os.path.exists(docked_ligand_path) and os.path.exists(crystal_ligand_path):
            try:
                if not docked_ligand_overlaps_with_crystal(docked_ligand_path,crystal_ligand_path):
                
                    self.PDB_2_npy(docked_ligand_path,coded,is_receptor=False)
            except:
                pass

    def remove_overlap(self, item, coded):
        '''
        if docked ligand doesn't overlap with crystal ligand
        convert docked ligand into npy format
        :param item: row of dataframe
        :param coded: bool
        :return:
        '''
        if coded:
            docked_ligand_path = item['code_destpath']
            crystal_ligand_path = item['code_crystal_destpath']
        else:
            docked_ligand_path = item['DestPath']
            crystal_ligand_path = item['crystal_destpath']
        if os.path.exists(docked_ligand_path) and os.path.exists(crystal_ligand_path):
            try:
                if docked_ligand_overlaps_with_crystal(docked_ligand_path, crystal_ligand_path):
                    print "remove"
                    os.popen('rm {}'.format(docked_ligand_path))
                else:
                    print "save"
            except:
                
                pass



    def process_PDB_to_npy(self,dataframe,coded):
        '''
        convert pdb data in /labeled_pdb and /unlabeled_pdb into npy form
        /labeled_npy and /unlabeled_npy


        directly convert /receptors , /crystal_ligands and /ligands into npy
        /docked_ligands need to remove the overlap result first


        :param dataframe: pandas.DataFrame : dataframe contain docked ligand
                                            and corresponding crystal ligand path
                          string : path of csv file
        :param coded: bool, if the result is coded
        :return:
        '''

        if type(dataframe) == str:
            try:
                dataframe = pd.read_csv(dataframe)
            except Exception as e:
                print e

        if hasattr(FLAGS,'cores'):
            self.process_num = FLAGS.cores

        print "running in {} processors".format(self.process_num)

        if FLAGS.orchestra_arrayjob:
            jobsize = FLAGS.orchestra_jobsize
            jobid = FLAGS.orchestra_jobid

            # using iloc to slice dataframe result doesn't contains end [start,end)
            linspace = np.linspace(0,len(dataframe),jobsize+1).astype(int)
            dataframe = dataframe.iloc[linspace[jobid-1]:linspace[jobid]]

        print "dataframe size ",len(dataframe)


        if coded:
            sourcePath = os.path.join(self.kaggleBasePath,'unlabeled_pdb')
        else:
            sourcePath = os.path.join(self.kaggleBasePath,'labeled_pdb')

        edge = np.linspace(0,len(dataframe),self.process_num+1).astype(int)
        process_list = [ multiprocessing.Process(target=self.process_convert,
                                                 args=(self.clashdetect,
                                                       dataframe,
                                                       coded,
                                                       range(edge[i],
                                                             edge[i+1])))
                         for i in range(self.process_num)]

        for p in process_list:
            print "process start: ",p
            p.start()


        if coded:
            for sourceFilePath in path_iterator(os.walk(os.path.join(sourcePath,self.receptorFolderName))):
                self.PDB_2_npy(sourceFilePath,coded,is_receptor=True)
            for sourceFilePath in list(set(dataframe['code_crystal_destpath'])):
                self.PDB_2_npy(sourceFilePath,coded,is_receptor=False)

        else:
            pass
            for sourceFilePath in path_iterator(os.walk(os.path.join(sourcePath,self.receptorFolderName))):
                self.PDB_2_npy(sourceFilePath,coded,is_receptor=True)

            for sourceFilePath in path_iterator(os.walk(os.path.join(sourcePath,self.crystalFolderName))):
                self.PDB_2_npy(sourceFilePath,coded,is_receptor=False)




        for p in process_list:
            print "process join: ",p
            p.join()




def parse_FLAG():
    try:
        opts,args = getopt.getopt(sys.argv[1:],None,["jobsize=","jobid=","cores="])
    except getopt.GetoptError as err:
        # print help information and exit:
        print str(err)  # will print something like "option -a not recognized"

        sys.exit(2)

    for name,value in opts:
        if name == '--jobsize':
            FLAGS.orchestra_jobsize = int(value)
            print "--jobsize ",value
        if name == '--jobid':
            FLAGS.orchestra_jobid = int(value)
            print "--jobid", value
        if name == '--cores':
            FLAGS.cores = int(value)

    if hasattr(FLAGS,"orchestra_jobsize") and hasattr(FLAGS,"orchestra_jobid"):
        FLAGS.orchestra_arrayjob = True

    print "orchestra job ",FLAGS.orchestra_arrayjob

    if hasattr(FLAGS,'cores'):
        print "cores num ",FLAGS.cores

def get_pdb():
    kaggle = kaggleDataset('jan_18')
    kaggle.database_from_csv('/home/xl198/remark/dec_17_small.csv')
    kaggle.convert('train_set.csv')
    kaggle.convert('train_receptor_crystal.csv', is_docked=False)
    kaggle.convert('test_set.csv', coded=True)
    kaggle.convert('test_receptor_crystal.csv', coded=True, is_docked=False)

def get_npy():
    kaggle = kaggleDataset('jan_13_big')
    kaggle.process_PDB_to_npy('/n/scratch2/xl198/data/jan_13_big/temp/train_docked_crystal_pair.csv',coded=False)
    kaggle.process_PDB_to_npy('/n/scratch2/xl198/data/jan_13_big/temp/test_docked_crystal_pair.csv', coded=True)


if __name__ == '__main__':
    parse_FLAG()
    kaggle = kaggleDataset('jan_18_big')
    kaggle.convert('train_docked_crystal_pair.csv',rm_overlap=True)
    kaggle.convert('test_docked_crystal_pair.csv',rm_overlap=True)
    #kaggle.database_from_csv('/home/xl198/remark/dec_17_nodude.csv')
    #kaggle.convert('train_set.csv')
    #kaggle.convert('train_receptor_crystal.csv', is_docked=False)
    #kaggle.convert('test_set.csv',)
    #kaggle.convert('test_receptor_crystal.csv', is_docked=False)



