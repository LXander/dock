import os,sys
import re
import shutil
from util import *

class obabelOperation:
    crystalSourceBase = '/n/scratch2/xl198/data/H/addH'
    receptorSourceBase = '/n/scratch2/xl198/data/H/data'
    pdbSourceBase ='/n/scratch2/xl198/data/pdbs'
    destBase = '/n/scratch2/xl198/data/'

    def __init__(self,baseFolder):
        '''
        Set dest folder and create it
        :param baseFolder: folder name for this convert eg. "jan_01"

        create folder '/n/scratch2/xl198/data/jan_01'
        create folder '/n/scratch2/xl198/data/jan_01/pdbs'
        '''

        self.destBaseFolder = os.path.join(self.destBase,baseFolder,'pdbs')
        create_parent_folder(self.destBaseFolder)
        create_folder(self.destBaseFolder)

    def get_source_file_path(self,sourcePath):
        '''
        return a iterator to retuen all the file's path under sourcePath
        :param sourcePath: e.g. '/n/scratch2/xl198/data'
        :return: file's full path
        '''
        for dirpath,dirname,filenames in os.walk(sourcePath):
            for filename in filenames:
                yield os.path.join(dirpath,filename)

    def get_source_dest_file_path(self,sourcePath,destPath):
        '''
        return source,dest path pair for all file under sourcePath
        :param sourcePath:
        :param destPath:
        :return:
        '''
        for sourceFilePath in self.get_source_file_path(sourcePath):
            destFilePath = re.sub(sourcePath,destPath,sourceFilePath)
            self.create_parent_folder(destFilePath)
            yield sourceFilePath,destFilePath

    def remove_hydrogen(self,sourcePath,destPath,offset=None,interval=1000):
        '''
        for all the file under sourcePath, remove hydrogen and then convert
        :param sourcePath:
        :param destPath:
        :param offset:
        :param interval:
        :return:
        '''
        count = 0
        for sourceFilePath,destFilePath in self.get_source_dest_file_path(sourcePath,destPath):
            count +=1
            cmd = 'obabel -ipdb {} -opdb -O {} -d'.format(sourceFilePath,destFilePath)
            if offset == None or (offset%interval == count%interval):
                if not os.path.exists(destFilePath):
                    #print cmd
                    os.popen(cmd)

    def convert_single_file_from_csv(self,csv_entry):
        PDB = csv_entry['PDBname']
        RS  = csv_entry['PDBResId']
        RES,Id = RS.split('_')

        sourceFilePath = os.path.join(self.pdbSourceBase,PDB,'_'.join([PDB,RES,'ligand','fast.pdb']))

        destFolder = os.path.join(self.destBaseFolder,'docked_ligands')

        create_folder(destFolder)

        destFilePath = os.path.join(destFolder,PDB,'_'.join([PDB,RES,'ligand','fast',Id+'.pdb']))

        create_parent_folder(destFilePath)

        cmd = 'obabel -ipdb %s -f %s -l %s -opdb -O %s'%(sourceFilePath,Id,Id,destFilePath)

        os.popen(cmd)



    def get_crystal(self,Id,destFolder):
        # get Crystal ligand based on ID: 1a2b_1234
        pass




