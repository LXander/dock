import os,sys
import re
import shutil

class obabelOperation:
    crystalSourceBase = '/n/scratch2/xl198/data/H/addH'
    receptorSourceBase = '/n/scratch2/xl198/data/H/data'
    pdbSourceBase ='/n/scratch2/xl198/data/pdbs'
    destBase = '/n/scratch2/xl198/data/'



    def __init__(self,baseFolder):
        self.destBaseFolder = os.path.join(self.destBase,baseFolder,'pdbs')
        self.create_folder(self.destBaseFolder)

    def get_source_file_path(self,sourcePath):
        for dirpath,dirname,filenames in os.walk(sourcePath):
            for filename in filenames:
                yield os.path.join(dirpath,filename)

    def get_source_dest_file_path(self,sourcePath,destPath):
        for sourceFilePath in self.get_source_file_path(sourcePath):
            destFilePath = re.sub(sourcePath,destPath,sourceFilePath)
            self.create_parent_folder(destFilePath)
            yield sourceFilePath,destFilePath

    def create_parent_folder(self,filePath):
        folderPath = os.path.dirname(filePath)
        if not os.path.exists(folderPath):
            os.mkdir(folderPath)

    def create_folder(self,folderPath):
        if not os.path.exists(folderPath):
            os.mkdir(folderPath)

    def remove_hydrogen(self,sourcePath,destPath,offset=None,interval=1000):
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

        self.create_folder(destFolder)

        destFilePath = os.path.join(destFolder,PDB,'_'.join([PDB,RES,'ligand','fast',Id+'.pdb']))

        self.create_parent_folder(destFilePath)

        cmd = 'obabel -ipdb %s -f %s -l %s -opdb -O %s'%(sourceFilePath,Id,Id,destFilePath)

        os.popen(cmd)

    def get_crystal(self,Id):
        # get Crystal ligand based on ID: 1a2b_1234
        pass




