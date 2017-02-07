import os, sys
import re
import numpy as np
import pandas as pd
import getopt

import threading
import multiprocessing

sys.path.append(os.path.dirname(os.path.dirname(sys.path[0])))
from util.createfolder import try_create_chain_folder, try_create_chain_parent_folder


class docking:
    def __init__(self, folderName):
        self.basePath = os.path.join(FLAGS.workplace, folderName)
        self.tempFolderPath = os.path.join(self.basePath, 'temp')
        self.dockingFolderPath = os.path.join(self.basePath, 'dock')
        self.thread_num = FLAGS.thread_num
        self.process_num = FLAGS.process_num

        try_create_chain_folder(self.tempFolderPath)
        try_create_chain_folder(self.dockingFolderPath)

    def write_dataframe(self, dataframe, filename):
        '''
        Write dataframe to tempFolderPath
        :param dataframe: pandas DataFrame
        :param filename:  name of output file without suffix
        :return:
        '''
        dataframe.to_csv(os.path.join(self.tempFolderPath, filename + '.csv'), index=False)

    def thread_convert(self, func, dataframe, index):
        '''
        job running in thread
        :param func: function running in single thread
        :param dataframe: pandas.DataFrame
        :param coded: bool
        :param index: list of int
        :return:
        '''
        for i in index:
            # print "dataframe size: {}, ix {}".format(len(dataframe),i)
            func(dataframe.iloc[i])

    def process_convert(self, func, dataframe, index):

        # linspace contain end value but range don't
        # so we use edge[i+1] to select value in index
        # end should be len(index)-1


        if len(index) < self.thread_num:
            for i in index:
                func(dataframe.iloc[i])
            return

        edge = np.linspace(0, len(index) - 1, self.thread_num + 1).astype(int)

        thread_list = [threading.Thread(target=self.thread_convert,
                                        args=(func,
                                              dataframe,
                                              range(index[edge[i]],
                                                    index[edge[i + 1]])))
                       for i in range(self.thread_num)]

        for t in thread_list:
            # print "thread start: ",t
            t.start()

        for t in thread_list:
            t.join()

    def convert(self, dataframe, convert_func):
        '''
        according the result of 'database_from_csv'
        running multiprocess to get result

        :param dataframe: pandas.DataFrame
                          string
        :param coded: bool
        :return:
        '''

        # convert_func = self.entry_convert if is_docked else self.receptor_crystal_convert



        if hasattr(FLAGS, 'cores'):
            self.process_num = FLAGS.cores

        # if get a str, read csv
        if type(dataframe) == str:
            try:
                dataframe = pd.read_csv(os.path.join(self.tempFolderPath, dataframe))
            except Exception as e:
                print e

        if hasattr(FLAGS, 'joboffset'):
            loc = FLAGS.joboffset * FLAGS.jobsize + FLAGS.jobid - 1
            print "convert offset single ", loc
            try:
                convert_func(dataframe.ix[loc])
            except:
                pass
            return None

        if hasattr(FLAGS, 'scan'):
            print "scan length ", len(dataframe)
            for i in range(len(dataframe)):
                print "convert single ", i

                convert_func(dataframe.ix[i])

            return None

        if FLAGS.arrayjob:
            jobsize = FLAGS.jobsize
            jobid = FLAGS.jobid

            # using iloc to slice dataframe result doesn't contains end [start,end)
            linspace = np.linspace(0, len(dataframe), jobsize + 1).astype(int)
            dataframe = dataframe.iloc[linspace[jobid - 1]:linspace[jobid]]

        print "dataframe size ", len(dataframe)

        # when there's not enough entry to comvert , decrease thread's num
        if len(dataframe) < self.process_num * self.thread_num:

            for i in range(len(dataframe)):
                convert_func(dataframe.iloc[i], coded)
            return
        edge = np.linspace(0, len(dataframe), self.process_num + 1).astype(int)
        process_list = [multiprocessing.Process(target=self.process_convert,
                                                args=(convert_func,
                                                      dataframe,
                                                      range(edge[i],
                                                            edge[i + 1])))
                        for i in range(self.process_num)]

        for p in process_list:
            print "process start: ", p
            p.start()

        for p in process_list:
            print "process end: ", p
            p.join()

    def convert_ritor(self,item):
        destPath = os.path.join(self.dockingFolderPath, 'rigor', 'docked_ligands')

        sourceFileName = item['Id'] + '_ligand_rigor.pdb'

        destFileName = item['Id'] + '_ligand_rigor_.pdb'

        sourceFilePath = os.path.join(FLAGS.rigor_sourcePath, item['PDBname'], sourceFileName)

        destFilePath = os.path.join(destPath, item['PDBname'], destFileName)

        try_create_chain_parent_folder(destFilePath)

        cmd = 'obabel -ipdb {} -opdb -O {} -m'.format(sourceFilePath, destFilePath)

        receptor_folder = os.path.join(self.dockingFolderPath, 'rigor', 'receptors')
        receptor_source = os.path.join(FLAGS.receptor_source, item['PDBname'], item['PDBname'] + '.pdb')
        receptor_dest = os.path.join(receptor_folder, item['PDBname'], item['PDBname'] + '.pdb')

        if not os.path.exists(receptor_dest):
            try_create_chain_parent_folder(receptor_dest)
            os.system('cp {} {}'.format(receptor_source, receptor_dest))

        crystal_folder = os.path.join(self.dockingFolderPath, 'rigor', 'crystal_ligands')
        crystal_source = os.path.join(FLAGS.crystal_source, item['PDBname'], item['Id'] + '_ligand.pdb')
        crystal_dest = os.path.join(crystal_folder, item['PDBname'], item['Id'] + '_ligand.pdb')

        if not os.path.exists(crystal_dest):
            try_create_chain_parent_folder(crystal_dest)
            os.system('cp {} {}'.format(crystal_source, crystal_dest))

        # print cmd
        #os.system(cmd)

    def convert_ritor_so(self,item):
        destPath = os.path.join(self.dockingFolderPath, 'rigor_so', 'docked_ligands')

        sourceFileName = item['Id'] + '_ligand_rigor_so.pdb'

        destFileName = item['Id'] + '_ligand_rigor_so_.pdb'

        sourceFilePath = os.path.join(FLAGS.rigor_so_sourcePath, item['PDBname'], sourceFileName)

        destFilePath = os.path.join(destPath, item['PDBname'], destFileName)

        try_create_chain_parent_folder(destFilePath)

        cmd = 'obabel -ipdb {} -opdb -O {} -m'.format(sourceFilePath, destFilePath)

        receptor_folder = os.path.join(self.dockingFolderPath, 'rigor_so', 'receptors')
        receptor_source = os.path.join(FLAGS.receptor_source, item['PDBname'], item['PDBname'] + '.pdb')
        receptor_dest = os.path.join(receptor_folder, item['PDBname'], item['PDBname'] + '.pdb')

        if not os.path.exists(receptor_dest):
            try_create_chain_parent_folder(receptor_dest)
            os.system('cp {} {}'.format(receptor_source, receptor_dest))

        crystal_folder = os.path.join(self.dockingFolderPath, 'rigor_so', 'crystal_ligands')
        crystal_source = os.path.join(FLAGS.crystal_source, item['PDBname'], item['Id'] + '_ligand.pdb')
        crystal_dest = os.path.join(crystal_folder, item['PDBname'], item['Id'] + '_ligand.pdb')

        if not os.path.exists(crystal_dest):
            try_create_chain_parent_folder(crystal_dest)
            os.system('cp {} {}'.format(crystal_source, crystal_dest))

        # print cmd
        os.system(cmd)

    def convert_function(self, item):
        destPath =os.path.join(self.dockingFolderPath,'fast','decked_ligands')

        sourceFileName = item['Id']+'_ligand_fast.pdb'

        destFileName = item['Id']+'_ligand_fast_.pdb'

        sourceFilePath = os.path.join(FLAGS.sourcePath,item['PDBname'],sourceFileName)

        destFilePath = os.path.join(destPath,item['PDBname'],destFileName)


        try_create_chain_parent_folder(destFilePath)

        cmd = 'obabel -ipdb {} -opdb -O {} -m'.format(sourceFilePath,destFilePath)

        receptor_folder = os.path.join(self.dockingFolderPath,'fast','receptors')
        receptor_source = os.path.join(FLAGS.receptor_source,item['PDBname'],item['PDBname']+'.pdb')
        receptor_dest = os.path.join(receptor_folder,item['PDBname'],item['PDBname']+'.pdb')

        if not os.path.exists(receptor_dest):
            try_create_chain_parent_folder(receptor_dest)
            os.system('cp {} {}'.format(receptor_source,receptor_dest))


        crystal_folder = os.path.join(self.dockingFolderPath,'fast','crystal_ligands')
        crystal_source = os.path.join(FLAGS.crystal_source,item['PDBname'],item['Id']+'_ligand.pdb')
        crystal_dest = os.path.join(crystal_folder,item['PDBname'],item['Id']+'_ligand.pdb')

        if not os.path.exists(crystal_dest):
            try_create_chain_parent_folder(crystal_dest)
            os.system('cp {} {}'.format(crystal_source,crystal_dest))

        #print cmd
        os.system(cmd)


class FLAGS:
    arrayjob = False
    workplace = '/n/scratch2/xl198/data'
    sourcePath = '/n/scratch2/xl198/data/pdbs'
    crystal_source = '/n/scratch2/xl198/data/H/addH'
    receptor_source ='/n/scratch2/xl198/data/H/data'

    rigor_sourcePath = '/n/scratch2/xl198/YI/rigor/final'
    rigor_so_sourcePath  = '/n/scratch2/xl198/YI/rigor_so/final'
    
    dockSourcePath = '/n/scratch2/xl198/dude/data/dude/pdbs/ligands'
    smina = '/home/xl198/program/smina/smina.static'
    thread_num = 1
    process_num = 1


def parse_FLAG():
    try:
        opts, args = getopt.getopt(sys.argv[1:], None, ["offset=", "jobsize=", "jobid=", "cores=", "scan"])
    except getopt.GetoptError as err:
        # print help information and exit:
        print str(err)  # will print something like "option -a not recognized"

        sys.exit(2)

    for name, value in opts:
        if name == '--offset':
            FLAGS.joboffset = int(value)
        if name == '--jobsize':
            FLAGS.jobsize = int(value)
            print "--jobsize ", value
        if name == '--jobid':
            FLAGS.jobid = int(value)
            print "--jobid", value
        if name == '--cores':
            FLAGS.cores = int(value)
        if name == '--scan':
            FLAGS.scan = True

    if hasattr(FLAGS, "jobsize") and hasattr(FLAGS, "jobid"):
        FLAGS.arrayjob = True

    print "orchestra job ", FLAGS.arrayjob

    if hasattr(FLAGS, 'cores'):
        print "cores num ", FLAGS.cores


if __name__ == '__main__':
    parse_FLAG()
    dockclass = docking('fusion')
    dockclass.convert('/n/scratch2/xl198/YI/rigor/selected.csv', dockclass.convert_ritor_so)
