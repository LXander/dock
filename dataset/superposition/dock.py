import os,sys
import re
import numpy as np
import pandas as pd
import getopt

import threading
import multiprocessing

sys.path.append(os.path.dirname(os.path.dirname(sys.path[0])))
from util.createfolder import try_create_chain_folder


class docking:

    def __init__(self,folderName):
        self.basePath = os.path.join(FLAGS.workplace,folderName)
        self.tempFolderPath = os.path.join(self.basePath,'temp')
        self.dockingFolderPath = os.path.join(self.basePath,'dock')
        self.thread_num = FLAGS.thread_num
        self.process_num = FLAGS.process_num

        try_create_chain_folder(self.tempFolderPath)
        try_create_chain_folder(self.dockingFolderPath)

    def write_dataframe(self,dataframe,filename):
        '''
        Write dataframe to tempFolderPath
        :param dataframe: pandas DataFrame
        :param filename:  name of output file without suffix
        :return:
        '''
        dataframe.to_csv(os.path.join(self.tempFolderPath,filename+'.csv'),index=False)

    def thread_convert(self,func,dataframe,index):
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

    def convert(self, dataframe,convert_func):
        '''
        according the result of 'database_from_csv'
        running multiprocess to get result

        :param dataframe: pandas.DataFrame
                          string
        :param coded: bool
        :return:
        '''

        #convert_func = self.entry_convert if is_docked else self.receptor_crystal_convert



        if hasattr(FLAGS, 'cores'):
            self.process_num = FLAGS.cores

        # if get a str, read csv
        if type(dataframe) == str:
            try:
                dataframe = pd.read_csv(os.path.join(self.tempFolderPath, dataframe))
            except Exception as e:
                print e

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

    def convert_function(self,item):
        dockDestPath = os.path.join(self.dockingFolderPath)
        receptor = item['receptor']
        ligand_box = item['ligand_box']
        ligand = item['ligand']
        dest_ligand = ligand.replace(FLAGS.dockSourcePath,dockDestPath)
        try_create_chain_folder(dest_ligand)
        cmd = "{} -r {} -l {} --autobox_ligand {} -o {} --num_modes=1000 --energy_range=100 --cpu=1 ".format(FLAGS.smina,receptor,ligand,ligand_box,dest_ligand)
        print cmd

class FLAGS:
    arrayjob = False
    workplace = '/n/scratch2/xl198/data'
    dockSourcePath = '/n/scratch2/xl198/dude/data/dude/pdbs/ligands/'
    smina = '/home/xl198/program/smina/smina.static'
    thread_num = 16
    process_num = 12

def parse_FLAG():
    try:
        opts,args = getopt.getopt(sys.argv[1:],None,["jobsize=","jobid=","cores="])
    except getopt.GetoptError as err:
        # print help information and exit:
        print str(err)  # will print something like "option -a not recognized"

        sys.exit(2)

    for name,value in opts:
        if name == '--jobsize':
            FLAGS.jobsize = int(value)
            print "--jobsize ",value
        if name == '--jobid':
            FLAGS.jobid = int(value)
            print "--jobid", value
        if name == '--cores':
            FLAGS.cores = int(value)

    if hasattr(FLAGS,"orchestra_jobsize") and hasattr(FLAGS,"orchestra_jobid"):
        FLAGS.arrayjob = True

    print "orchestra job ",FLAGS.orchestra_arrayjob

    if hasattr(FLAGS,'cores'):
        print "cores num ",FLAGS.cores


if __name__ == '__main__':
    parse_FLAG()
    dockclass = docking('Superposition')
    dockclass.convert('/n/scratch2/xl198/data/Superposition/forms/docking_pair.csv',dockclass.convert_function)