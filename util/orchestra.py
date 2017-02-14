import os,sys
import re
import numpy as np
import pandas as pd
import getopt

import threading
import multiprocessing

from createfolder import try_create_chain_folder,try_create_chain_parent_folder


class orchestra_job:


    def get(self, dataframe, begin, end=None):
        '''

        :param dataframe: pandas.DataFrame or list
        :param begin:
        :return: return index-th value from dataframe
        '''
        if isinstance(dataframe,pd.DataFrame):
            if end:
                return dataframe.iloc[begin:end]
            else:
                return dataframe.iloc[begin]
        elif isinstance(dataframe,list):
            if end:
                return dataframe[begin:end]
            else:
                return dataframe[begin]
        else:
            message = "doesn't support dataframe type ",type(dataframe)
            raise Exception(message)

    def thread_job(self,func,dataframe,index):
        '''
        run func for each item in dataframe

        :param func:
        :param dataframe: pandas.DataFrame or list
        :param index: list[int]
        :return:
        '''
        for i in index:
            func(self.get(dataframe,i))

    def process_job(self,func,dataframe,index):
        '''

        :param func:
        :param dataframe: pandas.DataFrame or list
        :param index: list[int]
        :return:
        '''

        if len(index) < self.thread_num:
            # if job num not enough, running job in process
            for i in index:
                func(self.get(dataframe,i))
        else:
            # create thread to run job
            edge = np.linspace(0, len(index) - 1, self.thread_num + 1).astype(int)

            thread_list = [threading.Thread(target=self.thread_job,
                                            args=(func,
                                                  dataframe,
                                                  range(index[edge[i]],
                                                        index[edge[i + 1]])))
                           for i in range(self.thread_num)]

            for t in thread_list:
                # print "thread start: ",t
                t.start()

            for t in thread_list:
                # print "thread end:",t
                t.join()

    def convert(self,dataframe,convert_func):
        '''

        :param dataframe:
        :param convert_func:
        :return:
        '''

        if hasattr(self, 'cores'):
            self.process_num = self.cores

        # if get a str, read csv
        if type(dataframe) == str:
            try:
                dataframe = pd.read_csv(os.path.join(self.formsFolder, dataframe))
            except Exception as e:
                print e

        if self.arrayjob:
            linspace = np.linspace(0, len(dataframe), self.jobsize + 1).astype(int)
            dataframe = self.get(dataframe,linspace[self.jobid - 1],linspace[self.jobid])

        if len(dataframe) < self.process_num * self.thread_num:
            for i in range(len(dataframe)):
                convert_func(self.get(dataframe,i))
        else:
            edge = np.linspace(0, len(dataframe), self.process_num + 1).astype(int)
            process_list = [multiprocessing.Process(target=self.process_job,
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

    def parse(self):
        try:
            opts, args = getopt.getopt(sys.argv[1:], None, ["jobsize=", "jobid=", "cores="])
        except getopt.GetoptError as err:
            # print help information and exit:
            print str(err)  # will print something like "option -a not recognized"

            sys.exit(2)

        for name, value in opts:
            if name == '--jobsize':
                self.jobsize = int(value)
                print "--jobsize ", value
            if name == '--jobid':
                self.jobid = int(value)
                print "--jobid", value
            if name == '--cores':
                self.cores = int(value)

        if hasattr(self, "jobsize") and hasattr(self, "jobid"):
            self.arrayjob = True

        print "orchestra job ", self.arrayjob

        if hasattr(self, 'cores'):
            print "cores num ", self.cores