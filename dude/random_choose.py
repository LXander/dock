import os,sys
import pandas as pd
import re
import random
import threading
import multiprocessing
import numpy as np
import getopt

test_receptor = ["aa2ar.pdb",
"ada17.pdb",
"adrb1.pdb",
"cp2c9.pdb",
"cp3a4.pdb",
"dyr.pdb",
"gcr.pdb",
"gria2.pdb",
"hdac2.pdb",
"kith.pdb",
"lkha4.pdb",
"mk10.pdb",
"mmp13.pdb",
"nram.pdb",
"pgh1.pdb",
"ppara.pdb",
"ppard.pdb",
"pparg.pdb",
"rxra.pdb",
"vgfr2.pdb"]

sys.path.append(os.path.dirname(sys.path[0])) 
from util.createfolder import create_chain_parent_folder,create_chain_folder,try_create_chain_folder


class select:
    dockedBasePath = '/n/scratch2/xl198/dude/data/dude/docked_dude/docked_ligands/'
    workplace = '/n/scratch2/xl198/dude/data'
    thread_num = 16
    process_num = 12

    def __init__(self,folderName):
        self.basePath = os.path.join(self.workplace,folderName)
        self.dataframeFolder = os.path.join(self.basePath,'dataframe')
        self.datasetFolder = os.path.join(self.basePath,'dataset')

        create_chain_folder(self.dataframeFolder)

    def write_dataframe(self,dataframe,filename):
        '''
        Write dataframe to tempFolderPath
        :param dataframe: pandas DataFrame
        :param filename:  name of output file without suffix
        :return:
        '''
        dataframe.to_csv(os.path.join(self.dataframeFolder,filename+'.csv'),index=False)

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

    def process_convert(self,func,dataframe,index):

        # linspace contain end value but range don't
        # so we use edge[i+1] to select value in index
        # end should be len(index)-1
        edge = np.linspace(0,len(index)-1,self.thread_num+1).astype(int)

        thread_list = [ threading.Thread(target=self.thread_convert,
                                         args=(func,
                                               dataframe,
                                               range(index[edge[i]],
                                                     index[edge[i+1]])))
                        for i in range(self.thread_num) ]

        for t in thread_list:
            #print "thread start: ",t
            t.start()

        for t in thread_list:
            t.join()

    def random_select(self,pool, num):
        if num>len(pool):
            raise Exception("Pool size {} smaller than chose num {}".format(len(pool),num))
        index = range(len(pool))
        random.shuffle(index)
        result = map(lambda i:pool[i],index[:num])
        return result

    def select_except_test(self,pool,test_pool):
        return list(set(pool)-set(test_pool))

    def select_file(self):


        #receptors = self.random_select(os.listdir(self.dockedBasePath), 20)
        receptors = self.select_except_test(os.listdir(self.dockedBasePath),test_receptor)
        actives_col   = []
        decoys_col    = []
        receptors_col = map(lambda receptor:os.path.join('receptors',receptor+'.pdb'),receptors)
        receptors_col = zip(receptors_col,receptors_col)

        for receptor in receptors:
            folder = os.path.join(self.dockedBasePath, receptor)
            actives = [ active for active in os.listdir(folder) if re.search("actives",active) ]
            decoys  = [ decoy  for decoy  in os.listdir(folder) if re.search("decoys", decoy ) ]
            if len(actives) == 0 or len(decoys) == 0:
                print "receptor {} doesn't have actives or decoys ligands".format(receptor)
                continue
            selected_decoys = self.random_select(decoys,len(actives))

            activesSource = map(lambda active:os.path.join(receptor,active),actives)
            activesDest   = map(lambda active:[active,os.path.join('actives',active)],activesSource)

            decoysSource  = map(lambda decoy:os.path.join(receptor,decoy),selected_decoys)
            decoysDest    = map(lambda decoy:[decoy,os.path.join('decoys',decoy)],decoysSource)

            actives_col +=activesDest
            decoys_col  +=decoysDest

        actives_df  = pd.DataFrame(actives_col,columns=['sourcePath','destPath'])
        decoys_df   = pd.DataFrame(decoys_col ,columns=['sourcePath','destPath'])
        receptors_df= pd.DataFrame(receptors_col, columns=['sourcePath','destPath'])

        self.write_dataframe(actives_df,'actives')
        self.write_dataframe(decoys_df,'decoys')
        self.write_dataframe(receptors_df,'receptors')

    def obabel_split(self,item,limit=400):
        source_file_path = os.path.join(self.dockedBasePath,item['sourcePath'])
        dest_file_path = os.path.join(self.datasetFolder,item['destPath'])

        dest_folder = dest_file_path.split('.')[0]
        dest_file = os.path.basename(dest_file_path)


        create_chain_folder(dest_folder)

        final_path = os.path.join(dest_folder,dest_file)

        file_base_name,file_suffix = final_path.split('.')

        for i in range(1,limit+1):
	    final_file_path = file_base_name+'_'+str(i)+'.'+file_suffix
            cmd = 'obabel -ipdb {} -f {} -l {} -opdb -O {}'.format(source_file_path, i, i, final_file_path)
            print i
            os.popen(cmd)
            content = os.popen("cat {} | wc -l".format(final_file_path)).read().strip()
            
            if content == '0':
                os.popen("rm {}".format(final_file_path))
                break

    def copy_file(self,item):
        '''
        smina only can use original receptor file
        prody need receptor which is convert by obabel
        '''
        sourceBasePath = '/n/scratch2/xl198/dude/data/dude/pdbs/receptors/'
        source_file_path = os.path.join(sourceBasePath,os.path.basename(item['sourcePath']))
        dest_file_path = os.path.join(self.datasetFolder,item['destPath'])

        create_chain_parent_folder(dest_file_path)

        cmd = 'cp {} {}'.format(source_file_path,dest_file_path)
        print cmd
        os.popen(cmd)

    def copy_ligand(self,item):
        '''
        copy original ligand directly to the /original folder
        :param item:  pd.DataFrame item
        :return:
        '''
        sourceBasePath ='/n/scratch2/xl198/dude/data/dude/pdbs/ligands'
        source_file_path = os.path.join(sourceBasePath,item['sourcePath'])
        dest_file_path = os.path.join(self.datasetFolder,'original',item['destPath'])
        
        try:
            create_chain_parent_folder(dest_file_path)
        except:
            pass
        cmd = "cp {} {}".format(source_file_path,dest_file_path)
        #print cmd
        os.popen(cmd)


    def convert(self,dataframe,is_docked=True,func = None):
        '''
        according the result of 'database_from_csv'
        running multiprocess to get result

        :param dataframe: pandas.DataFrame
                          string
        :param coded: bool
        :return:
        '''


        convert_func = self.obabel_split if is_docked else self.copy_file

        if func:
            convert_func = func

        # if get a str, read csv
        if type(dataframe) == str:
            try:
                dataframe = pd.read_csv(os.path.join(self.dataframeFolder,dataframe))
            except Exception as e:
                print e

        if not is_docked:
            for i in range(len(dataframe)):
                convert_func(dataframe.iloc[i]) 
           
        elif FLAGS.orchestra_arrayjob:
            jobsize = FLAGS.orchestra_jobsize
            jobid = FLAGS.orchestra_jobid
            linspace = np.linspace(0, len(dataframe), jobsize + 1).astype(int)
            #dataframe = dataframe.iloc[linspace[jobid - 1]:linspace[jobid]]
            index = range(linspace[jobid-1],linspace[jobid])
            for i in index:
                convert_func(dataframe.iloc[i])

        else:


            edge = np.linspace(0,len(dataframe),self.process_num+1).astype(int)
            process_list = [ multiprocessing.Process(target=self.process_convert,
                                                     args=(convert_func,
                                                           dataframe,
                                                           range(edge[i],
                                                                 edge[i+1])))
                             for i in range(self.process_num)]

            for p in process_list:
                #print "process start: ",p
                p.start()

            for p in process_list:
                p.join()


class FLAGS:
    orchestra_arrayjob = False

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



if __name__ == '__main__':
    parse_FLAG()
    sel = select('jan_22')
    #sel.select_file()
    #sel.convert('actives.csv',is_docked=True)
    #sel.convert('decoys.csv',is_docked=True)
    sel.convert('receptors.csv',is_docked=False)
    #sel.convert('decoys.csv',func=sel.copy_ligand)
