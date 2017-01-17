import os,sys
import pandas as pd
import re
import random

sys.path.append("..")
from util.createfolder import create_chain_parent_folder,create_chain_folder


class select:
    dockedBasePath = '/n/scratch2/xl198/dude/data/dude/docked_dude/docked_ligands/'
    workplace = ''
    thread_num = 16
    process_num = 12

    def __init__(self,folderName):
        self.basePath = os.path.join(self.workplace,folderName)
        self.dataframeFolder = os.path.join(self.basePath,'dataframe')

        create_chain_folder(self.dataframeFolder)

    def write_dataframe(self,dataframe,filename):
        '''
        Write dataframe to tempFolderPath
        :param dataframe: pandas DataFrame
        :param filename:  name of output file without suffix
        :return:
        '''
        dataframe.to_csv(os.path.join(self.dataframeFolder,filename+'.csv'),index=False)

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

    def random_select(self,pool, num):
        if num<len(pool):
            raise Exception("Pool size {} smaller than chose num {}".format(len(pool),num))
        index = range(len(pool))
        random.shuffle(index)
        result = map(lambda i:pool[i],index[:num])
        return result

    def select_file(self):


        receptors = self.random_select(os.listdir(self.dockedBasePath), 20)

        actives_col   = []
        decoys_col    = []
        receptors_col = map(lambda receptor:os.path.join('receptors',receptor+'.pdb'),receptors)
        receptors_col = zip(receptors_col,receptors_col)

        for receptor in receptors:
            folder = os.path.join(self.dockedBasePath, receptor)
            actives = [ active for active in os.listdir(folder) if re.search("actives",active) ]
            decoys  = [ decoy  for decoy  in os.listdir(folder) if re.search("decoys", decoy ) ]

            selected_decoys = self.random_select(decoys,len(actives))

            activesSource = map(lambda active:os.path.join(receptor,active),actives)
            activesDest   = map(lambda active:os.path.join('actives',receptor,active),actives)
            decoysSource  = map(lambda decoy:os.path.join(receptor,decoy),selected_decoys)
            decoysDest    = map(lambda decoy:os.path.join('decoys',receptor,decoy),decoys)

            actives_col +=zip(activesSource,activesDest)
            decoys_col  +=zip(decoysSource,decoysDest)


        actives_df  = pd.DataFrame(actives_col,columns=['sourcePath','destPath'])
        decoys_df   = pd.DataFrame(decoys_col ,columns=['sourcePath','destPath'])
        receptors_df= pd.DataFrame(receptors_col, columns=['sourcePath','destPath'])

        self.write_dataframe(actives_df,'actives')
        self.write_dataframe(decoys_df,'decoys')
        self.write_dataframe(receptors_df,'receptors')





