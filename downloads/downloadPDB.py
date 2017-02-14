import os,sys
import prody
import re
import getopt
import pandas as pd
import numpy  as np
import multiprocessing
import threading
from glob import glob

sys.path.append(os.path.dirname(sys.path[0]))
from util.createfolder import try_create_chain_parent_folder,try_create_chain_folder


class parseRCSB:

    def __init__(self):
        '''
        Create folder for kaggle dataset
        tempFolderPath used to store intermediate csv file
        kaggleBasePath used to store converted and coded data
        :param folerName:
        '''
        self.thread_num = FLAGS.thread_num
        self.process_num = FLAGS.process_num
        self.get_address = lambda PDB: 'https://files.rcsb.org/download/' + PDB + '.pdb'


    def error_log(self,content):
        with open(FLAGS.log_file,'a') as fout:
            fout.write(content)

    def repair(self,item):
        destfile = item.replace(FLAGS.splited_receptor_folder,FLAGS.repaired_receptor_folder)

        print destfile
        cmd = 'obminimize -cg -ff MMFF94 -h -n 500 -o pdb  {} > {}'.format(item,destfile)
        os.system(cmd)

    def docking(self,item):
        

    def downloads(self,item):
        download_address = self.get_address(item)
        if os.path.exists(os.path.join(FLAGS.rowdata_folder,item+'.pdb')):
            print item," exists"
            return None
        print 'download ',item
        os.system('wget -P {}  {}'.format(FLAGS.rowdata_folder,download_address))

        pdbname = item.lower()
        ligand_folder = os.path.join(FLAGS.splited_ligand_folder,pdbname)
        try_create_chain_folder(ligand_folder)

        try:
            parsed = prody.parsePDB(os.path.join(FLAGS.rowdata_folder,item+'.pdb'))
        except:
            self.error_log('can not parse {}.\n'.format(item))
            return None
        
        hetero = parsed.select('(hetero and not water) or resname ATP or resname ADP')
        receptor = parsed.select('protein or nucleic')
        if receptor is None:
            self.error_log("{} doesn't have receptor.\n".format(item))
            return None
        if hetero is None:
            self.error_log("{} doesn't have ligand.\n".format(item))
            return None
        ligand_flags = False
        for each in prody.HierView(hetero).iterResidues():
            if each.numAtoms() <= 10:
                continue
            else:
                ligand_flags = True
                ResId = each.getResindex()
                ligand_path = os.path.join(FLAGS.splited_ligand_folder,pdbname,"{}_{}_ligand.pdb".format(pdbname,ResId))
                try_create_chain_parent_folder(ligand_path)
                prody.writePDB(ligand_path,each)

        if ligand_flags:
            receptor_path = os.path.join(FLAGS.splited_receptor_folder,pdbname+'.pdb')
            prody.writePDB(receptor_path,receptor)
        else:
            self.error_log("{} doesn't convert, not ligand have more than 10 atoms.\n".format(item))

    def entry_convert(self, item, coded):
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

        cmd = 'obabel -ipdb {} -f {} -l {} -opdb -O {}'.format(source_file_path, Id, Id, dest_file_path)


        try_create_chain_parent_folder(dest_file_path)

        if not os.path.exists(dest_file_path):
            os.popen(cmd)

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
            func(dataframe[i])

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

    def convert(self):
        '''
        according the result of 'database_from_csv'
        running multiprocess to get result

        :param dataframe: pandas.DataFrame
                          string
        :param coded: bool
        :return:
        '''

        convert_func = self.downloads



        if hasattr(FLAGS, 'cores'):
            self.process_num = FLAGS.cores



        if hasattr(FLAGS,'arrayjob') and FLAGS.arrayjob:
            jobsize = FLAGS.jobsize
            jobid = FLAGS.jobid

            # using iloc to slice dataframe result doesn't contains end [start,end)
            linspace = np.linspace(0, len(FLAGS.pdb_list), jobsize + 1).astype(int)
            dataframe = FLAGS.pdb_list[linspace[jobid - 1]:linspace[jobid]]
        else:
            dataframe = FLAGS.pdb_list
        print "dataframe size ", len(dataframe)

        # when there's not enough entry to comvert , decrease thread's num
        if len(dataframe) < self.process_num * self.thread_num:
            for i in range(len(dataframe)):
                convert_func(dataframe[i])
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






class FLAGS:

    workplace = '/n/scratch2/xl198/data/rcsb'
    rowdata_folder = os.path.join(workplace,'row')
    splited_receptor_folder = os.path.join(workplace,'row_receptor')
    splited_ligand_folder = os.path.join(workplace,'ligands')
    repaired_receptor_folder = os.path.join(workplace,'repaired_receptor')
    log_file = 'error.log'
    thread_num = 1
    process_num = 1






def parse_FLAG():
    try:
        opts,args = getopt.getopt(sys.argv[1:],None,["offset=","jobsize=","jobid=","cores=","scan"])
    except getopt.GetoptError as err:
        # print help information and exit:
        print str(err)  # will print something like "option -a not recognized"

        sys.exit(2)

    for name,value in opts:
        if name == '--offset':
            FLAGS.joboffset = int(value)
        if name == '--jobsize':
            FLAGS.jobsize = int(value)
            print "--jobsize ",value
        if name == '--jobid':
            FLAGS.jobid = int(value)
            print "--jobid", value
        if name == '--cores':
            FLAGS.cores = int(value)
        if name == '--scan':
            FLAGS.scan = True

    if hasattr(FLAGS,"jobsize") and hasattr(FLAGS,"jobid"):
        FLAGS.arrayjob = True

    if hasattr(FLAGS,'arrayjob'):
        print "orchestra job ",FLAGS.arrayjob

    if hasattr(FLAGS,'cores'):
        print "cores num ",FLAGS.cores


    #content = open('/home/xl198/code/dock/downloads/target_PDB.txt').readline()
    #content = content.split(',')
    #content = map(lambda x:x.strip(),content)
    #FLAGS.pdb_list= content
    FLAGS.pdb_list = glob(os.path.join(FLAGS.splited_receptor_folder,"*.pdb"))


if __name__ == '__main__':
    parse_FLAG()
    parser = parseRCSB()
    parser.convert()
