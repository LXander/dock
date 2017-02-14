import os, sys
import prody
import re
import getopt
import pandas as pd
import numpy  as np
import multiprocessing
import threading
import tempfile
import random

class parseRCSB:
    '''
    Download pdb from rcsb and split it into receptor and ligand
    '''
    def __init__(self):

        self.thread_num = FLAGS.thread_num
        self.process_num = FLAGS.process_num
        self.log_file = FLAGS.log_file
        self.get_dataframe()

    def get_dataframe(self):

        content = open('target_PDB.txt').readline()
        content = content.split(',')
        content = map(lambda x: x.strip(), content)
        self.pdb_list = content

    def error_log(self, content):
        # write down error information
        with open(self.log_file, 'a') as fout:
            fout.write(content)


    def noise_scoring_function(self,ligand_name):

        random_scalar = lambda x:x*2 if random.randint(0,1) else x*0.5

        score_file_path = os.path.join(FLAGS.tempfolder,ligand_name+'.score')

        with open(score_file_path,'w') as fout:
            fout.write('{}    gauss(o=0,_w=0.5,_c=8)\n'.format(random_scalar(-0.035579)))
            fout.write('{}    gauss(o=3,_w=2,_c=8)\n'.format(random_scalar(-0.005156)))
            fout.write('{}     repulsion(o=0,_c=8)'.format(random_scalar(0.840245)))
            fout.write('{}    hydrophobic(g=0.5,_b=1.5,_c=8)'.format(random_scalar(-0.035069)))
            fout.write('{}    non_dir_h_bond(g=-0.7,_b=0,_c=8)'.format(random_scalar(-0.587439)))
            fout.write('{}        num_tors_div'.format(random_scalar(1.923)))

        return score_file_path

    def noise_docing(self,item):

        ligand_name = os.path.basename(item).split('.')[0]

        receptor = ligand_name.split('_')[0]
        receptor_path = os.path.join(FLAGS.receptor_folder,receptor+'.pdb')


        noise_ligand_path = item.replace(FLAGS.source_path,FLAGS.dest_path)
        noise_ligand_path = noise_ligand_path.replace('_ligand','_noise_ligand')

        scoring_file_path = self.noise_scoring_function(ligand_name)



        com = '{} -r {} -l {} --autobox_ligand {} -o {} --custom_scoring {} --num_modes=1000 --energy_range=100 --cpu=1 '\
            .format(FLAGS.smina,receptor_path,item,item,noise_ligand_path,scoring_file_path)






    def downloads(self, item):
        '''
        Download pdb from rcsb and split it into receptor and ligand
        :param item: 4 letter PDB ID '3EML'
        :return:
        '''

        # Download pdb to rowdata_folder
        download_address = 'https://files.rcsb.org/download/' + item + '.pdb'
        os.system('wget -P {} {}'.format(FLAGS.rowdata_folder,download_address))

        # create folder to store ligand
        pdbname = item.lower()
        ligand_folder = os.path.join(FLAGS.splited_ligand_folder, pdbname)
        if not os.path.exists(ligand_folder):
            os.mkdir(ligand_folder)

        # parse pdb
        try:
            parsed = prody.parsePDB(os.path.join(FLAGS.rowdata_folder, item + '.pdb'))
        except:
            self.error_log('can not parse {}.\n'.format(item))
            return None

        # select receptor and ligand
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
            if each.numAtoms() <= FLAGS.atom_num_threahold:
                # ignore ligand if atom num is less than threshold
                continue
            else:
                ligand_flags = True
                ResId = each.getResindex()
                ligand_path = os.path.join(FLAGS.splited_ligand_folder, pdbname,
                                           "{}_{}_ligand.pdb".format(pdbname, ResId))
                if not os.path.exists(os.path.dirname(ligand_path)):
                    os.mkdir(os.path.dirname(ligand_path))
                prody.writePDB(ligand_path, each)

        if ligand_flags:
            receptor_path = os.path.join(FLAGS.splited_receptor_folder, pdbname + '.pdb')
            prody.writePDB(receptor_path, receptor)
        else:
            self.error_log("{} doesn't convert, no ligand have more than 10 atoms.\n")

    def thread_convert(self, func, dataframe, index):
        for i in index:
            func(dataframe[i])

    def process_convert(self, func, dataframe, index):
        # linspace contain end value but range don't
        # so we use edge[i+1] to select value in index
        # end should be len(index)-1

        if len(index) < self.thread_num:
            for i in index:
                func(dataframe[i])
            return

        edge = np.linspace(0, len(index) - 1, self.thread_num + 1).astype(int)
        thread_list = [threading.Thread(target=self.thread_convert,
                                        args=(func,
                                              dataframe,
                                              range(index[edge[i]],
                                                    index[edge[i + 1]])))
                       for i in range(self.thread_num)]

        for t in thread_list:
            t.start()

        for t in thread_list:
            t.join()

    def convert(self):

        convert_func = self.downloads

        # when there's not enough entry to comvert , decrease thread's num
        if len(self.pdb_list) < self.process_num * self.thread_num:
            for i in range(len(self.pdb_list)):
                convert_func(self.pdb_list[i])
            return

        edge = np.linspace(0, len(self.pdb_list), self.process_num + 1).astype(int)
        process_list = [multiprocessing.Process(target=self.process_convert,
                                                args=(convert_func,
                                                      self.pdb_list,
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
    workplace = './rcsb'
    rowdata_folder = os.path.join(workplace, 'row')
    splited_receptor_folder = os.path.join(workplace, 'row_receptor')
    splited_ligand_folder = os.path.join(workplace, 'ligands')

    log_file = 'error.log'
    thread_num = 16
    process_num = 12
    atom_num_threahold = 10

def parse_FLAGS():
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

    print "orchestra job ",FLAGS.arrayjob

    if hasattr(FLAGS,'cores'):
        print "cores num ",FLAGS.cores

    FLAGS.tempfolder = tempfile.mkdtemp()


if __name__ == '__main__':
    parse_FLAGS()
    parser = parseRCSB()
    parser.convert()