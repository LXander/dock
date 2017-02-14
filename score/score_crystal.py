import os,sys
import re
import pandas as pd
from glob import glob
import getopt
sys.path.append(os.path.dirname(sys.path[0]))
from util.orchestra import orchestra_job
from util.createfolder import try_create_chain_parent_folder

smina = '/home/xl198/program/smina/smina.static'

class Score(orchestra_job):
    arrayjob = False
    workplace = '/n/scratch2/xl198/data/score'
    formsPath = os.path.join(workplace,'forms')
    dockSourcePath = '/n/scratch2/xl198/dude/data/dude/pdbs/ligands'
    smina = '/home/xl198/program/smina/smina.static'
    thread_num = 1
    process_num = 1


    def __init__(self):
        self.parse()

        try_create_chain_parent_folder(self.workplace)
        try_create_chain_parent_folder(self.formsPath)

    def scoring(self,crystal_path):
        file_name = os.path.basename(crystal_path)
        ligand_id = '_'.join(file_name.split('_')[:2])
        receptor = file_name.split('_')[0]
        receptor_path = os.path.dirname(crystal_path).replace('addH','data')
        receptor_file_path = os.path.join(receptor_path,receptor+'.pdb')

        cmd = '{} -r {} -l {} --score_only'.format(self.smina,receptor_file_path,crystal_path)
        print cmd
        result = os.popen(cmd)
        smina_score = re.search('Affinity: (-?\d+.\d+) \(kcal/mol\)',result).groups()[0]

        with open(os.path.join(self.formsPath,'score_'+str(self.jobid)+'.txt'),'a') as fout:
            fout.write(ligand_id + ',' + smina_score+'\n')

    def get_score(self):
        crystal_list = glob(os.path.join('/n/scratch2/xl198/data/H/addH', '*', '*.pdb'))
        self.convert(crystal_list,self.scoring)

if __name__ == '__main__':
    score = Score()
    score.get_score()