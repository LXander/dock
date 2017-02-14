import os,sys
import pandas
from glob import glob

smina = '/home/xl198/program/smina/smina.static'



def scoring(crystal_path):
    file_name = os.path.basename(crystal_path)
    receptor = file_name.split('_')[0]
    receptor_path = os.path.dirname(crystal_path).replace('addH','data')
    receptor_file_path = os.path.join(receptor_path,receptor+'.pdb')
    cmd = '{} -r {} -l {} --score_only'
    result = os.popen(cmd)
    print cmd

def get_score():
    crystal_list = glob(os.path.join('/n/scratch2/xl198/data/H/addH','*','*.pdb'))
    for crystal in crystal_list:
        scoring(crystal)

if __name__ == '__main__':
    get_score()