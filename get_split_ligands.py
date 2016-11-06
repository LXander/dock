import os,sys
import config
import re

'''
    This file is used to get seperated pdb file from Yi's folder
    mol2 file in organized directory and '@<TRIPOS>MOLECULE' indicate the
    number of ligands in file.

    To submit it with Jobarray, the base is needed cause Jobarray only allowed
    1000 jobId at a time

    1. get the path of input files
    2. create folder
    2. calculate the number of the ligands in files
    3. convert the top 10 and bottom 10
'''

# if only few ligand in the best and worst won't very difference
minimum_ligand_num = 100
# mol2 file organized in input_path
input_path = os.path.join(config.BASE_DATA,'result')
# select and convert pdb file to output_path
output_path = os.path.join(config.BASE_DATA,'select','ligands')


def convert(input_file):
    ligand_count = 0
    with open(input_file) as input:
        for line in input:
            if line == '@<TRIPOS>MOLECULE':
                ligand_count += 1

    dirname = os.path.dirname(input_file)
    receptor = os.path.basename(dirname)
    filename = os.path.basename(input_file).split('.')[0]

    ligand_output_path = os.path.join(output_path,receptor)
    if not os.path.exists(ligand_output_path):
        os.mkdir(ligand_output_path)

    if ligand_count > minimum_ligand_num:
        # convert top 10 ligand
        os.system('obabel -imol2 %s -f %d -l %d -opdb -O %s -d'%
                  (input_file,1,10,os.path.join(ligand_output_path,filename,'_top_.pdb')))
        # convert bottom 10 ligand
        os.system('obabel -imol2 %s -f %d -l %d -opdb -O %s -d' %
                  (input_file, ligand_count-10, ligand_count, os.path.join(ligand_output_path, filename, '_top_.pdb')))



def run(base,offset):

    # JobId start at 1
    cur = base * 1000 + offset -1



    file_path = None

    for dirname, dirnames, filenames in os.walk(input_path):
        if len(filenames) > cur:
            if cur>=0:
                file_path = os.path.join(dirname,filenames[cur])
            break

        else:
            cur -= len(filenames)

    if file_path and re.search('.mol$',file_path):
        # get the mol file we need
        convert(file_path)

def mian():
    args = sys.argv
    if len(args)>=3:
        base = int(args[1])
        offset = int(args[2])
        run(base,offset)
