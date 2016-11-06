import os, sys
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
minimum_ligand_num = 20
# mol2 file organized in input_path
input_path = os.path.join(config.BASE_DATA, 'result')
# select and convert pdb file to output_path
output_path = os.path.join(config.BASE_DATA, 'select', 'ligands')


def convert(input_file):
    ligand_count = 0
    with open(input_file) as input:
        for line in input:
            if line == '@<TRIPOS>MOLECULE\n':
                ligand_count += 1

    print "ligand_count %d" % ligand_count

    dirname = os.path.dirname(input_file)
    receptor = os.path.basename(dirname)
    filename = os.path.basename(input_file).split('.')[0]

    ligand_output_path = os.path.join(output_path, receptor)
    if not os.path.exists(ligand_output_path):
        os.mkdir(ligand_output_path)

    top = filename + '_top_.pdb'
    bottom = filename + '_bottom_.pdb'

    if ligand_count > minimum_ligand_num:
        # convert top 10 ligand
        tmp1 = os.path.join('/tmp', top)

        cmd1 = 'obabel -imol2 %s -f %d -l %d -opdb -O %s -d ' %\
              (input_file, 1, 10, tmp1)

        cmd2 = 'obabel -imol2 -m %s  -opdb -O %s -d ' %\
              (tmp1, os.path.join(ligand_output_path,top))

        os.system(cmd1)
        os.system(cmd2)


        # convert bottom 10 ligand
        tmp2 = os.path.join('/tmp', bottom)

        cmd3 = 'obabel -imol2 %s -f %d -l %d -opdb -O %s -d ' % \
              (input_file, ligand_count - 9, ligand_count, tmp2)

        cmd4 = 'obabel -imol2 -m %s  -opdb -O %s -d ' % \
               (tmp2, os.path.join(ligand_output_path, bottom))


        #print cmd4
        os.system(cmd3)
        os.system(cmd4)


def run(base, offset):
    # JobId start at 1
    cur = base * 1000 + offset - 1

    file_path = None

    for dirname, dirnames, filenames in os.walk(input_path):
        if len(filenames) > cur:
            if cur >= 0:
                file_path = os.path.join(dirname, filenames[cur])
            break

        else:
            cur -= len(filenames)

    # print 'file_path '+file_path

    if file_path and re.search('.mol$', file_path):
        # get the mol file we need
        print "convert"
        convert(file_path)


def main():
    args = sys.argv
    if len(args) >= 3:
        base = int(args[1])
        offset = int(args[2])
        print 'base %d offset %d' % (base, offset)
        run(base, offset)


if __name__ == '__main__':
    main()
