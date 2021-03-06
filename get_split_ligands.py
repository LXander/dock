import os, sys
import config
import re
import shutil
import util

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
# if atoms num too small it won't be too different between each other
minimum_atom_num = 100
# mol2 file organized in input_path
input_path = os.path.join(config.BASE_DATA, 'result')
# select and convert pdb file to output_path
output_path = os.path.join(config.BASE_DATA, 'select', 'ligands')


def convert(input_file):
    # convert from mol2 file into split pdb format

    # count the number of the single ligand in one file
    ligand_count = 0
    with open(input_file) as input:
        for line in input:
            if line == '@<TRIPOS>MOLECULE\n':
                ligand_count += 1

    print "ligand_count %d" % ligand_count


    # count the number of the atoms in ligand
    flag = False
    atom_num = 0
    with open(input_file) as input:
        for line in input:
            if not flag and line == '@<TRIPOS>ATOM\n':
                flag = True
            elif line == '@<TRIPOS>BOND\n':
                break
            elif flag:
                atom_num += 1

    # prepare the output path and file name
    dirname = os.path.dirname(input_file)
    receptor = os.path.basename(dirname)
    filename = os.path.basename(input_file).split('.')[0]

    # creat the folder if dones't exists
    # ( [output_path_base]/[receptor_name]  )
    # eg otuput_path_base [/n/scratch2/xl198/data/result/liagnds/]  receptor_name = [1a3b]
    #    folder [/n/scratch2/xl198/data/result/ligands/1a3b]
    ligand_output_path = os.path.join(output_path, receptor)
    if not os.path.exists(ligand_output_path):
        os.mkdir(ligand_output_path)

    # name the top 10 and bottom 10 to distinguish them
    top = filename + '_top_.pdb'
    bottom = filename + '_bottom_.pdb'
    print ligand_count
    print atom_num
    if ligand_count >= minimum_ligand_num:
        # convert top 10 ligand
        tmp1 = os.path.join('/tmp', top)

        cmd1 = 'obabel -imol2 %s -f %d -l %d -opdb -O %s -d ' %\
              (input_file, 1, 10, tmp1)

        cmd2 = 'obabel -ipdb -m %s  -opdb -O %s -d ' %\
              (tmp1, os.path.join(ligand_output_path,top))
	#print(cmd1)
	#print(cmd2)
        os.system(cmd1)
        os.system(cmd2)


        # convert bottom 10 ligand
        tmp2 = os.path.join('/tmp', bottom)

        cmd3 = 'obabel -imol2 %s -f %d -l %d -opdb -O %s -d ' % \
              (input_file, ligand_count - 9, ligand_count, tmp2)

        cmd4 = 'obabel -ipdb -m %s  -opdb -O %s -d ' % \
               (tmp2, os.path.join(ligand_output_path, bottom))


        #print cmd4
        os.system(cmd3)
        os.system(cmd4)

def get_pdb_and_crystal(input_file):
    # source to place crystal ligand
    crystal_source = os.path.join(config.BASE_DATA, 'H', 'addH')
    # dest to place converted ligand
    crystal_dest = os.path.join(config.BASE_DATA, 'select', 'crystal')

    # source to place pdb
    pdb_source = os.path.join(config.BASE_DATA, 'H', 'data')
    # dest to place pdb
    pdb_dest = os.path.join(config.BASE_DATA, 'select', 'receptor')

    dirname = os.path.dirname(input_file)
    receptor = os.path.basename(dirname)
    filename = os.path.basename(input_file).split('.')[0][:-5]

    receptor_in = os.path.join(pdb_source,receptor,receptor+'.pdb')
    shutil.copy(receptor_in,pdb_dest)

    crystal_in = os.path.join(crystal_source,receptor,filename+'.pdb')
    crystal_path = os.path.join(crystal_dest,receptor)
    if not os.path.exists(crystal_path):
        os.mkdir(crystal_path)

    crystal_out = os.path.join(crystal_path,filename+'.pdb')

    cmd = 'obabel -ipdb %s -opdb -O %s -d'%(crystal_in,crystal_out)
    os.system(cmd)

def check_avail(file_path):
    ID = util.get_id(os.path.basename(file_path))
    avail_list = []
    with open(config.AVAIL) as fr:
        for line in fr:
            avail_list.append(line.strip('\n'))
    return ID in avail_list


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
        flag = check_avail(file_path)
        if flag:
            # get the mol file we need
            print "get"
            convert(file_path)
            get_pdb_and_crystal(file_path)


def main():
    args = sys.argv
    if len(args) >= 3:
        base = int(args[1])
        offset = int(args[2])
        print 'base %d offset %d' % (base, offset)
        run(base, offset)


if __name__ == '__main__':
    main()
