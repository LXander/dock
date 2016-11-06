import os, sys
import random
from datetime import datetime
import subprocess
import time


class pickupexample:
    def __init__(self, db_path=None, output_path=None):

        # set db_path
        if db_path == None:
            sys.stderr.write("must set db_path")
            exit(1)
        else:
            self.db_path = db_path

        # set output_path
        if output_path == None:
            sys.stderr.write("must set output_path")
            exit(1)
        else:
            self.output_path = output_path

    def log(self, content, debug):
        # write debug information to log
        if debug:
            sys.stderr.write(content + '\n')

    def get_picked_list(self, num, debug=False):
        '''

        :param num: int , the number of files picked from database
        :param debug: boolean , if write log
        :return: a list of picked filename and output filename used to convert
        '''
        random.seed(datetime.now())

        folders = os.listdir(self.db_path)
        folders_num = len(folders)
        self.log("Open data base at %s" % self.db_path, debug)
        self.log("%d folders in data base" % folders_num, debug)

        picked_file = []

        for i in range(num):
            # randomly select folder
            folder = folders[random.randint(0, folders_num - 1)]
            folder_path = os.path.join(self.db_path, folder)

            # randomly select mol2 file
            files = os.listdir(folder_path)
            if len(files) == 1:
                file_name = files[0]
            else:
                file_name = files[random.randint(0, len(files) - 1)]

            # a pair: ['input_path/abcd_1234_lignad.mol2','output_path/abcd_1234_ligand.pdb']
            picked_pair = [
                os.path.join(folder_path, file_name),
                os.path.join(self.output_path, file_name.split('.')[0] + '.pdb')
            ]

            picked_file.append(picked_pair)

            self.log('pickup file %s' % file_name, debug)

        return picked_file

    def get_all_picked_list(self, num, debug=False):
        '''

        :param num: int , the number of files picked from database
        :param debug: boolean , if write log
        :return: a list of picked filename and output filename used to convert
        '''
        random.seed(datetime.now())

        folders = os.listdir(self.db_path)
        folders_num = len(folders)
        self.log("Open data base at %s" % self.db_path, debug)
        self.log("%d folders in data base" % folders_num, debug)

        picked_file = []

        for i in range(num):
            # randomly select folder
            folder = folders[i]
            folder_path = os.path.join(self.db_path, folder)

            # randomly select mol2 file
            files = os.listdir(folder_path)
            for file_name in files:
                # a pair: ['input_path/abcd_1234_lignad.mol2','output_path/abcd_1234_ligand.pdb']
                picked_pair = [
                    os.path.join(folder_path, file_name),
                    os.path.join(self.output_path, file_name.split('.')[0] + '.pdb')
                ]

                picked_file.append(picked_pair)

                # self.log('pickup file %s' % file_name, debug)

        return picked_file

    def get_converted_list(self, picked_file, debug):

        self.log('run convert command:\nobabel -i mol2 -m %s -o pdb -O %s' % (picked_file[0][0], picked_file[0][1]),
                 debug)

        jobs = []
        '''
	# use obabel to convert mol2 into seperate pdb file
        for input,output in picked_file:
            job = subprocess.Popen('obabel -i mol2 -m %s -o pdb -O %s'%(input,output), shell=True,
                                      stdout=subprocess.PIPE)
            jobs.append(job)

        # waiting for all convert job finish
        for job in jobs:
            job.wait()
	'''


        for input, output in picked_file:
            subprocess.call('obabel -i mol2 -f 1 -l 1 %s -o pdb -O %s' % (input, output), shell=True)
            # subprocess.call('obabel -i mol2 -m %s -o pdb -O %s'%(input,output),shell=True)

        filename_list = os.listdir(self.output_path)

        self.log('convert finished, total output', debug)

        return filename_list

    def pickup(self, num, debug=False):
        picked_file = self.get_picked_list(num, debug)
        converted_file = self.get_converted_list(picked_file, debug)
        return converted_file

    def pickup_all(self, num, debug=False):
        picked_file = self.get_all_picked_list(num, debug)
        converted_file = self.get_converted_list(picked_file, debug)


def test():
    p = pickupexample('/n/scratch2/xl198/data/H/wp_fast', '/home/xl198/playground/bucket')
    start = time.time()
    name_list = p.pickup_all(100, debug=True)
    end = time.time()

    print "totla time : %s" % (end - start)


test()
