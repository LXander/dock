import numpy as np
import pandas as pd
import re
import os,sys

#
# This code is used to get the target ligand meet the demand we want, and convert it
# there are remark in experimental data and docking result in different form
# according to the remark we can get the one we want
#
# and since mol2 doesn't have the newline so it can't convert to pdb directly
# thus we need first write them and then convert
#
# when convert by obabel , slice by -f and -l will help to get the one we want
#
# another thing is their name, sometime we need to convert only one and some time
# multiple of them, but then using -m to split the result, it will append the suffix
# by number and then we should do something to make it more readable
# also it need to run obabel with different parameter so may write different function to
# make code easier
#
# and if we need to get some top of the ligand or bottom of it, so when naming the result,
# we have to take care to make it more readable
#
# it's running on Orchestra, to take most advantage of it, once we could submit at most 1,000
# job and when most of them finished, submit another batch of it. if require us to use base and
# offset to identify the target ligand need to convert. The code to write jobarray.sh and submit it
# and monitor the number of job running are wrote in another file
#
#
# The mainly part of this class is
# get_sourcefile_path
# then parse the remark to see if meet our reqiurement
#
#
#
#

class convert:
    def __init__(self):
        # used as the source database
        self.input_path = '/n/scratch2/yw174/result/fast'


    def get_id(self,filename):
        '''
        extract the ID of a file
        :param filename: eg '1a2b_1234.pdb'
        :return: eg '1a2b_1234'
        '''
        ID = '_'.join(filename.split('_')[0:2])
        return ID

    def count_liangd_num(self,input_file):
        '''
        count the number of ligands for given ligands in mol2 foramt.

        :param input_file: path of input file eg. '/n/scratch2/xl198/data/source/ligands/1a2b/1a2b_1234_liagnd.mol2'
        :return:  [ID,atom_num] eg [1a2b_1234,10]
        '''

        # count the number of the single ligand in one file
        ligand_count = 0
        with open(input_file) as input:
            for line in input:
                if line == '@<TRIPOS>MOLECULE\n':
                    ligand_count += 1


        ID =self.get_id(os.path.basename(input_file))

        return [ID, ligand_count]

    def count_atom_num(self,input_file):
        '''
        count the number of atoms for given ligand in mol2 foramt.

        :param input_file: path of input file eg. '/n/scratch2/xl198/data/source/ligands/1a2b/1a2b_1234_liagnd.mol2'
        :return: [filename,atom_num] eg [1a2b_1234,10]
        '''

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

        ID = self.get_id(os.path.basename(input_file))

        return [ID, atom_num]

    def standarlize_exp(self,data):
        '''
        only first two columns is strng, and the rest are float,
        and some of the rest would have multiple value, should split them
        it's hard to convert with operator like gt or lt so we discard them at this time
        :param data: list of the remark
        :return: standarlized data
        '''

        result = []
        for i in range(len(data)):
            if i <2:
                # they're string so keep it as original
                result.append(data[i])
            elif data[i] == '':
                # Nothing here
                result.append(None)
            else:

                m = re.match('\[(.*)\]',data[i])
                if m :
                    # it was a list
                    nums = m.group(1)
                    nums = nums.split(',')
                    nums = [float(n) for n in nums]
                    result.append(nums)
                else:
                    nums = data[i].split('|')
                    nums = [ float(num) for num in nums if not re.search('>|<',num)]
                    result.append(nums)

        return result





    def parse_exp_remark(self,remark):
        '''
        parse the remark of expriment data
        :param remark:
        :return: a dict for the result, since
        '''
        remark = remark.strip('\n')
        remark = remark.split('*@*')
        remark = [r.split(':') for r in remark[1:]]
        columns = [r[0] for r in remark]
        data = [r[1] for r in remark]
        data = self.standarlize_exp(data)


    def run(self,input_file_path):
        #check if it is the file we want to use
        flag = self.check(input_file_path)