import os,sys
import pandas as pd
import numpy as np
import tempfile
import re
from glob import glob
import mdtraj as md
import getopt




def tanimoto(crystalFolder,ligandPath):

    file_name= os.path.basename(ligandPath)
    temp_file_path = os.path.join(FLAGS.tempPath,file_name)
    obabel_cmd = 'obabel -ipdb {} -l 1 -f 1 -opdb -O {}'.format(ligandPath,temp_file_path)
    os.system(obabel_cmd)


    crystalList = glob(os.path.join(crystalFolder,'*.pdb'))

    records = []

    for crystal_ligand in crystalList:
        command = os.popen('babel -d {} {} -ofpt -xfFP4'.format(crystal_ligand,temp_file_path))
        ls = command.read()
        try:
            tanimoto_similarity = re.split('=|\n', ls)[2]




def convert(fast_file_path):
    fast_name = os.path.basename(fast_file_path)
    receptor = fast_name.split('_')[0]
    original_name = fast_name.replace('ligand_fast','ligand')
    rigor_name= fast_name.replace('ligand_fast','ligand_rigor')
    rigor_so_name = fast_name.replace('ligand_fast','ligand_rigor_so')


    crystalFolder = os.path.join(FLAGS.crystalPath,receptor)





def run():
    FLAGS.tempPath = tempfile.mkdtemp()

    fastFileList = glob(os.path.join(FLAGS.fastPath,'*','*[fast]*.pdb'))

    index = range(len(fastFileList))

    if hasattr(FLAGS,'arrayjob') and FLAGS.arrayjob:
        if hasattr(FLAGS,'offset'):
            index = [len(index)/FLAGS.jobsize*FLAGS.offset+FLAGS.jobid]
        else:
            index = range(FLAGS.jobid-1,len(index),FLAGS.jobsize)

    for i in index:
        convert(fastFileList[i])


class FLAGS:
    fastPath = ''

    tanimoto_cutoff = 0.75

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

    print "orchestra job ",FLAGS.arrayjob

    if hasattr(FLAGS,'cores'):
        print "cores num ",FLAGS.cores