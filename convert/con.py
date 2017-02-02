import os,sys
from glob import glob
import getopt
import re
sourceSuffix = '.pdb'
destSuffix = '.mol2'
sourcePath = '/n/scratch2/xl198/data/jan_18_small/kaggle/train_pdb/receptors'
destPath = '/n/scratch2/xl198/data/jan_18_small/kaggle/train_pdb/mol2_receptors'
sys.path.append(os.path.dirname(sys.path[0]))
from util.createfolder import try_create_chain_parent_folder

class FLAGS:
    orchestra_arrayjob = False

def convert(offset = None):
    files = glob(os.path.join(sourcePath),'*.pdb')
    index = range(offset-1,len(files),FLAGS.jobsize) if offset else range(0,len(files))
    for i in index:
        sourceFilePath = files[i]
        destFilePath = sourceFilePath.replace(sourcePath,destPath)
        destFilePath = destFilePath.replace('.pdb','.mol2')
        try_create_chain_parent_folder(destFilePath)
        cmd = 'obabel -ipdb {} -omol2 -O {}'.format(sourceFilePath,destFilePath)
        #os.popen(cmd)
        print "cmd"


def main():
    if FLAGS.orchestra_arrayjob:
        convert(FLAGS.jobid)
    else:
        convert()

def parse_FLAG():
    try:
        opts,args = getopt.getopt(sys.argv[1:],None,["jobsize=","jobid=","cores="])
    except getopt.GetoptError as err:
        # print help information and exit:
        print str(err)  # will print something like "option -a not recognized"

        sys.exit(2)

    for name,value in opts:
        if name == '--jobsize':
            FLAGS.orchestra_jobsize = int(value)
            print "--jobsize ",value
        if name == '--jobid':
            FLAGS.orchestra_jobid = int(value)
            print "--jobid", value
        if name == '--cores':
            FLAGS.cores = int(value)

    if hasattr(FLAGS,"orchestra_jobsize") and hasattr(FLAGS,"orchestra_jobid"):
        FLAGS.orchestra_arrayjob = True

    print "orchestra job ",FLAGS.orchestra_arrayjob

    if hasattr(FLAGS,'cores'):
        print "cores num ",FLAGS.cores

if __name__ == '__main__':
    parse_FLAG()
    main()