import config
import os,sys
import getopt
#import subprocess
import time

'''
    This file is used to get the result form Yi
    and then organized them into folder
    and insert newline into it otherwise it can't be read into obabel

'''



def convert(raw_file_path):
    print "convert"
    file_name = os.path.basename(raw_file_path)
    receptor_name = file_name.split('_')[0]
    output_path = os.path.join(FLAGS.destPath, receptor_name)
    print output_path
    print os.path.exists(output_path)
    if not os.path.exists(output_path):
        os.mkdir(output_path)

    input_file_path = os.path.join(FLAGS.mediaPath,receptor_name,file_name)
    output_file_path = os.path.join(output_path,file_name.split('.')[0] + '.pdb')
    os.system('obabel -i mol2 %s -o pdb -O %s'%(input_file_path,output_file_path))


def run(input_file_path):
    
    file_name = os.path.basename(input_file_path)
    receptor_name = file_name.split('_')[0]
    output_path = os.path.join(FLAGS.mediaPath,receptor_name)
    if not os.path.exists(output_path):
        os.mkdir(output_path)

    output_file_path = os.path.join(output_path,file_name)
    with open(input_file_path) as infile:
        with open(output_file_path,'w') as outfile:
            for line in infile:
                outfile.write(line)
                if line == '@<TRIPOS>MOLECULE\n':
                    outfile.write('\n')

def get_all(num = None):
    base = config.BASE_YI
    files = os.listdir(base)
    size = num if num != None else len(files)
    sys.stderr.write("Convert %s files\n"%size)
    sys.stderr.write("first file is %s\n"%(os.path.join(base,files[0])))
    for i in range(size):
        run(os.path.join(base,files[i]))
        sys.stderr.write("write %d/%d\n"%(i+1,size))

def old_run_convert(base,offset):

    base_path = FLAGS.sourceBase
    files = os.listdir(base_path)
    #print base
    #print offset
    index = base*1000+offset

    if len(files)>index:
        run(os.path.join(base_path,files[index]))
        convert(os.path.join(base_path,files[index]))

def run_convert():
    #base_path = config.BASE_YI
    base_path = FLAGS.sourceBase
    files = os.listdir(base_path)
    #print base
    #print offset
    #index = base*1000+offset
    indexes = range(FLAGS.jobid-1,len(files),FLAGS.jobsize)
    for index in indexes :
        print indexes
        run(os.path.join(base_path,files[index]))
        convert(os.path.join(base_path,files[index]))

class FLAGS:
    arrayjob = False
    sourceBase ="/n/scratch2/xl198/YI/rigor_so/rigor_so"
    mediaPath = "/n/scratch2/xl198/YI/rigor_so/media"
    destPath = "/n/scratch2/xl198/YI/rigor_so/final"

def parse_FLAG():
    try:
        opts,args = getopt.getopt(sys.argv[1:],None,["jobsize=","jobid=","cores="])
    except getopt.GetoptError as err:
        # print help information and exit:
        print str(err)  # will print something like "option -a not recognized"

        sys.exit(2)

    for name,value in opts:
        if name == '--jobsize':
            FLAGS.jobsize = int(value)
            print "--jobsize ",value
        if name == '--jobid':
            FLAGS.jobid = int(value)
            print "--jobid", value
        if name == '--cores':
            FLAGS.cores = int(value)

    if hasattr(FLAGS,"orchestra_jobsize") and hasattr(FLAGS,"orchestra_jobid"):
        FLAGS.arrayjob = True

    print "orchestra job ",FLAGS.arrayjob

    if hasattr(FLAGS,'cores'):
        print "cores num ",FLAGS.cores

def main():
    '''
    args = sys.argv
    base  = int(args[1])
    offset = int(args[2])
    '''
    parse_FLAG()

    run_convert()
    #sys.stderr.write("run convert %s"%(base*1000+offset))
    '''
    if len(args)<2:
        num = None
    else:
        num = int(args[1])

    print num
    get_all(num)
    '''

if __name__ == '__main__':
   main()
