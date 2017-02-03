import os,sys
from glob import glob
import re
import getopt

pattern = re.compile(r'{(?P<key>[^\{\}]*)\:(?P<value>[^\{\}]*)}')

full_columns = ['Autovina_gauss2','PDBResId','RMSF','center','Autovina_Hydrogen','Autovina_repulsion','PDBname','Contact Similarity','Resolution(A)','rotation','Autovina_gauss1','Autovina_Affinity(kcal/mol)','Autovina_hydrophobic']

selected_columns = ['PDBname','PDBResId','RMSF','Autovina_Affinity(kcal/mol)']

def parse_remarks(filename,filedir=None):
    if filedir is None:
        fileloc=filename
    else:
        fileloc = os.path.join(filedir,filename)

    remarks=[]
    with open(fileloc,'r') as o:
        filelines= o.readlines()
        for eachline in filelines:
            if 'Remark:' in eachline:
                #print(eachline)
                rdict ={}
                for m in pattern.finditer(eachline):
                    result =m.groupdict()
                    rdict[result['key']]=result['value']

                remarks = [rdict[key] for key in selected_columns]
                with open(os.path.join(FLAGS.tempPath, 'temp_' + str(FLAGS.jobid) + '.txt'), 'a') as fout:
                    fout.write(','.join(remarks)+'/n')


def run():
    filelist= glob(os.path.join(FLAGS.sourceBase,'*.mol'))
    filelist = list(filelist)
    indexes = range(FLAGS.jobid-1,len(filelist),FLAGS.jobsize)
    for index in indexes:
        parse_remarks(filelist[index])

def merge():
    filelist = glob(os.path.join(FLAGS.tempPath,'*.txt'))
    filelist = list(filelist)
    with open(os.path.join(os.path.dirname(FLAGS.tempPath),'forms.csv')) as fout:
        for f in filelist:
            with open(f) as fin:
                for line in fin:
                    fout.write(line)

class FLAGS:
    arrayjob = False
    sourceBase ="/n/scratch2/xl198/YI/rigor_so/rigor_so"
    tempPath = '/n/scratch2/xl198/YI/rigor_so/temp'


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

