import os, sys , getopt
import config
import time

def get_job_num():
    # get how many job are there running on Orchestra
    r = os.popen('bjobs | wc -l')
    num = int(r.readline().strip('\n'))
    return num

def create_jobarray(base):
    Jobarray_name = 'jobarray_base_' + str(base) + '.sh'
    Jobarray = os.path.join(config.JOBARRAY,Jobarray_name)
    with open(Jobarray,'w') as job:
        job.write('#!/bin/bash\n')
        job.write('#BSUB -n 1                #each  job run on 1 core\n')
        job.write('#BSUB -W 10            #job run 12 hour\n')
        job.write('#BSUB -J jobArray[1-1000] #job array list goes begin,begin+1,begin+2...end\n')
        job.write('#BSUB -o /dev/null        #lsf output file\n')
        job.write('#BSUB -e /dev/null       #lsf error file\n')
        job.write('#BSUB -q mini         #submit to "short" queue\n')
        job.write('export LC_CTYPE=en_US.UTF-8\n')
        job.write('export LC_ALL=en_US.UTF-8\n')
        job.write('LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/xl198/local/lib:/home/xl198/boost/boost_1_54_0/stage/lib\n')
        job.write('export OMP_NUM_THREADS=1\n')
        job.write('export LC_ALL="en_US.UTF-8"\n')
        job.write('source /home/xl198/venv/data/bin/activate\n')
        job.write('python %s %s ${LSB_JOBINDEX} \n'%(config.CONVERT, base))


    return Jobarray

def submit_jobarray(Jobarray, base):
    r = os.popen('bsub < %s'%(Jobarray))
    submit_info = r.readline()
    sys.stderr.write(submit_info)
    sys.stderr.write("submit job %s to %s\n\n" % (base*1000, (base+1)*1000))
    #os.remove(Jobarray)

def check_loop():
    base = config.BASE_YI
    files = os.listdir(base)
    total = len(files)

    sys.stderr.write("\nConvert mol2 into pdb\n")
    sys.stderr.write("total commands num : %s\n"%total)
    cur = 0
    base = 0
    #base = offset + 1000 if offset + 1000 <= docking_num else docking_num
    while(1):
        if base*1000 >= total:
            sys.stderr.write('\nFinish\n')
            exit(0)
        num = get_job_num()
        if num<500:
            Jobarray = create_jobarray(base)
            submit_jobarray(Jobarray,base)

            #end = offset + 1000 if offset + 1000 <= docking_num else docking_num

        time.sleep(6)

if __name__ == '__main__':
    check_loop()