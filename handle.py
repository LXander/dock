import os, sys , getopt
import config
import time

def get_job_num():
    # get how many job are there running on Orchestra
    r = os.popen('bjobs | wc -l')
    num = int(r.readline().strip('\n'))
    return num

def create_jobarray(start, end, offset, script_index, debug_flag):
    Jobarray_name = 'jobarray_' + str(offset) + '.sh'
    Jobarray = os.path.join(config.JOBARRAY,Jobarray_name)
    with open(Jobarray,'w') as job:
        job.write('#!/bin/bash\n')
        job.write('#BSUB -n 1                #each  job run on 1 core\n')
        job.write('#BSUB -W 12:00            #job run 12 hour\n')
        job.write('#BSUB -J jobArray[1-1000] #job array list goes begin,begin+1,begin+2...end\n')
        job.write('#BSUB -o '+config.JOBARRAY+'/%J.%I.out        #lsf output file\n')
        job.write('#BSUB -e '+config.JOBARRAY+'/%J.%I.err       #lsf error file\n')
        job.write('#BSUB -q short         #submit to "short" queue\n')
        job.write('export LC_CTYPE=en_US.UTF-8\n')
        job.write('export LC_ALL=en_US.UTF-8\n')
        job.write('LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/xl198/local/lib:/home/xl198/boost/boost_1_54_0/stage/lib\n')
        job.write('export OMP_NUM_THREADS=1\n')
        job.write('export LC_ALL="en_US.UTF-8"\n')
        job.write('source /home/xl198/venv/data/bin/activate\n')
        job.write('python %s %s %s ${LSB_JOBINDEX} %s\n'%(config.RUN,script_index, offset,'-d'if debug_flag else ''))


    return Jobarray

def submit_jobarray(Jobarray, start, end):
    r = os.popen('bsub < %s'%(Jobarray))
    submit_info = r.readline()
    sys.stderr.write(submit_info)
    sys.stderr.write("submit job %s to %s\n\n" % (start, end))
    #os.remove(Jobarray)

def extract_file(script_index):
    r = os.popen('python %s filter %s'%(config.FILTER,config.ROW_FILE[script_index]))
    for line in r.readlines():
        sys.stderr.write(line)

def check_loop(script_index,debug_flag):
    extract_file(script_index)
    docking_num = 0
    with open(config.COMMANDS_FILE[script_index]) as fr:
        for line in fr:
            docking_num += 1

    sys.stderr.write("\nopen script file : %s\n"%config.COMMANDS_FILE[script_index])
    sys.stderr.write("total commands num : %s\n"%docking_num)
    offset = 0
    start = 1
    end = offset + 1000 if offset + 1000 <= docking_num else docking_num
    while(1):
        if offset >= docking_num:
            sys.stderr.write('\nFinish\n')
            exit(0)
        num = get_job_num()
        if num<500:
            Jobarray = create_jobarray(start,end, offset,script_index,debug_flag)
            submit_jobarray(Jobarray,start+offset,end)
            offset = end
            end = offset + 1000 if offset + 1000 <= docking_num else docking_num

        time.sleep(600)


def main():
    # eg : python handle.py 1 -d
    try:
        opts, args = getopt.getopt(sys.argv[1:], "d", ["help", "output="])
    except getopt.GetoptError as err:
        # print help information and exit:
        print str(err)  # will print something like "option -a not recognized"

        sys.exit(2)

 
    script_index = int(args[0])
    
    # if -d display debug information
    debug_flag = False
    for o, a in opts:
        if o == '-d':
            debug_flag = True

    check_loop(script_index,debug_flag)


if __name__ == '__main__':
    main()



