import os, sys, getopt
import config
import time
import pandas as pd
from itertools import chain

'''
This code is used to monitor the job on Orchestra

1. Calculate the total number of job
2. Check the number of current job is less than 500 then submit new job
3. To submit job : create a bash file to submit jobarray
'''


def get_job_num():
    # get how many job are there running on Orchestra
    r = os.popen('bjobs | wc -l')
    num = int(r.readline().strip('\n'))
    return num


def create_jobarray(base):
    # create jobarray.sh to submit job on Orchestra
    Jobarray_name = 'jobarray_base_' + str(base) + '.sh'
    Jobarray = os.path.join(config.JOBARRAY, Jobarray_name)
    with open(Jobarray, 'w') as job:
        job.write('#!/bin/bash\n')
        job.write('#BSUB -n 1                #each  job run on 1 core\n')
        job.write('#BSUB -W 240:00            #job run 12 hour\n')
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
        job.write('python %s %s ${LSB_JOBINDEX} \n' % (os.path.join(config._BASE_SCRIPT, 'run_obabel.py'), base))

    return Jobarray


def submit_jobarray(Jobarray, base):
    r = os.popen('bsub < %s' % (Jobarray))
    submit_info = r.readline()
    sys.stderr.write(submit_info)
    sys.stderr.write("submit job %s to %s\n\n" % (base * 1000, (base + 1) * 1000))
    # os.remove(Jobarray)


def get_job_num_get_split():
    count = 0
    # mol2 file organized in input_path
    input_path = os.path.join(config.BASE_DATA, 'result')

    for dirname, dirnames, filenames in os.walk(input_path):
        count += len(filenames)

    return count


def get_job_num_new_line_inserter():
    input_path = config.BASE_YI
    return len(os.listdir(input_path))


def get_csv_size():
    df = pd.read_csv('/home/xl198/remark/dec_1.csv')
    return len(df)


def get_in_folder():
    input_path = '/n/scratch2/xl198/data/H/data'
    count = 0

    for dirname, dirnames, filenames in os.walk(input_path):
        count += len(filenames)

    return count


def get_database_size():
    count = 0
    database_path = ''
    source_ligands_folder = os.path.join(database_path, 'ligands')
    source_crystal_folder = os.path.join(database_path, 'crystal')
    source_receptor_folder = os.path.join(database_path, 'receptor')
    for dirpath, dirname, filenames in chain(os.walk(source_crystal_folder), os.walk(source_ligands_folder),
                                             os.walk(source_receptor_folder)):
        count += len(filenames)

    return count


def check_loop():
    total = get_csv_size()

    # sys.stderr.write("\nConvert mol2 into pdb\n")
    sys.stderr.write("total commands num : %s\n" % total)
    cur = 0
    base = 0
    # i = 1
    # base = offset + 1000 if offset + 1000 <= docking_num else docking_num
    while (1):
        # i = 0
        if base * 1000 >= total:
            sys.stderr.write('\nFinish\n')
            exit(0)

        num = get_job_num()
        if num < 500:
            Jobarray = create_jobarray(base)
            submit_jobarray(Jobarray, base)

            # end = offset + 1000 if offset + 1000 <= docking_num else docking_num

            base += 1

        time.sleep(6)


if __name__ == '__main__':
    check_loop()
