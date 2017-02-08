import os,sys
from glob import glob

def add_hydrogen(filePath):
    fileName = os.path.basename(filePath)
    receptorName=  fileName.split('.')[0]

    with open('repair.sh', 'w') as w:
        w.write('# !/bin/bash\n')
        w.write('# BSUB -n 2\n')
        w.write('# BSUB -W 100:00\n')
        w.write('# BSUB -J reapir_{}\n'.format(receptorName))
        w.write('# BSUB -o /home/yw174/job/addH/{}.out\n'.format(receptorName))
        w.write('# BSUB -e /home/yw174/job/addH/{}.err\n'.format(receptorName))
        w.write('# BSUB -q long\n')
        #w.write('export PATH =$PATH:/home/yw174/usr/babel/bin/\n')
        #w.write('cd /home/yw174/program/pdb_sth\n')
        cmd = 'obminimize -cg -ff MMFF94 -h -n 500 {} > {}'.format(filePath,filePath.replace("protein","repair"))
        w.write(cmd + '\n')
    os.system('bsub < repair.sh')

def add_all():
    for receptor in glob(os.path.join('/n/scratch2/xl198/dude/data/receptor/protein','*.pdb')):
        add_hydrogen(receptor)




def add_one():
    receptor = glob(os.path.join('/n/scratch2/xl198/dude/data/receptor/protein','*.pdb'))[0]
    add_hydrogen(receptor)



add_all()
