import os,sys
import pandas as pd
import numpy  as np
import tempfile
import re
from glob import glob
import prody


def merge():

    fast_ligand = os.path.join(FLAGS.fas_path,Id+'_fast.pdb')
    rigor_ligand = os.path.join(FLAGS.rigor_path,Id+'_rigor.pdb')
    ritor_so_ligand = os.path.join(FLAGS.rigor_so_path,Id+'_rigor_so.pdb')





def tanimoto(item):

    Id = item['inner_index']

    sourceFileName = item['ID'] + '_' +'ligand' + '_' + item['dataset'] + '.pdb'
    sourFilePath = os.path.join(FLAGS.sourcePath,item['PDBname'],sourceFileName)

    destFileName = item['ID'] + '_' + item['dataset']+ '_' + str(item['inner_index'])+'.pdb'
    destFilePath = os.path.join(FLAGS.tempPath,destFileName)

    cmd = 'obabel -ipdb {} -f {} -l {} -opdb -O {}'.format(sourFilePath, Id, Id, destFilePath)

    os.system(cmd)

    # Calculate tanimoto similarity for all crystal ligand for the receptor
    crystal_list = glob(os.path.join(FLAGS.crystalPath),item['PDBname'],"*.pdb")

    tanimoto_list =[]

    for crystal_ligand in crystal_list:
        command = os.popen('babel -d {} {} -ofpt -xfFP4'.format(crystal_ligand,destFilePath))
        ls =command.read()
        try:
            tanimoto_similarity = re.split('=|\n', ls)[2]
            tanimoto_list.append(tanimoto_similarity)
        except:
            sys.stderr.write(ls)
            tanimoto_list.append( FLAGS.tanimoto_cutoff+0.1)

    return min(tanimoto_list)



def fusion():
    #Create tempdir to store splited result
    FLAGS.tempPath =tempfile.mkdtemp()

    df = pd.read_csv(FLAGS.recordFile)
    if FLAGS.tanimoto not in df.columns:
        df[FLAGS.tanimoto] = np.nan

    null_index = df[df[FLAGS.tanimoto].isnull()].index
    if len(null_index) > FLAGS.batchsize:
        null_index = null_index[:FLAGS.batchsize]

    df[FLAGS.tanimoto][null_index] = df.ix[null_index].apply(tanimoto, axis=1)

    df.to_csv(FLAGS.recordFile,index=False)




class FLAGS:
    recordFile = '/n/scratch2/xl198/data/fusion/forms/filter_fast.csv'
    tanimoto = 'tanimoto'
    tanimoto_cutoff = 0.75
    sourcePath = '/n/scratch2/xl198/data/pdbs'
    crystalPath = '/n/scratch2/xl198/data/H/addH'
    batchsize = 10


fusion()


