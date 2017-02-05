import os,sys
import pandas as pd
import numpy  as np
import tempfile
import re


def tanimoto(item):

    Id = item['inner_index']

    sourceFileName = item['ID'] + '_' +'ligand' + '_' + item['dataset'] + '.pdb'
    sourFilePath = os.path.join(FLAGS.sourcePath,item['PDBname'],sourceFileName)

    destFileName = item['ID'] + '_' + item['dataset']+ '_' + item['inner_index']+'.pdb'
    destFilePath = os.path.join(FLAGS.tempPath,destFileName)

    cmd = 'obabel -ipdb {} -f {} -l {} -opdb -O {}'.format(sourFilePath, Id, Id, destFilePath)

    os.system(cmd)

    command = os.popen('babel -d {} {} -ofpt -xfFP4'.format(sourFilePath,destFilePath))
    ls =command.read()
    try:
        tanimoto_similarity = re.split('=|\n', ls)[2]
    except:
        sys.stderr.write(ls)
        return FLAGS.tanimoto_cutoff+0.1
    return tanimoto_similarity


def fusion():
    #Create tempdir to store splited result
    FLAGS.tempPath =tempfile.mkdtemp()

    df = pd.read_csv(FLAGS.recordFile)
    if FLAGS.tanimoto not in df.columns:
        df[FLAGS.tanimoto] = np.nan

    null_index = df[df[FLAGS.tanimoto].isnull()]
    if len(null_index) > FLAGS.batchsize:
        null_index = null_index[:FLAGS.batchsize]

    df[null_index, FLAGS.tanimoto] = df[null_index].apply(tanimoto, axis=1)

    df.to_csv(FLAGS.recordFile,index=False)




class FLAGS:
    recordFile = '/n/scratch2/xl198/data/fusion/forms/filter_fast.csv'
    tanimoto = 'tanimoto'
    tanimoto_cutoff = 0.75
    sourcePath = '/n/scratch2/xl198/data/pdbs'
    batchsize = 10


fusion()


