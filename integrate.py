import sys,os
import numpy as np
from itertools import chain
import pandas as pd

dest_folder = '/n/scratch2/xl198/data/integrate'
source_path = '/n/scratch2/xl198/data/kaggle/train_set_npy'
test_file ='/n/scratch2/xl198/data/integrate/result.npy'


def test(i):
    look = pd.read_csv(os.path.join(dest_folder,'record.csv'))
    npy = np.load(test_file)
    idx = look.ix[i]['ID']
    path = os.path.join(source_path,'receptors', idx + '.npy')
    start = look.ix[i]['start']
    end = look.ix[i]['end']
    a = npy[start:end]
    b = np.load(path)
    print np.all(a == b)

def make_ligand_index():
    if not os.path.exists(dest_folder):
        os.mkdir(dest_folder)
    npy_list = []
    ids = []
    rec = []
    offset = [0]
    for dirpath,dirnames,filenames in chain(os.walk(os.path.join(source_path,'crystal_ligands')),os.walk(os.path.join(source_path,'docked_ligands'))):
        for filename in filenames:
            #print os.path.join(dirpath,filename)
            npy = np.reshape(np.load(os.path.join(dirpath,filename)),(-1,4))
            ids.append(filename.split('.')[0])
            #ids.append(filename.split('.')[0])
            rec.append(dirpath.split('/')[-1])
            offset.append(npy.shape[0])
            npy_list.append(npy)


    one_file = np.vstack(npy_list)

    np.save(os.path.join(dest_folder,'train_set_ligands.npy'),one_file)
    start = np.array(offset[:-1])
    start = np.cumsum(start)
    end = start + np.array(offset[1:])
    
    start = list(start)
    end = list(end)
    data = []
    for i in range(len(ids)):
	data.append([ids[i],rec[i],start[i],end[i]])

    record = pd.DataFrame(data=data,columns=['ID','receptor','start','end'])
    record.to_csv(os.path.join(dest_folder,'train_set_ligands.csv'),index=False)


def make_receptor_index():
    if not os.path.exists(dest_folder):
        os.mkdir(dest_folder)
    npy_list = []
    ids = []
    rec = []
    offset = [0]
    for dirpath, dirnames, filenames in chain(os.walk(os.path.join(source_path, 'receptors'))):
        for filename in filenames:
            npy = np.reshape(np.load(os.path.join(dirpath, filename)), (-1, 4))
            #ids.append('_'.join(filename.split('.')[0].split('_')[:3]))
            ids.append(filename.split('.')[0])
            #rec.append(dirpath.split('/')[-1])
            offset.append(npy.shape[0])
            npy_list.append(npy)

    one_file = np.vstack(npy_list)

    np.save(os.path.join(dest_folder,'train_set_receptors.npy'),one_file)
    start = np.array(offset[:-1])
    start = np.cumsum(start)
    end = start + np.array(offset[1:])

    start = list(start)
    end = list(end)
    data = []
    for i in range(len(ids)):
        data.append([ids[i], start[i], end[i]])

    record = pd.DataFrame(data=data, columns=['ID', 'start', 'end'])
    record.to_csv(os.path.join(dest_folder, 'train_set_receptors.csv'), index=False)


make_receptor_index()
make_ligand_index()
