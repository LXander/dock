import sys,os
import numpy as np
from itertools import chain
import pandas as pd

dest_folder = '/n/scratch2/xl198/data/integrate'
source_path = '/n/scratch2/xl198/data/kaggle_dec_04/test_set_npy'

def make_index():
    if not os.path.exists(dest_folder):
        os.mkdir(dest_folder)
    npy_list = []
    ids = []
    offset = [0]
    for dirpath,dirnames,filenames in chain(os.walk(os.path.join(source_path,'receptors'))):
        for filename in filenames:
            npy = np.load(os.path.join(dirpath,filename))
            #ids.append('_'.join(filename.split('.')[0].split('_')[:3]))
            ids.append(filename.split('.')[0])
            offset.append(npy.shape[0])
            npy_list.append(npy)


    one_file = np.vstack(npy_list)

    one_file.tofile(os.path.join(dest_folder,'result'))
    start = np.array(offset[:-1])
    start = np.cumsum(start)
    end = start + np.array(offset[1:])
    record = pd.DataFrame(data = [ids,start,end],columns=['ID','start','end'])
    record.to_csv(os.path.join(dest_folder,'record.csv'),index=False)


make_index()