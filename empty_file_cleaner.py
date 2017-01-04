from __future__ import print_function 
import os,sys,re
import pandas as pd

def record():
    records = []
    for dirpath,dirname,filenames in os.walk(input_path):
        for filename in filenames:
            filepath = os.path.join(dirpath,filename)
            with open(filepath) as fin:
                count = len(fin.readlines())
      		if not count:
		    records.append(filepath)

    df = pd.DataFrame(data = records,columns = ['Path'])
    df.to_csv(os.path.join(input_path,'record.csv'), index = False)

def remove():
    df = pd.read_csv(os.path.join(input_path,'record.csv'))
    rm = lambda x:os.popen("rm "+x)
    #rm = lambda x:print("rm "+x)
    #df['Path'].apply(rm)
    for i in range(len(df)):
        rm(df.ix[i]['Path'])
    

input_path = '/n/scratch2/xl198/dude/data/jan_05'

record()
