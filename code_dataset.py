import pandas as pd
import numpy as np
import random
import os,sys,re
import shutil

crystal_folder_name = 'crystal'
ligands_folder_name = 'ligands'
receptor_folder_name = 'receptor'

database_folders = [crystal_folder_name,ligands_folder_name,receptor_folder_name]

def dump_dataframe(dataframe, filename):
    dataframe.to_csv(filename+'.csv',index=False)

def dump_answer_list(answer_list,filename):
    with open(filename+'.txt','w') as fr:
        for row in answer_list:
            line = row[0] + ','
            line += ' '.join(row[1:])
            line += '/n'
            fr.write(line)


def create_folder(folder):
    if not os.path.exists(folder):
        os.mkdir(folder)

def copy_original(source,dest,raw):
    append_dest = lambda x:os.path.join(dest,x)
    folders =map(append_dest,database_folders)
    map(create_folder,folders)

    files = raw['path']
    for f in files:
        create_folder(os.path.dirname(os.path.join(dest,f)))
        shutil.copy(os.path.join(source,f),os.path.join(dest,f))


def copy_coded(source,dest,lookup):
    append_dest = lambda x: os.path.join(dest, x)
    folders = map(append_dest, database_folders)
    map(create_folder, folders)


    def _copy(item):
        create_folder(os.path.dirname(os.path.join(dest, item['file'])))
        shutil.copy(os.path.join(source, item['file']),
                    os.path.join(os.path.dirname(os.path.join(dest, item['file'])),item['code']+'.pdb'))




def generate_database_index(database_path):
    crystal_folder = os.path.join(database_path,'crystal')
    receptor_folder = os.path.join(database_path,'receptor')
    ligands_folder = os.path.join(database_path,'ligands')

    raw = []
    for dirpath,dirname,filenames in os.walk(ligands_folder):
        for filename in filenames:
            row = []
            # ID abcd_123
            row.append('_'.join(filename.split('_')[:2]))
            # PDBname abcd
            row.append(filename.split('_')[0])
            # path
            row.append(os.path.join(ligands_folder_name,filename.split('_')[0],filename))
            # innder_index
            row.append(int(filename.split('.')[0].split('_')[-1]))

    database = pd.DataFrame(data = raw,columns = ['ID','PDBname','path','inner_index'])

    dump_dataframe(database, 'database_index')
    #return database




def code_test_set(test_set):
    wholedata =pd.read_csv("/home/xl198/remark/fast.csv")
    sorted = wholedata.merge(test_set,on=['ID','inner_index'])
    gou = sorted.groupby("ID")
    groups = gou.groups
    raw_answer_list = []
    for key in groups.key():
        PDBname = key.split('_')[0]
        row = []
        row.append(os.path.join(receptor_folder_name,PDBname+'.pdb'))
        row.append(os.path.join(crystal_folder_name,PDBname,key+'.pdb'))
        for value in groups[key]:
            row.append(sorted.ix[value]['path'])
        raw_answer_list.append(row)

    total=sum(map(len,raw_answer_list))
    reindex = range(total)
    random.shuffle(reindex)

    coded_answer_list = []
    cur = 0
    for raw_row in raw_answer_list:
        coded_row = []
        for i in range(len(raw_row)):
            coded_row.append('L'+raw_row[i] if i else 'P'+raw_row[i])

    dump_answer_list(raw_answer_list,'raw_solution')
    dump_answer_list(coded_answer_list,'solution')

    lookup = []
    for i in range(len(coded_answer_list)):
        for j in range(len(coded_answer_list[i])):
            lookup.append([coded_answer_list[i][j],raw_answer_list[i][j]])

    lookup_tabel = pd.DataFrame(data = lookup,columns = ['code','file'])

    dump_dataframe(lookup_tabel,'lookup')






def divide_by_PDBname(database,testset_ratio):
    '''
    select train set and test set by testset_ratio
    :param database: DataFrame
    :param testset_ratio:
    :return:
    '''

    pdbs = list(set(database['PDBname']))
    train_pdb = []
    test_pdb = []
    for pdb in pdbs:
        if random.random()<testset_ratio:
            test_pdb.append(pdb)
        else:
            train_pdb.append(pdb)

    train_set = database.merge(pd.DataFrame(data = train_pdb,columns=['PDBname']))
    test_set = database.merge(pd.DataFrame(data = test_pdb,columns=['PDBname']))

    dump_dataframe(train_set, 'train_set')
    dump_dataframe(test_set, 'test_set')




