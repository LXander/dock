import pandas as pd
import numpy as np
import random
import os, sys, re
import shutil
from functools import partial

# folder_name used to store different data in different folder
crystal_folder_name = 'crystal_ligands'  # crystal liagnd
ligands_folder_name = 'docked_ligands'  # docking result
receptor_folder_name = 'receptors'  # receptor
suffix = '.npy'  # we may need to get different format of file '.pdb'\'.npy'

# for the convenient to create folders for train set at one time.
database_folders = [crystal_folder_name, ligands_folder_name, receptor_folder_name]

# for the convenient  to create folders for test set at one time.
test_set_folders = ['ligands', 'receptors']


def dump_dataframe(dataframe, filename):
    # store dataframe as csv file, without index
    dataframe.to_csv(filename + '.csv', index=False)


def dump_answer_list(answer_list, filename):
    # store answer list in to kaggle format: [Id][comma][Expected1][blankspace][Expected2]...
    with open(filename + '.txt', 'w') as fr:
        for row in answer_list:
            line = row[0] + ','
            line += ' '.join(row[1:])
            line += '\n'
            fr.write(line)


def create_folder(folder):
    # If folder doesn't exists create it, if exists ignore
    if not os.path.exists(folder):
        os.mkdir(folder)


def copy_original(source, dest, raw):
    '''
    copy file with original name form source to dest
    file sotre in raw

    source: path of source
    dest  : path of dest
    raw   : dataframe
    '''

    # create folders
    append_dest = lambda x: os.path.join(dest, x)
    folders = map(append_dest, database_folders)
    map(create_folder, folders)

    # copy docking result ligands
    files = raw['path']
    for f in files:
        create_folder(os.path.dirname(os.path.join(dest, f)))
        shutil.copy(os.path.join(source, f), os.path.join(dest, f))

    # copy crystal ligands
    crystal = raw.apply(
        lambda item: crystal_folder_name + '/' + item['PDBname'] + '/' + item['ID'] + '_ligand' + suffix, axis=1)

    crystal = list(set(crystal))

    for rec in crystal:
        create_folder(os.path.dirname(os.path.join(dest, rec)))
        shutil.copy(os.path.join(source, rec), os.path.join(dest, rec))

    # copy receptor
    receptor = raw.apply(lambda item: receptor_folder_name + '/' + item['PDBname'] + suffix, axis=1)
    receptor = list(set(receptor))

    for rec in receptor:
        create_folder(os.path.join(dest, 'receptor'))
        shutil.copy(os.path.join(source, rec), os.path.join(dest, rec))


def _copy(item, source, dest):
    # copy data from source path to dest path
    # put receptor in /receptors
    # put ligands in /ligands
    #
    # receptor's code begin with 'P' and ligands' code begin with 'L'
    #
    # print os.path.join(source, item['file'])
    create_folder(os.path.basename(os.path.join(dest,item['dest'])))
    shutil.copy(os.path.join(source, item['source']),
                os.path.join(dest, item['dest']))




def copy_coded(source, dest, copy_tabel):
    '''
    copy coded file from source to dest
    and changed file name from original to coded in lookup

    source  : path of source
    dest    : path of dest

    lookup  : dataframe [file,code]
    '''

    append_dest = lambda x: os.path.join(dest, x)
    folders = map(append_dest, test_set_folders)
    map(create_folder, folders)

    copy_function = partial(_copy, source=source, dest=dest)
    copy_tabel.apply(copy_function, axis=1)


def generate_database_index(database_path):
    crystal_folder = os.path.join(database_path, crystal_folder_name)
    receptor_folder = os.path.join(database_path, receptor_folder_name)
    ligands_folder = os.path.join(database_path, ligands_folder_name)

    raw = []
    for dirpath, dirname, filenames in os.walk(ligands_folder):
        for filename in filenames:
            row = []
            # ID abcd_123
            row.append('_'.join(filename.split('_')[:2]))
            # PDBname abcd
            row.append(filename.split('_')[0])
            # path
            row.append(os.path.join(ligands_folder_name, filename.split('_')[0], filename))
            # innder_index
            row.append(int(filename.split('.')[0].split('_')[-1]))
            raw.append(row)

    database = pd.DataFrame(data=raw, columns=['ID', 'PDBname', 'path', 'inner_index'])

    dump_dataframe(database, 'database_index')
    # return database


def code_test_set(test_set):
    wholedata = pd.read_csv("/home/xl198/remark/fast.csv")
    sorted = wholedata.merge(test_set, on=['ID', 'inner_index'])
    gou = sorted.groupby("ID")
    groups = gou.groups
    raw_answer_list = []
    for key in groups.keys():
        PDBname = key.split('_')[0]
        row = []
        row.append(os.path.join(receptor_folder_name, PDBname + suffix))
        row.append(os.path.join(crystal_folder_name, PDBname, key + '_ligand' + suffix))
        for value in groups[key]:
            row.append(sorted.ix[value]['path'])
        raw_answer_list.append(row)

    total = sum(map(len, raw_answer_list))
    reindex = range(total)
    random.shuffle(reindex)

    coded_answer_list = []
    cur = 0
    for raw_row in raw_answer_list:
        coded_row = []
        for i in range(len(raw_row)):
            coded_row.append('L' + str(reindex[cur]) if i else 'P' + str(reindex[cur]))
            cur += 1
        coded_answer_list.append(coded_row)

    dump_answer_list(raw_answer_list, 'raw_solution')
    dump_answer_list(coded_answer_list, 'solution')

    lookup = []
    for i in range(len(coded_answer_list)):
        for j in range(len(coded_answer_list[i])):
            lookup.append([raw_answer_list[i][j], coded_answer_list[i][j]])

    lookup_tabel = pd.DataFrame(data=lookup, columns=['code', 'path'])

    dump_dataframe(lookup_tabel, 'lookup')


    copy_tabel = []
    for i in range(len(coded_answer_list)):
        for j in range(len(coded_answer_list[i])):
            if j==0:
                copy_tabel.append([raw_answer_list[i][j],
                                   os.path.join(receptor_folder_name,coded_answer_list[i][j]+suffix)])
            else:
                copy_tabel.append([raw_answer_list[i][j],
                                   os.path.join(ligands_folder_name,coded_answer_list[i][0],coded_answer_list[i][j]+suffix)])


    coded_copy_tabel = pd.DataFrame(data = copy_tabel,columns = ['source','dest'])


def divide_by_PDBname(database, testset_ratio):
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
        if random.random() < testset_ratio:
            test_pdb.append(pdb)
        else:
            train_pdb.append(pdb)

    train_set = database.merge(pd.DataFrame(data=train_pdb, columns=['PDBname']))
    test_set = database.merge(pd.DataFrame(data=test_pdb, columns=['PDBname']))

    dump_dataframe(train_set, 'train_set')
    dump_dataframe(test_set, 'test_set')


