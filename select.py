import os, sys
import config
#import prody
#import subprocess

'''
    select ligand which have a rmsd more than 6

    ## basic step

    1. convert the mol2 file to seperate pdb files
    2. statistic the number of pdb
    3. start at the last one calculate the rmsd
    4. record 20 of them which have rmsd > 6 and copy them to final folder

    ## note
    since the original pdb have different order
    first we should open it with smina and save it again
    after that smina will add Hydrogen to it
    so need to use obabel to delete them
    thus we get the native position to calculate rmsd
    1. generate the path of the original ligand
    2. open it with smina and save it to temp folder
    3. delete hydrogen with obabel



    ## Unfinished

'''


def run(source_file_path):
    assert os.path.exists(source_file_path)

    # split part of the path to generate all the path we use
    file_name = os.path.basename(source_file_path)
    temp_folder_name = file_name.split('.')[0]
    source_folder_path = os.path.dirname(source_file_path)
    folder_name = os.path.basename(source_folder_path)

    # prepare the native


    # generate all the path used as temp folder
    base = config.COMMANDS_SELECT['temp']
    temp_mid_folder = os.path.join(base,folder_name)
    if not os.path.exists(temp_mid_folder):
        os.mkdir(temp_mid_folder)
    temp_folder = os.path.join(temp_mid_folder,temp_folder_name)
    assert not os.path.exists(temp_folder)
    os.mkdir(temp_folder)

    #use obabel to comvert mol2 into seperate pdb
    subprocess.call("obabel -i mol2 -m %s -o pdb -O %s"%
                    (source_file_path,os.path.join(temp_folder,temp_folder_name+'_.pdb')),
                    shell=True)

    # get converted_file
    temp_files = os.listdir(temp_folder)
    temp_dict = {}
    for temp_file in temp_files:
        key = temp_file.split('.')[0].split('_')[-1]
        temp_dict[int(key)] = temp_file

    # calculate the rmsd from the last one
    keys = temp_dict.keys()
    keys.sort(reverse=True)






