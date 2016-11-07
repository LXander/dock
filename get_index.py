import config
import os,sys
import re

def index(save_to,ligands_path ):
    '''

    :param save_to: the path to save index.txt
    :param ligands_path: the path to get index
    :return:
    '''

    files_path = []
    relative_root = os.path.basename(ligands_path)
    with open(os.path.join(save_to,'index.txt'),'w') as result:
        for folders in os.listdir(ligands_path):
            for f in os.listdir(os.path.join(ligands_path,folders)):
                result.write(os.path.join(relative_root,folders,f)+'\n')


def main():
    save_to = os.path.join(config.BASE_DATA,'select')
    ligands_path= os.path.join(save_to,'ligands')
    index(save_to,ligands_path)

if __name__ == '__main__':
    main()
