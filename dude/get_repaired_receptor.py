import os,sys
import json

json_file = '/home/xl198/remark/res.json'
repair_path = '/n/scratch2/xl198/dude/data/receptor/repair'
repaired_path = '/n/scratch2/xl198/dude/data/receptor/repaired'
def get_dic():
    dic= {}
    receptor_json = json.load(open(json_file))
    for item in receptor_json:
        dic[item['receptor']] = item['target'].lower()

    return dic

def move(dic):

    receptor_list = map(lambda x:os.path.join(repair_path,x),os.listdir(repair_path))
    for receptor in receptor_list:
        if len(open(receptor).readlines()):
            receptor_name = os.path.basename(receptor).split('.')[0]
            target_path = os.path.join(repaired_path,dic[receptor_name]+'.pdb')
            os.system('cp {} {}'.format(receptor,target_path))

move(get_dic())
