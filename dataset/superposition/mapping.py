import os,sys
import pandas as pd
import random

def mapping():
    cry = pd.read_csv(FLAGS.crystalFilePath)
    dec = pd.read_csv(FLAGS.decoyFilePath)

    cry_gou = cry.groupby("atom_num")
    dec_gou = dec.groupby("atom_num")
    dec_gous = dec_gou.groups
    records = []

    for name,cry_group in cry_gou:
        print "atom num ",name
        if not name in dec_gous.keys():
            sys.stderr.write("atom num {} doesn't have corresponding decoys.".format(name))
            continue
        else:
            dec_index = dec_gous[name]
            dec_group = dec.iloc[dec_index]
            
            dec_len = len(dec_group)
            cry_len = len(cry_group)
            if dec_len<cry_len:
                original_index = range(dec_len)*(cry_len/dec_len+1)
                sys.stderr.write("atom num {} crystal ligand num : {} , decoys liangd num : {} \n".format(name,cry_len,dec_len))
            else:
                original_index = range(dec_len)

            random.shuffle(original_index)

            selected_index = original_index[:cry_len]

            cry_files = list(cry_group['filename'])
            dec_files = list(dec_group.iloc[selected_index]['filename'])
            pairs = zip(cry_files,dec_files)
            for pair in pairs:
                records.append(pair)

    df = pd.DataFrame(records,columns=['ligand_box','ligand'])
    df['ligand_box'] = df['ligand_box'].apply(lambda x:x.strip())
    df['ligand'] = df['ligand'].apply(lambda x:x.strip())
    df['receptor'] = df['ligand_box'].apply(lambda x:x.replace("addH",'data').split('_')[0]+'.pdb')
    df.to_csv(os.path.join(FLAGS.formsFolder,'docking_pair.csv'),index=False)






class FLAGS:
    crystalFilePath = '/n/scratch2/xl198/data/Superposition/forms/crystal_small.csv'
    decoyFilePath = '/n/scratch2/xl198/data/Superposition/forms/decoys_atom_num.csv'
    formsFolder = '/n/scratch2/xl198/data/Superposition/forms'

mapping()
