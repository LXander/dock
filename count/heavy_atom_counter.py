import os,sys
import prody
from glob import glob


def counter():
    ligand_list = glob(os.path.join('/n/scratch2/xl198/dude/data/dude/pdbs/ligands','*','*[decoys]*.pdb'))
    count = 0
    with open('atom_count.csv','a') as fout:
        with open('err.log','a') as errout:
            for ligand in ligand_list:
                try:
                    count +=1
                    
                    parsed = prody.parsePDB(ligand)
                    atom_num = parsed.select('not hydrogen').numAtoms()
                    fout.write('{},{}\n'.format(os.path.basename(ligand),atom_num))
                    print count
                except Exception as e:
                    errout.write(e)

counter()
