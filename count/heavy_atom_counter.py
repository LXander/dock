import os,sys
import prody
from glob import glob


def counter():
    ligand_list = glob(os.path.join('/n/scratch2/xl198/dude/data/dude/pdbs/ligands','*','*[decoys]*.pdb'))
    with open('atom_count.csv','a') as fout:
        with open('err.log','a') as errout:
            for ligand in ligand_list:
                try:
                    parsed = prody.parsePDB(ligand)
                    count = parsed.select('not hydrogen').numAtoms()
                    fout.write('{},{}\n'.format(os.path.basename(ligand),count))
                except Exception as e:
                    errout.write(e)

counter()
