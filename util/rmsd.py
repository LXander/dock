import os,sys
import prody
import numpy as np

'''
usage
python rmsd.py a.pdb b.pdb
'''

def eval_rmsd(ligand_a, ligand_b):
    try:
        parsed_a = prody.parsePDB(ligand_a)
        parsed_b = prody.parsePDB(ligand_b)

        if not np.all(parsed_a.getElements() == parsed_b.getElements()):
            raise Exception("{} {} doesn't have correspond atom.".format(ligand_a, ligand_b))

        diff = parsed_a.getCoords() - parsed_b.getCoords()
        rmsd = np.sqrt(np.mean(np.power(diff,2)))
        print rmsd
    except Exception as e:
        print e

if __name__ == '__main__':
    args = sys.argv
    if len(args)<3:
        print "not enough argument"
    else:
        eval_rmsd(args[1],args[2])
