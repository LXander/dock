import os,sys
import re

mol2_path = '/n/scratch2/xl198/data/result'
pdb_path  = '/n/scratch2/xl198/data/pdbs'

if not os.path.exists(pdb_path):
    os.mkdir(pdb_path)

walk = os.walk(mol2_path)

for dirpath,dirnames,filenames in walk:
    for filename in filenames:
        
        source = os.path.join(dirpath,filename)
        temp = re.sub('/result/','/pdbs/',source)
        dest = re.sub('.mol','.pdb',temp)

        destdir = os.path.dirname(dest)
        if not os.path.exists(destdir):
            os.mkdir(destdir)
        #print source
        #print dest
        os.system("obabel -imol2 {} -opdb -O {}".format(source,dest))
	
