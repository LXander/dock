import os,sys
import pubchempy as pcp

def query_similar(ligand_path):
    result = os.popen('obabel -ipdb {} -osmi'.format(ligand_path))
    smiles = result.read().split('\t')[0]
    cids = pcp.get_cids(smiles, 'smiles', searchtype = 'similarity', list_return = 'flat')

