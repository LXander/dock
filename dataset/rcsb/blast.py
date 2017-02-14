import os, sys
import re
import numpy as np
import pandas as pd
import getopt
import tempfile
import threading
import multiprocessing
from glob import glob

from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

sys.path.append(os.path.dirname(os.path.dirname(sys.path[0])))
from util.createfolder import try_create_chain_folder, try_create_chain_parent_folder
from util.orchestra import orchestra_job


class Blast(orchestra_job):
    mutex = threading.Lock()
    arrayjob = False
    workplace = '/n/scratch2/xl198/data'
    thread_num = 10
    process_num = 1

    def __init__(self):
        self.parse()

        self.tempFolderPath =  tempfile.mkdtemp()

        self.basePath = os.path.join(self.workplace, 'blast')
        self.formsPath = os.path.join(self.basePath,'forms')


        try_create_chain_folder(self.tempFolderPath)
        try_create_chain_folder(self.formsPath)

    def blast_func(self,pdb_file):
        pdb_name = os.path.basename(pdb_file).split('.')[0]
        row_pdb = os.path.join('/n/scratch2/xl198/data/rcsb/row',pdb_name.upper()+'.pdb')
        parsed_xml_name = os.path.basename(pdb_file).replace('.pdb','.xml')
        parsed_xml_path = os.path.join(self.tempFolderPath,parsed_xml_name)

        handle = open(row_pdb,'rU')
        pairs = []
        for record in SeqIO.parse(handle,'pdb-seqres'):
            pairs.append([record.id,record.seq])

        for seq_id,seq in pairs:
            seq_size= len(seq)
            result_handle = NCBIWWW.qblast('blastp', 'nr', str(seq))
            blast_records = NCBIXML.read(result_handle)
            for ali in blast_records.alignments:
                # [u'gi', u'209447557', u'pdb', u'3EML', u'A']
                ali_info = ali.hit_id.split('|')
                for hsp in ali.hsps:
                    records = [hsp.score,hsp.bits,hsp.expect,hsp.identities,hsp.positives]
                    # [pdbname,seq_size,gi,id,type,name,chain]
                    full_record = [pdb_name,seq_id,str(seq_size)] + ali_info + records
                    self.mutex.acquire()
                    with open(os.path.join(self.formsPath,'blast_'+str(self.jobid)+'.txt'),'a') as fout:
                        fout.write(','.join(full_record)+'\n')
                    self.mutex.releaes()

    def run_blast(self):
        receptor_list = glob(os.path.join('/n/scratch2/xl198/data/rcsb/row_receptor','*.pdb'))
        self.convert(receptor_list,self.blast_func)


if __name__ == '__main__':
    blast = Blast()
    blast.run_blast()
