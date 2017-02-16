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

    thread_num = 2
    process_num = 1

    def __init__(self):
        self.parse()

        self.tempFolderPath =  tempfile.mkdtemp()

        self.basePath = os.path.join(self.workplace, 'blast')
        self.formsPath = os.path.join(self.basePath,'forms')
        self.err_log_file =os.path.join(self.basePath,'err.log')
        self.mergedPath = os.path.join(self.basePath,'merged_forms')

        try_create_chain_folder(self.tempFolderPath)
        try_create_chain_folder(self.formsPath)

    def err_log(self,message):
        with open(self.err_log_file,'a') as fout:
            fout.write(message+'\n')

    def blast_func(self,pdb_name):

        row_pdb = os.path.join('/n/scratch2/xl198/data/rcsb/row',pdb_name.upper()+'.pdb')


        handle = open(row_pdb,'rU')
        pairs = []
        for record in SeqIO.parse(handle,'pdb-seqres'):
            pairs.append([record.id,record.seq])

        for seq_id,seq in pairs:
            seq_size= len(seq)
            try:
                result_handle = NCBIWWW.qblast('blastp', 'nr', str(seq))
                blast_records = NCBIXML.read(result_handle)
                for ali in blast_records.alignments:
                    # [u'gi', u'209447557', u'pdb', u'3EML', u'A']
                    ali_info = ali.hit_id.split('|')
                    for hsp in ali.hsps:
                        records = [str(hsp.score),str(hsp.bits),str(hsp.expect),str(hsp.identities),str(hsp.positives)]
                        # [pdbname,seq_size,gi,id,type,name,chain]
                        full_record = [pdb_name,seq_id,str(seq_size)] + ali_info + records
                        self.mutex.acquire()
                        with open(os.path.join(self.formsPath,'blast_'+str(self.jobid)+'.txt'),'a') as fout:
                            fout.write(','.join(full_record)+'\n')
                        self.mutex.release()
            except Exception as e:
                self.err_log(pdb_name)
                print e



    def run_blast(self):
        receptor_list = os.listdir('/n/scratch2/xl198/data/rcsb/row_receptor')
        receptor = map(lambda x:x.split('.')[0],receptor_list)

        print len(receptor)

        converted_receptor = open('/n/scratch2/xl198/data/blast/merged_forms/merged_list.txt').readlines()
        converted_receptor = map(lambda x:x.strip(),converted_receptor)

        print len(converted_receptor)

        rest_receptor = list(set(receptor) - set(converted_receptor))

        print len(rest_receptor)

        self.convert(rest_receptor,self.blast_func)


    def merge_blast_result(self):

        blast_result = os.listdir(self.formsPath)
        blast_result_path = map(lambda x:os.path.join(self.formsPath,x),blast_result)

        # append result in to blast
        merged_blast = os.path.join(self.mergedPath,'blast.txt')
        for blast in blast_result_path:
            with open(merged_blast,'a') as fout:
                for line in open(blast):
                    fout.write(line)

        # get blasted receptor name
        blasted_receptors = []
        for line in open(merged_blast):
            blasted_receptors.append(line.split(',')[0])

        # write down the blasted receptor into file
        blasted_receptors_log  = os.path.join(self.mergedPath,'merged_list.txt')
        with open(blasted_receptors_log,'w') as fout:
            for item in blasted_receptors:
                fout.write(item+'\n')

        # move current forms folder, create a new empty folder
        count = 1
        moved_folder = os.path.join(self.basePath,'forms_'+str(count))
        while(os.path.exists(moved_folder)):
            count+=1
            moved_folder = os.path.join(self.basePath, 'forms_' + str(count))

        os.system('mv {} {}'.format(self.formsPath,moved_folder))
        os.system('mkdir -p {}'.format(self.formsPath))




if __name__ == '__main__':
    blast = Blast()
    #blast.merge_blast_result()
    blast.run_blast()
