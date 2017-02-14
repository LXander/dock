import os,sys
import pandas as pd
import numpy as np
import tempfile
import re
from glob import glob
import mdtraj as md
import prody
import getopt

sys.path.append(os.path.dirname(os.path.dirname(sys.path[0])))
from util.createfolder import try_create_chain_parent_folder
from av4.av4_atomdict import atom_dictionary


def get_similar_crystal_file(crystalFolder,ligandPath):
    '''
    compare the ligand with all the crystal ligand for the target receptor
    return the path of the crystal ligand which have a high similarity

    :param crystalFolder: folder store crystal ligands
    :param ligandPath: path for docking result
    :return: list of path
    '''

    file_name= os.path.basename(ligandPath)
    temp_file_path = os.path.join(FLAGS.tempPath,file_name)
    print 'obabel convert ',file_name
    obabel_cmd = 'obabel -ipdb {} -l 1 -f 1 -opdb -O {}'.format(ligandPath,temp_file_path)
    os.system(obabel_cmd)


    crystalList = glob(os.path.join(crystalFolder,'*.pdb'))

    records = []

    for crystal_ligand in crystalList:
        command = os.popen('babel -d {} {} -ofpt -xfFP4'.format(crystal_ligand,temp_file_path))
        ls = command.read()
        try:
            tanimoto_similarity = re.split('=|\n', ls)[2]
            if tanimoto_similarity>FLAGS.tanimoto_cutoff:
                records.append(crystal_ligand)
        except:
            records.append(crystal_ligand)



    return records

def parsePDB(file_path):
    '''
    given path of the docking result
    :param file_path:
    :return:
    '''

    lig = prody.parsePDB(file_path)
    coords = lig.getCoordsets()
    elements = lig.getElements()



    remarks = [line for line in open(file_path) if line[:6] == 'REMARK']
    affinity = [float(re.search('(-?\d+.\d+)',reamrk).group()) for reamrk in remarks]

    if coords.shape[0] != len(affinity):
        message = "{} parse error, frame {}, affinity {}".format(file_path,coords.shape[0],len(affinity))
        raise(message)

    return coords,elements,affinity

def overlap_filter(crystal_list, ligand_coords, ligand_affinity):


    for crystal_ligand_path in crystal_list:
        crystal_ligand = prody.parsePDB(crystal_ligand_path)
        # the shape of crystal coords is [1,n,3]
        # coordinate read by mdtraj is 10 time less than
        # literally so we times it by 10

        # [n,3]
        crysta_coords = crystal_ligand.getCoords()

        # [x,y,1,1,3]
        exp_ligand_coord = np.expand_dims(np.expand_dims(ligand_coords, -2), -2)
        print "exp_ligand_coord ",exp_ligand_coord.shape
        # [x,y,1,n,3]
        diff = exp_ligand_coord - crysta_coords
        print "diff ",diff.shape
        # [x,y,n]
        distance = np.squeeze(np.sqrt(np.sum(np.square(diff),-1)))
        print "distance ",distance.shape
        # [x,y]
        atom_overlap =(np.sum((distance<FLAGS.clash_cutoff_A).astype(np.float32),-1)>0).astype(np.float32)
        print "atom_overlap ",atom_overlap.shape
        # [x]
        ligand_not_overlap = np.mean(np.squeeze(atom_overlap),-1)<FLAGS.clash_size_cutoff
        print "ligand_not_overlap num",np.sum(ligand_not_overlap)
        print crystal_ligand_path
        ligand_coords = ligand_coords[ligand_not_overlap]
        ligand_affinity = ligand_affinity[ligand_not_overlap]

    sorted_index = np.argsort(ligand_affinity)
    ligand_coords = ligand_coords[sorted_index]
    ligand_affinity = ligand_affinity[sorted_index]

    return ligand_coords, ligand_affinity

def save_av4(filepath,labels,elements,multiframe_coords):
    concatenated_coords = multiframe_coords[0]
    for coords in multiframe_coords[1:]:
        concatenated_coords = np.hstack((concatenated_coords,coords))

    labels = np.asarray(labels*100,dtype=np.int32)
    elements = np.asarray(elements,dtype=np.int32)
    concatenated_coords = np.asarray(concatenated_coords,dtype=np.float32)

    if not(int(len(concatenated_coords[:,0])==int(len(elements)))):
        raise Exception('Number of atom elements is not equal, elements num {}, coords size {}'.format(len(elements),len(concatenated_coords[:,0])))

    if not(int(len(concatenated_coords[0,:])==int(3*len(labels)))):
        raise Exception('Number labels is not equal to the number of coordinate frames')

    number_of_examples = np.array([len(labels)],dtype=np.int32)
    av4_record = number_of_examples.tobytes()
    av4_record += labels.tobytes()
    av4_record += elements.tobytes()
    av4_record += multiframe_coords.tobytes()
    f = open(filepath,'w')
    f.write(av4_record)
    f.close()


def convert(fast_path):
    file_name = os.path.basename(fast_path)

    dest_name = file_name.replace('.pdb', '.av4')
    file_id = re.search('(^[a-zA-Z0-9]{4}_\d+)', file_name).group()
    receptor = file_name.split('_')[0]

    crystalFolder = os.path.join(FLAGS.crystalPath, receptor)
    similar_crystal = get_similar_crystal_file(crystalFolder, fast_path)
    lig_coords,lig_elements,lig_affinity = parsePDB()

    filtered_coords, filtered_affinity = overlap_filter(similar_crystal,lig_coords,lig_affinity)

    av4_dest_path = os.path.join(FLAGS.super_dest_path,receptor,dest_name)

    try_create_chain_parent_folder(av4_dest_path)

    save_av4(av4_dest_path, filtered_affinity, lig_elements, filtered_coords)




def run():
    FLAGS.tempPath = tempfile.mkdtemp()

    fastFileList = glob(os.path.join(FLAGS.super_source, '*', '*.pdb'))

    index = range(len(fastFileList))

    if hasattr(FLAGS,'arrayjob') and FLAGS.arrayjob:
        if hasattr(FLAGS,'offset'):
            index = [len(index)/FLAGS.jobsize*FLAGS.offset+FLAGS.jobid]
        else:
            index = range(FLAGS.jobid-1,len(index),FLAGS.jobsize)

    for i in index:
        convert(fastFileList[i])

class FLAGS:
    super_source = '/n/scratch2/xl198/data/Superposition/dock'
    fast_path = '/n/scratch2/xl198/data/pdbs'
    rigor_path = '/n/scratch2/xl198/YI/rigor/final'
    rigor_so_path = '/n/scratch2/xl198/YI/rigor_so/final'
    mix_path = '/n/scratch2/xl198/data/fusion/dock'
    crystalPath = '/n/scratch2/xl198/data/H/addH'
    dataframe =pd.read_csv("/n/scratch2/xl198/data/fusion/forms/simple_mix.csv")
    super_dest_path = '/n/scratch2/xl198/data/Syperposition/dock_av4'
    tanimoto_cutoff = 0.75
    clash_cutoff_A = 4
    clash_size_cutoff = 0.9

def parse_FLAG():
    try:
        opts,args = getopt.getopt(sys.argv[1:],None,["offset=","jobsize=","jobid=","cores=","scan"])
    except getopt.GetoptError as err:
        # print help information and exit:
        print str(err)  # will print something like "option -a not recognized"

        sys.exit(2)

    for name,value in opts:
        if name == '--offset':
            FLAGS.joboffset = int(value)
        if name == '--jobsize':
            FLAGS.jobsize = int(value)
            print "--jobsize ",value
        if name == '--jobid':
            FLAGS.jobid = int(value)
            print "--jobid", value
        if name == '--cores':
            FLAGS.cores = int(value)
        if name == '--scan':
            FLAGS.scan = True

    if hasattr(FLAGS,"jobsize") and hasattr(FLAGS,"jobid"):
        FLAGS.arrayjob = True

    print "orchestra job ",FLAGS.arrayjob

    if hasattr(FLAGS,'cores'):
        print "cores num ",FLAGS.cores

if __name__ == '__main__':
    parse_FLAG()
    run()
