import os,sys
import pandas as pd
import numpy as np
import tempfile
import re
from glob import glob
import mdtraj as md
import getopt
sys.path.append(os.path.dirname(os.path.dirname(sys.path[0])))
from util.createfolder import try_create_chain_parent_folder


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

    traj = md.load(file_path)
    coords = traj.xyz
    atom_type = [atom.name for atom in traj.topology.atoms]

    file_name = os.path.basename(file_path).split('.')[0]
    file_id = re.search('([a-zA-Z0-9]{4}_\d+)',file_name).group()
    file_dataset = re.search('ligand_(.+)$', file_name).groups()[0]
    if file_id == None or file_dataset == None:
        message = "cannot parse Id or dataset for {}".format(file_path)
        raise(message)

    select_by_id = FLAGS.dataframe[FLAGS.dataframe['ID'] == file_id]
    select_by_dataset = select_by_id[select_by_id['dataset'] == file_dataset ]
    affinity = np.array(select_by_dataset['smina_score'])

    #remarks = [line for line in open(file_path) if line[:6]=='REMARK' ]
    #affinity = [float(re.search(r'-?\d+.\d+',remark).group()) for remark in remarks]

    if len(affinity) != traj.n_frames:
        message ='Error when parse {}, fram num {}, affinity {}\n'.format(file_path,traj.n_frames,len(affinity))
        sys.stderr.write(message)
        raise Exception(message)

    else:
        return coords,atom_type,affinity

def parseLigand(fast_file_path):
    '''
    get the ligand's coordinate of whole protein fast docking, whole rigorous docking
    and site only rigorous docking
    :param fast_file_path:
    :return: coords,atom_type,affinity
    '''
    fast_name = os.path.basename(fast_file_path)
    receptor = fast_name.split('_')[0]
    original_name = fast_name.replace('ligand_fast','ligand')
    rigor_name= fast_name.replace('ligand_fast','ligand_rigor')
    rigor_so_name = fast_name.replace('ligand_fast','ligand_rigor_so')

    rigor_file_path = os.path.join(FLAGS.rigor_path,receptor,rigor_name)
    rigor_so_file_path= os.path.join(FLAGS.rigor_so_path,receptor,rigor_so_name)

    coords,atom_type,affinity= parsePDB(fast_file_path)

    if os.path.exists(rigor_file_path):
        rigor_coords,rigor_atom_type,rigor_affinity = parsePDB(rigor_file_path)
        if rigor_atom_type == atom_type:
            coords =np.concatenate((coords,rigor_coords))
            affinity = np.concatenate((affinity,rigor_affinity))
        else:
            message = 'different atom type for {} and {}\n'.format(fast_name,rigor_name)
            sys.stderr.write(message)
            raise Exception(message)

    if os.path.exists(rigor_so_file_path):
        rigor_so_coords,rigor_so_atom_type,rigor_so_affinity = parsePDB(rigor_so_file_path)
        if rigor_so_atom_type == atom_type:
            coords = np.concatenate((coords,rigor_so_coords))
            affinity = np.concatenate((affinity,rigor_so_affinity))
        else:
            message = 'different atom type for {} and {}\n'.format(fast_name,rigor_so_name)
            sys.stderr.write(message)
            raise Exception(message)


    return coords,atom_type,affinity


def overlap_filter(crystal_list, ligand_coords, ligand_affinity):

    for crystal_ligand_path in crystal_list:
        crystal_ligand = md.load(crystal_ligand_path)
        # the shape of crystal coords is [1,n,3]
        crysta_coords = crystal_ligand.xyz

        # [x,y,1,1,3]
        exp_ligand_coord = np.expand_dims(np.expand_dims(ligand_coords, -2), -2)
        # [x,y,1,n,3]
        diff = exp_ligand_coord - crysta_coords
        # [x,y,n]
        distance = np.squeeze(np.sqrt(np.sum(np.square(diff),-1)))
        # [x,y]
        atom_overlap =(np.sum((distance<FLAGS.clash_cutoff_A).astype(np.float32),-1)>0).astype(np.float32)
        # [x]
        ligand_not_overlap = np.mean(np.squeeze(atom_overlap),-1)<FLAGS.clash_size_cutoff

        ligand_coords = ligand_coords[ligand_not_overlap]
        ligand_affinity = ligand_affinity[ligand_not_overlap]

    sorted_index = np.argsort(ligand_affinity)
    ligand_coords = ligand_coords[sorted_index]
    ligand_affinity = ligand_affinity[sorted_index]

    return ligand_coords, ligand_affinity

def save_av4(filepath,labels,elements,multiframe_coords):
    concatenated_coords = multiframe_coords[0]
    for coords in multiframe_coords[1:]:
        concatenated_coords = np.concatenate((concatenated_coords,coords))

    labels = np.asarray(labels*100,dtype=np.int32)
    elements = np.asarray(elements,dtype=np.int32)
    multiframe_coords = np.asarray(multiframe_coords,dtype=np.float32)

    if not(int(len(multiframe_coords[:,0])==int(len(elements)))):
        raise Exception('Number of atom elements is not equal to the ')

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

def convert(fast_file_path):

    fast_name = os.path.basename(fast_file_path)
    mix_name = fast_name.replace('ligand_fast','ligand_mix')
    mix_name = mix_name.replace('.pdb','.av4')
    receptor = fast_name.split('_')[0]

    crystalFolder = os.path.join(FLAGS.crystalPath, receptor)
    similar_crystal = get_similar_crystal_file(crystalFolder,fast_file_path)

    ligand_coords,ligand_atom_type,ligand_afinity = parseLigand(fast_file_path)

    print 'parsed {} , shape'.format(fast_name)
    print ligand_coords.shape

    filtered_coords,filtered_affinity = overlap_filter(similar_crystal,ligand_coords,ligand_afinity)


    mix_dest_path = os.path.join(FLAGS.mix_path,receptor,mix_name)

    try_create_chain_parent_folder(mix_dest_path)

    save_av4(mix_dest_path,filtered_affinity,ligand_atom_type,filtered_coords)

def run():
    FLAGS.tempPath = tempfile.mkdtemp()

    fastFileList = glob(os.path.join(FLAGS.fast_path, '*', '*[fast]*.pdb'))

    index = range(len(fastFileList))

    if hasattr(FLAGS,'arrayjob') and FLAGS.arrayjob:
        if hasattr(FLAGS,'offset'):
            index = [len(index)/FLAGS.jobsize*FLAGS.offset+FLAGS.jobid]
        else:
            index = range(FLAGS.jobid-1,len(index),FLAGS.jobsize)

    for i in index:
        convert(fastFileList[i])

class FLAGS:
    fast_path = '/n/scratch2/xl198/data/pdbs'
    rigor_path = '/n/scratch2/xl198/YI/rigor/final'
    rigor_so_path = '/n/scratch2/xl198/YI/rigor_so/final'
    mix_path = '/n/scratch2/xl198/data/fusion/dock'
    crystalPath = '/n/scratch2/xl198/data/H/addH'
    dataframe =pd.read_csv("/n/scratch2/xl198/data/fusion/forms/simple_mix.csv")

    tanimoto_cutoff = 0.75
    clash_cutoff_A = 4
    clash_size_cutoff = 0.3

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
