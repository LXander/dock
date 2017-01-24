import numpy as np
import time,os,re,sys
import prody
from av4_atomdict import atom_dictionary
import getopt



sys.path.append("..")
from util.createfolder import create_chain_parent_folder,create_chain_folder,try_create_chain_folder

class FLAGS:
    orchestra_arrayjob = False

class stats:
    start = time.time()
    ligands_parsed = 0
    ligands_failed = 0

# crawls the database with crystal ligands
# looks which of them have docked ligands - if can acquire (?) and have pdb(?) - if can acquire
# moves all of them in a single folder with protein together
# writes statistics
# TODO make exceptions more informative

def convert_database_to_av4(database_path,positives_folder=None,decoys_folder=None,receptors_folder=None):
    """Crawls the folder (receptors in this case) and saves every PDB it finds
    into .npy array with 1) coordinates 2) mapped to the atom name number """

    # make a directory where the av4 form of the output will be written
    output_path = str(database_path+'_av4')
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    def save_av4(filepath,labels,elements,multiframe_coords):
        labels = np.asarray(labels,dtype=np.int32)
        elements = np.asarray(elements,dtype=np.int32)
        multiframe_coords = np.asarray(multiframe_coords,dtype=np.float32)

        if not (int(len(multiframe_coords[:,0]) == int(len(elements)))):
            raise Exception('Number of atom elements is not equal to the number of coordinates')

        if multiframe_coords.ndim==2:
            if not int(len(labels))==1:
                raise Exception ('Number labels is not equal to the number of coordinate frames')
        else:
            if not (int(len(multiframe_coords[0, 0, :]) == int(len(labels)))):
                raise Exception('Number labels is not equal to the number of coordinate frames')

        number_of_examples = np.array([len(labels)], dtype=np.int32)
        av4_record = number_of_examples.tobytes()
        av4_record += labels.tobytes()
        av4_record += elements.tobytes()
        av4_record += multiframe_coords.tobytes()
        f = open(filepath + ".av4", 'w')
        f.write(av4_record)
        f.close()


    count = 0
    database_ligand_path = os.path.join(database_path,'ligands')
    database_receptor_path = os.path.join(database_path,'receptors')
    for receptor in os.listdir(database_ligand_path):
        for ligand_name in os.listdir(os.path.join(database_ligand_path,receptor)):
            count +=1
            destFile = os.path.join(output_path,receptor, ligand_name+".av4")
            if os.path.exists(destFile):
                continue
            if FLAGS.orchestra_arrayjob and FLAGS.orchestra_jobid%FLAGS.orchestra_jobsize == count%FLAGS.orchestra_jobsize:
                continue
            ligand_folder = os.path.join(database_ligand_path,receptor,ligand_name)

            splited_ligands = os.listdir(ligand_folder)

            if len(splited_ligands) == 0:
                with open(os.path.join(database_path,'empty.txt'),'a') as fout:
                    fout.write(ligand_folder+'\n')
                continue

            path_to_receptor = os.path.join(database_receptor_path,receptor+'.pdb')
            path_to_first_ligand = os.path.join(ligand_folder,splited_ligands[0])

            try:
                prody_receptor = prody.parsePDB(path_to_receptor)
                prody_first_ligand = prody.parsePDB(path_to_first_ligand)

                multiframe_ligand_coords = prody_first_ligand.getCoords()
                labels = np.array([0])


                # if have more than one ligands, write them as one multiframe ligand
                if len(splited_ligands)>1:
                    for rest_ligand in splited_ligands[1:]:
                        prody_rest = prody.parsePDB(os.path.join(ligand_folder,rest_ligand))

                        # see if decoy is same as the initial ligand
                        if not all(np.asarray(prody_rest.getElements()) == np.asarray(prody_first_ligand.getElements())):
                            raise Exception('attempting to add ligand with different order of atoms')

                        multiframe_ligand_coords = np.dstack((multiframe_ligand_coords, prody_rest.getCoords()))

                        labels = np.concatenate((labels, [0]))

            except Exception as e:
                print e
                stats.ligands_failed += 1
                print "ligands parsed:", stats.ligands_parsed, "ligands failed:", stats.ligands_failed
                continue

            stats.ligands_parsed += 1
            print "ligands parsed:", stats.ligands_parsed, "ligands failed:", stats.ligands_failed

            # create an output path to write binaries for protein and ligands
            path_to_pdb_subfolder = os.path.join(output_path,receptor)

            try_create_chain_folder(path_to_pdb_subfolder)



            # convert atomnames to tags and write the data to disk
            def atom_to_number(atomname):
                atomic_tag_number = atom_dictionary.ATM[atomname.lower()]
                return atomic_tag_number

            print prody_receptor.getElements()

            receptor_elements = map(atom_to_number, prody_receptor.getElements())
            ligand_elements = map(atom_to_number, prody_first_ligand.getElements())


            receptor_output_path = os.path.join(path_to_pdb_subfolder,receptor)
            save_av4(receptor_output_path, [0], receptor_elements, prody_receptor.getCoords())
            ligand_output_path = os.path.join(path_to_pdb_subfolder,ligand_name)
            save_av4(ligand_output_path, labels, ligand_elements, multiframe_ligand_coords)


def parse_FLAG():
    try:
        opts,args = getopt.getopt(sys.argv[1:],None,["jobsize=","jobid=","cores="])
    except getopt.GetoptError as err:
        # print help information and exit:
        print str(err)  # will print something like "option -a not recognized"

        sys.exit(2)

    for name,value in opts:
        if name == '--jobsize':
            FLAGS.orchestra_jobsize = int(value)
            print "--jobsize ",value
        if name == '--jobid':
            FLAGS.orchestra_jobid = int(value)
            print "--jobid", value
        if name == '--cores':
            FLAGS.cores = int(value)

    if hasattr(FLAGS,"orchestra_jobsize") and hasattr(FLAGS,"orchestra_jobid"):
        FLAGS.orchestra_arrayjob = True

    print "orchestra job ",FLAGS.orchestra_arrayjob

    if hasattr(FLAGS,'cores'):
        print "cores num ",FLAGS.cores

if __name__ == '__main__':
    parse_FLAG()
    convert_database_to_av4(database_path="/home/ubuntu/xiao/data/newkaggle/dude/test_set/test_set",positives_folder="crystal_ligands",decoys_folder="docked_ligands",receptors_folder='receptors')


