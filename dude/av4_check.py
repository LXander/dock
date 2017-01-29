import os,sys,time,re
from glob import glob
import tensorflow as tf


def decode_av4(serialized_record):
    # decode everything into int32
    tmp_decoded_record = tf.decode_raw(serialized_record, tf.int32)
    # first four bytes determine the number of frames
    number_of_frames = tf.slice(tmp_decoded_record, [0], [1])
    # labels are saved as in32 * number of frames in the record
    labels = tf.slice(tmp_decoded_record, [1], number_of_frames)
    # elements are saved as int32 and their number is == to the number of atoms
    number_of_atoms = ((tf.shape(tmp_decoded_record) - number_of_frames - 1) / (3 * number_of_frames + 1))
    elements = tf.slice(tmp_decoded_record, number_of_frames + 1, number_of_atoms)

    # coordinates are saved as a stack of X,Y,Z where the first(vertical) dimension
    # corresponds to the number of atoms
    # second (horizontal dimension) is x,y,z coordinate of every atom and is always 3
    # third (depth) dimension corresponds to the number of frames

    coords_shape = tf.concat(0, [number_of_atoms, [3], number_of_frames])
    tmp_coords = tf.slice(tmp_decoded_record, number_of_frames + number_of_atoms + 1,
                          tf.shape(tmp_decoded_record) - number_of_frames - number_of_atoms - 1)
    multiframe_coords = tf.bitcast(tf.reshape(tmp_coords, coords_shape), type=tf.float32)

    return labels, elements, multiframe_coords

def train():
    sess = tf.Session()
    #file_name = '/home/ubuntu/common/data/new_kaggle/train_small/labeled_av4/1j0k/1j0k_1180_ligand.av4'
    file_name = '/n/scratch2/xl198/dude/data/jan_22/actives_av4/aa2ar/aa2ar_actives_356.av4'
    serialized_ligand = tf.read_file(file_name)
    ligand_labels, ligand_elements, multiframe_ligand_coords = decode_av4(serialized_ligand)

    sess.run(tf.initialize_all_variables())


    labels = sess.run([ligand_labels])
    print labels

train()

