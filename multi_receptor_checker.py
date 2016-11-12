import numpy as np
# import tensorflow as tf
import os, re, commands


def udf_write_log(filename, string):
    f = open(filename, "a")
    f.write("log mssg:")
    f.write(string)
    f.write("\n")
    f.close()


input_path = '/n/scratch2/xl198/data/H/data/'
output_path = '/n/scratch2/xl198/data/test/preprocessed_database/'
output_error_path = '/n/scratch2/xl198/data/test/ligands_with_preprocess_error'
output_receptor_path = '/n/scratch2/xl198/data/test/preprocessed_rectprot_database'
output_receptor_error_path = '/n/scratch2/xl198/data/test/receptor_with_prepeocess_error'


def check_ligand():
    for dirname, dirnames, filenames in os.walk(input_path):

        for filename in filenames:

            if re.search('_ligand.pdb$', filename):
                # if there is _ligand.pdb
                input_ligand_name = filename
                output_ligand_name = re.sub('_ligand.pdb', '_native_ligand.pdb', input_ligand_name)
                ligand_output_folder = output_path + (dirname.split("/")[-1])
                # make folder for the protein
                if not os.path.exists(ligand_output_folder):
                    os.mkdir(ligand_output_folder)
                    # put the protein inside that folder and check if conversion went as expected
                system_convert_ligand = 'obabel -i pdb ' + dirname + "/" + input_ligand_name + ' -h -o pdb -O ' + ligand_output_folder + "/" + output_ligand_name
                std_out = str(commands.getstatusoutput(system_convert_ligand))
                # if problems
                if not std_out == "(0, \'1 molecule converted\')":
                    # make ligands that converted with errors detectable
                    # if the folder ERROR_$liganad_output_folder does not exist, create one
                    print "WARNING unexpected output:", std_out
                    error_folder = "/".join(ligand_output_folder.split("/")[:-1]) + "/ERROR_" + str(
                        ligand_output_folder.split("/")[-1])
                    if not os.path.exists(error_folder):
                        str(commands.getstatusoutput("mv " + ligand_output_folder + " " + error_folder))
                        print "moved folder to ", error_folder, "--------------"
                    # if exists, put the ligand there
                    else:
                        print "moved ligand to ", error_folder, "--------------"
                        str(commands.getstatusoutput(
                            "mv " + ligand_output_folder + "/" + output_ligand_name + " " + error_folder + "/"))
                        # remove directory
                        os.rmdir(ligand_output_folder)

                    # finally, write a log file
                    udf_write_log(error_folder + "/" + re.sub('\.pdb$', '', output_ligand_name) + ".log", std_out)




def check_receptor():
    for dirname, dirnames, filenames in os.walk(input_path):

        for filename in filenames:

            if not re.search('_ligand.pdb$', filename):
                # if there is not _ligand.pdb
                input_ligand_name = filename
                output_ligand_name = re.sub('.pdb', '_native.pdb', input_ligand_name)
                ligand_output_folder = output_receptor_path + (dirname.split("/")[-1])
                # make folder for the protein
                if not os.path.exists(ligand_output_folder):
                    os.mkdir(ligand_output_folder)
                    # put the protein inside that folder and check if conversion went as expected
                system_convert_ligand = 'obabel -i pdb ' + dirname + "/" + input_ligand_name + ' -h -o pdb -O ' + ligand_output_folder + "/" + output_ligand_name
                std_out = str(commands.getstatusoutput(system_convert_ligand))
                # if problems
                if not std_out == "(0, \'1 molecule converted\')":
                    # make ligands that converted with errors detectable
                    # if the folder ERROR_$liganad_output_folder does not exist, create one
                    print "WARNING unexpected output:", std_out
                    error_folder = "/".join(ligand_output_folder.split("/")[:-1]) + "/ERROR_" + str(
                        ligand_output_folder.split("/")[-1])
                    if not os.path.exists(error_folder):
                        str(commands.getstatusoutput("mv " + ligand_output_folder + " " + error_folder))
                        print "moved folder to ", error_folder, "--------------"
                    # if exists, put the ligand there
                    else:
                        print "moved ligand to ", error_folder, "--------------"
                        str(commands.getstatusoutput(
                            "mv " + ligand_output_folder + "/" + output_ligand_name + " " + error_folder + "/"))
                        # remove directory
                        os.rmdir(ligand_output_folder)

                    # finally, write a log file
                    udf_write_log(error_folder + "/" + re.sub('\.pdb$', '', output_ligand_name) + ".log", std_out)


check_receptor()