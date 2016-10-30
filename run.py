import getopt, sys
import os
import subprocess
import config
import util

def docking(script_index,index,debug_flag):
    script_file = config.COMMANDS_FILE[script_index]
    script_file = util.path_check(script_file)


    with  open(script_file) as fr:
        # read one line
        # jobindex start as 1 but list start at 0

        command = fr.readlines()[index-1]

        util.debug(debug_flag,"Running command %s",(command))

        subprocess.call(command, shell=True)

def main():
    # eg : python run.py 1 2 -d
    try:
        opts, args = getopt.getopt(sys.argv[1:], "d", ["help", "output="])
    except getopt.GetoptError as err:
        # print help information and exit:
        print str(err)  # will print something like "option -a not recognized"

        sys.exit(2)


    script_index = int(args[0])
    index = int(args[1])

    #if -d display debug information
    debug_flag = False
    for o,a in opts:
        if o=='-d':
            debug_flag = True


    docking(script_index, index, debug_flag)

if __name__ == '__main__':
    main()