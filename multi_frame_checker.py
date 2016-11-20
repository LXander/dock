import config
import os,sys
import re
import commands

output_path = '/n/scratch2/xl198/data/multi'


def check(file_path):
    filename = os.path.basename(file_path)

    system_convert_ligand = 'obabel -i pdb %s -h -o pdb -O %s' % (file_path, os.path.join('/tmp/', filename))
    std_out = str(commands.getstatusoutput(system_convert_ligand))
    if not std_out == "(0, \'1 molecule converted\')":
        str(commands.getstatusoutut("cp " + file_path + " " + output_path))
        print std_out
    else:
	print "cp " + file_path + " " + output_path

def run(base, offset):
    input_path = '/n/scratch2/xl198/data/H/data'
    # JobId start at 1
    cur = base * 1000 + offset - 1

    file_path = None

    for dirname, dirnames, filenames in os.walk(input_path):
        if len(filenames) > cur:
            if cur >= 0:
                file_path = os.path.join(dirname, filenames[cur])
            break

        else:
            cur -= len(filenames)

    # print 'file_path '+file_path


    if file_path:
        check(file_path)

def main():
    args = sys.argv
    base = int(args[1])
    offset = int(args[2])

    run(base, offset)
    #sys.stderr.write("run convert %s" % (base * 1000 + offset))
    '''
    if len(args)<2:
        num = None
    else:
        num = int(args[1])

    print num
    get_all(num)
    '''


if __name__ == '__main__':
    main()
