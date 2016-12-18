import config
import os,sys
#import subprocess
import time

'''
    This file is used to get the result form Yi
    and then organized them into folder
    and insert newline into it otherwise it can't be read into obabel

'''

def convert(raw_file_path):
    print "convert"
    file_name = os.path.basename(raw_file_path)
    receptor_name = file_name.split('_')[0]
    output_path = os.path.join(config.BASE_CONVERT2PDB, receptor_name)
    print output_path
    print os.path.exists(output_path)
    if not os.path.exists(output_path):
        os.mkdir(output_path)

    input_file_path = os.path.join(config.BASE_CONVERT,receptor_name,file_name)
    output_file_path = os.path.join(output_path,file_name.split('.')[0] + '.pdb')
    os.system('obabel -i mol2 %s -o pdb -O %s'%(input_file_path,output_file_path))


def run(input_file_path):
    
    file_name = os.path.basename(input_file_path)
    receptor_name = file_name.split('_')[0]
    output_path = os.path.join(config.BASE_CONVERT,receptor_name)
    if not os.path.exists(output_path):
        os.mkdir(output_path)

    output_file_path = os.path.join(output_path,file_name)
    with open(input_file_path) as infile:
        with open(output_file_path,'w') as outfile:
            for line in infile:
                outfile.write(line)
                if line == '@<TRIPOS>MOLECULE\n':
                    outfile.write('\n')

def get_all(num = None):
    base = config.BASE_YI
    files = os.listdir(base)
    size = num if num != None else len(files)
    sys.stderr.write("Convert %s files\n"%size)
    sys.stderr.write("first file is %s\n"%(os.path.join(base,files[0])))
    for i in range(size):
        run(os.path.join(base,files[i]))
        sys.stderr.write("write %d/%d\n"%(i+1,size))

def run_convert(base,offset):
    base_path = config.BASE_YI
    files = os.listdir(base_path)
    #print base
    #print offset
    index = base*1000+offset
  
    if len(files)>index:
        run(os.path.join(base_path,files[index]))
        convert(os.path.join(base_path,files[index]))
        
        

def main():
    args = sys.argv
    base  = int(args[1])
    offset = int(args[2])
    
    run_convert(base,offset)
    sys.stderr.write("run convert %s"%(base*1000+offset))
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
