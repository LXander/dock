import os,sys

def path_check(path):
    if not os.path.exists(path):
        print "%s not exists!"%path
        exit(1)
    else:
        return path

def debug(debug_flag,output,content):
    if debug_flag:
        print output%content
