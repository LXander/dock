import os
import sys
import re


def path_iterator(sourcePath):
    for dirpath,dirname,filenames in os.walk(sourcePath):
        for filename in filenames:
            yield os.path.join(dirpath,filename)

def source_dest_iterator(sourcePath,map_func):
    for sourceFilePath in path_iterator(sourcePath):
        destFilePath = map_func(sourceFilePath)
        yield sourceFilePath,destFilePath