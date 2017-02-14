import os,sys
import pandas as pd
import numpy as np

class newworkScore:
    basePath = ''

    def __init__(self, folderName):
        pass

class sminaScore:
    basePath = ''


    def __init__(self,folderName):
        self.outputFolder = os.path.join(self.basePath,folderName,'sminaScore')
        self.tempFolder = os.path.join(self.basePath,folderName,'temp')


