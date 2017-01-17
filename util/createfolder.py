import os,sys,re


def create_chain_parent_folder(filePath):
    dirPath = os.path.dirname(filePath)
    create_chain_folder(dirPath)

def create_chain_folder(folderPath):
    dirPath = os.path.dirname(folderPath)
    if not os.path.exists(dirPath):
        create_chain_folder(dirPath)
    if not os.path.exists(folderPath):
        os.mkdir(folderPath)

def create_folder(folderPath):
    # If folder doesn't exists create it, if exists ignore
    if not os.path.exists(folderPath):
        os.mkdir(folderPath)

def create_parent_folder(filePath):
    dirPath = os.path.dirname(filePath)
    create_folder(dirPath)

