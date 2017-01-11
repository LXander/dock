from __future__ import print_function

import os,sys
import re

import mdtraj as md
import numpy as np
import pandas as pd

import scipy.cluster.hierarchy
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import fcluster

class cluster:


    def get_file_path(self,sourcePath):
        for dirpath,dirname,filenames in os.walk(sourcePath):
            for filename in filenames:
                if not len(open(os.path.join(dirpath,filename)).readlines()):
                    print("empty ",filename)
                    continue   
                yield os.path.join(dirpath,filename)

    def hierarchy_cluster(self,filePath,clusterNumber = 10):
        fileName = os.path.basename(filePath)
        ligand = '_'.join(fileName.split('_')[:3])
        receptor = fileName.split('_')[0]

   

        traj = md.load(filePath)
        print("frame number ",traj.n_frames)
        distances = np.empty((traj.n_frames,traj.n_frames))
        for i in range(traj.n_frames):
            distances[i] = md.rmsd(traj,traj,i)

        reduced_distances = squareform(distances,checks=False)

        Z = scipy.cluster.hierarchy.linkage(reduced_distances ,metric='average')

        clusterResult = fcluster(Z,clusterNumber,criterion='maxclust').reshape((-1,1))
        #print(clusterResult.shape)
        #print(clusterResult)
        ligands = np.array([ligand]*traj.n_frames).reshape((-1,1))
        receptors = np.array([receptor]*traj.n_frames).reshape((-1,1))
        inner_index = np.arange(1,traj.n_frames+1).reshape((-1,1))
        #print(inner_index)
        #print( "{} {} {} {}".format(ligands.shape[0],receptors.shape[0],inner_index.shape[0],clusterResult.shape[0]))
        labeledResult = np.hstack((ligands,receptors,inner_index,clusterResult))


        return labeledResult



    def test(self,filePath):
        print("test")
        filePathes = self.get_file_path(filePath)
        path = filePathes.next()
        print(path)    
        clusterResult = self.hierarchy_cluster(path)

        print(clusterResult)

    def clusterAll(self,filePath,savePath):
        allClusterResult = reduce(lambda x,y:np.concatenate((x,y)),self.get_file_path(filePath))

        df = pd.DataFrame(data = allClusterResult,columns = ['ligand','recult'])

        df.to_csv(savePath,index=False)

    def writeDownCluster(self,filePath,savePath):
        for path in self.get_file_path(filePath):
            pass        

if __name__ == '__main__':
    c = cluster()
    print("start")
    #c.test("/n/scratch2/xl198/dude/data/dude/docked_dude/docked_ligands")
    #c.test("/n/scratch2/xl198/dude/data/dude_400/pdbs/docked_ligands")
    
