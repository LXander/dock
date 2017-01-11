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
                yield os.path.join(dirpath,filename)

    def hierarchy_cluster(self,filePath,clusterNumber = 10):
        fileName = os.path.basename(filePath)
        ligand = '_'.join(fileName.split('_')[:3])
        receptor = fileName.split('_')[0]



        traj = md.load(filePath)
        distances = np.empty((traj.n_frames,traj.n_frames))
        for i in range(traj.n_frames):
            distances[i] = md.rmsd(traj,traj,i)

        reduced_distances = squareform(distances,checks=False)

        Z = scipy.cluster.hierarchy.linkage(reduced_distances ,methods='average')

        clusterResule = fcluster(Z,clusterNumber,criterion='maxcluster')


        ligands = np.array([ligand]*traj.n_frames).reshape((-1,1))
        receptors = np.array([receptor]*traj.n_frames).reshape((-1,1))
        inner_index = np.arange(1,traj.n_frames+1)
        labeledResult = np.vstack((ligands,receptors,inner_index,clusterResule)).transpose()


        return labeledResult



    def test(self,filePath):
        filePathes = self.get_file_path(filePath)

        clusterResult = self.hierarchy_cluster(filePathes.next())

        print(clusterResult)

    def clusterAll(self,filePath,savePath):
        allClusterResult = reduce(lambda x,y:np.concatenate((x,y)),self.get_file_path(filePath))

        df = pd.DataFrame(data = allClusterResult,columns = ['ligand','recult'])

        df.to_csv(savePath,index=False)

if __name__ == '__main__':
    c = cluster()
    c.test("")