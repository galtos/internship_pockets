# -*- coding: utf-8 -*-
"""
Created on Wed Jun  3 09:26:05 2020

@author: Guillaume
"""
import os
import pandas as pd
import numpy as np
from math import *
import random
import time
#
Intercept = 6.61543
hydrophobic_kyte = -0.52657#vec1[56]-vec2[56]
p_aromatic_residues = -0.26214#vec1[69]-vec2[69]
p_charged_residues = -0.42272#vec1[6]-vec2[6]
p_negative_residues = 0.33210#vec1[23]-vec2[23]
p_positive_residues = 0.59123#vec1[10]-vec2[10]
charge = -0.73920#vec1[24]-vec2[24]
p_Nlys_atom = -0.26117#vec1[15]-vec2[15]
p_Ocoo_atom = -0.16479#vec1[67]-vec2[67]
p_small_residues = -0.16480#vec1[7]-vec2[7]
p_C_atom = 0.19381#vec1[64]-vec2[64]
p_nitrogen_atom = -0.36986#vec1[9]-vec2[9]
RADIUS_HULL = -3.33198#vec1[31]-vec2[31]
SURFACE_HULL = -3.57145#vec1[73]-vec2[73] 
DIAMETER_HULL = 1.10726#vec1[49]-vec2[49]
VOLUME_HULL = 2.46588#vec1[28]-vec2[28]
RADIUS_CYLINDER = 1.07161#vec1[3]-vec2[3]
C_RESIDUES = -0.65951#vec1[63]-vec2[63]
PSI = -0.18916#vec1[58]-vec2[58]
PCI = 0.16018#vec1[12]-vec2[12]
INERTIA_2 = -0.32137#vec1[18]-vec2[18]
INERTIA_3 = -0.15445#vec1[70]-vec2[70]
A = -0.16218#vec1[32]-vec2[32]
C = -0.34735#vec1[34]-vec2[34]
D = -0.15076#vec1[36]-vec2[36]
G = -0.43486#vec1[37]-vec2[37]
F = -0.18449#vec1[38]-vec2[38]
I = -0.30728#vec1[39]-vec2[39]
H = -0.27257#vec1[40]-vec2[40]
K = -0.27611#vec1[41]-vec2[41]
M = -0.27610#vec1[42]-vec2[42]
L = -0.34013#vec1[43]-vec2[43]
Q = -0.11223#vec1[45]-vec2[45]
P = -0.11730#vec1[46]-vec2[46]
S = -0.30554#vec1[47]-vec2[47]
R = -0.22685#vec1[48]-vec2[48]
T = -0.30320#vec1[50]-vec2[50]
W = -0.15142#vec1[51]-vec2[51]
Y = -0.11486#vec1[53]-vec2[53]

def descriptors_similarity(vec1, vec2):
    sim = Intercept + \
    hydrophobic_kyte * sqrt((vec1[56]-vec2[56])**2) + \
    p_aromatic_residues * sqrt((vec1[69]-vec2[69])**2) + \
    p_charged_residues * sqrt((vec1[6]-vec2[6])**2) + \
    p_negative_residues *sqrt((vec1[23]-vec2[23])**2) + \
    p_positive_residues *sqrt((vec1[10]-vec2[10])**2) + \
    charge *sqrt((vec1[24]-vec2[24])**2) + \
    p_Nlys_atom *sqrt((vec1[15]-vec2[15])**2) + \
    p_Ocoo_atom *sqrt((vec1[67]-vec2[67])**2) + \
    p_small_residues *sqrt((vec1[7]-vec2[7])**2) + \
    p_C_atom *sqrt((vec1[64]-vec2[64])**2) + \
    p_nitrogen_atom *sqrt((vec1[9]-vec2[9])**2) + \
    RADIUS_HULL *sqrt((vec1[31]-vec2[31])**2) + \
    SURFACE_HULL *sqrt((vec1[73]-vec2[73])**2) + \
    DIAMETER_HULL *sqrt((vec1[49]-vec2[49])**2) + \
    VOLUME_HULL *sqrt((vec1[28]-vec2[28])**2) + \
    RADIUS_CYLINDER *sqrt((vec1[3]-vec2[3])**2) + \
    C_RESIDUES *sqrt((vec1[63]-vec2[63])**2) + \
    PSI *sqrt((vec1[58]-vec2[58])**2) + \
    PCI *sqrt((vec1[12]-vec2[12])**2) + \
    INERTIA_2 *sqrt((vec1[18]-vec2[18])**2) + \
    INERTIA_3 *sqrt((vec1[70]-vec2[70])**2) + \
    A *sqrt((vec1[32]-vec2[32])**2) + \
    C *sqrt((vec1[34]-vec2[34])**2) + \
    D *sqrt((vec1[36]-vec2[36])**2) + \
    G *sqrt((vec1[37]-vec2[37])**2) + \
    F *sqrt((vec1[38]-vec2[38])**2) + \
    I *sqrt((vec1[39]-vec2[39])**2) + \
    H *sqrt((vec1[40]-vec2[40])**2) + \
    K *sqrt((vec1[41]-vec2[41])**2) + \
    M *sqrt((vec1[42]-vec2[42])**2) + \
    L *sqrt((vec1[43]-vec2[43])**2) + \
    Q *sqrt((vec1[45]-vec2[45])**2) + \
    P *sqrt((vec1[46]-vec2[46])**2) + \
    S *sqrt((vec1[47]-vec2[47])**2) + \
    R *sqrt((vec1[48]-vec2[48])**2) + \
    T *sqrt((vec1[50]-vec2[50])**2) + \
    W *sqrt((vec1[51]-vec2[51])**2) + \
    Y *sqrt((vec1[53]-vec2[53])**2)
    return 1-exp(sim)/(1+exp(sim))



path_data_ph="../data/dt_72clean-50.csv"
dt = pd.read_csv(path_data_ph, sep=',', index_col = 0)
#dt.index.values
names_dt = dt.index.values

Nseeds = 30
path_data_clusters = "../results/kmeans_results_reglog/kmeans_reglog_seeds"+str(Nseeds)+"_clusters.txt"
dt_clusters = pd.read_csv(path_data_clusters, sep=',', index_col = 0, header = None)

#
dt.loc[names_dt[0],:]

count=0
for i in range(1):
    print(i)
    name_cluster_i = []
    for cluster in range(len(dt_clusters)):
        if dt_clusters.values[cluster] == i:
            name_cluster_i.append(dt_clusters.index.values[cluster])
    mat_dist= pd.DataFrame(0, index = name_cluster_i, columns=name_cluster_i)
    for j in range(20):#len(name_cluster_i)):
        for k in range(j):
            mat_dist.loc[name_cluster_i[j],name_cluster_i[k]] = descriptors_similarity(dt.loc[name_cluster_i[j],:],dt.loc[name_cluster_i[k],:])
    mat_dist.to_csv("../results/kmeans_results_reglog/mat_dist_30/"+"mat_dist_"+str(i)+".csv")

    
    
    
    
    
    
    
    


