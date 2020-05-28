# -*- coding: utf-8 -*-
"""
Created on Tue May 19 14:22:00 2020

@author: Guillaume
"""

import sys

for arg in sys.argv:
    print(arg)
    
print(sys.argv[0])

file_SO = sys.argv[0]
pocket_holo = sys.argv[1]
mat_SO = []
with open(file_SO,"r") as filin:
    for ligne in filin:
        ligne = ligne.strip('\n')
        mat_SO.append(ligne)

mat_SO = [i.split('\t') for i in mat_SO]

place_pocket = -1
for i in range(len(mat_SO[0])):
    if(pocket_holo == mat_SO[0][i]):
        place_pocket = i

if(place_pocket == -1):
    print("----ERROR----")
    print(pocket_holo)
    print("-------------")

place_pocket = place_pocket + 1 #decalage attention
max_SO = -1
place_fpocket = -1
for i in range(1,len(mat_SO)):
    if(float(mat_SO[i][place_pocket]) > max_SO):
        max_SO = float(mat_SO[i][place_pocket])
        place_fpocket = mat_SO[i][0]

print(place_fpocket)












