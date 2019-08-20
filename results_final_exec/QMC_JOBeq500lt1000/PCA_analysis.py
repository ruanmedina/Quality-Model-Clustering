# -*- coding: utf-8 -*-
"""
Created on Sun Aug 18 13:43:53 2019

@author: medina
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#df_ident_cluster = pd.DataFrame(columns=['prot', 'identidade', 'AffinityPropagation_Data', 'DBSCAN_Data', 'KMeans_Data', 'MeanShift_Data', 'SpectralClustering_Data', 'Ward_Data'])


dt = pd.read_csv("StandardScaler_pca_results.txt")

dt_comp_analy = pd.DataFrame(columns=['prot', 'Max_coor1', 'Max_coor2', 'Rep_var'])

variancia = {}

for l in range(len(dt)):
    
    mc1 = np.array(dt.columns)[1:7][np.array(abs(dt.iloc[l][1:7])).argmax()]
    mc2 = np.array(dt.columns)[7:13][np.array(abs(dt.iloc[l][7:13])).argmax()]
    
    variancia.update({dt['pbd'][l] : sum(dt.iloc[l][13:15])})
    
    d = pd.DataFrame(columns=['prot', 'Max_coor1', 'Max_coor2', 'Rep_var'])
    d['prot'] = [dt['pbd'][l]]
    d['Max_coor1'] = mc1.split('_')[0]    
    d['Max_coor2'] = mc2.split('_')[0]
    d['Rep_var'] = sum(dt.iloc[l][13:15])
    
    dt_comp_analy = dt_comp_analy.append(d, ignore_index=True)
    

import pickle
dt_comp_analy.to_pickle('QMC_JOBeq500it1000_freq_PCA_dataFrame.pickle')    


scores = ['molpdf', 'DOPE', 'DOPEHR', 'NDOPE', 'outlier', 'allowed']
N = len(scores)

freq_coor1 = list()
for s in scores:
    freq = 0 
    for val in dt_comp_analy['Max_coor1']:
        if val[1:] == s: freq = freq + 1
    print(freq)
    freq_coor1.append(freq)

freq_coor2 = list()
for s in scores:
    freq = 0 
    for val in dt_comp_analy['Max_coor2']:
        if val[1:] == s: freq = freq + 1
    freq_coor2.append(freq)
    
print(dt_comp_analy.set_index(['Max_coor1', 'Max_coor2', 'Rep_var']).count(level="Max_coor1"))
print(dt_comp_analy.set_index(['Max_coor1', 'Max_coor2', 'Rep_var']).count(level="Max_coor2"))

ind = np.arange(N)
width = 0.35

p1 = plt.bar(ind, freq_coor1, width)
p2 = plt.bar(ind, freq_coor2, width,
             bottom=freq_coor1)

plt.ylabel('Frequency of most significant on PCA')
plt.title('Frequency of atributes in the list of most significant on PCA')
plt.xticks(ind, scores)
#plt.yticks(np.arange(0, 81, 10))
plt.legend((p1[0], p2[0]), ('Coor1', 'Coor2'))


    