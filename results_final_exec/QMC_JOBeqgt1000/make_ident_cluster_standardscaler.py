#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 17:11:16 2019

@author: medina
"""

folders = ["1yq2", "DEBAIXO_1f4h", "DEBAIXO_6abu", "DEBAIXO_6ji0", "DECIMA_1f4h", "DECIMA_6abu", "DECIMA_6ji0"] 

list_files = ["AffinityPropagation_Data.csv", "DBSCAN_Data.csv", "KMeans_Data.csv", "MeanShift_Data.csv", "SpectralClustering_Data.csv", "Ward_Data.csv"]

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df_ident_cluster = pd.DataFrame(columns=['prot', 'identidade', 'AffinityPropagation_Data', 'DBSCAN_Data', 'KMeans_Data', 'MeanShift_Data', 'SpectralClustering_Data', 'Ward_Data'])

file_bas = './artigo.csv'
d2 = pd.read_csv(file_bas)

NORM_METHOD = 'StandardScaler'
cluster_date_dir = 'Clusters_Data_'+NORM_METHOD

for f in folders:
    
    print(f)
    print(str(list(d2[d2['ID']==f]['Identity'])[0]))
    identity = int(str(list(d2[d2['ID']==f]['Identity'])[0]).split('%')[0].split('(')[-1])
    
    d_aux = pd.DataFrame(columns=['prot', 'identidade', 'AffinityPropagation_Data', 'DBSCAN_Data', 'KMeans_Data', 'MeanShift_Data', 'SpectralClustering_Data', 'Ward_Data'])
    d_aux['prot'] = [f]
    d_aux['identidade'] = [identity]
    
    for clt in list_files:
        
        file_cls = './'+f+'/'+cluster_date_dir+'/'+clt
        print("reading file: " + str(file_cls))
        
        meth = clt.split('.')[0]
        d = pd.read_csv(file_cls)
        n_meth = list(d["Cluster"])[-1]
        
        d_aux[meth] = [int(n_meth) + 1]
        
    df_ident_cluster = df_ident_cluster.append(d_aux, ignore_index=True)
    
    
    
    

X = np.array(df_ident_cluster).T

plt.figure(figsize=(6 * 2 + 3, 9.5))
plt.subplots_adjust(left=.02, right=.98, bottom=.1, top=.96, wspace=.2, hspace=.2)
plot_num = 1

max_n = max(X[2:].reshape(-1)) + 1
min_n = min(X[2:].reshape(-1)) - 1

for clt in range(len(list_files)):
        
        # plot
        plt.subplot(2, 3, plot_num)
        plt.ylim(min_n, max_n)
        plt.xlabel("Percentual de Identidade (%)")
        plt.ylabel("Número de agrupamentos")
        plt.title(list_files[clt].split('.')[0], size=18)
        plt.scatter(X[1], X[2+clt])

        plot_num += 1

plt.savefig(NORM_METHOD+'_dispersion_graph_identity_n-clusters.png')
