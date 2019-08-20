#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 17:11:16 2019

@author: medina
"""

#Correto
folders = ['153l_', '1a17_', '1a1x_', '1a2pA', '1a2zA', '1a3aA', '1a3c_', '1a3k_', '1a44_', '1a53_', '1a6f_', '1a6jA', '1a6m_', '1a73A', '1a79A', '1a7vA', '1a8l_', '1a9nA', '1aa7A', '1aac_', '1abv_', '1ac6A', '1acf_', '1acz_', '1ad2_', '1ae9A', '1agi_', '1agjA', '1agrE', '1ahsA', '1aihA', '1aisA', '1aisB', '1ak4C', '1ak7_', '1akjD', '1akp_', '1al01', '1al3_', '1allA', '1alu_', '1alvA', '1aly_', '1am2_', '1amf_', '1amm_', '1amx_', '1aoa_', '1aoeA', '1aohA', '1aol_', '1ap8_', '1aqb_', '1aqcA', '1aqe_', '1aqzA', '1ash_', '1ast_', '1at0_', '1at3A', '1atzA', '1au1A', '1auiB', '1auoA', '1auyA', '1auz_', '1avaC', '1avgI', '1avpA', '1avqA', '1awcB', '1ax8_', '1axiB']

#Teste
#Deu erro em 1a79A, 1a9nA, 1abv_, 1aisA, 1ak7_, 1akp_, 1am2_
folders = ['153l_', '1a17_', '1a1x_', '1a2pA', '1a2zA', '1a3aA', '1a3c_', '1a3k_', '1a44_', '1a53_', '1a6f_', '1a6jA', '1a6m_', '1a73A', '1a7vA', '1a8l_', '1aa7A', '1aac_', '1ac6A', '1acf_', '1acz_', '1ad2_', '1ae9A', '1agi_', '1agjA', '1agrE', '1ahsA', '1aihA', '1aisB', '1ak4C', '1akjD', '1al01', '1al3_', '1allA', '1alu_', '1alvA', '1aly_', '1amf_', '1amm_', '1amx_', '1aoa_', '1aoeA', '1aohA', '1aol_', '1ap8_', '1aqb_', '1aqcA', '1aqe_', '1aqzA', '1ash_', '1ast_', '1at0_', '1at3A', '1atzA', '1au1A', '1auiB', '1auoA', '1auyA', '1auz_', '1avaC', '1avgI', '1avpA', '1avqA', '1awcB', '1ax8_', '1axiB']

list_files = ["AffinityPropagation_Data.csv", "DBSCAN_Data.csv", "KMeans_Data.csv", "MeanShift_Data.csv", "SpectralClustering_Data.csv", "Ward_Data.csv"]

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df_ident_cluster = pd.DataFrame(columns=['prot', 'identidade', 'AffinityPropagation_Data', 'DBSCAN_Data', 'KMeans_Data', 'MeanShift_Data', 'SpectralClustering_Data', 'Ward_Data'])

file_bas = './artigo.csv'
d2 = pd.read_csv(file_bas)

NORM_METHOD = 'MinMax'
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
