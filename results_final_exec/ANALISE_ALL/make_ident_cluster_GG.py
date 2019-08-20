#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 17:11:16 2019

@author: medina
"""

classes = ['QMC_JOBeq500lt1000', 'QMC_JOBeqgt1000']
folders = [
            ['1a6cA', 'BAIXO_1a4sA', 'BAIXO_1a8i_', 'BAIXO_1aco_', 'BAIXO_1aozA', 'CIMA_1a4sA', 'CIMA_1a8i_', 'CIMA_1aco_', 'CIMA_1aozA'], 
            ["1yq2", "DEBAIXO_1f4h", "DEBAIXO_6abu", "DEBAIXO_6ji0", "DECIMA_1f4h", "DECIMA_6abu", "DECIMA_6ji0"]  
          ]
list_files = ["AffinityPropagation_Data.csv", "DBSCAN_Data.csv", "KMeans_Data.csv", "MeanShift_Data.csv", "SpectralClustering_Data.csv", "Ward_Data.csv"]

NORM_METHOD = 'StandardScaler'
cluster_date_dir = 'Clusters_Data_'+NORM_METHOD

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df_ident_cluster = pd.DataFrame(columns=['prot', 'classe', 'identidade', 'AffinityPropagation_Data', 'DBSCAN_Data', 'KMeans_Data', 'MeanShift_Data', 'SpectralClustering_Data', 'Ward_Data'])

for cl in classes:
    
    n_class = classes.index(cl)

    file_bas = './../'+cl+'/artigo.csv'
    d2 = pd.read_csv(file_bas)
    
    for f in folders[n_class]:
        
        print(f)
        
        try:
            print(str(list(d2[d2['ID']==f]['Identity'])[0]))
            identity = int(str(list(d2[d2['ID']==f]['Identity'])[0]).split('%')[0].split('(')[-1])
        except:
            continue
        
        d_aux = pd.DataFrame(columns=['prot', 'classe', 'identidade', 'AffinityPropagation_Data', 'DBSCAN_Data', 'KMeans_Data', 'MeanShift_Data', 'SpectralClustering_Data', 'Ward_Data'])
        d_aux['prot'] = [f]
        d_aux['classe'] = n_class
        d_aux['identidade'] = [identity]
        
        for clt in list_files:
            
            file_cls = './../'+cl+'/'+f+'/'+cluster_date_dir+'/'+clt
            print("reading file: " + str(file_cls))
            
            meth = clt.split('.')[0]
            try:
                d = pd.read_csv(file_cls)
                n_meth = list(d["Cluster"])[-1]
                d_aux[meth] = [int(n_meth) + 1]
            except:
                continue
            
            
        df_ident_cluster = df_ident_cluster.append(d_aux, ignore_index=True)
    
    
df_ident_cluster = df_ident_cluster.dropna()

plt.figure(figsize=(6 * 2 + 3, 9.5))
plt.subplots_adjust(left=.035, right=.98, bottom=.1, top=.96, wspace=.2, hspace=.3)

X = np.array(df_ident_cluster).T
max_n = max(X[3:].reshape(-1)) + 1
min_n = min(X[3:].reshape(-1)) - 1

form = '.*dsp^'

for i in range (len(classes)):
    
    Y = df_ident_cluster[df_ident_cluster['classe']==i]
    X = np.array(Y).T
    
    print(X)
    
    plot_num = 1
    
    # colors used after on the plot
    #theColors='bgrcmykbgrcmykbgrcmykbgrcmyk'
    #colors = np.array([x for x in theColors])
    #colors = np.hstack([colors] * 20)
    
    for clt in range(len(list_files)):
        
            
            
            # plot
            plt.subplot(2, 3, plot_num)
            plt.ylim(min_n, max_n)
            plt.xlabel("Percentual de Identidade (%)")
            plt.ylabel("Número de agrupamentos")
            title = list_files[clt].split('.')[0].split('_')[0]
            plt.title(title, size=18)
            plt.scatter(X[2], X[3+clt], marker=form[i], label=classes[i])
            if title == 'DBSCAN':
                plt.legend()

            for j, name in enumerate(folders[i]):
                
                txt = name.split('_')
                if txt[-1] == '':
                    txt[-1] = txt[-2]
                if txt[0] == 'BAIXO' or txt[0] == 'DEBAIXO':
                    txt[-1] = txt[-1].upper()
                    
                print(name + "    " + txt[-1])
                
                plt.annotate(txt[-1], (list(X[2])[j], list(X[3+clt])[j]+(0.5*np.random.randn())))
            
            plot_num += 1


plt.savefig(NORM_METHOD+'_dispersion_graph_identity_n-clusters_GG.png')
