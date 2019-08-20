#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 17:11:16 2019

@author: medina
"""

classes = ['QMC_JOBlt100', 'QMC_JOBeq100lt250', 'QMC_JOBeq250lt500', 'QMC_JOBeq500lt1000', 'QMC_JOBeqgt1000']
folders = [
            ["1a0aA", "1a4pA", "1a8o_", "1af8_", "1afp_", "1agqA", "1ail_", "1am9A", "1apj_", "1atx_", "1awj_", "1a1iA", "1a32_", "1adr_", "1afj_", "1agg_", "1aho_", "1akhA", "1aoy_", "1atb_", "1aw0_"],
            ['153l_', '1a17_', '1a1x_', '1a2pA', '1a2zA', '1a3aA', '1a3c_', '1a3k_', '1a44_', '1a6f_', '1a6jA', '1a6m_', '1a73A', '1a7vA', '1a8l_', '1aa7A', '1aac_', '1ac6A', '1acf_', '1acz_', '1ad2_', '1ae9A', '1agi_', '1agjA', '1ahsA', '1aihA', '1aisB', '1ak4C', '1akjD', '1al01', '1al3_', '1allA', '1alu_', '1alvA', '1aly_', '1amf_', '1amm_', '1amx_', '1aoa_', '1aoeA', '1aohA', '1aol_', '1ap8_', '1aqb_', '1aqcA', '1aqe_', '1aqzA', '1ash_', '1ast_', '1at0_', '1at3A', '1atzA', '1au1A', '1auiB', '1auoA', '1auyA', '1auz_', '1avaC', '1avgI', '1avpA', '1avqA', '1awcB', '1ax8_', '1axiB'],
            ['12asA', '16pk_', '19hcA', '1a02N', '1a0cA', '1a0tP', '1a12A', '1a28A', '1a2oA', '1a3qA', '1a48_', '1a4iA', '1a4mA', '1a4yA', '1a59_', '1a7j_', '1a8d_', '1a8e_', '1a8q_', '1a8y_', '1a99A', '1ac5_', '1ad1A', '1ad3A', '1afrA', '1air_', '1ajsA', '1ak0_', '1ak5_', '1ako_', '1amp_', '1aop_', '1aq0A', '1aquA', '1arb_', '1aru_', '1auiA', '1axn_'],
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
            
            plot_num += 1


plt.savefig(NORM_METHOD+'_dispersion_graph_identity_n-clusters_all.png')
