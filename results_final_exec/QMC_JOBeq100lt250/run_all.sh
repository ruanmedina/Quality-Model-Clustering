#!/bin/bash

folders="153l_
1a17_
1a1x_
1a2pA
1a2zA
1a3aA
1a3c_
1a3k_
1a44_
1a53_
1a6f_
1a6jA
1a6m_
1a73A
1a7vA
1a8l_
1aa7A
1aac_
1ac6A
1acf_
1acz_
1ad2_
1ae9A
1agi_
1agjA
1agrE
1ahsA
1aihA
1aisB
1ak4C
1akjD
1al01
1al3_
1allA
1alu_
1alvA
1aly_
1amf_
1amm_
1amx_
1aoa_
1aoeA
1aohA
1aol_
1ap8_
1aqb_
1aqcA
1aqe_
1aqzA
1ash_
1ast_
1at0_
1at3A
1atzA
1au1A
1auiB
1auoA
1auyA
1auz_
1avaC
1avgI
1avpA
1avqA
1awcB
1ax8_
1axiB" 

for f in ${folders} ; do
    echo "\n\n ######## RUNNING ${f} ########## \n\n"
    cp ./qmc.py ./${f}/
    cp ./pdfgrep.py ./${f}/
    mv ./${f}/*.out ./${f}/analysis.out
    cd ./${f}/
    python3 qmc.py
    cd ..
done

python3 make_ident_cluster_standardscaler.py
python3 make_ident_cluster_minmax.py
#python3 PCA_analysis.py
#mv freq_PCA_dataFrame.pickle 100_freq_PCA_dataFrame.pickle 
