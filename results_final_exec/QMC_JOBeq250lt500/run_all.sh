#!/bin/bash

folders="12asA 1a02N 1a12A 1a3qA 1a4mA 1a7j_ 1a8q_ 1abrB 1ad3A 1ajsA 1ako_ 1aq0A 1aru_ 16pk_ 1a0cA 1a28A 1a48_ 1a4yA 1a8d_ 1a8y_ 1ac5_ 1afrA 1ak0_ 1amp_ 1aquA 1auiA 19hcA 1a0tP 1a2oA 1a4iA 1a59_ 1a8e_ 1a99A 1ad1A 1air_ 1ak5_ 1aop_ 1arb_ 1axn_" 

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
