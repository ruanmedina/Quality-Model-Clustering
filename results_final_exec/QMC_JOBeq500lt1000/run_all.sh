#!/bin/bash

folders="1a6cA
BAIXO_1a4sA
BAIXO_1a8i_
BAIXO_1aco_
BAIXO_1aozA
CIMA_1a4sA
CIMA_1a8i_
CIMA_1aco_
CIMA_1aozA"

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
