#!/bin/bash

folders="1yq2
DEBAIXO_1f4h
DEBAIXO_6abu
DEBAIXO_6ji0
DECIMA_1f4h
DECIMA_6abu
DECIMA_6ji0"

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
