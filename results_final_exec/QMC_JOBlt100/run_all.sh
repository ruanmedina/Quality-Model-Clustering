#!/bin/bash

folders="1a0aA 1a1w_ 1a4pA 1a8o_ 1af8_ 1afp_ 1agqA 1ail_ 1am9A 1apj_ 1atx_ 1awj_ 1a1iA 1a32_ 1a6s_ 1adr_ 1afj_ 1agg_ 1aho_ 1akhA 1aoy_ 1atb_ 1aw0_" 

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
