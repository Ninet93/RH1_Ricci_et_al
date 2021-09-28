#!/bin/bash

##########################################################################################
### Virginie Ricci, 2021
##########################################################################################


DATE=`date '+%d-%m-%Y %H:%M:%S'`
echo 'Start of the job' $DATE
##########################################################################################

module purge

##########################################################################################

# Don't forget to edit the files directory

# Topology test with PAUP
./paup4a168_osx
execute $CDS_file_MA_filtered
gettrees file=MrBayes_IQTree.nexus # concatenation of all generated phylogenetic trees
outgroup RefSeq
lset nst=6 rmatrix=estimate basefreq=estimate rates=gamma shape=estimate pinvar=estimate nthreads=auto
lscores all / khtest shtest autest
describe all / diagnose

##########################################################################################
DATE=`date '+%d-%m-%Y %H:%M:%S'`
echo 'End of the job' $DATE
