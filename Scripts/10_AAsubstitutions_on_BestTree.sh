#!/bin/bash

##########################################################################################
### Virginie Ricci, 2021
##########################################################################################


DATE=`date '+%d-%m-%Y %H:%M:%S'`
echo 'Start of the job' $DATE
##########################################################################################

module purge
module load R/3.5.0-goolf-1.7.20

##########################################################################################

# Don't forget to edit the files directory

# Mapping of AA on gene trees with PAUP
./paup4a168_osx
execute $CDS_file_MA_filtered
gettrees file=BestTree.nexus
outgroup RefSeq
describe / apo=y chg=yes diagnose=yes brlens=yes

Rscript PAUP_R_trees.R

##########################################################################################
DATE=`date '+%d-%m-%Y %H:%M:%S'`
echo 'End of the job' $DATE
