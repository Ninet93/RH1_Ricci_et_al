#!/bin/bash

##########################################################################################
### Virginie Ricci, 2021
##########################################################################################


DATE=`date '+%d-%m-%Y %H:%M:%S'`
echo 'Start of the job' $DATE
##########################################################################################

module purge
module load MrBayes/3.2.2-goolf-1.4.10-mpi
module load beagle-lib/2.1.2-goolf-1.7.20-Java-1.8.0_92
module load IQ-TREE/2.0-rc1-foss-2018b

##########################################################################################

# Don't forget to edit the files directory

# MrBayes for CDS and AA multiple alignment
mpirun -np 8 mb < MrBayes_CDS.txt > MrBayes_CDS.log
mpirun -np 8 mb < MrBayes_AA.txt > MrBayes_AA.log

# IQTree for CDS and AA multpile alignment
iqtree -s $CDS_file_MA_filtered --msub nuclear -T 8 -o $RefSeq -bb 1000 --nstep 1000 --merit AIC -m GTR+I+G --allnni
iqtree -s $AA_file_MA_filtered --msub nuclear -T 8 -o $RefSeq -bb 1000 --nstep 1000 --merit AIC -m JTT+I+G+F --allnni --seqtype AA

##########################################################################################
DATE=`date '+%d-%m-%Y %H:%M:%S'`
echo 'End of the job' $DATE
