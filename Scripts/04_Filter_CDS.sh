#!/bin/bash

##########################################################################################
### Virginie Ricci, 2021
##########################################################################################


DATE=`date '+%d-%m-%Y %H:%M:%S'`
echo 'Start of the job' $DATE
##########################################################################################

module purge
module load Python/3.5.2-goolf-1.7.20
module load Biopython/1.72-goolf-1.7.20-Python-3.5.2

##########################################################################################

# Don't forget to edit the files directory

python Filter_CDS.py $CDS_file_MA $path
# outputs: $CDS_file_MA_filtered and $AA_file_MA_filtered

##########################################################################################
DATE=`date '+%d-%m-%Y %H:%M:%S'`
echo 'End of the job' $DATE
