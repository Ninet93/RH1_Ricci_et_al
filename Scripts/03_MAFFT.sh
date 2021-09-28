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
module load MAFFT/7.310-goolf-1.7.20-with-extensions

##########################################################################################

# Don't forget to edit the files directory

# Multiple alignment of all consensus mRNA plus the RefSeq mRNA
mafft --auto --thread 8 $mRNA_file > $mRNA_file_MA

# Extraction of all consensus CDS plus the RefSeq CDS (as no gap was found in the mRNA multiple alignment, no need to perform MAFFT again)
python Extract_CDS.py $mRNA_file_MA Data/Coo_mRNA_Exon_CDS_RefSeq.txt $path
# output: $CDS_file_MA


##########################################################################################
DATE=`date '+%d-%m-%Y %H:%M:%S'`
echo 'End of the job' $DATE
