#!/bin/bash

##########################################################################################
### Virginie Ricci, 2021
##########################################################################################


DATE=`date '+%d-%m-%Y %H:%M:%S'`
echo 'Start of the job' $DATE
##########################################################################################

module purge
module load SAMtools/1.7-goolf-1.7.20

##########################################################################################

# Don't forget to edit the files directory

ID=$(sed -n ${SLURM_ARRAY_TASK_ID}p Data/Dataset.txt | cut -f4) # To run as an array from 3 to 518

# Extract mRNA BAM region
samtools view -Sb ${BWA}_${ID}.sorted.bam "NC_031984.2:14827754-14829360" > ${BWA}_${ID}_RH1.sorted.bam
samtools index ${BWA}_${ID}_RH1.sorted.bam

samtools view -Sb ${BWA}_${ID}.sorted.bam "NC_031984.2:14387282-14393197" > ${BWA}_${ID}_exoRH1.sorted.bam
samtools index ${BWA}_${ID}_exoRH1.sorted.bam

# Next steps on Geneious (see Materials and Methods)

##########################################################################################
DATE=`date '+%d-%m-%Y %H:%M:%S'`
echo 'End of the job' $DATE
