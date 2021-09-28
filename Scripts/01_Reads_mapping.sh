#!/bin/bash

##########################################################################################
### Virginie Ricci, 2021
##########################################################################################


DATE=`date '+%d-%m-%Y %H:%M:%S'`
echo 'Start of the job' $DATE
##########################################################################################

module purge
module load BWA/0.7.17-goolf-1.7.20
module load SAMtools/1.7-goolf-1.7.20

##########################################################################################

# Don't forget to edit the files directory

ID=$(sed -n ${SLURM_ARRAY_TASK_ID}p Data/Dataset.txt | cut -f4) # To run as an array from 3 to 518

# Input files from Ronco et al. 2021
ID_1=${ID}_R1.fastq
ID_2=${ID}_R2.fastq

# Indexing of RefSeq (GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.fna from NCBI)
bwa index $RefSeq

# Reads mapping of ID to RefSeq
bwa mem -t 8 $ref $ID_1 $ID_2 > ${BWA}_${ID}.sam

# Convert .sam to .bam
samtools view -Sb ${BWA}_${ID}.sam > ${BWA}_${ID}.bam

# Sort .bam
samtools sort -o ${BWA}_${ID}.sorted.bam ${BWA}_${ID}.bam

# Index .bam
samtools index ${BWA}_${ID}.sorted.bam

##########################################################################################
DATE=`date '+%d-%m-%Y %H:%M:%S'`
echo 'End of the job' $DATE
