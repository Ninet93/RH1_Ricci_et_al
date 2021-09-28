#!/bin/bash

##########################################################################################
### Virginie Ricci, 2021
##########################################################################################


DATE=`date '+%d-%m-%Y %H:%M:%S'`
echo 'Start of the job' $DATE
##########################################################################################

module purge
module load SAMtools/1.7-goolf-1.7.20
module load R/3.5.0-goolf-1.7.20

##########################################################################################

# Don't forget to edit the files directory

ID=$(sed -n ${SLURM_ARRAY_TASK_ID}p Data/Dataset.txt | cut -f4) # To run as an array from 3 to 518

# Mean and median reads coverage in the entire RefSeq genome
samtools depth ${BWA}_${ID}.bam > ${BWA}_${ID}.depth

mean=$(cut -f4 ${BWA}_${ID}.depth | datamash mean 1)
median=$(cut -f4 ${BWA}_${ID}.depth | datamash median 1)

echo ${ID}$'\t'${mean} >> Mean_ReadsCoverage.txt
echo ${ID}$'\t'${median} >> Median_ReadsCoverage.txt



# Mean and median reads coverage in reference CDS
Rscript Reads_coverage.R


##########################################################################################
DATE=`date '+%d-%m-%Y %H:%M:%S'`
echo 'End of the job' $DATE
