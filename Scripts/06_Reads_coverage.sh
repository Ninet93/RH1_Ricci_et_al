#!/bin/bash

##########################################################################################
### Virginie Ricci, 2021
##########################################################################################


DATE=`date '+%d-%m-%Y %H:%M:%S'`
echo 'Start of the job' $DATE
##########################################################################################
##########################################################################################

# Don't forget to edit the files directory

ID=$(sed -n ${SLURM_ARRAY_TASK_ID}p Data/Dataset.txt | cut -f4) # To run as an array from 3 to 518



# Reads coverage in the entire RefSeq genome
module purge
module load BEDTools/2.27.1-foss-2018b

bedtools genomecov -ibam ${BWA}_${ID}.bam | grep "^genome" | head 100 > ${BWA}_${ID}_genomecov.txt
	


# Reads coverage along positions in the RH1/exoRH1 reference CDS
module purge
module load SAMtools/1.7-goolf-1.7.20

samtools depth -aa -r "NC_031984.2:14827848-14828912" ${BWA}_${ID}_RH1.bam > ${BWA}_${ID}_RH1_CDS.bam.depth

samtools depth -aa -r "NC_031984.2:14387849-14388209" ${BWA}_${ID}_exoRH1.bam > ${BWA}_${ID}_exoRH1_CDS.bam.depth
samtools depth -aa -r "NC_031984.2:14388439-14388607" ${BWA}_${ID}_exoRH1.bam >> ${BWA}_${ID}_exoRH1_CDS.bam.depth
samtools depth -aa -r "NC_031984.2:14388789-14388954" ${BWA}_${ID}_exoRH1.bam >> ${BWA}_${ID}_exoRH1_CDS.bam.depth
samtools depth -aa -r "NC_031984.2:14390612-14390851" ${BWA}_${ID}_exoRH1.bam >> ${BWA}_${ID}_exoRH1_CDS.bam.depth
samtools depth -aa -r "NC_031984.2:14391687-14391806" ${BWA}_${ID}_exoRH1.bam >> ${BWA}_${ID}_exoRH1_CDS.bam.depth



module purge
module load R/3.5.0-goolf-1.7.20

# Mean and median reads coverage in reference CDS
Rscript Reads_coverage.R
output: Coverage_CDS_Overall.txt


##########################################################################################
DATE=`date '+%d-%m-%Y %H:%M:%S'`
echo 'End of the job' $DATE
