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
module load BLAST+/2.2.28

##########################################################################################

# Don't forget to edit the files directory

ID=$(sed -n ${SLURM_ARRAY_TASK_ID}p Data/Dataset.txt | cut -f4) # To run as an array from 3 to 518

OPSIN=RH1
#OPSIN=exoRH1

# Make BLAST database (Illumina genome assemblies from Ronco et al. 2021)
NCBI-Toolkit/21.0.0-goolf-1.7.20/bin/makeblastdb -in ${ID}.scf.fasta -out ${ID}_db -dbtype 'nucl'

# Extract individual CDS
python Extract_individual_CDS.py $CDS_file_MA_filtered $ID $path_input $path_output $OPSIN

# BLASTN of individual CDS and genome assembly
blastn -db ${ID}_db -query ${OPSIN}_${ID}.fasta -out ${OPSIN}_${ID}.blastn -outfmt '7 std stitle' -num_threads 8 -evalue 0.00001 -parse_deflines


##########################################################################################
DATE=`date '+%d-%m-%Y %H:%M:%S'`
echo 'End of the job' $DATE
