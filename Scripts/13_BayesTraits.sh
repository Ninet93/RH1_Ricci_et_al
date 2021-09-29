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

# Prepare input files and perform depth-related substitutions analysis
Rscript BayesTraits.R

##########################################################################################
DATE=`date '+%d-%m-%Y %H:%M:%S'`
echo 'End of the job' $DATE
