#!/bin/bash

##########################################################################################
### Virginie Ricci, 2021
##########################################################################################


DATE=`date '+%d-%m-%Y %H:%M:%S'`
echo 'Start of the job' $DATE
##########################################################################################

module purge
module load R/3.5.0-goolf-1.7.20
module load BayesTraits/3.0.2-Linux-Threaded

##########################################################################################

# Don't forget to edit the files directory

# Prepare input files and perform depth-related substitutions analysis
Rscript BayesTraits_input.R

# Run BayesTraits analysis
pos=$(sed -n ${SLURM_ARRAY_TASK_ID}p Data/BayesTrait_ListFiles.txt)

mod='dependent_MC'
#mod='independent_MC'

BayesTraitsV3 Data/b1_shallow_deep_species.tre Data/{pos}.txt < Scripts/parfile_${mod}.txt > Data/${pos}_${mod}_Log.txt

# Analyse BayesTraits output files
Rscript BayesTraits.R

##########################################################################################
DATE=`date '+%d-%m-%Y %H:%M:%S'`
echo 'End of the job' $DATE
