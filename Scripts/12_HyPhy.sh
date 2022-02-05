#!/bin/bash

##########################################################################################
### Virginie Ricci, 2021
##########################################################################################


DATE=`date '+%d-%m-%Y %H:%M:%S'`
echo 'Start of the job' $DATE
##########################################################################################

module purge
module load HyPhy/2.3.13-foss-2016b

##########################################################################################

# Don't forget to edit the files directory

method='FEL'
method='FUBAR'
method='aBSREL_shallow'
method='aBSREL_deep'

# aBSREL needs a tree with only shallow and deep species

HYPHYMP < HyPhy_${method}.txt > HyPhy_${method}.out

##########################################################################################
DATE=`date '+%d-%m-%Y %H:%M:%S'`
echo 'End of the job' $DATE
