#!/bin/bash

##########################################################################################
### Virginie Ricci, 2021
##########################################################################################


DATE=`date '+%d-%m-%Y %H:%M:%S'`
echo 'Start of the job' $DATE
##########################################################################################

module purge
module load PAML/4.9e-goolf-1.7.20

##########################################################################################

# Don't forget to edit the files directory

# Don't forget to set the foreground and background branches in the tree file

Model=M1a
#Model=M2a
#Model=M7
#Model=M8
#Model=branchsite_H0
#Model=branchsite_HA

# Branch site model needs a multiple alignment and a tree with only shallow and deep species

codeml CodeML_${Model}.ctl

##########################################################################################
DATE=`date '+%d-%m-%Y %H:%M:%S'`
echo 'End of the job' $DATE
