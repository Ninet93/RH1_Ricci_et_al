#!/bin/bash

##########################################################################################
### Virginie Ricci, 2021
##########################################################################################


DATE=`date '+%d-%m-%Y %H:%M:%S'`
echo 'Start of the job' $DATE
##########################################################################################

module purge
module load jModelTest/2.1.10-Java-1.8

##########################################################################################

# Don't forget to edit the files directory

# Nucleotide substitution model test with jModelTest
java -jar $EBROOTJMODELTEST/jModelTest.jar -tr 8 -d $CDS_file_MA_filtered -g 4 -i -f -AIC -BIC -a -s 11 -AICc -DT -p -w > 

# Nucleotide substitution model test with PAUP
./paup4a168_osx
execute $CDS_file_MA_filtered
outgroup haplo0
automodel AIC=yes AICc=yes BIC=yes DT=yes

# Amino acid substitution model test with ProtTest was performed locally (prottest-3.4.2.jar)

##########################################################################################
DATE=`date '+%d-%m-%Y %H:%M:%S'`
echo 'End of the job' $DATE
