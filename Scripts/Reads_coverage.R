library(dplyr)
library(BoSSA)
library(phytools)
library(seqinr)
library(stringr)
library(castor)
library(ggplot2)
library(zoo)
library(gridExtra)

# Dataset
Dataset = read.table('Data/Dataset.txt', sep='\t', header=T)

IDs = Dataset$ID

# Coordinates of RefSeq genes
COO_Onil = read.table('Data/Coo_mRNA_Exon_CDS_RefSeq.csv', sep=',', header = T)

OPSIN='RH1'
#OPSIN='exoRH1'

CDS_start = COO_Onil[(COO_Onil$Opsin == OPSIN) & (COO_Onil$Type == 'CDS'), ]$From_mRNA_start; CDS_start
CDS_stop = COO_Onil[(COO_Onil$Opsin == OPSIN) & (COO_Onil$Type == 'CDS'), ]$From_mRNA_stop; CDS_stop

# Reads coverage in the entire RefSeq genome and in the RefSeq CDS
Mean_median_CDS_GW = data.frame(matrix(ncol=5, nrow=0))
names(Mean_median_CDS_GW) = c('ID', 'CDS_mean', 'CDS_median', 'GW_mean', 'GW_median')

for (ID in IDs){
  RC_ID_all = read.table('${BWA}_', ID '_genomecov.txt', sep='\t', header=T)
  names(RC_ID_all) = c('Scf', 'CovDepth', 'Bp', 'Scf_bp', 'Ratio_Bp_Scf_bp')
  
  RC_ID = RC_ID_all[(RC_ID_all$CovDepth >= 1) & (RC_ID_all$CovDepth <= 50),]
  
  mean_GW = round(sum((RC_ID$CovDepth * RC_ID$Ratio_Bp_Scf_bp)),3)
  
  half_len_gen = sum(RC_ID$Bp)/2; half_len_gen
  half_len_gen_plus_one = half_len_gen + 1
  
  sum_bp=0
  for (i in seq(1, nrow(RC_ID))){
  
    if ((unique(RC_ID$Scf_bp) %% 2) == 0){
    
      sum_bp = sum_bp + RC_ID[i,]$Bp
    
      if (sum_bp == half_len_gen){
        median_GW=(RC_ID[i,]$CovDepth+RC_ID[i+1,]$CovDepth)/2
        break
      
      } else if (sum_bp > half_len_gen){
        median_GW=RC_ID[i,]$CovDepth
        break
      }
    }
  }
  
  RC_ID_CDS = read.table('${BWA}_', ID, '_', OPSIN, '_CDS.bam.depth', sep='\t', header=F)
  names(RC_ID_CDS) = c('Scf', 'Position', 'Cov')
  
  #CDS_start == min(RC_ID_CDS$Position)
  #CDS_stop == max(RC_ID_CDS$Position)
  
  mean_CDS = mean(RC_ID_CDS$Cov)
  median_CDS = median(RC_ID_CDS$Cov)
  
  tmp = data.frame(ID=ID, CDS_mean=mean_CDS, CDS_median=median_CDS, GW_mean=mean_GW, GW_median=median_GW)
  Mean_median_CDS_GW = rbind(Mean_median_CDS_GW, tmp)

}

Mean_median_CDS_GW$CDS_mean = round(Mean_median_CDS_GW$CDS_mean, 2)
Mean_median_CDS_GW$CDS_median = round(Mean_median_CDS_GW$CDS_median, 2)
Mean_median_CDS_GW$GW_mean = round(Mean_median_CDS_GW$GW_mean, 2)
Mean_median_CDS_GW$GW_mean = round(Mean_median_CDS_GW$GW_mean, 2)
Mean_median_CDS_GW$CDS_GW_mean = round(Mean_median_CDS_GW$CDS_mean / Mean_median_CDS_GW$GW_mean, 2)
Mean_median_CDS_GW$CDS_GW_medan = round(Mean_median_CDS_GW$CDS_median / Mean_median_CDS_GW$GW_median, 2)


write.table(Mean_median_CDS_GW, 'Coverage_CDS_Overall.txt', sep='\t', row.names = F, quote = F)
