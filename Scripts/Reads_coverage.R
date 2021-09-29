library(dplyr)
library(BoSSA)
library(phytools)
library(seqinr)
library(stringr)
library(castor)
library(ggplot2)
library(zoo)
library(gridExtra)

# Mean reads coverage in the entire RefSeq genome
MeanCov = read.table('Mean_ReadsCoverage.txt', sep='\t', header=F); names(MeanCov) = c('ID', 'Overall_MeanCov')

# Median reads coverage in the entire RefSeq genome
MedianCov = read.table('Median_ReadsCoverage.txt', sep='\t', header=F); names(MedianCov) = c('ID', 'Overall_MedianCov')

# Sequencing files
Sequencing = read.table('Sequencing_files.txt', sep='\t', header = F); names(Sequencing) = c('ID', 'Nb_SeqRun')

# Coordinates of RefSeq genes
COO_Onil = read.table('Coo_mRNA_Exon_CDS_RefSeq.csv', sep=',', header = T)

OPSIN='RH1'
#OPSIN='exoRH1'

CDS_start = COO_Onil[(COO_Onil$Opsin == OPSIN) & (COO_Onil$Type == 'CDS'), ]$From_mRNA_start; CDS_start
CDS_stop = COO_Onil[(COO_Onil$Opsin == OPSIN) & (COO_Onil$Type == 'CDS'), ]$From_mRNA_stop; CDS_stop

# Reads coverage in the RefSeq mRNA (concatenation of Geneious reads depth outputs for each ID)
Coverage_tmp = read.table('ReadsCoverage.txt', sep='\t', header=F); names(Coverage_tmp) = c('ID', 'Gene', 'Position', 'Coverage')

min_to_max_coverage = seq(min(Coverage_tmp$Position), max(Coverage_tmp$Position))

# Add lines of coordinates with 0 as coverage (absent Geneious reads depth outputs)
for (id in unique(Coverage_tmp$ID)){
  sub_id = Coverage_tmp[Coverage_tmp$ID == id,]
  pos_wo_cov = data.frame(Position=setdiff(min_to_max_coverage, sub_id$Position))
  if (nrow(pos_wo_cov) != 0){
    pos_wo_cov$Species_ID = unique(sub_id$Species_ID)
    pos_wo_cov$ID = unique(sub_id$ID)
    pos_wo_cov$Gene = unique(sub_id$Gene)
    pos_wo_cov$Coverage = 0
    Coverage_tmp = rbind(Coverage_tmp, pos_wo_cov)
  }
}

Coverage = Coverage_tmp[order(Coverage_tmp$Species_ID, Coverage_tmp$ID, Coverage_tmp$Position),]

# Extract CDS reads coverage
Coverage_CDS = Coverage[(Coverage$Position >= CDS_start) & (Coverage$Position <= CDS_stop),]
Coverage_CDS$Position_CDS = Coverage_CDS$Position - CDS_start + 1

# Get mean and median reads coverage in the RefSeq CDS
Coverage_CDS$CDS_MeanCov = ''
Coverage_CDS$CDS_MedianCov = ''


for (sp_ID in unique(Coverage$Species_ID)){
  IDs = unique(Coverage[Coverage$Species_ID == sp_ID,]$ID)

  for (id in IDs){
    Coverage_CDS[Coverage_CDS$ID == id,]$CDS_MeanCov = as.numeric(mean(Coverage_CDS[Coverage_CDS$ID == id,]$Coverage))
    Coverage_CDS[Coverage_CDS$ID == id,]$CDS_MedianCov = as.numeric(median(Coverage_CDS[Coverage_CDS$ID == id,]$Coverage))
  }
}


Coverage_CDS_Overall_tmp0 = left_join(Coverage_CDS, MeanCov)
Coverage_CDS_Overall_tmp = left_join(Coverage_CDS_Overall_tmp0, MedianCov)
Coverage_CDS_Overall = left_join(Coverage_CDS_Overall_tmp, Sequencing)

Coverage_CDS_Overall_final = unique(Coverage_CDS_Overall[c('ID', 'Nb_SeqRun', 'CDS_MeanCov', 'Overall_MeanCov, 'Ratio_Mean_CDS_Overall', 'CDS_MedianCov', 'Overall_MedianCov', 'Ratio_Median_CDS_Overall')])

Coverage_CDS_Overall_final$Ratio_Mean_CDS_Overall = (Coverage_CDS_Overall_final$CDS_MeanCov / Coverage_CDS_Overall_final$Overall_MeanCov) / Coverage_CDS_Overall_final$Nb_SeqRun
Coverage_CDS_Overall_final$Ratio_Median_CDS_Overall = (Coverage_CDS_Overall_final$CDS_MedianCov / Coverage_CDS_Overall_final$Overall_MedianCov) / Coverage_CDS_Overall_final$Nb_SeqRun

write.table(Coverage_CDS_Overall_final, 'Coverage_CDS_Overall.txt', sep='\t', row.names = F, quote = F)
