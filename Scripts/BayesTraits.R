library(ape)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(phangorn)
library(phytools)
library(seqinr)
library(stringr)

# Dataset
Dataset = read.csv('Data/Dataset.txt', sep='\t', header=T)
names(Dataset) = c('Species_abbr', 'Tribe', 'Sex', 'ID', 'Depth', 'CDS_unique_ID', 'CDS_tip_label', 'AA_unique_ID', 'CDS_ReadsCov', 'OverallRefSeq_ReadsCov', 'Ratio_CDS_Overall')
Dataset$Species_abbr = as.character(Dataset$Species_abbr)
Dataset$Species_ID = substr(Dataset$Species_abbr, 0, 6)
Dataset$ID = as.character(Dataset$ID)
Dataset$Tribe = as.character(Dataset$Tribe)
Dataset$Sex = as.character(Dataset$Sex)
Dataset$Depth = as.character(Dataset$Depth)
Dataset$CDS_unique_ID = as.character(Dataset$CDS_unique_ID)
Dataset$CDS_tip_label = as.character(Dataset$CDS_tip_label)
Dataset$AA_unique_ID = as.character(Dataset$AA_unique_ID)
Dataset$CDS_ReadsCov = as.character(Dataset$CDS_ReadsCov)
Dataset$OverallRefSeq_ReadsCov = as.character(Dataset$OverallRefSeq_ReadsCov)
Dataset$Ratio_CDS_Overall = as.character(Dataset$Ratio_CDS_Overall)


# Species tree from Ronco et al. 2021
phy_init = read.nexus('Data/b1.tre')


# CDS multiple alignment
CDS_MA = read.alignment(CDS_file_MA_filtered, format='fasta')
RefSeq_CDS = CDS_MA$seq[grep('haplo0', CDS_MA$nam)][[1]]


# Nucleotide and Amino acid tables with ambiguous characters
Nu_table = read.csv('Data/Nucleotide_table.txt', sep='\t', header=T)
AA_table = read.csv('Data/AminoAcid_table.txt', sep='\t', header=T)


# Get dataframe of codons
Codon_df = data.frame(matrix(ncol=CDS_MA$nb, nrow=nchar(RefSeq_CDS)/3)); names(Codon_df) = CDS_MA$nam
for (h in CDS_MA$nam){
  Codon_df[[h]] = strsplit(gsub("(.{3})", "\\1 ",toupper(CDS_MA$seq[grep(h, CDS_MA$nam)][[1]])), ' ')[[1]]
}


# Get dataframe of codons with ambiguous characters
Codon_ambiguous_df = data.frame(matrix(ncol=5, nrow=0))
names(Codon_ambiguous_df) = c('CDS_unique_ID', 'Position_AA', 'Initial_codon', 'Possible_codon', 'Possible_AA')
for (h in names(Codon_df)){
  for (p in seq(1, nrow(Codon_df))){
    codon_tmp = Codon_df[h][p,]
    codon_translated = translate(s2c(codon_tmp), ambiguous=T)
    
    if (codon_translated == 'X'){
      First=strsplit(codon_tmp, '')[[1]][1]
      Second=strsplit(codon_tmp, '')[[1]][2]
      Third=strsplit(codon_tmp, '')[[1]][3]
      for (nuu in c(First, Second, Third)){
        if (nuu != 'A' && nuu != 'T' && nuu != 'C' && nuu != 'G'){
          possible_nuu = Nu_table[Nu_table$Nu_code == nuu,]$Nu
          
          for (p_nuu in possible_nuu){
            p_codon = str_replace(codon_tmp, nuu, p_nuu)
            p_AA = translate(s2c(p_codon), ambiguous=T)
            tmp_df = data.frame(CDS_unique_ID=h, Position_AA=p, Initial_codon=codon_tmp, Possible_codon=p_codon, Possible_AA=p_AA)
            Codon_ambiguous_df = rbind(Codon_ambiguous_df, tmp_df)
          }
        }
      }
    }
  }
}

Codon_ambiguous_df$CDS_unique_ID = as.character(Codon_ambiguous_df$CDS_unique_ID)
Codon_ambiguous_df$Initial_codon  = as.character(Codon_ambiguous_df$Initial_codon)
Codon_ambiguous_df$Possible_codon = as.character(Codon_ambiguous_df$Possible_codon)
Codon_ambiguous_df$Possible_AA = as.character(Codon_ambiguous_df$Possible_AA)

# Get codons translated as 'X' and remove them from Codon_ambiguous_df
Codon_ambiguous_df_wo_X = Codon_ambiguous_df[!Codon_ambiguous_df$Possible_AA == 'X',]

# Manual correction - visual inspection with Geneious
Corrections = read.csv('Data/RH1_codons_manualcorrection.txt', sep='\t', header=T)
Corrections$CDS_unique_ID = as.character(Corrections$CDS_unique_ID)
Corrections$Initial_codon = as.character(Corrections$Initial_codon)
Corrections$Possible_codon = as.character(Corrections$Possible_codon)

Codon_ambiguous_df_final = full_join(Codon_ambiguous_df_wo_X, Corrections)

for (nr in rownames(Codon_ambiguous_df_final[is.na(Codon_ambiguous_df_final$Possible_AA),])){
  p_codon = Codon_ambiguous_df_final[nr,]$Possible_codon
  Codon_ambiguous_df_final[nr,]$Possible_AA = translate(s2c(p_codon), ambiguous=T)
}

# Get dataframe of AA with ambiguous characters
AA_ambiguous_df = data.frame(matrix(ncol=ncol(Codon_df), nrow=nrow(Codon_df))); names(AA_ambiguous_df) = names(Codon_df)
for (posi in seq(1, nrow(Codon_df))){
  tmp = Codon_df[posi, ]
  
  for (h in names(Codon_df)){
    tmp_h = tmp[h][[1]]

    check = Codon_ambiguous_df_final[(Codon_ambiguous_df_final$CDS_unique_ID == h) & (Codon_ambiguous_df_final$Position_AA=posi) & (Codon_ambiguous_df_final$Initial_codon == tmp_h),]

    if (nrow(check) != 0){
      aa_list = paste0(sort(unique(check$Possible_AA)), collapse='')
      aa = paste0('{', aa_list, '}')

      AA_ambiguous_df[posi, ][h] = aa


    }else{
      AA_ambiguous_df[posi, ][h] = translate(s2c(tmp_h), ambiguous=T)
    }
  }
}


# Extract species living in shallow and deep waters
Dataset_shallow = subset(Dataset, Depth == 'shallow')
Dataset_deep = subset(Dataset, Depth == 'deep')

Dataset_shallow_deep = rbind(Dataset_shallow, Dataset_deep)
Dataset_shallow_deep = Dataset_shallow_deep[!is.na(Dataset_shallow_deep$CDS_unique_ID),]
Dataset_shallow_deep_species = unique(Dataset_shallow_deep$Species_ID)

# Filter the species tree - shallow/deep species
phy_shallow_deep = drop.tip(phy_init, phy_init$tip.label[! phy_init$tip.label %in% Dataset_shallow_deep_species])
write.nexus(phy_shallow_deep, file='Data/b1_shallow_deep_species.tre')

# Filter the dataframe of AA with ambiguous characters - shallow/deep species
AA_ambiguous_df_shallow_deep = AA_ambiguous_df[c(unique(Dataset_shallow_deep$CDS_unique_ID))]


# Get dataframe of AA multi-allelic states
AA_ambiguous_df_shallow_deep_allelic = data.frame(matrix(0L, ncol=length(unique(Dataset_shallow_deep$CDS_tip_label)), nrow=nrow(AA_ambiguous_df_shallow_deep)))
names(AA_ambiguous_df_shallow_deep_allelic) = unique(Dataset_shallow_deep$CDS_tip_label)
for (l in seq(1, nrow(AA_ambiguous_df_shallow_deep_allelic))){
  for (h in names(AA_ambiguous_df)){
    tip_labels = Dataset_shallow_deep[Dataset_shallow_deep$CDS_unique_ID == h,]$CDS_tip_label
    
    h_char = AA_ambiguous_df[l, ][h][[1]]
    
    for (t in tip_labels){
      if (nchar(h_char) != 1){
        char_split = strsplit(h_char, '')[[1]]
        chars = char_split[c(-1, -length(char_split))]
        
        if(length(chars) == 1){
          AA_ambiguous_df_shallow_deep_allelic[l, ][t] = paste0(c(chars, chars), collapse='')
        }else{
          AA_ambiguous_df_shallow_deep_allelic[l, ][t] = paste0(chars, collapse='')
        }
      }else{
        AA_ambiguous_df_shallow_deep_allelic[l, ][t] = paste0(c(h_char, h_char), collapse='')
      }
    }
  }
}


# Get dataframe of only variable positions for AA multi-allelic states
AA_ambiguous_df_shallow_deep_allelic_var = AA_ambiguous_df_shallow_deep_allelic
Pos_AA = data.frame(matrix(0L, ncol=2, nrow=0))
names(Pos_AA) = c('Pos', 'AA')
variable_AAsites=c()
for (l in seq(1, nrow(AA_ambiguous_df_shallow_deep_allelic))){
 if (nrow(unique(t(AA_ambiguous_df_shallow_deep_allelic[l,]))) == 1){
   AA_ambiguous_df_shallow_deep_allelic_var[l,] = NA
   
 }else{
   possible_AA = unique(strsplit(paste(as.character(unique(t(AA_ambiguous_df_shallow_deep_allelic[l,]))), collapse=''), '')[[1]])
   variable_AAsites = append(variable_AAsites, l)
   
   for (AA in possible_AA){
     tmp = data.frame(Pos=l, AA=AA)
     Pos_AA = rbind(Pos_AA, tmp)
   }
   
 }
}


# Get dataframe of only variable positions for AA multi-allelic states per species
AA_ambiguous_df_shallow_deep_allelic_var_stats = data.frame(matrix(0L, ncol=length(Dataset_shallow_deep_species)+2, nrow=nrow(Pos_AA)))
names(AA_ambiguous_df_shallow_deep_allelic_var_stats) = c('Pos', 'AA', Dataset_shallow_deep_species)
AA_ambiguous_df_shallow_deep_allelic_var_stats$AA = as.character(AA_ambiguous_df_shallow_deep_allelic_var_stats$AA)
Pos_AA$AA = as.character(Pos_AA$AA)
AA_ambiguous_df_shallow_deep_allelic_var_stats = left_join(Pos_AA, AA_ambiguous_df_shallow_deep_allelic_var_stats)
AA_sp=c()

for (l in variable_AAsites){
  AA_sp=c()
  for (sp in Dataset_shallow_deep_species){
    AA_sp=c()
    sp_col = grep(sp, names(AA_ambiguous_df_shallow_deep_allelic_var))
    nb_sp = length(sp_col)
    
    for (col in sp_col){
      cols = AA_ambiguous_df_shallow_deep_allelic_var[col]
      cols_l = cols[l, ]
      
      AA_sp = append(AA_sp, strsplit(cols_l, '')[[1]])
    }
    
    for (AA in AA_sp){
      count_AA_tmp = as.data.frame(table(AA_sp))
      count_AA = as.integer(count_AA_tmp[count_AA_tmp$AA_sp == AA,]['Freq'])/length(AA_sp)
      row_l_AA = which(AA_ambiguous_df_shallow_deep_allelic_var_stats['Pos'] == l & AA_ambiguous_df_shallow_deep_allelic_var_stats['AA'] == AA)
      
      AA_ambiguous_df_shallow_deep_allelic_var_stats[row_l_AA,][sp] = count_AA
      
    }
  }
}

ncol_df = ncol(AA_ambiguous_df_shallow_deep_allelic_var_stats)
AA_ambiguous_df_shallow_deep_allelic_var_stats[3:ncol_df] = sapply(AA_ambiguous_df_shallow_deep_allelic_var_stats[3:ncol_df], function(x) as.numeric(x))
AA_ambiguous_df_shallow_deep_allelic_var_stats[3:ncol_df] = sapply(AA_ambiguous_df_shallow_deep_allelic_var_stats[3:ncol_df], function(x) { x[is.na(x)] = 0; x })
AA_ambiguous_df_shallow_deep_allelic_var_stats[3:ncol_df] = sapply(AA_ambiguous_df_shallow_deep_allelic_var_stats[3:ncol_df], function(x) format(round(as.numeric(x), 2), nsmall = 2))


# Get list of dataframes - per variable positions and per AA
list_per_pos_AA_stats=list()
list_per_pos_AA_stats_binary=list()
count=0
for (pos in variable_AAsites){
  count=count+1
  tmp_pos = AA_ambiguous_df_shallow_deep_allelic_var_stats[AA_ambiguous_df_shallow_deep_allelic_var_stats$Pos == pos,]
  Pos_AA = as.data.frame(t(tmp_pos[,-c(1,2)]))
  
  names(Pos_AA) = tmp_pos$AA
  Pos_AA$Species_ID = rownames(Pos_AA)
  rownames(Pos_AA) = NULL
  Pos_AA = left_join(Pos_AA, unique(Dataset[c('Species_ID', 'Depth')]))
  Pos_AA$Pos = pos
  
  Pos_AA_final = Pos_AA[c(length(Pos_AA), length(Pos_AA)-2, length(Pos_AA)-1,  seq(1, length(tmp_pos$AA)))]
  
  list_per_pos_AA_stats[[count]] = Pos_AA_final
  
  Pos_AA_final$Depth_val=1
  Pos_AA_final[Pos_AA_final$Depth == 'shallow', ]$Depth_val = 0
  
  for (n_AA in seq(4, 4+length(tmp_pos$AA)-1)){
    Pos_AA_final[,n_AA] = as.numeric(as.character(Pos_AA_final[,n_AA]))
    Pos_AA_final[,n_AA][Pos_AA_final[,n_AA] != 0] = 1
    Pos_AA_final_binary = Pos_AA_final[,c(2, n_AA, length(Pos_AA_final))]
    list_per_pos_AA_stats_binary[[count]] = Pos_AA_final_binary
    
    out_f = paste0(names(Pos_AA_final_binary)[2], '_', pos)
    write.table(Pos_AA_final_binary, file=paste0('Data/', out_f, '.txt'), quote = F, sep = '\t', col.names = F, row.names = F)
  }
  
  
  
  
}

length(variable_AAsites) == length(list_per_pos_AA_stats)

###
