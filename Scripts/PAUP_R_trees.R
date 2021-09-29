library(ape)
library(dplyr)
library(phangorn)
library(phytools)
library(seqinr)
library(stringr)

BestTree='BestTree_ML'
BestTree='BestTree_Bayes'

tree = ladderize(read.nexus(paste0(BestTree, '.nexus')))

total_nodes_tips = length(tree$tip.label) + tree$Nnode
nb_tips = seq(1, length(tree$tip.label))
nb_nodes = seq(length(tree$tip.label)+1, total_nodes_tips)

PAUP_tree = paste0('PAUP_', BestTree, '_nodes_AA.txt')

### PAUP nodes - PAUP_*_nodes_AA.txt generated with PAUP describe command
PAUP_nodes = read.csv(PAUP_tree, sep='\t', header=T) 

names(PAUP_nodes)[which(names(PAUP_nodes) == 'Node')] = 'PAUP_node'
names(PAUP_nodes)[which(names(PAUP_nodes) == 'Connected_node')] = 'PAUP_connected_node'

PAUP_nodes$PAUP_node = as.character(PAUP_nodes$PAUP_node)

PAUP_nodes_R = data.frame(matrix(ncol=5, nrow=0))
names(PAUP_nodes_R) = c('NodeTips', 'PAUP_node', 'PAUP_connected_node', 'R_node', 'R_connected_node')

for (no in nb_tips){
  t_tips = tree$tip.label[no]

  Paup_tmp = PAUP_nodes[PAUP_nodes$NodeTips == t_tips,]

  Paup_NodeTips = Paup_tmp$NodeTips
  Paup_node = Paup_tmp$PAUP_node
  Paup_connode = Paup_tmp$PAUP_connected_node
  R_node = no
  R_connode = tree$edge[which(tree$edge[,2] == no),][1]

  to_add = data.frame(t(c(Paup_NodeTips, Paup_node, Paup_connode, R_node, R_connode)))
  names(to_add) = names(PAUP_nodes_R)
  PAUP_nodes_R = rbind(PAUP_nodes_R, to_add)
}


while (nrow(PAUP_nodes_R) != nrow(PAUP_nodes)){
  for (no in setdiff(c(nb_tips, nb_nodes), PAUP_nodes_R$R_node)){
    Paup_tmp = PAUP_nodes_R[PAUP_nodes_R$R_connected_node == no,]
      
    if (dim(Paup_tmp)[1] != 0) {
      edge_tmp = tree$edge[which(tree$edge[,2] == no),]
        
      Paup_NodeTips = 'Node'
      Paup_node = unique(Paup_tmp$PAUP_connected_node)
      Paup_connode = PAUP_nodes[PAUP_nodes$PAUP_node == Paup_node,]$PAUP_connected_node
      R_node = no
      R_connode = edge_tmp[1]
        
      if (Paup_connode == 'root'){
        R_connode = 'root'
      }
        
      to_add = data.frame(t(c(Paup_NodeTips, Paup_node, Paup_connode, R_node, R_connode)))
      names(to_add) = names(PAUP_nodes_R)
      PAUP_nodes_R = rbind(PAUP_nodes_R, to_add)
    } 
  }
  PAUP_nodes_R = unique(PAUP_nodes_R)
}

### PAUP apo - PAUP_*_apo_AA.txt generated with PAUP describe command
PAUP_apo_f = str_replace(PAUP_tree, 'nodes', 'apo')
PAUP_apo = read.csv(PAUP_apo_f, sep='\t', header=T)
names(PAUP_apo)[which(names(PAUP_apo) == 'Node')] = 'PAUP_node'
names(PAUP_apo)[which(names(PAUP_apo) == 'Connected_node')] = 'PAUP_connected_node'
  
PAUP_apo$PAUP_node = str_replace(PAUP_apo$PAUP_node, 'node_', '')
PAUP_apo$PAUP_connected_node = str_replace(PAUP_apo$PAUP_connected_node, 'node_', '')
  
PAUP_apo$Changes = paste0(PAUP_apo$Change_start, PAUP_apo$Character, PAUP_apo$Change_stop)
  
PAUP_apo$NodeTips = 'Node'
  
for (no in PAUP_apo$PAUP_node){
  if (length(grep('_', no)) != 0){
    lines_no = which(PAUP_apo$PAUP_node == no)
      
    for (l_no in lines_no){
      PAUP_apo[l_no,]$NodeTips = PAUP_apo[l_no,]$PAUP_node
      
      PAUP_apo[l_no,]$PAUP_node = PAUP_nodes[PAUP_nodes$NodeTips == no,]$PAUP_node        
    }
  }
}

PAUP_nodes_R_apo = left_join(PAUP_nodes_R, PAUP_apo)

### R tree with PAUP info
tree_edge = as.data.frame(tree$edge)
names(tree_edge) = c('R_connected_node', 'R_node')
tree_edge$R_connected_node = as.character(tree_edge$R_connected_node)
tree_edge$R_node = as.character(tree_edge$R_node)

PAUP_nodes_R_apo_edge = left_join(tree_edge, PAUP_nodes_R_apo)

PAUP_nodes_R_apo_edge_label = data.frame(matrix(ncol=6, nrow=0))
names(PAUP_nodes_R_apo_edge_label) = names(PAUP_nodes_R_apo_edge)[c(1:5, 13)]

for (pn in PAUP_nodes_R_apo_edge$R_connected_node){
  pn_df = PAUP_nodes_R_apo_edge[PAUP_nodes_R_apo_edge$R_connected_node == pn,]
    
  pn_df_cn = unique(pn_df$R_node)
    
  for (cn in pn_df_cn){
    pn_df_cn_df = pn_df[pn_df$R_node == cn,]
      
    pn_df_cn_df_changes = pn_df_cn_df$Changes
      
    pn_df_cn_df_changes_list = paste(pn_df_cn_df_changes, collapse = '_')
      
    R_connected_node = unique(pn_df_cn_df$R_connected_node)
    R_node = unique(pn_df_cn_df$R_node)
    NodeTips = unique(pn_df_cn_df$NodeTips)
    PAUP_node = unique(pn_df_cn_df$PAUP_node)
    PAUP_connected_node = unique(pn_df_cn_df$PAUP_connected_node)
    Changes = pn_df_cn_df_changes_list
      
    to_add = data.frame(t(c(R_connected_node, R_node, NodeTips, PAUP_node, PAUP_connected_node, Changes)))
    names(to_add) = names(PAUP_nodes_R_apo_edge)[c(1:5, 13)]
    PAUP_nodes_R_apo_edge_label = rbind(PAUP_nodes_R_apo_edge_label, to_add)
  }
}
PAUP_nodes_R_apo_edge_label = unique(PAUP_nodes_R_apo_edge_label)
  
PAUP_nodes_R_apo_edge_label_order = left_join(tree_edge, PAUP_nodes_R_apo_edge_label); head(PAUP_nodes_R_apo_edge_label_order); dim(PAUP_nodes_R_apo_edge_label_order)
PAUP_nodes_R_apo_edge_label_order[PAUP_nodes_R_apo_edge_label_order$Changes == 'NA',]$Changes = ''

R_node_order = data.frame(R_node=str_sort(PAUP_nodes_R$R_node, numeric=T))
PAUP_nodes_R_order = left_join(R_node_order, PAUP_nodes_R, by='R_node'); head(PAUP_nodes_R_order); dim(PAUP_nodes_R_order)
PAUP_nodes_R_order_nodes = PAUP_nodes_R_order[PAUP_nodes_R_order$NodeTips == 'Node', ]



