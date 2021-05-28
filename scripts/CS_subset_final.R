library(scales)
library(ggplot2)
library(tidyr)
library(dplyr)
library(seqinr)
library(stringr)
library(Biostrings)



#The function to get cumulative scores and other parametres for each repetition:
##amino acid composition including gaps' composition;
##amyloidogenicy;
##amino acid conservatism
CS_subset<-function(alignment =dir(pattern  ='ali_noname_R.fa'), cumscore = dir(pattern = 'CumScores.tsv'), taxonomy =dir(pattern = 'result_table_all.tsv'), id = dir(pattern='*_nog_orthologs.txt'), cs_output ='CS_subset.tsv', aa_stat ='AA_plot.pdf',aa_stat_wNA='AA_plot_wNA.pdf', summary = 'Summary_stat.tsv', cs_stat_comp = 'CS_comp_plot.pdf', summary_wNA = 'Summary_wNA.tsv'){
  #This function count noticed parametres
  ## alignment - all previously gotten alignment for each repetition are taken by this function;
  ## cumsore - gotten by ArchCandy umulative scores for each sequence in the default orthologs' dataset (one for all repetitions)
  ## taxonomy - right taxonomy gotten from NCBI db for the orthologs' subset (one for all repetitions)
  ## id - the list of ID's for the analyzing protein (one for all repetitions)
  ## cs_output - cumulative scores for each repetitions;
  ## aa_stat - plots of analyzed parametres for each position in each repetition without NAs
  ## aa_stat_wNA - plots of analyzed parametres for each position in each repetition with NAs
  ## summary - summary of analyzed parametres for subsets withot NAs
  ## CS_comp_plot - comparison plots of amyloidogenicy with and without gaps for each repetition
  ## summary_wNA - summary of parametres for subsets with NAs
  
  ali<-read.fasta(alignment)
  ali2<-ali[unique(names(ali))]
  gap_ind<-list()
  for (i in (1:length(names(ali2)))){
    gap_ind[[i]]<-which(ali2[[i]]=='-')
  }
  aa<-c('n','d', 'q', 'e', 'k', 'r','i', 'v', 'l', 'f', 'c', 'm', 'a', 'w', 'g', 't', 's', 'y', 'p', 'h')
  aa_ind<-list()
  for (i in (1:length(names(ali2)))){
    aa_ind[[i]]<-which(ali2[[i]] %in% aa)
  }
  
  CumScore <- read.table(cumscore, sep = '\t')
  tax_table <- read.table(taxonomy, header = T, sep = '\t')
  id_table <- read.table(id, header = F, sep = '\t')
  names(id_table) <- c('ID', 'Gene', 'species', 'spID', 'aliace')
  
  tax_table2 <- merge(tax_table, id_table[,c('ID', 'species', 'spID')], by.x = 'species', by.y = 'species', all.y = T)
  tax_table2 <-tax_table2[order(tax_table2$kingdom, tax_table2$phylum, tax_table2$class, tax_table2$order),]
  
  
  CumScore_w_tax <- merge(y = CumScore, x = tax_table2, by.y = 'V1', by.x = 'ID',all.y = T)
  cumscor_test<-CumScore_w_tax%>% separate(ID, into=c('ID_sp', 'ID_seq'), sep = '\\.', extra = 'merge')
  cumscor<-cumscor_test[cumscor_test$ID_seq %in% names(ali2),]
  cumscor<-cumscor[match(names(ali2), cumscor$ID_seq),]

   
  rearrangeCumScore <- function(CumScore_w_NA, alignment){
    CumScore1 = CumScore_w_NA[!(is.na(CumScore_w_NA))]
    output = numeric(length(alignment))
    n = 1
    for(i in 1:length(alignment)){
      if(alignment[i] != '-'){
        output[i] = CumScore1[n]
        n = n + 1} 
      else{output[i] = NA}}
    return(output)
  }
  cumscore_subset<-list()
  for (i in 1:length(ali2)){
    cumscore_subset[[i]]<-rearrangeCumScore(cumscor[i,9:ncol(cumscor)], ali2[[i]])
  }
  #Prepare table for output table:
  test_table<-data.frame(cumscore_subset)
  if (ncol(test_table)<nrow(test_table)){
    test_table<-as.data.frame(t(test_table))
  }
  cumscor_out<-data.frame(ID = cumscor$ID_seq, species = cumscor$species, test_table)
  rownames(cumscor_out)<-rownames(cumscor)
  write.table(cumscor_out, file = cs_output, sep = '\t', row.names = T, col.names = T, dec = ',')  
  
  names(cumscor)[2]<-'ID'
  CumScore_w_tax2<-data.frame(cumscor[,1:8], cumscor_out)
  
  CumScore_w_tax2<-CumScore_w_tax2[,c(1,2,4:7, 3, 11:ncol(CumScore_w_tax2))]
  
  #AA_statistics with amyloidogenicy:
  ali_new<-readAAStringSet(alignment, format ='fasta')
  ali_new<-ali_new[unique(ali_new@ranges@NAMES)]
  ali_new<-AAMultipleAlignment(ali_new)
  
  #Count diff aa statistic:
  consensus_m<-consensusMatrix(ali_new, as.prob = T)
  gaps_c<-consensus_m['-',]
  consensus_m['-',][consensus_m['-',]==1]<-NA
  AA_cons <- apply(X = as.data.frame(consensus_m), MARGIN = 2, FUN = max)
  polar_amyl_aa <- c('N', 'Q')
  aromatic_aa <- c('F', 'W', 'Y')
  hydrophobic_aa <-c('A','L', 'M','I','V')
  
  polar_c<-apply(X=as.data.frame(consensus_m[polar_amyl_aa,]), MARGIN =2, FUN = sum)
  hydrophobic_c<-apply(X=as.data.frame(consensus_m[hydrophobic_aa,]), MARGIN =2, FUN = sum)
  aromatic_c<-apply(X=as.data.frame(consensus_m[aromatic_aa,]), MARGIN =2, FUN = sum)
  polar_c[is.na(AA_cons)]<-NA
  aromatic_c[is.na(AA_cons)]<-NA
  hydrophobic_c[is.na(AA_cons)]<-NA
  
  #Count amyloidogenic statistic:
  CS_cons <- apply(X = CumScore_w_tax2[,8:ncol(CumScore_w_tax2)],
                   MARGIN = 2,
                   FUN = function(x){sum(x>0,na.rm = T)/length(x[!is.na(x)])})
  
  #Amyloidogenicy stat for the ful:
  CS_cons_wNA <- apply(X = CumScore_w_tax2[,8:ncol(CumScore_w_tax2)],
                       MARGIN = 2,
                       FUN = function(x){sum(x>0,na.rm = T)/length(x)})
  
  #Group table for visualization:
  domain_stat<-data.frame(Gaps = gaps_c, AA_conservatism = AA_cons, Polar_NQ = polar_c, Hydrophobic = hydrophobic_c, Aromatic = aromatic_c, Amyloidogenic = CS_cons, position = 1:ncol(consensus_m))%>%pivot_longer(cols = -c('position'))
  domain_stat_wNA<-data.frame(Gaps = gaps_c, AA_conservatism = AA_cons, Polar_NQ = polar_c, Hydrophobic = hydrophobic_c, Aromatic = aromatic_c, Amyloidogenic = CS_cons_wNA, position = 1:ncol(consensus_m))%>%pivot_longer(cols = -c('position'))
  comparison_cs<-data.frame(Gaps = gaps_c, Amyloidogenicy_wo_NA = CS_cons, Amyloidogenicy_full = CS_cons_wNA, position = 1:ncol(consensus_m))%>%pivot_longer(cols = -c('position'))
  #Output table with all stat by each position:
  write.table(domain_stat, row.names=F, col.names = T, sep='\t', file = summary)
  write.table(domain_stat_wNA, col.names = T, row.names = F, sep = '\t', file = summary_wNA)
  #Visualization of parameteres counted for alignment positions excluding NAs:
  ggplot(domain_stat,mapping = aes(x = position, y = name, fill = value)) + geom_tile(height=.9) + 
    scale_fill_viridis_c(name = "",limits = c(0.0,1.0),oob=squish) +  theme_bw()
  ggsave(filename = aa_stat,device = 'pdf', width = 20, height = 4, units = 'cm')
  
  #Visualization of parameteres counted for alignment positions including NAs:
  ggplot(domain_stat_wNA,mapping = aes(x = position, y = name, fill = value)) + geom_tile(height=.9) + 
    scale_fill_viridis_c(name = "",limits = c(0.0,1.0),oob=squish) +  theme_bw()
  ggsave(filename = aa_stat_wNA,device = 'pdf', width = 20, height = 4, units = 'cm')
  
  #Visualization of amyloidogenicy comparison for alignment positions excluding and including NAs and gaps:
  ggplot(comparison_cs,mapping = aes(x = position, y = name, fill = value)) + geom_tile(height=.9) + 
    scale_fill_viridis_c(name = "",limits = c(0.0,1.0),oob=squish) +  theme_bw()
  ggsave(filename = cs_stat_comp,device = 'pdf', width = 20, height = 4, units = 'cm')
  
}
