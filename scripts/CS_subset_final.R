library(scales)
library(ggplot2)
library(tidyr)
library(dplyr)
library(seqinr)
library(stringr)




#The function to get cumulative scores and other parametres for each repetition:
##amino acid composition including gaps' composition;
##amyloidogenicy;
##amino acid conservatism
CS_subset<-function(alignment =dir(pattern  ='ali_noname_R.fa'), cumscore = dir(pattern = 'CumScores.tsv'), taxonomy =dir(pattern = 'result_table_all.tsv'), id = dir(pattern='*_nog_orthologs.txt'), cs_output ='CS_subset.tsv', PA_output ='Report_table.tsv', aa_stat ='AA_plot.pdf', summary = 'Summary_stat.tsv', cs_stat_comp = 'CS_comp_plot.pdf', summary_wNA = 'Summary_wNA.tsv', comp = 'Comparison_cs.tsv'){
## alignment - all previously gotten alignment for each repetition are taken by this function;
## cumsore - gotten by ArchCandy umulative scores for each sequence in the default orthologs' dataset (one for all repetitions)
## taxonomy - right taxonomy gotten from NCBI db for the orthologs' subset (one for all repetitions)
## id - the list of ID's for the analyzing protein (one for all repetitions)
## cs_output - cumulative scores for each repetitions;
## PA_output - number of amyloidogenic sequences for big taxonomy groups of each repetition
## aa_stat - plots of analyzed parametres for each position in each repetition
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
  #tax_table2[!complete.cases(tax_table2),]
  
  
  CumScore_w_tax <- merge(y = CumScore, x = tax_table2, by.y = 'V1', by.x = 'ID',all.y = T)
  cumscor_test<-CumScore_w_tax%>% separate(ID, into=c('ID_sp', 'ID_seq'), sep = '\\.', extra = 'merge')
  cumscor<-cumscor_test[cumscor_test$ID_seq %in% names(ali2),]
  cumscor<-cumscor[match(names(ali2), cumscor$ID_seq),]
  #cumscor2<-apply(cumscor_test[cumscor_test$ID_seq %in% names(ali2),9:ncol(cumscor_test)],1, function(x) x[!is.na(x)])
  
  
  
  
  
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
  
  #Visualization:
  #Was the problem with concatenating by merge - row multiplication inspite of equilibrity of dfs' nrows and restriction of the concatinating at all
  #CumScore_w_tax2<-merge(cumscor_out, cumscor[,1:8], by.x ='ID' , by.y = 'ID_seq', all.x = T)
  names(cumscor)[2]<-'ID'
  CumScore_w_tax2<-data.frame(cumscor[,1:8], cumscor_out)
  
  CumScore_w_tax2<-CumScore_w_tax2[,c(1,2,4:7, 3, 11:ncol(CumScore_w_tax2))]
  CumScore_w_tax_long<-CumScore_w_tax2 %>% pivot_longer(cols = -c('ID_sp','ID','kingdom','phylum','class','order', 'species'), names_to = 'position', values_to = 'score') %>%
    mutate(position = as.numeric(sub('X|V','', position))) 
  
  ggplot(CumScore_w_tax_long,mapping = aes(x = position, y = paste(class, order, ID), fill = score > 0)) + geom_tile() + theme_bw() +
    scale_fill_manual(values = c('TRUE' = 'red', 'FALSE' = 'darkblue'))
  
  
  #Summary statistics across taxons:
  CumScore_summary <- CumScore_w_tax2[,1:7]
  CumScore_summary$PA <- apply(CumScore_w_tax2[,8:ncol(CumScore_w_tax2)], 1, function(x)sum(x,na.rm = T)>0)
  
  
  PA_kingdom <- CumScore_summary %>% group_by(kingdom) %>% 
    summarise(N_kingdom = length(PA), R_kingdom = round(sum(PA)/length(PA),2))
  PA_phylum <- CumScore_summary %>% group_by(phylum) %>% 
    summarise(N_phylum = length(PA), R_phylum = round(sum(PA)/length(PA),2))
  PA_class <- CumScore_summary %>% group_by(class) %>% 
    summarise(N_class = length(PA), R_class = round(sum(PA)/length(PA),2))
  PA_order <- CumScore_summary %>% group_by(order) %>% 
    summarise(N_order = length(PA), R_order = round(sum(PA)/length(PA),2))
  
  PA_all <- merge(CumScore_summary, PA_kingdom)
  PA_all <- merge(PA_all, PA_phylum)
  PA_all <- merge(PA_all, PA_class)
  PA_all <- merge(PA_all, PA_order)
  
  PA_all_output <- distinct(PA_all[,c('kingdom', 'N_kingdom', 'R_kingdom',
                                      'phylum','N_phylum', 'R_phylum',
                                      'class','N_class', 'R_class',
                                      'order','N_order', 'R_order')])
  PA_all_output <- PA_all_output[order(PA_all_output$kingdom, 
                                       PA_all_output$phylum, 
                                       PA_all_output$class, 
                                       PA_all_output$order),]
  write.table(PA_all_output, file = PA_output, quote = F, row.names = F, dec = ',', sep = '\t')
  
  
  
  #AA_statistics with amyloidogenicy:
  ali_new<-readAAStringSet(alignment, format ='fasta')
  ali_new<-ali_new[unique(ali_new@ranges@NAMES)]
  ali_new<-AAMultipleAlignment(ali_new)
  
  #Count diff aa statistic:
  consensus_m<-consensusMatrix(ali_new, as.prob = T)
  gaps_c<-consensus_m['-',]
  consensus_m['-',][consensus_m['-',]==1]<-NA
  AA_cons <- apply(X = as.data.frame(consensus_m), MARGIN = 2, FUN = max)
  #hydrophilic<-toupper(c('n','d', 'q', 'e', 'k', 'r'))
  #hydrophobic<-toupper(c('i', 'v', 'l', 'f', 'c', 'm', 'a', 'w'))
  #neutral<- toupper(c('g', 't', 's', 'y', 'p', 'h'))
  #çàìåíèòü ãðóïïû è èç ïîëÿðíûõ îñòàâèòü òîëüêî NQ
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
  write.table(comparison_cs, col.names = F, row.names = F, sep = '\t', file = comp)
  ggplot(domain_stat,mapping = aes(x = position, y = name, fill = value)) + geom_tile(height=.9) + 
    scale_fill_viridis_c(name = "",limits = c(0.0,1.0),oob=squish) +  theme_bw()
  ggsave(filename = aa_stat,device = 'pdf', width = 20, height = 4, units = 'cm')
  
  ggplot(comparison_cs,mapping = aes(x = position, y = name, fill = value)) + geom_tile(height=.9) + 
    scale_fill_viridis_c(name = "",limits = c(0.0,1.0),oob=squish) +  theme_bw()
  ggsave(filename = cs_stat_comp,device = 'pdf', width = 20, height = 4, units = 'cm')
  
  
  #Summary table with aa profile and cum. score stat:
  #ali_df<-data.frame(names = ali_new@unmasked@ranges@NAMES, str_split_fixed(paste(ali_new),'',nchar(paste(ali_new))))
  #ali_df<-ali_df[order(match(ali_df$name,CumScore_w_tax2$ID)),]
  #names(ali_df)[1]<-'ID'
  #ali_df_w_tax<-merge(ali_df, CumScore_w_tax2[,1:7], by.x = 'ID', all.x = T, sort = F)
  #names(CumScore_w_tax2)[8:ncol(CumScore_w_tax2)]<-gsub('V','X', names(CumScore_w_tax2)[8:ncol(CumScore_w_tax2)])
  #united<-rbind(ali_df_w_tax, CumScore_w_tax2)
  #united<-united[order(united$ID),]
  
  #write.table(united, file = cs_aa_out, sep = '\t', row.names = F, col.names = T, dec = ',')
}

#Testing:
#Nsp1:
#CS_subset(alignment = 'NSP1_ali_noname_R.fa',
#cumscore = 'CumScores.tsv',
#taxonomy = 'result_table_all.tsv',
#id = '38HQQ_nog_orthologs.txt',
#cs_output = 'CS_subset_nsp1.tsv',
#PA_output = 'Report_table_nsp1.tsv',
#aa_stat = 'Nsp1_aa_stat.pdf',
#cs_aa_out = 'Nsp1_cs_aa.tsv')


#Nup100:
#CS_subset(alignment = 'NUP100_ali_noname_R.fa',
#          cumscore = 'CumScores.tsv',
 #         taxonomy = 'result_table_all.tsv',
  #        id = '38I30_nog_orthologs.txt',
   #       cs_output = 'CS_subset_nup100.tsv',
    #      PA_output = 'Report_table_nup100.tsv',
     #     aa_stat = 'Nup100_aa_stat.pdf',
      #    cs_aa_out = 'Nup100_cs_aa.tsv',
       #   summary ='Nup100_summary.tsv')

#Ð£Ð±Ñ€Ð°Ñ‚ÑŒ Ð² Ð³Ñ€Ð°Ñ„Ð¸ÐºÐµ ÐºÑ€Ð°ÑÐ½Ñ‹Ðµ Ñ€Ð°Ð¼ÐºÐ¸ + ÑÐ¾Ñ…Ñ€Ð°Ð½ÑÑ‚ÑŒ Ð½Ðµ Ð² Ð°Ð²Ñ‚Ð¾Ð¼Ð°Ñ‚Ð¸Ñ‡ÐµÑÐºÐ¾Ð¼ Ñ€Ð°Ð·Ñ€ÐµÑˆÐµÐ½Ð¸Ð¸
