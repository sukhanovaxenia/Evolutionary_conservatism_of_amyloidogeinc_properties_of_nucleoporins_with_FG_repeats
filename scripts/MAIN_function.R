#This sript contains three subfunctions and one SUPER-function which allows to get 10 different subsets of alignments and count for each of them following parametres:
##amino acid composition including gaps' composition;
##amyloidogenicy;
##amino acid conservatism

#Function to automatically get 10 different subsets of multiple alignments with MUSCLE algorithm:
muscle_refine<-function(fastaFile =dir(pattern = '*opisto_ortho.fa'), file=dir(pattern = 'tax_report.txt'), tree = dir(pattern = '*_tax_ID.tsv'), filtered= 'filt.fa', filtered_noname = 'filt_noname.fa', outputpath = 'ali_R.fa', outputpath_noname = 'ali_noname_R.fa'){
  #fastaFile - default multifasta with orthologs sequences (*_opisto_ortho.fa)
  #file - result_table_all.tsv file gained by Danilov Lavrentyii's python script
  #tree - full taxonomy table
  #input (optional) - optional multifile file where sequence IDs are replaced with species' names
  #filtered - filename for filtered multiple fasta (OBLIGATORY)
  #outputpath - file name for filtered multiple alignment with sequence labeles 'species's name.sequence ID'
  #outputpath_noname - file name for filtered multiple alignment with sequence labeles 'sequence ID'
  
  #Load multifasta file for replacing IDs on 'species's name.sequence ID':
  fastaFile <- readAAStringSet(fastaFile)
  seq_name = names(fastaFile)
  sequence = paste(fastaFile)
  df <- data.frame(seq_name, sequence)
  df$seq_name<-as.character(df$seq_name)
  df<-separate(df, seq_name, c('seq_ID', 'seq_name'),sep='[.]', extra='merge',fill='right')
  df<-df[order(df$seq_ID),]
  #Import tax_report.txt file
  ##ATTENTION: all \t should be removed from tax_report.txt file before import, otherwise it would cancel merging
  tax_rep<-read.csv(file, header=T, sep='|')
  tax_rep$code<-NULL
  tax_rep$primary.taxid<-NULL
  tax_rep$taxname<-as.character(tax_rep$taxname)
  tax_rep<-tax_rep[order(match(tax_rep$taxid,df$seq_ID)),]
  #Add species' names to the df with sequences merging by sequences' IDs:
  df<-merge(df, tax_rep,by.x='seq_ID', by.y = 'taxid')
  df<-df[!duplicated(df$seq_name),]
  names(df)[4]<-'species'
  print(str(df))
  #СCreating the second the same df for the following work:
  df2<-data.frame(name = df$species, ID = df$seq_name, sequence = df$sequence)
  df2$name<-as.character(df2$name)
  #Load the file with full taxonomy: kingdom, phylum, class, order, species
  tree<-read.table(tree, sep='\t', header = T)
  #Count frequencies of each class - the number of sequences in each of classes:
  freq<-as.data.frame(table(tree$class))
  #Add tha classes' frequencies into the full taxonomy table:
  new<-merge(tree, freq, by.x = 'class', by.y = 'Var1')
  #Merge taxonomy with sequences:
  full<-merge(new, df2[,c(1,3)], by.x = 'species', by.y = 'name')
  #Remove duplications appeared after merging:
  full[,-c(7,8)]<-lapply(full[,-c(7,8)], factor)
  full2<-full[unique(full$seq_ID),]
  full2[,c(1,2,3,4,5,6)]<-lapply(full2[,c(1,2,3,4,5,6)], factor)
  #Subsets classes with more than 3 and less than 10 sequences:
  full_3_10<-full2 %>% group_by(class) %>% filter(Freq >3 & Freq <=10)
  #Subset only classes with more than 10 sequences:
  full_out_10<- full2 %>% group_by(class) %>% filter(Freq >10)
  #Choose randomly 10 sequences for classes's subsets with more than 10 sequences:
  full_out_10<-ddply(full_out_10, .(class), function(x) x[sample(nrow(x),10, replace=T),])
  #Get the final df of filtered fasta
  full_ed<-rbind.fill(full_3_10, full_out_10)
  full_ed<-full_ed[duplicated(full_ed$seq_ID)==F,]
  #Prepare values for filtered multifasta export:
  sequences<-as.list(as.character(full_ed$sequence))
  name<-as.list(as.character(paste(full_ed$species, full_ed$seq_ID, sep=".")))
  write.fasta(sequences,name,file.out = filtered,open='w',as.string = T)
  #Loading previously filtered multifasta
  ali<-readAAStringSet(filepath = filtered, format = 'fasta')
  #Draft muscle alignment
  draft_muscle<-muscle(ali)
  #Refine alignment
  draft_muscle<-AAStringSet(draft_muscle)
  muscle_refine<-muscle(draft_muscle, refine = T)
  #Экспот результата
  muscle_refine<-AAStringSet(muscle_refine)
  writeXStringSet(muscle_refine, filepath = outputpath, format = 'fasta')
  
  #If you want to save alignment with only seq_ID:
  name2<-as.list(as.character(full_ed$seq_ID))
  write.fasta(sequences, name2, file.out=filtered_noname, open = 'w', as.string = T)
  ali_noname<-readAAStringSet(filepath = filtered_noname, format = 'fasta')
  draft_muscle_noname<-AAStringSet(muscle(ali_noname))
  muscle_refine_noname<-AAStringSet(muscle(draft_muscle_noname, refine = T))
  writeXStringSet(muscle_refine_noname, filepath = outputpath_noname, format = 'fasta')
}

#The function to gain amino acid composition:
aa_stat<-function(fasta = dir(pattern = 'ali_R.fa'), stat = 'ali_stat.tsv', plot = 'ali_aa_stat.pdf'){
  #fasta - alignment file, gained with R function with filtering (filter plus alignment)
  #stat - filename of the table with aa types proportions through alignment (for each position)
  #plot - filename of the proportions' stat graphic
  ali<-read.fasta(file = fasta)
  #Prepare a splitted by aa position alignment in a df format:
  fasta<-data.frame(names(ali), sequence = paste(ali))
  fasta2<-data.frame(names = fasta[,1], str_split_fixed(fasta$sequence, ',', nchar(as.character(fasta[,2]))))
  fasta2[,-1]<-fasta2[,-1]%>% mutate_all(as.character)
  fasta2<-data.frame(fasta2[,1], lapply(fasta2[,-1], function(x) gsub('"','',x)))
  fasta2[,2]<-data.frame(lapply(fasta2[,2], function(x) gsub('c','',x)))
  fasta2[,2]<-data.frame(lapply(fasta2[,2], function(x) gsub('\\(','',x)))
  fasta2[,2:ncol(fasta2)]<-data.frame(lapply(fasta2[,2:ncol(fasta2)], function(x) gsub(' ','',x)))
  fasta2[,2:ncol(fasta2)]<-data.frame(lapply(fasta2[,2:ncol(fasta2)], function(x) gsub('\n','',x)))
  fasta2[,2:ncol(fasta2)]<-data.frame(lapply(fasta2[,2:ncol(fasta2)], function(x) gsub('\\)','',x)))
  fasta2[,-1]<-fasta2[,-1]%>% mutate_all(as.character)
  fasta2<-fasta2[,sapply(fasta2, function(x) all(x!=''))]
  
  #Get the summary on aa contain and amount of each aa:  
  aa_contain<-unique(do.call('c',fasta2[,-1]))
  aa_table<-sapply(aa_contain, function(x)rowSums(t(fasta2[,-1])==x))
  aa_table<-t(aa_table)
  rownames(aa_table)<-aa_contain
  aa_table<-as.data.frame(aa_table)
  aa_table<-aa_table[1:nrow(aa_table)-1,]
  
  #Summary of aa types amount based on alignment content:
  aa_sum<-as.data.frame(matrix(nrow = 4, ncol = ncol(fasta2[,-1])))
  rownames(aa_sum)<-c('hydrophilic', 'hydrophobic', 'neutral', 'gap')
  aa_sum[1,1:ncol(aa_sum)]<-colSums(subset(aa_table, row.names(aa_table) %in% c('n','d', 'q', 'e', 'k', 'r')))
  aa_sum[2,1:ncol(aa_sum)]<-colSums(subset(aa_table, row.names(aa_table) %in% c('i', 'v', 'l', 'f', 'c', 'm', 'a', 'w')))
  aa_sum[3,1:ncol(aa_sum)]<-colSums(subset(aa_table, row.names(aa_table) %in% c('g', 't', 's', 'y', 'p', 'h')))
  aa_sum[4,1:ncol(aa_sum)]<-colSums(subset(aa_table, row.names(aa_table) %in% '-'))
  
  #Modify data to export only statistics (proportions of aa types) of the true length alignment:
  aa_sum<-as.data.frame(prop.table(as.matrix(aa_sum), 2))
  aa_sum<-data.frame(aa = rownames(aa_sum), aa_sum)
  aa_sum$aa<-as.character(aa_sum$aa)
  aa_sum3<-as.data.frame(apply(aa_sum, 2, function(x) x[!is.infinite(x)]))
  aa_sum3<-aa_sum3[,apply(aa_sum3, 2, function(x) all(!is.na(x)))]
  write.table(aa_sum3, file = stat, row.names = F, col.names = T, sep = '\t')
  
  #Modify data for following visualisation - tripping positions and proportions to a variable:
  aa_sum2<-pivot_longer(aa_sum, cols = -c(1), names_to = 'position', values_to = 'fraction', values_drop_na = T)
  aa_pos<-c(aa_sum2$position)
  aa_pos<-gsub('V', '',aa_pos)
  aa_sum2$position<-as.numeric(aa_pos)
  ggplot(aa_sum2)+geom_tile(aes(x = position, y = aa, fill = fraction))+theme_bw()+scale_fill_viridis_c(name = "",limits = c(0.0,1.0),oob=squish)
}

#The function to get cumulative scores for each repetition:
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

#The MAIN-function, including all three previous. No arguments should be told. It is launched in the directory for each protein separetly. 
##Previous functions should not be launched as all of them are in this function.
set.seed(10)
MAIN_fun<-function(x = dir(x, all.files = T), ..deparseOpts){
  i=1
  while (i<=10) { #заменить на while
    print(i)
    ali = paste('ali_noname_R',i, '.fa', sep ='')
    filt = paste('filt_noname_R',i,'.fa', sep='')
    stat_aa = paste('ali_stat',i,'.tsv', sep = '')
    plot_aa = paste('ali_aa_stat', i, '.pdf', sep = '')
    cs_out = paste('CS_subset', i, '.tsv', sep = '')
    aa_stplot = paste('AA_plot', i, '.pdf', sep= '')
    aa_stplot_na = paste('AA_plot_wNA',i,'.pdf',sep='')
    cs_comp_plot = paste('CS_comp_plot',i,'.pdf',sep='')
    sum_st = paste('Summary_stat', i, '.tsv', sep = '')
    sumNA = paste('Summary_wNA', i, '.tsv', sep = '')
    muscle_refine(outputpath_noname = ali, filtered_noname = filt)
    aa_stat(fasta =ali, stat = stat_aa, plot = plot_aa )
    CS_subset(alignment = ali, cs_output = cs_out, aa_stat = aa_stplot, aa_stat_wNA=aa_stplot_na, cs_stat_comp = cs_comp_plot, summary = sum_st, summary_wNA = sumNA)
    i = i+1  
  }
}

