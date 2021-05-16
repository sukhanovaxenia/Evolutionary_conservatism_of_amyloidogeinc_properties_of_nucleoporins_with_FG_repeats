library(scales)
library(ggplot2)
library(tidyr)
library(dplyr)
library(seqinr)
library(stringr)


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
 #pdf(plot)
 ggplot(aa_sum2)+geom_tile(aes(x = position, y = aa, fill = fraction))+theme_bw()+scale_fill_viridis_c(name = "",limits = c(0.0,1.0),oob=squish)
 #dev.off()
 ##pdf(plot) and dev.off() as argument 'plot' are for the case is exporting by a fuction is preferable (honestly, I got empty pdf files after function call)
}

#Test commands:
##NSP1:
#aa_stat(fasta = 'NSP1_ali_R.fa', stat = 'nsp1_ali_stat.tsv', plot = 'nsp1_aa_stat.pdf')
##NUP100:
#aa_stat(fasta = 'NUP100_ali_R.fa', stat = 'nup100_ali_stat.tsv', plot = 'nup100_aa_stat.pdf')


