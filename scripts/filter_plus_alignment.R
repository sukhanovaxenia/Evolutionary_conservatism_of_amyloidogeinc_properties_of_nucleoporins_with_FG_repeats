#Filtering and alignment:
library(muscle)
library(plyr)
library(seqinr)
library(tidyr)
library(dplyr)
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
