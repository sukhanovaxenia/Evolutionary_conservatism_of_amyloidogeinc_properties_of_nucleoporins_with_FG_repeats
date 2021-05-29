
ID_to_sp<-function(fasta,tax,out){
  #fasta-input alignment file
  #tax - input tax_report.txt table with deleted tabs
  #out-new fasta with changed ids
fastaFile <- readAAStringSet(fasta)
df <- data.frame(seq_name=names(fastaFile), sequence=paste(fastaFile))
df$seq_name<-as.character(df$seq_name)
df<-separate(df, seq_name, c('seq_ID', 'seq_name'),sep='[.]', extra='merge',fill='right')
df<-df[order(df$seq_ID),]

#Import tax_report.txt file 
##ATTENTION: all \t should be removed before import
tax_rep<-read.table(tax, header=T, sep='|')
tax_rep$code<-NULL
tax_rep$primary.taxid<-NULL
tax_rep$taxname<-as.character(tax_rep$taxname)
#Add species names to sequences merging by IDs:
df<-merge(df, tax_rep,by.x='seq_ID', by.y = 'taxid')
df<-df[!duplicated(df$seq_name),]
names(df)[4]<-'species'

#Creating multifaste with replaced by species names IDs:
library(seqinr)

write.fasta(sequences=as.list(as.character(df$sequence)),names=as.list(as.character(paste(df$seq_name, df$species, sep="."))),file.out=out,open='w',as.string = T)
}
