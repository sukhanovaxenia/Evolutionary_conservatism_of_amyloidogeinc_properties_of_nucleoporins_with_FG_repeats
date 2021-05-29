#For extracting IDs of taxonomic groups:

get_ID<-function(input, output=out){
  ##input - path to the input file *_opisto_ortho.fa from EggNOG
  ##output - path to the file with taxonomy IDs
out='taxid.txt'
fastaFile <- readAAStringSet(file=input)
df <- data.frame(seq_name=names(fastaFile))
df$seq_name<-as.character(df$seq_name)
df<-separate(df, seq_name, c('seq_ID', 'seq_name'),sep='[.]', extra='merge',fill='right')
df<-df[order(df$seq_ID),]
#Export file with taxonomy IDs:
write.table(df$seq_ID, file=out, quote=F, row.names=F, col.names=F)
}
