

tax_ID<-function(input,tax_list,tax_table, tax_w_ID){
  
  #On the following site: https://www.ncbi.nlm.nih.gov/Taxonomy/TaxIdentifier/tax_identifier.cgi
  ##load created taxid.txt  (в скрипте get_ID) list and export from the source tax_report.txt file with id and species' list
  
  #Import the default multifasta of orthologs *_opisto_ortho.fa:
  
  fastaFile <- readAAStringSet(file=input)
  seq_name = names(fastaFile)
  sequence = paste(fastaFile)
  df <- data.frame(seq_name, sequence)
  df$seq_name<-as.character(df$seq_name)
  df<-separate(df, seq_name, c('seq_ID', 'seq_name'),sep='[.]', extra='merge',fill='right')
  df<-df[order(df$seq_ID),]  
  
  #Import taxid.txt and tax_report.txt
  ##ATTENTION: all \t should be removed from tax_report.txt
  tax_rep<-read.csv(file=tax_list, header=T, sep='|')
  tax_rep$code<-NULL
  tax_rep$primary.taxid<-NULL
  tax_rep$taxname<-as.character(tax_rep$taxname)
  tax_rep$taxid<-as.character(tax_rep$taxid)
  tax_rep<-tax_rep[order(match(tax_rep$taxid, df$seq_ID)),]
  #Merge sequences and species names by IDs:
  df<-merge(df, tax_rep,by.x='seq_ID', by.y = 'taxid')
  df<-df[!duplicated(df$seq_name),]
  names(df)[4]<-'species'
  prot<-data.frame(seq_ID=df$seq_name, species=df$species)
 
  #Import taxonomy summary achieved by Danilov Lavrentyii's python scropt - result_table_all.tsv:
  
  tax<-read.csv(file=tax_table, header = F, row.names = NULL, fill = T, sep='\t')
  if (ncol(tax)==8){
    tax[,8]<-NULL
    colnames(tax)<-c('kingdom', 'phylum', 'class', 'order', 'species')
  } else {
    colnames(tax)<-c('kingdom', 'phylum', 'class', 'order','species')
  }
  tax<-tax[-1,]
  
  #Merge taxonomy summary table with sequences table by species' names:
  tax$species<-as.character(tax$species)
  prot$species<-as.character(prot$species)
  new<-merge(tax, prot, by.x='species', by.y='species', all.y=T, all.x=T)
  new<-new[,c(2,3,4,5,6,1)]
  write.table(new,file=tax_w_ID, sep='\t', quote = F, row.names = F)
}
