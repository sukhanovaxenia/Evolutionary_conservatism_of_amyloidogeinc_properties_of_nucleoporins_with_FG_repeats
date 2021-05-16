#library(dplyr)
#library(plyr)
#library(tidyr)
#install.packages("BiocManager")
#BiocManager::install("Biostrings")


fastaFile <- readAAStringSet(file)
seq_name = names(fastaFile)
sequence = paste(fastaFile)
df <- data.frame(seq_name, sequence)
df$seq_name<-as.character(df$seq_name)
df<-separate(df, seq_name, c('seq_ID', 'seq_name'),sep='[.]', extra='merge',fill='right')
df<-df[order(df$seq_ID),]
df$species<-NA

#загружаем полученный tax_report и оставляем только колонки taxid и taxname
##ATTENTION: предварительно перед загрузкой надо удалить \t в .txt файле, так как это помешает соединить таблицы
tax_rep<-read.csv(file, header=T, sep='|')
tax_rep$code<-NULL
tax_rep$primary.taxid<-NULL
tax_rep$taxname<-as.character(tax_rep$taxname)
#Добавляем в таблицу с последовательностями названия видов, соотнося id
for (i in 1:nrow(df)){
  if (df[i,1]==tax_rep[i,1]){
    df[i,4]<-tax_rep[i,2]
  }
}

#Создаем файл для выравнивания (опиционально):
sequences<-as.list(as.character(df$sequence))
name<-as.list(as.character(paste(df$species, df$seq_name, sep=".")))
library(seqinr)

write.fasta(sequences,name,file.out,open='w',as.string = T)
