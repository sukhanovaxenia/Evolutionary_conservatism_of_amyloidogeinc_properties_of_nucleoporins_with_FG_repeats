#Для получения списка ID:
#library(dplyr)
#library(plyr)
#library(tidyr)
#install.packages("BiocManager")
#BiocManager::install("Biostrings")
get_ID<-function(input, output){
#загрузка исходного инпута с базы EggNOG
  ##input - путь к файлу с ортологами с EggNOG в формате .fa
  ##output - пусть к файлу со списком ID

fastaFile <- readAAStringSet(file=input)
seq_name = names(fastaFile)
sequence = paste(fastaFile)
df <- data.frame(seq_name, sequence)
df$seq_name<-as.character(df$seq_name)
df<-separate(df, seq_name, c('seq_ID', 'seq_name'),sep='[.]', extra='merge',fill='right')
df<-df[order(df$seq_ID),]
df$species<-NA



#Выгружаем в отдельный файл список ID
id<-data.frame(df$seq_ID, row.names=NULL)
write.table(id, file=output, quote=F, row.names=F, col.names=F)
}