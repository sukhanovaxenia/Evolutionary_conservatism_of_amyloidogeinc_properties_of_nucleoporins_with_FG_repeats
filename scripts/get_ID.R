#Для получения списка ID:
#library(dplyr)
#library(plyr)
#library(tidyr)
#library(Biostrings)
get_ID<-function(input, output=out){
#загрузка исходного инпута с базы EggNOG
  ##input - путь к файлу с ортологами с EggNOG в формате .fa
  ##output - пусть к файлу со списком ID
out='taxid.txt'
fastaFile <- readAAStringSet(file=input)
df <- data.frame(seq_name=names(fastaFile))
df$seq_name<-as.character(df$seq_name)
df<-separate(df, seq_name, c('seq_ID', 'seq_name'),sep='[.]', extra='merge',fill='right')
df<-df[order(df$seq_ID),]

#Выгружаем в отдельный файл список ID
write.table(df$seq_ID, file=out, quote=F, row.names=F, col.names=F)
}
