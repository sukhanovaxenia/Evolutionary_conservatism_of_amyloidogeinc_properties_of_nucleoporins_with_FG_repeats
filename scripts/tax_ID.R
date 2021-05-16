
#library(dplyr)
#library(plyr)
#library(tidyr)
#install.packages("BiocManager")
#BiocManager::install("Biostrings")

tax_ID<-function(input,tax_list,tax_table, tax_w_ID){
  
  #На следующем сайте: https://www.ncbi.nlm.nih.gov/Taxonomy/TaxIdentifier/tax_identifier.cgi
  ##загружаем полученный на прошлом этапе (в скрипте get_ID) список id и выгружаем с сайта .txt файл с id и списком видов

#Загружаем исходный файл с ортологами

fastaFile <- readAAStringSet(file=input)
seq_name = names(fastaFile)
sequence = paste(fastaFile)
df <- data.frame(seq_name, sequence)
df$seq_name<-as.character(df$seq_name)
df<-separate(df, seq_name, c('seq_ID', 'seq_name'),sep='[.]', extra='merge',fill='right')
df<-df[order(df$seq_ID),]
df$species<-NA



#загружаем полученный tax_report и оставляем только колонки taxid и taxname
##ATTENTION: предварительно перед загрузкой надо удалить \t в .txt файле, так как это помешает соединить таблицы
tax_rep<-read.csv(file=tax_list, header=T, sep='|')
tax_rep$code<-NULL
tax_rep$primary.taxid<-NULL
tax_rep$taxname<-as.character(tax_rep$taxname)
#Добавляем в таблицу с последовательностями названия видов, соотнося id
for (i in 1:nrow(df)){
  if (df[i,1]==tax_rep[i,1]){
    df[i,4]<-tax_rep[i,2]
  }
}
prot<-data.frame(seq_ID=df$seq_name, species=df$species)
#Создаем файл для выравнивания (опиционально):
#sequences<-as.list(as.character(df$sequence))
#name<-as.list(as.character(paste(df$species, df$seq_name, sep="\t")))
#library(seqinr)

#write.fasta(sequences,name,file.out='C:/Users/sukha/OneDrive/Рабочий стол/c mac/Documents/Лаборатория/New_amyloids_coaggregation/Nucleoporines/Tree NSP1, NUP100/NSP1/NSP1_opisto_name3.fasta',open='w',as.string = T)

#Загружаем файл таксономии, полученный скриптом Лаврентия и выравненный по таксономическим группам вручную в excel :

tax<-read.csv(file=tax_table, header = F, row.names = NULL, fill = T, sep='\t')
#Немного корректируем баги, если все же не стали добавлять названия столбцов вручную:
#Если в файле от скрипта Лаврентия в таблице есть 8-ой столбец (ID), то его следует указать в следующем векторе присваивания названий - в противном случае не произойдет присваивание и будет ошибка
if (ncol(tax)==8){
  tax[,8]<-NULL
  colnames(tax)<-c('kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species')
} else {
  colnames(tax)<-c('kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species')
}


#Соединяем две таблицы через dplyr:
tax$species<-as.character(tax$species)
prot$species<-as.character(prot$species)
new<-merge(tax, prot, by.x='species', by.y='species', all.y=T, all.x=T)
new<-new[,c(2,3,4,5,6,7,1,8,9)]
write.table(new,file=tax_w_ID, sep='\t', quote = F, row.names = F)
}

#tax_ID(input = 'C:/Users/sukha/OneDrive/Рабочий стол/c mac/Documents/Лаборатория/New_amyloids_coaggregation/Nucleoporines/Tree NSP1, NUP100/NUP100/NUP100_opisto_orthologs.fa', tax_list = 'C:/Users/sukha/OneDrive/Рабочий стол/c mac/Documents/Лаборатория/New_amyloids_coaggregation/Nucleoporines/Tree NSP1, NUP100/NUP100/tax_report.txt', tax_table = 'C:/Users/sukha/OneDrive/Рабочий стол/c mac/Documents/Лаборатория/New_amyloids_coaggregation/Nucleoporines/Tree NSP1, NUP100/NUP100/NUP100_result_table.tsv', tax_w_ID = 'C:/Users/sukha/OneDrive/Рабочий стол/c mac/Documents/Лаборатория/New_amyloids_coaggregation/Nucleoporines/Tree NSP1, NUP100/NUP100/NUP100_tax_ID.tsv')

#tax_ID(input = 'C:/Users/sukha/OneDrive/Рабочий стол/c mac/Documents/Лаборатория/New_amyloids_coaggregation/Nucleoporines/Tree NSP1, NUP100/NSP1/NSP1_opisto_ortho.fa', tax_list = 'C:/Users/sukha/OneDrive/Рабочий стол/c mac/Documents/Лаборатория/New_amyloids_coaggregation/Nucleoporines/Tree NSP1, NUP100/NSP1/tax_report.txt', tax_table = 'C:/Users/sukha/OneDrive/Рабочий стол/c mac/Documents/Лаборатория/New_amyloids_coaggregation/Nucleoporines/Tree NSP1, NUP100/NSP1/NSP1_result_table2.txt', tax_w_ID = 'C:/Users/sukha/OneDrive/Рабочий стол/c mac/Documents/Лаборатория/New_amyloids_coaggregation/Nucleoporines/Tree NSP1, NUP100/NSP1/NSP1_tax_ID.tsv')
