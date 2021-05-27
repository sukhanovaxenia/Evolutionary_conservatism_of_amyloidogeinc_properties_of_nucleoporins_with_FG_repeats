#Filtering and alignment:
library(muscle)
library(plyr)
library(seqinr)
library(tidyr)
library(dplyr)
muscle_refine<-function(fastaFile =dir(pattern = '*opisto_ortho.fa'), file=dir(pattern = 'tax_report.txt'), tree = dir(pattern = '*_tax_ID.tsv'), filtered= 'filt.fa', filtered_noname = 'filt_noname.fa', outputpath = 'ali_R.fa', outputpath_noname = 'ali_noname_R.fa'){
  #fastaFile - фаста с набором последовательностей для белка до замены ID
  #file - файл с полученным с помощью скрипта Лаврентия таксономическим списком
  #tree - файл с полной таблицей таксономии
  #input (optional) - название файла для сохранения фасты с замененными на названия видов ID
  #filtered - название файла для сохранения отфильтрованной фасты и последующей её подгрузкой (ОБАЗЕТЕЛЕН)
  #outputpath -название для сохранения файла выравнивания
  
  #Загружаем набор последовательностей для замены ID на названия видов
  fastaFile <- readAAStringSet(fastaFile)
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
  
  #Создаем файл для фильтрации (опционально):
  #sequences<-as.list(as.character(df$sequence))
  #name<-as.list(as.character(paste(df$species, df$seq_name, sep=".")))
  #library(seqinr)
  #write.fasta(sequences,name,file.out = input,open='w',as.string = T)
  #Загружаем созданный ранее файл в формате 
  
  #input<-readAAStringSet(file=input, format = 'fasta') #path to the file with sequences with changed names
  
  #Создаем второй df для дальнейшей работы:
  df2<-data.frame(name = df$species, ID = df$seq_name, sequence = df$sequence)
  df2$name<-as.character(df2$name)
  #Загружаем файл полной таксономии: царство, фила, класс, семейство, род, вид
  tree<-read.table(tree, sep='\t', header = T)
  #Считаем частоту встречаемость каждого класса - количество последовательностей в классе
  freq<-as.data.frame(table(tree$class))
  #Добавляем информацию об объеме класса в полную таблицу таксономии
  new<-merge(tree, freq, by.x = 'class', by.y = 'Var1')
  #new$species<-paste(new$species, new$seq_ID, sep='|')
  #new<-new[,-8]
  #Объединяем таксономию с последовательностями
  full<-merge(new, df2[,c(1,3)], by.x = 'species', by.y = 'name')
  #Убираем потворяющиеся значения, возникшие при объединении
  full[,-c(9,10)]<-lapply(full[,-c(9,10)], factor)
  full2<-full[unique(full$seq_ID),]
  full2[,c(1,2,3,4,5,6,7)]<-lapply(full2[,c(1,2,3,4,5,6,7)], factor)
  #Сабсетим только классы, в которых от 3 до 10 (включительно) последовательностей
  full_3_10<-full2 %>% group_by(class) %>% filter(Freq >3 & Freq <=10)
  #Сабсетим классы, в которых более 10 последовательносетй:
  full_out_10<- full2 %>% group_by(class) %>% filter(Freq >10)
  #Выбираем случайно 10 последовательностей для классов большими наборами (>10):
  full_out_10<-ddply(full_out_10, .(class), function(x) x[sample(nrow(x),10, replace=T),])
  #Получаем итоговый df
  full_ed<-rbind.fill(full_3_10, full_out_10)
  full_ed<-full_ed[duplicated(full_ed$seq_ID)==F,]
  #Готовим переменные для создания фасты с последующей ее подгрузкой:
  sequences<-as.list(as.character(full_ed$sequence))
  name<-as.list(as.character(paste(full_ed$species, full_ed$seq_ID, sep=".")))
  write.fasta(sequences,name,file.out = filtered,open='w',as.string = T)
  #Загружаем ранее созданную фасту для выравнивания
  ali<-readAAStringSet(filepath = filtered, format = 'fasta')
  #Сырое выравнивание
  draft_muscle<-muscle(ali)
  #Refine выравнивание
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


#Test for NSP1:
#muscle_refine(fastaFile = 'NSP1_opisto_ortho.fa',
#file = "C:/Users/sukha/OneDrive/Документы/Лаборатория/New_amyloids_coaggregation/Nucleoporines/Tree NSP1, NUP100/NSP1/tax_report.txt", 
#tree = "C:/Users/sukha/OneDrive/Документы/Лаборатория/New_amyloids_coaggregation/Nucleoporines/Tree NSP1, NUP100/NSP1/NSP1_tax_ID.tsv",
#filtered = 'NSP1_filt.fa',
#filtered_noname = 'NSP1_filt_noname.fa',
#outputpath = 'NSP1_ali_R.fa',
#outputpath_noname = 'NSP1_ali_noname_R.fa')
