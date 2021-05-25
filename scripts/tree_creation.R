#Parsing ID according to the taxonomy gotten by Lanvrentyi's script and filled up
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
  tax_rep$taxid<-as.character(tax_rep$taxid)
  tax_rep<-tax_rep[order(match(tax_rep$taxid, df$seq_ID)),]
  #Добавляем в таблицу с последовательностями названия видов, соотнося id
  for (i in 1:nrow(df)){
    if (df[i,1]==tax_rep[i,1]){
      df[i,4]<-tax_rep[i,2]
    }
  }
  names(df)[4]<-'species'
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
    colnames(tax)<-c('kingdom', 'phylum', 'class', 'order', 'species')
  } else {
    colnames(tax)<-c('kingdom', 'phylum', 'class', 'order','species')
  }
  tax<-tax[-1,]
  
  #Соединяем две таблицы через dplyr:
  tax$species<-as.character(tax$species)
  prot$species<-as.character(prot$species)
  new<-merge(tax, prot, by.x='species', by.y='species', all.y=T, all.x=T)
  new<-new[,c(2,3,4,5,6,1)]
  write.table(new,file=tax_w_ID, sep='\t', quote = F, row.names = F)
}
#Building trees for phylums and histograms of PA in classes:
table_tree<-function(taxid){
library("viridis") 
# make a new variable of color IDs
library(RColorBrewer)
ColorSpace <- plasma(11 + 5)[1:11]
cumscores<-read.table(file="CumScores.tsv", header=F, sep="\t")
cumscores1<-cumscores
cumscores1$sum<-rowSums(cumscores[,2:ncol(cumscores)], na.rm = T)
cumscores_merged<-cumscores1[,c(1,ncol(cumscores1))]
cumscores_merged<-separate(cumscores_merged, col=V1, into=c("species_ID", "protein_ID"), extra="merge")

tax1 <- read.csv(file=taxid, header = T, row.names=NULL, fill=T, sep='\t')
a1<-merge(cumscores_merged, tax1, by.x="protein_ID", by.y="seq_ID", sort=F)
##Добавление переменной-кодификатора: является/нет амилоидогенным белком (sum>0):
#b<-mutate(b, amy=ifelse(b$sum>0, T, F))
#или
b<-mutate(a1,amy=ifelse(a1$sum>0, T, F))
#Меняем порядок видов согласно таблице из первого варианта
b<-with(b, b[order(species),])
b<-b[,c(2,1,3,9,8,7,6,5,4)]
#Обновляем нумерацию строк
row.names(b)<-NULL

tax1_brief<-b[,c(9,8,7,4,2)]


class_length<-aggregate(b$amy, by = list(b$class), FUN = length)
class_amyloid<-aggregate(b$amy, by = list(b$class), FUN = sum)
class_prop<-mutate(class_amyloid, proportion = x/class_length$x)
class_prop$x<-NULL

phylum_length<-aggregate(b$amy, by = list(b$phylum), FUN = length)
phylum_amyloid<-aggregate(b$amy, by = list(b$phylum), FUN = sum)
phylum_prop<-mutate(phylum_amyloid, proportion = x/phylum_length$x)
phylum_prop$x<-NULL

kingdom_length<-aggregate(b$amy, by = list(b$kingdom), FUN = length)
kingdom_amyloid<-aggregate(b$amy, by = list(b$kingdom), FUN = sum)
kingdom_prop<-mutate(kingdom_amyloid, proportion = x/kingdom_length$x)
kingdom_prop$x<-NULL

tax1_brief<-merge(tax1_brief,merge(class_length,class_prop,by="Group.1"), by.x="class", by.y="Group.1")
tax1_brief<-merge(tax1_brief,merge(phylum_length,phylum_prop,by="Group.1"), by.x="phylum", by.y="Group.1")
tax1_brief<-merge(tax1_brief,merge(kingdom_length,kingdom_prop,by="Group.1"), by.x="kingdom", by.y="Group.1")


colnames(tax1_brief)[6:11]<-c("N protein class", "Proportion of amyloids class", "N proteins phylum", "Proportion of amyloids phylum", "N proteins kingdom", "Proportion of amyloids kingdom")

write.table(tax1_brief, file="table_base.tsv",quote=F, sep='\t', row.names=F)


tax1_brief_ID<-merge(b[,1:2], tax1_brief, by.x="protein_ID", by.y="protein_ID")
tax1_brief_ID<-unite(tax1_brief_ID, col = "ID", species_ID, protein_ID, sep = ".")
write.table(tax1_brief_ID, file="table_base_ID.tsv",quote=F, sep='\t', row.names=F)

library(ape)
#Загружаем таблицу со всеми данными
tax<-read.table(file='table_base_ID.tsv', sep='\t', header = T)
#построить гистограммы для доли амилоидогенных
hist_classes<-hist(tax$Proportion.of.amyloids.class, breaks = 10, main = "Гистограмма доли потенциальных амилоидов разных классов", xlab = 'Доля потенциальных амилоидов (ПА)', ylab = 'Количество классов', las = 1)

# As the separating of the input table was unavailable (the was the wrong result and the separator was recognized between each letter - don't know why), so wrote the new table from the input file:
tax1<-read.table(file='table_base_ID.tsv', sep='\t', header = T)
# chose right variables

phylum<-tax1[,c("kingdom","phylum", "N.proteins.phylum", "Proportion.of.amyloids.phylum")]
phylum<-unique(phylum)
phylum$Proportion.of.amyloids.phylum<-round(phylum$Proportion.of.amyloids.phylum, 2)
phylum$color<-ColorSpace[as.integer(phylum$Proportion.of.amyloids.phylum * 10 + 1)]
text_tree_phylum<-''
for(k in unique(phylum$kingdom)){
  # collect phylum names for certain kingdom
  phylum1 <- unique(phylum$phylum[phylum$kingdom == k])
  
  kingdom_string <- paste(phylum1, collapse = ',')
  # remove last comma
  #kingdom_string <- substr(kingdom_string,1,nchar(kingdom_string)-1)
  # add the name of the kingodm
  kingdom_string <- paste0('(',kingdom_string,')',k)
  # append string of the tree
  text_tree_phylum <- paste0(text_tree_phylum,kingdom_string,',')
}

# remove the last comma, add brackets and  semicolon
text_tree_phylum <- substr(text_tree_phylum,1,nchar(text_tree_phylum)-1)
text_tree_phylum <- paste0('(',text_tree_phylum,')',';')

##in ape:
tree_phylum <- read.tree(text = text_tree_phylum)

ordered_phylum <- phylum[match(tree_phylum$tip.label, phylum$phylum),]

plot(tree_phylum,tip.color = ordered_phylum$color, main = 'Phylum', show.node.label = F, edge.width=0.2,node.depth=2, cex = 0.8)
legend(x = "bottomleft", legend = paste(seq(0,0.9,0.1), seq(0.1,1.0,0.1),sep = '-'), col = ColorSpace, pch = 19,
       title = 'Доля ПА', cex = 0.5)
}

#Example:
#tax_ID(input = 'Nup42_opisto_ortho.fa', tax_list = 'tax_report.txt', tax_table = 'result_table_all.tsv', tax_w_ID = 'Nup42_tax_ID.tsv')
#table_tree(taxid = 'Nup42_tax_id.tsv')

