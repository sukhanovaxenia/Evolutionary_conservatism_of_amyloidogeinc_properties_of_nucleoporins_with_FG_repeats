setwd("C:/Users/sukha/OneDrive/Рабочий стол/c mac/Documents/Лаборатория/New_amyloids_coaggregation/Nucleoporines/Tree NSP1, NUP100/NSP1/tree_vis_test")

#Данный скрипт построен на примере NSP1 и использовании соответствующих директории, названий файлов и путей для их записи

library(ape)
library(tidyr)
library(plyr)
library(dplyr)
#Загружаем дерево, полученное на выравнивании исходной мультифасты с EggNOG (БЕЗ добавления списка вида)
tree <- read.tree(file = 'NSP1_opisto_ortho_tree.nwk')
#Загружаем исходную с EggNOG мультифасту
fastaFile <- readAAStringSet(file='C:/Users/sukha/OneDrive/Рабочий стол/c mac/Documents/Лаборатория/New_amyloids_coaggregation/Nucleoporines/Tree NSP1, NUP100/NSP1/NSP1_opisto_ortho.fa')
seq_name = names(fastaFile)
sequence = paste(fastaFile)
gene_list <- data.frame(seq_name, sequence)
gene_list$seq_name<-as.character(gene_list$seq_name)
gene_list<-separate(gene_list, seq_name, c('seq_ID', 'seq_name'),sep='[.]', extra='merge',fill='right')
gene_list<-gene_list[order(gene_list$seq_ID),]
gene_list$species<-NA

#Загружаем загруженную с ncbi таблицу ID-вид для добавления правильных названий видов в таблицу от мультифасты:
tax_rep<-read.csv(file='C:/Users/sukha/OneDrive/Рабочий стол/c mac/Documents/Лаборатория/New_amyloids_coaggregation/Nucleoporines/Tree NSP1, NUP100/NSP1/tax_report.txt', header=T, sep='|')
tax_rep$code<-NULL
tax_rep$primary.taxid<-NULL
tax_rep$taxname<-as.character(tax_rep$taxname)
#Добавляем в таблицу с последовательностями названия видов, соотнося с id
for (i in 1:nrow(gene_list)){
       if (gene_list[i,1]==tax_rep[i,1]){
           gene_list[i,4]<-tax_rep[i,2]
      }
}

#Опционально - сохранение таблицы вид - ID - последовательность для последующего использования:
write.table(nsp1, file='C:/Users/sukha/OneDrive/Рабочий стол/c mac/Documents/Лаборатория/New_amyloids_coaggregation/Nucleoporines/Tree NSP1, NUP100/NSP1/NSP1_table2.tsv', quote = F, row.names = F, sep= '\t' )

#Загружаем файл, полученный после добавления к таксономической таблице ID белков
taxa_list<-read.table(file = 'NSP1_tax_ID.tsv', sep = '\t', as.is = T, header = T)

#Соединяем две таблицы по сходному столбцу
gene_taxa_list <- merge(x = gene_list, y = taxa_list, by.x = c('seq_name'), by.y=c('seq_ID'))

#Так как в дереве ID вида и ID белка пишутся вместе (даже при наличии разделителя в фасте) - необходимо получить аналогичный формат в таблице со списками видовых названий
gene_taxa_list$species_ID<-paste0(gene_taxa_list$seq_ID, '.', gene_taxa_list$seq_name)

#Сортируем согласно последовательности ID в таблице с видами для последующей замены)
gene_taxa_list <- gene_taxa_list[match(tree$tip.label, gene_taxa_list$species_ID),]

#Заменяем tip с ID на видовые названия
tree$tip.label <- gene_taxa_list$species.y

#Построение дерева (нельзя менять ширину ветвей, пришлось увеличить кол-во цветов для учитывания всех классов)
plot(tree, cex = 0.2, 
     tip.color = rainbow(60)[as.numeric(as.factor(gene_taxa_list$class))], 
     use.edge.length = F, no.margin = F, edge.width=0.1)

#Добавить легенду по классам
legend(x="topleft",legend=sort(unique(as.factor(gene_taxa_list2$class))), fill=unique(rainbow(60)), pt.cex=0.05, cex=0.3, pt.lwd=0.05, bty='n', ncol=2, y.intersp=0.8)
#Лучше сохранять в svg и мелкие изменения делать в редакторе векторной графики inscape/addube illustrator 