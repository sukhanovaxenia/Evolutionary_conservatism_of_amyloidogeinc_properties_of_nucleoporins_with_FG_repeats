setwd('/media/DATA/Department/laboratory/students/Xeniya/nucleoporins/Nsp1_ncbi_trees')
This script allows to build trees for phylum and also separetly for fungi and metazoa - in this research the latter was used for Nsp1.
library(ape)
#Загружаем таблицу со всеми данными
tax<-read.table(file='table_base_ID.tsv', sep='\t', header = T)
#построить гистограммы для доли амилоидогенных
hist_classes<-hist(tax$Proportion.of.amyloids.class, breaks = 10, main = "Гистограмма доли потенциальных амилоидов разных классов", xlab = 'Доля потенциальных амилоидов (ПА)', ylab = 'Количество классов', las = 1)

# As the separating of the input table was unavailable (the was the wrong result and the separator was recognized between each letter - don't know why), so wrote the new table from the input file:
tax1<-read.table(file='table_base_ID.tsv', sep='\t', header = T)
# chose right variables
tax_fungi<-subset(tax1, kingdom == 'Fungi', select = c(kingdom, phylum, class, N.protein.class, round(Proportion.of.amyloids.class,2)))
#removed duplications in classes
tax_fungi<-unique(tax_fungi)

#install.packages('viridis)  ##the package for a color pallete
library("viridis") 
# make a new variable of color IDs
library(RColorBrewer)
ColorSpace <- plasma(11 + 5)[1:11]
tax_fungi$color <- ColorSpace[as.integer(tax_fungi$Proportion.of.amyloids.class * 10 + 1)]

text_tree_fungi <- '' # create empty string for the tree
for(k in unique(tax_fungi$kingdom)){
  # collect phylum names for certain kingdom
  phylum <- unique(tax_fungi$phylum[tax_fungi$kingdom == k])
  kingdom_string <- '' # create empty string for the kingdom
  for(p in phylum){
    # collect classes names for certain phylum
    classes <- unique(tax_fungi$class[tax_fungi$phylum == p])
    # transform vector of names to the string
    phylum_string <- paste(classes, collapse = ',')
    # remove the comma at the end of the string
    #phylum_string <- substr(phylum_string,1,nchar(phylum_string)-1)
    # add the name of the phylum
    phylum_string <- paste0('(',phylum_string,')',p)
    # add the string to the string describing kingdom
    kingdom_string <- paste0(kingdom_string, phylum_string,',')
  }
  # remove last comma
  kingdom_string <- substr(kingdom_string,1,nchar(kingdom_string)-1)
  # add the name of the kingodm
  kingdom_string <- paste0('(',kingdom_string,')',k)
  # append string of the tree
  text_tree_fungi <- paste0(text_tree_fungi,kingdom_string,',')
}

# remove the last comma, add brackets and  semicolon
text_tree_fungi <- substr(text_tree_fungi,1,nchar(text_tree_fungi)-1)
text_tree_fungi <- paste0('(',text_tree_fungi,')',';')

## in ape:
tree_fungi <- read.tree(text = text_tree_fungi)

ordered_tax_fungi <- tax_fungi[match(tree_fungi$tip.label, tax_fungi$class),]

plot(tree_fungi, tip.color = ordered_tax_fungi$color, main = 'Fungi', show.node.label = T, edge.width=0.2, node.depth=2, node.pos=1, cex=0.8)
legend(x = "bottomright", legend = paste(seq(0,0.9,0.1), seq(0.1,1.0,0.1),sep = '-'), col = ColorSpace, pch = 19,
       title = 'Доля ПА')


# chose right variables
tax_metazoa<-subset(tax1, kingdom == 'Metazoa', select = c(kingdom, phylum, class, N.protein.class, round(Proportion.of.amyloids.class,2)))
#removed duplications in classes
tax_metazoa<-unique(tax_metazoa)

#install.packages('viridis)  ##the package for a color pallete
library("viridis")
library('RColorBrewer')
ColorSpace_m <- plasma(11 + 5)[1:11]
tax_metazoa$color <- ColorSpace_m[as.integer(tax_metazoa$Proportion.of.amyloids.class * 10 + 1)]

# make a new variable of color IDs

text_tree_metazoa <- '' # create empty string for the tree
for(k in unique(tax_metazoa$kingdom)){
  # collect phylum names for certain kingdom
  phylum <- unique(tax_metazoa$phylum[tax_metazoa$kingdom == k])
  kingdom_string <- '' # create empty string for the kingdom
  for(p in phylum){
    # collect classes names for certain phylum
    classes <- unique(tax_metazoa$class[tax_metazoa$phylum == p])
    # transform vector of names to the string
    phylum_string <- paste(classes, collapse = ',')
    # remove the comma at the end of the string
    #phylum_string <- substr(phylum_string,1,nchar(phylum_string)-1)
    # add the name of the phylum
    phylum_string <- paste0('(',phylum_string,')',p)
    # add the string to the string describing kingdom
    kingdom_string <- paste0(kingdom_string, phylum_string,',')
  }
  # remove last comma
  kingdom_string <- substr(kingdom_string,1,nchar(kingdom_string)-1)
  # add the name of the kingodm
  kingdom_string <- paste0('(',kingdom_string,')',k)
  # append string of the tree
  text_tree_metazoa <- paste0(text_tree_metazoa,kingdom_string,',')
}

# remove the last comma, add brackets and  semicolon
text_tree_metazoa <- substr(text_tree_metazoa,1,nchar(text_tree_metazoa)-1)
text_tree_metazoa <- paste0('(',text_tree_metazoa,')',';')

##in ape:
tree_metazoa <- read.tree(text = text_tree_metazoa)

ordered_tax_metazoa <- tax_metazoa[match(tree_metazoa$tip.label, tax_metazoa$class),]

plot(tree_metazoa, tip.color = ordered_tax_metazoa$color, main = 'Metazoa', show.node.label = T, edge.width=0.2, node.depth=2, node.pos=1, cex = 0.8)
legend(x = "bottomright", legend = paste(seq(0,0.9,0.1), seq(0.1,1.0,0.1),sep = '-'), col = ColorSpace, pch = 19,
       title = 'Доля ПА')

##Tree for phylums:
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

tax$class <- paste(tax$class, tax$N.protein.class, round(tax$Proportion.of.amyloids.class,2),
                   sep = '|')
tax$phylum<-paste(tax$phylum, tax$N.proteins.phylum, round(tax$Proportion.of.amyloids.phylum,2),
                  sep = '|')
str(tax)

table_for_tree <- tax[,2:4]
str(table_for_tree)
table_for_tree <- distinct(table_for_tree)

text_tree <- '' # create empty string for the tree
for(k in unique(table_for_tree$kingdom)){
  # collect phylum names for certain kingdom
  phylum <- unique(table_for_tree$phylum[table_for_tree$kingdom == k])
  kingdom_string <- '' # create empty string for the kingdom
  for(p in phylum){
    # collect classes names for certain phylum
    classes <- unique(table_for_tree$class[table_for_tree$phylum == p])
    # transform vector of names to the string
    phylum_string <- paste(classes, collapse = ',')
    # remove the comma at the end of the string
    #phylum_string <- substr(phylum_string,1,nchar(phylum_string)-1)
    # add the name of the phylum
    phylum_string <- paste0('(',phylum_string,')',p)
    # add the string to the string describing kingdom
    kingdom_string <- paste0(kingdom_string, phylum_string,',')
  }
  # remove last comma
  kingdom_string <- substr(kingdom_string,1,nchar(kingdom_string)-1)
  # add the name of the kingodm
  kingdom_string <- paste0('(',kingdom_string,')',k)
  # append string of the tree
  text_tree <- paste0(text_tree,kingdom_string,',')
}

# remove the last comma, add brackets and  semicolon
text_tree <- substr(text_tree,1,nchar(text_tree)-1)
text_tree <- paste0('(',text_tree,')',';')
tree<-read.tree(text=text_tree)
plot(tree, cex = 0.2, 
     tip.color = rainbow(77)[as.numeric(as.factor(tree$node.label))], #When trying to mark according to the node.labes have a trouble in the senc
     use.edge.length = F, no.margin = F, edge.width=0.1)

tree2<-tree$node.label[-1]
p<-ggtree(tree,branch.length = 'none', size = 0.25, layout = 'daylight')
p + geom_nodelab(geom='label', size = 2) + geom_tiplab(geom='text', color = rainbow(78)[as.numeric(as.factor(tree$tip.label))], size = 2)


## FUNGI
#class_color<-as.vector(tax_fungi$color)

#set the list of labels for color setting
#branches<-list('0.641' =1, '0.848' = 2, '0.971' = 3, '0.0' = 4, '0.77' = 5, '1.0' = 6, '0.75' = 7, '1.0' = 8, '0.5' = 9, '1.0' = 10, '0.75' = 11, '0.67' = 12, '0.67' = 13, '0.0' = 14, '0.0' = 15, '0.0' = 16, '1.0' = 17, '1.0' = 18, '0.0' = 19, '0.0' = 20, '1.0' = 21)
#tree_fungi<-groupOTU(tree_fungi, branches)

# circular view
#p_fungi<-ggtree(tree_fungi,branch.length = 'none', size = 0.25, layout = 'fan')
#p_fungi + geom_tiplab(aes(color=group), size = 2) + scale_color_manual(values = class_color) + geom_nodelab(geom = 'label', size = 2)

# bird view
#p_fungi<-ggtree(tree_fungi,branch.length = 'none', size = 0.25, layout = 'daylight')
#p_fungi + geom_tiplab(aes(color=group), size = 2) + scale_color_manual(values = class_color) + geom_nodelab(geom = 'label', size = 2)

#default phylogram view
#p_fungi<-ggtree(tree_fungi,branch.length = 'none', size = 0.25) 
#p_fungi + geom_tiplab(aes(color=group), size = 2) + scale_color_manual(values = class_color) + geom_nodelab(geom = 'label', size = 2)


## METAZOA
#tree_metazoa<-read.tree(text=text_tree_metazoa)
#class_color_m<-as.vector(tax_metazoa$color)

#set the list of labels for color setting
#branches_m<-list()
#for (i in 1:nrow(tax_metazoa)){
#  branches_m[i]<-as.character(round(tax_metazoa[i,5],2))
#  i+1
#}
#tree_metazoa<-groupOTU(tree_metazoa, branches_m)

# circular view
#p_metazoa<-ggtree(tree_metazoa,branch.length = 'none', size = 0.25, layout = 'fan')
#p_metazoa + geom_tiplab(aes(color=group), size = 2) + scale_color_manual(values = class_color_m) + geom_nodelab(geom = 'label', size = 2)

# bird view
#p_metazoa<-ggtree(tree_metazoa,branch.length = 'none', size = 0.25, layout = 'daylight')
#p_metazoa + geom_tiplab(aes(color=group), size = 2) + scale_color_manual(values = class_color_m) + geom_nodelab(geom = 'label', size = 2)

#default phylogram view
#p_metazoa<-ggtree(tree_metazoa,branch.length = 'none', size = 0.25) 
#p_metazoa + geom_tiplab(aes(color=group), size = 2) + scale_color_manual(values = class_color_m) + geom_nodelab(geom = 'label', size = 2)
