#Parsing ID according to the taxonomy gotten by Lanvrentyi's script and filled up
#Building trees for phylums and histograms of PA in classes:
table_tree<-function(taxid){
# make a new variable of color IDs
ColorSpace <- plasma(11 + 5)[1:11]
cumscores<-read.table(file="CumScores.tsv", header=F, sep="\t")
cumscores1<-cumscores
cumscores1$sum<-rowSums(cumscores[,2:ncol(cumscores)], na.rm = T)
cumscores_merged<-cumscores1[,c(1,ncol(cumscores1))]
cumscores_merged<-separate(cumscores_merged, col=V1, into=c("species_ID", "protein_ID"), extra="merge")

tax1 <- read.csv(file=taxid, header = T, row.names=NULL, fill=T, sep='\t')
a1<-merge(cumscores_merged, tax1, by.x="protein_ID", by.y="seq_ID", sort=F)
##Add new coding variable - is(T)/not(F) amyloidogenic protein (sum>0):
b<-mutate(a1,amy=ifelse(a1$sum>0, T, F))
#Order species according the default variant:
b<-with(b, b[order(species),])
b<-b[,c(2,1,3,9,8,7,6,5,4)]
#Update row numeration
row.names(b)<-NULL

tax1_brief<-b[,c(9,8,7,4,2)]

#Count parametres for taxonomic groups: kingdom, phylum, class:
##number of proteins
##number of potential amyloids (PA)
##proportion of PA
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

#Merging into the full table:  
tax1_brief<-merge(tax1_brief,merge(class_length,class_prop,by="Group.1"), by.x="class", by.y="Group.1")
tax1_brief<-merge(tax1_brief,merge(phylum_length,phylum_prop,by="Group.1"), by.x="phylum", by.y="Group.1")
tax1_brief<-merge(tax1_brief,merge(kingdom_length,kingdom_prop,by="Group.1"), by.x="kingdom", by.y="Group.1")

#Rename colums and export without sequence IDs:
colnames(tax1_brief)[6:11]<-c("N protein class", "Proportion of amyloids class", "N proteins phylum", "Proportion of amyloids phylum", "N proteins kingdom", "Proportion of amyloids kingdom")

write.table(tax1_brief, file="table_base.tsv",quote=F, sep='\t', row.names=F)

#Export table with sequence IDs:
tax1_brief_ID<-merge(b[,1:2], tax1_brief, by.x="protein_ID", by.y="protein_ID")
tax1_brief_ID<-unite(tax1_brief_ID, col = "ID", species_ID, protein_ID, sep = ".")
write.table(tax1_brief_ID, file="table_base_ID.tsv",quote=F, sep='\t', row.names=F)

#Building trees:

#Load table with parametres for taxonomic groups and sequences' IDs (table_base_ID.tsv):
tax<-read.table(file='table_base_ID.tsv', sep='\t', header = T)
#Building hist of PA's proportions for classes:
hist_classes<-hist(tax$Proportion.of.amyloids.class, breaks = 10, main = "Гистограмма доли потенциальных амилоидов разных классов", xlab = 'Доля потенциальных амилоидов (ПА)', ylab = 'Количество классов', las = 1)

# As the separating of the input table was unavailable (the was the wrong result and the separator was recognized between each letter - don't know why), so wrote the new table from the input file:
tax1<-read.table(file='table_base_ID.tsv', sep='\t', header = T)
#Choose right variables
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

#Remove the last comma, add brackets and  semicolon
text_tree_phylum <- substr(text_tree_phylum,1,nchar(text_tree_phylum)-1)
text_tree_phylum <- paste0('(',text_tree_phylum,')',';')

##Build trees in ape:
tree_phylum <- read.tree(text = text_tree_phylum)

ordered_phylum <- phylum[match(tree_phylum$tip.label, phylum$phylum),]

plot(tree_phylum,tip.color = ordered_phylum$color, main = 'Phylum', show.node.label = F, edge.width=0.2,node.depth=2, cex = 0.8)
legend(x = "bottomleft", legend = paste(seq(0,0.9,0.1), seq(0.1,1.0,0.1),sep = '-'), col = ColorSpace, pch = 19,
       title = 'Доля ПА', cex = 0.5)
}
