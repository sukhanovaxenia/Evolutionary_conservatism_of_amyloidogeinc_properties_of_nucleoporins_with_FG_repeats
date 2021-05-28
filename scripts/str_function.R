library(tidyr)
library(dplyr)
library(ggplot2)
library(muscle)
#This script allows to recoordinate found structured domains' coordinates for input protein on the each repetition alignment.
##It's needed to "normalize" structured domains' positions on gotten earlier parameteres of the datasets.

#Recoordinate structure positions on the default alignment of the repetition. Only datasets which include species's ID for which structure domains were found:
##new coordinates are counted for the starting ones and the endings separatly (coorinates are in the file Protein_str_dom.txt, where R* - is the number of repetitions)
## alignment - the default alignment of the exact repetition; 
## pat = the NCBI ID of the species for which structure were found (ID is written in the file Protein_str_dom.txt); 
## coords - vector of coordinates for startings and endings separetly.
position_to_column<- function(alignment, pat, coord){
  ali<-readAAStringSet(alignment, format = 'fasta')
  string <- as.character(ali[grepl(pat, ali@ranges@NAMES)])
  print(string)
  columns = c()
  column = 1
  tmp_coord = 1
  for(i in coord){
    while(tmp_coord != i & column < nchar(string)){
      symbol = substr(x = string, start = column, stop = column)
      if(symbol != '-'){tmp_coord = tmp_coord + 1}
      column = column + 1
      print(c(symbol, tmp_coord, column))
    }
    columns = c(columns, column)
  }
  return(columns)
}
#Achieving new AA_plots for filled with structure domains summary data. Coordinates should be found by yourself:
## sumfile - the file with parametres by positions (Summary_wNA*.tsv, where * - the number of the summary for the repetition with included ID)
## plotfile = the name of the new plot in the format AA_plot_str*.pdf
## start_coord - vector of starting coordinates of domains
## end_coord - vector of ending coords of domains
str_plot<-function(sumfile,plotfile,start_coord,end_coord, ..deparseOpts){
  stat<-sumfile
  aa_stplot = plotfile
  domain_stat<-read.table(file = stat, sep = '\t', header = T)
  ggplot(domain_stat,mapping = aes(x = position, y = name, fill = value)) + geom_tile(height=.9) + 
    scale_fill_viridis_c(name = "",limits = c(0.0,1.0),oob=squish) +  theme_bw()+
    annotate('rect',
             xmin=start_coord,
             xmax=end_coord, 
             ymin = 0.5, ymax=6.5, fill=NA, color='red')
  ggsave(filename = aa_stplot,device = 'pdf', width = 20, height = 4, units = 'cm')
}


