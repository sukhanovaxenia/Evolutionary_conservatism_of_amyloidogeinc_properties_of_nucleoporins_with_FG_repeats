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

#Testing for Nup1:
##R1:
position_to_column('ali_noname_R1.fa','YOR098C',coord = c(1,316,322))
xmin=c(1,444,450)
position_to_column('ali_noname_R1.fa','YOR098C',coord = c(340,355,1076))
xmax=c(484,504,1529)

##R2,R3,R5,R6,R7,R8,R10:
#is missed as YOR098C doesn't exist in this alignment
##R4:
position_to_column('ali_noname_R4.fa','YOR098C',coord = c(1,316,322))
xmin=c(1,390,396)
position_to_column('ali_noname_R4.fa','YOR098C',coord = c(340,355,1076))
xmax=c(422,440,1491)

##R9:
position_to_column('ali_noname_R9.fa','YOR098C',coord = c(1,316,322))
xmin=c(1,458,487)
position_to_column('ali_noname_R9.fa','YOR098C',coord = c(340,355,1076))
xmax=c(505,540,1497)

#Testing for Nup57:
#R6(start+end):
position_to_column('ali_noname_R6.fa','XP_006691614.1',coord = c(74,265))
xmin=c(452,761)
position_to_column('ali_noname_R6.fa','XP_006691614.1',coord = c(317,319))
xmax = c(842,844)

#Testing for NUP58:
#R9(start+end)
position_to_column('ali_noname_R9.fa','ENSRNOP00000017204',coord = c(327,415))
xmin=c(618)
xmax=c(706)

#Testing for NUP2:
#R2:
position_to_column('ali_noname_R2.fa','YLR335W',coord = c(1,51))
xmin = c(1)
xmax = c(155)
#R7:
position_to_column('ali_noname_R7.fa','YLR335W',coord = c(1,51))
xmin = c(1)
xmax = c(184)

#Testing for Nup42:
#R1:
position_to_column('ali_noname_R1.fa','YDR192C',coord = c(397,430))
xmin = 818
xmax = 860

#R3:
position_to_column('ali_noname_R3.fa','YDR192C',coord = c(397,430))
xmin = 824
xmax = 872

#R4:
position_to_column('ali_noname_R4.fa','YDR192C',coord = c(397,430))
xmin = 789
xmax = 831

#R5:
position_to_column('ali_noname_R5.fa','YDR192C',coord = c(397,430))
xmin = 792
xmax = 840

#Testing for Nup49:
#R1:
position_to_column('ali_noname_R1.fa','XP_006693639.1',coord = c(246,470))
xmin = 703
xmax = 1046

#R9:
position_to_column('ali_noname_R9.fa','XP_006693639.1',coord = c(246,470))
xmin = 544
xmax =926

#R10:
position_to_column('ali_noname_R10.fa','XP_006693639.1',coord = c(246,470))
xmin = 489
xmax = 867

#Testing for Nup50:
#R1:
position_to_column('ali_noname_R1.fa','ENSP00000345895',coord = c(1,351))
xmin = c(1,759)
position_to_column('ali_noname_R1.fa','ENSP00000345895',coord = c(109,468))
xmax = c(225,881)

#R6:
position_to_column('ali_noname_R6.fa','ENSP00000345895',coord = c(1,351))
xmin = c(1,796)
position_to_column('ali_noname_R6.fa','ENSP00000345895',coord = c(109,468))
xmax = c(376,923)

#Testing for Nup54:
#R6:
position_to_column('ali_noname_R6.fa','ENSP00000264883',coord = c(1,453))
xmin = c(1,773)
position_to_column('ali_noname_R6.fa','ENSP00000264883',coord = c(491,507))
xmax = c(832,862)

#Testing for Nup153:
#R2:
position_to_column('ali_noname_R2.fa','ENSP00000444029',coord = c(722,773,851,1407))
xmin = c(1401,1493,1628,2512)
position_to_column('ali_noname_R2.fa','ENSP00000444029',coord = c(761,822,890,1429))
xmax = c(1481,1599,1721,2541)

#Testing for Nup159:
##No matches with IDs

#Testing for Nsp1:
#R4:
position_to_column('ali_noname_R4.fa','YJL041W',coord = c(496,608))
#642 802
xmin = c(642)
xmax = c(802)

#R8:
position_to_column('ali_noname_R8.fa','YJL041W',coord = c(496,608))
#689 880
xmin = c(689)
xmax = c(880)
