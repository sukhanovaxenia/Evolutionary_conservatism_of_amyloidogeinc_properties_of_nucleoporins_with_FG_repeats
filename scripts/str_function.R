library(tidyr)
library(dplyr)
library(ggplot2)
library(muscle)
#This script allows to recoordinate found structured domains' coordinates for input protein on the each repetition alignment.
##It's needed to "normalize" structured domains' positions on gotten earlier parameteres of the datasets.

#Get alignment for each sample data with structure domains - realigning structure domains on each repetition dataset of the protein:
muscle_str<-function(x=dir(x, all.files = T), ..deparseOpts){
  i=1
  while (i<=10){
    print(i)
    fasta<-paste('ali_noname_R',i,'_str.fa',sep='')
    ali<-paste('str_R',i,'_ali.fa',sep='')
    structure<-readAAStringSet(file = fasta)
    mu<-muscle(structure)
    mu<-AAStringSet(mu)
    writeXStringSet(mu, filepath = ali)
    i=i+1
  }
}
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
position_to_column('ali_noname_R1.fa','YOR098C',coord = c(54,126,244,253,286,304,339,355,414,454,508,515,542,604,617,645,668,694,705,736,768,791,813,833,862,961,974,999,1028,1067,1099,1130,1158,1201,1221,1256,1289,1391,1418,1470,1532,1557,1616,1636,1670,1690))
xmin=c(49,124,144,170,245,359,416,505,581,624,689,735,767,796,837,903,947,965,995,1021,1067,1116,1197,1229,1287,1325,1336,1356,1416,1481,1502)
xmax=c(103,139,162,223,315,408,474,558,621,658,710,751,784,833,846,943,960,983,1011,1044,1102,1177,1212,1260,1304,1332,1342,1401,1440,1497,1523)

##R2,R3,R5,R6,R7,R8,R10:
#is missed as YOR098C doesn't exist in this alignment
##R4:
position_to_column('ali_noname_R4.fa','YOR098C', coord = c(10,140,350,404,514,591,667,709,751,792,840,958,1051,1079,1140,1183,1225,1249,1275,1320,1387,1426,1455))
xmin=c(23,188,435,497,660,780,889,989,1063,1154,1202,1337,1447)
position_to_column('ali_noname_R4.fa','YOR098C', coord = c(60,315,381,496,579,618,686,740,774,824,912,1008,1065,1088,1161,1205,1238,1264,1311,1359,1400,1438,1468))
xmax=c(88,389,474,619,760,825,936,1038,1086,1186,1288,1390,1478)
##R9:
position_to_column('ali_noname_R9.fa','YOR098C', coord = c(8,127,220,374,432,497,556,709,882,973,1099,1155,1201,1261,1300,1443))
xmin=c(17,219,344,572,661,736,829,1052,1281,1374)
position_to_column('ali_noname_R9.fa','YOR098C', coord = c(56,202,345,404,466,514,612,846,961,1059,1115,1180,1233,1281,1415,1465))
xmax=c(133,326,510,602,695,762,911,1245,1362,1480) 

#Testing for Nup57:
#R6(start+end):
position_to_column('ali_noname_R6.fa','XP_006691614.1',coord = c(771,778,807,854))
#882 882 882 882
##It means that no str domains fit the alignment

#Testing for NUP58:
#R9(start+end)
position_to_column('ali_noname_R9.fa','ENSRNOP00000017204',coord = c(257,260,601,689))
#548 551 887 887
xmin=c(548)
xmax=c(551)

#Testing for NUP2:
#R2:
position_to_column('ali_noname_R2.fa','YLR335W',coord = c(48,134,150,191,809))
#152  597  644  809 2045
xmin = c(152,597,644,809)
position_to_column('ali_noname_R2.fa','YLR335W',coord = c(60,142,171,196,810))
#196  633  710  814 2045
xmax = c(196,633,710,814)
#R7:
position_to_column('ali_noname_R7.fa','YLR335W',coord = c(29,69,1641))
#136  227 1965
xmin = c(136,227)
position_to_column('ali_noname_R7.fa','YLR335W',coord = c(46,107,1642))
#159  268 1965
xmax = c(159,268)

#Testing for Nup42:
#R1:
position_to_column('ali_noname_R1.fa','YDR192C',coord = c(419,795,824))
#849 870 870
xmin = 849
position_to_column('ali_noname_R1.fa','YDR192C',coord = c(422,816,837))
#852 870 870
xmax = 852

#R3:
position_to_column('ali_noname_R3.fa','YDR192C',coord = c(613,833,851,866))
#881 881 881 881
xmin = NULL
position_to_column('ali_noname_R3.fa','YDR192C',coord = c(617,846,858,879))
#881 881 881 881
xmax = NULL

#R4:
position_to_column('ali_noname_R4.fa','YDR192C',coord = c(403,780,807))
#795 841 841
xmin = 795
position_to_column('ali_noname_R4.fa','YDR192C',coord = c(411,798,819))
#803 841 841
xmax = 803

#R5:
position_to_column('ali_noname_R5.fa','YDR192C',coord = c(544,795,810))
#850 850 850
xmin = NULL
position_to_column('ali_noname_R5.fa','YDR192C',coord = c(548,799,839))
#850 850 850
xmax = NULL

#Testing for Nup49:
#R1:
position_to_column('ali_noname_R1.fa','XP_006693639.1',coord = c(604,752,802,868,919))
#1106 1106 1106 1106 1106
xmin = NULL
position_to_column('ali_noname_R1.fa','XP_006693639.1',coord = c(702,791,839,912,931))
#1106 1106 1106 1106 1106
xmax = NULL

#R9:
position_to_column('ali_noname_R9.fa','XP_006693639.1',coord = c(522,605,665,705,736,771,835,850,905))
#993 993 993 993 993 993 993 993
xmin = NULL
position_to_column('ali_noname_R9.fa','XP_006693639.1',coord = c(579,644,683,715,758,793,842,895,908))
#993 993 993 993 993 993 993 993
xmax = NULL

#R10:
position_to_column('ali_noname_R10.fa','XP_006693639.1',coord = c(529,612,652,682,698,725,795,850,872,879,897))
#935 935 935 935 935 935 935 935 935 935 935
xmin = NULL
position_to_column('ali_noname_R10.fa','XP_006693639.1',coord = c(586,647,658,687,713,768,821,867,875,889,917))
#935 935 935 935 935 935 935 935 935 935 935
xmax = NULL

#Testing for Nup50:
#R1:
position_to_column('ali_noname_R1.fa','ENSP00000345895',coord = c(37,117,131,163,193,226,234,267,726))
#80 243 261 344 477 567 581 629 886
xmin = c(80,243,261,344,477,567,581,629)
position_to_column('ali_noname_R1.fa','ENSP00000345895',coord = c(94,121,133,169,195,227,249,287,848))
#198 247 263 350 479 568 609 661 886
xmax = c(198,247,263,350,479,568,609,661)

#R6:
position_to_column('ali_noname_R6.fa','ENSP00000345895',coord = c(73,226,241,279,297,309,353,390,799,896))
#270 605 636 695 739 751 798 836 928 928
xmin = c(270,605,636,695,739,751,798,836)
position_to_column('ali_noname_R6.fa','ENSP00000345895',coord = c(131,235,247,288,298,322,361,408,890,926))
#426 617 642 722 740 765 806 854 928 928
xmax = c(426,617,642,722,740,765,806,854)

#Testing for Nup54:
#R6:
position_to_column('ali_noname_R6.fa','ENSP00000264883',coord = c(113,161,197,211,236,259,291,321,373,403,429,563,628,579,820,893,935))
#348 419 462 478 514 541 577 612 665 697 742 876 876 876 876 876 876
xmin = c(348,419,462,478,541,577,612,665,697,742)
position_to_column('ali_noname_R6.fa','ENSP00000264883',coord = c(121,173,207,212,240,277,298,362,398,413,553,621,667,803,873,920,950))
#362 433 474 479 518 561 584 653 692 707 876 876 876 876 876 876 876
xmax = c(362,433,474,479,518,561,584,653,681,707)

#Testing for Nup153:
position_to_column('ali_noname_R2.fa','ENSP00000444029',coord = c(911,1405,1465,1508,1536,1654,2467))
#1742 2510 2590 2641 2641 2641 2641
xmin = c(1742, 2510, 2590)
position_to_column('ali_noname_R2.fa','ENSP00000444029',coord = c(917,1436,1479,1520,1570,1695,2490))
#1757 2550 2606 2641 2641 2641 2641
xmax = c(1757, 2550, 2606)

#Testing for Nup159:
##No matches with IDs

#Testing for Nsp1:
#R4:
position_to_column('ali_noname_R4.fa','YJL041W',coord = c(143,292,468,496,562,578,712,884,1171,1213,245))
#154 381 614 642 745 772 962
xmin = c(154,381,614,642,745,772,962)
position_to_column('ali_noname_R4.fa','YJL041W',coord = c(146,321,477,515,574,583,717,894,1178,1217,1250))
#157 431 623 674 768 777 967
xmax = c(157,431,623,674,768,777,967)

#R8:
position_to_column('ali_noname_R8.fa','YJL041W',coord = c(126,195,238,285,370,620,816,1050,1093,1120))
#156 235 334 429 529 892
xmin = c(156,235,334,429,529,892)
position_to_column('ali_noname_R8.fa','YJL041W',coord = c(135,220,250,292,381,640,825,1054,1098,1125))
#165 260 369 436 541 912
xmax = c(165,260,369,436,541,912)
