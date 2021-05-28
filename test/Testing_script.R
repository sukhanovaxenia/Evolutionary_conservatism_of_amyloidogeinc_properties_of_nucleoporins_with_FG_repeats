library(dplyr)
library(plyr)
library(tidyr)
library(Biostrings)

#The example for Nup2:
#Step one - get IDs which then will be used for tax_report.txt and phyliptree.phy achieving:
source(file = '../scripts/get_ID.R', encoding = 'UTF8')
get_ID(input = 'Nup2_opisto_ortho.fa')

#Step two - getting the file with matched taxonomy and sequence IDs:
source(file = '../scripts/tax_ID.R', encoding = 'UTF8')
tax_ID(input = 'Nup2_opisto_ortho.fa', tax_list = 'tax_report.txt', 
       tax_table = 'result_table_all.tsv', tax_w_ID = 'Nup2_tax_ID.tsv')

#Start the analysis for 10 times:
source(file = '../scripts/MAIN_function.R', encoding = 'UTF8')
MAIN_fun()

#Recalculate the structure domains positions and moderate plots of the exact subsets:
source(file = '../scripts/str_function.R', encoding = 'UTF8')
#R2: please check if 'YLR335W' in 'ali_noname_R2.fa'
position_to_column('ali_noname_R2.fa','YLR335W',coord = c(1,51))
# the function two new coordinates of the structured region
str_plot(sumfile = 'Summary_wNA2.tsv', start_coord = 1, end_coord = 210, plotfile = 'AA_plot_str2.pdf')

#Get taxonomy summary and build phylum tree for the whole multifasta file:
source(file = '../scripts/tree_creation.R', encoding = 'UTF8')
table_tree('Nup2_tax_ID.tsv')



