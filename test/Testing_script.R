#The example for Nup2:
#Step one - get IDs which then will be used for tax_report.txt and phyliptree.phy achieving:
source(file = 'C:/Users/sukha/OneDrive/Документы/AmyloEvo/scripts/get_ID.R', encoding = 'UTF8')
get_ID(input = 'Nup2_opisto_ortho.fa')

#Step two - getting the file with matched taxonomy and sequence IDs:
source(file = 'C:/Users/sukha/OneDrive/Документы/AmyloEvo/scripts/tax_ID.R', encoding = 'UTF8')
tax_ID(input = 'Nup2_opisto_ortho.fa', tax_list = 'tax_report.txt', tax_table = 'result_table_all.tsv', tax_w_ID = 'Nup2_tax_ID.tsv')

#Start the analysis for 10 times:
source(file = 'C:/Users/sukha/OneDrive/Документы/AmyloEvo/scripts/MAIN_function.R', encoding = 'UTF8')
MAIN_fun()

#Recoordinate the structure domains positions and moderate plots of the exact subsets:
source(file = 'C:/Users/sukha/OneDrive/Документы/AmyloEvo/scripts/str_function.R', encoding = 'UTF8')
#R2:
position_to_column('ali_noname_R2.fa','YLR335W',coord = c(1,51))
#xmin = c(1)
#xmax = c(155)
str_plot(sumfile = 'Summary_wNA2.tsv', start_coord = 1, end_coord = 155, plotfile = 'AA_plot_str2.pdf')

#R7:
position_to_column('ali_noname_R7.fa','YLR335W',coord = c(1,51))
#xmin = c(1)
#xmax = c(184)
str_plot(sumfile = 'Summary_wNA7.tsv', start_coord = 1, end_coord = 184, plotfile = 'AA_plot_str7.pdf')

#Get taxonomy summary and build phylum tree for the whole multifasta file:
source(file = 'C:/Users/sukha/OneDrive/Документы/AmyloEvo/scripts/tree_creation.R', encoding = 'UTF8')
table_tree('Nup2_tax_ID.tsv')


