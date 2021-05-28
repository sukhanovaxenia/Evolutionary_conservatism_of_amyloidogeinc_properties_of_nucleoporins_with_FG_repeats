#The example for Nup1:
#Step one - get IDs which then will be used for tax_report.txt and phyliptree.phy achieving:
source(file = 'C:/Users/sukha/OneDrive/Документы/AmyloEvo/scripts/get_ID.R', encoding = 'UTF8')
get_ID(input = 'Nup1_opisto_ortho.fa')

#Step two - getting the file with matched taxonomy and sequence IDs:
source(file = 'C:/Users/sukha/OneDrive/Документы/AmyloEvo/scripts/tax_ID.R', encoding = 'UTF8')
tax_ID(input = 'Nup1_opisto_ortho.fa', tax_list = 'tax_report.txt', tax_table = 'result_table_all.tsv', tax_w_ID = 'Nup1_tax_ID.tsv')

#Start the analysis for 10 times:
source(file = 'C:/Users/sukha/OneDrive/Документы/AmyloEvo/scripts/MAIN_function.R', encoding = 'UTF8')
MAIN_fun()

#Recoordinate the structure domains positions and moderate plots of the exact subsets:
source(file = 'C:/Users/sukha/OneDrive/Документы/AmyloEvo/scripts/str_function.R', encoding = 'UTF8')
#R1
position_to_column('ali_noname_R1.fa','YOR098C',coord = c(1,316,322))
#xmin=c(1,444,450)
position_to_column('ali_noname_R1.fa','YOR098C',coord = c(340,355,1076))
#xmax=c(484,504,1529)
str_plot(sumfile = 'Summary_wNA1.tsv', start_coord = c(1,444, 450), end_coord = c(484, 504, 1529), plotfile = 'AA_plot_str1.pdf')


#R4:
position_to_column('ali_noname_R4.fa','YOR098C',coord = c(1,316,322))
#xmin=c(1,390,396)
position_to_column('ali_noname_R4.fa','YOR098C',coord = c(340,355,1076))
#xmax=c(422,440,1491)
str_plot(sumfile = 'Summary_wNA4.tsv', start_coord = c(1,390,396), end_coord = c(422, 440, 1491), plotfile = 'AA_plot_str4.pdf')

#R9:
position_to_column('ali_noname_R9.fa','YOR098C',coord = c(1,316,322))
#xmin=c(1,458,487)
position_to_column('ali_noname_R9.fa','YOR098C',coord = c(340,355,1076))
#xmax=c(505,540,1497)
str_plot(sumfile = 'Summary_wNA9.tsv', start_coord = c(1,458,487), end_coord = c(505, 540, 1497), plotfile = 'AA_plot_str9.pdf')

#Get taxonomy summary and build phylum tree for the whole multifasta file:
source(file = 'C:/Users/sukha/OneDrive/Документы/AmyloEvo/scripts/tree_creation.R', encoding = 'UTF8')
table_tree('Nup1_tax_ID.tsv')


