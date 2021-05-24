# Evolutionary_conservatism_of_amyloidogeinc_properties_of_nucleoporins_with_FR_repeats
This repository includes used in the research work .R scripts and testing files.

Scripts:
get_ID - this script allows to get IDs of all orthologs sequences in the input alignment, the output format - .txt
Parse_NCBI_taxonomy_v2 - this python script is used to get the taxonomy table. !NB: is should be filled up by hand if for some specieces there is no taxonomy position according to the NCBI database
*phyliptree.phy can be loaded using taxid.txt (fisrt script) from the link:https://www.ncbi.nlm.nih.gov/Taxonomy/CommonTree/wwwcmt.cgi
**!NB: only kingdom, phylum, class, order, species should be chosen
tax_ID - this script is used to match full filled taxonomy (filling is made by hand) with sequences' IDs (without taxonomy IDs)
tax_ID_short - the shorter and faster script as the previous one
*tax_report.txt is loaded using taxid.txt from the resource: https://www.ncbi.nlm.nih.gov/Taxonomy/TaxIdentifier/tax_identifier.cgi
ID_to_species_short_script - is used to change sequences' IDs into species names
*this script isn't obligatory, however was necessary for the manual check of the sequences' quality and tree creation
filter_plus_alignment - this script includes the function to filter our sequences according to criterias and make automatical MUSCLE alignment
ali_aa_statistics_plot_R - this script allows to get the alignments' statistic(amino acid and gap's composition) and illustrate them with ggplot2 package
CS_subset_final - this scripts is used to count amino acid and gap composition(another algorithm), amyloidogenicy, amino acid conservaism; illustrate the parametres distributions across the alignment; build parametres' tables
Super_function - this script unites three previous functions and launch them without any arguments 10 times by your directory (necessary files should be in this directory: CumScores.tsv, *_opisto_ortho.fa, *_nog_orthologs.txt, result_table_all.tsv, tax_report.txt, taxid.txt, phyliptree.phy, *_tax_ID.tsv
str_function - this script allows to reccordinate structure domain coordinates according to each sample alignment and build new plots with painted structure domains
tree_creation - this script is used to build phylums trees and get taxonomy statistics of potential amyloids for kingdom, phylum and class
