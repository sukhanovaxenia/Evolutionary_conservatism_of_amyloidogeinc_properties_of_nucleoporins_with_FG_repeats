# Evolutionary_conservatism_of_amyloidogeinc_properties_of_nucleoporins_with_FG_repeats
This repository includes used in the research work .R scripts and testing files.
All scripts were designed for R version 4.0.5 (2021-03-31) and RStudio version 1.4.1106.

Scripts (in the order of the analysis):

1. Step one: get_ID - this script allows to get IDs of all orthologs sequences in the input alignment, the output format - .txt
*input file - *_opisto_ortho.fa - the multiple alignment of the protein's orthologs from Opistokhonta group
*output file - 'taxid.txt' - the .txt list of species' taxonomy IDs

2. Step two: tax_ID - this script is used to match full filled taxonomy (filling is made by hand) with sequences' IDs (without taxonomy IDs)

*input files: input - '*_opisto_ortho.fa', tax_list - 'tax_report.txt', tax_table - 'result_table_all.tsv'

*output file - '*_tax_ID.tsv', consisting of 6 columns: kingdom, phylum, class, order, seq_ID (unique sequence ID) and species

**tax_report.txt is loaded using taxid.txt from the resource: https://www.ncbi.nlm.nih.gov/Taxonomy/TaxIdentifier/tax_identifier.cgi

**result_table_all.tsv is gained by a special private python script, written by Danilov Lavrentyii, and full filled by hand according to the NCBI Taxonomy database

3. ID_to_species_short_script - is used to change sequences' IDs into species names

**this script isn't obligatory, however was necessary for the manual check of the sequences' quality and tree creation

*fasta - '*_opisto_ortho.fa', the default multiple alignment 

*tax - 'tax_report.txt', loaded from https://www.ncbi.nlm.nih.gov/Taxonomy/TaxIdentifier/tax_identifier.cgi NCBI species list of sequences of the analyzing alignment

*out - '.fa', the new multiple alignment file where all ids are changed to species names

4. filter_plus_alignment - this script includes the function to filter our sequences according to criterias and make automatical MUSCLE alignment

inputs:

*fastaFile - '*_opisto_ortho.fa', fasta with multiple alignment of orthologs' sequences
 
*file - 'tax_report.txt', the list of species and species IDs included into the alignment

*tree - '*_tax_ID.tsv', file of the full taxonomy: kingdom, phylum, class, order, seq_ID, species

outputs:

*filtered - 'filt.fa', filtered multiple alignment file

*outputpath - 'ali_R.fa', MUSCLE alignment filtered file where each sequence is named by 'seq_ID.species'

*outputpath_noname - 'ali_noname_R.fa', MUSCLE alignment filtered file with only seq_IDs

5. ali_aa_statistics_plot_R - this script allows to get the alignments' statistic(amino acid and gap's composition) and illustrate them with ggplot2 package

input:

*fasta - 'ali_R.fa', alignment file, gained with R function with filtering (filter plus alignment)

outputs:

*stat - 'ali_stat.tsv', filename of the table with aa types proportions through alignment (for each position)

*plot - 'ali_aa_stat.pdf', filename of the proportions' stat graphic

6. CS_subset_final - this scripts is used to count amino acid and gap composition(another algorithm), amyloidogenicy, amino acid conservaism; illustrate the parametres distributions across the alignment; build parametres' tables

inputs:

*alignment - 'ali_noname_R.fa', all previously gotten alignment for each repetition are taken by this function;

*cumsore - 'CumScores.tsv', gotten by ArchCandy umulative scores for each sequence in the default orthologs' dataset (one for all repetitions)

*taxonomy - 'result_table_all.tsv', right taxonomy gotten from NCBI db for the orthologs' subset (one for all repetitions)

*id - '*_nog_orthologs.txt', the list of ID's for the analyzing protein (one for all repetitions)

outputs:

*cs_output - 'CS_subset.tsv', cumulative scores for each repetitions;

*PA_output - 'Report_table.tsv', number of amyloidogenic sequences for big taxonomy groups of each repetition

*aa_stat - 'AA_plot.pdf', plots of analyzed parametres for each position in each repetition

*aa_stat_wNA - 'AA_plot_wNA.pdf', plots of analyzed parametres for each position in each repetition, ncluding NAs

*summary - 'Summary_stat.tsv', summary of analyzed parametres for subsets withot NAs

*CS_comp_plot - 'CS_comp_plot.pdf', comparison plots of amyloidogenicy with and without gaps for each repetition

*summary_wNA - 'Summary_wNA.tsv', summary of parametres for subsets with NAs  

*comp - 'Comparison_cs.tsv', comparison of the statistics with both included NA and not postitions

7. Step three: MAIN_function - this script unites three previous functions and launch them without any arguments 10 times by your directory (necessary files should be in this directory: 
CumScores.tsv, *_opisto_ortho.fa, *_nog_orthologs.txt, result_table_all.tsv, tax_report.txt, taxid.txt, phyliptree.phy, *_tax_ID.tsv

inputs:

*'CumScores.tsv', gotten by ArchCandy umulative scores for each sequence in the default orthologs' dataset (one for all repetitions)

*'*_opisto_ortho.fa', fasta with multiple alignment of orthologs' sequences

*'*_nog_orthologs.txt', the list of ID's for the analyzing protein (one for all repetitions)

*'result_table_all.tsv', right taxonomy gotten from NCBI db for the orthologs' subset (one for all repetitions)

*'tax_report.txt', the list of species and species IDs included into the alignment

*taxid.txt, the .txt list of species' taxonomy IDs

*'phyliptree.phy', the .phy taxonomy common tree, including kingdom, phylum, class, order and species (one for all repetitions)

*'*_tax_ID.tsv', file of the full taxonomy: kingdom, phylum, class, order, seq_ID, species

outputs:

*'ali_noname_R*.fa - filtered aligment with only IDs and sequences

*'filt_noname_R*.fa' - filtered multifasta

*'ali_stat*.tsv' - amino acid statistics

*'ali_aa_stat*.pdf' - file with the proportions' stat graphic

*'CS_subset*.tsv' - CumScores parametres for each subset

*'Report_table*.tsv' - number of amyloidogenic sequences for big taxonomy groups of each repetition

*'AA_plot*.pdf' - plots of gained parametres (not including NAs) for all alignments' positions

*'AA_plot_wNA*.pdf' - plots of gained parametres (including NAs) for all alignments' positions

*CS_comp_plot*.pdf' - comparison of amyloidogenicy behaviour depending on the gaps distribution for subsets without anf with NAs

*'Summary_stat*.tsv' - summary of all parametres (not including NAs)

*'Summary_wNA*.isv' - summary of all parametres (including NAs)


8. Step four: str_function - this script allows to reccordinate structure domain coordinates according to each sample alignment and build new plots with painted structure domains

 1) position_to_column (function, the output is shown in the console):
 
   inputs:
   
   *alignment - 'ali_noname_R*.fa',the default alignment of the exact repetition; 

   *pat = e.g. 'YOR098C' (for  Nup1), the NCBI ID of the species for which structure were found (ID is written in the file *_str_dom.txt, where * - one of nucleoporins); 

   *coords - e.g. c(1,316,322) (for Nup1), vector of coordinates for startings and endings separetly (positions are in the file *s_str_dom.txt.
 
 2) str_plot (function)
   
   inputs:
   
   *sumfile - 'Summary_wNA*.tsv', where * - the number of the summary for the repetition with included ID; the file with parametres by positions (Summary_wNA*.tsv, where * - the number of the summary for the repetition with included ID)
   
   *start_coord - vector of starting coordinates of domains
   
   *end_coord - vector of ending coords of domains
   
   output:
   
   *plotfile = 'AA_plot_str*.pdf',the name of the new plot;

9. Step five: tree_creation - this script is used to build phylums trees and get taxonomy statistics of potential amyloids for kingdom, phylum and class
  1) tax_ID (function) is used to gain kingdom, phylum and class statistic of following parametres: number of proteins, proportion of PA
   
   inputs:
   
   *input - '*_opisto_ortho.fa', default multifasta of orthologs
   
   *tax_list - 'tax_report.txt', the list of species and species IDs included into the alignment
   
   *tax_table - 'result_table_all.tsv', right taxonomy gotten from NCBI db for the orthologs' subset (one for all repetitions)
   
   output:
   
   *tax_w_ID - '*_tax_ID.tsv, output file of the full taxonomy: kingdom, phylum, class, order, seq_ID, species
  
  2) table_tree (function):
  
  input:
  
  taxid - '*_tax_ID.tsv, input file of the full taxonomy: kingdom, phylum, class, order, seq_ID, species
  
  otputs:
  
  'table_base.tsv' - table with taxoomy statistics of PA without sequences' IDs
  
  'table_base_ID.tsv' - table with taxoomy statistics of PA with sequences' IDs
  
