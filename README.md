# Evolutionary_conservatism_of_amyloidogeinc_properties_of_nucleoporins_with_FG_repeats
This repository includes used in the research work .R scripts and testing files.
All scripts were designed for R version 4.0.5 (2021-03-31) and RStudio version 1.4.1106.

Scripts (in the order of the analysis):

PART 1. Rates of potential amyloids among different taxonomic groups.

1.1. Data preprocessing

1.1.1. Extraction of sequences IDs from fasta file, downloaded from the EggNOG database. 

get_ID - this script allows to get IDs of all orthologs sequences in the input alignment, the output format - .txt
*input file - *_opisto_ortho.fa - the multiple alignment of the protein's orthologs from Opistokhonta group
*output file - 'taxid.txt' - the .txt list of species' taxonomy IDs

1.1.2. Addtion of taxonomy (kingdom, phylum, class) to the IDs

tax_ID - this script is used to match full filled taxonomy (filling is made by hand) with sequences' IDs (without taxonomy IDs)
*input files: input - '*_opisto_ortho.fa', tax_list - 'tax_report.txt', tax_table - 'result_table_all.tsv'
*output file - '*_tax_ID.tsv', consisting of 6 columns: kingdom, phylum, class, order, seq_ID (unique sequence ID) and species
**tax_report.txt is loaded using taxid.txt from the resource: https://www.ncbi.nlm.nih.gov/Taxonomy/TaxIdentifier/tax_identifier.cgi
**result_table_all.tsv is table obtained by parsing of newik tree, downloaded from NCBI Taxonomy database, which was used for all Nups

1.1.3. Optional. Replacement of species IDs in fasta file with species names.
ID_to_species_short_script - is used to change sequences' IDs into species names
**this script isn't obligatory, however was necessary for the manual check of the sequences' quality and tree creation
*fasta - '*_opisto_ortho.fa', the default multiple alignment 
*tax - 'tax_report.txt', loaded from https://www.ncbi.nlm.nih.gov/Taxonomy/TaxIdentifier/tax_identifier.cgi NCBI species list of sequences of the analyzing alignment
*out - '.fa', the new multiple alignment file where all ids are changed to species names

1.2. Calculation rates of potential amyloids among different taxonomic groups.

1.2.1. tree_creation - this script is used to build phylums trees and get taxonomy statistics of potential amyloids for kingdom, phylum and class
input:
taxid - '*_tax_ID.tsv, input file of the full taxonomy: kingdom, phylum, class, order, seq_ID, species
outputs:
'table_base.tsv' - table with taxonomy statistics of PA without sequences' IDs
the phylogenetic tree with colorcoded rates of potential amyloids


PART 2. Search of regions with conservative amyloid properties

2.1. The MAIN function allows to take a ramdom subset of sequenses, perform two-step multiple alignment, rearrange the amyloidogenic score (Cumulative Score from ArchCandy) and calculate different statistics for each position of the alignment.
The MAIN function launch three subfunctions listed below for 10 times.

2.1.1. muscle_refine --- this function performs filtering and two-step aligning of our multifasta 

input:

fastaFile - default multifasta with orthologs sequences (*_opisto_ortho.fa)

file - result_table_all.tsv file gained by Danilov Lavrentyii's python script

tree - full taxonomy table

input (optional) - optional multifile file where sequence IDs are replaced with species' names

filtered - filename for filtered multiple fasta (OBLIGATORY)

output:

outputpath - file name for filtered multiple alignment with sequence labeles 'species's name.sequence ID'

outputpath_noname - file name for filtered multiple alignment with sequence labeles 'sequence ID'

2.1.2. aa_stat --- this function allows to get the amino acid composition:

input:

fasta - alignment file, gained with R function with filtering (filter plus alignment)

output:

stat - filename of the table with aa types proportions through alignment (for each position)

plot - filename of the proportions' stat graphic

2.1.3. CS_subset --- this function count subsetted cumulative scores, aa composition and amyloidogenicy parameteres and builds summary and comparitive plots

input:

alignment - all previously gotten alignment for each repetition are taken by this function;

cumsore - gotten by ArchCandy umulative scores for each sequence in the default orthologs' dataset (one for all repetitions)

taxonomy - right taxonomy gotten from NCBI db for the orthologs' subset (one for all repetitions)

id - the list of ID's for the analyzing protein (one for all repetitions)

output:

cs_output - cumulative scores for each repetitions;

aa_stat - plots of analyzed parametres for each position in each repetition without NAs

aa_stat_wNA - plots of analyzed parametres for each position in each repetition with NAs

summary - summary of analyzed parametres for subsets withot NAs

CS_comp_plot - comparison plots of amyloidogenicy with and without gaps for each repetition

summary_wNA - summary of parametres for subsets with NAs


MAIN_function - this script unites three previous functions and launch them without any arguments 10 times by your directory (necessary files should be in this directory: 
CumScores.tsv, *_opisto_ortho.fa, *_nog_orthologs.txt, result_table_all.tsv, tax_report.txt, taxid.txt, phyliptree.phy, *_tax_ID.tsv

inputs:

'CumScores.tsv', gotten by ArchCandy cumulative scores for each sequence in the default orthologs' dataset (one for all repetitions)

'*_opisto_ortho.fa', fasta with multiple alignment of orthologs' sequences

'*_nog_orthologs.txt', the list of ID's for the analyzing protein (one for all repetitions)

'result_table_all.tsv', right taxonomy gotten from NCBI db for the orthologs' subset (one for all repetitions)

'tax_report.txt', the list of species and species IDs included into the alignment

'taxid.txt', the .txt list of species' taxonomy IDs

'phyliptree.phy', the .phy taxonomy common tree, including kingdom, phylum, class, order and species (one for all repetitions)

'*_tax_ID.tsv', file of the full taxonomy: kingdom, phylum, class, order, seq_ID, species

outputs:

'ali_noname_R*.fa - filtered aligment with only IDs and sequences

'filt_noname_R*.fa' - filtered multifasta

'ali_stat*.tsv' - amino acid statistics

'ali_aa_stat*.pdf' - file with the proportions' stat graphic

'CS_subset*.tsv' - CumScores parametres for each subset

'AA_plot*.pdf' - plots of gained parametres (not including NAs) for all alignments' positions

'AA_plot_wNA*.pdf' - plots of gained parametres (including NAs) for all alignments' positions

CS_comp_plot*.pdf' - comparison of amyloidogenicy behaviour depending on the gaps distribution for subsets without anf with NAs

'Summary_stat*.tsv' - summary of all parametres (not including NAs)

'Summary_wNA*.isv' - summary of all parametres (including NAs)


2.2. Mapping of structured regions on alignment

2.2.1. str_function - this script allows to reccordinate structure domain coordinates according to each sample alignment and build new plots with painted structure domains

position_to_column (function, the output is shown in the console):
inputs:
*alignment - 'ali_noname_R*.fa',the default alignment of the exact repetition; 
*pat = e.g. 'YOR098C' (for  Nup1), the NCBI ID of the species for which structure were found (ID is written in the file *_str_dom.txt, where * - one of nucleoporins); 
*coords - e.g. c(1,316,322) (for Nup1), vector of coordinates for startings and endings separetly (positions are in the file *s_str_dom.txt.

2.2.2. str_plot (function)
inputs:
*sumfile - 'Summary_wNA*.tsv', where * - the number of the summary for the repetition with included ID; the file with parametres by positions (Summary_wNA*.tsv, where * - the number of the summary for the repetition with included ID)
*start_coord - vector of starting coordinates of domains
*end_coord - vector of ending coords of domains
output:
*plotfile = 'AA_plot_str*.pdf',the name of the new plot;
  
  
========================================================================

The list of used packages and versions:

dplyr - 1.0.5

plyr - 1.8.6

tidyr - 1.1.3

Biostrings - 2.58.0

scales - 1.1.1

ggplot2 - 3.3.3

stringr - 1.4.0

muscle - 3.32.0

seqinr - 4.2.5

ape - 5.5

viridis - 0.6.1

RColorBrewer - 1.1.2

