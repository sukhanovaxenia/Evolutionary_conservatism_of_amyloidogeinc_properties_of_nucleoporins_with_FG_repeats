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

1.1.2. Optional. Replacement of species IDs in fasta file with species names.
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

2.1. The MAIN function allow to take a ramdom subset of sequenses, perform two-step multiple alignment, rearrange the amyloidogenic score (Cumulative Score from ArchCandy) and calculate different statistics for each position of the alignment.
The MAIN function launch three subfunctions listed below for 10 times.

2.1.1. muscle_refine --- DESCRIPTION in the script filter_plus_alignment.R

2.1.2. aa_stat --- DESCRIPTION in the script ali_aa_statistics_plot_R.R

2.1.3. CS_subset --- DESCRIPTION in the script CS_subset_final.R

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
1) position_to_column (function, the output is shown in the console):
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
   
