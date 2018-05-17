**StructProt_EnrichBoundaries_script.R** 

This script assesses the enrichment of CTCF, RAD21 and SMC3 ChIP-seq peaks at TAD boundaries for a given TAD partition.
See methods section in Zufferey & Tavernari et al. for details.

_Input_

* TadsFile: Tab-separated file containing the list of TADs. Each line (no header) should represent a TAD, with genomic coordinates (chr, start, end)
* chr: Chromosome (e.g. '6')
* resolution_kb: Resolution in kb (e.g. '10')
* OutFolder: Folder where results should be saved
* ChrSizes_file: Tab-separated file containing the chromosome sizes. Each line should correspond to a chromosome (e.g. chr1   249250621)
* ProteinPeaks_Folder: Folder containing the three lists of peaks, one for each protein, in .bed format. Each list should be named as 'protein'\_peaks.bed (e.g. CTCF_peaks.bed)

_Output_

* Enrich_StructProt\_\*.pdf: the plot of the average number of peaks every 5 kb around TAD boundaries for the three proteins
* Enrich_StructProt\_\*.txt: the table used to generate the plot, with all the values for each protein
* StructProteins\_\*.txt: a table containing, for each protein, the ratio of boundaries tagged, the fold change of the peak and its empirical p-value

_Requires_

R
