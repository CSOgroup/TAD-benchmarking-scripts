#!/bin/bash

#### HistMod_script.sh
# Script to assess the share of TADs showing a significant H3K27me3/H3K36me3 log10-ratio for a given TAD partition.
# Input: see following arguments.
# Output:
# - Results_df_*.txt: the number of TADs, the binsize used for permutation and the share of TADs showing a significant H3K27me3/H3K36me3 log10-ratio
# - Density_lr_*.pdf: The density plot of H3K27me3/H3K36me3 log10-ratio in real and permuted scenarios 
#
# Refer to Zufferey & Tavernari et al. for details.
# 
# Contact daniele.tavernari@unil.ch for enquiries

###### INPUT ######
TadsFile="TopDom_final_domains.txt" # Tab-separated file containing the list of TADs. Each line (no header) should represent a TAD, with genomic coordinates (chr, start, end)
h27_fc="HistoneSignals/ENCFF594HSG_chr6.bedGraph" # BedGraph file containing the fold change vs control of ChIP-seq tracks for H3K27me3 mark
h36_fc="HistoneSignals/ENCFF662QFK_chr6.bedGraph" # BedGraph file containing the fold change vs control of ChIP-seq tracks for H3K36me3 mark
chr=6 # Chromosome
OutFolder="Outdir2/" # Folder where results should be saved
share=0.1 # Share of average TAD size to use as bin size for permutation test
nshuf=10 # Number of shufflings for permutation test
fdr_thresh=0.1 # FDR threshold
###################

binsizes_file=${OutFolder}"/binsizedf_chr"${chr}"_share"${share}".txt"
./HistMod_computeBinsizes_script.R ${TadsFile} ${binsizes_file} ${share}

while read binsize; do
  	binnedFile=${OutFolder}"/chr"${chr}_"binsize"${binsize}".bed"
	./HistMod_binChr.py -o ${binnedFile} -c ${chr} -b ${binsize}

	# H3K27me3
	filein_fc=${h27_fc}
	filein_fc_binned_intersections=${OutFolder}"h27_fc_intersections.bedGraph"
	h27_binned=${OutFolder}"h27_fc_binsize"${binsize}".bedGraph"
	bedtools intersect -wao -a ${binnedFile} -b ${filein_fc} > ${filein_fc_binned_intersections}
	./HistMod_aggregateBdg.R ${filein_fc_binned_intersections} ${h27_binned} ${binsize}

	# H3K36me3
	filein_fc=${h36_fc}
	filein_fc_binned_intersections=${OutFolder}"h36_fc_intersections.bedGraph"
	h36_binned=${OutFolder}"h36_fc_binsize"${binsize}".bedGraph"
	bedtools intersect -wao -a ${binnedFile} -b ${filein_fc} > ${filein_fc_binned_intersections}
	./HistMod_aggregateBdg.R ${filein_fc_binned_intersections} ${h36_binned} ${binsize}
done <${binsizes_file}

./HistMod_FDRclassic_script.R ${TadsFile} ${share} ${nshuf} ${fdr_thresh} ${chr} ${OutFolder} ${h27_binned} ${h36_binned}
