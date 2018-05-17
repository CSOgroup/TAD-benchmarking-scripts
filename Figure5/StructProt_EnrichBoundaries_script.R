#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
TadsFile = args[1] # Tab-separated file containing the list of TADs. Each line (no header) should represent a TAD, with genomic coordinates (chr, start, end)
chr = as.numeric(args[2]) # Chromosome (e.g. '6')
resolution_kb = as.numeric(args[3]) # Resolution in kb (e.g. '10')
OutFolder = args[4] # Folder where results should be saved
ChrSizes_file = args[5] # Tab-separated file containing the chromosome sizes. Each line should correspond to a chromosome (e.g. chr1   249250621)
ProteinPeaks_Folder = args[6] # Folder containing the three lists of peaks, one for each protein, in .bed format. Each list should be named as 'protein'_peaks.bed (e.g. CTCF_peaks.bed)

#### StructProt_EnrichBoundaries_script.R
# Script to assess the enrichment of CTCF, RAD21 and SMC3 ChIP-seq peaks at TAD boundaries for a given TAD partition.
# Input: see command line arguments.
# Output:
# - Enrich_StructProt_*.pdf: the plot of the average number of peaks every 5 kb around TAD boundaries for the three proteins
# - Enrich_StructProt_*.txt: the table used to generate the plot, with all the values for each protein
# - StructProteins_*.txt: a table containing, for each protein, the ratio of boundaries tagged, the fold change of the peak and its empirical p-value
#
# Refer to Zufferey & Tavernari et al. for details.
# 
# Contact daniele.tavernari@unil.ch for enquiries

####### START ######
start_time <- Sys.time()
####################

###### INPUT ######

OutFolder = paste0(OutFolder,"/")
color1 = "darkgoldenrod2" 
color2 = "royalblue4"
color3 = "orangered3"

res_step = 5
proteins = c("CTCF", "RAD21", "SMC3")

ChrSizes = read.table(file = ChrSizes_file, quote = '', stringsAsFactors=FALSE, sep = "\t", header = FALSE)
ChrSize = as.numeric(ChrSizes[ChrSizes[,1]==paste0('chr',as.character(chr)),2])

###################

############## ENRICHMENT AT BOUNDARY SITES ##############
Final_Ratios_df <- data.frame(TadsFile = TadsFile, resolution_kb = NA, protein = NA, domains_ratio = NA)
Final_EnrichPeaks_df <- as.data.frame(matrix(nrow = 0, ncol = (as.integer(1000/res_step)+3) ))
colnames(Final_EnrichPeaks_df) <- c("TadsFile", "resolution_kb", "protein", as.character(c(1:as.integer(1000/res_step))))

Tads_chr <- read.table(TadsFile, quote = '', stringsAsFactors=FALSE, sep = "\t", header = FALSE)
Allprots_enrichs <- data.frame(row.names = c(1:as.integer(1000/res_step)))

for (protein in proteins)
{
   Peaks <- read.table(paste0(ProteinPeaks_Folder, protein,'_peaks.bed'), quote = '', stringsAsFactors=FALSE, sep = "\t", header = FALSE)
   Peaks_chr <- Peaks[Peaks[,1]==paste0("chr",as.character(chr)),c(2,3)]
   Enrich_df = data.frame(matrix(vector(), 0, as.integer(1000/res_step)))

   Tad_edges <- unique(c(Tads_chr[,2], (Tads_chr[,3]+1)))
   for (tad_edge in Tad_edges)
   {
      if (((tad_edge-500*1000)<0) | ((tad_edge+500*1000)>ChrSize)) next
      Tad_surroundings <- seq(from = tad_edge-500*1000, to = tad_edge+500*1000, by = res_step*1000)
      Peaks_tad <- Peaks_chr[Peaks_chr[,2]>=(tad_edge-500*1000) & Peaks_chr[,1]<=(tad_edge+500*1000),]
      Vec_tad <- vector(mode = "numeric", length = as.integer(1000/res_step))
      for (pt in 1:dim(Peaks_tad)[1])
      {
         query <- as.numeric(Peaks_tad[pt,])
         query <- query[(query>Tad_surroundings[1]) & (query<Tad_surroundings[length(Tad_surroundings)]) ]
         a <- hist(query, breaks = Tad_surroundings, plot = FALSE)
         Vec_tad <- Vec_tad + as.numeric(a$counts>0)
      }
      Enrich_df <- rbind(Enrich_df, Vec_tad)
   }
   Ratio_tagged <- sum(as.numeric(rowSums(Enrich_df[,c((length(Vec_tad)/2-resolution_kb/res_step+1):(length(Vec_tad)/2+resolution_kb/res_step))])>0))/dim(Enrich_df)[1]
   Allprots_enrichs[,protein] <- colMeans(Enrich_df)
   Ratios_temp_df <- data.frame(TadsFile = TadsFile, resolution_kb = resolution_kb, protein = protein, domains_ratio = Ratio_tagged)
   Final_Ratios_df <- rbind(Final_Ratios_df, Ratios_temp_df)
   
}

pdf(paste0(OutFolder, "Enrich_StructProt_chr",chr,"_res",resolution_kb,"kb.pdf"))
ymax <- max(Allprots_enrichs) 
ymax = max(ymax, 0.14)
ymin <- min(Allprots_enrichs) 
plot_df = data.frame(row.names = c(1:as.integer(1000/res_step)), x = c(1:as.integer(1000/res_step)), CTCF = Allprots_enrichs[,"CTCF"], RAD21 = Allprots_enrichs[,"RAD21"], SMC3 = Allprots_enrichs[,"SMC3"] )
write.table(plot_df, file = paste0(OutFolder, "Enrich_StructProt_chr",chr,"_res",resolution_kb,"kb.txt"), row.names = F, col.names = T, quote = F )
plot(x=c(1:as.integer(1000/res_step)), y=Allprots_enrichs[,"CTCF"], type = "l", col = color1, main = paste0("chr",chr,", resolution ",resolution_kb,"kb"), xlab = "Distance from TAD boundary", ylab = paste0("Average # of peaks / ", as.character(res_step)," kb"), xaxt='n', cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
axis(1, at=seq(from = 0, to = as.integer(1000/res_step), by = as.integer(1000/res_step)/4), labels=c("-500 kb","-250 kb", "0", "+250 kb", "+500 kb"), cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
lines(x=c(1:as.integer(1000/res_step)), y=Allprots_enrichs[,"RAD21"], col = color2)
lines(x=c(1:as.integer(1000/res_step)), y=Allprots_enrichs[,"SMC3"], col = color3)
legend(x="topright", legend = proteins, col = c(color1, color2, color3), lty = 1, cex = 1.5)
dev.off()

Final_EnrichPeaks_temp_annotation_df <- data.frame(TadsFile = TadsFile, resolution_kb = resolution_kb, protein = proteins)
Final_EnrichPeaks_temp_values_df <- as.data.frame(t(Allprots_enrichs))
Final_EnrichPeaks_temp_df <- cbind(Final_EnrichPeaks_temp_annotation_df, Final_EnrichPeaks_temp_values_df)
Final_EnrichPeaks_df <- rbind(Final_EnrichPeaks_df, Final_EnrichPeaks_temp_df)

##########################################################

################# QUANTIFYING ENRICHMENT #################

mean_fc_over_bg <- function( profile, na.rm = TRUE, ratio_from_peak = 0.01, ratio_from_end = 0.6 ){
   if (na.rm) profile <- profile[!is.na(profile)]
   peak_val <- mean(profile[ (as.integer(length(profile)/2)-as.integer(length(profile)*ratio_from_peak)) : (as.integer(length(profile))/2+as.integer(length(profile)*ratio_from_peak)) ])
   bg <- mean(profile[c( 1:(as.integer(length(profile)/2)-as.integer(length(profile)*(1-ratio_from_end))), (as.integer(length(profile))/2+as.integer(length(profile)*(1-ratio_from_end))):length(profile) )])
   
   fc_over_bg <- peak_val/bg - 1
   return(fc_over_bg)
}

mean_pval_of_peak <- function( profile, plot_dist = FALSE, na.rm = TRUE, ratio_from_peak = 0.01, ratio_from_end = 0.6 ){
   if (na.rm) profile <- profile[!is.na(profile)]
   peak_val <- mean(profile[ (as.integer(length(profile)/2)-as.integer(length(profile)*ratio_from_peak)) : (as.integer(length(profile))/2+as.integer(length(profile)*ratio_from_peak)) ])
   bg <- profile[c( 1:(as.integer(length(profile)/2)-as.integer(length(profile)*(1-ratio_from_end))), (as.integer(length(profile))/2+as.integer(length(profile)*(1-ratio_from_end))):length(profile) )]

   if (plot_dist)
   {
      hist(bg, 30, xlim = c(min(bg), max(peak_val, max(bg) )),  main = "Distribution of background and value of peak")
      abline(v = peak_val, col = "red")
   }
   
   pval_of_peak <- pnorm( peak_val, mean = mean(bg), sd = sd(bg), lower.tail = FALSE )
   return(pval_of_peak)
}

Final_Ratios_df$fc_over_bg <- NA
Final_Ratios_df$pval_of_peak <- NA

for (protein in proteins)
{
   profile <- Final_EnrichPeaks_df[ ((Final_EnrichPeaks_df$TadsFile==TadsFile) & (Final_EnrichPeaks_df$resolution_kb==resolution_kb) & (Final_EnrichPeaks_df$protein==protein) ) , c(4:203)]

   fc_over_bg_val <- mean_fc_over_bg(profile)
   pval_of_peak_val <- mean_pval_of_peak(profile, plot_dist = FALSE)

   Final_Ratios_df$fc_over_bg[ Final_Ratios_df$TadsFile==TadsFile & Final_Ratios_df$resolution_kb==resolution_kb & Final_Ratios_df$protein==protein ] <- fc_over_bg_val
   Final_Ratios_df$pval_of_peak[ Final_Ratios_df$TadsFile==TadsFile & Final_Ratios_df$resolution_kb==resolution_kb & Final_Ratios_df$protein==protein ] <- pval_of_peak_val
}
Final_Ratios_df = Final_Ratios_df[!is.na(Final_Ratios_df$protein),]
write.table(Final_Ratios_df, file = paste0( OutFolder, "StructProteins_chr",chr,"_res",resolution_kb,"kb.txt" ), quote = FALSE, sep = "\t", row.names = FALSE)

##########################################################

####### END #######
end_time <- Sys.time()
time_taken <- end_time - start_time
print(time_taken)
###################