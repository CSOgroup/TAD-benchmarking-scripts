#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
TadsFile = args[1]
outfile = args[2]
share = as.numeric(args[3])

binsize_df = data.frame(row.names = TadsFile, binsize = NA)

Tads <- read.table(TadsFile, quote = '', stringsAsFactors=FALSE, sep = "\t", header = FALSE)
colnames(Tads) = c("chr", "start", "end")
Tads[,"size"] = Tads[,"end"]-Tads[,"start"]
binsize_df[TadsFile,"binsize"] = round(mean(Tads$size)*share)

write.table(binsize_df, outfile, sep = "\t", row.names = F, col.names = F )
##