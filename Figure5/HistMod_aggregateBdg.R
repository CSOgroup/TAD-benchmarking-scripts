#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
intersection_file = args[1]
fileOut = args[2]
binsize = as.numeric(args[3])


is = read.table(file = intersection_file)
is = is[,c("V1", "V2", "V3", "V8","V9")]
colnames(is) = c("chr", "start", "end", "fc","bp")
is[,"fc_bp"] = as.numeric(as.character(is$fc))*as.numeric(as.character(is$bp))
is[is.na(is$fc_bp),'fc_bp'] = 0.0
is$fc = NULL
is$bp = NULL

is_aggr_df = aggregate(fc_bp~chr+start+end, data = is, FUN = sum)
is_aggr_df$fc_bp = as.numeric(is_aggr_df$fc_bp)/binsize

write.table(is_aggr_df, file = fileOut, row.names = F, col.names = F, sep="\t", quote = F)
