#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
TadsFile = args[1]
share = as.numeric(args[2])
nshuf = as.numeric(args[3])
fdr_thresh = as.numeric(args[4])
chr = as.numeric(args[5])
OutFolder = args[6]
h27_binned = args[7]
h36_binned = args[8]

binsOverlapping = TRUE

Results_df = data.frame(row.names = TadsFile, nTads = NA, binsize = NA, shareTadsSignif_lr_BH = NA )

Tads <- read.table(TadsFile, quote = '', stringsAsFactors=FALSE, sep = "\t", header = FALSE)
colnames(Tads) = c("chr", "start", "end")

# computing binsize
Tads[,"size"] = Tads[,"end"]-Tads[,"start"]
binsize = round(mean(Tads$size)*share)
Tads$size = NULL

h27_binned = read.table(h27_binned, quote = '', stringsAsFactors=FALSE, sep = "\t", header = FALSE)
h36_binned = read.table(h36_binned, quote = '', stringsAsFactors=FALSE, sep = "\t", header = FALSE)

colnames(h27_binned) = c("chr", "start", "end", "reads_h27")
colnames(h36_binned) = c("chr", "start", "end", "reads_h36")
h_binned = merge(h27_binned, h36_binned, by = c("chr", "start", "end") ) # "bin"     "chr.x"   "start.x" "end.x"   "reads.x" "chr.y"   "start.y" "end.y"   "reads.y"
h_binned = h_binned[((h_binned$reads_h27>0) & (h_binned$reads_h36>0)),]

x27 = h_binned$reads_h27
x36 = h_binned$reads_h36
sign = as.numeric(x27>x36)
sign[sign==0] = -1
h_binned$lr = log10(h_binned$reads_h27/h_binned$reads_h36)

Tads$avg_lr = 0
Results_df[TadsFile,"nTads"] = nrow(Tads)
for (tad in rownames(Tads))
{
   # with these conditions: also bins overlapping the boundaries
   if (binsOverlapping)
   {
      condition = (Tads[tad,"start" ]<=h_binned[,"start"] & Tads[tad,"end" ]>=h_binned[,"end"]) | (Tads[tad,"start" ]<h_binned[,"start"] & Tads[tad,"end" ]>h_binned[,"start"]) | (Tads[tad,"start" ]<h_binned[,"end"] & Tads[tad,"end" ]>h_binned[,"end"]) | (Tads[tad,"start" ]>=h_binned[,"start"] & Tads[tad,"end" ]<=h_binned[,"end"])
   }
   h_binned_tad = h_binned[ condition ,]
   if (nrow(h_binned_tad)==0)
   {
      Tads[tad, "avg_lr"] = 0
      next
   }
   h_binned_tad = h_binned_tad[order(h_binned_tad$start),]
   rownames(h_binned_tad) = 1:nrow(h_binned_tad)
   if (h_binned_tad["1","start"] < Tads[tad,"start" ]) {h_binned_tad[1,"start"] = Tads[tad,"start" ]}
   if (h_binned_tad[as.character(nrow(h_binned_tad)),"end"] > Tads[tad, "end" ]) {h_binned_tad[as.character(nrow(h_binned_tad)),"end"] = Tads[tad, "end" ]}
   Tads[tad, "avg_lr"] = sum(h_binned_tad[,"lr"]*( h_binned_tad[,"end"]-h_binned_tad[,"start"] ))/( h_binned_tad[nrow(h_binned_tad),"end"]-h_binned_tad[1,"start"] )
}

# Shuffle bins
Tads_ws = Tads
diff_each_lr = vector(mode = "numeric", length = nshuf)
for (s in 1:nshuf)
{
   h_binned_shuffled = h_binned
   shuf_pos = sample(rownames(h_binned_shuffled), size = length(rownames(h_binned_shuffled)))
   h_binned_shuffled[,"lr"] = h_binned_shuffled[shuf_pos,"lr"]
   
   Tads_ws[,paste0("avg_lr_shuf",s)] = 0
   for (tad in rownames(Tads))
   {
      # with these conditions: also bins overlapping the boundaries
      if (binsOverlapping)
      {
         condition = (Tads[tad,"start" ]<=h_binned_shuffled[,"start"] & Tads[tad,"end" ]>=h_binned_shuffled[,"end"]) | (Tads[tad,"start" ]<h_binned_shuffled[,"start"] & Tads[tad,"end" ]>h_binned_shuffled[,"start"]) | (Tads[tad,"start" ]<h_binned_shuffled[,"end"] & Tads[tad,"end" ]>h_binned_shuffled[,"end"]) | (Tads[tad,"start" ]>=h_binned_shuffled[,"start"] & Tads[tad,"end" ]<=h_binned_shuffled[,"end"])
      }
      h_binned_tad = h_binned_shuffled[ condition ,]
      if (nrow(h_binned_tad)==0)
      {
         Tads_ws[tad, paste0("avg_lr_shuf",s)] = 0
         next
      }
      h_binned_tad = h_binned_tad[order(h_binned_tad$start),]
      rownames(h_binned_tad) = 1:nrow(h_binned_tad)
      if (h_binned_tad[1,"start"] < Tads[tad,"start" ]) {h_binned_tad[1,"start"] = Tads[tad,"start" ]}
      if (h_binned_tad[nrow(h_binned_tad),"end"] > Tads[tad, "end" ]) {h_binned_tad[nrow(h_binned_tad),"end"] = Tads[tad, "end" ]}
      Tads_ws[tad, paste0("avg_lr_shuf",s)] = sum(h_binned_tad[,"lr"]*( h_binned_tad[,"end"]-h_binned_tad[,"start"] ))/( h_binned_tad[nrow(h_binned_tad),"end"]-h_binned_tad[1,"start"] )
      
   }
   
   s = s+1
   
}

Tads_ws_lr = Tads_ws[,c("start", "end", colnames(Tads_ws)[grepl("*lr*",colnames(Tads_ws))])]


# BH approach, lr
distr_shuf_lr = as.numeric(unlist(Tads_ws_lr[,paste0("avg_lr_shuf",c(1:nshuf))]))
median_shuf_lr = median(distr_shuf_lr)
Tads_ws_lr$pval = NA
Tads_ws_lr$qval_BH = NA
Tads_ws_lr$qval_Bonferroni = NA
for (tad in rownames(Tads_ws_lr))
{
   query = Tads_ws_lr[tad, "avg_lr"]
   if (query>=median_shuf_lr){ Tads_ws_lr[tad, "pval"] = sum(distr_shuf_lr>=query)/length(distr_shuf_lr) }
   if (query<median_shuf_lr){ Tads_ws_lr[tad, "pval"] = sum(distr_shuf_lr<=query)/length(distr_shuf_lr) }
}
Tads_ws_lr$qval_BH = p.adjust(Tads_ws_lr$pval, method = "BH")
Tads_ws_lr$qval_Bonferroni = p.adjust(Tads_ws_lr$pval, method = "bonferroni")


distr_real_lr = Tads_ws$avg_lr
distr_shuf_lr = as.numeric(unlist(Tads_ws[,paste0("avg_lr_shuf",c(1:nshuf))]))
pdf(paste0(OutFolder, "Density_lr_binsize",binsize,"_Real_vs_Shuf.pdf"))
d_shuf <- density(distr_shuf_lr, na.rm = T)
d_real <- density(distr_real_lr, na.rm = T)
plot(d_shuf, main = "Density plot of lr", xlab = "lr", ylab = "Density", col = "black")
lines(d_real, col = "red")
legend(x="topright", legend = c("With random permutations", "Real"), col = c("black","red"), lty = 1)
dev.off()

Results_df[TadsFile,"binsize"] = binsize
Results_df[TadsFile,"nTads"] = nrow(Tads)
Results_df[TadsFile,"shareTadsSignif_lr_BH"] = sum(Tads_ws_lr$qval_BH<fdr_thresh)/nrow(Tads)

write.table(Results_df, paste0(OutFolder, "Results_df_shuf",nshuf,"_share",share,"_fdrthresh",fdr_thresh,".txt"), sep = "\t", row.names = F, col.names = T )
