options(scipen=100)

startTime <- Sys.time()

source(paste0(< path to R file with missing data information (callers that fail); data frame with caller/resol/norm columns>))  # load missing_data_DT
source(paste0(<path to R file with color information>)) # => imports the colors, also the label sizes for the plots

returnNA <- TRUE 

library(optparse)
library(foreach)
library(doMC)
library(dplyr)

registerDoMC(10)

my_callers <- all_callers
my_resolutions <- c("10","50", "100", "250", "1000")
my_norm <- c("ICE", "LGF")

# hard-coded settings used for our analyses
all_resolutions <- as.character(sort(as.numeric(as.character(all_resolutions))))
all_resolutions <- "10"
all_norm <- "ICE"

subLab <- "(ordered from left to right by increasing % of uniqueness)"

outFold <- file.path(<path to output folder>)
system(paste0("mkdir -p ", outFold))

plotType <- "pdf"
myHeight <- ifelse(plotType == "png", 480,7)
myWidth <- ifelse(plotType == "png", 600, 10)

chg_mar <- c(3,0,0,5)

cexLeg <- 0.9

chrFile <- <path to file containing chromosome size information (table without header, size information in [1,2] table cell)
chrSize <- read.delim(chrFile, header=F)[1,2]

allCallers_DT <- foreach(curr_caller = all_callers, .combine="rbind") %dopar% {
  cat(paste0("... start: ", curr_caller, "\n"))

  all_DT <- foreach(curr_norm = all_norm, .combine = 'rbind') %do% {
    # ensure all files in the same order
    dt_norm <- foreach(res=all_resolutions, .combine='rbind') %do% {
      if(any(as.character(missing_data_DT$caller) == as.character(curr_caller) &
             as.character(missing_data_DT$resol) == as.character(res) &
             as.character(missing_data_DT$norm) == as.character(curr_norm))){
        if(returnNA) {
          return(data.frame(chromo = NA, start = NA, end = NA, resol = res, norm = curr_norm))
        } else{
          return(NULL)
        }
      }
      myfile <- <retrieve partition for the current caller, normalization, resolution [BED format chromo/start/end without header]>
      tmp_dt <- read.delim(myfile, header=FALSE, stringsAsFactors = FALSE)
      colnames(tmp_dt) <- c("chromo", "start", "end")
      tmp_dt$resol <- rep(res, nrow(tmp_dt))
      tmp_dt$norm <- rep(curr_norm, nrow(tmp_dt))
      tmp_dt
    }
    dt_norm
  }
  if(is.null(all_DT))
    return(NULL)
  if(nrow(all_DT) < 1)
    return(NULL)
  if(nrow(na.omit(all_DT)) < 1)
    return(NULL)
  all_DT <- as.data.frame(all_DT)
  all_DT$start <- as.numeric(as.character(all_DT$start))
  all_DT$end <- as.numeric(as.character(all_DT$end))
  all_DT$size <- all_DT$end - all_DT$start + 1
  all_DT$resol <- factor(all_DT$resol, levels= all_resolutions)
  all_DT$caller <- curr_caller

  all_DT
}

# hard-coded: do it for bin size 10kb
# fig. 4c and 4d: 2* bin size tolerance radius
binSize <- 10000
tolRad <- 2*binSize
domains_DT$start_end <- paste0( domains_DT$chromo, "_", domains_DT$start, "_", domains_DT$end)

domainsMatch_DT <- foreach(i = seq_len(nrow(domains_DT)), .combine='rbind') %dopar% {
  curr_start <- domains_DT$start[i]
  curr_end <- domains_DT$end[i]
  start_start_match <- abs(domains_DT$start - curr_start) <= tolRad
  start_end_match <- abs(domains_DT$end - curr_start) <= tolRad
  end_end_match <- abs(domains_DT$end - curr_end) <= tolRad
  end_start_match <- abs(domains_DT$start - curr_end) <= tolRad
  nStrictConserv <- length(unique(domains_DT$caller[start_start_match & end_end_match]))
  nStartConserv <- length(unique(domains_DT$caller[start_start_match]))
  nEndConserv <- length(unique(domains_DT$caller[end_end_match]))
  nStartBDconserv <- length(unique(domains_DT$caller[start_start_match | start_end_match]))
  nEndBDconserv <- length(unique(domains_DT$caller[end_end_match | end_start_match]))
  data.frame(nStrictConserv=nStrictConserv,
             nStartConserv=nStartConserv,
             nEndConserv=nEndConserv,
             nStartBDconserv=nStartBDconserv,
             nEndBDconserv=nEndBDconserv,
             stringsAsFactors = F
  )
}

stopifnot(nrow(domainsMatch_DT) == nrow(domains_DT))
allDomains_match_DT <- cbind(domains_DT, domainsMatch_DT)

####################################################################
# get conservation by caller
nCallers <- length(unique(allDomains_match_DT$caller))

conservByCaller_DT <- do.call(rbind, by(allDomains_match_DT, allDomains_match_DT$caller, function(x) {
  tmpBoundaries <- c(x$nStartBDconserv, x$nEndBDconserv)
  tad_dt <- data.frame( caller = unique(x$caller),
            type ="domain",            
            nbrConserv = names(table(x$nStrictConserv)),
             nbrRegions =  as.numeric(table(x$nStrictConserv)), 
            stringsAsFactors = FALSE)            
  boundary_dt <- data.frame( caller = unique(x$caller),
                       type ="boundary",            
                       nbrConserv =  names(table(tmpBoundaries)),
                       nbrRegions =  as.numeric(table(tmpBoundaries)), 
                       stringsAsFactors = FALSE) 
  rbind(tad_dt, boundary_dt)  
})
)
rownames(conservByCaller_DT) <- NULL

rowOrder_init <- c("1", "2-5", "6-10", "11-15", ">15")

#################################################################### plot TAD conservation
byCaller_TAD_DT <- conservByCaller_DT[conservByCaller_DT$type == "domain",]
byCaller_TAD_DT$type <- NULL
# split the conservation for different # of callers
byCaller_TAD_DT$nbrConserv <- as.numeric(as.character(byCaller_TAD_DT$nbrConserv))
byCaller_TAD_DT$nbrConserv <- ifelse(byCaller_TAD_DT$nbrConserv == 1, rowOrder_init[1], 
                                     ifelse(byCaller_TAD_DT$nbrConserv <= 5, rowOrder_init[2], 
                                        ifelse(byCaller_TAD_DT$nbrConserv <= 10, rowOrder_init[3],
                                            ifelse(byCaller_TAD_DT$nbrConserv <= 15, rowOrder_init[4],
                                                   ifelse(byCaller_TAD_DT$nbrConserv >15, rowOrder_init[5],NA)))))
stopifnot(!is.na(byCaller_TAD_DT$nbrConserv))
rowOrder <- rowOrder_init
rowOrder <- rowOrder[rowOrder %in% byCaller_TAD_DT$nbrConserv]

stopifnot(is.numeric(byCaller_TAD_DT$nbrRegions))
byCaller_TAD_DT <- aggregate(nbrRegions ~., data=byCaller_TAD_DT, FUN=sum)

byCaller_TAD_DT <- reshape(byCaller_TAD_DT, timevar="caller", idvar="nbrConserv", direction="wide")
byCaller_TAD_DT[is.na(byCaller_TAD_DT)] <- 0
colnames(byCaller_TAD_DT) <- gsub("nbrRegions.", "", colnames(byCaller_TAD_DT))
rownames(byCaller_TAD_DT) <- byCaller_TAD_DT$nbrConserv
byCaller_TAD_DT$nbrConserv <- NULL
byCaller_TAD_DT <- byCaller_TAD_DT[rowOrder,]

tmpDT <- byCaller_TAD_DT[rowOrder[1],]/colSums(byCaller_TAD_DT, na.rm=T)
myOrder <- names(sort(tmpDT, decreasing=FALSE))
byCaller_TAD_DT <- byCaller_TAD_DT[, myOrder]

nCatego <- nrow(byCaller_TAD_DT)

if(nCatego <= 11) {
  myCols <- rev(brewer.pal(n = nCatego, name = "Spectral"))
} else {
  myCols <- colorRampPalette(c("blue", "red"))( nCatego )
}

totDomains <- colSums(byCaller_TAD_DT, na.rm=T)

# plot the ratio
byCaller_TAD_DT_ratio <- byCaller_TAD_DT
byCaller_TAD_DT_ratio <- apply(byCaller_TAD_DT_ratio, 2, function(x) x/sum(x) * 100)

outFile <- file.path(<path to outFile>)
pdf(outFile, width=myWidth, height=myHeight)
par(mar = par()$mar + chg_mar)
xpos <- barplot(as.matrix(byCaller_TAD_DT_ratio), las=2, col=myCols, ylim=c(0, 100+10),
                cex.axis = 0.7,
                axes = FALSE, ylab="% of TADs")
title("Domain conservation (ratio)")
mtext(side=3, subLab, font=3, cex =0.8)
axis(2, at = seq(0,100, 20), labels=seq(0,100, 20), las=1)
par(xpd=TRUE)
legend(x=max(xpos)+1, y=1,xjust=0, yjust=0, 
       legend=rownames(byCaller_TAD_DT_ratio), cex=cexLeg, fill=myCols, bty="n", title="conserv. in\n# callers")
text(x = xpos, y = 100, label = totDomains, cex=0.6, pos=3)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


#################################################################### plot boundary conservation
### => do the same for the boundaries
byCaller_BD_DT <- conservByCaller_DT[conservByCaller_DT$type == "boundary",]
byCaller_BD_DT$type <- NULL
byCaller_BD_DT$nbrConserv <- as.numeric(as.character(byCaller_BD_DT$nbrConserv))
byCaller_BD_DT$nbrConserv <- ifelse(byCaller_BD_DT$nbrConserv == 1, rowOrder_init[1], 
                                    ifelse(byCaller_BD_DT$nbrConserv <= 5, rowOrder_init[2], 
                                         ifelse(byCaller_BD_DT$nbrConserv <= 10, rowOrder_init[3],
                                            ifelse(byCaller_BD_DT$nbrConserv <= 15, rowOrder_init[4],
                                                   ifelse(byCaller_BD_DT$nbrConserv >15, rowOrder_init[5],NA)))))
stopifnot(!is.na(byCaller_BD_DT$nbrConserv))
rowOrder <- rowOrder_init
rowOrder <- rowOrder[rowOrder %in% byCaller_BD_DT$nbrConserv]

byCaller_BD_DT <- aggregate(nbrRegions ~., data=byCaller_BD_DT, FUN=sum)

byCaller_BD_DT <- reshape(byCaller_BD_DT, timevar="caller", idvar="nbrConserv", direction="wide")
byCaller_BD_DT[is.na(byCaller_BD_DT)] <- 0
colnames(byCaller_BD_DT) <- gsub("nbrRegions.", "", colnames(byCaller_BD_DT))
rownames(byCaller_BD_DT) <- byCaller_BD_DT$nbrConserv
byCaller_BD_DT$nbrConserv <- NULL
byCaller_BD_DT <- byCaller_BD_DT[rowOrder,]


tmpDT <- byCaller_BD_DT[rowOrder[1],]/colSums(byCaller_BD_DT, na.rm=T)
myOrder <- names(sort(tmpDT, decreasing=F))
byCaller_BD_DT <- byCaller_BD_DT[, myOrder]

nCatego <- nrow(byCaller_BD_DT)

if(nCatego <= 11) {
  myCols <- rev(brewer.pal(n = nCatego, name = "Spectral"))
} else {
  myCols <- colorRampPalette(c("blue", "red"))( nCatego )
}

totBound <- colSums(byCaller_BD_DT, na.rm=T)

byCaller_BD_DT_ratio <- byCaller_BD_DT
byCaller_BD_DT_ratio <- apply(byCaller_BD_DT_ratio, 2, function(x) x/sum(x) * 100)

outFile <- file.path(<path to outFile>)
pdf(outFile, width=myWidth, height=myHeight)
par(mar = par()$mar + chg_mar)
xpos <- barplot(as.matrix(byCaller_BD_DT_ratio), las=2, col=myCols, ylim=c(0, 100+10),
                cex.axis = 0.7,
                axes = FALSE, ylab="% of boundaries")
title("Boundary conservation (ratio)")
mtext(side=3, subLab, font=3, cex =0.8)
axis(2, at = seq(0,100, 20), labels=seq(0,100, 20), las=1)
par(xpd=TRUE)
legend(x=max(xpos)+1, y=1,xjust=0, yjust=0, legend=rownames(byCaller_BD_DT_ratio), cex=cexLeg, fill=myCols, bty="n", title="conserv. in\n# callers")
text(x = xpos, y = 100, label = totBound, cex=0.6, pos=3)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


#######
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))
