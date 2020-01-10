#' Calculate Measure of Concordance (MoC); adapted from formula given in Pfitzner et al. 2009
#'
#' @param file1 Path to file for 1st partition (3-column BED format chromo-start-end; no header)
#' @param file2 Path to file for 2nd partition (3-column BED format chromo-start-end; no header)
#' @param chrSize Chromo size [bp]
#' @param nCpu Nbr of cpu available [use of foreach]
#' @param fillInter Consider inter-TAD gaps as TADs ? [in our analyses: TRUE]
#' @param meanWithInter Compute MoC for the inter-TADs separately and take the mean of TADs and interTADs MoCs [in our analyses: FALSE]
#' @param correctClust Add penalty term for the number of domains ? [in our analyses: FALSE]
#' @param binSize Bin size [bp]
#' @param noMinusOne Do not substract 1 to the MoC ? [in our analyses: FALSE]
#' @param fillInterGapZero Set score to 0 if comparing 1 TAD with interTAD ? [in our analyses: FALSE]

# requirement: doMC, foreach

get_MoC <- function(file1, file2, chrSize, nCpu = 1, fillInter=TRUE, meanWithInter = FALSE, correctClust = FALSE, binSize = NA, noMinusOne = FALSE, fillInterGapZero = FALSE) {

  library(foreach)
  library(doMC)

  if(correctClust)
    if(is.na(binSize))
      stop("should provide binSize for correctClust!\n")

  if(fillInterGapZero)
    fillInter <- TRUE
  
  registerDoMC(cores=nCpu)
  
  if(fillInter & meanWithInter)
    stop("not meaningful")
  
  if (file.info(as.character(file1))$size == 0) {
    set1DT <- data.frame(chromo=character(0), start = numeric(0), end =numeric(0))
  } else {
    set1DT <- read.delim(file1, header = FALSE, stringsAsFactors = FALSE)
    colnames(set1DT) <- c("chromo", "start", "end")
  }
  if (file.info(as.character(file2))$size == 0) {
    set2DT <- data.frame(chromo=character(0), start = numeric(0), end =numeric(0))
  } else {
    set2DT <- read.delim(file2, header=FALSE, stringsAsFactors = FALSE)
    colnames(set2DT) <- c("chromo", "start", "end")
  }
  if(fillInterGapZero){
    set1DT_nofilled <- set1DT
    set2DT_nofilled <- set2DT
  }
  
  if(fillInter) {
    set1DT <- fill_part(set1DT, chrSize)
    set2DT <- fill_part(set2DT, chrSize)
  }
  
  if(fillInterGapZero){
    if(nrow(set1DT) > 0)
    set1DT$partType <- unlist(sapply(1:nrow(set1DT), function(x)
                ifelse(any(set1DT_nofilled$start == set1DT$start[x] & set1DT_nofilled$end == set1DT$end[x]), "domain", "gap")))
    if(nrow(set2DT) > 0)
    set2DT$partType <- unlist(sapply(1:nrow(set2DT), function(x)
      ifelse(any(set2DT_nofilled$start == set2DT$start[x] & set2DT_nofilled$end == set2DT$end[x]), "domain", "gap")))
  }
  
  if(nrow(set1DT) > 0)
    rownames(set1DT) <- paste0("P", 1:nrow(set1DT))
  if(nrow(set2DT) > 0)
    rownames(set2DT) <- paste0("Q", 1:nrow(set2DT))
  
  if(nrow(set1DT) == 0 & nrow(set2DT) > 0){
    MoC_score <- 0
  } else if(nrow(set1DT) > 0 & nrow(set2DT) == 0){
    MoC_score <- 0
  } else if( (nrow(set1DT) == nrow(set2DT))  & ( all(set1DT$start == set2DT$start) & all(set1DT$end == set2DT$end) ) ) {
      MoC_score <- 1
  } else if( (nrow(set1DT) == nrow(set2DT)) & (nrow(set1DT) == 1) ) {
      MoC_score <- 1
  } else {
    all_fragmentDT <- foreach(i = 1:nrow(set1DT), .combine='rbind') %dopar% {
      ref_start <- set1DT$start[i]
      ref_end <- set1DT$end[i]
      # if exact match
      all_matches <- which((set2DT$start == ref_start & set2DT$end == ref_end)  |
                             # nested
                             (set2DT$start >= ref_start & set2DT$end <= ref_end) |
                             # overlap left
                             (set2DT$start <= ref_start & set2DT$end >= ref_start) |
                             # overlap right
                             (set2DT$start <= ref_end & set2DT$end >= ref_end))
      all_matches_2 <- which(set2DT$end >= ref_start & set2DT$start <= ref_end)
      stopifnot(all(all_matches == all_matches_2))
      
      fragmentDT <- foreach(i_match = all_matches, .combine='rbind') %dopar% {
        # "fragment" size
        tmp_range <- c(set2DT$start[i_match] :set2DT$end[i_match])
        tmp_range <- tmp_range[tmp_range >= ref_start & tmp_range <= ref_end]
        frag_overlap <- tmp_range[length(tmp_range)] - tmp_range[1] + 1
        c(rownames(set1DT)[i], rownames(set2DT)[i_match], frag_overlap)
      }
      fragmentDT
    }
  
    # there is no intersect -> 0
    if(is.null(all_fragmentDT)) {
      MoC_score <- 0
    } else {
      # through matrix otherwise drop if nrow=1
      all_fragmentDT <- as.data.frame(matrix(all_fragmentDT, ncol=3), stringsAsFactors=F)
      rownames(all_fragmentDT) <- NULL
      colnames(all_fragmentDT) <- c("set1", "set2", "intersect_size")
      all_fragmentDT$intersect_size <- as.numeric(as.character(all_fragmentDT$intersect_size))
      
      MoC_score <- foreach(i = 1:nrow(all_fragmentDT), .combine = 'sum') %dopar% {
        d1 <- all_fragmentDT$set1[i]
        d2 <- all_fragmentDT$set2[i]
        # get the size of the domain from P  
        Pi <- set1DT[d1, "end"]  - set1DT[d1, "start"] + 1
        # get the size of the domain from Q 
        Qi <- set2DT[d2, "end"]  - set2DT[d2, "start"] + 1
        Fij <- all_fragmentDT$intersect_size[i]
        moc <- ((Fij*Fij)/(Pi * Qi))
        if(fillInterGapZero) { 
           if(set1DT[all_fragmentDT$set1[i], "partType" ] != set2DT[all_fragmentDT$set2[i], "partType" ])
              moc <- 0          
        }
        moc
      }
      # if wanted, correct for the number of clusters
      if(correctClust) {
        # penalize number of clusters / number tot of bins
        penaltyTerm <- nrow(all_fragmentDT) / (chrSize/binSize)
        MoC_score <- MoC_score + (1 - penaltyTerm) 
      }
    }
    
	if(! noMinusOne)
      MoC_score <- MoC_score - 1

    MoC_score <- MoC_score/(sqrt(nrow(set1DT) * nrow(set2DT)) - 1) 
  }
  if(meanWithInter) {
    bd_only_set1DT <- bd_only(set1DT, chrSize = chrSize)
    file1_BD_only <- sub("_final_domains.txt", "_final_domains_BD_only.txt", file1)
    write.table(bd_only_set1DT, file = file1_BD_only, col.names = F, row.names = F, quote=F, sep="\t")
    
    bd_only_set2DT <- bd_only(set2DT, chrSize = chrSize)
    file2_BD_only <- sub("_final_domains.txt", "_final_domains_BD_only.txt", file2)
    write.table(bd_only_set2DT, file = file2_BD_only, col.names = F, row.names = F, quote=F, sep="\t")
    
      MoC_score_inter <- get_MoC(file1_BD_only, file2_BD_only,chrSize, fillInter=fillInter, meanWithInter = FALSE )
    MoC_score <- (MoC_score + MoC_score_inter)/2
  }
  MoC_score
}


################
fill_part <- function(fillDT, chrSize) {
  stopifnot(all(colnames(fillDT) == c("chromo", "start", "end")))
  
  if(nrow(fillDT) == 0) {
    fillDT <- data.frame(chromo = "chr6",
                                  start = 1,
                                  end = chrSize)
    return(fillDT)
  }
  
  # add the intra-TAD also as domains
  fillDT_BD <- foreach(x = 1:nrow(fillDT), .combine='rbind') %dopar% {
    if(x == 1) {
      if(fillDT$start[1] == 1) {
        tmpDT <- data.frame(chromo = fillDT$chromo[1],
                            start = fillDT$start[1],
                            end = fillDT$end[1])
      } else{
        tmpDT <- data.frame(chromo = c(fillDT$chromo[1],fillDT$chromo[1]),
                            start = c(1, fillDT$start[1]),
                            end = c(fillDT$start[1]-1,fillDT$end[1] ))
      }
    } else{
      if(fillDT$start[x] > fillDT$end[x-1] + 1) {
        tmpDT <- data.frame(chromo = c(fillDT$chromo[x], fillDT$chromo[x]),
                            start = c(fillDT$end[x-1]+1, fillDT$start[x]),
                            end = c(fillDT$start[x] -1, fillDT$end[x]))
      } else{
        tmpDT <- data.frame(chromo = fillDT$chromo[x],
                            start = fillDT$start[x],
                            end = fillDT$end[x])
      }
    }
  }
  fillDT_BD <- as.data.frame(fillDT_BD)
  fillDT_BD$start <- as.numeric(as.character(fillDT_BD$start))
  fillDT_BD$end <- as.numeric(as.character(fillDT_BD$end))
  lastrow <- nrow(fillDT_BD)
  if(fillDT_BD$end[lastrow] < chrSize ){
    endDT <- data.frame(chromo = fillDT_BD$chromo[lastrow],
                        start = fillDT_BD$end[lastrow] + 1,
                        end = chrSize)
    fillDT <- rbind(fillDT_BD, endDT)
  } else{
    fillDT <- fillDT_BD
  }
  
  stopifnot(all( fillDT$end > fillDT$start  ))
  stopifnot(!any(duplicated(fillDT$start)))
  stopifnot(!any(duplicated(fillDT$end)))
  return(fillDT)
}


