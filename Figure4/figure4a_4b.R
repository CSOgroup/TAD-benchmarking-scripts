options(scipen=100)

startTime <- Sys.time()

source(paste0(< path to R file with missing data information (callers that fail); data frame with caller/resol/norm columns>))  # load missing_data_DT
source(paste0(<path to R file with color information>)) # => imports the colors, also the label sizes for the plots

returnNA <- TRUE # if NA -> in the plot but nothing, if FALSE, not appear in the plot at all

library(foreach)
library(doMC)
library(tools)
library(lattice)

registerDoMC(10)

my_callers <- <callers to consider>
my_resolutions <- <which resolutions to plot [kb]>
my_norm <- <which normalization method [ICE or LGF]>

all_resolutions <- my_resolutions
all_norm <- my_norm

all_resolutions <- as.character(sort(as.numeric(as.character(all_resolutions))))

#****************** some hard-coded settings used for our analyses
all_resolutions <- "10"
all_norm <- "ICE"
binSize <- 10000
stopifnot(binSize == (as.numeric(all_resolutions) * 1000))
# which tolerance radius to iterate over
allTolRadNbrs <- c(0:5)

rowOrder_init <- c("1", "2-5", "6-10", "11-15", ">15")

# set TRUE or FALSE if the data are already prepared
# (if FALSE -> Rdata directly loaded for plotting)
buildTableDomains <- TRUE
buildTableMatch <- TRUE
buildAllTables <- TRUE

stopifnot(as.character(binSize/1000) == all_resolutions)
stopifnot(length(all_norm) == 1)

#******************************************************************
foldPrefix <- paste0(setDir, "/mnt/ed2/shared/TADcompare/pipeline/12_07_Cmap")
mypatt <- paste0(paste0(foldPrefix, "_", all_resolutions, "k_"), collapse="|")

subLab <- "(ordered from left to right by increasing % of uniqueness)"

outFold <- file.path(<path output folder>)
system(paste0("mkdir -p ", outFold))

plotType <- "svg"
myHeight <- ifelse(plotType == "png", 480,7)
myWidth <- ifelse(plotType == "png", 600, 10)

myWidthFacet <- 12
myHeightFacet <- 8

myHeightJoy <- 14
myWidthJoy <- 10

myHeightLattice <- 10
myWidthLattice <- 12

chrFile <- <path to file containing chromosome size information (table without header, size information in [1,2] table cell)
chrSize <- read.delim(chrFile, header=F)[1,2]

if(buildAllTables) {
  #******************************************* START ITERATING OVER THE TOLERANCE RADIUS TO BUILD THE TABLE
  for(tolRadNbr in allTolRadNbrs) {
    tolRad <- tolRadNbr * binSize
    cat(paste0("**START tolRadNbr = ", tolRadNbr, " => tolRad = ", tolRad, "\n"))
    ########################################################################################################################################  BUILD domains_DT
    if(buildTableDomains) {
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
                      myfile <- <path to file with partition for the current caller, normalization, resolution [BED format chromo/start/end without header]>
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
                outFile <- file.path(<path to output file>)
                allCallers_DT$tolRadNbr <- tolRadNbr
                save(allCallers_DT, file = outFile)
    } else {
              outFile <- file.path(<path to output file>)
    }
    ########################################################################################################################################  BUILD domainsMatch_DT
    domains_DT <- eval(parse(text = load(outFile)))
    assign(paste0("allDomains_DT_", tolRadNbr), domains_DT)
    stopifnot(as.character(domains_DT$resol) == all_resolutions)
    stopifnot(as.character(domains_DT$norm) == all_norm)
    domains_DT$start_end <- paste0( domains_DT$chromo, "_", domains_DT$start, "_", domains_DT$end)
        
    if(buildTableMatch) {
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
              outFile <- file.path(<path to output file>)
              allDomains_match_DT$tolRadNbr <- tolRadNbr
              save(allDomains_match_DT, file = outFile)
    } else{
       outFile <- file.path(<path to output file>)
    }
    assign(paste0("allDomains_match_DT_", tolRadNbr), allDomains_match_DT)
    
    ########################################################################################################################################  BUILD byCaller_TAD_DT_ratio
    allDomains_match_DT <- eval(parse(text = load(outFile)))
    nCallers <- length(unique(allDomains_match_DT$caller))
    conservByCaller_DT <- do.call(rbind, by(allDomains_match_DT, allDomains_match_DT$caller, function(x) {
      tmpBoundaries <- c(x$nStartBDconserv, x$nEndBDconserv)
      tad_dt <- data.frame( caller = unique(x$caller),
                type ="domain",            
                nbrConserv = names(table(x$nStrictConserv)),
                 nbrRegions =  as.numeric(table(x$nStrictConserv)), 
                stringsAsFactors = F)            
      boundary_dt <- data.frame( caller = unique(x$caller),
                           type ="boundary",            
                           nbrConserv =  names(table(tmpBoundaries)),
                           nbrRegions =  as.numeric(table(tmpBoundaries)), 
                           stringsAsFactors = F) 
      rbind(tad_dt, boundary_dt)  
    })
    )
    rownames(conservByCaller_DT) <- NULL
    
    byCaller_TAD_DT <- conservByCaller_DT[conservByCaller_DT$type == "domain",]
    byCaller_TAD_DT$type <- NULL
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
    myOrder <- names(sort(tmpDT, decreasing=F))
    byCaller_TAD_DT <- byCaller_TAD_DT[, myOrder]
    nCatego <- nrow(byCaller_TAD_DT)
    
    if(nCatego <= 11) {
      myCols <- rev(brewer.pal(n = nCatego, name = "Spectral"))
    } else {
      myCols <- colorRampPalette(c("blue", "red"))( nCatego )
    }
    
    totDomains <- colSums(byCaller_TAD_DT, na.rm=T)
    byCaller_TAD_DT_ratio <- byCaller_TAD_DT
    byCaller_TAD_DT_ratio <- apply(byCaller_TAD_DT_ratio, 2, function(x) x/sum(x) * 100)
    
    outFile <- file.path(<path to output file>)
    save(byCaller_TAD_DT_ratio, file = outFile)
    byCaller_TAD_DT_ratio <- eval(parse(text = load(outFile)))
    
    byCaller_TAD_DT_ratioDF <- as.data.frame(byCaller_TAD_DT_ratio)
    byCaller_TAD_DT_ratioDF$catego <- rownames(byCaller_TAD_DT_ratioDF)
    byCaller_TAD_DT_ratioDF$tolRadNbr <- tolRadNbr
    rownames(byCaller_TAD_DT_ratioDF) <- NULL
    assign(paste0("byCaller_TAD_DT_ratioDF_", tolRadNbr), byCaller_TAD_DT_ratioDF)
    
    ########################################################################################################################################  BUILD byCaller_BD_DT_ratio
    ####### boundary conservation
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
    # tmpDT[is.na(tmpDT)] <- 0
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
    
    outFile <- file.path(<path to output file>)
    save(byCaller_BD_DT_ratio, file = outFile)
    assign(paste0("byCaller_BD_DT_ratio_", tolRadNbr), byCaller_BD_DT_ratio)
    
    byCaller_BD_DT_ratio <- eval(parse(text = load(outFile)))
    byCaller_BD_DT_ratioDF <- as.data.frame(byCaller_BD_DT_ratio)
    byCaller_BD_DT_ratioDF$catego <- rownames(byCaller_BD_DT_ratioDF)
    byCaller_BD_DT_ratioDF$tolRadNbr <- tolRadNbr
    rownames(byCaller_BD_DT_ratioDF) <- NULL
    assign(paste0("byCaller_BD_DT_ratioDF_", tolRadNbr), byCaller_BD_DT_ratioDF)
    
    ########################################################################################################################################  BUILD strictDomainConserv_DT 
    allDomains_match_DT$start_end <- paste0(allDomains_match_DT$chromo, "_", allDomains_match_DT$start, "_", allDomains_match_DT$end)
    
    # for this part, do not take duplicates !
    allDomainsUnique_match_DT <- allDomains_match_DT[!duplicated(allDomains_match_DT$start_end),]
    rm(allDomains_match_DT)
    # caller does not make sense anymore
    allDomainsUnique_match_DT$caller <- NULL
    
    # plot how many conserved ---> for the domains
    strictDomainConserv_DT <- data.frame(nbrConserv = names(table(allDomainsUnique_match_DT$nStrictConserv)),
                                         nbrDomains =  as.numeric(table(allDomainsUnique_match_DT$nStrictConserv)),
                                         stringsAsFactors = F)
    strictDomainConserv_DT$nbrDomains_cumsum <- cumsum(strictDomainConserv_DT$nbrDomains)
    
    outFile <- file.path(<path to output file>)
    strictDomainConserv_DT$tolRadNbr <- tolRadNbr
    save(strictDomainConserv_DT, file = outFile)
    assign(paste0("strictDomainConserv_DT_", tolRadNbr), strictDomainConserv_DT)
    strictDomainConserv_DT <- eval(parse(text = load(outFile)))
    
    ########################################################################################################################################  SAVE boundConserv_DT 
    allBoundaries <- c(allDomainsUnique_match_DT$nStartBDconserv, allDomainsUnique_match_DT$nEndBDconserv)
    
    boundConserv_DT <- data.frame(nbrConserv = names(table(allBoundaries)),
                                  nbrBoundaries =  as.numeric(table(allBoundaries)),
                                  stringsAsFactors = F)
    
    boundConserv_DT$nbrBoundaries_cumsum <- cumsum(boundConserv_DT$nbrBoundaries)
    
    nTot <- sum(boundConserv_DT$nbrBoundaries)
    
    outFile <- file.path(outFold, paste0("boundConserv_DT_", all_resolutions, "_", all_norm, "_tolRadNbr", tolRadNbr, ".Rdata"))
    boundConserv_DT$tolRadNbr <- tolRadNbr
    save(boundConserv_DT, file = outFile)
    
    assign(paste0("boundConserv_DT_", tolRadNbr), boundConserv_DT)
    
    boundConserv_DT <- eval(parse(text = load(outFile)))
  } # end iterating over the tolRadNbrs  
    
  #************************************************************************************************************************* MERGE THE DATA FRAMES FOR ALL TOLERANCE RADIUS
  # get the name of the data frames
  all_TADratio_DT_names <- paste0("byCaller_TAD_DT_ratioDF_", allTolRadNbrs)
  # create a list of data frames
  all_TADratio_DT_list <- lapply(all_TADratio_DT_names, get)
  # combine into a single dataframe
  all_TADratio_DT <- do.call(rbind, all_TADratio_DT_list) 
  
  outFile <- file.path(outFold, paste0("all_TADratio_DT.Rdata"))
  save(all_TADratio_DT, file = outFile)
  
  all_BDratio_DT_names <- paste0("byCaller_BD_DT_ratioDF_", allTolRadNbrs)
  all_BDratio_DT_list <- lapply(all_BDratio_DT_names, get)
  all_BDratio_DT <- do.call(rbind, all_BDratio_DT_list) 
  
  outFile <- file.path(outFold, paste0("all_BDratio_DT.Rdata"))
  save(all_BDratio_DT, file = outFile)
  
  all_domainConserv_DT_names <- paste0("strictDomainConserv_DT_", allTolRadNbrs)
  all_domainConserv_DT_list <- lapply(all_domainConserv_DT_names, get)
  all_domainConserv_DT <- do.call(rbind, all_domainConserv_DT_list)
  
  outFile <- file.path(<path to outFile>)
  save(all_domainConserv_DT, file = outFile)
  
  all_BDconserv_DT_names <- paste0("boundConserv_DT_", allTolRadNbrs)
  all_BDconserv_DT_list <- lapply(all_BDconserv_DT_names, get)
  all_BDconserv_DT <- do.call(rbind, all_BDconserv_DT_list) 
  
  outFile <- file.path(<path to output file>)
  save(all_BDconserv_DT, file = outFile)
} # end - if buildAllTables

#*************************************************************************************************************************
#************************************************************************************************************************* PLOT ACROSS TOLERANCE RADIUS
#*************************************************************************************************************************
all_TADratio_DT <- eval(parse(text = load(file.path(<path to input file>))))
all_BDratio_DT <- eval(parse(text = load(file.path(<path to input file>))))
all_domainConserv_DT <- eval(parse(text = load(file.path(<path to input file>))))
all_BDconserv_DT <- eval(parse(text = load(file.path(<path to input file>))))

all_domainConserv_DT$nbrConserv <- as.numeric(as.character(all_domainConserv_DT$nbrConserv))
all_domainConserv_DT$nbrDomains <- as.numeric(as.character(all_domainConserv_DT$nbrDomains))
all_domainConserv_DT$nbrDomains_cumsum <- as.numeric(as.character(all_domainConserv_DT$nbrDomains_cumsum))
all_domainConserv_DT$tolRadNbr <- as.numeric(as.character(all_domainConserv_DT$tolRadNbr))

all_BDconserv_DT$nbrConserv <- as.numeric(as.character(all_BDconserv_DT$nbrConserv))
all_BDconserv_DT$nbrBoundaries <- as.numeric(as.character(all_BDconserv_DT$nbrBoundaries))
all_BDconserv_DT$nbrBoundaries_cumsum <- as.numeric(as.character(all_BDconserv_DT$nbrBoundaries_cumsum))
all_BDconserv_DT$tolRadNbr <- as.numeric(as.character(all_BDconserv_DT$tolRadNbr))

allToPlots <- c("boundaries", "domains") 

for(toPlot in allToPlots) {
  if(toPlot == "boundaries") {
    conserv_DT <-  all_BDconserv_DT
    ratio_DT <- all_BDratio_DT
  } else if (toPlot == "domains") {
    conserv_DT <- all_domainConserv_DT
    ratio_DT <- all_TADratio_DT
  } else {
    stop("...error...\n")
  }
  
  mytit <- paste0(toTitleCase(toPlot), ": # callers conservation across tolerance radius")
  subTit <- paste0(all_norm, "; bin size = ", all_resolutions, " kb")
  
  myColours <- brewer.pal(6,"Blues")
  ## Create your own list with
  my.settings <- list(
    superpose.polygon=list(col=myColours[2:5], border="transparent"),
    strip.background=list(col=myColours[6]),
    strip.border=list(col="black")
  )
  tolRadLev <- as.character(sort(unique(as.numeric(as.character(conserv_DT$tolRadNbr)))))
  conserv_DT$tolRadNbr <- factor(conserv_DT$tolRadNbr, levels = tolRadLev )
  conserv_DT$tolRadNbr_txt <- paste0("#*binSize = ", as.character(conserv_DT$tolRadNbr))
  conserv_DT$tolRadNbr_txt <- factor(conserv_DT$tolRadNbr_txt, levels = paste0("#*binSize = ", tolRadLev ))
  
  outFile <- file.path(outFold,paste0(toPlot, "_nbrConservHist_acrossTolRad_", all_norm, "_", all_resolutions, "_lattice.", plotType))
  p <- barchart(
               as.formula(paste0("nbr", toTitleCase(toPlot), " ~ nbrConserv | tolRadNbr_txt")), 
                data=conserv_DT,
                horiz=F,
                main=mytit,
                # scales=list(alternating=1, y = "free"),              
                scales=list(alternating=1),  
                xlab = "# callers", 
                ylab = paste0("# ", toPlot),
                panel=function(x,y,...){
                  panel.grid(h=-1, v=0); 
                  panel.barchart(x,y,...)
                },
                par.settings = my.settings,
                par.strip.text=list(col="white", font=2),
                origin = 0,
                col = myColours[2],
				layout = c(length(allTolRadNbrs), 1))
  
  # trellis.device(device = plotType, filename = outFile)
  do.call(plotType, list(file = outFile, height = myHeight*0.8, width = myWidth*3))
  print(p)
  foo <- dev.off()
} # end-for boundaries or domains



#######
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))
