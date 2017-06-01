### Jinliang Yang
### April 24th, 2015


getmrkres <- function(dir="slurm-scripts/wholeset/"){
  
  files <- list.files(path = dir, pattern="\\.mrkRes")
  ## file line of the shell file:
  message(sprintf("[ %s ] mrkRes files detected!", length(files)))
  
  resout <- data.frame()
  for(i in 1:length(files)){
    myfile <- paste(dir, files[i], sep="")
    tab5rows <- read.table(myfile, header=TRUE, nrows=5)
    classes <- sapply(tab5rows, class)
    res <- read.table(myfile, header=TRUE, colClasses=classes)
    message(sprintf("###>>> input [ %s ] rows for [ %s ]", nrow(res), myfile))
    
    res$chr <- as.numeric(as.character(gsub("_.*", "", res$ibdid)))
    #res$chr <- as.numeric(as.character(sub("chr", "", res$chr)))
    #res <- subset(res, !is.na(chr))
    res$pos <- as.numeric(as.character(gsub(".*_", "", res$ibdid)))
    
    res <- subset(res, GenVar > 0.1)
    if(nrow(res) > 0){
      message(sprintf("###>>> [ %s ] IBD detected with GenVar > 0.1!", nrow(res)))
      res$file <- myfile
      resout <- rbind(resout, res)
    }
  }
  
  resout$trait <- gsub("_.*", "", resout$file)
  resout$mode <- gsub(".*_", "", resout$file)
  resout$tsf <- gsub("_1.*", "", resout$file)
  resout$tsf <- gsub(".*_", "", resout$tsf)
 
  return(resout)
}

###### mrkres ####################################
mrkres <- getmrkres(dir="slurm-scripts/wholeset/")

write.table(mrkres, "cache/mrkres_GenVar_wholeset.csv", sep=",", row.names=FALSE, quote=FALSE)


mrkres <- read.csv("cache/mrkres_GenVar_wholeset.csv")
mrkres <- mrkres[order(mrkres$GenVar, decreasing=TRUE),]
