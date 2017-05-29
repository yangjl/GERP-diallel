# Jinliang Yang
# updated: April/24/2015
# Purpose: quick plot of GWAS results

getBayes <- function(inputfile="/home/NSF-SAM-GS/GenSel/SAM_run41000.mrkRes1", mfcutoff=0){
  
  tab5rows <- read.table(inputfile, header=TRUE, nrows=5)
  classes <- sapply(tab5rows, class)
  res <- read.table(inputfile, header=TRUE, colClasses=classes)
  message(sprintf("###>>> input [ %s ]", nrow(res)))
  
  res$chr <- as.numeric(as.character(gsub("_.*", "", res$ibdid)))
  #res$chr <- as.numeric(as.character(sub("chr", "", res$chr)))
  #res <- subset(res, !is.na(chr))
  res$pos <- as.numeric(as.character(gsub(".*_", "", res$ibdid)))
  res <- subset(res, ModelFreq > mfcutoff)
  message(sprintf("removing >mfcutff, remainning [ %s ]", nrow(res)))
  return(res)
}


source("lib/quickMHTplot.R")

#### pai=0.995
bayes1 <- getBayes(inputfile="slurm-scripts/wholeset/gy_perse_1a2.mrkRes1", mfcutoff=0)
#bayes1s <- subset(bayes1, ModelFreq > 0.005)
#pdf("largedata/lgraphs/asi_gs.pdf", width=10, height=4)
quickMHTplot(res=bayes1, main="ASI, pai=0.9999, chain length=101000", cex=.9, pch=16, 
             col=rep(c("slateblue", "cyan4"), 5), 
             GAP=5e+06, ylab="model frequency", yaxis=NULL,
             col2plot="ModelFreq")


