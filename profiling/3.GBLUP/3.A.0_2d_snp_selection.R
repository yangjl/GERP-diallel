### Jinliang Yang
### date: 10-03-2016
### purpose: 2D SNP selection for GWAS


# Note: snp11m is filtered data MAF > 0.1 and 
data_cleaning <- function(){
  ### SNP data
  library("data.table")
  snp11m <- fread("largedata/SNP/allsnps_11m_gerpv2.csv", header=TRUE, sep=",") 
  snp11m <- as.data.frame(snp11m)
  snp11m[is.na(snp11m$RS), ]$RS <- 0
  
  ### cM
  map <- fread("largedata/SNP/allsnps_11m_genetic.map")
  map$AGPv2_pos <- round(map$AGPv2_pos/1000000, 0)
  map <- as.data.frame(map)
  names(map) <- c("snpid", "chr", "genetic")
  
  ### Exonic bp per cM => quantiles
  excm <- read.csv("cache/exonic_cM.csv")
  excm$qt0 <- 9
  excm[excm$exonbp <= quantile(excm$exonbp)[2], ]$qt0 <- 1
  excm[excm$exonbp > quantile(excm$exonbp)[2] & excm$exonbp <= quantile(excm$exonbp)[3], ]$qt0 <- 2
  excm[excm$exonbp > quantile(excm$exonbp)[3] & excm$exonbp <= quantile(excm$exonbp)[4], ]$qt0 <- 3
  excm[excm$exonbp > quantile(excm$exonbp)[4], ]$qt0 <- 4
  
  map$cM <- paste(map$chr, map$genetic, sep="_")
  map2 <- merge(map[, -2], excm[, c("cM", "exonbp", "qt0")], by="cM")
  
  snpdf <- merge(map2, snp11m, by="snpid")
  # round to nearest .05 or .10
  snpdf$frq <- round(snpdf$MAF*2, 1) / 2
  
  return(snpdf)
  
}


#############################################################################
##### generate same number of random SNPs
getRandomSNP <- function(snpdf, verbose=FALSE, outfile="largedata/SNP/randomsnp/rsnp1.csv"){
  
  ###
  sub <- subset(snpdf, RS > 0)
  myr <- subset(snpdf, RS <= 0)
  out <- data.frame()
  a <- c(0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50)
  for(q in 1:4){
    for (f in a){
      num <- nrow(subset(sub, qt0 == q & frq == f))
      subr <- subset(myr, qt0 == q & frq == f)
      
      if(verbose){
        message(sprintf("### [q= %s], [f=%s]: target num: [%s] and pools: [%s]",
                        q, f, num, nrow(subr)))
      }
      
      idx <- sample(1:nrow(subr), num, replace=FALSE)
      out <- rbind(out, subr[idx, ])
    }
  }
  
  #write.table(sub, "largedata/SNP/randomsnp/rsnp0.csv", sep=",", row.names=FALSE, quote=FALSE)
  write.table(out, outfile, sep=",", row.names=FALSE, quote=FALSE)
}

###################################
set.seed(1234567)
snpdf <- data_cleaning()

for(i in 10:10){
  getRandomSNP(snpdf, verbose=FALSE, 
               outfile= paste0("largedata/SNP/randomsnp/rsnp", i, ".csv"))
}
