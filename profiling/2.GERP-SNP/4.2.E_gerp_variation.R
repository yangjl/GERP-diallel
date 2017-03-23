### Jinliang Yang
### 9/3/2015
### purpose: find the variation of the GERP

library("plyr")
library("data.table", lib="~/bin/Rlib/")
#load data in the local machine for plotting
ob <- load("largedata/lcache/snpnzRS.RData")
#snpnz

### genotyping information
snp11m <- fread("largedata/SNP/allsnps_11m.dsf5")
snp11m <- as.data.table(snp11m)

##############################################################
source("lib/use_p2g.R")

get_map <- function(){
  gerp_b0 <- as.data.frame(subset(snpnz, RS >0)) #506898
  gerp_b0 <- merge(gerp_b0[, 1:2], snp11m, by="snpid")[, c("snpid", "chr", "pos")]
  names(gerp_b0)[3] <- "Physical"
  
  ####
  test_p2g(chr=9)
  map <- use_p2g(df=gerp_b0)
  message(sprintf("###>>> Got the genetic position [ %s ] SNP using GAM function", nrow(map)))
  
  return(map)
}

map <- get_map()


gerp_b0 <- as.data.frame(subset(snpnz, RS >0)) #506898
gerpsnp <- merge(gerp_b0[, 1:2], snp11m, by="snpid")
gerpsnp <- merge(map[, c(1,4)], gerpsnp, by.x="marker", by.y="snpid")
gerpsnp <- gerpsnp[order(gerpsnp$chr, gerpsnp$pos), ]

write.table(gerpsnp, "largedata/GERPv2/gerpsnp_506898.csv", sep=",", row.names=FALSE, quote=FALSE)
#note gerpsnp <- subset(gerpsnp, B73 != "N")
