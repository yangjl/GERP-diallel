### Jinliang
### Jan. 8th, 2015

### circular shuffling
CirShuffling <- function(gerp=gerp, SN=1000000, times=10, outfile="allsnps_11m_gerpv2"){
  for(i in 1:times){
    gerp$RS <- c(gerp$RS[(SN + 1):nrow(gerp)], gerp$RS[1:SN])
    outfile0 <- paste(outfile, "_cs", i, ".csv", sep="")
    write.table(gerp, outfile0, sep=",", row.names=FALSE, quote=FALSE)
    message(sprintf("###>>> output [ %s ] !!!", outfile0))
  }
}

#########
library(data.table, lib="~/bin/Rlib/")
gerp <- fread("largedata/SNP/allsnps_11m_gerpv2_tidy.csv", sep=",")
gerp <- subset(gerp, RS>0) #506898      5
write.table(gerp, "largedata/newGERPv2/allgeno/gerpv2_b0_cs0.csv", sep=",",
            row.names=FALSE, quote=FALSE)

CirShuffling(gerp=gerp, SN=50000, times=10, outfile="largedata/newGERPv2/allgeno/gerpv2_b0")

## turn all GERP=> 1 no information
gerp1 <- gerp
gerp1$RS <- 1
write.table(gerp1, "largedata/newGERPv2/allgeno/gerpv2_b0_cs999.csv", sep=",",
            row.names=FALSE, quote=FALSE)







#################### circular shuffling
CirShuffling_method2 <- function(allgerp, mygerp, SN=1000000, times=10, outfile="allsnps_11m_gerpv2"){
  for(i in 1:times){
    allgerp$RS <- c(allgerp$RS[(SN + 1):nrow(allgerp)], allgerp$RS[1:SN])
    outfile0 <- paste(outfile, "_cs", i, ".csv", sep="")
    gerp <- subset(allgerp, snpid %in% mygerp$snpid)
    write.table(gerp, outfile0, sep=",", row.names=FALSE, quote=FALSE)
    message(sprintf(">>> output [ %s ] !!!", outfile0))
  }
}

library(data.table, lib="~/bin/Rlib/")
allgerp <- fread("largedata/SNP/allsnps_11m_gerpv2_tidy.csv", sep=",")
allgerp <- as.data.frame(allgerp)

gerp0 <- subset(allgerp, RS > 0) #167455      5
CirShuffling_method2(allgerp, mygerp=gerp0, SN=50000, times=10, outfile="largedata/SNP/geno_b0_cs/gerpv2_b1")


gerp1 <- subset(allgerp, RS > 1) #167455      5
CirShuffling_method2(allgerp, mygerp=gerp1, SN=50000, times=10, outfile="largedata/SNP/geno_b1_cs/gerpv2_b1")


gerp2 <- subset(allgerp, RS > 2) #28977     5
CirShuffling_method2(allgerp, mygerp=gerp2, SN=50000, times=10, outfile="largedata/SNP/geno_b2_cs/gerpv2_b2")

