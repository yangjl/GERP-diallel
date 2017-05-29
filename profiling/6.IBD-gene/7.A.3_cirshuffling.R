### Jinliang
### Jan. 8th, 2015

### circular shuffling
CirShuffling <- function(gerp=gerp, SN=1000000, times=10, outfile="allsnps_11m_gerpv2"){
  
  outfile0 <- paste(outfile, "_cs", 0, ".csv", sep="")
  write.table(gerp, outfile0, sep=",", row.names=FALSE, quote=FALSE)
  message(sprintf(">>> output [ %s ] !!!", outfile0))
  
  for(i in 1:times){
    gerp$RS <- c(gerp$RS[(SN + 1):nrow(gerp)], gerp$RS[1:SN])
    outfile0 <- paste(outfile, "_cs", i, ".csv", sep="")
    write.table(gerp, outfile0, sep=",", row.names=FALSE, quote=FALSE)
    message(sprintf(">>> output [ %s ] !!!", outfile0))
  }
}

#########
library(data.table, lib="~/bin/Rlib/")
gerp <- fread("largedata/SNP/gerp11m_in_gene_b0.csv", sep=",")
# 313821

########## update GERP 
genegeno <- read.csv("largedata/newGERPv2/genegeno/genegeno_cs0.csv")

v3 <- read.csv("largedata/gerpsnp_v3_345176_del.csv")
gv3 <- merge(genegeno, v3[, c("marker", "RS")], by.x="snpid", by.y="marker")
dim(gv3)
# 221960     6


CirShuffling(gerp=gerp, SN=50000, times=10, outfile="largedata/newGERPv2/genegeno/genegeno")


CirShuffling(gerp=gerp, SN=50000, times=10, outfile="largedata/SNP/gene_bph_cs/gene_bph")


## turn all GERP=> 1 no information
gerp1 <- gerp
gerp1$RS <- 1
write.table(gerp1, "largedata/newGERPv2/genegeno/genegeno_cs999.csv", sep=",",
            row.names=FALSE, quote=FALSE)
