### Jinliang Yang
### 1.5.2014
### purpose: run the gerpIBD program


###############
gerpfile <- list.files(path="largedata/newGERPv2/allgeno", pattern="csv$", full.names=TRUE)
inputdf <- data.frame(
  d="largedata/IBD/allsnps_11m_IBD.bed", 
  s="largedata/SNP/allsnps2_newgerp2_50k.dsf7", #deine deleterious alleles
  g=gerpfile, 
  f="largedata/newGERPv2/allgeno/k5.txt",
  out=gsub(".csv", "", gerpfile),
  l=0,
  t="all"
)

library(farmeR)
run_gerpIBD(inputdf, email="yangjl0930@gmail.com", runinfo = c(TRUE, "med", 2) )

############### trait-specific k
gerpfile <- list.files(path="largedata/newGERPv2/allgeno", pattern="csv$", full.names=TRUE)
kfile <- list.files(path="largedata/snpeff/perse", pattern="_k.txt$", full.names=TRUE)
inputdf1 <- data.frame(
  d="largedata/IBD/allsnps_11m_IBD.bed", 
  s="largedata/SNP/allsnps2_newgerp2_50k.dsf7", #deine deleterious alleles
  g= rep(gerpfile, times = 7),
  f= rep(kfile, each=12),
  out=gsub("allgeno", "allgeno_k", gsub(".csv", "", gerpfile)),
  l=0,
  t="k"
)

inputdf1$out <- paste0(inputdf1$out, "_", gsub(".*/|_k.txt", "", inputdf1$f), "_perse")

#library(farmeR)
#run_gerpIBD(inputdf, email="yangjl0930@gmail.com", runinfo = c(TRUE, "med", 2) )


############### trait-specific k for BPH
gerpfile <- list.files(path="largedata/newGERPv2/allgeno", pattern="csv$", full.names=TRUE)
kfile <- list.files(path="largedata/snpeff/BPH", pattern="_k.txt$", full.names=TRUE)
inputdf2 <- data.frame(
  d="largedata/IBD/allsnps_11m_IBD.bed", 
  s="largedata/SNP/allsnps2_newgerp2_50k.dsf7", #deine deleterious alleles
  g= rep(gerpfile, times = 7),
  f= rep(kfile, each=12),
  out=gsub("allgeno", "allgeno_kbph", gsub(".csv", "", gerpfile)),
  l=0,
  t="k"
)

inputdf2$out <- paste0(inputdf2$out, "_", gsub(".*/|_k.txt", "", inputdf2$f), "_bph")

library(farmeR)
run_gerpIBD(inputdf=rbind(inputdf1, inputdf2), email="yangjl0930@gmail.com", runinfo = c(TRUE, "med", 2) )


library(farmeR)
run_gerpIBD(inputdf=inputdf2, email="yangjl0930@gmail.com", runinfo = c(TRUE, "serial", 4) )





##############################################
