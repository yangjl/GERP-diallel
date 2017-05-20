### Jinliang Yang
### 05-20-2017, I love you.
### purpose: response with reviewer 2, using random k values


library(data.table)

dsf <- fread("largedata/SNP/allsnps2_newgerp2_50k.dsf7", data.table=FALSE)

### 100kb bin
dsf$bin <- paste(dsf$chr, round(dsf$pos/100000, 0), sep="_")
length(unique(dsf$bin))
k1 <- dsf[!duplicated(dsf$bin),]
dim(k1)
## 15950
write.table(k1[, -ncol(k1)], "largedata/sgeno_100k/allsnps2_newgerp2_100k.dsf7", sep="\t", row.names=FALSE, quote=FALSE)



########## shuffling k=d/a 10 times for each trait
set.seed(1234567)
kfile <- list.files(path="largedata/snpeff/perse", pattern="_k.txt$", full.names=TRUE)

inputdf1 <- data.frame(
  d="largedata/IBD/allsnps_11m_IBD.bed", 
  s="largedata/sgeno_100k/allsnps2_newgerp2_100k.dsf7", #define deleterious alleles
  g= "largedata/newGERPv2/allgeno/gerpv2_b0_cs0.csv",
  f= kfile, # degree of dominance, k values
  out= gsub("snpeff\\/perse", "sgeno_100k", gsub(".txt", "", kfile)),
  l=0,
  t="k"
)

library(farmeR)
run_gerpIBD(inputdf1, email="yangjl0930@gmail.com", runinfo = c(FALSE, "med", 2) )
###>>> In this path: cd /home/jolyang/Documents/Github/GERP-diallel
###>>> RUN: sbatch -p bigmemm --mem 8G --ntasks=1 --time=8:00:00 slurm-script/run_gerpibd_array.sh




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
