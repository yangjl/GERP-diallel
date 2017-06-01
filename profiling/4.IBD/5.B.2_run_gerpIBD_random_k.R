### Jinliang Yang
### 05-20-2017, I love you.
### purpose: response with reviewer 2, using random k values


library(data.table)
shuffling_k <- function(kfile, stimes){
    
    # kfile <- list.files(path="largedata/snpeff/perse", pattern="_k.txt$", full.names=TRUE)
    # stimes: number of shuffling times
    
    for(i in 1:length(kfile)){
        kvalue <- fread(kfile[i], header=TRUE, data.table=FALSE)
        
        outfile <- gsub("_k", paste0("_k_stime", 0), kfile[i])
        write.table(kvalue, outfile, sep="\t", row.names=FALSE, quote=FALSE)
        
        for(j in 1:stimes){
            out <- kvalue
            out$h <- sample(kvalue$h, nrow(kvalue), replace = FALSE)
            outfile <- gsub("_k", paste0("_k_stime", j), kfile[i])
            write.table(out, outfile, sep="\t", row.names=FALSE, quote=FALSE)
        }
    }
}

########## shuffling k=d/a 10 times for each trait
set.seed(1234567)
kfile <- list.files(path="largedata/snpeff/perse", pattern="_k.txt$", full.names=TRUE)
shuffling_k(kfile, stimes=10)


############### trait-specific k
kout <- list.files(path="largedata/snpeff/perse", pattern="_k_stime", full.names=TRUE)
#gerpfile <- list.files(path="largedata/newGERPv2/allgeno", pattern="cs0.csv$", full.names=TRUE)

inputdf1 <- data.frame(
  d="largedata/IBD/allsnps_11m_IBD.bed", 
  s="largedata/SNP/allsnps2_newgerp2_50k.dsf7", #define deleterious alleles
  g= "largedata/newGERPv2/allgeno/gerpv2_b0_cs0.csv",
  f= kout, # degree of dominance, k values
  out=gsub("snpeff\\/perse", "sgeno", gsub(".txt", "", kout)),
  l=0,
  t="k"
)

library(farmeR)
run_gerpIBD(inputdf1, email="yangjl0930@gmail.com", runinfo = c(FALSE, "bigmemm", 1) )
###>>> In this path: cd /home/jolyang/Documents/Github/GERP-diallel
###>>> RUN: sbatch -p bigmemm --mem 8G --ntasks=1 --time=8:00:00 slurm-script/run_gerpibd_array.sh

############### trait-specific adk
kout <- list.files(path="largedata/snpeff/perse", pattern="_k_stime0", full.names=TRUE)
#gerpfile <- list.files(path="largedata/newGERPv2/allgeno", pattern="cs0.csv$", full.names=TRUE)

inputdf1 <- data.frame(
    d="largedata/IBD/allsnps_11m_IBD.bed", 
    s="largedata/SNP/allsnps2_newgerp2_50k.dsf7", #define deleterious alleles
    g= "largedata/newGERPv2/allgeno/gerpv2_b0_cs0.csv",
    f= kout, # degree of dominance, k values
    out=gsub("snpeff\\/perse", "sgeno", gsub(".txt", "", kout)),
    l=0,
    t="adk"
)

library(farmeR)
run_gerpIBD(inputdf1, email="yangjl0930@gmail.com", runinfo = c(FALSE, "bigmemm", 1) )
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
