## Jinliang Yang
## Oct. 13th, 2014
## phenotypic data of offpvp diallel

#setwd("~/Documents/Github/pvpDiallel/")
trait2gs <- function(outdir="largedata/pheno/wholeset/"){
  
  trait <- read.csv("data/trait_matrix_updated_BPH.csv")
  trait$Hyb <- paste(trait$P1, trait$P2, sep="x")
  
  ti <- c("ASI", "DTP", "DTS", "EHT",  "GY", "PHT",  "TW")
  ts <- c("valHyb", "BPHmax", "pBPHmax", "BPHmin", "pBPHmin", "MPH", "pMPH", "pBPH", "BPH")
  
  for(i in 1:7){
    pheno <- subset(trait, trait==ti[i])
    
    ### ===> trait per se
    pheno1 <- pheno[, c(1,4,5)]
    names(pheno1) <- c("Genotype", ti[i], "Fix")
    pheno1$Fix <- 1
    
    out <- paste(outdir, tolower(ti[i]), "_perse.txt", sep="")
    write.table(pheno1, out, sep="\t", row.names=FALSE, quote=FALSE) 
    message(sprintf(">>> output [ %s ] per se pheno to [ %s ]", ti[i], out)) 
    
    ### ==> other traits
    for(j in 2:9){
      pheno1 <- pheno[, c("Hyb", ts[j], "valP1")]
      names(pheno1) <- c("Genotype", ti[i], "Fix")
      pheno1$Fix <- 1
      
      out <- paste(outdir, tolower(ti[i]), "_", ts[j], ".txt", sep="")
      write.table(pheno1, out, sep="\t", row.names=FALSE, quote=FALSE) 
      message(sprintf(">>> output [ %s ] | [ %s ] pheno to [ %s ]", ti[i], ts[j], out)) 
      
    }
  } 
}

trait2gs(outdir="largedata/pheno/wholeset/")
