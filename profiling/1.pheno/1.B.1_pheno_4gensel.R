## Jinliang Yang
## 03-19-2017
## phenotypic data of offpvp diallel

#setwd("~/Documents/Github/pvpDiallel/")
trait2gs <- function(traitfile = "data/hyb_heterosis.csv", outdir="largedata/pheno/wholeset/"){
  
  trait <- read.csv(traitfile)
  trait$Hyb <- paste(trait$Par1, trait$Par2, sep="x")
  dir.create(outdir, showWarnings = FALSE)
  
  ti <- c("ASI", "DTP", "DTS", "EHT",  "GY", "PHT",  "TW")
  #ts <- c("valHyb", "BPHmax", "pBPHmax", "BPHmin", "pBPHmin", "MPH", "pMPH", "pBPH", "BPH")
  ts <- c("valHyb", "pMPH", "MPH")
  for(i in 1:7){
    pheno <- subset(trait, trait==ti[i])
    
    ### ===> trait per se
    pheno1 <- pheno[, c(1,4,5)]
    names(pheno1) <- c("Genotype", ti[i], "Fix")
    pheno1$Fix <- 1
    
    out <- paste(outdir, "/", tolower(ti[i]), "_perse.txt", sep="")
    write.table(pheno1, out, sep="\t", row.names=FALSE, quote=FALSE) 
    message(sprintf(">>> output [ %s ] per se pheno to [ %s ]", ti[i], out)) 
    
    ### ==> other traits
    for(j in 2:3){
      pheno1 <- pheno[, c("Hyb", ts[j], "valPar1")]
      names(pheno1) <- c("Genotype", ti[i], "Fix")
      pheno1$Fix <- 1
      
      out <- paste(outdir, "/", tolower(ti[i]), "_", ts[j], ".txt", sep="")
      write.table(pheno1, out, sep="\t", row.names=FALSE, quote=FALSE) 
      message(sprintf(">>> output [ %s ] | [ %s ] pheno to [ %s ]", ti[i], ts[j], out)) 
      
    }
  } 
}

trait2gs(traitfile = "data/hyb_heterosis.csv", outdir="largedata/pheno/wholeset/")
