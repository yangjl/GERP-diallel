## Jinliang Yang
## Oct. 13th, 2014
## phenotypic data of offpvp diallel

#setwd("~/Documents/Github/pvpDiallel/")
trait2gsCV <- function(){
  
  trait <- read.csv("data/trait_matrix_updated_pBPH.csv")
  trait$Hyb <- paste(trait$P1, trait$P2, sep="x")
  
  ti <- c("ASI", "DTP", "DTS", "EHT",  "GY", "PHT",  "TW")
  parent <- unique(c(as.character(trait$P1), as.character(trait$P2)))
  
  for(i in 1:7){
    pheno <- subset(trait, trait==ti[i])
   
    for(j in 1:length(parent)){
      #### training phenotype
      trainpheno <- subset(pheno, P1!=parent[j] & P2!=parent[j])
      pheno1 <- trainpheno[, c(1,4,5)]
      names(pheno1) <- c("Genotype", ti[i], "Fix")
      pheno1$Fix <- 1
      #### validation phenotype
      testpheno <- subset(pheno, P1==parent[j] | P2==parent[j])
      pheno2 <- testpheno[, c(1,4,5)]
      names(pheno2) <- c("Genotype", ti[i], "Fix")
      pheno2$Fix <- 1
      #### output phenotype
      out1 <- paste("largedata/pheno/CV/", tolower(ti[i]), "_train", j, ".txt", sep="")
      write.table(pheno1, out1, sep="\t", row.names=FALSE, quote=FALSE) 
      out2 <- paste("largedata/pheno/CV/", tolower(ti[i]), "_test", j, ".txt", sep="")
      write.table(pheno2, out2, sep="\t", row.names=FALSE, quote=FALSE) 
      message(sprintf(">>> output pheno [ %s ], set [%s] to [ %s ] and [ %s ]", ti[i],j, out1, out2))
    }
  } 
}

####
system("mkdir largedata/pheno/CV/")
seed(1234567)
trait2gsCV()
