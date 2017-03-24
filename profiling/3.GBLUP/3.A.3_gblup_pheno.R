### Jinliang Yang
### Sept 4th, 2015
### get degree of dominance



### get the pedigree info
get_pheno_mx <- function(){
  files <- list.files(path="largedata/pheno/wholeset/", pattern="txt", full.names=TRUE)
  trait <- read.table(files[1], header=TRUE)
  names(trait)[2] <-  gsub(".*//|.txt", "", files[1]) 
  trait <- trait[, -3]
  for(i in 2:length(files)){
    tem <- read.table(files[i], header=TRUE)
    names(tem)[2] <-  gsub(".*//|.txt", "", files[i]) 
    tem <- tem[, -3]
    trait <- merge(trait, tem, by="Genotype")
  }
  names(trait)[1] <- "Sample_ID"
  write.table(trait, "largedata/pheno/wholeset/trait_mx.dat", sep="\t", row.names=FALSE, quote=FALSE)
  return(trait)
}

trait <- get_pheno_mx()
idx <- which(names(trait) =="gy_perse")
