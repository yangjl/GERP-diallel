# Jinliang Yang
# Jan 12th, 2015
# harvest the results of model training with gerp and random SNPs

harvestCV <- function(dir="slurm-scripts/", fileptn="\\.ghatREL1"){
  
  files <- list.files(path = dir, pattern=fileptn, full.names=TRUE)
  ## file line of the shell file:
  message(sprintf("[ %s ] files detected!", length(files)))
  
  res <- unlist(lapply(1:length(files), function(i){
    #genotype gHat DTP Fix  meanBias PEV=Var(g/y)   R^2 
    ghat <- read.table(files[i], skip=1, header=FALSE)
    if(i %% 10000 ==0) {
      # Print on the screen some message
      message(sprintf("###>>> finished reading [ %s ] files!", i))
    }
    return(cor(ghat$V2, ghat$V3))
  }))
  
  resout <- data.frame(file=files, r=res)
  return(resout)
}

########
SplitName <- function(infile=resout){
  
  infile$file <- as.character(infile$file)
  infile$file <- gsub(".*/", "", infile$file)
  
  infile$cs <- gsub("_.*", "", infile$file)
  infile$trait <- gsub("_cv.*", "", infile$file)
  infile$trait <- gsub("_.2", "", infile$trait)
  infile$trait <- gsub(".*_", "", infile$trait)
  
  infile$mode <- gsub("_cv.*", "", infile$file) 
  infile$mode <- gsub(".*_", "", infile$mode)
  infile$cv <- paste0("cv", gsub(".*_cv|_.*", "", infile$file))
  infile$sp <- paste0("sp", gsub(".*_sp|\\..*", "", infile$file))
  infile$rel <- gsub(".*\\.", "", infile$file)
  
  infile$type <- infile$cs
  infile$type <- gsub("cs0", "real", infile$type)
  infile$type <- gsub("cs.*", "random", infile$type)
  
  print(table(infile$trait))
  return(infile)
}

#### extract with real data
collect_res <- function(dir="slurm-scripts/cv_b2/"){
  res1 <- harvestCV(dir=dir, fileptn="\\.ghatREL")
  res1 <- SplitName(infile=res1) #885
  print(table(res1$trait))
  
  #rand1 <- subset(rand1, trait != "asi")
  #rand1$trait <- tolower(rand1$trait)
  return(res1)
  
}

library(plyr, lib="~/bin/Rlib/")
################################################################
#7*3*5*100*11 = [1] 115500

gene_perse0 <- collect_res(dir="largedata/SNP/gene_perse_cs/")
write.table(gene_perse0, "cache/gene_perse_2016.csv", sep=",", row.names=FALSE, quote=FALSE)

###############################
gene_bph0 <- collect_res(dir="largedata/SNP/gene_bph_cs/")

write.table(gene_bph0, "cache/gene_bph_2016.csv", sep=",", row.names=FALSE, quote=FALSE)

