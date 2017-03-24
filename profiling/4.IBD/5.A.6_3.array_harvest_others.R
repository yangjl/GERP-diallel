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
  
  infile$cs <- gsub(".*_cs|_.*", "", infile$file)
  infile$trait <- gsub("_.*", "", infile$file)
  infile$mode <- gsub(".*_|.ghat.*", "", infile$file) 
  
  infile$cv <- gsub(".*_train|_sp.*", "", infile$file)
  infile$sp <- gsub(".*_sp|_cs.*", "", infile$file)
  
  infile$type <- "cs"
  infile[infile$cs==0,]$type <- "real"
  infile[infile$cs==999,]$type <- "null"
  
  print(table(infile$trait))
  return(infile)
}

#### extract with real data
collect_res <- function(dir="largedata/newGERPv2/allgeno_perse_d/"){
  res1 <- harvestCV(dir=dir, fileptn="\\.ghatREL")
  res2 <- SplitName(infile=res1) #885
  #rand1 <- subset(rand1, trait != "asi")
  #rand1$trait <- tolower(rand1$trait)
  return(res2)
}

######### trait perse
res_a <- collect_res(dir="largedata/newGERPv2/allgeno_perse_a/") ### done
write.table(res_a[, c("file", "trait", "r", "cs", "mode", "cv", "sp", "type")],
            "largedata/newGERPv2/res_a2_perse_42000.csv", sep=",", row.names=FALSE, quote=FALSE)

res_d <- collect_res(dir="largedata/newGERPv2/allgeno_perse_d/") ### done
write.table(res_d[, c("file", "trait", "r", "cs", "mode", "cv", "sp", "type")],
            "largedata/newGERPv2/res_d2_perse_42000.csv", sep=",", row.names=FALSE, quote=FALSE)

res_k5 <- collect_res(dir="largedata/newGERPv2/allgeno_perse_k5/") ### done
write.table(res_k5[, c("file", "trait", "r", "cs", "mode", "cv", "sp", "type")],
            "largedata/newGERPv2/res_k5_perse_42000.csv", sep=",", row.names=FALSE, quote=FALSE)

res_k <- collect_res(dir="largedata/newGERPv2/allgeno_perse_k/") ### done
write.table(res_k[, c("file", "trait", "r", "cs", "mode", "cv", "sp", "type")],
            "largedata/newGERPv2/res_realk_perse_42000.csv", sep=",", row.names=FALSE, quote=FALSE)


######### trait bph
res_a <- collect_res(dir="largedata/newGERPv2/allgeno_bph_a/") ### done
write.table(res_a[, c("file", "trait", "r", "cs", "mode", "cv", "sp", "type")],
            "largedata/newGERPv2/res_a2_bph_42000.csv", sep=",", row.names=FALSE, quote=FALSE)

res_d <- collect_res(dir="largedata/newGERPv2/allgeno_bph_d/") ### done
write.table(res_d[, c("file", "trait", "r", "cs", "mode", "cv", "sp", "type")],
            "largedata/newGERPv2/res_d2_bph_42000.csv", sep=",", row.names=FALSE, quote=FALSE)

res_k5 <- collect_res(dir="largedata/newGERPv2/allgeno_bph_k5/") ### done
write.table(res_k5[, c("file", "trait", "r", "cs", "mode", "cv", "sp", "type")],
            "largedata/newGERPv2/res_k5_bph_42000.csv", sep=",", row.names=FALSE, quote=FALSE)

res_k <- collect_res(dir="largedata/newGERPv2/allgeno_bph_k/") ### done
write.table(res_k[, c("file", "trait", "r", "cs", "mode", "cv", "sp", "type")],
            "largedata/newGERPv2/res_realk_bph_42000.csv", sep=",", row.names=FALSE, quote=FALSE)

