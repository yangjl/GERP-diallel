### Jinliang Yang
### April 24th, 2015


getghat <- function(dir="slurm-scripts/wholeset/"){
  
  files <- list.files(path = dir, pattern="\\.ghatREL")
  ## file line of the shell file:
  message(sprintf("[ %s ] ghatREL files detected!", length(files)))
  
  resout <- data.frame()
  for(i in 1:length(files)){
    myfile <- paste(dir, files[i], sep="")
    #genotype gHat DTP Fix  meanBias PEV=Var(g/y)   R^2 
    ghat <- read.table(myfile, skip=1, header=FALSE)
    names(ghat) <- c("genotype", "ghat")
    ghat$file <- files[i]
    resout <- rbind(resout, ghat[, c("file", "genotype", "ghat")])
  }
  
  file2 <- list.files(path = dir, pattern="\\.cgrRes")
  ## file line of the shell file:
  message(sprintf("[ %s ] cgrRes files detected!", length(files)))
  
  res2 <- data.frame()
  for(i in 1:length(file2)){
    myfile2 <- paste(dir, file2[i], sep="")
    #genotype gHat DTP Fix  meanBias PEV=Var(g/y)   R^2 
    cgr <- read.table(myfile2, skip=2, header=FALSE)
    cgr$file <- file2[i]
    res2 <- rbind(res2, cgr)
  }
  
  resout$file <- gsub("\\.ghatREL.*", "", resout$file)
  res2$file <- gsub("\\.cgrRes.*", "", res2$file)
  resout <- merge(resout, res2, by="file")
  resout$ghat <- resout$ghat + resout$V2
  return(resout)
}

###### ghat ####################################
ghat <- getghat(dir="slurm-scripts/gerpall_wholeset/")
ghat$trait <- gsub("_.*", "", ghat$file)
ghat$mode <- gsub(".*_", "", ghat$file)
ghat$tsf <- gsub("_..$", "", ghat$file)
ghat$tsf <- gsub(".*_", "", ghat$tsf)

write.table(ghat, "cache/gerpall_ghat.csv", sep=",", row.names=FALSE, quote=FALSE)

