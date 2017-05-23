# Jinliang Yang
# 05-22-2017
# harvest the results of model training with gerp and random SNPs

harvest_gsout <- function(dir="largedata/sgeno", fileptn="\\.out"){
    files <- list.files(path = dir, pattern = fileptn, full.names=TRUE)
    ## file line of the shell file:
    message(sprintf("[ %s ] files detected!", length(files)))
    
    out <- data.frame(file=files, genvar=-9, resvar=-9, totvar=-9, h2=-9)
    for(i in 1:nrow(out)){
        con <- file(as.character(out$file[i]), "r")
        text <- readLines(con)
        val1 <- grep("Residual Variance", text, value=TRUE)
        val1 <- gsub(".* = ", "", val1)
        
        val2 <- grep("Genetic  Variance", text, value=TRUE)
        val2 <- gsub(".* = ", "", val2)
        
        val3 <- grep("Total Variance", text, value=TRUE)
        val3 <- gsub(".* = ", "", val3)
        
        val4 <- grep("Proportion of Variance a/c Markers", text, value=TRUE)
        val4 <- gsub(".* = ", "", val4)
        
        if(is.null(val1)) val1 = -9
        if(is.null(val2)) val2 = -9
        if(is.null(val3)) val3 = -9
        if(is.null(val4)) val4 = -9
        
        out$resvar[i] <- val1
        out$genvar[i] <- val2
        out$totvar[i] <- val3
        out$h2[i] <- val4
        close(con)
    }
  return(out)
}


########
SplitName <- function(infile=resout){
  
  infile$file <- as.character(infile$file)
  infile$file <- gsub(".*\\/", "", infile$file)
  
  infile$stime <- gsub(".*_stime|_.*", "", infile$file)
  infile$trait <- gsub("_.*", "", infile$file)
  return(infile)
}

res1 <- harvest_gsout(dir="largedata/sgeno", fileptn="\\.out")
res2 <- SplitName(infile=res1) #885

write.table(res2, "largedata/var_explained_shuffling.csv", 
            sep=",", row.names=FALSE, quote=FALSE)

#### 100-kb
res1 <- harvest_gsout(dir="slurm-script/subgeno", fileptn="\\.out")
res2 <- SplitName(infile=res1) #885
res2$trait <- gsub("ws_|_h2.*", "", res2$file)
write.table(res2, "largedata/var_explained_100kb.csv", 
            sep=",", row.names=FALSE, quote=FALSE)
