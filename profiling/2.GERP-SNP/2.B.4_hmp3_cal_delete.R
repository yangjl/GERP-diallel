### Jinliang Yang
### deterimine the deleterious alleles


##############
library(tidyr)
library(plyr)
library(pbapply)
library(data.table)

cal_del_allele <- function(chri = 10, outbase="largedata/Alignment/hmp3_GERPv3_major"){
  
  ### reading data
  res1 <- fread(paste0("largedata/Alignment/hmp3_chr", chri, "_gerpv3.txt"), header=FALSE, data.table=FALSE)
  
  res0 <- subset(res1, V2 != "-" & V2 != "N")
  res0$sp <- gsub("-.*", "", res0$V1)
  res0$V1 <- gsub(".*-", "", res0$V1)
  
  snpid <- table(res0$V1)
  t2 <- snpid[snpid > 1]
  res2 <- subset(res0, V1 %in% names(t2))
  res2 <- spread(res2, sp, V2)
  
  ### calculate maf
  message(sprintf("###>>> calculating major allele for chr [%s] ... ", chri))
  
  pboptions(type="txt", char="-")
  outfile <- paste0(outbase, "_chr", chri, ".csv")
  write.table(data.frame(V1="snpid", V2="major", V3="majc", V4="minor", V5="minc"), outfile, sep=",", 
              row.names=FALSE, col.names=FALSE, quote=FALSE, append =TRUE)
  ma <- pbapply(res2, 1, function(x){
    n <- table(t(x[-1]))
    #idx <- which(names(n) == "-")
    #if(length(idx) > 0) n <- n[-idx]
    if(length(n) == 1){
      n <- sort(n, decreasing=TRUE)
      out <- data.frame(snpid=x[1], major=names(n)[1], majc=n[1], minor="N", minc=0)
    }else if(length(n) >=2){
      n <- sort(n, decreasing=TRUE)
      out <- data.frame(snpid=x[1], major=names(n)[1], majc=n[1], minor=names(n)[2], minc=n[2])
    }
    write.table(out, outfile, sep=",", row.names=FALSE, col.names=FALSE, quote=FALSE, append =TRUE)
  })
}

##get command line args
options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)

myi <- as.numeric(as.character(args[1]))

out <- cal_del_allele(chri = myi, outbase="largedata/Alignment/hmp3_GERPv3_major")


