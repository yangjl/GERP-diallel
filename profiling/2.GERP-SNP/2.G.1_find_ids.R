### Jinliang Yang
### 10/26/2016
### purpose: find lines in hmp3

##### get landrace lines
get_founder <- function(){
  pvp <- read.csv("largedata/GERPv2/gerpsnp_506898.csv", nrow=5)
  ids <- names(pvp)[10:21]
  
  line <- read.table("~/dbcenter/HapMap/HapMap3/vcf.header")
  newline <- toupper(as.vector(t(line)))
  idx <- grep("BKN", newline)
  bkn <- line[idx]
  
  idx2 <- which(newline %in% ids)
  p <- line[idx2]
  
  write.table(data.frame(v1=c(bkn, p)), "~/dbcenter/HapMap/HapMap3/bkn_pvp_samples.txt", 
              sep="\n", row.names=FALSE, quote=FALSE, col.names=FALSE)
}
######
get_founder()











