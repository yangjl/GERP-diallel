### Jinliang Yang
### 09/08/2015
### updated: 03-23-2017



##########################################
getk <- function(files, outpwd="largedata/snpeff"){
  output <- data.frame()
  for(i in 1:length(files)){
    h1 <- fread(files[i], header=TRUE, data.table=FALSE)
    names(h1) <- c("snpid","chr","pos","Effect_A","Effect_D","Effect_A2","Effect_D2","h2_mrk_A", 
                   "h2_mrk_D","H2_mrk","h2_mrk_A_p","h2_mrk_D_p","H2_mrk_p","log10_h2_mrk_A","log10_h2_mrk_D","log10_H2_mrk")
    #h1 <- subset(h1, Effect_A !=0)
    #h1 <- subset(h1, H2_mrk > 1e-7 & Effect_A !=0)
    
    trait <- gsub(".*\\/|_snpeff_.*", "", files[i])
    h1$k <- h1$Effect_D/h1$Effect_A
    
    #out <- h1[, c("snpid", "k")]
    out <- h1[!is.na(h1$k),]

    if(sum(out$k > 1) > 0){
      if(sum(out$k > 2) > 0){
        out[out$k > 2, ]$k <- 2.1
      }
      #out[out$k > 1, ]$k <- rescale(out[out$k > 1, ]$k, c(1, 2))
    }
    if(sum(out$k < -1) > 0){
      if(sum(out$k < -2) > 0){
        out[out$k < -2, ]$k <- -2.1
      }
      #out[out$k < -1, ]$k <- rescale(out[out$k < -1, ]$k, c(-2, -1))
    }
    
    message(sprintf("###>>> trait [ %s ], snp # [ %s ], k ranged [ %s - %s ]", trait, nrow(out), max(out$k), min(out$k)))
    write.table(out, paste0(outpwd, "/", trait, "_k.txt"), sep="\t", row.names=FALSE, quote=FALSE)
  }
}
##############
library(data.table)
files1 <- list.files(path="largedata/GBLUP", pattern="snpe$", full.names=TRUE)
getk(files=files1, outpwd="largedata/lcache/")

file2 <- list.files(path="largedata/snpeff/BPH", pattern="snpe$", full.names=TRUE)
getk(file2, outpwd="largedata/snpeff/BPH")

