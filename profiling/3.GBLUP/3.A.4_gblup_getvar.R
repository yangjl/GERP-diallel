### Jinliang Yang
### Sept 19th, 2015

library("data.table")
library("plyr")
##########################################
getvar <- function(gblup_efile="largedata/snpeff/rsnp1/gy_perse_snpeff_ce.snpe",
                 genofile="largedata/SNP/randomsnp/rsnp1.csv", nf=1){
  
  h1 <- fread(gblup_efile, header=TRUE, data.table=FALSE)
  #h1 <- as.data.frame(h1)
  names(h1) <- c("snpid","chr","pos","Effect_A","Effect_D","Effect_A2","Effect_D2","h2_mrk_A", 
                 "h2_mrk_D","H2_mrk","h2_mrk_A_p","h2_mrk_D_p","H2_mrk_p","log10_h2_mrk_A","log10_h2_mrk_D","log10_H2_mrk")
  #h1 <- subset(h1, Effect_A !=0)
  h1 <- subset(h1, H2_mrk > nf*mean(h1$H2_mrk))
  message(sprintf("###>>> remaining SNPs [ %s ]", nrow(h1)))
  
  ### get the parental genotype info
  geno <- fread(genofile, data.table=FALSE)
  geno <- geno[, c("snpid", "frq", "qt0", "genetic", "exonbp")]
  
  hgeno <- merge(geno, h1, by="snpid")
  tab1 <- ddply(hgeno, .(frq), summarise,
               totvar=sum(H2_mrk))
  
  tab2 <- ddply(hgeno, .(frq), nrow)
  tab <- merge(tab1, tab2, by="frq")
  tab$nvar <- tab$totvar/tab$V1*10000
  
  return(tab)
}

###################################
res0 <- getvar(gblup_efile="largedata/snpeff/rsnp0/gy_perse_snpeff_ce.snpe",
               genofile="largedata/SNP/randomsnp/rsnp0.csv", nf=1)
write.table(res0, "cache/rsnp_var0.csv", sep=",", row.names=FALSE, quote=FALSE)

res1 <- getvar(gblup_efile="largedata/snpeff/rsnp1/gy_perse_snpeff_ce.snpe",
              genofile="largedata/SNP/randomsnp/rsnp1.csv", nf=1)
write.table(res1, "cache/rsnp_var1.csv", sep=",", row.names=FALSE, quote=FALSE)

### collect all data
out <- data.frame()
for(i in 0:10){
  res1 <- getvar(gblup_efile=paste0("largedata/snpeff/rsnp", i, "/gy_perse_snpeff_ce.snpe"),
                 genofile=paste0("largedata/SNP/randomsnp/rsnp", i, ".csv"), nf=1)
  res1$sample <- i
  out <- rbind(out, res1)
}
write.table(out, "cache/rsnp_var_nf1.csv", sep=",", row.names=FALSE, quote=FALSE)

### collect all data
out5 <- data.frame()
for(i in 0:10){
  res1 <- getvar(gblup_efile=paste0("largedata/snpeff/rsnp", i, "/gy_perse_snpeff_ce.snpe"),
                 genofile=paste0("largedata/SNP/randomsnp/rsnp", i, ".csv"), nf=5)
  res1$sample <- i
  out5 <- rbind(out5, res1)
}
write.table(out5, "cache/rsnp_var_nf5.csv", sep=",", row.names=FALSE, quote=FALSE)




# heterosis trait: BPHmax ----------------------------
### collect all data
out1 <- data.frame()
for(i in 0:10){
  res1 <- getvar(gblup_efile=paste0("largedata/snpeff/rsnp", i, "/gy_BPHmax_snpeff_ce.snpe"),
                 genofile=paste0("largedata/SNP/randomsnp/rsnp", i, ".csv"), nf=1)
  res1$sample <- i
  out1 <- rbind(out1, res1)
}
write.table(out1, "cache/rsnp_var_BPH_nf1.csv", sep=",", row.names=FALSE, quote=FALSE)

### collect all data
out5 <- data.frame()
for(i in 0:10){
  res1 <- getvar(gblup_efile=paste0("largedata/snpeff/rsnp", i, "/gy_BPHmax_snpeff_ce.snpe"),
                 genofile=paste0("largedata/SNP/randomsnp/rsnp", i, ".csv"), nf=5)
  res1$sample <- i
  out5 <- rbind(out5, res1)
}
write.table(out5, "cache/rsnp_var_BPH_nf5.csv", sep=",", row.names=FALSE, quote=FALSE)

### collect all data
out4 <- data.frame()
for(i in 0:10){
  res1 <- getvar(gblup_efile=paste0("largedata/snpeff/rsnp", i, "/gy_BPHmax_snpeff_ce.snpe"),
                 genofile=paste0("largedata/SNP/randomsnp/rsnp", i, ".csv"), nf=4)
  res1$sample <- i
  out4 <- rbind(out4, res1)
}
write.table(out4, "cache/rsnp_var_BPH_nf4.csv", sep=",", row.names=FALSE, quote=FALSE)

### collect all data
out3 <- data.frame()
for(i in 0:10){
  res1 <- getvar(gblup_efile=paste0("largedata/snpeff/rsnp", i, "/gy_BPHmax_snpeff_ce.snpe"),
                 genofile=paste0("largedata/SNP/randomsnp/rsnp", i, ".csv"), nf=3)
  res1$sample <- i
  out3 <- rbind(out3, res1)
}
write.table(out3, "cache/rsnp_var_BPH_nf3.csv", sep=",", row.names=FALSE, quote=FALSE)

### collect all data
out0 <- data.frame()
for(i in 0:10){
  res1 <- getvar(gblup_efile=paste0("largedata/snpeff/rsnp", i, "/gy_BPHmax_snpeff_ce.snpe"),
                 genofile=paste0("largedata/SNP/randomsnp/rsnp", i, ".csv"), nf=0)
  res1$sample <- i
  out0 <- rbind(out0, res1)
}
write.table(out0, "cache/rsnp_var_BPH_nf0.csv", sep=",", row.names=FALSE, quote=FALSE)
