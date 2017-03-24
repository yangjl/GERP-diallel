### Jinliang Yang
### beanplot

#http://www.jstatsoft.org/v28/c01/paper

runttest <- function(resfile="cache/gerpall_perse.csv"){
  resin <- read.csv(resfile)
  print(table(resin$trait))
  print(head(resin))
  
  resin$trait <- as.character(resin$trait)
  myt <- unique(resin$trait)
  resin$cs <- gsub("cs.*", "cs", resin$cs)

  #### >>>>
  res <- data.frame()
  for(ti in myt){
    test <- t.test(subset(resin, type=="real" & trait== ti)$r, 
                   subset(resin, type=="cs" & trait == ti)$r, alternative ="greater")
    tem <- data.frame(trait= ti, pval=test$p.value,
                      r_real=mean(subset(resin, type=="real" & trait== ti )$r),
                      r_cs=mean(subset(resin, type=="cs" & trait == ti )$r))
    res <- rbind(res, tem)
  }
  res$file <- resfile
  return(res)
}

#####
res1 <- runttest(resfile="largedata/newGERPv2/res_a2_perse_42000.csv")
res2 <- runttest(resfile="largedata/newGERPv2/res_d2_perse_42000.csv")

res2 <- runttest(resfile="largedata/newGERPv2/res_d2_bph_42000.csv")

res3 <- runttest(resfile="cache/gerpall_pBPHmax.csv")

res4 <- runttest(resfile="cache/gerpall_BPHmin.csv")
res5 <- runttest(resfile="cache/gerpall_pBPHmin.csv")
res6 <- runttest(resfile="cache/gerpall_MPH.csv")
res7 <- runttest(resfile="cache/gerpall_pMPH.csv")



res7
pval <- rbind(res1, res2, res3, res4, res5, res6, res7)
pval <- subset(pval, mode %in% c("a2", "d2"))
write.table(pval, "cache/pval_gerpall.csv", sep=",", row.names=FALSE, quote=FALSE)



##################################################################################################################









