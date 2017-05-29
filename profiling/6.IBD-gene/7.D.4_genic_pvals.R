### Jinliang Yang
### beanplot

#http://www.jstatsoft.org/v28/c01/paper

runttest2 <- function(resfile="cache/genic_perse.csv"){
  resin <- read.csv(resfile)
  print(table(resin$trait))
  print(head(resin))
  #### >>>>
  res <- data.frame()
  myt <- unique(resin$trait)
  resin$cs <- gsub("cs0", "real", resin$cs)
  resin$cs <- gsub("cs.*", "cs", resin$cs)
  mode <- unique(resin$mode)
  for(ti in myt){
    for(modei in mode){
      test <- t.test(subset(resin, cs=="real" & trait== ti & mode == modei)$r, 
                     subset(resin, cs=="cs" & trait == ti & mode == modei)$r, alternative ="greater")
      tem <- data.frame(trait= ti, pval=test$p.value, mode=modei,
                        r_real=mean(subset(resin, cs=="real" & trait== ti & mode == modei)$r),
                        r_cs=mean(subset(resin, cs=="cs" & trait == ti & mode == modei)$r))
      res <- rbind(res, tem)
    }
    
  }
  res$file <- resfile
  return(res)
}

#####
res1 <- runttest2(resfile="cache/genic_perse.csv")
res2 <- runttest2(resfile="cache/genic_BPHmax.csv")
res3 <- runttest2(resfile="cache/genic_pBPHmax.csv")
res4 <- runttest2(resfile="cache/genic_BPHmin.csv")
res5 <- runttest2(resfile="cache/genic_pBPHmin.csv")
res6 <- runttest2(resfile="cache/genic_MPH.csv")
res7 <- runttest2(resfile="cache/genic_pMPH.csv")

pval <- rbind(res1, res2, res3, res4, res5, res6, res7)
#pval <- subset(pval, mode %in% c("a2", "d2"))
write.table(pval, "cache/pval_genic.csv", sep=",", row.names=FALSE, quote=FALSE)

### format data
pval <- read.csv("cache/pval_genic.csv")
idx <- grep("MPH", pval$file)
pval2 <- pval[-idx, ]

pval2$FDR <- p.adjust(pval2$pval, method="fdr" )
subset(pval2, FDR < 0.05)

pval3 <- pval2
pval3$trait <- toupper(pval3$trait)
pval3$mode <- gsub("a2", "additive", pval3$mode)
pval3$mode <- gsub("d2", "dominant", pval3$mode)
pval3$file <- gsub(".*_", "", pval3$file)
pval3$file <- gsub("\\.csv", "", pval3$file)
names(pval3) <- c("phenotype", "p-value", "mode", "r_real", "r_cs", "trait", "FDR")

pval4 <- subset(pval3, trait %in% c("perse", "BPHmax", "BPHmin"))
write.table(pval4, "manuscript/Figure_Table/Table_S4_genicsnps_FDR.csv")



