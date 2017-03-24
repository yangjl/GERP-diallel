### Jinliang Yang
### test the accuracy differences

###################################################################################
library(plyr, lib="~/bin/Rlib/")
runttest <- function(res0){
  
  #res0$R2 <- res0$r^2
  res0$trait <- as.character(res0$trait)
  res0$type <- as.character(res0$type)
  res0$cs <- as.character(paste0("cs", res0$cs))
  res0$cv <- as.character(paste0("cv", res0$cv))
  res0$sp <- as.character(paste0("sp", res0$sp))
  
  tab <- ddply(res0, .(trait, cs, sp), summarise,
                r = mean(r) )
  
  #tab$trait <- as.character(tab$trait)
  #tab$type <- as.character(tab$type)
  
  #tab <- res0
  tab$type <- "cs"
  tab[tab$cs == "cs0",]$type <- "real"
  tab[tab$cs == "cs999",]$type <- "null"
  
  myt <- unique(tab$trait)
  res <- data.frame()
  for(i in 1:length(myt)){
    real <- subset(tab, type == "real" & trait== myt[i])
    
    rand <- subset(tab, type == "cs" & trait == myt[i])
    nl <- subset(tab, type == "null" & trait == myt[i])
    
    message(sprintf("###>>> real [ %s ], null [ %s ] and random [ %s ]", 
                    nrow(real), nrow(nl), nrow(rand)))
    test1 <- t.test(real$r, rand$r, alternative="greater")
    test2 <- try(t.test(real$r, nl$r, alternative="greater"))
    tem <- data.frame(trait = myt[i], pval1=test1$p.value, pval2=test2$p.value,
                      r_real = mean(subset(res0, type == "real" & trait== myt[i] )$r),
                      r_cs = mean(subset(res0, type == "cs" & trait == myt[i] )$r),
                      r_null = mean(subset(res0, type == "null" & trait == myt[i] )$r))
    res <- rbind(res, tem)
  }
  return(res)
}


res_a <- read.csv("largedata/newGERPv2/res_a2_perse_42000.csv")
t1 <- runttest(res0=res_a)

res_d <- read.csv("largedata/newGERPv2/res_d2_perse_42000.csv")
t2 <- runttest(res0=res_d)

res_k5 <- read.csv("largedata/newGERPv2/res_k5_perse_42000.csv")
t3 <- runttest(res0=res_k5)

res_k <- read.csv("largedata/newGERPv2/res_realk_perse_42000.csv")
t4 <- runttest(res0=res_k)
t4$diff <- t4$r_real - t4$r_cs

######
res4 <- read.csv("largedata/newGERPv2/res_d2_bph_42000.csv")
t4 <- runttest(res0=res4)

res5 <- read.csv("largedata/newGERPv2/res_k5_bph_42000.csv")
t5 <- runttest(res0=res5)


