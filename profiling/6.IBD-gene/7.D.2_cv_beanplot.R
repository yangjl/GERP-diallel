### Jinliang Yang
### beanplot

library("beanplot")
library(ggplot2)
library(plyr)
#http://www.jstatsoft.org/v28/c01/paper
mybean <- function(res0, mymode="a2", ...){
  res0$mode <- as.character(res0$mode)
  res1 <- subset(res0, mode == mymode)
  res1$trait <- toupper(res1$trait)
  #print(nrow(res1))
  par(lend = 1, mai = c(0.8, 0.8, 0.5, 0.5))
  res1$type <- factor(res1$type, levels = c("real", "random"))
  res1$trait <- factor(res1$trait, levels = toupper(c("dts", "dtp", "tw", "asi","pht", "eht", "gy")))
  beanplot(m ~ type + trait, data = res1, kernel="cosine", ll = 0.04, cex=1.5, side = "both", bw=0.02,
           border = NA, col = list( c("#008080", "#008080"), c("grey", "grey")), ...)
  #legend("bottomleft", fill = c("black", "grey"),
  #       legend = c("Group 2", "Group 1"))
  #return(res0)
  
}

############################################################

res01 <- read.csv("cache/gene_perse_2016.csv")
res1 <- ddply(res01, .(trait, mode, sp, type), summarise,
              r = mean(r),
              m = median(r))
mean(subset(res1, mode=="a2")$r) #0.8352665
mean(subset(res1, mode=="d2")$r) #0.7495015
mean(subset(res1, mode=="h2")$r) #0.7741429

mean(subset(res1, mode=="a2" & trait == "gy" & type=="real")$r)
mean(subset(res1, mode=="a2" & trait == "gy" & type=="random")$r)

mean(subset(res1, mode=="d2" & trait == "gy" & type=="real")$r)
mean(subset(res1, mode=="d2" & trait == "gy" & type=="random")$r)
mean(subset(res1, mode=="h2" & trait == "gy" & type=="real")$r)
mean(subset(res1, mode=="h2" & trait == "gy" & type=="random")$r)

mean(subset(res1, mode=="h2" & trait == "tw" & type=="real")$r) - mean(subset(res1, mode=="h2" & trait == "tw" & type=="random")$r)
mean(subset(res1, mode=="h2" & trait == "asi" & type=="real")$r) - mean(subset(res1, mode=="h2" & trait == "asi" & type=="random")$r)
mean(subset(res1, mode=="h2" & trait == "pht" & type=="real")$r) - mean(subset(res1, mode=="h2" & trait == "pht" & type=="random")$r)

mean(subset(res1, mode=="h2" & trait == "gy" & type=="real")$r) - mean(subset(res1, mode=="h2" & trait == "gy" & type=="random")$r)





res02 <- read.csv("cache/gene_bph_2016.csv")
res2 <- ddply(res02, .(type, trait, mode, sp), summarise,
              r = mean(r),
              m = median(r))

res02 <- read.csv("cache/g0_k_bph_2016.csv")
res2 <- ddply(res02, .(trait, mode, sp, type), summarise,
              r = mean(r),
              m = median(r))
mean(subset(res2, mode=="a2")$r)
mean(subset(res2, mode=="d2")$r)
mean(subset(res2, mode=="h2")$r)


mean(subset(res2, mode=="a2" & trait == "gy" & type=="real")$r)
mean(subset(res2, mode=="a2" & trait == "gy" & type=="random")$r)

mean(subset(res2, mode=="d2" & trait == "gy" & type=="real")$r)
mean(subset(res2, mode=="d2" & trait == "gy" & type=="random")$r)
mean(subset(res2, mode=="h2" & trait == "gy" & type=="real")$r)
mean(subset(res2, mode=="h2" & trait == "gy" & type=="random")$r)


############

#######
pdf("graphs/SFig4_6plots.pdf", height=10, width=15)
par(mfrow=c(2,3))
mybean(res1, mymode = "a2", ylim=c(0, 1), main="Additive", ylab="Cross-validation Accuracy")
mybean(res1, mymode = "d2", ylim=c(0, 1), main="Dominance", ylab="Cross-validation Accuracy")
mybean(res1, mymode = "h2", ylim=c(0, 1), main="Incomplete Dominance", ylab="Cross-validation Accuracy")

mybean(res2, mymode = "a2", ylim=c(0, 1), main="Additive", ylab="Cross-validation Accuracy")
mybean(res2, mymode = "d2", ylim=c(0, 1), main="Dominance", ylab="Cross-validation Accuracy")
mybean(res2, mymode = "h2", ylim=c(0, 1), main="Incomplete Dominance", ylab="Cross-validation Accuracy")
dev.off()




###################################################################################
myt <- c( "dtp", "dts", "tw", "asi", "pht", "eht",  "gy")
runttest <- function(res0, mymode="h2", mytrait=myt[7]){
  
  res0$R2 <- res0$r^2
  real <- subset(res0, type == "real" & trait== mytrait & mode == mymode)
  rand <- subset(res0, type == "random" & trait == mytrait & mode == mymode)
  
  message(sprintf("###>>> real [ %s ] and random [ %s ]", nrow(real), nrow(rand)))
  test <- t.test(real$r, rand$r, alternative ="greater")
  print(sprintf("### Trait [ %s ]", mytrait))
  print(test$p.value)
  #return(rbind(out1, out2[, -1]))
}


for(i in 1:7){
  runttest(res0=res1, mymode="h2", mytrait=myt[i])
}

for(i in 1:7){
  runttest(res0=res2, mymode="a2", mytrait=myt[i])
}








