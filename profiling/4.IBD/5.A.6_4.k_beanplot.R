### Jinliang Yang
### beanplot

library(beanplot)
library(ggplot2)
library(plyr)
#http://www.jstatsoft.org/v28/c01/paper
mybean <- function(inputdf, mymode="a2", ...){
  res1 <- ddply(inputdf, .(trait, mode, sp, type), summarise,
                r = mean(r),
                m = median(r))
  res1 <- subset(res1, type != "null")
  
  res1$mode <- as.character(res1$mode)
  res1 <- subset(res1, mode == mymode)
  res1$trait <- toupper(res1$trait)
  #print(nrow(res1))
  par(lend = 1, mai = c(0.8, 0.8, 0.5, 0.5))
  res1$type <- factor(res1$type, levels = c("real", "cs"))
  res1$trait <- factor(res1$trait, levels = toupper(c("dts", "dtp", "tw", "asi","pht", "eht", "gy")))
  beanplot(m ~ type + trait, data = res1, kernel="cosine", ll = 0.04, cex=1.5, side = "both", bw=0.02,
           border = NA, col = list(c("#d41243", "#d41243"), c("grey", "grey")), ...)
  #legend("bottomleft", fill = c("black", "grey"),
  #       legend = c("Group 2", "Group 1"))
  #return(res0)
  
  out <- data.frame()
  myt <- toupper(c( "dtp", "dts", "tw", "asi", "pht", "eht",  "gy"))
  for(myti in myt){
    real <- subset(res1, type == "real" & trait == myti)
    rand <- subset(res1, type == "cs" & trait == myti)
    
    message(sprintf("###>>> real [ %s ] and random [ %s ]", nrow(real), nrow(rand)))
    test <- t.test(real$r, rand$r, alternative ="greater")
    
    tem <- data.frame(trait=myti, realmean=mean(real$r), csmean=mean(rand$r), pval=test$p.value)
    out <- rbind(out, tem)
  }

  return(out)
}

############################################################

res01 <- read.csv("largedata/newGERPv2/res_realk_perse_42000.csv")
res02 <- read.csv("largedata/newGERPv2/res_realk_bph_42000.csv")

pdf("graphs/Fig3_perse_BPH_2plots.pdf", height=5, width=10)
par(mfrow=c(1,2))
res1 <- mybean(inputdf=res01, mymode = "h2", ylim=c(0, 1), main="Trait per se", ylab="Cross-validation Accuracy")
res2 <- mybean(inputdf=res02, mymode = "h2", ylim=c(0, 1), main="Heterosis (BPH)", ylab="Cross-validation Accuracy")
dev.off()



############
res3 <- read.csv("largedata/newGERPv2/res_a2_perse_42000.csv")
res4 <- read.csv("largedata/newGERPv2/res_a2_bph_42000.csv")
res5 <- read.csv("largedata/newGERPv2/res_d2_perse_42000.csv")
res6 <- read.csv("largedata/newGERPv2/res_d2_bph_42000.csv")
res7 <- read.csv("largedata/newGERPv2/res_k5_perse_42000.csv")
res8 <- read.csv("largedata/newGERPv2/res_k5_bph_42000.csv")

pdf("graphs/FigS5_BPH_4plots.pdf", height=10, width=10)
par(mfrow=c(2,2))

out3 <- mybean(res3, mymode = "a2", ylim=c(0, 1), main="Trait per se with additive", ylab="Cross-validation Accuracy")
out5 <- mybean(res5, mymode = "d2", ylim=c(0, 1), main="Trait per se with dominance", ylab="Cross-validation Accuracy")
#out7 <- mybean(res7, mymode = "h2", ylim=c(0, 1), main="Incomplete Dominance", ylab="Cross-validation Accuracy")

out4 <- mybean(res4, mymode = "a2", ylim=c(0, 1), main="Heterosis with additive", ylab="Cross-validation Accuracy")
out6 <- mybean(res6, mymode = "d2", ylim=c(0, 1), main="Heterosis with dominance", ylab="Cross-validation Accuracy")
#out8 <- mybean(res8, mymode = "h2", ylim=c(0, 1), main="Incomplete Dominance", ylab="Cross-validation Accuracy")
dev.off()

#######
pdf("graphs/SFig3_perse_3plots.pdf", height=4, width=12)
par(mfrow=c(1,3))

mybean(res1, mymode = "a2", ylim=c(0, 1), main="Additive", ylab="Cross-validation Accuracy")
mybean(res1, mymode = "d2", ylim=c(0, 1), main="Dominance", ylab="Cross-validation Accuracy")
mybean(res1, mymode = "h2", ylim=c(0, 1), main="Incomplete Dominance", ylab="Cross-validation Accuracy")
dev.off()


############
pdf("graphs/Figure5_2plot2.pdf", height=5, width=10)
par(mfrow=c(1,2))
mybean(res1, mymode = "h2", ylim=c(0, 1), main="Trait per se", ylab="Cross-validation Accuracy")
mybean(res2, mymode = "h2", ylim=c(0, 1), main="Heterosis", ylab="Cross-validation Accuracy")
dev.off()


############ For BAPG meeting
my_BAPG_bean <- function(res0, mymode="a2", ...){
  res0$mode <- as.character(res0$mode)
  res1 <- subset(res0, mode == mymode)
  res1$trait <- toupper(res1$trait)
  #print(nrow(res1))
  par(lend = 1, mai = c(0.8, 0.8, 0.5, 0.5))
  res1$type <- factor(res1$type, levels = c("real", "random"))
  res1$trait <- factor(res1$trait, levels = toupper(c("dtp", "pht", "gy")))
  beanplot(m ~ type + trait, data = res1, kernel="cosine", ll = 0.04, cex=1.5, side = "both", bw=0.02,
           border = NA, col = list(c("#d41243", "#d41243"), c("grey", "grey")), ...)
  #legend("bottomleft", fill = c("black", "grey"),
  #       legend = c("Group 2", "Group 1"))
  #return(res0)
  
}
res1 <- subset(res1, trait %in% c("dtp", "pht", "gy"))
res2 <- subset(res2, trait %in% c("dtp", "pht", "gy"))
pdf("graphs/Figure3_BAPG.pdf", height=5, width=10)
par(mfrow=c(1,2))
my_BAPG_bean(res1, mymode = "h2", ylim=c(0, 1), main="Trait per se", ylab="Cross-validation Accuracy")
my_BAPG_bean(res2, mymode = "h2", ylim=c(0, 1), main="Heterosis", ylab="Cross-validation Accuracy")
dev.off()








res1 <- subset(res0, mode == "d2")
par(lend = 1, mai = c(0.8, 0.8, 0.5, 0.5))
res1$type <- factor(res1$type, levels = c("real", "random"))
res1$trait <- factor(res1$trait, levels = c("tw", "dtp", "dts", "pht", "eht", "asi", "gy"))
beanplot(R2 ~ type + trait, data = res1, ll = 0.04, cex=1.5,
         main = "Dominant GERP Score", ylab = "cross-validation accuracy", side = "both",
         border = NA, col = list(c("blue", "red"), c("grey", "black")))
#legend("bottomleft", fill = c("black", "grey"),
#       legend = c("Group 2", "Group 1"))

res1 <- subset(res0, mode == "h2")
par(lend = 1, mai = c(0.8, 0.8, 0.5, 0.5))

out$trait <- factor(out$trait, levels = c("tw", "dtp", "dts", "pht", "eht", "asi", "gy"))

out$trait <- "gy"
out$type <- factor(out$type, levels = c("real", "random"))
beanplot(R2 ~ type + trait, data = out, ll = 0.04, cex=1.5,
         main = "h GERP Score", ylab = "cross-validation accuracy", side = "both",
         border = NA, col = list(c("blue", "red"), c("grey", "black")))









# A basic box with the conditions colored
library(ggplot2, lib="~/bin/Rlib/")

res0 <- read.csv("cache/g2_k_perse.csv")
res0 <- subset(res0, type=="real")
out2 <- ddply(res0, .(sp, trait, mode), summarise,
              r = mean(r),
              m = median(r))
out2$trait <- toupper(out2$trait)
out2$trait <- factor(out2$trait, levels = toupper(c("asi", "dts", "dtp", "tw", "pht", "eht", "gy")))

p2 <- ggplot(data=out2) +
  geom_boxplot(aes(x= mode, y=r, fill= factor(mode, labels=c("A", "D", "k")) ) ) +
  #guides(fill=FALSE) +
  labs(y=NULL, fill="") + theme_bw() +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  xlab("Trait per se") + ylab("Prediction Accuracy (r)") + facet_grid(~ trait) 

pdf("graphs/Figure5_a.pdf", height=5, width=10)
p2
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
  runttest(res0=res2, mymode="d2", mytrait=myt[i])
}

for(i in 1:7){
  runttest(res0=res1, mymode="h2", mytrait=myt[i])
}

for(i in 1:7){
  runttest(res0=res2, mymode="h2", mytrait=myt[i])
}



