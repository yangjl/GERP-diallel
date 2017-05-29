### Jinliang Yang
### beanplot

#http://www.jstatsoft.org/v28/c01/paper




GetBeanPlot <- function(res0=res0, m="a2", factors =c("real", "random"), ...){
  res1 <- subset(res0, mode %in% m)
  par(lend = 1, mai = c(0.8, 0.9, 0.5, 0.5))
  res1$type <- factor(res1$type, levels = factors)
  res1$trait <- toupper(res1$trait)
  res1$trait <- factor(res1$trait, levels = toupper(c("asi", "dtp", "dts", "eht", "gy", "pht", "tw")))
  beanplot(r ~ type + trait, data = res1, ll = 0.04, cex=1.5, ylim=c(0,1.2),
           ylab = "Accuracy (r)", side = "both",
           border = NA, col = list(c("blue", "red"), c("grey", "black")),
           ...)
  #legend("bottomleft", fill = c("black", "grey"),
  #       legend = c("Group 2", "Group 1"))
}


#### trait per se:
library("beanplot")
res0 <- read.csv("cache/cv_results0.csv")
par(mfrow=c(2,2))
GetBeanPlot(res0=res0, m="a2", factors =c("real", "random"), main="Trait per se with additive model")
GetBeanPlot(res0=res0, m="d2", factors =c("real", "random"), main="Trait per se with dominant model")
    

#### pHPH
phph <- read.csv("cache/genic_cv_all_pHPH.csv")
phph$type <- "A"
phph[phph$cs %in% "cs0",]$type <- "real"
phph[phph$cs != "cs0", ]$type <- "random"


GetBeanPlot(res0=phph, m="a2", factors =c("real", "random"), main="pHPH with additive model")
GetBeanPlot(res0=phph, m="d2", factors =c("real", "random"), main="pHPH with dominant model")


######
runttest <- function(res0=res0, mymode="d2", mytrait=myt[1]){
  res <- data.frame()
  for(i in 1:length(mytrait)){
    test <- t.test(subset(res0, type=="real" & trait== mytrait[i] & mode == mymode)$r, 
                   subset(res0, type=="random" & trait == mytrait[i] & mode == mymode)$r, alternative ="greater")
    tem <- data.frame(trait=mytrait[i], pval=test$p.value)
    res <- rbind(res, tem)
  }
  return(res)
}

#####
myt <- c("asi", "dtp", "dts", "eht", "gy", "pht", "tw")
runttest(res0=res0, mymode="a2", mytrait=myt)
runttest(res0=res0, mymode="d2", mytrait=myt)
runttest(res0=phph, mymode="a2", mytrait=myt)
runttest(res0=phph, mymode="d2", mytrait=myt)














