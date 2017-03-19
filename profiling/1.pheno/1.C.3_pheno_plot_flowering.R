## Jinliang Yang
## Oct. 13th, 2014
## phenotypic data of offpvp diallel

library(magrittr)
#setwd("~/Documents/Github/pvpDiallel/")
trait <- read.csv("data/trait_matrix.csv")



hybplot <- function(trait, myt="GY"){
  
  dts <- subset(trait, trait == myt)
  
  par(mfrow=c(1,1))
  dts <- dts[order(dts$valHyb, decreasing=T),]
  
  xrange <- range(c(dts$valP1, dts$valP2, dts$valHyb))
  plot(dts$valP1, 1:66, xlim= xrange, main="hybrids (red) + inbreds (black) ", 
       xlab= myt, ylab="Hybrids")
  abline(h=1:66, col="grey" )
  points(dts$valP2, 1:66)
  points(dts$valHyb, 1:66, col="red")
  
}

hybplot(trait, myt="PHT")

pht <- subset(trait, trait == "PHT")
hist(pht$valHyb)
shapiro.test(pht$valHyb)
res <- c()
for(i in 1:nrow(pht)){
  tem <- max(pht$valP1[i], pht$valP2[i])
  res <- c(res, tem)
}
cor(res, pht$BPHmax)


allinbred <- c(dts$valP1, dts$valP2)
inbred <- allinbred[!duplicated(allinbred)]

allhy <- dts$valHyb

par(mfrow=c(2,1))
hist(allhy, xlim=c(50, 80), main="Days to silking, Hybrids (N=66)", xlab="Hybrids")
hist(inbred, xlim=c(50, 80), main="Days to silking, Inbreds (N=12)", xlab="Inbreds")

##################
dtp <- subset(trait, trait == "DTP")
allinbred <- c(dtp$valP1, dtp$valP2)
inbred <- allinbred[!duplicated(allinbred)]

allhy <- dtp$valHyb

par(mfrow=c(2,1))
hist(allhy, xlim=c(50, 80), main="Days to pollen, Hybrids (N=66)", xlab="Hybrids")
hist(inbred, xlim=c(50, 80), main="Days to pollen, Inbreds (N=12)", xlab="Inbreds")

##################
dtp <- subset(trait, trait == "GY")
allinbred <- c(dtp$valP1, dtp$valP2)
inbred <- allinbred[!duplicated(allinbred)]

allhy <- dtp$valHyb[!duplicated(dtp$valHyb)]

par(mfrow=c(2,1))
hist(allhy, xlim=c(30, 205), main="GY, Hybrids (N=66)", xlab="Hybrids")
hist(inbred, xlim=c(30, 205), main="GY, Inbreds (N=12)", xlab="Inbreds")






