### Jinliang Yang
### Sept 19th, 2015

##>>>>>
library("beanplot")
library("dplyr")

avar5 <- read.csv("cache/rsnp_var_nf1.csv")

avar5$frq1 <- round(avar5$frq, 1)
res <- ddply(avar5, .(frq1, sample), summarise,
             totvar = sum(totvar),
             numsnp = sum(V1))
res$nvar <- res$totvar/res$numsnp
write.table(res, "cache/persnp_var_perse.csv", sep=",", row.names=FALSE, quote=FALSE)
### Note: 0.15 is missing from the analysis for the random samples
### therefore, concatenate 0.1 and 0.15 into a single catergory
myd <- subset(res, sample != 0)
myd0 <- subset(res, sample == 0)

#### for trait perse
beanplot(nvar ~ frq1, data = myd, ll = 0.04, cex=1.5, border = NA,
         what=c(0, 1, 0, 1), col=c(c("grey", "black")), cex=1.3, 
         xlab="Allele Frequency", ylab="Variance Explained" )
for(i in 1:nrow(myd0)){
  lines(x=c(i-0.3, i+0.3), y=c(myd0$nvar[i], myd0$nvar[i]), lwd=2, col="red")
}



##############
beanplot(nvar ~ frq, data = myd, ll = 0.04, cex=1.5, border = NA,
         what=c(0, 1, 0, 1), col=c(c("grey", "black")) )
for(i in 1:nrow(myd0)){
  lines(x=c(i-0.3, i+0.3), y=c(myd0$nvar[i], myd0$nvar[i]), lwd=2, col="red")
}





avar1 <- read.csv("cache/rsnp_var_nf1.csv")

### Note: 0.15 is missing from the analysis for the random samples
### therefore, concatenate 0.1 and 0.15 into a single catergory
myd <- subset(avar1, sample != 0)
#myd[myd$frq==0.10, ]$frq <- 0.15

myd0 <- subset(avar1, sample == 0)

#### for trait perse
beanplot(nvar ~ frq, data = myd, ll = 0.04, cex=1.5, border = NA,
         what=c(0, 1, 0, 1), col=c(c("grey", "black")) )
for(i in 1:nrow(myd0)){
  lines(x=c(i-0.3, i+0.3), y=c(myd0$nvar[i], myd0$nvar[i]), lwd=2, col="red")
}



avar5 <- read.csv("cache/rsnp_var_BPH_nf5.csv")

### Note: 0.15 is missing from the analysis for the random samples
### therefore, concatenate 0.1 and 0.15 into a single catergory
myd <- subset(avar5, sample != 0)
#myd[myd$frq==0.10, ]$frq <- 0.15

myd0 <- subset(avar5, sample == 0)
#myd0[myd0$frq==0.15, ]$totvar <- sum(myd0[myd0$frq<=0.15, ]$totvar)
#myd0 <- myd0[-1,]

#### for trait perse
beanplot(nvar ~ frq, data = myd, ll = 0.04, cex=1.5, border = NA,
         what=c(0, 1, 0, 1), col=c(c("grey", "black")) )
for(i in 1:nrow(myd0)){
  lines(x=c(i-0.3, i+0.3), y=c(myd0$nvar[i], myd0$nvar[i]), lwd=2, col="red")
}





avar1 <- read.csv("cache/rsnp_var_BPH_nf5.csv")

### Note: 0.15 is missing from the analysis for the random samples
### therefore, concatenate 0.1 and 0.15 into a single catergory
myd <- subset(avar1, sample != 0)
#myd[myd$frq==0.10, ]$frq <- 0.15

myd0 <- subset(avar1, sample == 0)

#### for trait perse
beanplot(nvar ~ frq, data = myd, ll = 0.04, cex=1.5, border = NA,  what=c(0, 1, 0, 1), col=c(c("grey", "black")) )
for(i in 1:nrow(myd0)){
  lines(x=c(i-0.3, i+0.3), y=c(myd0$nvar[i], myd0$nvar[i]), lwd=2, col="red")
}

