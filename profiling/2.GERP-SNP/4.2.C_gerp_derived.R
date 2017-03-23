### Jinliang Yang
### August 19th, 2015

library("plyr")
library("ggplot2")
#load data in the local machine for plotting

hmp3 <- read.table("data/smaller.txt", header=T)

snptab <- ddply(hmp3, .(GERP), summarise,
                frqmean= mean(derived.freq),
                frqsd = sd(derived.freq))

plotReg <- function(x, y, ...){
  #lm1 <- lm(y~x)
  lo <- loess(y~x)
  #p_conf1 <- predict(lm1, interval="confidence")
  #p_pred1 <- predict(lm1,interval="prediction")
  plot(x, y, ...)
  lines(predict(lo), col='red', lwd=2)
  #matlines(x,p_conf1[,c("lwr","upr")], col="grey", lty=2, lwd=2, type="b", pch="+")
}

plotReg(x=as.numeric(as.character(snptab$GERP)), y=snptab$frqmean,
        pch=16, col="cornflowerblue", xlab="MAF", ylab="Avg. GERP", main="GERP vs. MAF")

qplot(x=as.numeric(as.character(snptab$GERP)), y=snptab$frqmean, geom='smooth', span =0.75,
      xlab="GERP Score", ylab="Derived Frequency")

qplot(x=hmp3$GERP, y=hmp3$derived.freq, geom='smooth', span =0.75,
      xlab="GERP Score", ylab="Derived Frequency")

qplot(x=hmp3$GERP, y=hmp3$derived.freq, geom='smooth', span =0.9)

##############################################################
pdf("graphs/Fig2_ab.pdf", width=10, height=5)

par(mfrow=c(1,2))
hist(snpnz$RS, breaks=50, main="Distribution of GERP (N=1.2M)", xlab="GERP score", 
     col="antiquewhite3", cex.lab=1)
abline(v=0, lty=2, col="grey", lwd=3)
nrow(subset(snpnz, RS>2))
nrow(subset(snpnz, RS<0))

dev.off()