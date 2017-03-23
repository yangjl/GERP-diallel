### Jinliang Yang
### August 19th, 2015

library("plyr")
library("ggplot2")
library("data.table")
#load data in the local machine for plotting

hmp3 <- fread("largedata/del_allele_frq_chr1.txt", header=T)
hmp3

hmp3[, GERP1 := round(RS, 1)] 
hmp3[, GERP2 := round(RS, 2)]
hmp3[, GERP3 := round(RS, 3)] 


snptab1 <- hmp3[, .(meandaf = mean(daf)), by= GERP1] 
snptab2 <- hmp3[, .(meandaf = mean(daf)), by= GERP2] 
snptab3 <- hmp3[, .(meandaf = mean(daf)), by= GERP3]

write.table(snptab1, "cache/daf_gerp1_chr1.csv", sep=",", row.names=FALSE, quote=FALSE)
write.table(snptab2, "cache/daf_gerp2_chr1.csv", sep=",", row.names=FALSE, quote=FALSE)
write.table(snptab3, "cache/daf_gerp3_chr1.csv", sep=",", row.names=FALSE, quote=FALSE)

snptab <- read.csv("cache/daf_gerp2_chr1.csv")
plotReg <- function(snptab, ...){
  #lm1 <- lm(y~x)
  snptab <- snptab[order(snptab$GERP2),]
  plx <- predict(loess(snptab$meandaf ~ snptab$GERP2), se=T)
  
  x <- snptab$GERP2
  y <- snptab$meandaf
  #p_conf1 <- predict(lm1, interval="confidence")
  #p_pred1 <- predict(lm1,interval="prediction")
  plot(x, y, ...)
  #lines(predict(lm1), col='red', lwd=2)
  #matlines(x,p_conf1[,c("lwr","upr")], col="grey", lty=2, lwd=2, type="b", pch="+")
  
  lines(snptab$GERP2, plx$fit, col="cornflowerblue", lwd=2)
  lines(snptab$GERP2, plx$fit - qt(0.975,plx$df)*plx$se, lty=2, lwd=2, col="black")
  lines(snptab$GERP2, plx$fit + qt(0.975,plx$df)*plx$se, lty=2, lwd=2, col="black")
}

plotReg(snptab,
        pch=16, col="antiquewhite3", xlab="GERP Score", ylab="Deleterious Allele Frequency", 
        main="")


#############################################
pdf("graphs/Fig1_c_new.pdf", width=5, height=5)
plotReg(snptab,
        pch=16, col="antiquewhite3", xlab="GERP Score", ylab="Deleterious Allele Frequency", 
        main="")
dev.off()



qplot(x=as.numeric(as.character(snptab$GERP2)), y=snptab$meandaf, geom='smooth', span =0.75,
      xlab="GERP Score", ylab="Derived Frequency")

qplot(x=hmp3$GERP, y=hmp3$derived.freq, geom='smooth', span =0.75,
      xlab="GERP Score", ylab="Derived Frequency")

qplot(x=hmp3$GERP, y=hmp3$derived.freq, geom='smooth', span =0.9)

##############################################################
