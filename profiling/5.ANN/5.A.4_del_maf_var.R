### Jinliang Yang
### August 19th, 2015

library("plyr")
library("ggplot2")
library("data.table")
#load data in the local machine for plotting

hmp3 <- fread("largedata/del_allele_frq_chr1.txt", header=T)
dim(hmp3)

### get GERP > 0 SNPs
hmp3 <- hmp3[RS > 0]

hmp3[, dafbin1 := round(daf, 1)] 
hmp3[, dafbin2 := round(daf, 2)]
hmp3[, dafbin3 := round(daf, 3)] 


snptab1 <- hmp3[, .(vargerp = var(RS)), by= dafbin1] 
snptab2 <- hmp3[, .(vargerp = var(RS)), by= dafbin2] 
snptab3 <- hmp3[, .(vargerp = var(RS)), by= dafbin3]

write.table(snptab1, "cache/vargerp_daf_bin1.csv", sep=",", row.names=FALSE, quote=FALSE)
write.table(snptab2, "cache/vargerp_daf_bin2.csv", sep=",", row.names=FALSE, quote=FALSE)
write.table(snptab3, "cache/vargerp_daf_bin3.csv", sep=",", row.names=FALSE, quote=FALSE)


plotReg <- function(snptab, fitline=TRUE, ...){
  #lm1 <- lm(y~x)
  snptab <- snptab[order(snptab[,1]),]
  x <- snptab[,1]
  y <- snptab[,2]
  
  plx <- predict(loess(y ~ x), se=T)
  
  #p_conf1 <- predict(lm1, interval="confidence")
  #p_pred1 <- predict(lm1,interval="prediction")
  plot(x, y, ...)
  #lines(predict(lm1), col='red', lwd=2)
  #matlines(x,p_conf1[,c("lwr","upr")], col="grey", lty=2, lwd=2, type="b", pch="+")
  if(fitline){
      lines(x, plx$fit, col="cornflowerblue", lwd=2)
      lines(x, plx$fit - qt(0.975, plx$df)*plx$se, lty=2, lwd=2, col="black")
      lines(x, plx$fit + qt(0.975, plx$df)*plx$se, lty=2, lwd=2, col="black")
  }
 
}

snptab <- read.csv("cache/vargerp_daf_bin3.csv")
#snptab <- subset(snptab, !is.na(vardaf) & snptab[,1] > 0 & snptab[,1] < 5)
#snptab <- subset(snptab, GERP3 > 0)
plotReg(snptab, fitline=TRUE,
        pch=16, col="antiquewhite3", xlab="Deleterious Allele frequency", ylab="Variance of GERP Score", 
        main="")


#############################################
pdf("graphs/SFigure_var_maf.pdf", width=5, height=5)
plotReg(snptab,
        pch=16, col="antiquewhite3", xlab="GERP Score", ylab="Deleterious Allele Frequency", 
        main="")
dev.off()



