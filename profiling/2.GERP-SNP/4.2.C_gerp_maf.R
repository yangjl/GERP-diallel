### Jinliang Yang
### 11/23/2014
### purpose: find the relationship of RS with SNP frq

library("plyr")
library("ggplot2")
#load data in the local machine for plotting
ob <- load("largedata/lcache/snpnzRS.RData")
snpnz$MAF2 <- as.factor(round(snpnz$MAF, 3))


snptab <- ddply(subset(snpnz, RS >0 ), .(MAF2), summarise,
                rsmean= mean(RS),
                rssd = sd(RS))

plotReg <- function(x, y, ...){
  lm1 <- lm(y~x)
  p_conf1 <- predict(lm1,interval="confidence")
  #p_pred1 <- predict(lm1,interval="prediction")
  plot(x, y, ...)
  abline(lm1, lwd=2, col="red") ## fit
  matlines(x,p_conf1[,c("lwr","upr")], col="grey", lty=2, lwd=2, type="b", pch="+")
}

##############################################################
pdf("graphs/Fig2_ab.pdf", width=10, height=5)

par(mfrow=c(1,2))
hist(snpnz$RS, breaks=50, main="Distribution of GERP (N=1.2M)", xlab="GERP score", 
     col="antiquewhite3", cex.lab=1)
abline(v=0, lty=2, col="grey", lwd=3)
nrow(subset(snpnz, RS>2))
nrow(subset(snpnz, RS<0))
plotReg(x=as.numeric(as.character(snptab$MAF2)), y=snptab$rsmean,
        pch=16, col="cornflowerblue", xlab="MAF", ylab="Avg. GERP", main="GERP vs. MAF")

dev.off()

#############################################
pdf("graphs/Fig1_c_new.pdf", width=5, height=5)
plotReg(x=as.numeric(as.character(snptab$MAF2)), y=snptab$rsmean,
        pch=16, col="cornflowerblue", xlab="MAF", ylab="Mean GERP Score", main="")
dev.off()

#matlines(d$x,p_pred1[,c("lwr","upr")],col=2,lty=2,type="b",pch=1)
#matlines(nd$x,p_conf2[,c("lwr","upr")],col=4,lty=1,type="b",pch="+")
#matlines(nd$x,p_pred2[,c("lwr","upr")],col=4,lty=2,type="b",pch=1)




###############################
hetplot=function(a,b){  
  bob=data.frame(a,b)
  #namex=deparse(substitute(a));
  #namey=deparse(substitute(b));
  the_plot=  ggplot(bob)+geom_point(aes(y=b,x=a), cex=2, shape = 3) + 
    theme_bw() +
    geom_smooth(aes(y=b,x=a), method="lm") +   
    geom_text(aes(-Inf,Inf,label=paste("r2=",round(summary(lm(b~a))$r.squared,2))),vjust=2,hjust=-0.5) +
    
    theme(plot.title = element_text(lineheight=3, face="bold", color="black", size=29)) +
    ggtitle("Simulated Maize Lines") +
    labs(x="Number of SNPs", y="Log10 (L1 / L2)");
  print(the_plot)
}
hetplot(a=as.numeric(as.character(snptab$MAF2)), b=snptab$rsmean)



### GERP data
ob <- load("largedata/lcache/4.1.A_gerpdis.RData")
chrall <- rbind(chr1, chr2, chr3, chr4, chr5, chr6, chr7, chr8, chr9, chr10)

chrall <- data.table(chrall)

### GERP>0 sites
chrall0 <- subset(chrall, RS > 0) # 86006888        4
chrall0[, snpid:=paste0(chr, "_", pos)]

### merging by setkey
setkey(snp11m, snpid)
setkey(chrall, snpid)

snpgs <- snp11m[chrall, roll=T]
system.time(new.T <- zea.grin.T[amesWyears.T, roll =T])


