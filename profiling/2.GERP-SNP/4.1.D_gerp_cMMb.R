### Jinliang Yang
### GERP in cM/Mb
### Feb 3th, 2016


### server and local
gerp <- read.csv("largedata/GERPv2/gerpsnp_506898.csv")
gerp <- gerp[, 1:8]
write.table(gerp, "cache/gerpsnp_506898_gp.csv", sep=",", row.names=FALSE, quote=FALSE)


#####
gerp <- read.csv("cache/gerpsnp_506898_gp.csv")
library(plyr)
BINSIZE = 5000000
gerp$bin <- paste(gerp$chr, round(gerp$pos/BINSIZE, 0), sep="_")
res <- ddply(gerp, .(bin), summarise,
            mgerp = mean(RS),
            gen = 1000000*(max(genetic) - min(genetic))/(max(pos) - min(pos)) )
range(res$gen)
#####
mygerp <- merge(gerp, res, by="bin")
quantile(mygerp$gen)
write.table(mygerp, "cache/gerpsnp_506898_gp_withcm_mb.csv", sep=",", row.names=FALSE, quote=FALSE)

#######
t.test(subset(res, gen < 0.5)$mgerp, subset(res, gen > 0.5)$mgerp)

res$sq <- 0
res[res$gen < 0.5, ]$sq <- 1
res[res$gen >= 0.5, ]$sq <- 2

table(res$sq)

library("beanplot")
pdf("graphs/Fig1d.pdf", width=5, height=5)
beanplot(mgerp ~ sq, data = res, kernel="cosine", ll = 0.04, cex=1.5, side = "no", cut=3,
         border = NA, col= list(c("#cd5b45", "antiquewhite3", "antiquewhite3", "antiquewhite3")), 
         xaxt="n", ylab="mean ERP Score", xlab="Recombination rate in quantile intervals")
axis(side =1, at =1:2, labels =c("< 0.5", " >= 0.5"))
dev.off()



#################
gerp <- read.csv("cache/gerpsnp_506898_gp.csv")


library(plyr)
BINSIZE = 5000000
gerp$bin <- paste(gerp$chr, round(gerp$pos/BINSIZE, 0), sep="_")
res <- ddply(gerp, .(bin), summarise,
             mgerp = mean(RS),
             gen = 5000000*(max(genetic) - min(genetic))/(max(pos) - min(pos)) )
range(res$gen)

res <- subset(res, !is.na(gen))
t.test(subset(res, gen < 0.5)$mgerp, subset(res, gen > 0.5)$mgerp)
write.table(res, "cache/mgerp_cm.csv", sep=",", row.names = FALSE, quote = FALSE)

cutoff <- quantile(res$gen, c(0.25, 0.5, 0.75))
res$sq <- 0
res[res$gen < cutoff[1], ]$sq <- 1
res[res$gen >= cutoff[1] & res$gen < cutoff[2], ]$sq <- 2
res[res$gen >= cutoff[2] & res$gen < cutoff[3], ]$sq <- 3
res[res$gen >= cutoff[3], ]$sq <- 4

pdf("graphs/Fig1d.pdf", width=5, height=5)
#par(mar = c(7, 4, 4, 2) + 0.1)
beanplot(mgerp ~ sq, data = res, kernel="cosine", ll = 0.04, cex=1.5, side = "no", cut=10, ylim=c(0.5, 1.5),
         border = NA, col=list("#cd5b45", "antiquewhite3", "antiquewhite3", "antiquewhite3"), 
         xaxt="n", ylab="GERP Score", xlab="Recombination rate in quantile intervals")
#axis(side =1, at =1:4, las=3, labels =c("[0, 25%)", "[25%, 50%)", "[50%, 75%)", "[75%, 100%)"))
axis(side =1, at =1:4, labels =c("<25%", "25-50%", "50-75%", " >75%"), cex.lab=0.25)

dev.off()

t.test(subset(res, sq == 1)$mgerp, subset(res, sq ==2)$mgerp)
t.test(subset(res, sq == 1)$mgerp, subset(res, sq ==3)$mgerp)
t.test(subset(res, sq == 1)$mgerp, subset(res, sq ==4)$mgerp)

t.test(subset(res, sq == 2)$mgerp, subset(res, sq ==3)$mgerp)
t.test(subset(res, sq == 2)$mgerp, subset(res, sq ==4)$mgerp)

t.test(subset(res, sq == 3)$mgerp, subset(res, sq ==4)$mgerp)
