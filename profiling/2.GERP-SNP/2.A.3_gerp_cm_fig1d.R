### Jinliang Yang
### GERP in cM/Mb
### Feb 3th, 2016


library("data.table")

v23 <- fread("largedata/SNP/gerpsnp_v3.bed", data.table=FALSE)
v23$agpv3 <- paste(v23$V1, v23$V3, sep="_")

### v2 
gerp <- fread("cache/gerpsnp_506898_gp.csv", data.table=FALSE)
gerp <- merge(gerp, v23[, c("V4", "agpv3")], by.x="marker", by.y="V4")

## read in v3 RS
gerp3 <- data.frame()
for(i in 1:10){
    rs3 <- fread(paste0("largedata/GERPv3/AGPv3_gerp_chr", i, ".csv"), data.table=FALSE)
    rs3$marker <- paste(rs3$chrom, rs3$pos, sep="_")
    sub <- subset(rs3, marker %in% gerp$agpv3)
    gerp3 <- rbind(gerp3, sub)
}
dim(gerp3)
# 410183      5

names(gerp)[3] <- "RSv2"
gerp23 <- merge(gerp, gerp3[, c("marker", "RS")], by.x="agpv3", by.y="marker")

write.table(gerp23, "largedata/gerp23.csv", sep=",", row.names=FALSE, quote=FALSE)


####################################
library(plyr)
BINSIZE = 1000000
gerp23$bin <- paste(gerp23$chr, round(gerp23$pos/BINSIZE, 0), sep="_")
res <- ddply(gerp23, .(bin), summarise,
             mgerp = mean(RS),
             gen = BINSIZE*(max(genetic) - min(genetic))/(max(pos) - min(pos)) )
range(res$gen)

res <- subset(res, !is.na(gen))
t.test(subset(res, gen < 0.5)$mgerp, subset(res, gen > 0.5)$mgerp)
write.table(res, "cache/mgerp_cm.csv", sep=",", row.names = FALSE, quote = FALSE)



##########################
plotmgerp <- function(mgerp, getpdf, outfile, ...){
    
    #mgerp <- read.csv("cache/mgerp_cm.csv")
    
    mgerp <- mgerp[order(mgerp$gen),]
    plx <- predict(loess(mgerp$mgerp ~ mgerp$gen), se=T)
    x <- mgerp$gen
    y <- mgerp$mgerp
    
    plot(x, y, ...)
    lines(mgerp$gen, plx$fit, col="cornflowerblue", lwd=2)
    lines(mgerp$gen, plx$fit - qt(0.975,plx$df)*plx$se, lty=2, lwd=2, col="black")
    lines(mgerp$gen, plx$fit + qt(0.975,plx$df)*plx$se, lty=2, lwd=2, col="black")
    
    if(getpdf == TRUE){
        pdf(outfile, width=wt, height=ht)
        par(mar=c(5,5,4,2))
        plot(x, y, ...)
        lines(mgerp$gen, plx$fit, col="cornflowerblue", lwd=2)
        lines(mgerp$gen, plx$fit - qt(0.975,plx$df)*plx$se, lty=2, lwd=2, col="black")
        lines(mgerp$gen, plx$fit + qt(0.975,plx$df)*plx$se, lty=2, lwd=2, col="black")
        
        dev.off()
    }
}

mgerp <- read.csv("cache/mgerp_cm.csv")
par(font= 2, font.lab= 2, font.axis= 2)
fs <- 1.6 # times bigger than default
ht= 6; wt= 6 #figure height and weight
getpdf <- TRUE # get figures in seperated pdf [TRUE] or not [FALSE]
plotmgerp(mgerp, getpdf=TRUE, outfile="graphs/Fig1d.pdf",
        pch=16, col="antiquewhite3", xlab="Recombination Rate (cM/Mb)", ylab="GERP Score", 
        main="", cex.axis=fs, cex.lab=fs)

cor.test(mgerp$mgerp, mgerp$gen)

# get figures in seperated pdf [TRUE] or not [FALSE]
plotmgerp(mgerp=subset(mgerp, gen < 0.5), getpdf=FALSE, outfile="graphs/Fig1d.pdf",
          pch=16, col="antiquewhite3", xlab="Recombination Rate (cM/Mb)", ylab="GERP Score", 
          main="", cex.axis=fs, cex.lab=fs)



res <- read.csv("cache/mgerp_cm.csv")
cutoff <- quantile(res$gen, c(0.25, 0.5, 0.75))
res$sq <- 0
res[res$gen < cutoff[1], ]$sq <- 1
res[res$gen >= cutoff[1] & res$gen < cutoff[2], ]$sq <- 2
res[res$gen >= cutoff[2] & res$gen < cutoff[3], ]$sq <- 3
res[res$gen >= cutoff[3], ]$sq <- 4

library("beanplot")
pdf("graphs/Fig1d.pdf", width=5, height=5)
#par(mar = c(7, 4, 4, 2) + 0.1)
beanplot(mgerp ~ sq, data = res, kernel="cosine", ll = 0.04, cex=1.5, side = "no", cut=10, ylim=c(0.5, 5),
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
