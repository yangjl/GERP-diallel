## Jinliang Yang
## Oct. 13th, 2014
## phenotypic data of offpvp diallel

#setwd("~/Documents/Github/pvpDiallel/")
H2 <- read.csv("data/DIalleleHeritability.csv")

h <- read.csv("cache/loh_pBPHmax_median.csv")

H2 <- merge(H2, h, by.x="Traits", by.y="trait")
H2 <- H2[order(H2$h),]


pdf("graphs/Fig1a_v2.pdf", width=5, height=5)
barplot(H2[,4], ylim=c(0, 1), col="antiquewhite3", names.arg = H2$Traits, ylab="Heritability", 
        cex.axis=1.3, cex.names=1.3, cex.lab=1.3)
box()
dev.off()


