### Jinliang Yang
### chr plot of the avg RS in bins


#######################
gerPlot <- function(binsize=1000000, tab=tab1m, ...){
  
  source("~/Documents/Github/zmSNPtools/Rcodes/rescale.R")
    
  tab <- tab[order(tab$chr, tab$pos), ]
  #### read chr length
  cl <- read.table("data/ZmB73_RefGen_v2.length", header=FALSE)
  names(cl) <- c("chrom", "BP")
  plot(c(0, max(cl$BP)/binsize), c(10,110), type="n", 
       xlab="Physical position (kb)", ylab="", yaxt="n", bty="n", ...)

  axis(side=2, tick =FALSE, las=1, at=c(105, 95, 85, 75, 65, 55, 45, 35, 25, 15), 
       labels=paste("Chr", 1:10, sep=""))
  #### chromosome
  for (i in 1:10){
    lines(c(0, cl[i,]$BP/binsize), c(105-10*(i-1), 105-10*(i-1)), lwd=2, col="grey")
    #lines (c(centromere[i,]$Start,centromere[i,]$End),
    #       c(105-10*i, 105-10*i),lwd=5, col="tomato") 
  }
  ### core plot
  cols <- rep(c("slateblue", "cyan4"), 5)
  for (chri in 1:10){
    mytab <- subset(tab, chr == chri)
    mytab$binrs <- rescale(mytab$binrs, c(-5, 5))
    points(mytab$pos, 105 - 10*(chri-1) + mytab$binrs, pch=19, cex=0.6, col=cols[chri])
  } 
}
###
ob <- load(file="cache/4.1.B_gerpbins.RData")


ob <- load(file="cache/4.1.B_gerpbins_big0.RData")

#tab <- rbind(tab1, tab2, tab3, tab4, tab5, tab6, tab7, tab8, tab9, tab10)
pdf("graphs/Figure_Sn.gerp1m.pdf", width=7, height=7)
gerPlot(binsize=1000000, tab=tab1m, main="1-Mb bin")
dev.off()

gerPlot(binsize=100000, tab=tab100k, main="100-kb bin")

pdf("manuscript/SI/Figure_Sn.gerp10k.pdf", width=7, height=7)
gerPlot(binsize=10000, tab=tab10k, main="10-kb bin")
dev.off()

gerPlot(binsize=1000, tab=tab1k, main="1-kb bin")





