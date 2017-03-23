### Jinliang Yang
### 9/3/2015
### purpose: find the deleterious per line

gerpsnp <- read.csv("largedata/GERPv2/gerpsnp_506898.csv")
gerpsnp <- subset(gerpsnp, B73 != "N")


#### compute variation of GERP per cM
gerp_var <- function(gerpsnp, binsize=1){
  df <- gerpsnp
  df$B73 <- as.character(df$B73)
  for(i in 11:ncol(df)){
    df[,i] <- as.character(df[,i])
    df[df[,i] != df$B73 & df[,i] != "N", i] <- df[df[,i] != df$B73 & df[,i] != "N", ]$RS
    df[df[,i] == df$B73 | df[,i] == "N", i] <- 0  
  }
  
  ### compute the sum of GERP in genetic bins
  df$bin <- paste(df$chr, round(df$genetic/binsize, 0), sep="_")
  mydf <- df[, 11:ncol(df)]
  mydf[, 1:11] <- apply(mydf[, 1:11], 2, as.numeric)
  
  binsum1 <- binsum2 <- data.frame()
  for(bini in unique(mydf$bin)){
    tem <- subset(mydf, bin==bini)
    out2 <- out1 <- tem[1, ]
    out1[, -12] <- apply(tem[, -12], 2, sum)
    out2[, -12] <- apply(tem[, -12], 2, function(x) sum(x != 0))
    binsum1 <- rbind(binsum1, out1)
    binsum2 <- rbind(binsum2, out2)
  }
  
  #binsum1$var <- apply(as.matrix(binsum[, 1:11]), 1, var)
  binsum1$mean <- apply(as.matrix(binsum1[, 1:11]), 1, mean)
  binsum2$mean <- apply(as.matrix(binsum2[, 1:11]), 1, mean)
  return(list(binsum1, binsum2))
}

####
library("plyr", lib="~/bin/Rlib/")
out <- gerp_var(gerpsnp, binsize=1)

gerpmean <- out[[1]]
write.table(gerpmean, "data/gerpmean_cm.csv", sep=",", row.names=FALSE, quote=FALSE )

gerpcount <- out[[2]]
write.table(gerpcount, "data/gerpcount_cm.csv", sep=",", row.names=FALSE, quote=FALSE )


###############################################################################################


gerp_mc_plot <- function(tab=gerpmean){
  
  source("~/Documents/Github/zmSNPtools/Rcodes/rescale.R")
  
  tab <- tab[order(tab$chr, tab$pos), ]
  plot(c(0, max(tab$pos)), c(10,110), type="n", 
       xlab="Genetic Distance (cM)", ylab="", yaxt="n", bty="n")
  
  axis(side=2, tick =FALSE, las=1, at=c(100, 90, 80, 70, 60, 50, 40, 30, 20, 10), 
       labels=paste("Chr", 1:10, sep=""), line=-1.5)
  #### chromosome
  for (i in 1:10){
    lines(c(0, max(subset(tab, chr==i)$pos)), c(100-10*(i-1), 100-10*(i-1)), lwd=2, col="grey")
    #lines (c(centromere[i,]$Start,centromere[i,]$End),
    #       c(105-10*i, 105-10*i),lwd=5, col="tomato") 
  }
  
  cent <- read.csv("data/centromere_AGPv2.csv")
  ### core plot
  cols <- rep(c("slateblue", "cyan4"), 5)
  for (chri in 1:10){
    mytab <- subset(tab, chr == chri & pos > 0)
    subc <- subset(cent, chr == chri)
    mytab$mean <- rescale(mytab$mean, c(0, 9))
    rect(xleft=subc$cm_start, ybottom= 110 - 10*chri, xright=subc$cm_end, ytop=109 - 10*(chri-1),
         col="red", border="red")
    points(mytab$pos, 100 - 10*(chri-1) + mytab$mean, pch=19, cex=0.5, col=cols[chri])
    
  } 
}

#######################
gerpmean <- read.csv("data/gerpmean_cm.csv")
gerpmean$chr <- as.numeric(as.character(gsub("_.*", "", gerpmean$bin)))
gerpmean$pos <- as.numeric(as.character(gsub(".*_", "", gerpmean$bin)))


#tab <- rbind(tab1, tab2, tab3, tab4, tab5, tab6, tab7, tab8, tab9, tab10)
pdf("graphs/suppl_del_chrs.pdf", width=10, height=5)
gerp_mc_plot(tab=gerpmean)
dev.off()
