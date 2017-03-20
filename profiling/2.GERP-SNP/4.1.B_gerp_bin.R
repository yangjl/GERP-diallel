### Jinliang Yang
### exploring genome-wide feature of GERP score

library("plyr")
#library("data.table")
GerpBin <- function(chr=chr10, BINSIZE=1000000){

  message(sprintf("Chr: tot length [ %s ], >0 [ %s ], <0 [ %s ]", nrow(chr),
                  nrow(subset(chr, RS>0)), nrow(subset(chr, RS<0)) ))
  
  #### get the bin average RS
  chr$bin <- paste(chr$chr, round(chr$pos/BINSIZE, 0), sep="_")
  tab <- ddply(chr, .(bin), summarise,
               binrs = mean(RS))
  #tab$binrs <- tab$binrs/BINSIZE
  
  tab$chr <- as.numeric(as.character(gsub("_.*", "", tab$bin)))
  tab$pos <- as.numeric(as.character(gsub(".*_", "", tab$bin)))
  return(tab)
  
}

### server and local
ob <- load("largedata/lcache/4.1.A_gerpdis.RData")
chrall <- rbind(chr1, chr2, chr3, chr4, chr5, chr6, chr7, chr8, chr9, chr10)
dim(chrall)
#[1] 130,896,913         4

### binsize = 1000bp
tab1k <- GerpBin(chr=chrall, BINSIZE=1000)
tab10k <- GerpBin(chr=chrall, BINSIZE=10000)
tab100k <- GerpBin(chr=chrall, BINSIZE=100000)
tab1m  <- GerpBin(chr=chrall, BINSIZE=1000000)

save(file="cache/4.1.B_gerpbins.RData", list=c("tab1k", "tab10k", "tab100k", "tab1m"))
#tab1k <- rbbind(tab1, tab2, tab3, tab4, tab5, tab6, tab7, tab8, tab9, tab10)
#plot(x=tab$pos, y=tab$binrs, pch=19, cex=0.3, main="Chr10", xlab="Chr", ylab="Avg RS")


#############################
chrbig <- subset(chrall, RS > 0)
dim(chrbig)
#[1] 86006888        4
### binsize = 1000bp
tab1k <- GerpBin(chr=chrbig, BINSIZE=1000)
tab10k <- GerpBin(chr=chrbig, BINSIZE=10000)
tab100k <- GerpBin(chr=chrbig, BINSIZE=100000)
tab1m  <- GerpBin(chr=chrbig, BINSIZE=1000000)
save(file="cache/4.1.B_gerpbins_big0.RData", list=c("tab1k", "tab10k", "tab100k", "tab1m"))

