### Jinliang Yang
### 11/20/2014
### check how many SNPs have >0 GERP

### SNP data
library("data.table")
library("plyr")
snp11m <- fread("largedata/SNP/allsnps_11m_gerpv2_tidy.csv", header=TRUE, sep=",") 
# Read 11051839 rows and 19 (of 19) columns from 0.626 GB file in 00:04:53
snp11m$bin1m <- paste(snp11m$chr, round(snp11m$pos/1000000), sep="_" )
snp11m$bin100k <- paste(snp11m$chr, round(snp11m$pos/100000), sep="_" )
snp11m$bin10k <- paste(snp11m$chr, round(snp11m$pos/10000), sep="_" )
snp11m$bin1k <- paste(snp11m$chr, round(snp11m$pos/1000), sep="_" )



cat3 <- subset(snp11m, RS > 2) #28977
cat2 <- subset(snp11m, RS > 1 & RS <= 2) #138478
cat1 <- subset(snp11m, RS > 0 & RS <= 1) #339443
cat0 <- subset(snp11m, is.na(RS) ) #9800436
catn1 <- subset(snp11m, RS > -1 & RS < 0) #393913
catn2 <- subset(snp11m, RS > -2 & RS <= -1) #160300
catn3 <- subset(snp11m, RS > -3 & RS <= -2) #128744
catn4 <- subset(snp11m, RS > -4 & RS <= -3) #15547
catn5 <- subset(snp11m, RS <= -4) #46001








for (i in 1:length(unique(cat2$bin1m)) ){
  idx <- cat2
}











snpnon <- snp11m[!is.na(snp11m$RS),] #1251403 21
snpnz <- subset(snpnon, select = c("snpid", "RS", "MAF", "missing"))
save(list="snpnz", file="largedata/lcache/snpnzRS.RData")

snpb <- subset(snpnon, RS >0 )
#system("mkdir largedata/SNP/inmarker")
write.table(subset(snpb, select="snpid"), "largedata/SNP/inmarker/snp_gerp_bg0.txt", sep="\t", 
            row.names=FALSE, col.names=FALSE, quote=FALSE)

snpb2 <- subset(snpnon, RS >2 )
#system("mkdir largedata/SNP/inmarker")
write.table(subset(snpb2, select="snpid"), "largedata/SNP/inmarker/snp_gerp_bg2.txt", sep="\t", 
            row.names=FALSE, col.names=FALSE, quote=FALSE)

snps4 <- subset(snpnon, RS < -4 )
#system("mkdir largedata/SNP/inmarker")
write.table(subset(snps4, select="snpid"), "largedata/SNP/inmarker/snp_gerp_sm4.txt", sep="\t", 
            row.names=FALSE, col.names=FALSE, quote=FALSE)




#############################################################################

##### generate same number of random SNPs
getRandomSNP <- function(fromdf=snp11m, num=nrow(snpb), nrep=10, basenm="snprandom_set", writeto="largedata/SNP/inmarker/"){
  
  ###
  for(i in 1:nrep){
    idx <- sample(1:nrow(fromdf), num, replace=FALSE)
    tempdf <- snp11m[idx,]
    tempout <- paste0(writeto, basenm, i, ".txt")
    write.table(subset(tempdf, select="snpid"), tempout, sep="\t", 
                row.names=FALSE, col.names=FALSE, quote=FALSE)
  }
  
  message(sprintf("###==> [%s] times of random [%s] SNPs were put in [%s]!", nrep, num, writeto))
  
}

###
set.seed(1234567)
getRandomSNP(fromdf=snp11m, num=nrow(snpb), nrep=10, basenm="snprandom_set", writeto="largedata/SNP/inmarker/")
getRandomSNP(fromdf=snp11m, num=nrow(snpb2), nrep=10, basenm="snp29k_set", writeto="largedata/SNP/inmarker/")
getRandomSNP(fromdf=snp11m, num=nrow(snps4), nrep=10, basenm="snps4_set", writeto="largedata/SNP/inmarker/")





