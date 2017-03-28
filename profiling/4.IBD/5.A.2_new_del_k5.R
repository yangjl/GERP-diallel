### Jinliang Yang
### change multiple alignment Major allele as beneficial allele rather than B73 allele
### udpated: 03-27-2017

del <- read.csv("largedata/Alignment/conserved_alleles_AGPv2.csv")
idx <- which(as.character(del$Zea) != as.character(del$major))
length(idx) #68195
nrow(del) #506839

del <- fread("largedata/gerpsnp_v3_345176_del.csv", data.table=FALSE)
### 345176     21
names(del)[2] <- "ben"
############
library("data.table")
dsf <- fread("largedata/SNP/allsnps_11m.dsf5", header=TRUE)
dsf <- as.data.frame(dsf)

dsf0 <- subset(dsf, snpid %in% del$marker)


dsf1 <- del
dsf1 <- dsf1[, c(3,6,7,2,1,8,4, 10:ncol(dsf1)) ]
#names(dsf1)[1:7] <- c("snpid","chr", "pos", "major","minor", "MAF","missing")
names(dsf1)[1:7] <- c("snpid","chr", "pos", "major","snpidv3", "MAF", "RSv3")
write.table(dsf1, "largedata/SNP/allsnps_gerpv3_snpidv2_345k.dsf7", sep="\t", row.names=FALSE, quote=FALSE)


###### standard k
kfile <- dsf1[, c("snpid", "RSv3")]
names(kfile)[1:2] <- c("snpid", "RS")
#names(kfile)[2] <- "h"
#kfile$h <- 0.5
write.table(kfile, "largedata/GERPv3/gerp_snpidv2_RSv3.txt", sep="\t", row.names=FALSE, quote=FALSE)


###### standard k
kfile <- dsf1[, 1:2]
names(kfile)[2] <- "h"
kfile$h <- 0.5
write.table(kfile, "largedata/newGERPv2/allgeno/k5.txt", sep="\t", row.names=FALSE, quote=FALSE)

