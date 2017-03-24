### Jinliang Yang
### change multiple alignment Major allele as beneficial allele rather than B73 allele


del <- read.csv("largedata/Alignment/conserved_alleles_AGPv2.csv")
idx <- which(as.character(del$Zea) != as.character(del$major))
length(idx) #68195
nrow(del) #506839

library("data.table", lib="~/bin/Rlib")
dsf <- fread("largedata/SNP/allsnps_11m.dsf5", header=TRUE)
dsf <- as.data.frame(dsf)

dsf0 <- subset(dsf, snpid %in% del$snpid)


test <- merge(del, dsf0[, 1:8], by="snpid")
idx <- which(test$Zea != test$B73 & test$B73 != "N")

####
names(del)[2] <- "ben"
dsf1 <- merge(del[,1:2], dsf0, by="snpid")
nrow(dsf1[dsf1$major != dsf1$ben & dsf1$ben != "N", ]) #133030
nrow(dsf1[dsf1$B73 != dsf1$ben & dsf1$ben != "N", ]) #75677

dsf1[dsf1$major != dsf1$ben & dsf1$ben != "N", ]$major <- dsf1[dsf1$major != dsf1$ben & dsf1$ben != "N", ]$minor
write.table(dsf1[, -2], "largedata/SNP/allsnps_newgerp2_50k.dsf7", sep="\t", row.names=FALSE, quote=FALSE)

dsf1$major <- dsf1$B73
write.table(dsf1[, -2], "largedata/SNP/allsnps2_newgerp2_50k.dsf7", sep="\t", row.names=FALSE, quote=FALSE)


###### standard k
kfile <- dsf1[, 1:2]
names(kfile)[2] <- "h"
kfile$h <- 0.5
write.table(kfile, "largedata/newGERPv2/allgeno/k5.txt", sep="\t", row.names=FALSE, quote=FALSE)

