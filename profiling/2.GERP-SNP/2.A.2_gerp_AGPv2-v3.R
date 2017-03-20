### Jinliang Yang
### 03-19-2017
### check how many SNPs have >0 GERP

### SNP data
library("data.table")
library("plyr")
snp11m <- fread("largedata/SNP/allsnps_11m.bed3", header=FALSE, data.table=FALSE) 
# Read 11051839 rows and 19 (of 19) columns from 0.626 GB file in 00:04:53

snp11m$V1 <- gsub("chr", "", snp11m$V1)
for(i in 1:10){
    sub <- subset(snp11m, V1 == i)
    #sub$snpid <- paste(sub$V1, sub$V3, sep="_")
    write.table(sub, paste0("largedata/SNP/chr", i, ".bed"), row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
}

# convert to AGPv4 http://plants.ensembl.org/Zea_mays/Tools/AssemblyConverter?db=core

gerp <- read.csv("cache/gerpsnp_506898_gp.csv")

snpnon <- snp11m[!is.na(snp11m$RS),] #1251403 21
snpnz <- subset(snpnon, select = c("snpid", "RS", "MAF", "missing"))
save(list="snpnz", file="largedata/lcache/snpnzRS.RData")