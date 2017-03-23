### Jinliang Yang
### 03-19-2017
### check how many SNPs have >0 GERP


gsnp <- read.csv("cache/gerpsnp_506898_gp.csv")

bed <- gsnp[, 1:2]
bed$chr <- gsub("_.*", "", bed$marker)
bed$start <- gsub(".*_", "", bed$marker)
bed$end <- bed$start

bed$start <- as.numeric(as.character(bed$start)) - 1

bed <- bed[, c("chr", "start", "end", "marker")]
write.table(bed, paste0("largedata/SNP/gerpsnp_v2.bed"), 
            row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
# convert to v2->v4 and v4->v3 AGPv4 
# http://plants.ensembl.org/Zea_mays/Tools/AssemblyConverter?db=core

