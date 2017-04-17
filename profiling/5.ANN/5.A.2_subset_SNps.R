### Jinliang Yang
### 04/16/2017
### subset SNPs from VCF file


library(data.table)
del <- fread("largedata/gerpsnp_v3_345176_del.csv", data.table=FALSE)

r <- del[, c("snpid", "chr", "pos")]
r$chr <- gsub("_.*", "", r$snpid)
r$pos <- gsub(".*_", "", r$snpid)
write.table(r[, c("chr", "pos")], "/home/jolyang/dbcenter/HapMap/HapMap3//gerpsnp_v3_chr_pos.txt", 
            sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

####
# cp gerpsnp_v3_chr_pos.txt /home/jolyang/dbcenter/HapMap/HapMap3/

