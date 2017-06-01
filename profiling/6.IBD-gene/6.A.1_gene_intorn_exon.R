### Jinliang Yang
### Feb 22nd, 2015
### select positive and negative GERP SNPs in gene, exon and intron regions.
#=>> largedata/SNP/allsnps_11m_gerpv2_tidy.csv

###################
source("lib/get_fgs_intron_exon.R")
############
fgs <- get_fgs_intron_exon(gff="~/dbcenter/AGP/AGPv2/ZmB73_5b_FGS.gff",
                           info="~/dbcenter/AGP/AGPv2/ZmB73_5b_FGS_info.txt")

#  [ 994386 ] lines of FGSv2 were loaded!
#  [ 39656 ] canonical trascripts were loaded!
#  [ 39656 ] genes' [ 39656 ]  canonical trascripts were found!
#  [ 39656 ] genes' [ 189308 ]  canonical exons were found!
#  [ 29225 ] genes' [ 149652 ]  canonical introns were found!

### BED5: chrom, start, end, name, and score
genebed <- fgs[[1]][, c("seqname", "start", "end", "geneid", "feature")]
genebed$start <- genebed$start - 1
write.table(genebed, "largedata/SNP/AGPv2_canonical_gene.bed5", 
            row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")


exonbed <- fgs[[2]][, c("seqname", "start", "end", "txid", "feature")]
exonbed$start <- exonbed$start - 1
write.table(exonbed, "largedata/SNP/AGPv2_canonical_exon.bed5", 
            row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")


intronbed <- fgs[[3]][, c("seqname", "start", "end", "txid", "feature")]
intronbed$start <- intronbed$start - 1
write.table(intronbed, "largedata/SNP/AGPv2_canonical_intron.bed5", 
            row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

