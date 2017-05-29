### Jinliang Yang
### Feb 22nd, 2015
### compute GERPs in features


########
gerp11m <- read.csv("largedata/SNP/allsnps_11m_gerpv2_tidy.csv")

bed5 <- gerp11m
bed5$end <- bed5$pos
bed5 <- bed5[, c("chr", "pos", "end", "snpid", "RS")]
names(bed5)[2] <- "start"
bed5$start <- bed5$start - 1

write.table(bed5, "largedata/SNP/allsnps_11m.bed5", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)

#########################################################################################################
### run the following sh
cd largedata/SNP
bedtools intersect -a allsnps_11m.bed5 -b AGPv2_canonical_gene.bed5 -wa > allsnps_11m_in_gene.bed
bedtools intersect -a allsnps_11m.bed5 -b AGPv2_canonical_exon.bed5 -wa > allsnps_11m_in_exon.bed
bedtools intersect -a allsnps_11m.bed5 -b AGPv2_canonical_intron.bed5 -wa > allsnps_11m_in_intron.bed

#########################################################################################################


get_gerp_in_feature <- function(gerp11m=gerp11m){
  gene <- read.table("largedata/SNP/allsnps_11m_in_gene.bed", header=FALSE)
  dim(subset(gene, V5 > 0)) #316983
  dim(subset(gene, V5 < 0)) #491505
  
  exon <- read.table("largedata/SNP/allsnps_11m_in_exon.bed", header=FALSE)
  dim(subset(exon, V5 >0)) #163741
  dim(subset(exon, V5 <0)) #218844
  
  intron <- read.table("largedata/SNP/allsnps_11m_in_intron.bed", header=FALSE)
  dim(subset(intron, V5 >0)) #148355
  dim(subset(intron, V5 <0)) #285472
  
  
  #### sample the same number of >0 and <0 GERP in each feature
  sub1 <- subset(gerp11m, snpid %in% subset(gene, V5 > 0)$V4)
  #> dim(sub1)
  #[1] 313821      5
  write.table(sub1, "largedata/SNP/gerp11m_in_gene_b0.csv", sep=",", quote=FALSE, row.names=FALSE)
  
  sub2 <- subset(gerp11m, snpid %in% subset(gene, V5 < 0)$V4)
  idx2 <- sample(1:nrow(sub2), nrow(sub1))
  write.table(sub2[idx2, ], "largedata/SNP/gerp11m_in_gene_s0.csv", sep=",", quote=FALSE, row.names=FALSE)
  
  sub3 <- subset(gerp11m, snpid %in% subset(exon, V5 > 0)$V4)
  write.table(sub3, "largedata/SNP/gerp11m_in_exon_b0.csv", sep=",", quote=FALSE, row.names=FALSE)
  
  sub4 <- subset(gerp11m, snpid %in% subset(exon, V5 < 0)$V4)
  idx4 <- sample(1:nrow(sub4), nrow(sub3))
  write.table(sub2[idx4, ], "largedata/SNP/gerp11m_in_exon_s0.csv", sep=",", quote=FALSE, row.names=FALSE)
  
  sub5 <- subset(gerp11m, snpid %in% subset(intron, V5 > 0)$V4)
  write.table(sub3, "largedata/SNP/gerp11m_in_intron_b0.csv", sep=",", quote=FALSE, row.names=FALSE)
  
  sub6 <- subset(gerp11m, snpid %in% subset(intron, V5 < 0)$V4)
  idx6 <- sample(1:nrow(sub6), nrow(sub5))
  write.table(sub2[idx6, ], "largedata/SNP/gerp11m_in_intron_s0.csv", sep=",", quote=FALSE, row.names=FALSE)
  
  
}

get_gerp_in_feature(gerp11m=gerp11m)