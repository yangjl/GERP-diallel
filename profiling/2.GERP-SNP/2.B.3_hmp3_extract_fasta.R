### Jinliang Yang
### 10-26-2016
### extract SNP from fasta file

library("data.table")
gerpb0_hmp3 <- function(){
  for(i in 1:10){
    ### get the parental genotype info

    ### bcftools query -S bkn_pvp_samples.txt -f '%CHROM\t%POS\t%REF\t%ALT[\t%TGT]\n' 
    ### chr8_bisnp_n35.vcf.gz > chr8_gt_n35.txt
    
    geno <- fread(paste0("~/dbcenter/HapMap/HapMap3/chr", i, "_gt_n35.txt"), data.table=FALSE)
    ids <- read.table("~/dbcenter/HapMap/HapMap3/bkn_pvp_samples.txt")
    
    names(geno) <- c("chr", "pos", "ref", "alt", as.character(ids$V1))
    bed5 <- geno[, c("chr", "pos", "pos", "ref", "B73")]
    bed5$B73 <- gsub("/.*", "", bed5$B73)
    
    names(bed5) <- c("chrom", "start", "end", "ref", "score")
    bed5$start <- bed5$start - 1
    
    #------- get AGPv3 GERP > 0 sites
    gerpv3 <- fread(paste0("largedata/GERPv3/AGPv3_gerp_chr", i, ".csv"), data.table=FALSE)
    bed <- subset(bed5, end %in% gerpv3$pos)
    bed$start <- format(bed$start, scientific = FALSE)
    bed$end <- format(bed$end, scientific = FALSE)
    
    message(sprintf("###>>> chr [%s], SNPs [%s], GERP>0 [%s], merged [%s]",
                    i, nrow(bed5), nrow(gerpv3), nrow(bed)))
    
    write.table(bed, paste0("largedata/Alignment/GERPv3_b0_chr", i, ".bed"), sep="\t",
                row.names=FALSE, quote=FALSE, col.names=FALSE)
  }
}
########
gerpb0_hmp3()


### -------------
pre_bed5 <- function(){
  for(i in 1:10){
    gfile <- paste0("largedata/Alignment/GERPv3_b0_chr", i, ".bed")
    v3 <- fread(gfile, data.table=FALSE)
    names(v3) <- c("chr", "start", "end", "name", "B73")
    ### pull out seq for different species
    sps <- c("Zea", "Coelorachis","Vossia","Sorghum","Oryza","Setaria",
             "Brachypodium","Hordeum","Musa","Populus","Vitis","Arabidopsis","Panicum")
    v3$B73 <- v3$name
    v3$name <- paste(v3$chr, v3$end, sep="_")
    
    
    bed5 <- data.frame(chr= rep(sps, times=nrow(v3)), start=rep(v3$start, each=length(sps)),
                       end = rep(v3$end, each=length(sps)), name= rep(v3$name, each=length(sps)),
                       score= rep(v3$B73, each=length(sps)))
    bed5$name <- paste(bed5$chr, bed5$name, sep="-")
    write.table(bed5, paste0("largedata/Alignment/GERPv3_chr", i, ".bed"), sep="\t", 
                row.names=FALSE, quote=FALSE, col.names=FALSE)
  }
}

###########
pre_bed5()

# bedtools getfasta -name -tab -fi roast.chrom.1.msa.in -bed GERPv3_chr1.bed -fo hmp3_chr1_gerpv3.txt
# bedtools getfasta -name -tab -fi roast.chrom.2.msa.in -bed GERPv3_chr2.bed -fo hmp3_chr2_gerpv3.txt
# bedtools getfasta -name -tab -fi roast.chrom.3.msa.in -bed GERPv3_chr3.bed -fo hmp3_chr3_gerpv3.txt
# bedtools getfasta -name -tab -fi roast.chrom.4.msa.in -bed GERPv3_chr4.bed -fo hmp3_chr4_gerpv3.txt
# bedtools getfasta -name -tab -fi roast.chrom.5.msa.in -bed GERPv3_chr5.bed -fo hmp3_chr5_gerpv3.txt
# bedtools getfasta -name -tab -fi roast.chrom.6.msa.in -bed GERPv3_chr6.bed -fo hmp3_chr6_gerpv3.txt
# bedtools getfasta -name -tab -fi roast.chrom.7.msa.in -bed GERPv3_chr7.bed -fo hmp3_chr7_gerpv3.txt
# bedtools getfasta -name -tab -fi roast.chrom.8.msa.in -bed GERPv3_chr8.bed -fo hmp3_chr8_gerpv3.txt
# bedtools getfasta -name -tab -fi roast.chrom.9.msa.in -bed GERPv3_chr9.bed -fo hmp3_chr9_gerpv3.txt

# bedtools getfasta -name -tab -fi roast.chrom.10.msa.in -bed GERPv3_chr10.bed -fo hmp3_chr10_gerpv3.txt


### results checking: comparing with B73
res1 <- fread("largedata/Alignment/hmp3_chr10_gerpv3.txt", data.table=FALSE, header=FALSE)
res1$sp <- gsub("-.*", "", res1$V1)
res <- subset(res1, sp %in% "Zea")

geno <- fread("largedata/Alignment/GERPv3_chr10.bed", data.table=FALSE, header=FALSE)
test <- merge(res, geno[, c("V4", "V5")], by.x="V1", by.y="V4", sort=FALSE)

test$V2 <- as.character(test$V2)
test$V5 <- as.character(test$V5)
idx <- which(test$V2 != test$V5 & test$V5 != "N")
length(idx)/nrow(test)
### 0


