### Jinliang Yang
### August 19th, 2015
### get the genetic load in landrace and pvp diallel lines

library("data.table")
### get the parental genotype info
get_del_count <- function(){
  #### v2 and v3 coordinates conversion
  v2 <- read.table("largedata/Alignment/GERP_b0_AGPv2.bed", header=FALSE)
  v3 <- read.table("data/output_GERP_b0_AGPv2.bed")
  #sum(v3$V4 %in% v2$V4)/nrow(v2)
  #[1] 0.9999231
  names(v3) <- c("chr", "start", "end", "v2id", "B73")
  #end: v3 coordinates
  #v2id: v2 chr_pos
  v3$v3id <- paste(v3$chr, v3$end, sep="_")
  message(sprintf("###>> number of same id between v2 and v3: [ %s ]", sum(with(v3, v2id == v3id))))
  
  ### 
  pvpgeno <- fread("largedata/GERPv2/gerpsnp_506898.csv", data.table=FALSE, header=TRUE)
  pvpv3 <- merge(v3[, c("v2id", "v3id")], pvpgeno, by.x="v2id", by.y="marker")
  
  outf <- data.frame()
  for(i in 1:10){
    pchr <- subset(pvpv3, chr == i)
    del <- fread(paste0("largedata/del_allele_frq_chr", i, ".txt"), header=TRUE)
    #note minor alleles here are obtained from alignment, should (or non-major) be deleterious 
    out <- merge(del, pchr, by.x="snpid", by.y="v3id")
    message(sprintf("###>>> chr [ %s ]: [ %s ] percentage of overlapped (or [ %s ] / [ %s ])",
                    i, 100*round(nrow(out)/nrow(pchr), 3), nrow(out), nrow(pchr)))
    outf <- rbind(outf, out)
  }
  
  return(as.data.frame(outf))
  
}

##### GET interested sites
outf <- get_del_count()
bed <- outf[, c("snpid", "chr", "pos")]
bed$chr <- gsub("_.*", "", bed$snpid)
bed$pos <- gsub(".*_", "", bed$snpid)
bed$start <- as.numeric(bed$pos) -1
bed$end <- as.numeric(bed$pos)
bed$start <- format(bed$start, scientific =FALSE)
bed$end <- format(bed$end, scientific =FALSE)

write.table(bed[, c("chr", "start", "end")], "~//dbcenter/HapMap/HapMap3/pvp_sites_v3.bed", 
            sep="\t", row.names=FALSE, quote=FALSE, col.names=FALSE)


### GET LOAD
idx1 <- which(names(outf) == "B73")
idx2 <- ncol(outf)
outk <- data.frame()
for(j in idx1:idx2){
  #outf[, j] <- as.character()
  tem <- data.frame(parent=names(outf)[j], 
                    #num_del=sum(outf$major.x != outf[, j] & outf$major.x != "N" & outf[,j] != "N"),
                    num_del=sum(outf$minor.x == outf[, j] & outf[,j] != "N"))
  outk <- rbind(outk, tem)
}




##### get landrace lines
line <- read.table("~/dbcenter/HapMap/HapMap3/vcf.header")
line <- as.vector(t(line))
idx <- grep("BKN", line)
bkn <- line[idx]
write.table(data.frame(v1=bkn), "~/dbcenter/HapMap/HapMap3/bkn_samples.txt", 
            sep="\t", row.names=FALSE, quote=FALSE, col.names=FALSE)

#[1] "BKN009" "BKN010" "BKN011" "BKN014" "BKN015" "BKN016" "BKN017" "BKN018"
#[9] "BKN019" "BKN020" "BKN022" "BKN023" "BKN025" "BKN026" "BKN027" "BKN029"
#[17] "BKN030" "BKN031" "BKN032" "BKN033" "BKN034" "BKN035" "BKN040"


### RUNNINg
# bcftools query -R pvp_sites_v3.bed -S bkn_samples.txt -f 
cmd <- paste("bcftools query -R pvp_sites_v3.bed -S bkn_samples.txt -f", 
             "'%CHROM\t%POS[\t%SAMPLE]\n' merged_flt_c10.vcf.gz")

library("farmeR")
cmd <- c("cd /home/jolyang/dbcenter/HapMap/HapMap3")
for(i in 1:9){
  ext <- paste0("pigz -d -p 16 merged_flt_c", i, ".vcf.gz")
  bgzip <- paste0("bgzip merged_flt_c", i, ".vcf -@ 16")
  idx <- paste0("tabix -p vcf merged_flt_c", i, ".vcf.gz")
  tmp <- paste0("bcftools query -R pvp_sites_v3.bed -S bkn_samples.txt",
                " -f '%CHROM\t%POS\t[%SAMPLE\t]\n'",
                " merged_flt_c", i, ".vcf.gz > chr", i, "_frq.txt")
  tmp2 <- paste0("bcftools query -R pvp_sites_v3.bed -S bkn_samples.txt",
                " -f '%CHROM\t%POS\t[%GT\t]\n'",
                " merged_flt_c", i, ".vcf.gz > chr", i, "_bkn_gt.txt")
  cmd <- c(cmd, ext, bgzip, idx)
}

set_farm_job(slurmsh = "slurm-script/getbzip.sh",
             shcode = cmd, wd = NULL, jobid = "bzip",
             email = "yangjl0930@gmail.com", runinfo = c(TRUE, "bigmemm", 16))


# bcftools query -R pvp_sites_v3.bed -S bkn_samples.txt -f 
cmd <- paste("bcftools query -S bkn_samples.txt -f", 
             "'%CHROM\t%POS\t%REF\t%ALT[\t%TGT]\n' merged_flt_c10.vcf.gz > out.txt")


chr10 <- fread("~/dbcenter/HapMap/HapMap3/genotype_chr10.txt")
