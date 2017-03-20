### Jinliang Yang
### August 21th, 2016


library("farmeR")
cmd <- c("cd /home/jolyang/dbcenter/HapMap/HapMap3")
for(i in 1:10){
  tmp <- paste0("bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/DP\t%INFO/NZ\t%INFO/MAF\n'",
               " merged_flt_c", i, ".vcf.gz > chr", i, "_frq.txt")
  cmd <- c(cmd, tmp)
}

set_farm_job(slurmsh = "slurm-script/getmaf.sh",
             shcode = cmd, wd = NULL, jobid = "getmaf",
             email = "yangjl0930@gmail.com", runinfo = c(TRUE, "bigmemh", "8"))

###### main scripts
library("data.table")

gp <- readgerp()
write.table(gp, "largedata/GERPv3_hmp3.csv", sep=",", row.names=FALSE, quote=FALSE)

##############################
readgerp <- function(){
  
  res <- data.frame()
  #### read the table fast
  for(i in 1:10){
    
    ### get allele freq and keep only biallelic sites
    frq <- fread(paste0("~/dbcenter/HapMap/HapMap3/chr", i, "_frq.txt"))
    names(frq) <- c("chr", "pos", "ref", "alt", "dp", "nz", "maf")
    frq$nch <- nchar(frq$alt)
    frq <- subset(frq, nch == 1)

    ### read v3 sites
    gerpfile= paste0("/group/jrigrp/gerp/GERPv3/roast.chrom.", i, ".msa.in.rates.full")
    chr <- fread(gerpfile, header=FALSE)
    names(chr) <- c("N", "RS")
    chr$pos <- 1:nrow(chr)
    
    ###
    if( max(chr$pos) > max(frq$pos) ){
      message(sprintf("###>>> chr [ %s ] OK!", i))
    }else{
      message(sprintf("###>>> GERP last base: [ %s ] and frq last base: [ %s ]", max(chr$pos), max(frq$pos)))
    }
    
    chr <- subset(chr, pos %in% frq$pos)
    achr <- merge(chr, frq, by="pos")
    achr <- as.data.frame(achr)
    achr <- achr[, c(-7, -10)]
    
    res <- rbind(res, achr)
  }
  return(res)
}
