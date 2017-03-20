### Jinliang Yang
### 10/26/2016
### purpose: find GT of the given SNPs in hmp3

library("farmeR")

for(i in 1:9){
  cmd <- c("cd /home/jolyang/dbcenter/HapMap/HapMap3")
  #ext <- paste0("pigz -d -p 16 merged_flt_c", i, ".vcf.gz")
  #bgzip <- paste0("bgzip merged_flt_c", i, ".vcf -@ 16")
  #idx <- paste0("tabix -p vcf merged_flt_c", i, ".vcf.gz")
  
  ## annotate VCF header becasue there is an error FORMAT=>INFO for AD
  b0 <- paste0("bcftools annotate -h change_AD.txt merged_flt_c", i, ".vcf.gz -Oz -o merged_flt_ad_c", i, ".vcf.gz")
  ## keep only the biallelic SNPs for N=35 samples
  b1 <- paste0("bcftools view merged_flt_ad_c", i, ".vcf.gz ", 
               "-m2 -M2 -v snps -S bkn_pvp_samples.txt -Oz -o ", "chr", i, "_bisnp_n35.vcf.gz") 
  ## convert vcf.gz into SNP txt
  b2 <- paste0("bcftools query -S bkn_pvp_samples.txt",
               " -f '%CHROM\t%POS\t%REF\t%ALT[\t%TGT]\n'",
               " chr", i, "_bisnp_n35.vcf.gz > chr", i, "_gt_n35.txt")
  
  shid <- paste0("slurm-script/run_", i, ".sh")
  cmd <- c(cmd, b0, b1, b2)
  cat(cmd, file=shid, sep="\n", append=FALSE)
}

##########
shcode <- "sh slurm-script/run_$SLURM_ARRAY_TASK_ID.sh"

set_array_job(shid="slurm-script/run.sh", shcode=shcode,
              arrayjobs="1-9", wd=NULL, jobid="getvcf", email="yangjl0930@gmail.com",
              run = c(FALSE, "bigmemh", 1, "8G"))
###>>> In this path: cd /home/jolyang/Documents/Github/pvpDiallel
###>>> RUN: sbatch -p bigmemh --mem 8G --ntasks=1 --time=5:00:00 slurm-script/run.sh







