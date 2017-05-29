# Jinliang Yang
# April 21, 2015
# run the wholeset with GERP SNPs

source("lib/RunWhoeSet_genic.R")
####

mysh <- RunWholeSet_genic(inppwd="slurm-scripts/wholeset/")
setUpslurm(slurmsh="slurm-scripts/wholeset_testrun.sh",
           codesh=mysh,
           wd=NULL, jobid="testrun", email=TRUE)
  
###>>> In this path: cd /home/jolyang/Documents/Github/pvpDiallel
###>>> [ note: --ntasks=INT, number of cup ]
###>>> [ note: --mem=16000, 16G memory ]
###>>> RUN: sbatch -p bigmemh --ntasks 2 slurm-scripts/wholeset_testrun.sh



source("~/Documents/Github/zmSNPtools/Rcodes/collect_gsout.R")
res <- collect_gsout(dir = "slurm-scripts/genic_wholeset/", fileptn ="out")

main_res <- function(res = res){
  res$trait <- gsub("_.*", "", res$file)
  res$transf <- gsub("_1.*", "", res$file)
  res$transf <- gsub(".*_", "", res$transf)
  res$mode <- gsub(".*_1", "", res$file)
  res$mode <- gsub("\\.out1", "", res$mode)
  return(res)
}

#####
res <- main_res(res=res)
write.table(res, "cache/gs_wholeset_h2.csv", sep=",", row.names=FALSE, quote=FALSE)

gs <- read.csv("cache/gs_wholeset_h2.csv")


