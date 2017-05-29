# Jinliang Yang
# April 21, 2015
# run the wholeset with GERP SNPs

source("lib/RunWhoeSet_genic.R")
####
source("~/Documents/Github/zmSNPtools/Rcodes/setUpslurm.R")

gs2 <- read.csv("cache/gs_wholeset_h2.csv")
mysh <- RunWholeSet_genic(inppwd="slurm-scripts/genic_wholeset/", priors=gs2)

for(i in 1:10){
  setUpslurm(slurmsh= paste0("slurm-scripts/genic_wholeset_run", i, ".sh"),
             codesh= mysh[(7*(i-1)+1) : (7*i)],
             wd=NULL, jobid= paste0("genic_run", i), email="yangjl0930@gmail.com")
  
}

###>>> In this path: cd /home/jolyang/Documents/Github/pvpDiallel
###>>> [ note: --ntasks=INT, number of cup ]
###>>> [ note: --mem=16000, 16G memory ]
###>>> RUN: sbatch -p bigmemh --ntasks=1 slurm-scripts/genic_wholeset_run1.sh


source("~/Documents/Github/zmSNPtools/Rcodes/collect_gsout.R")
res <- collect_gsout(dir = "slurm-scripts/wholeset", fileptn ="out")

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


