# Jinliang Yang
# updated: Jan 25th, 2016
# run the wholeset with all GERP SNPs

source("lib/slurm4GenSel.R")
source("~/Documents/Github/zmSNPtools/Rcodes/set_arrayjob.R")
setup_newbin_array <- function(
  ### note: it is for 7 traits with 3 modes for one random shuffling or real data
  genobase="largedata/SNP/geno_b0_cs/gerpv2_b0_cs0", 
  jobdir="slurm-scripts/get_newbin", inpbase= "cs0",
  ptype="pBPHmax",
  jobbase="run_newbin_job", jobid =1){
  
  ### prior information
  wd <- getwd()
  #test run of the 66 diallel of trait per se with additive model
  ti <- tolower(c("ASI", "DTP", "DTS", "EHT",  "GY", "PHT",  "TW"))
  res <- c(0.38, 0.46, 0.46, 15, 88, 41, 0.64)
  gen <- c(0.18, 5.1, 6.0, 123, 65, 377, 0.82)
  
  dir.create(jobdir, showWarnings = FALSE)
  shcommand <- c()
  for(myti in 1:7){
    for(mode in c("a2", "d2", "h2")){
      ### the first one use gs
      myinp <- paste0(jobdir, "/", inpbase, "_", ti[myti], "_", mode,"_ws", ".inp")
      GenSel_inp(
        inp= myinp, pi=0.999,
        findsale ="no",
        geno=paste0(wd, "/", genobase, "_", ti[myti], "_", mode, ".gs.newbin"), 
        pheno=paste0(wd, "/largedata/pheno/wholeset/", tolower(ti[myti]), "_", ptype, ".txt"),
        chainLength=41000, burnin=1000, 
        varGenotypic = gen[myti], 
        varResidual = res[myti]
      )
      shcommand <- c(shcommand, paste("GenSel4R", myinp))
    }
  }
  
  
  
  #################
  jobstart = jobid
  for(i in 1:length(shcommand)){
    cat(shcommand[i], file=paste0(jobdir, "/", jobbase, jobid, ".sh"), sep="\n", append=FALSE)
    jobid <- jobid + 1
  }
  jobend <- jobid -1
  message(sprintf("###>>> setup array jobs: [ %s - %s]", jobstart, jobend))
  set_arrayjob(shid=paste0(jobdir, "/", jobbase, ".sh"),
               shcode=paste0("sh ", jobdir, "/", jobbase, "$SLURM_ARRAY_TASK_ID.sh"),
               arrayjobs= paste0("1-", jobend),
               wd=NULL, jobid=jobbase, email="yangjl0930@gmail.com")
  
}
#newbin_array_7traits_3modes(genobase="largedata/SNP/geno_b0_cs/gerpv2_b0_cs0",
setup_newbin_array(
  ### note: it is for 7 traits with 3 modes for one random shuffling or real data
  genobase="largedata/SNP/geno_b0_cs/gerpv2_b0_cs0", 
  ptype="perse",
  jobdir="slurm-scripts/gwas_b0", inpbase= "ws",
  jobbase="run_ws", jobid =1)
  
#### BPHmax
setup_newbin_array(
  genobase="largedata/SNP/bph_b0_cs/gerpv2_b0_cs0", 
  ptype="BPHmax",
  jobdir="slurm-scripts/gwas_bph_b0", inpbase= "ws_bph",
  jobbase="run_bph_ws", jobid =1)

######################################################################################
main_res <- function(res = res, ptype="perse"){
  res$file <- gsub("ws_", "", res$file)
  res$file <- gsub("_ws\\.out1", "", res$file)
  
  res$trait <- gsub("_.*", "", res$file)
  res$mode <- gsub(".*_", "", res$file)
  res$ptype <- ptype
  return(res)
}

##### perse results
source("~/Documents/Github/zmSNPtools/Rcodes/collect_gsout.R")
res1 <- collect_gsout(dir = "slurm-scripts/gwas_b0", fileptn ="out")

res1 <- main_res(res=res1, ptype="perse")
write.table(res1, "cache/gerpsnp_wholeset_perse.csv", sep=",", row.names=FALSE, quote=FALSE)

### hph results
res2 <- collect_gsout(dir = "slurm-scripts/gwas_bph_b0", fileptn ="out")

res2 <- main_res(res=res2, ptype="bph")
res2$file <- gsub("bph_", "", res2$file)
res2$trait <- gsub("_.*", "", res2$file)
write.table(res2, "cache/gerpsnp_wholeset_bph.csv", sep=",", row.names=FALSE, quote=FALSE)


