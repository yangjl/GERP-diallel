### Jinliang Yang
### Jan. 9th, 2014

source("lib/cv_array_jobs.R")
############# GERP > 0, BPH #############
### 70 (10cs x 7 traits) array jobs, =>gerpid x 3modes
for(i in 1:11){
  setup_gerpibd_array_7traits(
    outdir="largedata/SNP/bph_b0_cs", jobbase="gerpid_bph_job", jobid = 7*(i-1)+1,
    kfile_path="largedata/snpeff/BPH/",
    genobase= paste0("largedata/SNP/bph_b0_cs/gerpv2_b0_cs", i-1))
}

## note: it is for 7 traits with 3 modes for one random shuffling or real data
for(i in 1:11){
  setup_newbin_array(
    genobase= paste0("largedata/SNP/bph_b0_cs/gerpv2_b0_cs", i-1), jobid=21*(i-1)+1,
    phenobase="largedata/pheno/CV5fold_BPHmax",
    jobdir="largedata/SNP/bph_b0_cs/get_newbin", inpbase= paste0("cs",i -1),
    jobbase="run_newbin_job")
}

##### gensel: 10 sp x (7traits x 5 cv x 3 modes)
for(i in 1:11){
  setup_gensel_array(
    outdir="largedata/SNP/bph_b0_cs/", jobbase="bph_gs_job", jobid=100*(i-1)+1,
    inpbase= paste0("largedata/SNP/bph_b0_cs/cs", i-1),
    phenobase="largedata/pheno/CV5fold_BPHmax",
    genobase= paste0("largedata/SNP/bph_b0_cs/gerpv2_b0_cs", i-1))
}


###>>> setup gensel array jobs: [ 1001 - 1100]
###>>> In this path: cd /home/jolyang/Documents/Github/pvpDiallel
###>>> [ note: --ntasks=INT, number of cup ]
###>>> [ note: --mem=16000, 16G memory ]
###>>> RUN: sbatch -p serial largedata/SNP/bph_b0_cs//bph_gs_job.sh









