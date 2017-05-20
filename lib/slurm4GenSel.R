# Jinliang Yang
# Purpose: slurm script for GenSel running
# date: Nov.2nd.2014
# location: farm


slurm4GenSel <- function(sh="largedata/GenSel/CL_test.sh", 
                      sbatho="/home/jolyang/Documents/pvpDiallel/slurm-log/testout-%j.txt",
                      sbathe="/home/jolyang/Documents/pvpDiallel/slurm-log/error-%j.txt",
                      sbathJ="jobid",
                      
                      pi=0.995, findsale ="no",
                      geno="/Users/yangjl/Documents/GWAS2_KRN/SNP/merged/geno_chr", 
                      pheno="/Users/yangjl/Documents/Heterosis_GWAS/pheno2011/reports/cd_GenSel_fullset.txt",
                      map="/Users/yangjl/Documents/linkage.map",
                      chainLength=1000, burnin=100, varGenotypic=1.4, varResidual=2
                      ){
  #####
  wd <- getwd()
  inp <- gsub("sh$", "inp", sh);
  log <- gsub("sh$", "log", sh)
  
  ##### parameters pass to GenSel input
  GenSel_inp(inp=inp, geno=geno, pheno=pheno, map=map, pi=pi, findsale ="no", 
             chainLength=chainLength, burnin=burnin, varGenotypic=varGenotypic, varResidual=varResidual)
  
  #### parameters pass to slurm script
  cat(paste("#!/bin/bash"),
      #-D sets your project directory.
      #-o sets where standard output (of your batch script) goes.
      #-e sets where standard error (of your batch script) goes.
      #-J sets the job name.
      paste("#SBATCH -D", wd, sep=" "),
      paste("#SBATCH -o", sbatho, sep=" "),
      paste("#SBATCH -e", sbathe, sep=" "),
      paste("#SBATCH -J", sbathJ, sep=" "),
      "set -e",
      "set -u",
      "",
      paste("GenSel4R", inp, ">", log),
      paste("python /home/jolyang/bin/send_email.py -s", inp),
      
      file=sh, sep="\n", append=FALSE);
  
      message(paste("###>>> In this path: cd ", wd, sep=""), "\n",
              paste("###>>> note --ntask=x, 8GB of memory per CPU"),"\n",
              paste("###>>> RUN: sbatch -p serial --mem 16000", sh))
}

############################
GenSel_inp <- function(inp="CL_test.inp", pi=0.995, findsale ="no",
                       geno="/Users/yangjl/Documents/GWAS2_KRN/SNP/merged/geno_chr", 
                       pheno="/Users/yangjl/Documents/Heterosis_GWAS/pheno2011/reports/cd_GenSel_fullset.txt",
                       chainLength=1000, burnin=100, varGenotypic=1.4, varResidual=2){
  
    cat(paste("// gensel input file written", Sys.time(), sep=" "), 
        
        "analysisType Bayes",
        "bayesType BayesC",
        paste("chainLength", chainLength, sep=" "),
        paste("burnin", burnin=burnin, sep=" "),
        paste("probFixed", pi, sep=" "),
        
        paste("varGenotypic",  varGenotypic, sep=" "),
        paste("varResidual",  varResidual, sep=" "),
        "nuRes 10",
        "degreesFreedomEffectVar 4",
        "outputFreq 100",
        "seed 1234",
        "mcmcSamples yes",
        "plotPosteriors no",
        paste("FindScale", findsale),
        "modelSequence no",
        "isCategorical no",
        #"linkageMap AGPv2",
        "addMapInfoToMarkers no",
        "windowBV no",
        "",
        "// markerFileName",
        paste("markerFileName", geno, sep=" "), 
        "",
        "// phenotypeFileName",
        paste("phenotypeFileName", pheno, sep=" "),
        "",
        "// mapOrderFileName",
        #paste("mapOrderFileName", map, sep=" "),
        
        file=inp, sep="\n"
    )	
}
