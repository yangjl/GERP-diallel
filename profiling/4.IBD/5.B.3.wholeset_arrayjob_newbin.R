### Jinliang Yang
### Jan. 9th, 2014
### updated: 4/7/2016

library(farmeR)

### convert geno into newbin
gsfiles <- list.files(path="/home/jolyang/Documents/Github/GERP-diallel/largedata/sgeno_100k", 
                      pattern="gs$", full.names=TRUE)

inputdf <- data.frame(
  pi=0.995, geno=gsfiles,
  trainpheno="/home/jolyang/Documents/Github/pvpDiallel/largedata/pheno/CV5fold/gy_train1_sp1.txt",
  testpheno="/home/jolyang/Documents/Github/pvpDiallel/largedata/pheno/CV5fold/gy_test1_sp1.txt",
  chainLength=10, burnin=1, varGenotypic=1.4, varResidual=2,
  out = gsub(".*/|.gs", "", gsfiles)
)

run_GenSel4(inputdf, inpdir="largedata/newGERPv2/allgeno", cmdno = 1,
            email="yangjl0930@gmail.com", runinfo = c(FALSE, "bigmemh", 1) )
###>>> In this path: cd /home/jolyang/Documents/Github/GERP-diallel
###>>> RUN: sbatch -p bigmemh --mem 8G --ntasks=1 --time=8:00:00 slurm-script/run_gensel_array.sh



library(farmeR)

### convert geno into newbin
gsfiles <- list.files(path="/home/jolyang/Documents/Github/GERP-diallel/largedata/sgeno", 
                      pattern="gs$", full.names=TRUE)

inputdf <- data.frame(
    pi=0.995, geno=gsfiles,
    trainpheno="/home/jolyang/Documents/Github/pvpDiallel/largedata/pheno/CV5fold/gy_train1_sp1.txt",
    testpheno="/home/jolyang/Documents/Github/pvpDiallel/largedata/pheno/CV5fold/gy_test1_sp1.txt",
    chainLength=10, burnin=1, varGenotypic=1.4, varResidual=2,
    out = gsub(".*/|.gs", "", gsfiles)
)

run_GenSel4(inputdf, inpdir="largedata/sgeno/", cmdno = 1,
            email="yangjl0930@gmail.com", runinfo = c(FALSE, "bigmemh", 1) )
###>>> In this path: cd /home/jolyang/Documents/Github/GERP-diallel
###>>> RUN: sbatch -p bigmemh --mem 8G --ntasks=1 --time=8:00:00 slurm-script/run_gensel_array.sh






