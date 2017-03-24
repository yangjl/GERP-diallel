### Jinliang Yang
### Jan. 9th, 2014
### updated: 4/7/2016

library(farmeR)

### convert geno into newbin
gsfiles <- list.files(path="/home/jolyang/Documents/Github/pvpDiallel/largedata/newGERPv2/allgeno", 
                      pattern="gs$", full.names=TRUE)

inputdf <- data.frame(
  pi=0.995, geno=gsfiles,
  trainpheno="/home/jolyang/Documents/Github/pvpDiallel/largedata/pheno/CV5fold/gy_train1_sp1.txt",
  testpheno="/home/jolyang/Documents/Github/pvpDiallel/largedata/pheno/CV5fold/gy_test1_sp1.txt",
  chainLength=10, burnin=1, varGenotypic=1.4, varResidual=2,
  out = gsub(".*/|.gs", "", gsfiles)
)

run_GenSel4(inputdf, inpdir="largedata/newGERPv2/allgeno", cmdno = 1,
            email="yangjl0930@gmail.com", runinfo = c(TRUE, "med", 1) )

###################################
### convert geno into newbin
gsfiles <- list.files(path="/home/jolyang/Documents/Github/pvpDiallel/largedata/newGERPv2/allgeno_k", 
                      pattern="gs$", full.names=TRUE)

inputdf <- data.frame(
  pi=0.995, geno=gsfiles,
  trainpheno="/home/jolyang/Documents/Github/pvpDiallel/largedata/pheno/CV5fold/gy_train1_sp1.txt",
  testpheno="/home/jolyang/Documents/Github/pvpDiallel/largedata/pheno/CV5fold/gy_test1_sp1.txt",
  chainLength=10, burnin=1, varGenotypic=1.4, varResidual=2,
  out = gsub(".*/|.gs", "", gsfiles)
)

run_GenSel4(inputdf, inpdir="largedata/newGERPv2/allgeno_k", cmdno = 21,
            email="yangjl0930@gmail.com", runinfo = c(TRUE, "bigmemh", 1) )

###################################
### convert geno into newbin
gsfiles <- list.files(path="/home/jolyang/Documents/Github/pvpDiallel/largedata/newGERPv2/allgeno_kbph", 
                      pattern="gs$", full.names=TRUE)

inputdf <- data.frame(
  pi=0.995, geno=gsfiles,
  trainpheno="/home/jolyang/Documents/Github/pvpDiallel/largedata/pheno/CV5fold/gy_train1_sp1.txt",
  testpheno="/home/jolyang/Documents/Github/pvpDiallel/largedata/pheno/CV5fold/gy_test1_sp1.txt",
  chainLength=10, burnin=1, varGenotypic=1.4, varResidual=2,
  out = gsub(".*/|.gs", "", gsfiles)
)

run_GenSel4(inputdf, inpdir="largedata/newGERPv2/allgeno_k", cmdno = 21,
            email="yangjl0930@gmail.com", runinfo = c(TRUE, "serial", 1) )
