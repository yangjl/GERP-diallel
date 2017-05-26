### Jinliang Yang
### purpose: random sampled k=d/a
### updated: 05-22-2017

### Get inputdf for each mode
get_inputdf <- function(mygeno=a2, phenopwd){
  #7traits x 3 modes x 12 permuation x  5 fold 100 rand x
  inputdf <- data.frame()
  
  #test run of the 66 diallel of trait per se with additive model
  ti <- tolower(c("ASI", "DTP", "DTS", "EHT",  "GY", "PHT",  "TW"))
  res <- c(0.38, 0.46, 0.46, 15, 88, 41, 0.64)
  gen <- c(0.18, 5.1, 6.0, 123, 65, 377, 0.82)
  
  for(i in 1:7){ #traits
    idx <- grep(ti[i], mygeno)
      
    tem <- data.frame(
        pi = 0.995, geno = mygeno[idx],
        pheno = paste0(phenopwd, ti[i], "_perse.txt"),
        chainLength=41000, burnin=1000, varGenotypic=gen[i], varResidual=res[i]
    ) 
    inputdf <- rbind(inputdf, tem)
  }
  
  return(inputdf)
}

library(farmeR)
####### For trait perse
phenopwd <- "/home/jolyang/Documents/Github/GERP-diallel/largedata/pheno/wholeset/"
h2 <- list.files(path="/home/jolyang/Documents/Github/GERP-diallel/largedata/sgeno", 
                 pattern="h2.gs.newbin$", full.names=TRUE)
### add
inputdf1 <- get_inputdf(mygeno=h2, phenopwd)
inputdf1$out <- gsub(".*\\/|.gs.newbin", "", inputdf1$geno)

run_GenSel4(inputdf=inputdf1, cv=FALSE, inpdir="largedata/sgeno", cmd=1,
            email="yangjl0930@gmail.com", runinfo = c(FALSE, "serial", 1) )

###>>> In this path: cd /home/jolyang/Documents/Github/GERP-diallel
###>>> RUN: sbatch -p bigmemm --mem 6G --ntasks=1 --time=8:00:00 slurm-script/run_gensel_array.sh






