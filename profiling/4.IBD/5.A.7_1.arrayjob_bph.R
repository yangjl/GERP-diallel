### Jinliang Yang
### Jan. 9th, 2014

### Get inputdf for each mode
get_inputdf <- function(mygeno=a2, phenopwd){
  #7traits x 3 modes x 12 permuation x  5 fold 100 rand x
  inputdf <- data.frame()
  
  #test run of the 66 diallel of trait per se with additive model
  ti <- tolower(c("ASI", "DTP", "DTS", "EHT",  "GY", "PHT",  "TW"))
  res <- c(0.38, 0.46, 0.46, 15, 88, 41, 0.64)
  gen <- c(0.18, 5.1, 6.0, 123, 65, 377, 0.82)
  
  for(i in 2:3){ #traits
    for(p in 1:12){ #gerp
      for(f in 1:5){ #5-fold
        for(r in 1:100){ #pheno
          tem <- data.frame(
            pi = 0.995, geno = mygeno[p],
            trainpheno = paste0(phenopwd, ti[i], "_train", f, "_sp", r, ".txt"),
            testpheno = paste0(phenopwd, ti[i], "_test", f, "_sp", r, ".txt"),
            chainLength=11000, burnin=1000, varGenotypic=gen[i], varResidual=res[i]
          ) 
          inputdf <- rbind(inputdf, tem)
        }
      }
    }
  }
  return(inputdf)
}

library(farmeR)
####### For trait perse
phenopwd <- "/home/jolyang/Documents/Github/GERP-diallel/largedata/pheno/CV5fold_MPH/"
a2 <- list.files(path="/home/jolyang/Documents/Github/GERP-diallel/largedata/newGERPv2/allgeno", 
                 pattern="a2.gs.newbin$", full.names=TRUE)
d2 <- list.files(path="/home/jolyang/Documents/Github/GERP-diallel/largedata/newGERPv2/allgeno", 
                 pattern="d2.gs.newbin$", full.names=TRUE)
h2 <- list.files(path="/home/jolyang/Documents/Github/GERP-diallel/largedata/newGERPv2/allgeno", 
                 pattern="h2.gs.newbin$", full.names=TRUE)
### add
inputdf1 <- get_inputdf(mygeno=a2, phenopwd)

inputdf1 <- read.csv("largedata/newGERPv2/inputdf_a2_bph_42000.csv")
inputdf1$geno <- gsub("pvpDiallel", "GERP-diallel", inputdf1$geno)
inputdf1$trainpheno <- gsub("pvpDiallel", "GERP-diallel", inputdf1$trainpheno)
inputdf1$testpheno <- gsub("pvpDiallel", "GERP-diallel", inputdf1$testpheno)
inputdf1$trainpheno <- gsub("CV5fold_BPH", "CV5fold_MPH", inputdf1$trainpheno)
inputdf1$testpheno <- gsub("CV5fold_BPH", "CV5fold_MPH", inputdf1$testpheno)

inputdf1$out <- paste0(gsub(".*/|.txt", "", inputdf1$trainpheno), 
                       gsub(".*gerpv2_b0|.gs.newbin", "", inputdf1$geno))

write.table(inputdf1, "largedata/newGERPv2/inputdf_a2_mph_42000.csv", sep=",", quote=FALSE)

run_GenSel4(inputdf=inputdf1, inpdir="largedata/newGERPv2/allgeno_mph_a", cmdno=100,
            email="yangjl0930@gmail.com", runinfo = c(FALSE, "serial", 1) )
###>>> In this path: cd /home/jolyang/Documents/Github/GERP-diallel
###>>> RUN: sbatch -p serial --mem 1500 --time=80:00:00 --ntasks=1 slurm-script/run_gensel_array.sh

### dom
inputdf2 <- get_inputdf(mygeno=d2, phenopwd)

inputdf2 <- read.csv("largedata/newGERPv2/inputdf_d2_bph_42000.csv")
inputdf2$geno <- gsub("pvpDiallel", "GERP-diallel", inputdf2$geno)
inputdf2$trainpheno <- gsub("pvpDiallel", "GERP-diallel", inputdf2$trainpheno)
inputdf2$testpheno <- gsub("pvpDiallel", "GERP-diallel", inputdf2$testpheno)
inputdf2$trainpheno <- gsub("CV5fold_BPH", "CV5fold_MPH", inputdf2$trainpheno)
inputdf2$testpheno <- gsub("CV5fold_BPH", "CV5fold_MPH", inputdf2$testpheno)

inputdf2$out <- paste0(gsub(".*/|.txt", "", inputdf2$trainpheno), 
                       gsub(".*gerpv2_b0|.gs.newbin", "", inputdf2$geno))

write.table(inputdf2, "largedata/newGERPv2/inputdf_d2_mph_42000.csv", sep=",", quote=FALSE)


inputdf2 <- read.csv("largedata/newGERPv2/inputdf_d2_mph_42000.csv")
run_GenSel4(inputdf=inputdf2, inpdir="largedata/newGERPv2/allgeno_mph_d", cmdno=100,
            shid = "slurm-script/run_gensel_d2_array.sh",
            email="yangjl0930@gmail.com", runinfo = c(FALSE, "bigmemm", 1) )
###>>> In this path: cd /home/jolyang/Documents/Github/pvpDiallel
###>>> RUN: sbatch -p serial --mem 1500 --ntasks=1 --time=80:00:00 slurm-script/run_gensel_d2_array.sh

### k=0.5
inputdf3 <- get_inputdf(mygeno=h2, phenopwd)

inputdf3 <- read.csv("largedata/newGERPv2/inputdf_h2_bph_42000.csv")
inputdf3$geno <- gsub("pvpDiallel", "GERP-diallel", inputdf3$geno)
inputdf3$trainpheno <- gsub("pvpDiallel", "GERP-diallel", inputdf3$trainpheno)
inputdf3$testpheno <- gsub("pvpDiallel", "GERP-diallel", inputdf3$testpheno)
inputdf3$trainpheno <- gsub("CV5fold_BPHmax", "CV5fold_MPH", inputdf3$trainpheno)
inputdf3$testpheno <- gsub("CV5fold_BPHmax", "CV5fold_MPH", inputdf3$testpheno)

inputdf3$out <- paste0(gsub(".*/|.txt", "", inputdf3$trainpheno), 
                       gsub(".*gerpv2_b0|.gs.newbin", "", inputdf3$geno))

write.table(inputdf3, "largedata/newGERPv2/inputdf_h2_mph_42000.csv", sep=",", quote=FALSE)

inputdf3 <- read.csv("largedata/newGERPv2/inputdf_h2_mph_42000.csv")
run_GenSel4(inputdf=inputdf3, inpdir="largedata/newGERPv2/allgeno_mph_k5", cmdno=100,
            shid = "slurm-script/run_gensel_h2_array.sh",
            email="yangjl0930@gmail.com", runinfo = c(FALSE, "bigmemm", 1) )

###>>> In this path: cd /home/jolyang/Documents/Github/pvpDiallel
###>>> RUN: sbatch -p bigmemm --mem 5G --ntasks=1 --time=80:00:00 slurm-script/run_gensel_h2_array.sh







