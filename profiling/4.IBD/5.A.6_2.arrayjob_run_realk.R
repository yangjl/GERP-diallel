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
  
  for(i in 1:7){ #traits
    for(p in 1:12){ #gerp
      for(f in 1:5){ #5-fold
        for(r in 1:100){ #pheno
          tem <- data.frame(
            pi = 0.995, geno = mygeno[p, ti[i]],
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
### real k for trait per se
phenopwd <- "/home/jolyang/Documents/Github/pvpDiallel/largedata/pheno/CV5fold/"

h2 <- list.files(path="/home/jolyang/Documents/Github/pvpDiallel/largedata/newGERPv2/allgeno_k", 
                 pattern="h2.gs.newbin$", full.names=TRUE)

df <- data.frame(asi=h2[grep("asi", h2)], dtp=h2[grep("dtp", h2)], dts=h2[grep("dts", h2)], 
                 eht=h2[grep("eht", h2)], gy=h2[grep("gy", h2)], pht=h2[grep("pht", h2)], 
                 tw=h2[grep("tw", h2)])

inputdf3 <- get_inputdf(mygeno=df, phenopwd)
inputdf3$out <- paste0(gsub(".*/|.txt", "", inputdf3$trainpheno), 
                       gsub(".*gerpv2_b0|.gs.newbin", "", inputdf3$geno))

write.table(inputdf, "largedata/newGERPv2/inputdf_realk_perse_42000.csv", sep=",", quote=FALSE)

inputdf <- read.csv("largedata/newGERPv2/inputdf_realk_perse_42000.csv")
run_GenSel4(inputdf=inputdf, inpdir="largedata/newGERPv2/allgeno_perse_k", cmdno=500,
            shid = "slurm-script/run_gensel_realk_array.sh",
            email="yangjl0930@gmail.com", runinfo = c(TRUE, "serial", 1) )
###>>> In this path: cd /home/jolyang/Documents/Github/pvpDiallel
###>>> RUN: sbatch -p med --mem 2600 --ntasks=1 slurm-script/run_gensel_h2_array.sh


### real k for trait bph
phenopwd <- "/home/jolyang/Documents/Github/pvpDiallel/largedata/pheno/CV5fold_BPHmax/"
h2 <- list.files(path="/home/jolyang/Documents/Github/pvpDiallel/largedata/newGERPv2/allgeno_kbph", 
                 pattern="h2.gs.newbin$", full.names=TRUE)

df <- data.frame(asi=h2[grep("asi", h2)], dtp=h2[grep("dtp", h2)], dts=h2[grep("dts", h2)], 
                 eht=h2[grep("eht", h2)], gy=h2[grep("gy", h2)], pht=h2[grep("pht", h2)], 
                 tw=h2[grep("tw", h2)])

inputdf2 <- get_inputdf(mygeno=df, phenopwd)
inputdf2$out <- paste0(gsub(".*/|.txt", "", inputdf2$trainpheno), 
                       gsub(".*gerpv2_b0|.gs.newbin", "", inputdf2$geno))

write.table(inputdf2, "largedata/newGERPv2/inputdf_realk_bph_42000.csv", sep=",", quote=FALSE)

inputdf2 <- read.csv("largedata/newGERPv2/inputdf_realk_bph_42000.csv")
run_GenSel4(inputdf=inputdf2, inpdir="largedata/newGERPv2/allgeno_bph_k", cmdno=500,
            shid = "slurm-script/run_gensel_realk_array.sh",
            email="yangjl0930@gmail.com", runinfo = c(TRUE, "serial", 1) )
###>>> In this path: cd /home/jolyang/Documents/Github/pvpDiallel
###>>> RUN: sbatch -p med --mem 2600 --ntasks=1 slurm-script/run_gensel_h2_array.sh






