## Jinliang Yang
## Oct. 13th, 2014
## phenotypic data of offpvp diallel

#setwd("~/Documents/Github/pvpDiallel/")
trait5fold <- function(times=10, typep="valHyb", outdir="largedata/pheno/CV5fold/"){
  
  trait <- read.csv("data/trait_matrix_updated_BPH.csv")
  ### => Table_S1.trait_matrix.csv
  trait$Hyb <- paste(trait$P1, trait$P2, sep="x")
  ti <- c("ASI", "DTP", "DTS", "EHT",  "GY", "PHT",  "TW")
  parent <- unique(c(as.character(trait$P1), as.character(trait$P2)))
  
  dir.create(outdir, showWarnings = FALSE)
  
  #### sample number of randomization
  for(j in 1:times){
    #### split the sample into 5-fold
    idx <- sample(1:66, 66, replace=FALSE)
    idx1 <- idx[1:13]
    idx2 <- idx[14:26]
    idx3 <- idx[27:39]
    idx4 <- idx[40:52]
    idx5 <- idx[53:66]
    
    for(i in 1:7){
      pheno <- subset(trait, trait==ti[i])
      pheno <- pheno[, c("Hyb", typep, "valP1")]
      names(pheno) <- c("Genotype", ti[i], "Fix")
      pheno$Fix <- 1
      
      #### 1 training and validation phenotype
      out_t1 <- paste(outdir, tolower(ti[i]), "_train1_sp", j, ".txt", sep="")
      write.table(pheno[c(idx2, idx3, idx4, idx5), ], out_t1, sep="\t", row.names=FALSE, quote=FALSE) 
      out_v1 <- paste(outdir, tolower(ti[i]), "_test1_sp", j, ".txt", sep="")
      write.table(pheno[idx1, ], out_v1, sep="\t", row.names=FALSE, quote=FALSE)
      #### 2 training and validation phenotype
      out_t2 <- paste(outdir, tolower(ti[i]), "_train2_sp", j, ".txt", sep="")
      write.table(pheno[c(idx1, idx3, idx4, idx5), ], out_t2, sep="\t", row.names=FALSE, quote=FALSE) 
      out_v2 <- paste(outdir, tolower(ti[i]), "_test2_sp", j, ".txt", sep="")
      write.table(pheno[idx2, ], out_v2, sep="\t", row.names=FALSE, quote=FALSE)
      #### 3 training and validation phenotype
      out_t3 <- paste(outdir, tolower(ti[i]), "_train3_sp", j, ".txt", sep="")
      write.table(pheno[c(idx1, idx2, idx4, idx5), ], out_t3, sep="\t", row.names=FALSE, quote=FALSE) 
      out_v3 <- paste(outdir, tolower(ti[i]), "_test3_sp", j, ".txt", sep="")
      write.table(pheno[idx3, ], out_v3, sep="\t", row.names=FALSE, quote=FALSE)
      #### 4 training and validation phenotype
      out_t4 <- paste(outdir, tolower(ti[i]), "_train4_sp", j, ".txt", sep="")
      write.table(pheno[c(idx1, idx2, idx3, idx5), ], out_t4, sep="\t", row.names=FALSE, quote=FALSE) 
      out_v4 <- paste(outdir, tolower(ti[i]), "_test4_sp", j, ".txt", sep="")
      write.table(pheno[idx4, ], out_v4, sep="\t", row.names=FALSE, quote=FALSE)
      #### 5 training and validation phenotype
      out_t5 <- paste(outdir, tolower(ti[i]), "_train5_sp", j, ".txt", sep="")
      write.table(pheno[c(idx1, idx2, idx3, idx4), ], out_t5, sep="\t", row.names=FALSE, quote=FALSE) 
      out_v5 <- paste(outdir, tolower(ti[i]), "_test5_sp", j, ".txt", sep="")
      write.table(pheno[idx5, ], out_v5, sep="\t", row.names=FALSE, quote=FALSE)
    } 
    message(sprintf(">>> output random set [%s] !", j))
  }  
}

#### trait per se
set.seed(12345)
# total file 7 traits x 5fold x 10 sample =350 x 2(test file and validation file)
trait5fold(times=100, typep="valHyb", outdir="largedata/pheno/CV5fold/")

#### BPHmax
set.seed(12345)
# total file 7 traits x 5fold x 10 sample =350 x 2(test file and validation file)
trait5fold(times=100, typep="BPHmax", outdir="largedata/pheno/CV5fold_BPHmax/") 

#### pBPHmax
set.seed(12345)
# total file 7 traits x 5fold x 10 sample =350 x 2(test file and validation file)
trait5fold(times=100, typep="pBPHmax", outdir="largedata/pheno/CV5fold_pBPHmax/") 

#### BPHmin
set.seed(12345)
# total file 7 traits x 5fold x 10 sample =350 x 2(test file and validation file)
trait5fold(times=100, typep="BPHmin", outdir="largedata/pheno/CV5fold_BPHmin/")

#### pBPHmin
set.seed(12345)
# total file 7 traits x 5fold x 10 sample =350 x 2(test file and validation file)
trait5fold(times=100, typep="pBPHmin", outdir="largedata/pheno/CV5fold_pBPHmin/") 

#### MPH
set.seed(12345)
# total file 7 traits x 5fold x 10 sample =350 x 2(test file and validation file)
trait5fold(times=100, typep="pBPHmax", outdir="largedata/pheno/CV5fold_MPH/") 

#### pMPH
set.seed(12345)
# total file 7 traits x 5fold x 10 sample =350 x 2(test file and validation file)
trait5fold(times=100, typep="pBPHmax", outdir="largedata/pheno/CV5fold_pMPH/") 

#### pBPH
set.seed(12345)
# total file 7 traits x 5fold x 10 sample =350 x 2(test file and validation file)
trait5fold(times=100, typep="pBPH", outdir="largedata/pheno/CV5fold_pBPH/") 

#### BPH
set.seed(12345)
# total file 7 traits x 5fold x 10 sample =350 x 2(test file and validation file)
trait5fold(times=100, typep="BPH", outdir="largedata/pheno/CV5fold_BPH/") 







