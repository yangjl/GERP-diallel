### Jinliang Yang
### 9/4/2015
### purpose: using p2g to find the genetic map of Diallel pop


test_p2g <- function(chr=9){
  source("~/Documents/Github/zmSNPtools/Rcodes/p2g.R")
  train <- read.csv("~/Documents/Github/zmSNPtools/shareData/ISU_integrated_IBM_AGPv2_training.csv")
  train$Physical <- as.numeric(as.character(gsub(",", "", train$Physical)))
  ### test the proformance of the function using the training data
  
  predict <- subset(train, Chr == chr)[, -3]
  #predict <- data.frame(marker=c("A", "B", "C"), chr=1, Physical=c(467, 41238467, 164555621)) 
  out <- p2g(predict, train)
  res <- merge(train, out, by.x="Marker", by.y="marker")
  
  message(sprintf("###>>> p2g testing results: correlation R2 = [ %s ]", cor(res$Genetic, res$genetic)))
}

use_p2g <- function(){
  source("~/Documents/Github/zmSNPtools/Rcodes/p2g.R")
  train <- read.csv("~/Documents/Github/zmSNPtools/shareData/ISU_integrated_IBM_AGPv2_training.csv")
  train$Physical <- as.numeric(as.character(gsub(",", "", train$Physical)))
  ### test the proformance of the function using the training data
  
  predict <- subset(train, Chr == chr)[, -3]
  #predict <- data.frame(marker=c("A", "B", "C"), chr=1, Physical=c(467, 41238467, 164555621)) 
  out <- p2g(predict, train)
  res <- merge(train, out, by.x="Marker", by.y="marker")
  
  message(sprintf("###>>> p2g testing results: correlation R2 = [ %s ]", cor(res$Genetic, res$genetic)))
}

