### Jinliang Yang
### plot the relationship of gerp and h


get_k <- function(files){
  
  geno <- read.csv("largedata/GERPv2/gerpsnp_506898.csv")
  geno <- geno[, 1:5]
  
  output <- list()
  for(i in 1:length(files)){
    h1 <- read.table(files[i], header=TRUE)
    names(h1) <- c("snpid","chr","pos","Effect_A","Effect_D","Effect_A2","Effect_D2","h2_mrk_A", 
                   "h2_mrk_D","H2_mrk","h2_mrk_A_p","h2_mrk_D_p","H2_mrk_p","log10_h2_mrk_A","log10_h2_mrk_D","log10_H2_mrk")
    
    h1$k <- h1$Effect_D/h1$Effect_A
    
    trait <- gsub(".*/|_.*", "", files[i])
    output[[trait]] <- merge(geno[, 1:3], h1,  by.y="snpid", by.x="marker")
    gh[[i]] <- 
    message(sprintf("###>>> processing file: [ %s ]", files[i]))
  }
  return(output)
}

cor_k_others <- function(k, RScutoff=0){
  out <- data.frame()
  for(i in 1:7){
    k[[i]] <- subset(k[[i]], RS > RScutoff)
    tst1 <- cor.test(k[[i]]$Effect_A, k[[i]]$RS)
    tst2 <- cor.test(k[[i]]$Effect_D, k[[i]]$RS)
    tst3 <- cor.test(k[[i]]$k, k[[i]]$RS)
    tst4 <- cor.test(k[[i]]$h2_mrk_A, k[[i]]$RS)
    tst5 <- cor.test(k[[i]]$h2_mrk_D, k[[i]]$RS)
    tst6 <- cor.test(k[[i]]$H2_mrk, k[[i]]$RS)
    
    
    out1 <- data.frame(trait=names(k)[i], test="A_RS", p=tst1$p.value, cor=tst1$estimate)
    out2 <- data.frame(trait=names(k)[i], test="D_RS", p=tst2$p.value, cor=tst2$estimate)
    out3 <- data.frame(trait=names(k)[i], test="K_RS", p=tst3$p.value, cor=tst3$estimate)
    out4 <- data.frame(trait=names(k)[i], test="h2_A_RS", p=tst4$p.value, cor=tst4$estimate)
    out5 <- data.frame(trait=names(k)[i], test="h2_D_RS", p=tst5$p.value, cor=tst5$estimate)
    out6 <- data.frame(trait=names(k)[i], test="h2_RS", p=tst6$p.value, cor=tst6$estimate)
    
    out <- rbind(out, out1, out2, out3, out4, out5, out6)
    #print()
  }
  return(out)
}

####
files <- list.files(path="largedata/snpeff", pattern="snpe$", full.names=TRUE)
k <- get_k(files)

kcor0 <- cor_k_others(k, RScutoff=0)
kcor1 <- cor_k_others(k, RScutoff=1)
kcor2 <- cor_k_others(k, RScutoff=2)
# RS bigger not better

save(file="largedata/lcache/k_gerp.RData", list=c("k", "kcor0", "kcor1", "kcor2"))

