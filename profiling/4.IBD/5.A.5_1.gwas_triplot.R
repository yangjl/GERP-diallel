# Jinliang Yang
# 9/21, 2014
# purpose: get the prdiction results of triploid genomes


pred_triploid <- function(basepwd="slurm-scripts/gwas_b1/gerpIBD_k_b1_"){
  
  ti <- tolower(c("ASI", "DTP", "DTS", "EHT",  "GY", "PHT",  "TW"))
  out <- list()
  for(i in 1:7){
    cgr <- read.table(paste0(basepwd, ti[i], ".cgrRes1"), skip=2)
    
    ghat <- read.table(paste0(basepwd, ti[i], ".ghatREL1"), skip=1)
    #Animal gHat ASI Fix  meanBias# PEV=Var(g/y)#   R^2# 
    ghat$V2 <- cgr$V2 + ghat$V2
    ghat_a2b <- read.table(paste0(basepwd, ti[i], "_predict_a2b.ghat1"), header=TRUE)
    ghat_a2b$gHat <- ghat_a2b$gHat + cgr$V2
    ghat_ab2 <- read.table(paste0(basepwd, ti[i], "_predict_ab2.ghat1"), header=TRUE)
    ghat_ab2$gHat <- ghat_ab2$gHat + cgr$V2
    
    message(sprintf("###>>> trait [ %s ]:", ti[i]))
    
    alld <- merge(ghat[, 1:2], ghat_a2b[, 1:2], by.x="V1", by.y="Animal_ID")
    alld <- merge(alld, ghat_ab2[, 1:2], by.x="V1", by.y="Animal_ID")
    
    alld$max <- apply(alld[, 3:4], 1, max)
    alld$min <- apply(alld[, 3:4], 1, min)
    names(alld) <- c("plantid", "AB", "AAB", "ABB", "max", "min")
    out[[ti[i]]] <- alld
  }
  return(out)
}


k_b1 <- pred_triploid(basepwd="slurm-scripts/gwas_b1/gerpIBD_k_b1_")
k_b2 <- pred_triploid(basepwd="slurm-scripts/gwas_b2/gerpIBD_k_b2_")
k_b0 <- pred_triploid(basepwd="slurm-scripts/gwas_b0/gerpIBD_k_")
save(file="cache/triploid.RData", list=c("k_b1", "k_b2", "k_b0"))

###################################################
ob <- load("cache/triploid.RData")

plot_tri <- function(kdat=k_b2, mytrait="GY", getmax=TRUE, ...){
  
  trait <- read.csv("data/trait_matrix.csv")
  trait$plantid <- paste(trait$P1, trait$P2, sep="x")
  dat <- kdat[[tolower(mytrait)]]
  
  myt <- subset(trait, trait == toupper(ti[j]))
  dat <- merge(dat, myt[, c("plantid", "valP1", "valP2")])
  
  if(getmax){
    dat$BP <- apply(dat[, 7:8], 1, max)
  }else{
    dat$BP <- apply(dat[, 7:8], 1, min)
  }
  
  
  val <- c(dat$AB, dat$AAB, dat$ABB, dat$BP)
  plot(c(range(val)), c(1,66), type="n", 
       xlab="Breeding Value", ylab="", yaxt="n", bty="n", ...)
  
  #axis(side=2, tick =FALSE, las=1, at=c(105, 95, 85, 75, 65, 55, 45, 35, 25, 15), 
  #     labels=paste("Chr", 1:10, sep=""))
  #### chromosome
  dat <- dat[order(dat$AB),]
  for (i in 1:nrow(dat)){
    bv <- c(dat$AB[i], dat$AAB[i], dat$ABB[i], dat$BP[i])
    lines(c(range(bv)), c(67-i, 67-i), lwd=2, col="grey")
    points(dat$AB[i], 67-i, pch=16, col="red")
    points(dat$AAB[i], 67-i, pch=15, col="blue")
    points(dat$ABB[i], 67-i, pch=17, col="darkgreen")
    points(dat$BP[i], 67-i, pch=18, col="black")
    #lines (c(centromere[i,]$Start,centromere[i,]$End),
    #       c(105-10*i, 105-10*i),lwd=5, col="tomato")
  } 
  print(t.test(dat$AB, dat$max, paired=TRUE))
  print(t.test(dat$AB, dat$min, paired=TRUE))
  print(t.test(dat$max, dat$min, paired=TRUE))
  #t.test(dat$AB, dat$min, paired=TRUE)
}

ti <- tolower(c("ASI", "DTP", "DTS", "EHT",  "GY", "PHT",  "TW"))

pdf("graphs/Fig5_triploid.pdf", width=5, height=8)
plot_tri(kdat=k_b0, mytrait="GY", main="Grain Yield")
dev.off()

pdf("graphs/Supp_Fig5_triploid.pdf", width=5, height=8)
plot_tri(kdat=k_b0, mytrait="ASI",getmax=FALSE, main="Anthesis-Silking Interval")
plot_tri(kdat=k_b0, mytrait="DTP", getmax=TRUE,main="Days to 50% Pollen Shed")
plot_tri(kdat=k_b0, mytrait="DTS", getmax=TRUE,main="Days to 50% Silking")
plot_tri(kdat=k_b0, mytrait="EHT", getmax=TRUE,main="Ear Height")
plot_tri(kdat=k_b0, mytrait="PHT",getmax=TRUE, main="Plant Height")
plot_tri(kdat=k_b0, mytrait="TW",getmax=TRUE, main="Test Weight")
dev.off()


#######
ab <- density(dat$AB)
aab <- density(dat$AAB)
abb <- density(dat$ABB)

plot(aab, xlim=c(150, 200))
lines(abb)

