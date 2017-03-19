# Jinliang Yang
# Octo 27th, 2014
# purpose: cor plot of 7 phenotypic traits

get_pheno <- function(trait=trait, pheno="valHyb"){
  ti <- c("ASI", "DTS","DTP", "TW", "PHT", "EHT", "GY")
  myp <- subset(trait, trait==ti[1])
  myp <- myp[, c("Hyb", pheno)]
  names(myp)[2] <- ti[1]
  for(i in 2:length(ti)){
    tem <- subset(trait, trait==ti[i])
    tem <- tem[, c("Hyb", pheno)]
    names(tem)[2] <- ti[i]
    myp <- merge(myp, tem, by="Hyb")
  }
  return(myp)
}

#################
source("lib/Correlation_plot.R")
trait <- read.csv("data/hyb_heterosis.csv")


### get percentage of MPH
myp3 <- get_pheno(trait=trait, pheno="valHyb")
pdf("graphs/Fig_S10.pdf", width=7, height=7)
pairs(myp3[, 2:8], text.panel = diag, upper.panel=panel.smooth,
      lower.panel=panel.cor, gap=0, main="", pch=19, col="grey", lwd=2)
dev.off()


cor.test(myp3$DTS, myp3$DTP, method = "spearman")
cor.test(myp3$PHT, myp3$EHT, method = "spearman")
cor.test(myp3$GY, myp3$DTP, method = "spearman")
cor.test(myp3$GY, myp3$DTS, method = "spearman")

#pdf("~/Documents/Heterosis_GWAS/HGWAS_proj/reports/S.F2_cor.pdf", height=8, width=8)
#pairs(pheno[, 4:10], text.panel = diag, upper.panel=panel.smooth, 
#      lower.panel=panel.cor, gap=0, main="", pch=19, col="grey", lwd=2)
#dev.off()
