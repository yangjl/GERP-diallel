### Jinliang Yang
### beanplot of the final results with all the SNPs


#########################################
##>>>>>
library("beanplot")
add_bean_plot <- function(resdf = HPH, mymode="a2", ...){
  
  res0 <- resdf
  
  res0$cs <- gsub("cs0", "real", res0$cs)
  res0$cs <- gsub("cs.*", "cs", res0$cs)
  res1 <- subset(res0, mode == mymode)
  par(lend = 1, mai = c(0.8, 0.8, 0.5, 0.5))
  
  res1$cs <- factor(res1$cs, levels = c("real", "cs"))
  res1$trait <- toupper(res1$trait)
  res1$trait <- factor(res1$trait, levels = toupper(c("asi", "dtp", "dts", "eht", "gy", "pht", "tw")) )
  beanplot(r ~ cs + trait, data = res1, ll = 0.04, cex=1.5,
           side = "both", border = NA, col = list(c("blue", "red"), c("grey", "black")) , ...)
  #legend("bottomleft", fill = c("black", "grey"),
  #       legend = c("Group 2", "Group 1"))
  
}

##### get the data frames
bph_max <- read.csv("cache/gerpall_BPHmax.csv")
bph_max <- subset(bph_max, trait != "asi")
bph_min <- read.csv("cache/gerpall_BPHmin.csv")

BPH <- rbind(bph_max, bph_min)
mean(subset(BPH, cs == "real" & mode == "a2")$r) #[1] 0.49
mean(subset(BPH, cs == "real" & mode == "d2")$r) #[1] 0.42


pbph_max <- read.csv("cache/gerpall_pBPHmax.csv")
pbph_max <- subset(pbph_max, trait != "asi")
pbph_min <- read.csv("cache/gerpall_pBPHmin.csv")
mean(subset(pBPH, cs == "real" & mode == "a2")$r) #[1] 0.29
mean(subset(pBPH, cs == "real" & mode == "d2")$r) #[1] 0.24

pBPH <- rbind(pbph_max, pbph_min)
perse <- read.csv("cache/gerpall_perse.csv")
mean(subset(perse, cs == "real" & mode == "a2")$r) #[1] 0.81
mean(subset(perse, cs == "real" & mode == "d2")$r) #[1] 0.70

##################################################################################################################

pdf("graph/Figure_gerpall.pdf", width=10, height=8)
par(mfrow=c(2,2))
add_bean_plot(resdf = perse, mymode="a2", 
              main = "Traits per se using additive model", ylab = "Accuracy (r)")
add_bean_plot(resdf = perse, mymode="d2", 
              main = "Traits per se using dominant model", ylab = "Accuracy (r)")

add_bean_plot(resdf = BPH, mymode="a2", 
              main = "BPH using additive model", ylab = "Accuracy (r)")
add_bean_plot(resdf = BPH, mymode="d2", 
              main = "BPH using dominant model", ylab = "Accuracy (r)")

#add_bean_plot(resdf = pHPH, mymode="a2", 
#              main = "pBPH with additive model", ylab = "CV Accuracy (r)")
#add_bean_plot(resdf = pHPH, mymode="d2", 
#              main = "pBPH with dominant model", ylab = "CV Accuracy (r)")

dev.off()


#####
pval <- read.csv("cache/pval_gerpall.csv")
subset(pval, pval <= 0.05)

idx <- grep("MPH", pval$file)
pval2 <- pval[-idx, ]

pval2$FDR <- p.adjust(pval2$pval, method="fdr" )
subset(pval2, FDR < 0.05)

###>>>
pval3 <- pval2
pval3$trait <- toupper(pval3$trait)
pval3$mode <- gsub("a2", "additive", pval3$mode)
pval3$mode <- gsub("d2", "dominant", pval3$mode)
pval3$file <- gsub(".*_", "", pval3$file)
pval3$file <- gsub("\\.csv", "", pval3$file)
names(pval3) <- c("phenotype", "p-value", "mode", "r_real", "r_cs", "trait","FDR")

###>>>
pval4 <- subset(pval3, trait %in% c("perse", "BPHmax", "BPHmin"))
write.table(pval4, "manuscript/Figure_Table/Table_S3_allsnps_FDR.csv", sep=",", )







