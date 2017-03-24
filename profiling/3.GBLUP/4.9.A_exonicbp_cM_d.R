### Jinliang Yang
### August/25/2016
### (exonic bp/cM) vs d


get_exonicbp_cM <- function(gff_file="~/dbcenter/AGP/AGPv2/ZmB73_5b_FGS_exon.gff"){
  gff <- read.table(gff_file, header=FALSE)
  gff$tx <- gsub("Parent=|;.*", "", gff$V9)
  gff$geneid <- gsub("_.*", "", gff$tx)
  message(sprintf("### Processing [ %s ] transcripts for [ %s ] genes", 
                  length(unique(gff$tx)), length(unique(gff$geneid))))
  gff$exonlen <- gff$V5 - gff$V4
  
  ###### GERP physical and cM
  gerp <- read.csv("cache/gerpsnp_506898_gp.csv")
  gerp$bingen <- paste(gerp$chr, round(gerp$genetic, 0), sep="_")
  
  out <- data.frame()
  for(gbin in unique(gerp$bingen)){
    chunk <- subset(gerp, bingen %in% gbin)
    chri <- gsub("_.*", "", gbin)
    pos <- range(chunk$pos)
    
    bp <- sum(subset(gff, V1 %in% chri & V4 > pos[1] & V4 < pos[2])$exonlen)
    
    tem <- data.frame(chr=chri, start=pos[1], end=pos[2], cM=gbin, exonbp = bp)
    message(sprintf("### processing [ %s ] ...", gbin))
    out <- rbind(out, tem)
  }
  
  out$genetic <- as.numeric( gsub(".*_", "", out$cM))
  out$midpos <- out$start + out$exonbp
  return(out)
}

########
out <- get_exonicbp_cM(gff_file="~/dbcenter/AGP/AGPv2/ZmB73_5b_FGS_exon.gff")
write.table(out, "cache/exonic_cM.csv", sep=",", row.names=FALSE, quote=FALSE)
out <- read.csv("cache/exonic_cM.csv")


######## k values 5x variance explained!
res2 <- read.csv("cache/kval_perse_5x.csv")
gerp <- read.csv("cache/gerpsnp_506898_gp.csv")
gerp$bingen <- paste(gerp$chr, round(gerp$genetic, 0), sep="_")

res3 <- merge(res2, gerp, by.x="snpid", by.y="marker")

res4 <- merge(res3, out, by.x="bingen", by.y="cM")


t <- subset(res4, trait == "GY")
r <- cor(t$exonbp, abs(t$Effect_D))
plot(t$exonbp, abs(t$Effect_D), main= paste("GY cor=", round(r,5)))






#################
#cols <- wes_palette(7, name = "Zissou", type = "continuous")
cols <- c("#f6546a", "#daa520", "#00ff00", "#66cdaa", "#3b5998", "#8a2be2", "#ff00ff")
theme_set(theme_grey(base_size = 18)) 

#lty1 <- getlty(df=out, eff="effa", cutoff=0.05)$l
p1 <- ggplot(res4, aes(x=round(exonbp/1000,0), y=abs(Effect_D))) +
  facet_grid(~trait, scales = "free") +
  theme_bw() +
  xlab("Exonic Kb") +
  ylab("abs(dominance effect)") +
  guides(colour=FALSE, linetype=FALSE) +
  
  geom_point() +
  theme(axis.text.y = element_text(angle = 90, hjust = 1))
p1

#lty1 <- getlty(df=out, eff="effa", cutoff=0.05)$l
p2 <- ggplot(res4, aes(x=round(exonbp/1000,0), y=abs(Effect_A))) +
  facet_grid(~trait, scales = "free") +
  theme_bw() +
  xlab("Exonic Kb") +
  ylab("abs(Additive effect)") +
  guides(colour=FALSE, linetype=FALSE) +
  
  geom_point() +
  theme(axis.text.y = element_text(angle = 90, hjust = 1))
p2
