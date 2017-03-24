### Jinliang Yang
### Sept 4th, 2015
### updated: 03-23-2017
### impute genotype for hybrids



################
geno_for_k <- function(geno, ped, RS_cutoff=0, outfile="test"){
  
  geno <- subset(geno, RS > RS_cutoff)
  message(sprintf("###>>> [ %s ] SNPs for computing!", nrow(geno)))
  res <- geno[, 1:5]
  for(i in 1:nrow(ped)){
    p1 <- ped$P1[i]
    p2 <- ped$P2[i]
    # note: need to get the reference allele!
    geno$tem1 <- 3
    geno[geno[, p1] == geno[, "B73"] & geno[, "B73"] != "N", ]$tem1 <- 1
    if(nrow(geno[geno[, p1] != geno[, "B73"] & geno[, p1]!="N" & geno[, "B73"] != "N", ]) > 0){
      geno[geno[, p1] != geno[, "B73"] & geno[, p1]!="N" & geno[, "B73"] != "N", ]$tem1 <- 0
    }
    
    geno$tem2 <- 3
    geno[geno[, p2] == geno[, "B73"] & geno[, "B73"] != "N", ]$tem2 <- 1
    geno[geno[, p2] != geno[, "B73"] & geno[, p2]!="N" & geno[, "B73"] != "N", ]$tem2 <- 0
    
    res$out <- geno$tem1 + geno$tem2
    names(res)[ncol(res)] <- ped$Genotype[i]
    message(sprintf("###>>> impute genotype for [ %sth ] hybrid [ %s ]", i, ped$Genotype[i]))
  }
  
  #### output chr by chr
  for(chri in 1:10){
    chr <- subset(res, chr==chri)
    chrmx <- chr[, -1:-5]
    chrmx[chrmx>3] <- 3
    tchr <- as.data.frame(t(chrmx))
    names(tchr) <- chr$marker
    outchr <- cbind(data.frame(Sample_ID=ped$Genotype), tchr)
    write.table(outchr, paste0(outfile, chri, ".txt"), sep="\t", 
                row.names=FALSE, quote=FALSE)
  }
  message(sprintf("###>>> DONE!")) 
} 


### get the parental genotype info
geno <- read.csv("largedata/gerpsnp_506898.csv")
v23 <- read.csv("largedata/gerp23.csv")

geno <- merge(v23[, c("marker", "agpv3", "RS")], geno[, -3], by="marker")
write.table(geno, "largedata/gerpsnp_v3_410183.csv", sep=",", row.names=FALSE, quote=FALSE)


### get the pedigree info
ped <- read.table("largedata/pheno/wholeset/asi_perse.txt", header=TRUE)
ped$P1 <- gsub("x.*", "", ped$Genotype)
ped$P2 <- gsub(".*x", "", ped$Genotype)

### output map file
map <- geno[, c("marker", "chr", "pos")]
map <- map[order(map$chr, map$pos), ]
write.table(map, "largedata/SNP/genotype_h.map", sep="\t", row.names=FALSE, quote=FALSE)

####
geno_for_k(geno, ped, RS_cutoff=0, outfile="largedata/SNP/genotype_h")

