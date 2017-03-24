### Jinliang Yang
### Sept 4th, 2016
### get degree of dominance from random SNP


################
geno_for_k <- function(geno, outfile="test"){
  
  ### get the pedigree info
  ped <- read.table("largedata/pheno/wholeset/asi_perse.txt", header=TRUE)
  ped$P1 <- gsub("x.*", "", ped$Genotype)
  ped$P2 <- gsub(".*x", "", ped$Genotype)
  
  #geno <- subset(geno, RS > RS_cutoff)
  message(sprintf("###>>> [ %s ] SNPs for computing!", nrow(geno)))
  res <- geno[, c("snpid", "chr", "pos", "frq", "qt0", "genetic", "exonbp")]
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
    chrmx <- chr[, -1:-7]
    chrmx[chrmx>3] <- 3
    tchr <- as.data.frame(t(chrmx))
    names(tchr) <- chr$snpid
    outchr <- cbind(data.frame(Sample_ID=ped$Genotype), tchr)
    write.table(outchr, paste0(outfile, chri, ".txt"), sep="\t", 
                row.names=FALSE, quote=FALSE)
  }
  message(sprintf("###>>> DONE!")) 
} 

### get the parental genotype info
for(i in 0:10){
  geno <- read.csv(paste0("largedata/SNP/randomsnp/rsnp", i, ".csv"))
  
  ### output map file
  map <- geno[, c("snpid", "chr", "pos")]
  map <- map[order(map$chr, map$pos), ]
  write.table(map, paste0("largedata/SNP/randomsnp/rsnp", i, ".map"), sep="\t", 
              row.names=FALSE, quote=FALSE)
  ### run the function
  geno_for_k(geno, outfile= paste0("largedata/SNP/randomsnp/rsnp", i, "_chr"))
}

