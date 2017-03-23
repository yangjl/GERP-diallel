# Jinliang Yang
# updated: 12/12/2014
# Purpose: select SNPs from GERP conserved region


GERPelem2bed6 <- function(indf=gerpelemt, chr="chr1"){

  elem <- read.table(infile, header=FALSE)
  names(elem) <- c("chrom", "start",  "end",  "score",  "strand", "length", "name", "V8")
  elem <- elem[, c("chrom", "start", "end", "name", "score", "strand")]
  elem$chrom <- chr;
  #elem$start <- elem$start -1;
  elem$name <- paste(elem$chrom, elem$start, sep="_")
  return(elem)
}



gerp <- read.table("largedata/GERPv2/gerpelemt.bed6", header=FALSE)

### GERP command
bedtools intersect -a largedata/GERPv2/gerpelemt.bed6 -b largedata/SNP/allsnps_11m.bed3 
> largedata/SNP/allsnps_11m_gerpelemt.bed

