### Jinliang Yang
### updated: 03-19-2017
### exploring genome-wide feature of GERP score

library("data.table")
readgerp <- function(){
  #### read the table fast
  out <- data.frame()
  for(chri in 1:10){
    gerpfile <- paste0("/group/jrigrp/gerp/GERPv3/roast.chrom.", chri, ".msa.in.rates.full")
    
    chr <- fread(gerpfile, header=FALSE, data.table=FALSE)
    names(chr) <- c("N", "RS")
    chr$chrom <- chri
    chr$pos <- 1:nrow(chr)
    
    
    message(sprintf("Chr [%s]: >0 [ %s ], <0 [ %s ]", chri, 
                    nrow(subset(chr, RS>0)), nrow(subset(chr, RS<0)) ))
    
    chr <- subset(chr, RS > 0)
    
    outfile <- paste0("largedata/GERPv3/AGPv3_gerp_chr", chri, ".csv")
    write.table(chr, outfile, sep=",", row.names=FALSE, quote=FALSE)
    
  }
}

####-----------------
allgerp <- readgerp()
# Read 301433381 rows and 2 (of 2) columns from 1.449 GB file in 00:02:11
# Chr [1]: >0 [ 35672919 ], <0 [ 15361221 ]
# Read 237865861 rows and 2 (of 2) columns from 1.136 GB file in 00:00:50
# Chr [2]: >0 [ 27149438 ], <0 [ 11830592 ]
# Read 232221667 rows and 2 (of 2) columns from 1.099 GB file in 00:00:44
# Chr [3]: >0 [ 25265529 ], <0 [ 11159023 ]
# Read 242024971 rows and 2 (of 2) columns from 1.135 GB file in 00:01:09
# Chr [4]: >0 [ 24828834 ], <0 [ 11350803 ]
# Read 217906509 rows and 2 (of 2) columns from 1.039 GB file in 00:01:13
# Chr [5]: >0 [ 25150515 ], <0 [ 10616256 ]
# Read 169354735 rows and 2 (of 2) columns from 0.804 GB file in 00:00:39
# Chr [6]: >0 [ 18722827 ], <0 [ 8220868 ]
