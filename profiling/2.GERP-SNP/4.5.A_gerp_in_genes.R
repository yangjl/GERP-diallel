### Jinliang Yang
### Feb 9th, 2015

library(data.table)
gerp0 <- fread("largedata/GERPv2/gerp130m.csv", sep=",")
#[1] 130896913         4
gerp <- subset(gerp0, RS > 0)
#[1] 86006888        4

gerp$start <- gerp$pos -1
sub <- subset(gerp, select = c("chr", "start", "pos", "RS"))

options(scipen=999)
write.table(sub, "~/Documents/Github/Misc/largedata/1.gc/gerpv2.bed4", sep="\t",
            row.names=FALSE, col.names=FALSE, quote=FALSE)


sub2 <- fread("~/Documents/Github/Misc/largedata/1.gc/gerpv2.bed4", header=FALSE)

######################## Filtering the duplicated annotation #################
### Note: transfter to the Misc project
################################################################################

### using the following to sumarise, by supper slow
sapply(1:10, 
       function(i){
         tem <- subset(gerp, chr == gene$seqname[i] & pos >= gene$start[i] & pos <= gene$end[i])
         #res <- data.frame(geneid = gene$attribute[i], gerpsum= sum(tem))
         message(sprintf("###>>> [ %s ] genes", i))
         return(sum(tem$RS));
       })
