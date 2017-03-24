### Jinliang Yang
### 2016-08-22

library("farmeR")

set_array_job(shid="slurm-script/getmajor.sh",
             shcode='R --no-save --args ${SLURM_ARRAY_TASK_ID} < profiling/2.SNP/2.G.4_hmp3_cal_delete.R',
             arrayjobs="1-9", wd=NULL, jobid="getmajor", 
             email="yangjl0930@gmail.com",
             run=c(FALSE, "bigmemm", 2, "16G"))

###>>> In this path: cd /home/jolyang/Documents/Github/pvpDiallel
###>>> RUN: sbatch -p bigmemh --mem 8G --ntasks=1 --time=12:00:00 slurm-script/getmajor.sh




library("data.table")

df <- data.frame()
for(i in 1:10){
    t <- fread(paste0("largedata/Alignment/hmp3_GERPv3_major_chr", i, ".csv"), data.table=FALSE)
    df <- rbind(df, t)
}

gsnpv3 <- fread("largedata/gerpsnp_v3_410183.csv", data.table=FALSE)
df2 <- merge(df[, c("snpid", "major")], gsnpv3[, -7:-8], by.x="snpid", by.y="agpv3")
names(df2)[2] <- "del"
write.table(df2, "largedata/gerpsnp_v3_345176_del.csv", sep=",", row.names=FALSE, quote=FALSE)

sum(df2$del == df2$B73)
#[1] 299838
sum(df2$del != df2$B73)
#[1] 45338
