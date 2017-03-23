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