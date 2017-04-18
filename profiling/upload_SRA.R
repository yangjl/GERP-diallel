
library("farmeR")

cmd <- paste0("ascp -i ~/bin/aspera.openssh -Tr -l100m -k1 ",
              "-d /group/jrigrp4/diallel_fq/fastq ",
              "subasp@upload.ncbi.nlm.nih.gov:uploads/yangjl0930@gmail.com_7MfdNojq/GERP-diallel/")

set_farm_job(slurmsh = "slurm-script/upload_sra.sh",
             shcode = cmd, wd = NULL, jobid = "upload",
             email = "yangjl0930@gmail.com", runinfo = c(FALSE, "bigmemh", "1", "5G", "80:00:00"))

###>>> In this path: cd /home/jolyang/Documents/Github/GERP-diallel
###>>> RUN: sbatch -p bigmemh --mem 5G --ntasks=1 --time=80:00:00 slurm-script/upload_sra.sh

## ssh farm -p 2022