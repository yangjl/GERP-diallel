#!/bin/bash
#SBATCH -D /home/jolyang/Documents/Github/GERP-diallel
#SBATCH -o /home/jolyang/Documents/Github/GERP-diallel/slurm-log/testout-%j.txt
#SBATCH -e /home/jolyang/Documents/Github/GERP-diallel/slurm-log/err-%j.txt
#SBATCH -J upload
#SBATCH --mail-user=yangjl0930@gmail.com
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL #email if fails
set -e
set -u

ascp -i ~/bin/aspera.openssh -QT -l100m -k1 -d /group/jrigrp4/diallel_fq/fastq subasp@upload.ncbi.nlm.nih.gov:uploads/yangjl0930@gmail.com_7MfdNojq/GERP-diallel/
