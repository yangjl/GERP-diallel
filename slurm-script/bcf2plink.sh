#!/bin/bash
#SBATCH -D /home/jolyang/Documents/Github/GERP-diallel
#SBATCH -o /home/jolyang/Documents/Github/GERP-diallel/slurm-log/testout-%j.txt
#SBATCH -e /home/jolyang/Documents/Github/GERP-diallel/slurm-log/err-%j.txt
#SBATCH -J bcftools
#SBATCH --mail-user=yangjl0930@gmail.com
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL #email if fails
set -e
set -u

cd /home/jolyang/dbcenter/AllZeaGBS/
tabix -p vcf AllZeaGBSv2.7_publicSamples_imputedV3b_agpv3_sorted.vcf.gz
bcftools view -S /home/jolyang/Documents/Github/GERP-diallel/cache/id11.txt --no-update -m2 -M2 -v snps AllZeaGBSv2.7_publicSamples_imputedV3b_agpv3_sorted.vcf.gz -Oz -o AllZeaGBSv2.7_elite11_imputedV3b.vcf.gz
