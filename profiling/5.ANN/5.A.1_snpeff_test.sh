
### Jinliang Yang
### 04/13/2017

### snpEff

# step1: build a genome for it
java -jar snpEff.jar build -gff3 AGPv3 -v > err.txt


# step2:

java -Xmx32g -jar snpEff.jar AGPv3 ~/dbcenter/HapMap/HapMap3/chr10_bisnp_n35.vcf.gz -csvStats > test.ann.vcf

java -Xmx32g -jar snpEff.jar AGPv3 -i bed /home/jolyang/Documents/Github/GERP-diallel/largedata/SNP/gerpsnp_v3.bed -csvStats bed_sum.csv > test.ann.bed


java -Xmx4g -jar snpEff.jar Zea_mays examples/test.chr22.vcf > test.chr22.ann.vcf
