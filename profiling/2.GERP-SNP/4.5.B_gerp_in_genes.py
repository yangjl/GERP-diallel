### Jinliang Yang
### ipython
### Feb 9th, 2015

import pandas as pd
import numpy as np

### subset of the gerp
gerp = pd.read_csv("largedata/GERPv2/gerp130m.csv")
len(gerp) #130896913
gerp = gerp[gerp['RS']>0]
len(gerp) #86006888

######################## get v2 gene annotation #################
fgsv2 = pd.read_csv("/home/jolyang/dbcenter/AGP/AGPv2/ZmB73_5b_FGS.gff", sep="\t", header=None)
fgsv2.columns = ["seqname", "source", "feature", "start", "end", "score",
                  "strand", "frame", "attribute"]
gene = fgsv2[fgsv2['feature'] == 'gene']


######## loop through the gene
for index, row in gene.iterrows():
  
  tem = gerp[(gerp['chr'] == row["seqname"]) & (gerp["pos"].astype(int) >= row["start"].astype(int)) & \
  (gerp["pos"].astype(int) <= row["end"].astype(int))]
  
  print index, row

 tem = gerp[(gerp['chr'].astype(str) == row["seqname"].astype(str))]
