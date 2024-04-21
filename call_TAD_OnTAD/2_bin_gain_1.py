import os
import pandas as pd
import numpy as np
def savetxt(filename,x):
    np.savetxt(filename,x,delimiter = '\t',fmt='%s') 

chrname = []
for i in range(1,23):
    chrname.append("chr"+str(i))
chrname.append("chrX")

dirs = ["GM12878","HMEC","HUVEC","IMR90","K562","NHEK"]
res = ["5","10","25","50"]

for d in dirs:
    for r in res:
        for c in chrname:
            f1 = pd.read_csv("Input/"+d+"/raw/sparse/"+r+"000/"+c+"_raw_"+r+"kb_sparse.txt", sep = "\t",header=None)
            f1[0] = f1[0] + 1
            f1[1] = f1[1] + 1

            f1[0] = f1[0].astype(int)
            f1[1] = f1[1].astype(int)
            f1[2] = f1[2].astype(int)
        
            savetxt("Input/"+d+"/raw/sparse/"+r+"000/"+c+"_raw_"+r+"kb_sparse_gain_1.txt", f1)
            
    