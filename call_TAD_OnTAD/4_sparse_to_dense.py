import numpy as np
import os 

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
        os.system("mkdir -p Input/"+d+"/ICE/dense/"+r+"000")
        for c in chrname:
            os.system("python sparseToDense.py -b /mnt/f/project1/12_revise/TAD_calling/0.contact_matrices/hic_bed/"+r+"kb/"+c+"_"+r+"kb.bed -o Input/"+d+"/ICE/dense/"+r+"000/"+c+"_ice_"+r+"kb_dense.matrix Input/"+d+"/ICE/sparse/"+r+"000/"+c+"_ice_"+r+"kb_sparse.txt")


