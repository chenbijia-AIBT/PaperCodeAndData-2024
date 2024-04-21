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
        os.system("mkdir -p Input/"+d+"/ICE/sparse/"+r+"000")
        for c in chrname:
            os.system("ice --results_filename Input/"+d+"/ICE/sparse/"+r+"000/"+c+"_ice_"+r+"kb_sparse.txt --filter_low_counts_perc 0.02 --filter_high_counts_perc 0 --max_iter 100 --eps 0.1 --remove-all-zeros-loci --output-bias 1 --verbose 1 Input/"+d+"/raw/sparse/"+r+"000/"+c+"_raw_"+r+"kb_sparse_gain_1.txt")
    