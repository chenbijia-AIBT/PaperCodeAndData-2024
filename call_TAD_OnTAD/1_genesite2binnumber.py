import numpy as np
import pandas as pd
import os
def savetxt(filename,x):
    np.savetxt(filename,x,delimiter = '\t',fmt='%s') 

dict1 = []
for i in range(1,23):
    dict1.append("chr"+str(i))
dict1.append("chrX")

dict2 = ['GM12878','HMEC','HUVEC','IMR90','K562','NHEK']
res = ["5","10","25","50"]

for i in dict2:
    for r in res:
        for item in dict1:
            f = pd.read_csv("0.GSE63525_data/" + i + "/"+r+"kb_resolution_intrachromosomal/" + item + "/MAPQGE30/"+ item +"_"+r+"kb.RAWobserved", sep = "\t", header=None)

            f[0] = f[0]/(int(r)*1000)   
            f[1] = f[1]/(int(r)*1000)

            f[0] = f[0].astype(int)
            f[1] = f[1].astype(int)
            f[2] = f[2].astype(int)
        
            isExists=os.path.exists("Input/" + i + "/raw/sparse/"+r+"000")
            if not isExists:
                os.makedirs("Input/" + i + "/raw/sparse/"+r+"000")
        
            savetxt("Input/" + i + "/raw/sparse/"+r+"000/" + item + "_raw_"+r+"kb_sparse.txt", f)










