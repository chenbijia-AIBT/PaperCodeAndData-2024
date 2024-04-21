import pandas as pd
import numpy as np
def savetxt(filename,x):
	np.savetxt(filename,x,delimiter = '\t',fmt='%s')
    
f = pd.read_table('./3_scan/fimo.tsv',sep='\t',header=None,low_memory=False)
df = f.drop(columns=[1])
dfDropCol = df.drop([0])
        
savetxt("./4_out_fimo/fimo.txt",dfDropCol)
