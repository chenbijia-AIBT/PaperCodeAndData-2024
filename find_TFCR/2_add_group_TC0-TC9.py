import pandas as pd
 
data = pd.read_csv("TFCR_annotation.txt",sep="\t")
new_data = data.sort_values("X31.889.1",ascending=True) ###选择data的第七列进行排序，即TFCR的复杂度
TFCR = new_data.values.tolist()

#根据输入文件的行数设定分组数（目的是将该文件按照score/复杂度尽量平均分成10组），K值是文件长度/10
K=12154      
for i in range(len(TFCR)):
    if i <= K-1:
        TFCR[i].append("TC0")
    elif i <= 2*K-1:
          TFCR[i].append("TC1")
    elif i <= 3*K-1:
          TFCR[i].append("TC2")
    elif i <= 4*K-1:
          TFCR[i].append("TC3")
    elif i <= 5*K-1:
          TFCR[i].append("TC4")  
    elif i <= 6*K-1:
          TFCR[i].append("TC5")
    elif i <= 7*K:
          TFCR[i].append("TC6")
    elif i <= 8*K-1:
          TFCR[i].append("TC7")
    elif i <= 9*K-1:
          TFCR[i].append("TC8")
    elif i > 9*K-1:
          TFCR[i].append("TC9")
          
result = sorted(TFCR)
result1 = pd.DataFrame(result)
index_a = data.columns
index_b = index_a.insert(19, "group")
result1.columns = index_b
result1.to_csv("TFCR_annotation_group.txt",index=False,header=True,sep="\t")
