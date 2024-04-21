import numpy as np
from collections import Counter
import os

def savetxt(filename,x):
    np.savetxt(filename,x,delimiter = '\t',fmt='%s')
    
chrname = []
for i in range(1,23):
    chrname.append(str(i))
chrname.append("X")

dirs = ["GM12878","HMEC","HUVEC","IMR90","K562","NHEK"]
for s in dirs: 
    
################5kb#####################################
    os.system("mkdir -p 3.TAD_boundary/"+s+"/5000")
    result = []
    for chrom in chrname:
        left = []
        right = []
        f = open("1.OnTAD_output/"+s+"/5000/OnTAD_ICE_pen0.1_max_400_min_6_chr"+chrom+".tad")
        lines=f.readlines()
        for i in range(1,len(lines)):
            line = lines[i].strip()
            sep_line = line.split("\t")
            left.append(sep_line[0])
            right.append(sep_line[1])
        f.close()
        left_count = Counter(left)
        right_count = Counter(right)
        left_name = list(left_count)
        right_name = list(right_count)
    
        boundary = []
        for j in left_name:
            if j not in right_name:
                boundary.append([int(j),int(left_count[j])])
            else:
                if int(left_count[j])>int(right_count[j]):
                    boundary.append([int(j),int(left_count[j])])
                else:
                    boundary.append([int(j),int(right_count[j])])
                    
        for k in right_name:
            if k not in left_name:
                boundary.append([int(k),int(right_count[k])])
                
        for h in range(len(boundary)):
            if int(boundary[h][1]) < 3: 
                result.append(["chr"+chrom,(int(boundary[h][0])-2)*5000,(int(boundary[h][0])-1)*5000,int(boundary[h][1])])
            else:
                result.append(["chr"+chrom,(int(boundary[h][0])-2)*5000,(int(boundary[h][0])-1)*5000,"3+"])
                
    savetxt("3.TAD_boundary/"+s+"/5000/boundary_level_3",result)

####################10kb#####################################
    os.system("mkdir -p 3.TAD_boundary/"+s+"/10000")
    result = []
    for chrom in chrname:
        left = []
        right = []
        f = open("1.OnTAD_output/"+s+"/10000/OnTAD_ICE_pen0.1_max_200_min_3_chr"+chrom+".tad")
        lines=f.readlines()
        for i in range(1,len(lines)):
            line = lines[i].strip()
            sep_line = line.split("\t")
            left.append(sep_line[0])
            right.append(sep_line[1])
        f.close()
        left_count = Counter(left)
        right_count = Counter(right)
        left_name = list(left_count)
        right_name = list(right_count)
        
        boundary = []
        for j in left_name:
            if j not in right_name:
                boundary.append([int(j),int(left_count[j])])
            else:
                if int(left_count[j])>int(right_count[j]):
                    boundary.append([int(j),int(left_count[j])])
                else:
                    boundary.append([int(j),int(right_count[j])])
    
        for k in right_name:
            if k not in left_name:
                boundary.append([int(k),int(right_count[k])])
                
        for h in range(len(boundary)):
            if int(boundary[h][1]) < 3: 
                result.append(["chr"+chrom,(int(boundary[h][0])-2)*10000,(int(boundary[h][0])-1)*10000,int(boundary[h][1])])
            else:
                result.append(["chr"+chrom,(int(boundary[h][0])-2)*10000,(int(boundary[h][0])-1)*10000,"3+"])
                
    savetxt("3.TAD_boundary/"+s+"/10000/boundary_level_3",result)
    
####################25kb#####################################
    os.system("mkdir -p 3.TAD_boundary/"+s+"/25000")
    result = []
    for chrom in chrname:
        left = []
        right = []
        f = open("1.OnTAD_output/"+s+"/25000/OnTAD_ICE_pen0.1_max_80_min_1.2_chr"+chrom+".tad")
        lines=f.readlines()
        for i in range(1,len(lines)):
            line = lines[i].strip()
            sep_line = line.split("\t")
            left.append(sep_line[0])
            right.append(sep_line[1])
        f.close()
        left_count = Counter(left)
        right_count = Counter(right)
        left_name = list(left_count)
        right_name = list(right_count)
        
        boundary = []
        for j in left_name:
            if j not in right_name:
                boundary.append([int(j),int(left_count[j])])
            else:
                if int(left_count[j])>int(right_count[j]):
                    boundary.append([int(j),int(left_count[j])])
                else:
                    boundary.append([int(j),int(right_count[j])])
                    
        for k in right_name:
            if k not in left_name:
                boundary.append([int(k),int(right_count[k])])
                
        for h in range(len(boundary)):
            if int(boundary[h][1]) < 3: 
                result.append(["chr"+chrom,(int(boundary[h][0])-2)*25000,(int(boundary[h][0])-1)*25000,int(boundary[h][1])])
            else:
                result.append(["chr"+chrom,(int(boundary[h][0])-2)*25000,(int(boundary[h][0])-1)*25000,"3+"])
                
    savetxt("3.TAD_boundary/"+s+"/25000/boundary_level_3",result)
    
####################50kb#####################################
    os.system("mkdir -p 3.TAD_boundary/"+s+"/50000")
    result = []
    for chrom in chrname:
        left = []
        right = []
        f = open("1.OnTAD_output/"+s+"/50000/OnTAD_ICE_pen0.1_max_40_min_0.6_chr"+chrom+".tad")
        lines=f.readlines()
        for i in range(1,len(lines)):
            line = lines[i].strip()
            sep_line = line.split("\t")
            left.append(sep_line[0])
            right.append(sep_line[1])
        f.close()
        left_count = Counter(left)
        right_count = Counter(right)
        left_name = list(left_count)
        right_name = list(right_count)
        
        boundary = []
        for j in left_name:
            if j not in right_name:
                boundary.append([int(j),int(left_count[j])])
            else:
                if int(left_count[j])>int(right_count[j]):
                    boundary.append([int(j),int(left_count[j])])
                else:
                    boundary.append([int(j),int(right_count[j])])
                    
        for k in right_name:
            if k not in left_name:
                boundary.append([int(k),int(right_count[k])])
                
        for h in range(len(boundary)):
            if int(boundary[h][1]) < 3: 
                result.append(["chr"+chrom,(int(boundary[h][0])-2)*50000,(int(boundary[h][0])-1)*50000,int(boundary[h][1])])
            else:
                result.append(["chr"+chrom,(int(boundary[h][0])-2)*50000,(int(boundary[h][0])-1)*50000,"3+"])
                
    savetxt("3.TAD_boundary/"+s+"/50000/boundary_level_3",result)



