import numpy as np
import os

def savetxt(filename,x):
    np.savetxt(filename,x,delimiter = '\t',fmt='%s')
    
chrname = []
for i in range(1,23):
    chrname.append(str(i))
chrname.append("X")

dirs = ["GM12878","HMEC","HUVEC","IMR90","K562","NHEK","HCT116"]
for s in dirs: 

################5kb#####################################
    os.system("mkdir -p 2.TAD_region/"+s+"/5000")
    domain = []
    for chrom in chrname:
        f = open("1.OnTAD_output/"+s+"/5000/OnTAD_ICE_pen0.1_max_400_min_6_chr"+chrom+".tad")
        lines=f.readlines()
        for i in range(1,len(lines)):
            line = lines[i].strip()
            sep_line = line.split("\t")
            domain.append(["chr"+chrom,str((int(sep_line[0])-1)*5000),str((int(sep_line[1])-1)*5000),int(sep_line[2])])
        f.close()
    
    savetxt("2.TAD_region/"+s+"/5000/TAD_region.bed",domain)

################10kb#####################################
    os.system("mkdir -p 2.TAD_region/"+s+"/10000")
    domain = []
    for chrom in chrname:
        f = open("1.OnTAD_output/"+s+"/10000/OnTAD_ICE_pen0.1_max_200_min_3_chr"+chrom+".tad")
        lines=f.readlines()
        for i in range(1,len(lines)):
            line = lines[i].strip()
            sep_line = line.split("\t")
            domain.append(["chr"+chrom,str((int(sep_line[0])-1)*10000),str((int(sep_line[1])-1)*10000),int(sep_line[2])])
        f.close()
    
    savetxt("2.TAD_region/"+s+"/10000/TAD_region.bed",domain)

################25kb#####################################
    os.system("mkdir -p 2.TAD_region/"+s+"/25000")
    domain = []
    for chrom in chrname:
        f = open("1.OnTAD_output/"+s+"/25000/OnTAD_ICE_pen0.1_max_80_min_1.2_chr"+chrom+".tad")
        lines=f.readlines()
        for i in range(1,len(lines)):
            line = lines[i].strip()
            sep_line = line.split("\t")
            domain.append(["chr"+chrom,str((int(sep_line[0])-1)*25000),str((int(sep_line[1])-1)*25000),int(sep_line[2])])
        f.close()
    
    savetxt("2.TAD_region/"+s+"/25000/TAD_region.bed",domain)

################50kb#####################################
    os.system("mkdir -p 2.TAD_region/"+s+"/50000")
    domain = []
    for chrom in chrname:
        f = open("1.OnTAD_output/"+s+"/50000/OnTAD_ICE_pen0.1_max_40_min_0.6_chr"+chrom+".tad")
        lines=f.readlines()
        for i in range(1,len(lines)):
            line = lines[i].strip()
            sep_line = line.split("\t")
            domain.append(["chr"+chrom,str((int(sep_line[0])-1)*50000),str((int(sep_line[1])-1)*50000),int(sep_line[2])])
        f.close()
    
    savetxt("2.TAD_region/"+s+"/50000/TAD_region.bed",domain)






