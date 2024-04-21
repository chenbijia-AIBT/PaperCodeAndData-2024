library(data.table)
library(ggplot2)
library(tidyverse)
library(plyr)
library(ggpubr)
library(ggthemes)
library(scales)
####### HCT116 #############
# 为基因组区域分配边界level
df <- fread("F:/project1/9_new_results/0.genome/10000_abs.bed")
boundary <- fread("F:/project1/9_new_results/0.TAD/boundary/HCT116/boundary_level_all")
boundary$V5 <- gsub(pattern = "[3-9]",replacement = "3+",x = boundary$V4)
genome <- merge(df,boundary,by=c("V1","V2","V3"),all.x=T)
colnames(genome) <- c("V1","V2","V3","V4","contact","levels_contact")
genome$contact[is.na(genome$contact)] <- 0
genome$levels_contact[is.na(genome$levels_contact)] <- 0

## DNase-seq bigwig信号
bigwig <- fread("F:/project1/9_new_results/0.TFCR/HCT116/DNase-seq_data/bigwig/10000_readCounts.tab")
colnames(bigwig) <- c("V1","V2","V3","signal")
g1 <- merge(genome,bigwig,by=c("V1","V2","V3"),all.x=T)

#DSB数目
dsb <- fread("F:/project1/9_new_results/0.DSB/HCT116_NT_END.narrowPeak")
dsb$number <- 1
chroms.vector <- c(paste0("chr",c(1:22,"X","Y")))
g1_d <- lapply(chroms.vector, FUN=function(temp.chrom){
  temp.g1 <- g1[V1==temp.chrom]
  temp.dsb <- dsb[V1==temp.chrom]
  aa=temp.g1[, {
    temp.dsb[V2 <= g1.V3 & V3 >= g1.V2][,sum(number)]
  },list(V1,g1.V2=V2,g1.V3=V3,V4,contact,levels_contact,signal)]
})
g2 <- do.call("rbind", g1_d)  
colnames(g2) <- c("seqnames","start","end","bin","contact","levels_contact","signal","DSB_number")

#基因表达数据 
exp <- fread("F:/project1/9_new_results/0.TFCR/HCT116/gene_exp/gene_expression.tsv")
exp <- separate(exp,col="gene_id",into=c("gene_id","gene_id_1"),sep = "[.]")
exp <- exp[,c(1,8)]
protein_gene <- fread("F:/project1/4_gene_and_boundary/hg19.annotation/hg19_protein_gene_promoter_2kb")
protein_gene <- separate(protein_gene,col="gene_id",into=c("gene_id","gene_id_1"),sep = "[.]")
merge_data <- merge(protein_gene,exp,by="gene_id",all.x=T)
sum(is.na(merge_data$FPKM)) 
merge_data <- merge_data[,-5]
merge_data <- merge_data[,c(2,8,9,1,5:7,10)] 
gene <- as.data.table(merge_data)
g2_g <- lapply(chroms.vector, FUN=function(temp.chrom){
  temp.g2 <- g2[seqnames==temp.chrom]
  temp.gene <- gene[seqnames==temp.chrom]
  aa=temp.g2[, {
    temp.gene[promoter_start <= g2.end & promoter_end >= g2.start][,paste(gene_len,transcript_num,FPKM)]
  },list(seqnames,g2.start=start,g2.end=end,bin,contact,levels_contact,signal,DSB_number)]
})
g2_gene <- do.call("rbind", g2_g)  
g2_gene <- separate(g2_gene,col="V1",into=c("gene_len","transcript_num","gene_exp"),sep = "[ ]")
colnames(g2_gene)[colnames(g2_gene) %in% c("g2.start","g2.end")] <- c("start","end")
g2_gene <- transform(g2_gene, gene_exp = as.numeric(gene_exp))
g2_gene$gene_exp1 <- g2_gene$gene_exp
meanexp <- aggregate(gene_exp ~ seqnames+start+end, data = g2_gene, mean)
colnames(meanexp)[4] <- "meanexp"
maxexp <- aggregate(gene_exp1 ~ seqnames+start+end, data = g2_gene, max)
colnames(maxexp)[4] <- "maxexp"
g3 <- merge(g2,meanexp,by=c("seqnames","start","end"),all.x=T)
g3 <- merge(g3,maxexp,by=c("seqnames","start","end"),all.x=T)
g3 <- transform(g3, meanexp = as.numeric(meanexp), maxexp = as.numeric(maxexp))
g2_gene <- transform(g2_gene, gene_len = as.numeric(gene_len))
g2_gene$gene_len1 <- g2_gene$gene_len
meanlen <- aggregate(gene_len ~ seqnames+start+end, data = g2_gene, mean)
colnames(meanlen)[4] <- "meanlen"
maxlen <- aggregate(gene_len1 ~ seqnames+start+end, data = g2_gene, max)
colnames(maxlen)[4] <- "maxlen"
g3 <- merge(g3,meanlen,by=c("seqnames","start","end"),all.x=T)
g3 <- merge(g3,maxlen,by=c("seqnames","start","end"),all.x=T)
g2_gene <- transform(g2_gene, transcript_num = as.numeric(transcript_num))
g2_gene$transcript_num1 <- g2_gene$transcript_num
meantranscript <- aggregate(transcript_num ~ seqnames+start+end, data = g2_gene, mean)
colnames(meantranscript)[4] <- "meantranscript"
maxtranscript <- aggregate(transcript_num1 ~ seqnames+start+end, data = g2_gene, max)
colnames(maxtranscript)[4] <- "maxtranscript"
g3 <- merge(g3,meantranscript,by=c("seqnames","start","end"),all.x=T)
g3 <- merge(g3,maxtranscript,by=c("seqnames","start","end"),all.x=T)

#基因数目
gene_a <- subset(gene,FPKM>5)
gene_a$number <- 1
chroms.vector <- c(paste0("chr",c(1:22,"X","Y")))
g3_g <- lapply(chroms.vector, FUN=function(temp.chrom){
  temp.g3 <- g3[seqnames==temp.chrom]
  temp.gene_a <- gene_a[seqnames==temp.chrom]
  aa=temp.g3[, {
    temp.gene_a[promoter_start <= g3.end & promoter_end >= g3.start][,sum(number)]
  },list(seqnames,g3.start=start,g3.end=end,bin,contact,levels_contact,signal,DSB_number,meanexp,maxexp,meanlen,maxlen,meantranscript,maxtranscript)]
})
g62 <- do.call("rbind", g3_g)  
colnames(g62)[colnames(g62) %in% c("g3.start","g3.end","V1")] <- c("start","end","active_gene_number")
summary(g62)

#必需基因数目
e_gene <- fread("F:/project1/4_gene_and_boundary/essential_gene/essential_gene_promoter_2kb")
gene_e <- subset(e_gene,type==0)
gene_e$number <- 1
chroms.vector <- c(paste0("chr",c(1:22,"X","Y")))
g62_g <- lapply(chroms.vector, FUN=function(temp.chrom){
  temp.g62 <- g62[seqnames==temp.chrom]
  temp.gene_e <- gene_e[chromosome==temp.chrom]
  aa=temp.g62[, {
    temp.gene_e[promoter_start <= g62.end & promoter_end >= g62.start][,sum(number)]
  },list(seqnames,g62.start=start,g62.end=end,bin,contact,levels_contact,signal,DSB_number,meanexp,maxexp,meanlen,maxlen,meantranscript,maxtranscript,active_gene_number)]
})
g6 <- do.call("rbind", g62_g)  
colnames(g6)[colnames(g6) %in% c("g62.start","g62.end","V1")] <- c("start","end","essential_gene_number")

# 必需基因的一些性质
g6_g <- lapply(chroms.vector, FUN=function(temp.chrom){
  temp.g6 <- g6[seqnames==temp.chrom]
  temp.e_gene <- e_gene[chromosome==temp.chrom]
  aa=temp.g6[, {
    temp.e_gene[promoter_start <= g6.end & promoter_end >= g6.start][,paste(CS,`gene_age(12)`)]
  },list(seqnames,g6.start=start,g6.end=end,bin,contact,levels_contact,signal,DSB_number,meanexp,maxexp,meanlen,maxlen,meantranscript,maxtranscript,active_gene_number,essential_gene_number)]
})
g6_gene <- do.call("rbind", g6_g)  
g6_gene <- separate(g6_gene,col="V1",into=c("CS","gene_age"),sep = "[ ]")
colnames(g6_gene)[colnames(g6_gene) %in% c("g6.start","g6.end")] <- c("start","end")
summary(g6_gene)
g6_gene <- transform(g6_gene, CS = as.numeric(CS), gene_age = as.numeric(gene_age))
meanCS <- aggregate(CS ~ seqnames+start+end, data = g6_gene, mean)
colnames(meanCS)[4] <- "meanCS"
g61 <- merge(g6,meanCS,by=c("seqnames","start","end"),all.x=T)
meanage <- aggregate(gene_age ~ seqnames+start+end, data = g6_gene, mean)
colnames(meanage)[4] <- "meanage"
g62 <- merge(g61,meanage,by=c("seqnames","start","end"),all.x=T)

g7 <- g62
g7$levels_contact <- factor(g7$levels_contact, levels = c("0","1","2","3+"))
g71 <- subset(g7,signal==0)
g71$levels_signal <- 0
g72 <- subset(g7,signal!=0)
g72 <- g72[order(g72$signal)] 
g72$levels_signal <- as.character(as.numeric(cut(1:nrow(g72),breaks = 3)))
g7 <- rbind(g71,g72)
g7$levels_signal <- factor(g7$levels_signal, levels = c("0","1","2","3"))
g7 <- g7 %>% filter(!is.na(meanexp))

g7_mean <- aggregate(cbind(DSB_number) ~ levels_signal + levels_contact, 
                     data = g7, 
                     FUN = mean)
g7_agg <- aggregate(cbind(levels_signal,levels_contact) ~ levels_signal + levels_contact, 
                    data = g7, 
                    FUN = length)
g7_agg <- g7_agg[,1:3]
colnames(g7_agg)[3] <- "bin_count"
g7_mean <- left_join(g7_mean,g7_agg,by=names(g7_agg)[1:2])

##方块图（热图）
ggplot(g7_mean,aes(x = levels_signal,y = levels_contact,fill = DSB_number)) + 
  geom_tile() +
  geom_text(aes(label=bin_count),size=2.5)+     
  scale_fill_distiller(palette = "Spectral")+
  theme_few() +
  labs(x = "DNase-seq Signal",y = "Bin level",fill="DSB density")+
  theme(axis.text=element_text(size=9,face = "plain"),
        legend.text = element_text(size = 9, face = "plain"),    
        legend.key.size=unit(0.3,'cm'))

##气泡图
g7_mean <- aggregate(cbind(meanexp,maxexp,meanlen,maxlen,meantranscript,maxtranscript,active_gene_number,essential_gene_number,meanCS,meanage) ~ levels_signal + levels_contact, 
                     data = g7, 
                     FUN = mean)
g7_agg <- aggregate(cbind(levels_signal,levels_contact) ~ levels_signal + levels_contact, 
                    data = g7, 
                    FUN = length)
g7_agg <- g7_agg[,1:3]
colnames(g7_agg)[3] <- "bin_count"
g7_mean <- left_join(g7_mean,g7_agg,by=names(g7_agg)[1:2])
ggplot(g7_mean,aes(x=levels_signal,y=levels_contact,color=meanexp,size=active_gene_number))+
  geom_point()+
  scale_color_gradientn(colours = c("white", muted("blue")))+ 
  theme_few() +
  labs(x = "DNase-seq Signal",y = "Bin level",color="Gene exp.",size="Active gene density") +
  theme(axis.text=element_text(size=9,face = "plain"),
        legend.text = element_text(size = 9, face = "plain"),      
        legend.key.size=unit(0.3,'cm')) +
  guides(size=guide_legend(order = 1),colorbar=guide_legend(order=2))

ggplot(g7_mean,aes(x=levels_signal,y=levels_contact,color=meanCS,size=essential_gene_number))+
  geom_point()+
  scale_color_gradientn(colours = c("darkred", "orange", "white"))+ 
  theme_few() +
  labs(x = "DNase-seq Signal",y = "Bin level",color="CS value",size="Essential gene density") +
  theme(axis.text=element_text(size=9,face = "plain"),
        legend.text = element_text(size = 9, face = "plain"),      
        legend.key.size=unit(0.3,'cm'))+
  guides(size=guide_legend(order = 1),colorbar=guide_legend(order=2))
options(scipen=200)
ggplot(g7_mean,aes(x=levels_signal,y=levels_contact,color=meanage,size=meanlen))+
  geom_point()+
  scale_color_gradientn(colours = c("darkgreen", "white"))+ 
  theme_few() +
  labs(x = "DNase-seq Signal",y = "Bin level",color="Gene age",size="Gene length") +
  theme(axis.text=element_text(size=9,face = "plain"),
        legend.text = element_text(size = 9, face = "plain"),      
        legend.key.size=unit(0.3,'cm'))+
  guides(size=guide_legend(order = 1),colorbar=guide_legend(order=2))

##DSB的bigwig信号值
bigwig <- fread("F:/project1/9_new_results/0.genome/DSB_signal/HCT116/10000/readCounts.tab")
colnames(bigwig) <- c("seqnames","start","end","DSB_signal")
g63 <- merge(g7,bigwig,by=c("seqnames","start","end"),all.x=T)
g8 <- g63 

g8_mean <- aggregate(cbind(DSB_signal) ~ levels_signal + levels_contact, 
                     data = g8, 
                     FUN = mean)
g8_agg <- aggregate(cbind(levels_signal,levels_contact) ~ levels_signal + levels_contact, 
                    data = g8, 
                    FUN = length)
g8_agg <- g8_agg[,1:3]
colnames(g8_agg)[3] <- "bin_count"
g8_mean <- left_join(g8_mean,g8_agg,by=names(g8_agg)[1:2])
##方块图（热图）
ggplot(g8_mean,aes(x = levels_signal,y = levels_contact,fill = DSB_signal)) + 
  geom_tile() +
  geom_text(aes(label=bin_count),size=2.5)+   
  scale_fill_distiller(palette = "Spectral")+
  theme_few() +
  labs(x = "DNase-seq Signal",y = "Bin level",fill="DSB bigWig signal")+
  theme(axis.text=element_text(size=9,face = "plain"),
        legend.text = element_text(size = 9, face = "plain"),    
        legend.key.size=unit(0.3,'cm'))

####箱图
##根据bin中的DSB数目分类比较
g70 <- g62
g70$levels_DSB <- g70$DSB_number
g70$levels_DSB[which(g70$levels_DSB > 2)] <- "3+"
g70$levels_DSB <- factor(g70$levels_DSB, levels = c("0","1","2","3+"))
my_comparisons1 <- list(c("0", "1"),
                        c("1", "2"),
                        c("2", "3+"))

ggplot(g70, aes(x=levels_DSB, y=contact)) +
  stat_boxplot(geom ="errorbar",width=0.3) +
  geom_boxplot(aes(fill=levels_DSB),size=0.4,outlier.size=0.8,width=0.7,outlier.shape = NA) +
  scale_fill_manual(values = c("#26547b","#20707c","#1f9077","#2fac67")) + 
  theme_few() +
  labs(x = "DSB number", y = "Bin level", title = "") +  
  guides(fill = "none") +       
  theme(axis.title = element_text(size=12),
        axis.text.y=element_text(size=10),
        axis.text.x=element_text(size=11),
        plot.title = element_text(hjust = 0.5,size=12),
        axis.ticks.length=unit(0.1,'cm')) + 
  stat_compare_means(comparisons=my_comparisons1,
                     aes(label=..p.signif..), 
                     label.x = 1.5, 
                     tip.length = 0.002)
aggregate(contact ~ levels_DSB, data = g70, FUN = mean)
g_HCT116 <- g70[,c(5,19)]
g_HCT116$cell_line <- "HCT116"
#方差分析（ANOVA）
aov1 <- aov(contact ~ levels_DSB, data = g70)  
aov1
summary(aov1)

ggplot(g70, aes(x=levels_DSB, y=signal)) +
  stat_boxplot(geom ="errorbar",width=0.3) +
  geom_boxplot(aes(fill=levels_DSB),size=0.4,outlier.size=0.8,width=0.7,outlier.shape = NA) +
  scale_fill_manual(values = c("#26547b","#20707c","#1f9077","#2fac67")) + 
  theme_few() +
  labs(x = "DSB number", y = "DNase-seq Signal", title = "") +  
  guides(fill = "none") +  
  coord_cartesian(ylim = c(0, 0.4)) +          
  theme(axis.title = element_text(size=12),
        axis.text.y=element_text(size=10),
        axis.text.x=element_text(size=11),
        plot.title = element_text(hjust = 0.5,size=12),
        axis.ticks.length=unit(0.1,'cm')) + 
  stat_compare_means(comparisons=my_comparisons1,
                     aes(label=..p.signif..), 
                     label.x = 1.5, 
                     step.increase = 0.05,
                     tip.length = 0.001)

ggplot(g70, aes(x=levels_DSB, y=log10(meanexp+1))) +
  stat_boxplot(geom ="errorbar",width=0.3) +
  geom_boxplot(aes(fill=levels_DSB),size=0.4,outlier.size=0.8,width=0.7,outlier.shape = NA) +
  scale_fill_manual(values = c("#26547b","#20707c","#1f9077","#2fac67")) + 
  theme_few() +
  labs(x = "DSB number", y = "log10(Gene Expression+1)", title = "") +  
  guides(fill = "none") +  
  coord_cartesian(ylim = c(0, 5)) +          
  theme(axis.title = element_text(size=12),
        axis.text.y=element_text(size=10),
        axis.text.x=element_text(size=11),
        plot.title = element_text(hjust = 0.5,size=12),
        axis.ticks.length=unit(0.1,'cm')) + 
  stat_compare_means(comparisons=my_comparisons1,
                     aes(label=..p.signif..), 
                     label.x = 1.5, 
                     label.y = 4,
                     step.increase = 0.1,
                     tip.length = 0.01)

ggplot(g70, aes(x=levels_DSB, y=meanCS)) +
  stat_boxplot(geom ="errorbar",width=0.3) +
  geom_boxplot(aes(fill=levels_DSB),size=0.4,outlier.size=0.8,width=0.7,outlier.shape = NA) +
  scale_fill_manual(values = c("#26547b","#20707c","#1f9077","#2fac67")) + 
  theme_few() +
  labs(x = "DSB number", y = "CS value of genes", title = "") +  
  guides(fill = "none") +  
  coord_cartesian(ylim = c(-2, 2.3)) +          
  theme(axis.title = element_text(size=12),
        axis.text.y=element_text(size=10),
        axis.text.x=element_text(size=11),
        plot.title = element_text(hjust = 0.5,size=12),
        axis.ticks.length=unit(0.1,'cm')) + 
  stat_compare_means(comparisons=my_comparisons1,
                     aes(label=..p.signif..), 
                     label.y = 1.4,
                     step.increase = 0.05,
                     tip.length = 0.01)


### 比较bin level均值 #########
data_merge <- rbind(g_HCT116,g_NHEK_b,g_NHEK_c)
data_merge$levels_DSB <- factor(data_merge$levels_DSB, levels = c("0","1","2","3+"))
ggline(data_merge,x="levels_DSB",y="contact",
       color = "cell_line",
       shape = "cell_line",
       add = c("mean"),
       size = 0.8)+ 
  scale_color_manual(values = c("#ddb96f","#ad88a2","#6dadac")) + 
  labs(x = "DSB number", y = "Mean bin-level", title = "") +  
  theme_bw() +
  theme(axis.title = element_text(size=12),
        axis.text.y=element_text(size=12),
        axis.text.x=element_text(size=12),
        axis.ticks.length=unit(0.1,'cm'),
        legend.title = element_blank(),
        legend.position = c(0.72,0.2)) 
aggregate(contact ~ levels_DSB + cell_line, data = data_merge, mean)

#方差分析（ANOVA）
file_dsb <- list(g_HCT116,g_NHEK_b,g_NHEK_c)
file_aov <- NULL
for (f in file_dsb){
  aov1 <- aov(contact ~ levels_DSB, data = f)
  test <- summary(aov1)[[1]][1,]
  file_aov <- rbind(file_aov,test)
}


### 统计不同DSB数目的bin数目 #######
# 创建数据
df_count <- data.frame(cell_line=rep(c("HCT116","NHEK (BLESS)","NHEK (DSBCapture)"),each=4),  
                       DSB_number=rep(c("0","1","2","3+"),times=3),   
                       count=c(238208,52025,13510,5836,291344,17175,1026,34,242430,49755,14103,3291))
df_count$cell_line <- factor(df_count$cell_line,levels = c("NHEK (DSBCapture)","NHEK (BLESS)","HCT116"))
df_count$DSB_number <- factor(df_count$DSB_number,levels = c("3+","2","1","0"))

ggplot(df_count, aes(x = cell_line,y = count,fill = DSB_number))+
  geom_bar(stat="identity",width = 0.6,position = position_dodge(width=0.7),alpha = 0.9)+ 
  scale_fill_manual(values = c("#2fac67","#1f9077","#20707c","#26547b"))+
  geom_text(aes(label = count), position=position_dodge(width = 0.7),size = 3,hjust = 0)+
  labs(x=" ",y="Bin count",fill="DSB number",title="")+
  theme_classic()+
  theme(axis.ticks.length=unit(0.1,'cm'))+
  scale_y_break(c(53000,230000),scales = "free")+   
  coord_flip()+
  guides(fill = guide_legend(reverse = T))






