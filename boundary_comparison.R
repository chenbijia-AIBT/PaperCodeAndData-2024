library(data.table)
library(ggplot2)
library(tidyverse)
library(plyr)
library(ggthemes)
library(pacman)
library(ggpubr)
#all cell lines#############
#*不同level边界的比较 #########
# protein-coding genes
### **基因表达量 ##########
protein_gene <- fread("F:/project1/4_gene_and_boundary/hg19.annotation/hg19_protein_gene_promoter_2kb")
protein_gene <- separate(protein_gene,col="gene_id",into=c("gene_id","gene_id_1"),sep = "[.]")
gene_select <- protein_gene %>% dplyr::select(1,9,10,4,7,8)
filename <- c("GM12878","HMEC","HUVEC","IMR90","K562","NHEK","HCT116")
chroms.vector <- c(paste0("chr",c(1:22,"X","Y")))
file_final <- NULL
file_final1 <- NULL
for (f in filename){
  exp <- fread(paste0("F:/project1/9_new_results/0.TFCR/",f,"/gene_exp/gene_expression.tsv"))
  exp <- separate(exp,col="gene_id",into=c("gene_id","gene_id_1"),sep = "[.]")
  exp <- exp[,c(1,8)]
  gene <- as.data.table(left_join(gene_select,exp,by="gene_id"))
  b <- fread(paste0("F:/project1/9_new_results/0.TAD/boundary/",f,"/boundary_level_3"))
  g_b <- lapply(chroms.vector, FUN=function(temp.chrom){
    temp.gene <- gene[seqnames==temp.chrom]
    temp.b <- b[V1==temp.chrom]
    aa=temp.gene[, {
      temp.b[V2 <= gene.promoter_end & V3 >= gene.promoter_start][,paste(V4)]
    },list(seqnames,gene.promoter_start=promoter_start,gene.promoter_end=promoter_end,gene_id,gene_len,transcript_num,FPKM)]
  })
  gene_b <- do.call("rbind", g_b)
  data <- gene_b %>% group_by(seqnames,gene.promoter_start,gene.promoter_end) %>% filter(V1 == max(V1))  
  colnames(data)[colnames(data) %in% c("gene.promoter_start","gene.promoter_end","V1")] <- c("promoter_start","promoter_end","level")
  data1 <- merge(gene,data,all.x=T)
  data1$level[is.na(data1$level)] <- 0
  data1$cell_lines <- f
  file_final <- rbind(file_final,data1)  
  data$cell_lines <- f
  file_final1 <- rbind(file_final1,data)
}
file_final$level <- factor(file_final$level,levels=c("0","1","2","3+"))
file_final$level_group <- c("0"="0","1"="1+","2"="1+","3+"="1+")[as.character(file_final$level)]
file_final$level_group <- factor(file_final$level_group,levels=c("0","1+"))
file_final1$level <- factor(file_final1$level,levels=c("1","2","3+"))

file_final$FPKM[file_final$FPKM > 7000] <- 7000
file_final$log10 <- log10(file_final$FPKM+1)
stat_wilcox <- rstatix::wilcox_test(group_by(file_final, cell_lines),log10 ~ level_group)  #获取组间两两比较的 p 值
stat_wilcox <- rstatix::add_significance(stat_wilcox, "p")   #根据 p 值添加显著性标记 * 符号
stat_wilcox.test <-  rstatix::add_xy_position(stat_wilcox, x = "cell_lines")
ggplot(file_final, aes(x=cell_lines, y=log10)) +
  geom_boxplot(aes(fill=level_group),size=0.4,outlier.size=0.8,width=0.6,outlier.shape = NA) +
  scale_fill_manual(values = c("#f6f6f6","#8d5108")) +
  theme_few() +
  labs(x = "", y = "log10(Gene Expression+1)", title = "", fill = "level") +  
  scale_y_continuous(limits = c(0, 4.5)) +          
  theme(axis.title = element_text(size=11),
        axis.text.y=element_text(size=10),
        axis.text.x=element_text(size=10,angle = 45,vjust = 1,hjust = 1),
        axis.ticks.length=unit(0.1,'cm')) +
  stat_pvalue_manual(stat_wilcox.test,
                     label = 'p.signif',
                     tip.length = 0.001,
                     y.position = 4,
                     hide.ns = TRUE)
ggline(file_final,x="level",y="FPKM",
       color = "cell_lines",
       add = c("median"),
       size = 0.8)+ 
  scale_color_manual(values = c("#8dd3c7","#fdb461","#e78ac4","#fb8073","#80b1d3","#e5c494","#999999")) + 
  labs(x = "Boundary level", y = "Median Gene Expression", title = "") +  
  theme_bw() +
  theme(axis.title = element_text(size=12),
        axis.text.y=element_text(size=11),
        axis.text.x=element_text(size=11),
        axis.ticks.length=unit(0.1,'cm'),
        legend.title = element_blank(),
        legend.position = "right") 
aggregate(FPKM ~ level + cell_lines, data = file_final, mean)

file_final1$FPKM[file_final1$FPKM > 7000] <- 7000
file_final1$log10 <- log10(file_final1$FPKM+1)
stat_wilcox <- rstatix::wilcox_test(group_by(file_final1, cell_lines),log10 ~ level)  
stat_wilcox <- rstatix::add_significance(stat_wilcox, "p")   
stat_wilcox.test <-  rstatix::add_xy_position(stat_wilcox, x = "cell_lines")
ggplot(file_final1, aes(x=cell_lines, y=log10)) +
  geom_boxplot(aes(fill=level),size=0.4,outlier.size=0.8,width=0.6,outlier.shape = NA) +
  scale_fill_manual(values = c("#f3e3b7","#d6b66a","#af6e22")) + 
  theme_few() +
  labs(x = "", y = "log10(Gene Expression+1)", title = "") +      
  theme(axis.title = element_text(size=11),
        axis.text.y=element_text(size=10),
        axis.text.x=element_text(size=10,angle = 45,vjust = 1,hjust = 1),
        axis.ticks.length=unit(0.1,'cm')) + 
  stat_pvalue_manual(stat_wilcox.test, 
                     label = 'p.signif',
                     tip.length = 0.001,
                     hide.ns = TRUE)


###**基因长度########################
ggline(file_final,x="level",y="gene_len",
       color = "cell_lines",
       add = c("median"),
       size = 0.8)+ 
  scale_color_manual(values = c("#8dd3c7","#fdb461","#e78ac4","#fb8073","#80b1d3","#e5c494","#999999")) + 
  labs(x = "Boundary level", y = "Median Gene Length", title = "") +  
  theme_bw() +
  theme(axis.title = element_text(size=12),
        axis.text.y=element_text(size=11),
        axis.text.x=element_text(size=11),
        axis.ticks.length=unit(0.1,'cm'),
        legend.title = element_blank(),
        legend.position = "right") 

file_final$gene_len[file_final$gene_len > 200000] <- 200000
stat_wilcox <- rstatix::wilcox_test(group_by(file_final, cell_lines),gene_len ~ level_group) 
stat_wilcox <- rstatix::add_significance(stat_wilcox, "p")  
stat_wilcox.test <-  rstatix::add_xy_position(stat_wilcox, x = "cell_lines")
ggplot(file_final, aes(x=cell_lines, y=gene_len)) +
  geom_boxplot(aes(fill=level_group),size=0.4,outlier.size=0.8,width=0.6,outlier.shape = NA) +
  scale_fill_manual(values = c("#f6f6f6","#8d5108")) + 
  theme_few() +
  labs(x = "", y = "Gene Length", title = "", fill = "level")  +  
  scale_y_continuous(limits = c(0, 220000))+          
  theme(axis.title = element_text(size=11),
        axis.text.y=element_text(size=9),
        axis.text.x=element_text(size=10,angle = 45,vjust = 1,hjust = 1),
        axis.ticks.length=unit(0.1,'cm')) + 
  stat_pvalue_manual(stat_wilcox.test, 
                     label = 'p.signif',
                     tip.length = 0.001,
                     y.position = 200000,
                     hide.ns = TRUE)

file_final1$gene_len[file_final1$gene_len > 200000] <- 200000
stat_wilcox <- rstatix::wilcox_test(group_by(file_final1, cell_lines),gene_len ~ level)  
stat_wilcox <- rstatix::add_significance(stat_wilcox, "p")   
stat_wilcox.test <-  rstatix::add_xy_position(stat_wilcox, x = "cell_lines")
ggplot(file_final1, aes(x=cell_lines, y=gene_len)) +
  geom_boxplot(aes(fill=level),size=0.4,outlier.size=0.8,width=0.6,outlier.shape = NA) +
  scale_fill_manual(values = c("#f3e3b7","#d6b66a","#af6e22")) + 
  theme_few() +
  labs(x = "", y = "Gene Length", title = "") +       
  theme(axis.title = element_text(size=11),
        axis.text.y=element_text(size=9),
        axis.text.x=element_text(size=10,angle = 45,vjust = 1,hjust = 1),
        axis.ticks.length=unit(0.1,'cm')) + 
  stat_pvalue_manual(stat_wilcox.test, 
                     label = 'p.signif',
                     tip.length = 0.001,
                     step.increase = 0.05,
                     hide.ns = TRUE)

###**转录本数目#############################
ggline(file_final,x="level",y="transcript_num",
       color = "cell_lines",
       add = c("mean"),
       size = 0.8)+
  scale_color_manual(values = c("#8dd3c7","#fdb461","#e78ac4","#fb8073","#80b1d3","#e5c494","#999999")) +
  labs(x = "Boundary level", y = "Mean Transcript Num", title = "") +
  theme_bw() +
  theme(axis.title = element_text(size=12),
        axis.text.y=element_text(size=11),
        axis.text.x=element_text(size=11),
        axis.ticks.length=unit(0.1,'cm'),
        legend.title = element_blank(),
        legend.position = "right")

stat_wilcox <- rstatix::wilcox_test(group_by(file_final, cell_lines),transcript_num ~ level_group)  
stat_wilcox <- rstatix::add_significance(stat_wilcox, "p")   
stat_wilcox.test <-  rstatix::add_xy_position(stat_wilcox, x = "cell_lines")
ggplot(file_final, aes(x=cell_lines, y=transcript_num)) +
  geom_boxplot(aes(fill=level_group),size=0.4,outlier.size=0.8,width=0.6,outlier.shape = NA) +
  scale_fill_manual(values = c("#f6f6f6","#8d5108")) +
  theme_few() +
  labs(x = "", y = "Number of transcript types", title = "", fill = "level") +  
  scale_y_continuous(limits = c(0, 28)) +          
  theme(axis.title = element_text(size=11),
        axis.text.y=element_text(size=10),
        axis.text.x=element_text(size=10,angle = 45,vjust = 1,hjust = 1),
        axis.ticks.length=unit(0.1,'cm')) +
  stat_pvalue_manual(stat_wilcox.test,
                     label = 'p.signif',
                     tip.length = 0.001,
                     y.position = 25,
                     hide.ns = TRUE)

###**基因CS值（必需性）################################
e_gene <- fread("F:/project1/10_reorganize/0_data_source/TH_score/essential_genes_promoter2000_sorted.bed")
egene_select <- e_gene %>% dplyr::select(1,12,13,8,5,6,9,10,14)
colnames(egene_select) <- c("seqnames","promoter_start","promoter_end","gene_name","CS","type","transcript_num","gene_age","HK")
gene <- egene_select
filename <- c("GM12878","HMEC","HUVEC","IMR90","K562","NHEK","HCT116")
chroms.vector <- c(paste0("chr",c(1:22,"X","Y")))
file_final <- NULL
file_final1 <- NULL
for (f in filename){
  b <- fread(paste0("F:/project1/9_new_results/0.TAD/boundary/",f,"/boundary_level_all"))
  g_b <- lapply(chroms.vector, FUN=function(temp.chrom){
    temp.gene <- gene[seqnames==temp.chrom]
    temp.b <- b[V1==temp.chrom]
    aa=temp.gene[, {
      temp.b[V2 <= gene.promoter_end & V3 >= gene.promoter_start][,paste(V4)]
    },list(seqnames,gene.promoter_start=promoter_start,gene.promoter_end=promoter_end,gene_name,CS,type,transcript_num,gene_age,HK)]
  })
  gene_b <- do.call("rbind", g_b)
  data <- gene_b %>% group_by(seqnames,gene.promoter_start,gene.promoter_end) %>% filter(V1 == max(V1))  
  colnames(data)[colnames(data) %in% c("gene.promoter_start","gene.promoter_end","V1")] <- c("promoter_start","promoter_end","level")
  data1 <- merge(gene,data,all.x=T)
  data1$level[is.na(data1$level)] <- 0
  data1$cell_lines <- f
  file_final <- rbind(file_final,data1)
  data$cell_lines <- f
  file_final1 <- rbind(file_final1,data)
}
file_final$level[which(file_final$level>2)] <- "3+"
file_final$level <- factor(file_final$level,levels=c("0","1","2","3+"))
file_final$level_group <- c("0"="0","1"="1+","2"="1+","3+"="1+")[as.character(file_final$level)]
file_final$level_group <- factor(file_final$level_group,levels=c("0","1+"))
file_final1$level[which(file_final1$level>2)] <- "3+"
file_final1$level <- factor(file_final1$level,levels=c("1","2","3+"))

stat_wilcox <- rstatix::wilcox_test(group_by(file_final, cell_lines),CS ~ level_group)  
stat_wilcox <- rstatix::add_significance(stat_wilcox, "p")  
stat_wilcox.test <-  rstatix::add_xy_position(stat_wilcox, x = "cell_lines")
ggplot(file_final, aes(x=cell_lines, y=CS)) +
  geom_boxplot(aes(fill=level_group),size=0.4,outlier.size=0.8,width=0.6,outlier.shape = NA) +
  scale_fill_manual(values = c("#f6f6f6","#8d5108")) + 
  theme_few() +
  labs(x = "", y = "CRISPR score (CS) of genes", title = "", fill = "level")  +  
  scale_y_continuous(limits = c(-2, 2)) +          
  theme(axis.title = element_text(size=11),
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10,angle = 45,vjust = 1,hjust = 1),
        plot.title = element_text(hjust = 0.5,size=12),
        axis.ticks.length=unit(0.1,'cm')) + 
  stat_pvalue_manual(stat_wilcox.test, 
                     label = 'p.signif',
                     tip.length = 0.001,
                     hide.ns = TRUE)
ggline(file_final,x="level",y="CS",
       color = "cell_lines",
       add = c("median"),
       size = 0.8)+ 
  scale_color_manual(values = c("#8dd3c7","#fdb461","#e78ac4","#fb8073","#80b1d3","#e5c494","#999999")) + 
  labs(x = "Boundary level", y = "Median CRISPR score (CS)", title = "") +  
  theme_bw() +
  theme(axis.title = element_text(size=12),
        axis.text.y=element_text(size=11),
        axis.text.x=element_text(size=11),
        axis.ticks.length=unit(0.1,'cm'),
        legend.title = element_blank(),
        legend.position = "right") 

stat_wilcox <- rstatix::wilcox_test(group_by(file_final1, cell_lines),CS ~ level)  
stat_wilcox <- rstatix::add_significance(stat_wilcox, "p")   
stat_wilcox.test <-  rstatix::add_xy_position(stat_wilcox, x = "cell_lines")
ggplot(file_final1, aes(x=cell_lines, y=CS)) +
  geom_boxplot(aes(fill=level),size=0.4,outlier.size=0.8,width=0.6,outlier.shape = NA) +
  scale_fill_manual(values = c("#f3e3b7","#d6b66a","#af6e22")) + 
  theme_few() +
  labs(x = "", y = "CRISPR score (CS) of genes", title = "") +  
  scale_y_continuous(limits = c(-2, 2)) +          
  theme(axis.title = element_text(size=11),
        axis.text.y=element_text(size=10),
        axis.text.x=element_text(size=10,angle = 45,vjust = 1,hjust = 1),
        axis.ticks.length=unit(0.1,'cm')) + 
  stat_pvalue_manual(stat_wilcox.test, 
                     label = 'p.signif',
                     tip.length = 0.001,
                     hide.ns = TRUE)

## **基因年龄 #######################
ggline(file_final,x="level",y="gene_age",
       color = "cell_lines",
       add = c("mean"),
       size = 0.8)+ 
  scale_color_manual(values = c("#8dd3c7","#fdb461","#e78ac4","#fb8073","#80b1d3","#e5c494","#999999")) + 
  labs(x = "Boundary level", y = "Mean gene age", title = "") +  
  theme_bw() +
  theme(axis.title = element_text(size=12),
        axis.text.y=element_text(size=11),
        axis.text.x=element_text(size=11),
        axis.ticks.length=unit(0.1,'cm'),
        legend.title = element_blank(),
        legend.position = "right") 

#方差分析（ANOVA）
filename <- c("GM12878","HMEC","HUVEC","IMR90","K562","NHEK","HCT116")
file_aov <- NULL
for (f in filename){
  file_split <- file_final %>% filter(cell_lines==f)
  aov1 <- aov(gene_age ~ level, data = file_split)
  test <- summary(aov1)[[1]][1,]
  test$cell_lines <- f
  file_aov <- rbind(file_aov,test)
}


stat_wilcox <- rstatix::wilcox_test(group_by(file_final, cell_lines),gene_age ~ level_group)  
stat_wilcox <- rstatix::add_significance(stat_wilcox, "p")  
stat_wilcox.test <-  rstatix::add_xy_position(stat_wilcox, x = "cell_lines")
ggplot(file_final, aes(x=cell_lines, y=gene_age)) +
  geom_boxplot(aes(fill=level_group),size=0.4,outlier.size=0.8,width=0.6,outlier.shape = NA) +
  scale_fill_manual(values = c("#f6f6f6","#8d5108")) + 
  theme_few() +
  labs(x = "", y = "Gene age", title = "", fill = "level")  +  
  scale_y_continuous(limits = c(-2, 2)) +          
  theme(axis.title = element_text(size=11),
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10,angle = 45,vjust = 1,hjust = 1),
        plot.title = element_text(hjust = 0.5,size=12),
        axis.ticks.length=unit(0.1,'cm')) + 
  stat_pvalue_manual(stat_wilcox.test, 
                     label = 'p.signif',
                     tip.length = 0.001,
                     hide.ns = TRUE)

stat_wilcox <- rstatix::wilcox_test(group_by(file_final1, cell_lines),gene_age ~ level)  
stat_wilcox <- rstatix::add_significance(stat_wilcox, "p")   
stat_wilcox.test <-  rstatix::add_xy_position(stat_wilcox, x = "cell_lines")
ggplot(file_final1, aes(x=cell_lines, y=gene_age)) +
  geom_boxplot(aes(fill=level),size=0.4,outlier.size=0.8,width=0.6,outlier.shape = NA) +
  scale_fill_manual(values = c("#f3e3b7","#d6b66a","#af6e22")) + 
  theme_few() +
  labs(x = "", y = "Gene age", title = "") +  
  scale_y_continuous(limits = c(-2, 2)) +          
  theme(axis.title = element_text(size=11),
        axis.text.y=element_text(size=10),
        axis.text.x=element_text(size=10,angle = 45,vjust = 1,hjust = 1),
        axis.ticks.length=unit(0.1,'cm')) + 
  stat_pvalue_manual(stat_wilcox.test, 
                     label = 'p.signif',
                     tip.length = 0.001,
                     hide.ns = TRUE)


###**TFCR比较 ####################
filename <- c("GM12878","HMEC","HUVEC","IMR90","K562","NHEK","HCT116")
chroms.vector <- c(paste0("chr",c(1:22,"X","Y")))
file_final <- NULL
for (f in filename){
  tfcr <- fread(paste0("F:/project1/9_new_results/0.TFCR/",f,"/TFCR_annotation_group.txt"))
  colnames(tfcr)[7] <- "complexity"
  tfcr <- tfcr[,c(1:4,7,8,20)]
  b <- fread(paste0("F:/project1/9_new_results/0.TAD/boundary/",f,"/boundary_level_all"))
  t_b <- lapply(chroms.vector, FUN=function(temp.chrom){
    temp.tfcr <- tfcr[seqnames==temp.chrom]
    temp.b <- b[V1==temp.chrom]
    aa=temp.tfcr[, {
      temp.b[V2 <= tfcr.end & V3 >= tfcr.start][,paste(V4)]
    },list(seqnames,tfcr.start=start,tfcr.end=end,width,complexity,annotation,group)]
  })
  tfcr_b <- do.call("rbind", t_b)
  data <- tfcr_b %>% group_by(seqnames,tfcr.start,tfcr.end) %>% filter(V1 == max(V1))  
  colnames(data)[colnames(data) %in% c("tfcr.start","tfcr.end","V1")] <- c("start","end","level")
  data1 <- merge(tfcr,data,all.x=T)
  data1$level[is.na(data1$level)] <- 0
  data1$cell_lines <- f
  file_final <- rbind(file_final,data1)
}
file_final1 <- file_final
file_final1$level[which(file_final1$level > 2)] <- "3+"
file_final1$level <- factor(file_final1$level,levels=c("0","1","2","3+"))
## complexity
stat_wilcox <- rstatix::wilcox_test(group_by(file_final1, cell_lines),complexity ~ level) 
stat_wilcox <- rstatix::add_significance(stat_wilcox, "p")  
stat_wilcox.test <-  rstatix::add_xy_position(stat_wilcox, x = "cell_lines")
ggplot(file_final1, aes(x=cell_lines, y=complexity)) +
  geom_boxplot(aes(fill=level),size=0.4,outlier.size=0.8,width=0.6,outlier.shape = NA) +
  scale_fill_manual(values = c("#f6f6f6","#f3e3b7","#d6b66a","#af6e22")) + 
  theme_few() +
  labs(x = "", y = "TFCR complexity", title = "") +     
  theme(axis.title = element_text(size=12),
        axis.text.y=element_text(size=10),
        axis.text.x=element_text(size=11),
        plot.title = element_text(hjust = 0.5,size=12),
        axis.ticks.length=unit(0.1,'cm')) + 
  stat_pvalue_manual(stat_wilcox.test, 
                     label = 'p.signif',
                     tip.length = 0.001,
                     step.increase = 0.06,
                     y.position = 90,
                     hide.ns = TRUE)

ggline(file_final1,x="level",y="complexity",
       color = "cell_lines",
       add = c("mean"),
       size = 0.8)+ 
  scale_color_manual(values = c("#8dd3c7","#fdb461","#e78ac4","#fb8073","#80b1d3","#e5c494","#999999")) + 
  labs(x = "Boundary level", y = "Mean TC value", title = "") +  
  theme_bw() +
  theme(axis.title = element_text(size=12),
        axis.text.y=element_text(size=11),
        axis.text.x=element_text(size=11),
        axis.ticks.length=unit(0.1,'cm'),
        legend.title = element_blank(),
        legend.position = "right") 

## width
stat_wilcox <- rstatix::wilcox_test(group_by(file_final1, cell_lines),width ~ level) 
stat_wilcox <- rstatix::add_significance(stat_wilcox, "p")  
stat_wilcox.test <-  rstatix::add_xy_position(stat_wilcox, x = "cell_lines")
ggplot(file_final1, aes(x=cell_lines, y=width)) +
  geom_boxplot(aes(fill=level),size=0.4,outlier.size=0.8,width=0.6,outlier.shape = NA) +
  scale_fill_manual(values = c("#f6f6f6","#f3e3b7","#d6b66a","#af6e22")) + 
  theme_few() +
  labs(x = "", y = "TFCR width", title = "") +    
  theme(axis.title = element_text(size=12),
        axis.text.y=element_text(size=10),
        axis.text.x=element_text(size=11),
        plot.title = element_text(hjust = 0.5,size=12),
        axis.ticks.length=unit(0.1,'cm')) + 
  stat_pvalue_manual(stat_wilcox.test, 
                     label = 'p.signif',
                     tip.length = 0.001,
                     step.increase = 0.06,
                     hide.ns = TRUE)
ggline(file_final1,x="level",y="width",
       color = "cell_lines",
       add = c("mean"),
       size = 0.8)+ 
  scale_color_manual(values = c("#8dd3c7","#fdb461","#e78ac4","#fb8073","#80b1d3","#e5c494","#999999")) + 
  labs(x = "Boundary level", y = "Mean TFCR width", title = "") +  
  theme_bw() +
  theme(axis.title = element_text(size=12),
        axis.text.y=element_text(size=11),
        axis.text.x=element_text(size=11),
        axis.ticks.length=unit(0.1,'cm'),
        legend.title = element_blank(),
        legend.position = "right") 

## annotation
file_final1$annotation <- str_replace(file_final1$annotation, "Exon.*", "Exon")
file_final1$annotation <- str_replace(file_final1$annotation, "Intron.*", "Intron")
file_final1$annotation <- str_replace(file_final1$annotation, "Promoter.*", "Promoter")
file_final1$annotation <- str_replace(file_final1$annotation, "Downstream.*", "Downstream")
file_final1$annotation <- str_replace(file_final1$annotation, "Distal.*", "Distal")
file_final1$annotation <- str_replace(file_final1$annotation, ".*UTR", "UTR")
file_final1$annotation <- factor(file_final1$annotation,levels=c("Distal","Intron","Exon","UTR","Downstream","Promoter"))
file_final1$number <- 1
ggplot(file_final1,aes(level,number,fill=annotation))+
  geom_bar(stat="identity",position="fill")+ 
  scale_fill_manual(values = c("#61cc5f","#31b57b","#1e948c","#2c738d","#3c508b","#443b84"))+
  scale_y_continuous(expand = expansion(mult=c(0.01,0.02)),labels = scales::percent_format())+
  labs(x="Boundary level",y="Genome fraction",fill=" ",title="")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.1,'cm'))+
  guides(fill=guide_legend(title=NULL))+
  facet_wrap(. ~ cell_lines,ncol = 3,nrow = 3)

