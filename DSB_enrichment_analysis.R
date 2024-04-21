library(data.table)
library(ggplot2)
library(tidyverse)
library(plyr)
library(ggthemes)
library(pacman)
library(ggpubr)
###HCT116和NHEK中DSB与边界的overlap情况
filename <- c("HCT116","NHEK (BLESS)","NHEK (DSBCapture)")
chroms.vector <- c(paste0("chr",c(1:22,"X","Y")))
file_final <- NULL
for (f in filename){
  boundary <- fread(paste0("F:/project1/9_new_results/0.TAD/boundary/",f,"/boundary_level_3"))
  genome <- fread("F:/project1/9_new_results/0.genome/10000_abs.bed")
  a <- boundary[,1:3]
  b <- genome[,1:3]
  diff <- setDT(b)[!a, on = names(b)]
  diff_sample <- diff[sample(1:nrow(diff),nrow(boundary),replace = FALSE),]
  diff_sample$V4 <- 0
  boundary <- rbind(boundary,diff_sample)
  dsb <- fread(paste0("F:/project1/9_new_results/0.DSB/select/",f,".narrowPeak"))
  dsb$number <- 1
  b_d <- lapply(chroms.vector, FUN=function(temp.chrom){
    temp.boundary <- boundary[V1==temp.chrom]
    temp.dsb <- dsb[V1==temp.chrom]
    aa=temp.boundary[, {
      temp.dsb[V2 <= boundary.V3 & V3 >= boundary.V2][,sum(number)]
    },list(V1,boundary.V2=V2,boundary.V3=V3,V4)]
  })
  boundary_dsb <- do.call("rbind", b_d)  
  colnames(boundary_dsb) <- c("seqnames","start","end","level","DSB_number")
  boundary_dsb$cell_line <- f
  file_final <- rbind(file_final,boundary_dsb)
}
# 存在DSB的边界的比例的比较
file_final$DSB_group <- gsub(pattern = "0",replacement = "no DSB",x = file_final$DSB_number)
file_final$DSB_group[which(file_final$DSB_group != "no DSB")] <- "DSB"
file_final$number <- 1
file_final$DSB_group <- factor(file_final$DSB_group,levels=c("no DSB","DSB"))
ggplot(file_final,aes(level,number,fill=DSB_group))+
  geom_bar(stat="identity",position="fill")+ 
  scale_fill_manual(values = c("#4da0a0","#9b3a74"))+
  scale_y_continuous(expand = expansion(mult=c(0.01,0.02)))+
  labs(x="Boundary level",y="Fraction",fill=" ",title="")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.1,'cm'),
        axis.text=element_text(size=10,face = "plain"),
        axis.title = element_text(size=12),
        legend.title = element_blank(),            
        legend.text = element_text(size = 10, face = "plain"),     
        legend.position = "right",         
        legend.key.size=unit(0.4,"cm"))+
  guides(fill=guide_legend(title=NULL))+
  facet_wrap(. ~ cell_line,ncol = 3,nrow = 1)  

compare_means(DSB_number~level,file_final,group.by = "cell_line")

# 存在DSB数目的均值的比较
library(Rmisc)
test <- summarySE(file_final, measurevar = "DSB_number", groupvars = c("level", "cell_line"))  # summarySE函数提供了标准差、标准误以及 95% 的置信区间
dev.new()
ggplot(test, aes(x = cell_line, y = DSB_number, colour = level)) + 
  geom_errorbar(aes(ymin = DSB_number - ci, ymax = DSB_number + ci), 
                position = position_dodge(0.7), width = 0.2, size = 0.8) + 
  geom_point(position = position_dodge(0.7),size=2) +
  scale_color_manual(values = c("#d0d0d0","#f3e3b7","#d6b66a","#af6e22")) + 
  theme_bw() +
  labs(x = "", y = "Mean number of DSBs per boundary", title = "",color="level") +  
  scale_y_continuous(expand = c(0,0),limits = c(0,1.1)) +
  theme(axis.title = element_text(size=12),
        axis.text.y=element_text(size=10),
        axis.text.x=element_text(size=11),
        plot.title = element_text(hjust = 0.5,size=12),
        axis.ticks.length=unit(0.1,'cm')) 

## HCT116
# bigwig信号
bigwig <- fread("F:/project1/9_new_results/0.genome/DSB_signal/HCT116/10000/readCounts.tab")
colnames(genome) <- c("seqnames","start","end","signal")
colnames(bigwig) <- c("V1","V2","V3","signal")
boundary <- fread("F:/project1/9_new_results/0.TAD/boundary/HCT116/boundary_level")
df <- merge(boundary,bigwig,by=c("V1","V2","V3"),all.x=T)
my_comparisons <- list(c("1", "2"),
                       c("2", "3"),
                       c("3", "4+"))
ggplot(df, aes(x=V4, y=signal)) +
  stat_boxplot(geom ="errorbar",width=0.3) +
  geom_boxplot(aes(fill=V4),size=0.4,outlier.size=0.8,width=0.6,outlier.shape = NA) +
  scale_fill_manual(values = c("#f3e3b7","#d6b66a","#af6e22","#8d5108")) + 
  theme_few() +
  labs(x = "boundary level", y = "DSB BigWig signal", title = "") +  
  guides(fill = "none") +     
  scale_y_continuous (limits = c (0, 6))+ 
  theme(axis.title = element_text(size=12),
        axis.text.y=element_text(size=10),
        axis.text.x=element_text(size=11),
        plot.title = element_text(hjust = 0.5,size=12),
        axis.ticks.length=unit(0.1,'cm')) + 
  stat_compare_means(comparisons=my_comparisons,
                     aes(label=..p.signif..), label.x = 1.5, label.y = 5)          
aggregate(signal ~ V4, data = df, median)

#方差分析（ANOVA）
aov1 <- aov(signal ~ V4, df)  
aov1
summary(aov1)

df$V4 <- factor(df$V4,levels=c("1","2","3","4+"))
ggline(df,x="V4",y="signal",
       color = "V4",
       add = c("mean_se"),
       size = 0.8)+ 
  scale_color_manual(values = c("#f3e3b7","#d6b66a","#af6e22","#8d5108")) + 
  labs(x = "Boundary level", y = "DSB BigWig signal (mean.±se.)", title = "") +  
  theme_bw() +
  theme(axis.title = element_text(size=12),
        axis.text.y=element_text(size=11),
        axis.text.x=element_text(size=11),
        axis.ticks.length=unit(0.1,'cm'),
        legend.title = element_blank(),
        legend.position = "right") 


##* 边界分有无TFCR比较DSB #####################
filename <- c("HCT116","NHEK(BLESS)","NHEK(DSBCapture)")
chroms.vector <- c(paste0("chr",c(1:22,"X","Y")))
file_final <- NULL
file_percent <- NULL
for (f in filename){
  boundary <- fread(paste0("F:/project1/9_new_results/0.TAD/boundary/",f,"/boundary_level_3"))
  dsb <- fread(paste0("F:/project1/9_new_results/0.DSB/select/",f,".narrowPeak"))
  dsb$number <- 1
  b_d <- lapply(chroms.vector, FUN=function(temp.chrom){
    temp.boundary <- boundary[V1==temp.chrom]
    temp.dsb <- dsb[V1==temp.chrom]
    aa=temp.boundary[, {
      temp.dsb[V2 <= boundary.V3 & V3 >= boundary.V2][,sum(number)]
    },list(V1,boundary.V2=V2,boundary.V3=V3,V4)]
  })
  boundary_dsb <- do.call("rbind", b_d)  
  colnames(boundary_dsb) <- c("seqnames","start","end","level","DSB_number")
  tfcr <- fread(paste0("F:/project1/9_new_results/0.TFCR/",f,"/TFCR_annotation_group.txt"))
  b_t <- lapply(chroms.vector, FUN=function(temp.chrom){
    temp.boundary <- boundary_dsb[seqnames==temp.chrom]
    temp.tfcr <- tfcr[seqnames==temp.chrom]
    aa=temp.boundary[, {
      temp.tfcr[start <= boundary_dsb.end & end >= boundary_dsb.start][,paste(group)]
    },list(seqnames,boundary_dsb.start=start,boundary_dsb.end=end,level,DSB_number)]
  })
  df <- do.call("rbind", b_t)  
  colnames(df)[2:3] <- c("start","end")
  df$TC_group <- as.numeric(gsub(pattern = "TC",replacement = "",x = df$V1))
  data <- df %>% 
    group_by(seqnames,start,end) %>% 
    filter(TC_group == max(TC_group)) %>%   
    unique() 
  merge_data <- merge(boundary_dsb,data,
                      by=c("seqnames","start","end","level","DSB_number"),
                      all.x=T)
  colnames(merge_data)[6] <- "TFCR"
  merge_data$TC_group[is.na(merge_data$TC_group)] <- "boundary_no_TFCR"
  merge_data$cell_line <- f
  merge_data$TC_group[which(merge_data$TC_group == 9)] <- "boundary_TC9"
  merge_data$TC_group[which(merge_data$TC_group != "boundary_TC9" & merge_data$TC_group != "boundary_no_TFCR")] <- "boundary_TC0-8"
  merge_data$DSB_group <- gsub(pattern = "0",replacement = "no DSB",x = merge_data$DSB_number)
  merge_data$DSB_group[which(merge_data$DSB_group != "no DSB")] <- "DSB"
  file_final <- rbind(file_final,merge_data)
  a <- data.frame(table(merge_data$TC_group,merge_data$DSB_group))
  a <- ddply(a,.(Var1),transform,percent=Freq/sum(Freq)*100) 
  a$cell_line <- f
  file_percent <- rbind(file_percent,a)
}
file_final$TC_group <- factor(file_final$TC_group,levels=c("boundary_TC9","boundary_TC0-8","boundary_no_TFCR"))
file_final$DSB_group <- factor(file_final$DSB_group,levels=c("no DSB","DSB"))
file_final$number <- 1
ggplot(file_final,aes(TC_group,number,fill=DSB_group))+
  geom_bar(stat="identity",position="fill",width = 0.7)+ 
  scale_fill_manual(values = c("#4da0a0","#9b3a74"))+
  scale_y_continuous(expand = expansion(mult=c(0.01,0.02)))+
  labs(x="Boundary level",y="Fraction",fill=" ",title="")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.1,'cm'),
        axis.text=element_text(size=10,face = "plain"),
        axis.title = element_text(size=12),
        legend.title = element_blank(),            
        legend.text = element_text(size = 10, face = "plain"),     
        legend.position = "right",         
        legend.key.size=unit(0.4,"cm"))+
  guides(fill=guide_legend(title=NULL))+
  facet_wrap(. ~ cell_line,ncol = 3,nrow = 1)  

file_dsb <- subset(file_percent,Var2=="DSB")
file_dsb$Var1 <- factor(file_dsb$Var1,levels=c("boundary_TC9","boundary_TC0-8","boundary_no_TFCR"))
dev.new()
ggplot(file_dsb,aes(cell_line,percent,fill=Var1))+
  geom_bar(stat="identity",position="dodge",width = 0.8)+ 
  scale_fill_manual(values = c("#08519c","#73b2d7","#cfe1ef"))+
  scale_y_continuous(expand = c(0,0),limits = c(0,100))+
  labs(x=" ",y="Fraction of DSBs (%)",fill=" ",title="")+
  theme_classic()+
  theme(axis.ticks.length=unit(0.1,'cm'),
        axis.text=element_text(size=10,face = "plain"),
        axis.title = element_text(size=12),
        legend.title = element_blank(),            
        legend.text = element_text(size = 10, face = "plain"),     
        legend.position = "right",         
        legend.key.size=unit(0.4,"cm"))+
  guides(fill=guide_legend(title=NULL)) 


## * TFCR分level比较存在DSB的比例---------------------
filename <- c("HCT116","NHEK(BLESS)","NHEK(DSBCapture)")
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
  dsb <- fread(paste0("F:/project1/9_new_results/0.DSB/select/",f,".narrowPeak"))
  dsb$number <- 1
  t_d <- lapply(chroms.vector, FUN=function(temp.chrom){
    temp.tfcr <- data1[seqnames==temp.chrom]
    temp.dsb <- dsb[V1==temp.chrom]
    aa=temp.tfcr[, {
      temp.dsb[V2 <= data1.end & V3 >= data1.start][,sum(number)]
    },list(seqnames,data1.start=start,data1.end=end,width,complexity,annotation,group,level)]
  })
  result <- do.call("rbind", t_d) 
  colnames(result)[colnames(result) %in% c("V1")] <- "DSB_number"
  result$DSB_density <- result$DSB_number/result$width
  result$cell_line <- f
  file_final <- rbind(file_final,result)
}
file_final1 <- file_final
file_final$level[which(file_final$level > 2)] <- "3+"
file_final$level <- factor(file_final$level,levels=c("0","1","2","3+"))
file_final$DSB_group <- gsub(pattern = "0",replacement = "no DSB",x = file_final$DSB_number)
file_final$DSB_group[which(file_final$DSB_group != "no DSB")] <- "DSB"
file_final$number <- 1
file_final$DSB_group <- factor(file_final$DSB_group,levels=c("no DSB","DSB"))
ggplot(file_final,aes(level,number,fill=DSB_group))+
  geom_bar(stat="identity",position="fill")+ 
  scale_fill_manual(values = c("#4da0a0","#9b3a74"))+
  scale_y_continuous(expand = expansion(mult=c(0.01,0.02)))+
  labs(x="Boundary level",y="Fraction",fill=" ",title="")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.1,'cm'),
        axis.text=element_text(size=10,face = "plain"),
        axis.title = element_text(size=12),
        legend.title = element_blank(),            
        legend.text = element_text(size = 10, face = "plain"),     
        legend.position = "right",         
        legend.key.size=unit(0.4,"cm"))+
  guides(fill=guide_legend(title=NULL))+
  facet_wrap(. ~ cell_line,ncol = 3,nrow = 1)  

ggline(file_final1,x="level",y="DSB_density",
       color = "cell_line",
       add = c("mean"),
       size = 0.8)+ 
  scale_color_manual(values = c("#8dd3c7","#fdb461","#e78ac4")) + 
  labs(x = "Boundary level", y = "Mean DSB density", title = "") +  
  theme_bw() +
  theme(axis.title = element_text(size=12),
        axis.text.y=element_text(size=11),
        axis.text.x=element_text(size=11),
        axis.ticks.length=unit(0.1,'cm'),
        legend.title = element_blank(),
        legend.position = "right") 
