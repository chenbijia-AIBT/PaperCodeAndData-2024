library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
setwd("F:/project1/9_new_results/0.TFCR/HCT116/0.find_TFCR_1.4/")
peak <- readPeakFile("f_tfbed_short/fimo.txt")
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
peakAnno <- annotatePeak(peak,tssRegion = c(-2000,2000),TxDb = txdb,annoDb = "org.Hs.eg.db")
write.table(as.data.frame(peakAnno),"5_annotation/TFCR_annotation.txt",sep="\t",row.names = F,quote = F)

