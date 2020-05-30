#theory of super-recruitement
#import spread regions obtained from chase clustering
setwd("C:/Users/egj19/OneDrive - McGill University/Desktop/Analysis")
library("rtracklayer")
b.spread=import.gff("bt245_k27me3_spread_regions.gff")
d.spread=import.gff("dipg13_k27me3_spread_regions.gff")
b.spread=makeGRangesFromDataFrame(b.spread)
d.spread=makeGRangesFromDataFrame(d.spread)

spread=intersect(b.spread,d.spread)
export.bed(spread,"k27me3_spread_regions.bed")

#ring1b counts from these regions using seqmonk
setwd("C:/Users/egj19/OneDrive - McGill University/Desktop/Analysis/peakcalls")
data=read.table("seqmonk20.txt",header=T)
names(data)
data$BT245=data[5]/data[8]
data$BT245_KO=data[7]/data[6]
data$DIPG13=data[9]/data[10]
data$DIPG13_KO=data[11]/data[12]
data1=data[c(8,9,10,11)]
names(data1)=paste(c("K27M",'KO',"K27M","KO"))
#run Limma
data1 <- normalizeBetweenArrays(data1, method="cyclicloess")
lm_p <- NULL
y <- new("EList", list(E=data1))
design <- model.matrix(~factor(colnames(data1)))
fit <- lmFit(y, design)
fit <- eBayes(fit)
res <- topTable(fit, coef=2, sort.by="none", n=Inf)
res$row=rownames(res)
lm_p <- rbind (row, res$p.value)
results <- cbind (data, res$P.Value)

#correlation plots on regions with spread k27me3
setwd("C:/Users/egj19/OneDrive - McGill University/Desktop/Analysis/peakcalls")
data=read.table("Parental_allmarks_chip_input.txt",header=T)
data.spread=as.data.frame(data[c(1,2,3,4,5)])
data$RING1b_BT245=data$BT245.DMSO.2.Rx_XChIP_RING1B_condition3.sorted.dup.bam/data$BT245.DMSO.3.Rx_cells_ChIP1_Input_1.sorted.dup.bam
data$RING1b_DIPG13=data$DIPG13p14_RING1B.sorted.dup.bam/data$DIPG13.C14_Input.sorted.dup.bam
data$H3K27me3_BT245=(data$BT245.C24_H3K27me3_1.sorted.dup.bam/data$BT245.C24_Input_1.sorted.dup.bam)*1
data$H3K27me3_DIPG13=(data$DIPG13.C14_H3K27me3.sorted.dup.bam/data$DIPG13.C14_Input.sorted.dup.bam)*1
data$H3K27ac_BT245=(data$BT245.nko.c1p5.Rx_cells_ChIP1_H3K27ac_1.sorted.dup.bam/data$BT245.nko.c1p5.Rx_cells_ChIP1_Input_1.sorted.dup.bam)*1
data$H3K27ac_DIPG13=(data$DIPGXIIIp11_cells_ChIP1_H3K27ac_1.sorted.dup.bam/data$DIPG13p11_input.sorted.dup.bam)*1
data$H3K36me2_BT245=data$BT245.C24.Rx_cells_ChIP1_H3K36me2_1.sorted.dup.bam/data$BT245.C24_Input_1.sorted.dup.bam
data$H3K36me2_DIPG13=data$DIPGXIIIp11_cells_ChIP1_H3K36me2_1.sorted.dup.bam/data$DIPG13p11_input.sorted.dup.bam
data$H3K4me3_BT245=data$BT245.C22_cells_ChIP1_H3K4me3_1.sorted.dup.bam/data$BT245.C22_cells_ChIP1_input_1.sorted.dup.bam
data$H3K4me3_DIPG13=data$DIPGXIIIp11_cells_ChIP1_H3K4me3_1.sorted.dup.bam/data$DIPG13p11_input.sorted.dup.bam
data$H3K4me1_BT245=data$BT245.C22_cells_ChIP1_H3K4me1_1.sorted.dup.bam/data$BT245.C22_cells_ChIP1_input_1.sorted.dup.bam
data$H3K4me1_DIPG13=data$DIPGXIII.C12.Rx_cells_ChIP1_H3K4me1_1.sorted.dup.bam.x/data$DIPG13.C12_Input.sorted.dup.bam
data$H3K4me3_BT245_KO=data$BT245.ko.c2p8.Rx_XChIP_H3K4me3.sorted.dup.bam/data$BT245.ko.c2p3.Rx_cells_ChIP1_Input_1.sorted.dup.bam
data$H3K4me3_DIPGX13_KO=data$DIPGXIII.ko.c5p5.Rx_XChIP_H3K4me3.sorted.dup.bam/data$DIPG13.ko.c5p5_Input.sorted.dup.bam
data$CBX2_BT245=data$BT245.P47_XChIP_CBX2_c1.sorted.dup.bam/data$BT245.P47.RX_XCHIP_INPUT.sorted.dup.bam
data$CBX2_DIPG13=data$DIPG13_p73_XxChIP_CBX2_c2.sorted.dup.bam.x/data$DIPG13.P73_XxChIP_Input.sorted.dup.bam.x
data$H2AK119ub_BT245=data$BT245.P47.H2AK119ub.sorted.dup.bam/data$BT245.P47.RX_XCHIP_INPUT.sorted.dup.bam
data$H2AK119ub_DIPG13=data$DIPG13.P29_H2AK119UB.sorted.dup.bam.x/data$DIPG13.P29.RX_XCHIP_INPUT.sorted.dup.bam
data$SUZ12_BT245=data$BT245p11_SUZ12_1.sorted.dup.bam/data$DIPG13p11_input.sorted.dup.bam
data$SUZ12_DIPG13=data$DIPG13p14_SUZ12.sorted.dup.bam/data$DIPG13.pool_Input.sorted.dup.bam
data$H3K27me2_BT245=data$BT245.nko.c1p5.Rx_cells_ChIP1_H3K27me2.sorted.dup.bam/data$BT245.nko.c1p5_ChIP1_Input_1.sorted.dup.bam
data$H3K27me2_DIPG13=data$DIPGXIII.nko.c12p3.Rx_cells_ChIP1_H3K27me2_1.sorted.dup.bam/data$DIPGXIII.nko.c12p3.Rx_cells_ChIP1_Input_1.sorted.dup.bam
data$H3K27me1_BT245=data$BT245.nko.c1p5.Rx_cells_ChIP1_H3K27me1_1.sorted.dup.bam/data$BT245.nko.c1p5_ChIP1_Input_1.sorted.dup.bam
data$H3K27me1_DIPG13=data$DIPGXIII.nko.c12p3.Rx_cells_ChIP1_H3K27me1_1.sorted.dup.bam/data$DIPGXIII.nko.c12p3.Rx_cells_ChIP1_Input_1.sorted.dup.bam
data$H2AK119ub_G477=data$G477.P25_H2K119UB.sorted.dup.bam/data$G477.pool_Input.sorted.dup.bam
data$H2AK119ub_pcGBM2=data$pcGBM2.p27_H2AK119ub.sorted.dup.bam/data$pcGBM2.p27_Input.sorted.dup.bam
data$H2AK119ub_BT245_KO_2=data$BT245.C2P7.K27MKO_H2AK119UB.sorted.dup.bam.y/data$BT245.K27MKO.C2P7_INPUT.sorted.dup.bam.y
data$H2AK119ub_DIPG13_KO_2=data$DIPG13.KO.C5P5_H2K119UB.sorted.dup.bam/data$DIPG13.ko.c5p5_Input.sorted.dup.bam
data$H2AK119ub_DIPG13_KO_3=data$DIPG13.C2P3.K27MKO_H2AK119UB.sorted.dup.bam.y/data$DIPG13.K27MKO.C2P3_INPUT.sorted.dup.bam.y
data$H2AK119ub_DIPG13_2=data$DIPG13_P66_Rx_XChIP_H2AK119ub_80M.sorted.dup.bam/data$DIPG13_P66_Rx_XChIP_INPUT.sorted.dup.bam.y
data$RING1B=data$RING1b_BT245
data$H2AK119ub=data$H2AK119ub_BT245
data$CBX2=data$CBX2_BT245
data$SUZ12=data$SUZ12_BT245
data$H3K27me3=data$H3K27me3_BT245


data.spread=makeGRangesFromDataFrame(data.spread)
data1=subsetByOverlaps(data.spread,d.spread)
data1=as.data.frame(data1)
data1$probe=paste(data1$seqnames,data1$start,data1$end)
data$probe=paste(data$seqnames,data$start,data$end)
library(dplyr)
data1=left_join(data1,data,by="probe")
names(data1)
#load package
library(ggcorrplot)
# Compute a correlation matrix
corr <- round(cor(data1[c(124:128)]), 3)
head(corr[, 1:4])
#compute p-value matrix
p.mat <- cor_pmat(data1[c(39:43)])
head(p.mat[, 1:4])
# Get upper triangle of the correlation matrix
get_upper_tri <- function(corr){
  corr[lower.tri(corr)]<- NA
  return(corr)
}
upper_tri <- get_upper_tri(corr)
# Melt the correlation matrix
library(reshape2)
melted_cormat <- melt(upper_tri, na.rm = TRUE)
#reorder function
reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}
# Reorder the correlation matrix
cormat <- reorder_cormat(corr)
upper_tri <- get_upper_tri(corr)
# Melt the correlation matrix
melted_cormat <- melt(upper_tri, na.rm = TRUE)
#heatmap
library(ggplot2)
ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "dodgerblue2", high = "firebrick2", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1)) +
  ggtitle("BT245 Parental")+
  coord_fixed()+
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 4)+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    legend.direction = "horizontal",
    legend.justification = c(1, 0),
    legend.position = c(0.5, 0.8))+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))
ggsave("CorrelationMatrix_PRC_BT245_NormalizedToInput_Theory.tiff", units="in", width=8, height=7, dpi=1200, compression = 'lzw')

#ggplot correlation matrix for all marks in ko 
#import dataset with read values normalized to total read count and filtered for 3 reads
setwd("C:/Users/egj19/OneDrive - McGill University/Desktop/Analysis/peakcalls")
data2=read.table("KO_allmarks_chip_input.txt",header=T)
data2$RING1b_BT245=data2$BT245.ko.C2P8.Rx_XChIP_RING1B.sorted.dup.bam/data2$BT245.ko.c2p3.Rx_cells_ChIP1_Input_1.sorted.dup.bam
data2$RING1b_DIPG13=data2$DIPGXIII.ko.c5p6.7.Rx_XChIP_RING1B.sorted.dup.bam/data2$DIPGXIII.ko.c5p6.7.Rx_XChIP_Input.sorted.dup.bam
data2$H2AK119ub_BT245=data2$BT245.ko.C2P8_H2AK119UB.sorted.dup.bam/data2$BT245.ko.C2P8.Rx_XChIP_INPUT.sorted.dup.bam
data2$H2AK119ub_DIPG13=data2$DIPG13.ko.C5P33_H2AK119UB.sorted.dup.bam/data2$DIPG13.ko.C5P33.Rx_XChIP_INPUT.sorted.dup.bam
data2$H3K27me3_BT245=data2$BT245.ko.c4p6.Rx.C_cells_ChIP1_H3K27me3_1.sorted.dup.bam/data2$BT245.ko.c4p6.Rx_cells_ChIP1_Input_1.sorted.dup.bam
data2$H3K27me3_DIPG13=data2$DIPGXIII.ko.c10p5.Rx.C_cells_ChIP1_H3K27me3_1.sorted.dup.bam/data2$DIPGXIII.ko.c10p5.Rx_cells_ChIP1_Input_1.bam.sorted.dup.bam
data2$H3K27ac_BT245=data2$BT245.ko.c2p3.Rx_cells_ChIP1_H3K27ac_1.sorted.dup.bam/data2$BT245.ko.c2p3.Rx_cells_ChIP1_Input_1.sorted.dup.bam
data2$H3K27ac_DIPG13=data2$DIPGXIII.ko.c5p5.Rx_cells_ChIP1_H3K27ac_1.sorted.dup.bam.x/data2$DIPG13.ko.c5p5_Input.sorted.dup.bam
data2$H3K36me2_BT245=data2$BT245.ko.c4p6.Rx_cells_ChIP1_H3K36me2_1.sorted.dup.bam/data2$BT245.ko.c4p6.Rx_cells_ChIP1_Input_1.sorted.dup.bam
data2$H3K36me2_DIPG13=data2$DIPGXIII.ko.c10p5.Rx_cells_ChIP1_H3K36me2_1.sorted.dup.bam/data2$DIPGXIII.ko.c10p5.Rx_cells_ChIP1_Input_1.bam.sorted.dup.bam
data2$SUZ12_BT245=data2$BT245.ko.c4p9.Rx_cells_ChIP1_SUZ12_1.sorted.dup.bam/data2$BT245.ko.c4p9_Input_1.sorted.dup.bam.x
data2$SUZ12_DIPG13=data2$DIPGXIII.ko.c5p5.Rx_cells_ChIP1_SUZ12_1.sorted.dup.bam/data2$DIPG13.ko.c5p5_Input.sorted.dup.bam
data2$H3K27me2_BT245=data2$BT245.ko.c4p9.Rx_cells_ChIP1_H3K27me2_1.sorted.dup.bam/data2$BT245.ko.c4p9_Input_1.sorted.dup.bam.x
data2$H3K27me2_DIPG13=data2$DIPGXIII.ko.c5p5.Rx_cells_ChIP1_H3K27me2_1.sorted.dup.bam/data2$DIPGXIII.ko.c5p5.Rx_cells_ChIP1_Input_1.sorted.dup.bam
data2$H3K27me1_BT245=data2$BT245.ko.c4p9.Rx_cells_ChIP1_H3K27me1_1.sorted.dup.bam/data2$BT245.ko.c4p9_Input_1.sorted.dup.bam.y
data2$H3K27me1_DIPG13=data2$DIPGXIII.ko.c5p5.Rx_cells_ChIP1_H3K27me1_1.sorted.dup.bam/data2$DIPG13.ko.c5p5_Input.sorted.dup.bam
data2$H3K4me3_BT245=data2$BT245.ko.c2p8.Rx_XChIP_H3K4me3.sorted.dup.bam/data2$BT245.ko.c2p3.Rx_cells_ChIP1_Input_1.sorted.dup.bam
data2$H3K4me3_DIPG13=data2$DIPGXIII.ko.c5p5.Rx_XChIP_H3K4me3.sorted.dup.bam/data2$DIPG13.ko.c5p5_Input.sorted.dup.bam
data2$CBX2_BT245=data2$BT245.C2P8_XChIP_CBX2_c1.sorted.dup.bam/data2$BT245.ko.c2p3.Rx_cells_ChIP1_Input_1.sorted.dup.bam
data2$CBX2_DIPG13=data2$DIPG13.C5P21_XxChIP_CBX2_c2.sorted.dup.bam/data2$DIPG13.C5P21_XxChIP_input.sorted.dup.bam

data.spread=makeGRangesFromDataFrame(data.spread)
data3=subsetByOverlaps(data.spread,d.spread)
data3=as.data.frame(data3)
data3$probe=paste(data3$seqnames,data3$start,data3$end)
data2$probe=paste(data2$seqnames,data2$start,data2$end)


library(dplyr)
data3=left_join(data3,data2,by="probe")
data3=data3[complete.cases(data3),]
names(data3)

data3$RING1B=data3$RING1b_DIPG13
data3$H2AK119ub=data3$H2AK119ub_DIPG13
data3$CBX2=data3$CBX2_DIPG13
data3$SUZ12=data3$SUZ12_DIPG13
data3$H3K27me3=data3$H3K27me3_DIPG13

names(data3)
#load package
library(ggcorrplot)
# Compute a correlation matrix
corr <- round(cor(data3[c(67:71)]), 3)
head(corr[, 1:4])
#compute p-value matrix
p.mat <- cor_pmat(data3[c(32:35)])
head(p.mat[, 1:4])
# Get upper triangle of the correlation matrix
get_upper_tri <- function(corr){
  corr[lower.tri(corr)]<- NA
  return(corr)
}
upper_tri <- get_upper_tri(corr)
# Melt the correlation matrix
library(reshape2)
melted_cormat <- melt(upper_tri, na.rm = TRUE)
#reorder function
reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}
# Reorder the correlation matrix
cormat <- reorder_cormat(corr)
upper_tri <- get_upper_tri(corr)
# Melt the correlation matrix
melted_cormat <- melt(upper_tri, na.rm = TRUE)
#heatmap
library(ggplot2)
ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "dodgerblue2", high = "firebrick2", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1)) +
  ggtitle("DIPGXIII K27M KO")+
  coord_fixed()+
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 4)+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    legend.direction = "horizontal",
    legend.justification = c(1, 0),
    legend.position = c(0.5, 0.8))+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))
ggsave("CorrelationMatrix_PRC_DIPG13_KO_NormalizedToInput_Theory.tiff", units="in", width=8, height=7, dpi=1200, compression = 'lzw')

#are these regions associated with any DEG?
bt245.genes=read.table("seqmonk17.txt",header=T)
bt245.genes$probe=paste(bt245.genes$Chromosome,bt245.genes$Start,bt245.genes$End)
bt245.genes.spread=bt245.genes[c(2,3,4)]
bt245.genes.spread=makeGRangesFromDataFrame(bt245.genes.spread)
bt245.genes1=subsetByOverlaps(bt245.genes.spread,b.spread)
bt245.genes1=as.data.frame(bt245.genes1)
bt245.genes1$probe=paste(bt245.genes1$seqnames,bt245.genes1$start,bt245.genes1$end)
bt245.genes1=left_join(bt245.genes1,bt245.genes,by="probe")

#RING1b peaks within these regions
bt245.ring1b=subsetByOverlaps(BT245,b.spread)
bt245.suz12=subsetByOverlaps(b.suz12.k27m,spread)
a1=length(subsetByOverlaps(bt245.ring1b,bt245.suz12))/length(bt245.ring1b)

bt245.ko.ring1b=subsetByOverlaps(BT245.ko.1,b.spread)
bt245.ko.suz12=subsetByOverlaps(b.suz12.ko.2,b.spread)
a2=length(subsetByOverlaps(bt245.ko.ring1b,bt245.ko.suz12))/length(bt245.ko.ring1b)


#dipg13
#correlation plots on regions with spread k27me3
setwd("C:/Users/egj19/OneDrive - McGill University/Desktop/Analysis/peakcalls")
data=read.table("Parental_allmarks_chip_input.txt",header=T)
names(data)
data.spread=as.data.frame(data[c(2,3,4,5)])
data$RING1b=data$DIPG13p14_RING1B.sorted.dup.bam/data$DIPG13.C14_Input.sorted.dup.bam
data$SUZ12=data$DIPG13p14_SUZ12.sorted.dup.bam/data$DIPG13.pool_Input.sorted.dup.bam
data$H3K27me3=data$DIPGXIII.C12.Rx_cells_ChIP1_H3K27me3_1.sorted.dup.bam/data$DIPGXIII.C12.Rx_cells_ChIP1_Input_1.sorted.dup.bam
data$H2AK119ub=data$DIPG13.P29_H2AK119UB.sorted.dup.bam/data$DIPG13.P29.RX_XCHIP_INPUT.sorted.dup.bam

data.spread=makeGRangesFromDataFrame(data.spread)
data1=subsetByOverlaps(data.spread,d.spread)
data1=as.data.frame(data1)
data1$probe=paste(data1$seqnames,data1$start,data1$end)
data$probe=paste(data$Chromosome,data$Start,data$End)
library(dplyr)
data1=left_join(data1,data,by="probe")
names(data1)
library(ggcorrplot)
# Compute a correlation matrix
corr <- round(cor(data1[c(53:56)]), 3)
head(corr[, 1:4])
#compute p-value matrix
p.mat <- cor_pmat(data1[c(39:42)])
head(p.mat[, 1:4])
# Get upper triangle of the correlation matrix
get_upper_tri <- function(corr){
  corr[lower.tri(corr)]<- NA
  return(corr)
}
upper_tri <- get_upper_tri(corr)
# Melt the correlation matrix
library(reshape2)
melted_cormat <- melt(upper_tri, na.rm = TRUE)
#reorder function
reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}
# Reorder the correlation matrix
cormat <- reorder_cormat(corr)
upper_tri <- get_upper_tri(corr)
# Melt the correlation matrix
melted_cormat <- melt(upper_tri, na.rm = TRUE)
#heatmap
library(ggplot2)
ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ 
  ggtitle("DIPGXIII Parental")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1)) +
  coord_fixed()+
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 4)+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    legend.direction = "horizontal",
    legend.justification = c(1, 0),
    legend.position = c(0.5, 0.8))+
  
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))

ggsave("CorrelationPlot_DIPG13_Parental_Theory_normalizedtoinputs.tiff", units="in", width=5, height=4, dpi=600, compression = 'lzw')

#repeat for dipg13 ko
#import dataset with read values normalized to total read count and filtered for 3 reads
data2=read.table("seqmonk16.txt",header=T)
data2$RING1b=data2$DIPGXIII.ko.c5p6.7.Rx_XChIP_RING1B.sorted.dup.bam/data2$DIPGXIII.ko.c5p6.7.Rx_XChIP_Input.sorted.dup.bam
data2$SUZ12=data2$DIPGXIII.ko.c5p5.Rx_cells_ChIP1_SUZ12_1.sorted.dup.bam/data2$DIPG13.ko.c5p5_Input.sorted.dup.bam
data2$H3K27me3=data2$DIPGXIII.ko.c10p5.Rx.C_cells_ChIP1_H3K27me3_1.sorted.dup.bam/data2$DIPGXIII.ko.c10p5.Rx_cells_ChIP1_Input_1.bam.sorted.dup.bam
data2$H2AK119ub=data2$DIPG13.ko.C5P33_H2AK119UB.sorted.dup.bam/data2$DIPG13.ko.C5P33.Rx_XChIP_INPUT.sorted.dup.bam

data.spread=makeGRangesFromDataFrame(data.spread)
data4=subsetByOverlaps(data.spread,d.spread)
data4=as.data.frame(data4)
data4$probe=paste(data4$seqnames,data4$start,data4$end)
data2$probe=paste(data2$Chromosome,data2$Start,data2$End)

library(dplyr)
data4=left_join(data4,data2,by="probe")
data4=data4[complete.cases(data4),]
names(data4)


#load package
library(ggcorrplot)
# Compute a correlation matrix
corr <- round(cor(data4[c(32:35)]), 3)
head(corr[, 1:4])
#compute p-value matrix
p.mat <- cor_pmat(data3[c(32:35)])
head(p.mat[, 1:4])
# Get upper triangle of the correlation matrix
get_upper_tri <- function(corr){
  corr[lower.tri(corr)]<- NA
  return(corr)
}
upper_tri <- get_upper_tri(corr)
# Melt the correlation matrix
library(reshape2)
melted_cormat <- melt(upper_tri, na.rm = TRUE)
#reorder function
reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}
# Reorder the correlation matrix
cormat <- reorder_cormat(corr)
upper_tri <- get_upper_tri(corr)
# Melt the correlation matrix
melted_cormat <- melt(upper_tri, na.rm = TRUE)
#heatmap
library(ggplot2)
ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1)) +
  coord_fixed()+
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 4)+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    legend.direction = "horizontal",
    legend.justification = c(1, 0),
    legend.position = c(0.5, 0.8))+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))
ggsave("CorrelationPlot_DIPG13_KO_Theory_normalizedtoinputs.tiff", units="in", width=5, height=4, dpi=600, compression = 'lzw')


#Extract rign1b peaks at k27me3 cluster
setwd("C:/Users/egj19/OneDrive - McGill University/Desktop/Analysis")
b.chase=read.table("chase_bt245_k27me3_cluster.gff")
d.chase=read.table("chase_dipg13_k27me3_cluster.gff")
b.chase=makeGRangesFromDataFrame(b.chase,seqnames.field= "V1",start.field = "V2",end.field = "V3")
d.chase=makeGRangesFromDataFrame(d.chase,seqnames.field= "V1",start.field = "V2",end.field = "V3")

overlaps1=subsetByOverlaps(BT245,b.chase)
overlaps1
overlaps2=subsetByOverlaps(DIPG13,d.chase)
overlaps2

overlaps3=intersect(overlaps1,overlaps2)

export.bed(overlaps2,"test.bed")

#associate these peaks with genes
setwd("C:/Users/egj19/OneDrive - McGill University/Desktop/Analysis/peakcalls")

b.anno=read.csv("BT245-DMSO-2-Rx_XChIP_RING1B_condition3.annotated.csv")
b.anno$join=paste(b.anno$Chr,b.anno$Start,b.anno$End)
b.anno.1=makeGRangesFromDataFrame(b.anno)
b.anno.2=subsetByOverlaps(b.anno.1,overlaps1)
b.anno.2=as.data.frame(b.anno.2)
b.anno.2$join=paste(b.anno.2$seqnames,b.anno.2$start,b.anno.2$end)
b.anno.2=left_join(b.anno.2,b.anno)
b.anno.2=as.data.frame(b.anno.2[,c(1,2,3,4,10)])
colnames(b.anno.2)=c("Chromosome","Start","End","Width","Probe")
setwd("C:/Users/egj19/OneDrive - McGill University/Desktop/Analysis")
b.deg=read.table("BT245_DEG.txt",header=T)
b.anno.2=left_join(b.anno.2,b.deg,by="Probe")
b.anno.2=b.anno.2[complete.cases(b.anno.2),]

#141 of these genes are upregulated while 365 are downregulated
test=b.anno.2[b.anno.2$Log2_Fold_Change<c(0),]
test=makeGRangesFromDataFrame(test,seqnames.field = "Chromosome.x",start.field = "Start.x",end.field = "End.x")
export.bed(test,"test.bed")

#plot pichart of DEG numbers
library(ggplot2)
pie3 <- data.frame(Genes=c("Downregulated","Upregulated"),value=c(365,141))
library(dplyr)
# Basic piechart
pie3 <- pie3 %>% 
  arrange(desc(Genes)) %>%
  mutate(prop = value/ 506) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

ggplot(pie3, aes(x="", y=prop, fill=Genes)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  geom_text(aes(y = ypos, label = value), color = "white", size=5) +
  scale_fill_manual(values = c("red","#4E84C4"))+
  theme(legend.background = element_rect(size=2, linetype="solid"))
ggsave("Rplot.tiff", units="in", width=5, height=4, dpi=600, compression = 'lzw')

#volcano plot
names(b.anno.2)
threshold_OE <- b.anno.2$P.value<0.5
length(which(threshold_OE))
b.anno.2$threshold <- threshold_OE 

ggplot(b.anno.2) +
  geom_point(aes(x=Log2_Fold_Change, y=-log2(P.value), colour=threshold)) +
  ggtitle("RING1b") +
  xlab("log2 fold change") + 
  ylab("-log2 P-value") +
  #scale_y_continuous(limits = c(0,50)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))+theme_minimal()


#are there sites that lose ring1b due to k27me3 loss?
b.lostpeaks
b.lostk27me3.haifen=setdiff(b.k27me3.haifen.ko,subsetByOverlaps(b.k27me3.haifen.ko,intersect(b.k27me3.haifen.k27m,b.k27me3.haifen.ko)))
test=subsetByOverlaps(b.lostpeaks,b.lostk27me3.haifen)

export.bed(test,"test.bed")

#In Chase, how many bins of RING1b-K27me3 cluster overlaps
setwd("C:/Users/egj19/OneDrive - McGill University/Desktop/Analysis")
b245_k27me3_clusteronpeaks=import.gff("b245_k27me3_clusteronpeaks.gff")
bt245_k27ac_clusteronpeaks=import.gff("bt245_k27ac_clusteronpeaks.gff")
b245_KO_k27me3_clusteronpeaks=import.gff("b245_KO_k27me3_clusteronpeaks.gff")
b245_KO_k27ac_clusteronpeaks=import.gff("b245_KO_k27ac_clusteronpeaks.gff")
dipg13_k27ac_clusteronpeaks=import.gff("dipg13_k27ac_clusteronpeaks.gff")
dipg13_k27me3_clusteronpeaks=import.gff("dipg13_k27me3_clusteronpeaks.gff")
dipg13_KO_k27ac_clusteronpeaks=import.gff("dipg13_KO_k27ac_clusteronpeaks.gff")
dipg13_KO_k27me3_clusteronpeaks=import.gff("dipg13_KO_k27me3_clusteronpeaks.gff")

dipg13_KO_k27me3_clusteronpeaks=makeGRangesFromDataFrame(dipg13_KO_k27me3_clusteronpeaks)
export.bed(dipg13_KO_k27me3_clusteronpeaks,"dipg13_KO_k27me3_clusteronpeaks.bed")

#log2 figure
library(readxl)
cprc1.annotation=read_excel("cprc1.annotation.xlsx")
range=makeGRangesFromDataFrame(cprc1.annotation)
range1=makeGRangesFromDataFrame(data)
range2=makeGRangesFromDataFrame(data2)
range3=subsetByOverlaps(range1,range2)

hits1=as.list(findOverlaps(range,range3))
hits2=as.list(findOverlaps(range,range3))

cprc1.annotation$RING1b_BT245=as.numeric(lapply(extractList(data$RING1b_BT245,hits1),mean))
cprc1.annotation$RING1b_DIPG13=as.numeric(lapply(extractList(data$RING1b_DIPG13,hits1),mean))
cprc1.annotation$RING1b_BT245_KO=as.numeric(lapply(extractList(data2$RING1b_BT245_KO,hits2),mean))
cprc1.annotation$RING1b_DIPG13_KO=as.numeric(lapply(extractList(data2$RING1b_DIPG13_KO,hits2),mean))
#fix
cprc1.annotation=cprc1.annotation[complete.cases(cprc1.annotation),]

#combine log2 fold
library(dplyr)
setwd("C:/Users/egj19/OneDrive - McGill University/Desktop/Analysis/BEDs")
K27M_DEG=read.table("K27M_DEG.txt",header=T)
test=K27M_DEG[K27M_DEG$P.adj<0.1,]
cprc1.annotation=left_join(cprc1.annotation,test,by="Probe")
cprc1.annotation=cprc1.annotation[complete.cases(cprc1.annotation),]

cprc1.annotation$RING1b_K27M=(cprc1.annotation$RING1b_BT245+cprc1.annotation$RING1b_DIPG13)/2
cprc1.annotation$RING1b_KO=(cprc1.annotation$RING1b_BT245_KO+cprc1.annotation$RING1b_DIPG13_KO)/2
cprc1.annotation$RING1b_log2=log2(cprc1.annotation$RING1b_K27M/cprc1.annotation$RING1b_KO)


#plot scatter plot
library(ggplot2)
ggplot(cprc1.annotation, aes(x=RING1b_log2, y=Log2_FC)) +
  geom_point(size=2,color="#56B4E9")+geom_smooth(method=lm,  linetype="dashed",
                                                         color="#E69F00", fill="#999999")+theme_minimal()+
  xlab("Log2 RING1B RPKMs [K27M/K27M KO]")+ylab("Log2 Gene Expression Fold Change [K27M/K27M KO]")
ggsave("Scatterplot_RING1B_geneexpression_log2_DEG_cprc1.tiff", units="in", width=5, height=4, dpi=600, compression = 'lzw')


#pcgbm2 
test=read.table("seqmonk51.txt",header=T)
test$H3K27me3=test$pcGBM2.p31_H3K27me3.sorted.dup.bam/test$pcGBM2.p27.Rx_cells_ChIP1_Input_1.sorted.dup.bam
summary(test$H3K27me3)
test=test[test$H3K27me3>1.2,]
test=makeGRangesFromDataFrame(test)
#export.bed(test,"pcGBM2_H3K27me3_10Kb.bed")
