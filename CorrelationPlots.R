library(GenomicRanges)
library(dplyr)
setwd("C:/Users/egj19/OneDrive - McGill University/Desktop/Analysis/peakcalls")
#correlation plots using RPKMs normalized to inputs. Load Data that is clean
data=read.table("Parental_allmarks_chip_input.txt",header=T)
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
data=data[complete.cases(data),]
data=data[is.finite(data$H2AK119ub_G477),]
data=data[is.finite(data$H2AK119ub_pcGBM2),]
data=data[is.finite(data$CBX2_DIPG13),]
data=data[is.finite(data$H3K27me1_DIPG13),]
#correlation matrix of X Parental
data$RING1b=data$RING1b_BT245
data$SUZ12=data$SUZ12_BT245
data$CBX2=data$CBX2_BT245
data$H3K27me3=data$H3K27me3_BT245
data$H2AK119ub=data$H2AK119ub_BT245
data$H3K27me2=data$H3K27me2_BT245
data$H3K27me1=data$H3K27me1_BT245
data$H3K27ac=data$H3K27ac_BT245
data$H3K36me2=data$H3K36me2_BT245
data$H3K4me3=data$H3K4me3_BT245
data$H3K4me1=data$H3K4me1_BT245
#data$Pooled_H3K27me2_3=data$H3K27me3_DIPG13+data$H3K27me2_DIPG13

names(data)
#load package
library(ggcorrplot)
# Compute a correlation matrix
corr <- round(cor(data[c(118:128)]), 3)
head(corr[, 1:6])
#compute p-value matrix
p.mat <- cor_pmat(data[c(78:88)])
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
  scale_fill_gradient2(low = "royalblue2", high = "red1", mid = "snow1", 
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

ggsave("CorrelationMatrix_allmarks_BT245_NormalizedToInput.tiff", units="in", width=8, height=7, dpi=1200, compression = 'lzw')

#ggplot correlation matrix for all marks in ko 
#import dataset with read values normalized to total read count and filtered for 3 reads
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

data2$RING1b=data2$RING1b_DIPG13
data2$SUZ12=data2$SUZ12_DIPG13
data2$CBX2=data2$CBX2_DIPG13
data2$H3K27me3=data2$H3K27me3_DIPG13
data2$H2AK119ub=data2$H2AK119ub_DIPG13
data2$H3K27me2=data2$H3K27me2_DIPG13
data2$H3K27me1=data2$H3K27me1_DIPG13
data2$H3K27ac=data2$H3K27ac_DIPG13
data2$H3K36me2=data2$H3K36me2_DIPG13
data2$H3K4me3=data2$H3K4me3_DIPG13
#data2$Pooled_H3K27me2_3=data2$H3K27me2_BT245+data2$H3K27me3_BT245
names(data2)
#load package
library(ggcorrplot)
# Compute a correlation matrix
corr <- round(cor(data2[c(62:71)]), 3)
head(corr[, 1:6])
#compute p-value matrix
p.mat <- cor_pmat(data2[c(38:43)])
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
  scale_fill_gradient2(low = "royalblue2", high = "red1", mid = "snow1", 
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


ggsave("CorrelationMatrix_allmarks_DIPG13_K27MKO_NormalizedToInput.tiff", units="in", width=8, height=7, dpi=1200, compression = 'lzw')

#run limma on RING1b
library("limma")
ring1b_rpkms=right_join(data,data2,by="Probe")
rownames(ring1b_rpkms)=ring1b_rpkms[,1]
ring1b_rpkms=as.matrix(ring1b_rpkms[c(50,95)])
names(ring1b_rpkms)=paste(c("K27M","K27M KO"))

read.ilmn(ring1b_rpkms[1,],ring1b_rpkms[2,],path=NULL,ctrlpath = NULL)

test=aov(ring1b_rpkms[1,]~ring1b_rpkms[2,])

#boxplot values for H2AK119ub
median(log2(data$H2AK119ub_BT245*1.02130055+1))/1.058
median(log2(data$H2AK119ub_DIPG13*1.269379747+1))
median(log2(data2$H2AK119ub_BT245*0.948839873+1))/1.138
median(log2(data2$H2AK119ub_DIPG13*1.280060624+1))

df <- data.frame(Line=c("BT245 Parental","BT245 K27M KO","DIPGXIII Parental", "DIPGXIII K27M KO"),level=c(0.8767348,1,0.8111081,1),Condition=c("K27M","K27M KO","K27M","K27M KO"))


ggplot(data=df, aes(x=Line, y=level, fill=Condition)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_minimal()+
  ylab("Relative H2AK119ub1 level")+
  xlab("Cell Line")+
  scale_fill_manual(values=c("firebrick1","dodgerblue2","firebrick1","dodgerblue2"))

ggsave("H2AK119ub_levels_Barplot_CHIP.tiff", units="in", width=8, height=7, dpi=1200, compression = 'lzw')

#extract bins at promoters
library(GenomicRanges)
data.range=data[c(3,4,5)]
data.range=makeGRangesFromDataFrame(data.range)
#load annotated bins from UCSC
data.range.annotation=read.table("data.range.annotation.bed",header=F)
data.range.annotation=data.range.annotation[data.range.annotation$V4==c("1_Active_Promoter",'2_Weak_Promoter',"3_Poised_Promoter"),]
data.range.annotation=makeGRangesFromDataFrame(data.range.annotation,start.field = "V2",seqnames.field = "V1",end.field = "V3")
#create intersections
data.range.1=subsetByOverlaps(data.range,data.range.annotation)
data.range.1=as.data.frame(data.range.1)
data.range.1$link=paste(data.range.1$seqnames,data.range.1$start,data.range.1$end)
data$link=paste(data$Chromosome,data$Start,data$End)

Promoter_bins_Parental_allmarks=left_join(data.range.1,data)
Promoter_bins_Parental_allmarks=Promoter_bins_Parental_allmarks[complete.cases(Promoter_bins_Parental_allmarks),]
names(Promoter_bins_Parental_allmarks)
Promoter_bins_Parental_allmarks=Promoter_bins_Parental_allmarks[c(9,10,11,54:72)]
write.table(Promoter_bins_Parental_allmarks,"Promoter_bins_Parental_allmarks.txt",sep = " ")
