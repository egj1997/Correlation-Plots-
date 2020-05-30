#kmeans
library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(readxl)
#djihad heatmap
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
data$H3K4me1_DIPG13=data$DIPGXIII.C12.Rx_cells_ChIP1_H3K4me1_1.sorted.dup.bam/data$DIPGXIII.C12.Rx_cells_ChIP1_Input_1.sorted.dup.bam
data$CBX2_BT245=data$BT245.P47_XChIP_CBX2_c1.sorted.dup.bam/data$BT245.P47.RX_XCHIP_INPUT.sorted.dup.bam
data$H2AK119ub_BT245=data$BT245.P47.H2AK119ub.sorted.dup.bam/data$BT245.P47.RX_XCHIP_INPUT.sorted.dup.bam
data$H2AK119ub_DIPG13=data$DIPG13.P29_H2AK119UB.sorted.dup.bam/data$DIPG13.P29.RX_XCHIP_INPUT.sorted.dup.bam
data$SUZ12_BT245=data$BT245p11_SUZ12_1.sorted.dup.bam/data$DIPG13p11_input.sorted.dup.bam
data$SUZ12_DIPG13=data$DIPG13p14_SUZ12.sorted.dup.bam/data$DIPG13.pool_Input.sorted.dup.bam
data$H3K27me2_BT245=data$BT245.nko.c1p5.Rx_cells_ChIP1_H3K27me2.sorted.dup.bam/data$BT245.nko.c1p5_ChIP1_Input_1.sorted.dup.bam
data$H3K27me2_DIPG13=data$DIPGXIII.nko.c12p3.Rx_cells_ChIP1_H3K27me2_1.sorted.dup.bam/data$DIPGXIII.nko.c12p3.Rx_cells_ChIP1_Input_1.sorted.dup.bam
data$H3K27me1_BT245=data$BT245.nko.c1p5.Rx_cells_ChIP1_H3K27me1_1.sorted.dup.bam/data$BT245.nko.c1p5_ChIP1_Input_1.sorted.dup.bam
data$H3K27me1_DIPG13=data$DIPGXIII.nko.c12p3.Rx_cells_ChIP1_H3K27me1_1.sorted.dup.bam/data$DIPGXIII.nko.c12p3.Rx_cells_ChIP1_Input_1.sorted.dup.bam
#correlation matrix of X Parental
data$RING1b=data$RING1b_BT245
data$H2AK119ub=data$H2AK119ub_BT245
data$CBX2=data$CBX2_BT245
data$SUZ12=data$SUZ12_BT245
data$H3K27me3=data$H3K27me3_BT245
data$H3K27me2=data$H3K27me2_BT245
data$H3K27me1=data$H3K27me1_BT245
data$H3K27ac=data$H3K27ac_BT245
data$H3K4me1=data$H3K4me1_BT245
data$H3K4me3=data$H3K4me3_BT245
data$H3K36me2=data$H3K36me2_BT245
#Choose dataframe
names(data)
BT245.t=data[c(1:5,107:117)]
names(BT245.t)
#kmeans table
kmat=BT245.t[,c(6:16)]
#BT245 peaks
BT245=import.bed("BT245_RING1b_peaks_optimized.bed")
BT245_KO=import.bed("BT245_C2P8_RING1b_peaks_optimized.bed")
#find overlaps in table
BT245range=makeGRangesFromDataFrame(BT245.t,seqnames.field = "seqnames",end.field = "End",start.field="Start")
region=union(BT245,BT245_KO)
hits=subsetByOverlaps(BT245range,region)
hits=as.data.frame(hits)
hits$link=paste(hits$seqnames,hits$start,hits$end)
BT245.t$link=paste(BT245.t$seqnames,BT245.t$start,BT245.t$end)
hits=left_join(hits,BT245.t,by="link")
names(hits)
hits$Probe.x=rownames(hits)
#create kmeans table
kmat=hits[,c(12:22)]
rownames(kmat)=hits$Probe.x

cluster=kmeans(kmat, 2, iter.max = 10, nstart = 1,
       algorithm = c("Hartigan-Wong", "Lloyd", "Forgy",
                     "MacQueen"), trace=FALSE)
clusters.see=fitted(cluster,method = "classes")
clusters.see=as.data.frame(clusters.see)
clusters.see$Probe.x=rownames(clusters.see)
names(clusters.see)=paste(c("Cluster","Probe.x"))
#generate bed files from 5 clusters
clusters.see=left_join(clusters.see,hits,by="Probe.x")
cluster1=clusters.see[clusters.see$Cluster==1,]
cluster1=makeGRangesFromDataFrame(cluster1,seqnames.field = "seqnames.x",start.field = "start.x",end.field = "end.x")
cluster1=subsetByOverlaps(BT245,cluster1)
export.bed(cluster1,"cluster1.1.bed")

cluster2=clusters.see[clusters.see$Cluster==2,]
cluster2=makeGRangesFromDataFrame(cluster2,seqnames.field = "seqnames.x",start.field = "start.x",end.field = "end.x")
cluster2=subsetByOverlaps(BT245,cluster2)
export.bed(cluster2,"cluster2.1.bed")

cluster3=clusters.see[clusters.see$Cluster==3,]
cluster3=makeGRangesFromDataFrame(cluster3)
cluster3=subsetByOverlaps(BT245,cluster3)
export.bed(cluster3,"cluster3.bed")

cluster4=clusters.see[clusters.see$Cluster==4,]
cluster4=makeGRangesFromDataFrame(cluster4)
cluster4=subsetByOverlaps(BT245,cluster4)
export.bed(cluster4,"cluster4.bed")

cluster5=clusters.see[clusters.see$Cluster==5,]
cluster5=makeGRangesFromDataFrame(cluster5)
cluster5=subsetByOverlaps(BT245,cluster5)
export.bed(cluster5,"cluster5.bed")

#plot enrichmennt
data=read_excel("kmeans.xlsx")
names(data)
data=as.data.frame(data)
data=data[data$`Cell Line`=="BT245 K27M KO",]
data$RING1B=scale(data$RING1B)
data$CBX2=scale(data$CBX2)
data$H2AK119ub=scale(data$H2AK119ub)
data$H3K27me3=scale(data$H3K27me3)
data$H3K27me2=scale(data$H3K27me2)
data$H3K27me1=scale(data$H3K27me1)
data$H3K27ac=scale(data$H3K27ac)
data$H3K36me2=scale(data$H3K36me2)
data$H3K4me3=scale(data$H3K4me3)
data$H3K4me1=scale(data$H3K4me1)
data$SUZ12=scale(data$SUZ12)

library(dplyr)
library(reshape2)
library(ggplot2)
data=melt(data,id.vars=c("Cluster","Cell Line"))
ggplot(data = data, aes(variable, Cluster, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "dodgerblue1", high = "firebrick1", mid = "white", 
                       midpoint = 0, limit = c(-3.3,3.2), space = "Lab", 
                       name="Scaled\nEnrichment\nScore") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                   size = 12, hjust = 1))+
  ggtitle("BT245 K27M KO")+
  coord_fixed()+xlab("ChIP")+scale_y_continuous(breaks=c(1:6),labels = c(1:6))
ggsave("EnrichmentPlot_kmeansCluster_BT245_K27MKO_2.tiff", units="in", width=7, height=4, dpi=600, compression = 'lzw')

#RPKM levels of clusters
setwd("C:/Users/egj19/OneDrive - McGill University/Desktop/Analysis/Peakcalls")
cluster1.anno=read_excel("cluster1.anno.xlsx")
cluster2.anno=read_excel("cluster2.anno.xlsx")
cluster3.anno=read_excel("cluster3.anno.xlsx")
cluster4.anno=read_excel("cluster4.anno.xlsx")
cluster5.anno=read_excel("cluster5.anno.xlsx")
cluster6.anno=read_excel("cluster6.anno.xlsx")
setwd("C:/Users/egj19/OneDrive - McGill University/Desktop/Analysis/BEDs")

BT245_RPKM=read.table("BT245_RPKM.txt",header=T)
DIPG13_RPKM=read.table("DIPG13_RPKM.txt",header=T)
RPKM=left_join(BT245_RPKM,DIPG13_RPKM,by="Probe")
RPKM$K27M=(RPKM$BT245_Parental1+RPKM$BT245_Parental2+RPKM$BT245_Parental3+RPKM$BT245_Parental4+RPKM$BT245_Parental5+RPKM$BK.D13.DMSO.R2.bam+RPKM$BK.D13.DMSO.R3.bam+RPKM$BK.D13.DMSO.R1.bam)/8
RPKM$KO=(RPKM$C12.P2.bam+RPKM$C8.P2.bam+RPKM$BK.D13c5.DMSO.R1.bam+RPKM$BK.D13c5.DMSO.R2.bam+RPKM$C5.P2.bam+RPKM$C10.P2.bam+RPKM$BT245_Parental1+RPKM$BT245_Parental2+RPKM$BT245_Parental3+RPKM$BT245_Parental4+RPKM$BT245_Parental5)/11

cluster1.anno=left_join(cluster1.anno,RPKM,by="Probe")
cluster1.anno=cluster1.anno[complete.cases(cluster1.anno),]

cluster2.anno=left_join(cluster2.anno,RPKM,by="Probe")
cluster2.anno=cluster2.anno[complete.cases(cluster2.anno),]

cluster3.anno=left_join(cluster3.anno,RPKM,by="Probe")
cluster3.anno=cluster3.anno[complete.cases(cluster3.anno),]

cluster4.anno=left_join(cluster4.anno,RPKM,by="Probe")
cluster4.anno=cluster4.anno[complete.cases(cluster4.anno),]

cluster5.anno=left_join(cluster5.anno,RPKM,by="Probe")
cluster5.anno=cluster5.anno[complete.cases(cluster5.anno),]

cluster6.anno=left_join(cluster6.anno,RPKM,by="Probe")
cluster6.anno=cluster6.anno[complete.cases(cluster6.anno),]

padding  <- function(a,b) {
  zz <-rep(NA,length(X)-length(Y))
  YY  <- c(Y,zz)
  out <- data.frame(cbind(a,YY))
}

Y=cluster6.anno$KO
X=cluster6.anno$K27M
dd1 <- padding(X,Y)
library(reshape)
names(dd1)=c("K27M","Wild-type")
dd1=melt(dd1)
names(dd1)=c("Condition","RPKM")

library("ggsci")
library("ggplot2")
library("gridExtra")
ggplot(data = subset(dd1,!is.na(dd1$RPKM)), aes(x=Condition, y=RPKM,fill=Condition)) + 
  geom_boxplot()+  stat_boxplot()+
  labs(title="Cluster 6 RPKM",x="", y = "RPKMs")+ylim(0,50)+
  theme_classic()+scale_fill_manual(values=c("firebrick2","dodgerblue3"))
ggsave("Boxplot_Cluster6_RPKMs.tiff", units="in", width=7, height=4, dpi=600, compression = 'lzw')
