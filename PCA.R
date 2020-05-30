#create PCA matrix
setwd("C:/Users/egj19/OneDrive - McGill University/Desktop/Analysis/peakcalls")
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
data$H2AK119ub_DIPG13=data$DIPG13.P29_H2AK119UB.sorted.dup.bam.x/data$DIPG13.P29.RX_XCHIP_INPUT.sorted.dup.bam
data$SUZ12_BT245=data$BT245p11_SUZ12_1.sorted.dup.bam/data$DIPG13p11_input.sorted.dup.bam
data$SUZ12_DIPG13=data$DIPG13p14_SUZ12.sorted.dup.bam/data$DIPG13.pool_Input.sorted.dup.bam
data$H3K27me2_BT245=data$BT245.nko.c1p5.Rx_cells_ChIP1_H3K27me2.sorted.dup.bam/data$BT245.nko.c1p5_ChIP1_Input_1.sorted.dup.bam
data$H3K27me2_DIPG13=data$DIPGXIII.nko.c12p3.Rx_cells_ChIP1_H3K27me2_1.sorted.dup.bam/data$DIPGXIII.nko.c12p3.Rx_cells_ChIP1_Input_1.sorted.dup.bam
data$H3K27me1_BT245=data$BT245.nko.c1p5.Rx_cells_ChIP1_H3K27me1_1.sorted.dup.bam/data$BT245.nko.c1p5_ChIP1_Input_1.sorted.dup.bam
data$H3K27me1_DIPG13=data$DIPGXIII.nko.c12p3.Rx_cells_ChIP1_H3K27me2_1.sorted.dup.bam/data$DIPGXIII.nko.c12p3.Rx_cells_ChIP1_Input_1.sorted.dup.bam
data$H2AK119ub_G477=data$G477.P25_H2K119UB.sorted.dup.bam/data$G477.pool_Input.sorted.dup.bam
data$H2AK119ub_pcGBM2=data$pcGBM2.p27_H2AK119ub.sorted.dup.bam/data$pcGBM2.p27_Input.sorted.dup.bam
data$H2AK119ub_BT245_KO_2=data$BT245.C2P7.K27MKO_H2AK119UB.sorted.dup.bam.y/data$BT245.K27MKO.C2P7_INPUT.sorted.dup.bam.y
data$H2AK119ub_DIPG13_KO_2=data$DIPG13.KO.C5P5_H2K119UB.sorted.dup.bam/data$DIPG13.ko.c5p5_Input.sorted.dup.bam
data$H2AK119ub_DIPG13_KO_3=data$DIPG13.C2P3.K27MKO_H2AK119UB.sorted.dup.bam.y/data$DIPG13.K27MKO.C2P3_INPUT.sorted.dup.bam.y
data$H2AK119ub_DIPG13_2=data$DIPG13_P66_Rx_XChIP_H2AK119ub_80M.sorted.dup.bam/data$DIPG13_P66_Rx_XChIP_INPUT.sorted.dup.bam.y
data=data[complete.cases(data),]
data=data[is.finite(data$H2AK119ub_G477),]
data=data[is.finite(data$H2AK119ub_pcGBM2),]
data$probe=paste(data$seqnames,data$start,data$end)
setwd("C:/Users/egj19/OneDrive - McGill University/Desktop/Analysis/peakcalls")
data2=read.table("KO_allmarks_chip_input.txt",header=T)
data2$RING1b_BT245_KO=data2$BT245.ko.C2P8.Rx_XChIP_RING1B.sorted.dup.bam/data2$BT245.ko.C2P8.Rx_XChIP_INPUT.sorted.dup.bam
data2$RING1b_DIPG13_KO=data2$DIPGXIII.ko.c5p6.7.Rx_XChIP_RING1B.sorted.dup.bam/data2$DIPGXIII.ko.c5p6.7.Rx_XChIP_Input.sorted.dup.bam
data2$H2AK119ub_BT245_KO=data2$BT245.ko.C2P8_H2AK119UB.sorted.dup.bam/data2$BT245.ko.C2P8.Rx_XChIP_INPUT.sorted.dup.bam
data2$H2AK119ub_DIPG13_KO=data2$DIPG13.ko.C5P33_H2AK119UB.sorted.dup.bam/data2$DIPG13.ko.C5P33.Rx_XChIP_INPUT.sorted.dup.bam
data2$H3K27me3_BT245_KO=data2$BT245.ko.c4p6.Rx.C_cells_ChIP1_H3K27me3_1.sorted.dup.bam/data2$BT245.ko.c4p6.Rx_cells_ChIP1_Input_1.sorted.dup.bam
data2$H3K27me3_DIPG13_KO=data2$DIPGXIII.ko.c10p5.Rx.C_cells_ChIP1_H3K27me3_1.sorted.dup.bam/data2$DIPGXIII.ko.c10p5.Rx_cells_ChIP1_Input_1.bam.sorted.dup.bam
data2$H3K27ac_BT245_KO=data2$BT245.ko.c4p6.Rx_cells_ChIP1_H3K27ac_1.sorted.dup.bam/data2$BT245.ko.c4p6.Rx_cells_ChIP1_Input_1.sorted.dup.bam
data2$H3K27ac_DIPG13_KO=data2$DIPGXIII.ko.c10p5.Rx_cells_ChIP1_H3K27ac_1.sorted.dup.bam/data2$DIPGXIII.ko.c10p5.Rx_cells_ChIP1_Input_1.bam.sorted.dup.bam
data2$H3K36me2_BT245_KO=data2$BT245.ko.c4p6.Rx_cells_ChIP1_H3K36me2_1.sorted.dup.bam/data2$BT245.ko.c4p6.Rx_cells_ChIP1_Input_1.sorted.dup.bam
data2$H3K36me2_DIPG13_KO=data2$DIPGXIII.ko.c10p5.Rx_cells_ChIP1_H3K36me2_1.sorted.dup.bam/data2$DIPGXIII.ko.c10p5.Rx_cells_ChIP1_Input_1.bam.sorted.dup.bam
data2$SUZ12_BT245_KO=data2$BT245.ko.c4p9.Rx_cells_ChIP1_SUZ12_1.sorted.dup.bam/data2$BT245.ko.c4p9_Input_1.sorted.dup.bam.x
data2$SUZ12_DIPG13_KO=data2$DIPGXIII.ko.c5p5.Rx_cells_ChIP1_SUZ12_1.sorted.dup.bam/data2$DIPG13.ko.c5p5_Input.sorted.dup.bam
data2$H3K27me2_BT245_KO=data2$BT245.ko.c4p9.Rx_cells_ChIP1_H3K27me2_1.sorted.dup.bam/data2$BT245.ko.c4p9_Input_1.sorted.dup.bam.x
data2$H3K27me2_DIPG13_KO=data2$DIPGXIII.ko.c5p5.Rx_cells_ChIP1_H3K27me2_1.sorted.dup.bam/data2$DIPGXIII.ko.c5p5.Rx_cells_ChIP1_Input_1.sorted.dup.bam
data2$H3K27me1_BT245_KO=data2$BT245.ko.c4p9.Rx_cells_ChIP1_H3K27me1_1.sorted.dup.bam/data2$BT245.ko.c4p9_Input_1.sorted.dup.bam.y

merge=left_join(data,data2,by='probe')
merge=merge[complete.cases(merge),]
names(merge)
matrix=merge[,c(89:90,97:102,158:159)]
#matrix.2=merge[,c(1:3,77:78,158:159)]
#matrix.2=makeGRangesFromDataFrame(matrix.2,seqnames.field = "seqnames.x",start.field = "start.x",end.field = "end.x")
#matrix.2=subsetByOverlaps(matrix.2,cgi_10kb)
#matrix.2=as.data.frame(matrix.2)
#matrix.2$probe=paste(matrix.2$seqnames,matrix.2$start,matrix.2$end)
#matrix.2=left_join(matrix.2,merge[,c(62,77:78,158:159)],by="probe")
#matrix.2=matrix.2[complete.cases(matrix.2),]
#matrix.2=matrix.2[,c(7:10)]
rownames(matrix)=merge$probe
#rownames(matrix.2)=merge$probe
#names(matrix.2)=paste(c("BT245","DIPGXIII","BT245 K27M KO","DIPGXIII K27M KO"))
names(matrix)=paste(c("BT245 Parental","DIPGXIII Parental","G477","pcGBM2","BT245 K27M KO 2","DIPGXIII K27M KO 2","DIPGXIII K27M KO 3","DIPGXIII Parental 2","BT245 K27M KO","DIPGXIII K27M KO"))
matrix$sdev=apply(subset(matrix, select = c(1:10)), 1, sd, na.rm=TRUE)
#matrix.2$sdev=apply(subset(matrix.2, select = c(1:4)), 1, sd, na.rm=TRUE)
#matrix.2=matrix.2[order(-matrix.2$sdev),]
#matrix.2=matrix.2[1:10000,c(1:4)]
matrix=matrix[order(-matrix$sdev),]
matrix_1k=matrix[,c(1:4,5:8,9:10)]
#matrix_1k.melt=melt(matrix_1k)
library(factoextra)
matrix_1k=matrix_1k[,-c(5:8)]
res.pca <- prcomp(matrix_1k, scale = TRUE)

tiff("PCA_BiPlot_H2AK119ub_all_3.tiff",width=10,height=7,compression="lzw",res = 600,units="in")
fviz_pca_biplot(res.pca, label = "var",
                col.var = "black", # Variables color
                col.ind = "cos2",  # Individuals color
                gradient.cols = c("white", "#2E9FDF", "#FC4E07"),
                title=("PCA Biplot of the top 1000 variant bins"),
                
)

dev.off()

