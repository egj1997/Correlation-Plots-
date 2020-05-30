BiocManager::install("ComplexHeatmap")
library("ComplexHeatmap")
#prepare K27M dataset
setwd("C:/Users/egj19/OneDrive - McGill University/Desktop/Analysis/peakcalls")
#correlation plots using RPKMs normalized to inputs. Load Data that is clean
data=read.table("Parental_allmarks_chip_input.txt",header=T)

data$RING1b_BT245=data$BT245.DMSO.2.Rx_XChIP_RING1B_condition3.sorted.dup.bam/data$BT245.DMSO.3.Rx_cells_ChIP1_Input_1.sorted.dup.bam
data$RING1b_DIPG13=data$DIPG13p14_RING1B.sorted.dup.bam/data$DIPG13.pool_Input.sorted.dup.bam
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
#correlation matrix of X Parental
data$RING1b=data$RING1b_BT245
data$SUZ12=data$SUZ12_BT245
data$CBX2=data$CBX2_BT245
data$H3K27me3=data$H3K27me3_BT245
data$H2AK119ub=data$H2AK119ub_BT245
data$H3K27me2=data$H3K27me2_BT245
data$H3K27ac=data$H3K27ac_BT245
data$H3K36me2=data$H3K36me2_BT245
data$H3K4me3=data$H3K4me3_BT245
data$H3K4me1=data$H3K4me1_BT245

data1=as.matrix(data[,48:66])
data1=data1[order(-data$RING1b),]
data1=data1[1:5000,]
#prepare KO dataset
data2=read.table("seqmonk3.txt",header=T)
data2$RING1b=scale((data2$RING1b_BT245_KO+data2$RING1b_DIPG13_KO)/2)
data2$SUZ12=scale((data2$SUZ12_BT245_KO+data2$SUZ12_DIPG13_KO)/2)
data2$H3K27me3=scale((data2$H3K27me3_BT245_KO+data2$H3K27me3_DIPG13_KO)/2)
data2$H3K27ac=scale((data2$H3K27ac_BT245_KO+data2$H3K27ac_DIPG13_KO)/2)
data2$H3K36me2=scale((data2$H3K36me2_BT245_KO+data2$H3K36me2_DIPG13_KO)/2)
data2$H2AK119ub=scale((data2$H2AK119ub_BT245_KO+data2$H2AK119ub_DIPG13_KO)/2)
names(data2)
data3=as.matrix(data2[17:22])
data3=data3[order(-data3[,1]),]

#complexHeatmap
library(ComplexHeatmap)
Heatmap(data1, 
        name = "K27M", #title of legend
        column_title = "Histone Mark/ Transcription Factor", row_title = "10Kb Bins",
        row_names_gp = gpar(fontsize = 7) # Text size for row names
)

Heatmap(data3[1:1000,], 
        name = "KO", #title of legend
        column_title = "Histone Mark/ Transcription Factor", row_title = "10Kb Bins",
        row_names_gp = gpar(fontsize = 7) # Text size for row names
)

BiocManager::install("cluster")
library("cluster")
set.seed(2)
Heatmap(data1[1:5000,], name = "K27M_H3K27ac", k = 2)

Heatmap(data1[1:50,], name ="K27M", 
        split = data.frame(H3K27me3= data1$K27M_H3K27me3, H3K27ac=data1[5,]),
        row_names_gp = gpar(fontsize = 7))

#d3heatmap
library("d3heatmap")
d3heatmap(scale(data1[1:1000,]), colors = "RdYlBu",
          k_row = 4, # Number of groups in rows
          k_col = 2 # Number of groups in columns
)



#densitymap
densityHeatmap(data1[1:50,])
#heatmap2
install.packages("gplots")
library("gplots")
heatmap.2(df, scale = "none", col = bluered(100), 
          trace = "none", density.info = "none")