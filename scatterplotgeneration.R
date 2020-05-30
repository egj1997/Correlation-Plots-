library(readxl)
data <- read_excel("5kb.xlsx")
data[11]=log10(data[5]/data[10])
data[12]=log10(data[7]/data[6])
data[13]=log10(data[8]/data[9])

BiocManager::install("scatterplot3d")
library("scatterplot3d")
scatterplot3d(data[11], data[12], data[13], angle=20,col.axis="blue",
              col.grid="lightblue", main="Three-dimensional scatterplot", 
              pch=21, box=F, cex.symbols=2)

#create a scatterplot using matrix
scatterplot3d(logs[,1], y=logs[,2], z=logs[,3],
              main="3D Scatter Plot of 5Kb bins",
              xlab="RING1b", ylab="CBX2", zlab="K27me3", angle=60,axis=TRUE,grid=TRUE,box=TRUE,highlight.3d=TRUE)

logs=data[,c(11,12,13)]
logs[is.na(logs)] <- 0
try=logs[1:100,]
try=as.matrix(try)
logs=as.matrix(logs)
scatter3d(x = logs[,1], y = logs[,2], z = logs[,3], 
          surface=FALSE, grid = FALSE, ellipsoid = TRUE,
          surface.col = c("#999999", "#E69F00", "#56B4E9"))

#scatter using plotly
library(plotly)

colnames(logs)=c("RING1b","CBX2","K27me3")
logs=as.data.frame(logs)
p <- plot_ly(logs, x = logs[,1],y=logs[,2], z =logs[,3], colors = c('#BF382A', '#0C4B8E')) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'RING1b'),
                      yaxis = list(title = 'CBX2'),
                      zaxis = list(title = 'H3K27me3')))
p
chart_link = api_create(p, filename="scatter3d")
chart_link

plot(logs[,1],logs[,2])
ggplot(logs,aes(x=logs[,1],y=logs[,2]))+scatter.smooth()

#data
data=read.table("Parental_allmarks_chip_input.txt",header=T)
data$RING1b_BT245=log10(data$BT245.DMSO.2.Rx_XChIP_RING1B_condition3.sorted.dup.bam/data$BT245.DMSO.3.Rx_cells_ChIP1_Input_1.sorted.dup.bam)
data$RING1b_DIPG13=log10(data$DIPG13p14_RING1B.sorted.dup.bam/data$DIPG13.pool_Input.sorted.dup.bam)
data$H3K27me3_BT245=log10((data$BT245.C24_H3K27me3_1.sorted.dup.bam/data$BT245.C24_Input_1.sorted.dup.bam)*1)
data$H3K27me3_DIPG13=log10((data$DIPG13.C14_H3K27me3.sorted.dup.bam/data$DIPG13.C14_Input.sorted.dup.bam)*1)
data$H3K27ac_BT245=log10((data$BT245.nko.c1p5.Rx_cells_ChIP1_H3K27ac_1.sorted.dup.bam/data$BT245.nko.c1p5.Rx_cells_ChIP1_Input_1.sorted.dup.bam)*1)
data$H3K27ac_DIPG13=log10((data$DIPGXIIIp11_cells_ChIP1_H3K27ac_1.sorted.dup.bam/data$DIPG13p11_input.sorted.dup.bam)*1)
names(data)
logs2=data[,c(49,51,53)]
names(logs2)=c("RING1b","H3K27me3","H3K27Ac")
#data2

data=read.table("KO_allmarks_chip_input.txt",header=T)
data$RING1b_BT245=log10(data$BT245.ko.C2P8.Rx_XChIP_RING1B.sorted.dup.bam/data$BT245.ko.C2P8.Rx_XChIP_INPUT.sorted.dup.bam)
data$RING1b_DIPG13=log10(data$DIPGXIII.ko.c5p6.7.Rx_XChIP_RING1B.sorted.dup.bam/data$DIPGXIII.ko.c5p6.7.Rx_XChIP_Input.sorted.dup.bam)
data$H3K27me3_BT245=log10((data$BT245.ko.c4p6.Rx.C_cells_ChIP1_H3K27me3_1.sorted.dup.bam/data$BT245.ko.c4p6.Rx_cells_ChIP1_Input_1.sorted.dup.bam)*1)
data$H3K27me3_DIPG13=log10((data$DIPGXIII.ko.c10p5.Rx.C_cells_ChIP1_H3K27me3_1.sorted.dup.bam/data$DIPGXIII.ko.c10p5.Rx_cells_ChIP1_Input_1.bam.sorted.dup.bam)*1)
data$H3K27ac_BT245=log10((data$BT245.ko.c4p6.Rx_cells_ChIP1_H3K27ac_1.sorted.dup.bam/data$BT245.ko.c4p6.Rx_cells_ChIP1_Input_1.sorted.dup.bam)*1)
data$H3K27ac_DIPG13=log10((data$DIPGXIII.ko.c5p5.Rx_cells_ChIP1_H3K27ac_1.sorted.dup.bam/data$DIPGXIII.ko.c5p5.Rx_cells_ChIP1_Input_1.sorted.dup.bam)*1)
names(data)
logs2=data[,c(31,33,35)]
names(logs2)=c("RING1b","H3K27me3","H3K27Ac")
rownames(logs2) <- data[,1]
#sort by ring1b 
logs2.sort <- logs2[order(logs2$RING1b),]
#pick top 1%
library(dplyr)
logs2.1=top_frac(logs2.sort, 0.01, RING1b)
logs2.1=logs2.1[order(logs2.1$RING1b),]
#plot with matrix
library(scatterplot3d)
logs2.1=as.matrix(logs2.1)
scatterplot3d(logs2.1[,1], y=logs2.1[,2], z=logs2.1[,3],
              main="3D Scatter Plot of 5Kb bins",
              xlab="k27me3", ylab="ring1b", zlab="k27ac", angle=60,axis=TRUE,grid=TRUE,box=TRUE,highlight.3d=TRUE)

#interactive plot
p <- plot_ly(logs2.10[1:100,], x = logs2.10[1:100,1],y=logs2.10[1:100,2], z =logs2.10[1:100,3], colors = c('#BF382A', '#0C4B8E')) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'k27me3'),
                      yaxis = list(title = 'ring1b'),
                      zaxis = list(title = 'k27ac')))
p
chart_link = api_create(p, filename="scatter3d_k27ac_k27me3_ring1b")
chart_link


Sys.setenv("plotly_username"="egj1997")
Sys.setenv("plotly_api_key"="v3j00ThyiIiXH5oYV8x1")

#pca
library(factoextra)

logs2.pca <- prcomp(logs2, scale = TRUE)
fviz_eig(logs2.pca)


tiff("PCA_Arrow_BT245_parental.tiff",res=600,width = 10, height = 8, units = 'in')

fviz_pca_var(logs2.pca,geom = c("text","arrow"),
             col.ind = "cos2", # Color by the quality of representation
             col.var="blue",fill.var = "red",col.circle="red",
             label=c("var"),   # Avoid text overlapping
             title="PCA Biplot - BT245 Parental"
)
dev.off()
#PCA


tiff("PCA_scatter_DIPG13_KO.tiff",res=600,width = 10, height = 8, units = 'in')


fviz_pca_biplot(logs2.pca, label ="var",col.ind = "cos2", 
                gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"))+ggtitle('DIPGXIII K27M KO')
dev.off()

#load K27M data
data3=read.table("seqmonk5.txt",header=T)
data3[17]=log10(data3[5]/data3[6])
data3[18]=log10(data3[11]/data3[15])
data3[21]=log10(data3[13]/data3[14])
data3[22]=log10(data3[12]/data3[16])
data3[19]=log10(data3[7]/data3[8])
data3[20]=log10(data3[9]/data3[10])

data3[23]=(data3[17]+data3[18])/2
data3[24]=(data3[19]+data3[20])/2
data3[25]=(data3[22]+data3[21])/2

logs3=data3[,c(23,24,25)]
names(logs3)=c("H3K27me3","RING1b","H3K27Ac")
rownames(logs3) <- data3[,1]

#sort by ring1b 
logs3.sort <- logs3[order(logs3$RING1b),]
#pick top 1%
library(dplyr)
logs3.1=top_frac(logs3.sort, 0.01, RING1b)
logs3.1=logs3.1[order(logs3.1$RING1b),]

#pca for k27m top 1%
logs3.pca <- prcomp(logs3.1, scale = TRUE)
fviz_eig(logs3.pca)

fviz_pca_ind(logs3.pca, col.ind = "cos2", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             label=c("var"),
             legend.title = list(color = "Variability"),# Avoid text overlapping (slow if many points)
)

#k27m biplot top1%
fviz_pca_biplot(logs3.pca, label=c("var"),
                col.var = "#00AFBB", # Variables color
                col.ind = "#FC4E07"  # Individuals color
)



fviz_pca_biplot(logs3.pca, 
                # Fill individuals by groups
                geom.ind = "point",
                pointshape = 21,
                pointsize = 1.5,
                # Color variable by groups
                col.var = factor(c("RING1b","H3K27me3","H3K27ac")),
                legend.title = list(color = "Enrichment Trend"),
                repel = TRUE,
                palette=c("#00AFBB", "#E7B800", "#FC4E07")# Avoid label overplotting
)
+
  ggpubr::fill_palette("jco")+      # Indiviual fill color
  ggpubr::color_palette("npg")      # Variable colors


var_coord_func <- function(loadings, comp.sdev){
  loadings*comp.sdev
}


#load KO data
data4=read.table("seqmonk4.txt",header=T)
data4[16]=log10(data4[8]/data4[15])
data4[17]=log10(data4[10]/data4[12])
data4[20]=log10(data4[9]/data4[14])
data4[21]=log10(data4[5]/data4[11])
data4[18]=log10(data4[7]/data4[13])
data4[19]=log10(data4[6]/data4[12])

data4[22]=(data4[16]+data4[17])/2
data4[23]=(data4[18]+data4[19])/2
data4[24]=(data4[20]+data4[21])/2

logs4=data4[,c(22,23,24)]
names(logs4)=c("H3K27me3","RING1b","H3K27Ac")
rownames(logs4) <- data4[,1]

#sort by ring1b 
logs4.sort <- logs4[order(logs4$RING1b),]
#pick top 1%
library(dplyr)
logs4.1=top_frac(logs4.sort, 0.01, RING1b)
logs4.1=logs4.1[order(logs4.1$RING1b),]

logs7.pca=prcomp(logs3,scale=TRUE)
#pca for top 1%
library(factoextra)
logs4.pca <- prcomp(logs4.1, scale = TRUE)
fviz_eig(logs3.pca)

fviz_pca_ind(logs7.pca, col.ind = "cos2", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             label=c("var"),
             legend.title = list(color = "Variability"),# Avoid text overlapping (slow if many points)
)

#ko biplot top1%
fviz_pca_biplot(logs4.pca, label=c("var"),
                col.var = "#00AFBB", # Variables color
                col.ind = "#FC4E07"  # Individuals color
)



fviz_pca_biplot(logs4.pca, 
                # Fill individuals by groups
                geom.ind = "point",
                pointshape = 21,
                pointsize = 1.5,
                # Color variable by groups
                col.var = factor(c("RING1b","H3K27me3","H3K27ac")),
                legend.title = list(color = "Enrichment Trend"),
                repel = TRUE,
                palette=c("#00AFBB", "#E7B800", "#FC4E07")# Avoid label overplotting
)

#load K27M data
data3=read.table("seqmonk5.txt",header=T)
data3[17]=log10(data3[5]/data3[6])
data3[18]=log10(data3[11]/data3[15])
data3[21]=log10(data3[12]/data3[16])
data3[22]=log10(data3[13]/data3[14])
data3[19]=log10(data3[7]/data3[8])
data3[20]=log10(data3[9]/data3[10])

data3[23]=(data3[17]+data3[18])/2
data3[24]=(data3[19]+data3[20])/2
data3[25]=(data3[21]+data3[22])/2

logs3=data3[,c(23,24,25)]
names(logs3)=c("H3K27me3","RING1b","H3K27Ac")
rownames(logs3) <- data3[,1]

#sort by ring1b 
logs3.sort <- logs3[order(logs3$RING1b),]
#pick top 1%
library(dplyr)
logs3.1=top_frac(logs3.sort, 0.01, RING1b)
logs3.1=logs3.1[order(logs3.1$RING1b),]

#pca for k27m top 1%
logs3.pca <- prcomp(logs3.1, scale = TRUE)
fviz_eig(logs3.pca)

#visualize k27m data
fviz_pca_biplot(logs3.pca, label=c("var"),
                col.var = "#00AFBB", # Variables color
                col.ind = "#FC4E07"  # Individuals color
)
#find matched ranges of top 1%
library(GenomicRanges)
library(plyr)
k27m.ring1b.top1=logs3.sort[c(513395:518579),]
k27m.ring1b.top1[4]=rownames(k27m.ring1b.top1)
ko.ring1b.top1=logs4.sort[c(498906:504944),]
ko.ring1b.top1[4]=rownames(ko.ring1b.top1)


matched=inner_join(k27m.ring1b.top1,ko.ring1b.top1,by="V4")
names(matched)=c("H3K27me3_K27m","RING1b_K27M","H3K27ac_K27m","Bins","H3K27me3_KO","RING1b_KO","H3K27ac_KO")
logs5.pca <- prcomp(matched[c(1:3)], scale = TRUE)
logs6.pca<- prcomp(matched[c(5:7)], scale = TRUE)
fviz_eig(logs5.pca)
fviz_eig(logs6.pca)


fviz_pca_biplot(logs5.pca, label=c("var"),
                col.var = "#00AFBB", # Variables color
                col.ind = "#FC4E07"  # Individuals color
)

fviz_pca_biplot(logs6.pca, label=c("var"),
                col.var = "#00AFBB", # Variables color
                col.ind = "#FC4E07"  # Individuals color
)

#repeat biplot with suz12 in k27m:
data5=read.table("seqmonk6.txt",header=T)
names(data5)
data5[17]=log10(data5[5]/data5[6])
data5[18]=log10(data5[7]/data5[11])
data5[21]=log10(data5[8]/data5[12])
data5[22]=log10(data5[9]/data5[10])
data5[19]=log10(data5[13]/data5[16])
data5[20]=log10(data5[14]/data5[15])
data5[23]=(data5[17]+data5[18])/2
data5[24]=(data5[19]+data5[20])/2
data5[25]=(data5[21]+data5[22])/2
logs5=data5[,c(23,24,25)]
names(logs5)=c("H3K27me3","SUZ12","H3K27Ac")
rownames(logs5) <- data5[,1]
#sort by suz12 
logs5.sort <- logs5[order(logs5$SUZ12),]
#pick top 1%
library(dplyr)
logs5.1=top_frac(logs5.sort, 0.01, SUZ12)
logs5.1=logs5.1[order(logs5.1$SUZ12),]
#pca for k27m top 1%
logs5.pca <- prcomp(logs5.1, scale = TRUE)
fviz_eig(logs5.pca)
#Visualize Overlaps
fviz_pca_biplot(logs5.pca, label=c("var"),
                col.var = "#00AFBB", # Variables color
                col.ind = "#FC4E07"  # Individuals color
)

#repeat biplot with suz12 in k27m:
data6=read.table("seqmonk7.txt",header=T)
names(data6)
data6[17]=log10(data6[7]/data6[8])
data6[18]=log10(data6[14]/data6[15])
data6[21]=log10(data6[5]/data6[6])
data6[22]=log10(data6[12]/data6[11])
data6[19]=log10(data6[9]/data6[10])
data6[20]=log10(data6[13]/data6[11])
data6[23]=(data6[17]+data6[18])/2
data6[24]=(data6[19]+data6[20])/2
data6[25]=(data6[21]+data6[22])/2
logs6=data6[,c(23,24,25)]
names(logs6)=c("H3K27me3","SUZ12","H3K27Ac")
rownames(logs6) <- data6[,1]
#sort by suz12 
logs6.sort <- logs6[order(logs6$SUZ12),]
#pick top 1%
library(dplyr)
logs6.1=top_frac(logs6.sort, 0.01, SUZ12)
logs6.1=logs6.1[order(logs6.1$SUZ12),]
#pca for k27m top 1%
logs6.pca <- prcomp(logs6.1, scale = TRUE)
fviz_eig(logs6.pca)
#Visualize Overlaps
fviz_pca_biplot(logs6.pca, label=c("var"),
                col.var = "#00AFBB", # Variables color
                col.ind = "#FC4E07"  # Individuals color
)



#3d plot
install.packages("plot3D")
library("plot3D")
scatter3D(logs3.1$H3K27me3, logs3.1$RING1b, logs3.1$H3K27Ac, theta = 15, phi = 20)
# Compute Coordinates
#::::::::::::::::::::::::::::::::::::::::
loadings <- logs2.pca$rotation
sdev <- logs2.pca$sdev
var.coord <- t(apply(loadings, 1, var_coord_func, sdev)) 
head(var.coord)

#ggplot pca
library(ggfortify)
df <- logs2
autoplot(prcomp(df))
ggbiplot(df)+scale_color_gradient()

#ggplot correlation of ring1bpeaks and k27me3, k27ac
data3=read.table("dataset_ring1bpeaksk27mko_rx_normalizedtosize.txt",header=T)
logs3=data3[c(7,6,5,8,9)]
colnames(logs3)=c("H3K27me3","H3K27ac", "RING1b","SUZ12","H2AK119ub")
rownames(logs3)=data3[,1]
corlogs3 <- round(cor(logs3),4)
head(corlogs3)

library(reshape2)
melted_corlogs3 <- melt(corlogs3)
head(melted_corlogs3)

library(ggplot2)
ggheatmap <- ggplot(melted_corlogs3, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-0.5,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()

ggheatmap + 
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 6) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(2, 2.5),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))



#DIffBIND PACKAGE

library("DiffBind")
setwd("C:/Users/egj19/OneDrive - McGill University/Desktop/Analysis/peakcalls")
dat <- dba(sampleSheet = 'ring1b.csv')
dat <- dba.count(dat)

dat<- dba.contrast(dat, categories=DBA_CONDITION,minMembers=2)

dat_analysis=dba.analyze(dat)
dba.plotPCA(dat,label=DBA_CONDITION)

plot(dat,contrast=1)

corvals <- dba.plotHeatmap(dat)
dba.plotVenn(ring1b,dat$masks$Consensus)

#aggregate plot using soGGi

library("soGGi")

browseVignettes("soGGi")
library(ggplot2)
library(rtracklayer)
cgi=import.bed("CGI_centred10kb_hg19.bed")
a1 <- regionPlot("BigWig/Normalized/BT245-P47-H2AK119ub.RxNormalized.bw",cgi,format="bigwig")
a2<-regionPlot("BigWig/Normalized/BT245_RNF2KO_C10_H2AK119ub.RxNormalized.bw",cgi,format="bigwig")
plotRegion(a1)

c(a1,a2)
amerged=rbind(a1,a2)
amerged=list(a1,a2)
plotRegion(amerged,colourBy="Sample", groupBy="Sample", freeScale=TRUE)
assay()

#ggplot correlation matrix for all marks in k27m 
#import dataset with read values normalized to total read count and filtered for 3 reads
data=read.table("seqmonk2.txt",header=T)
data$K27M_RING1b=(data$RING1b_BT245+data$RING1b_DIPG13)/2
data$K27M_SUZ12=(data$SUZ12_BT245+data$SUZ12_DIPG13)/2
data$K27M_CBX2=(data$CBX2_BT245)
data$K27M_H3K27me3=(data$H3K27me3_BT245+data$H3K27me3_DIPG13)/2
data$K27M_H3K27ac=(data$H3K27ac_BT245+data$H3K27ac_DIPG13)/2
data$K27M_H3K36me2=(data$H3K36me2_BT245+data$H3K36me2_DIPG13)/2
data$K27M_H2AK119ub=(data$H2AK119ub_BT245+data$H2AK119ub_DIPG13)/2
data$K27M_H3K4me1=(data$H3K4me1_BT245+data$H3K4me1_DIPG13)/2
data$K27M_H3K4me3=(data$H3K4me3_BT245+data$H3K4me3_DIPG13)/2

colnames(data)[c(22,23,24,25,26,27,28,29,30)]=paste(c("RING1b","SUZ12","CBX2","H3K27me3","H3K27ac","H3K36me2","H2AK119ub","H3K4me1","H3K4me3"))
data1=data[order(-data$H3K27me3),]

#remove myc region
data1=data1[!data1$Chromosome=="8",]
data1=data1[1:1000,]
library(GenomicRanges)
findOverlaps(makeGRangesFromDataFrame(chr8.remove),makeGRangesListFromDataFrame(data))

#load package
library(ggcorrplot)
# Compute a correlation matrix
corr <- round(cor(data1[c(22,23,24,25,26,27,28,29,30)]), 3)
head(corr[, 1:6])
#compute p-value matrix
p.mat <- cor_pmat(data[c(22,23,24,25,26,27,28,29,30)])
head(p.mat[, 1:4])
#visualize matrix with ggcorrplot
ggcorrplot(corr, hc.order = TRUE, type = "lower",outline.col = "white",ggtheme = ggplot2::theme_gray,colors = c("#6D9EC1", "white", "#E46726"))
ggcorrplot(corr, method = "circle")           
#load ggplot2
library(ggplot2)
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
                       midpoint = 0.5, limit = c(0,1), space = "Lab", 
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
#ggplot correlation matrix for all marks in ko 
#import dataset with read values normalized to total read count and filtered for 3 reads
data=read.table("seqmonk3.txt",header=T)
data$RING1b=(data$RING1b_BT245_KO+data$RING1b_DIPG13_KO)/2
data$SUZ12=(data$SUZ12_BT245_KO+data$SUZ12_DIPG13_KO)/2
data$H3K27me3=(data$H3K27me3_BT245_KO+data$H3K27me3_DIPG13_KO)/2
data$H3K27ac=(data$H3K27ac_BT245_KO+data$H3K27ac_DIPG13_KO)/2
data$H3K36me2=(data$H3K36me2_BT245_KO+data$H3K36me2_DIPG13_KO)/2
data$H2AK119ub=(data$H2AK119ub_BT245_KO+data$H2AK119ub_DIPG13_KO)/2
# Compute a correlation matrix
corr <- round(cor(data[c(17,18,19,20,21,22)]), 3)
head(corr[, 1:6])
# Get upper triangle of the correlation matrix
get_upper_tri <- function(corr){
  corr[lower.tri(corr)]<- NA
  return(corr)
}
upper_tri <- get_upper_tri(corr)
# Melt the correlation matrix
library(reshape2)
melted_cormat <- melt(upper_tri, na.rm = TRUE)
#heatmap
ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0.5, limit = c(0,1), space = "Lab", 
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
