#djihad heatmap
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
data$RING1b_pcGBM2=data$pcGBM2p7_RING1B.sorted.dup.bam/data$pcGBM2p18_Input.sorted.dup.bam
data$RING1b_G477=data$G477p5_RING1B.sorted.dup.bam/data$G477p6_Input.sorted.dup.bam
data=data[complete.cases(data),]
data=data[is.finite(data$H2AK119ub_G477),]
data=data[is.finite(data$H2AK119ub_pcGBM2),]
data=data[is.finite(data$RING1b_G477),]
data=data[is.finite(data$RING1b_pcGBM2),]
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
names(merge)
BT245=data[c(2:5,78:88)]
RING1B=merge[c(1:5,80:81,107:108,162:163)]
RING1B=RING1B[complete.cases(RING1B),]
names(RING1B)=paste(c("Chromosome","Start","End","Width","Strand","BT245","DIPGXIII","pcGBM2","G477","BT245 KO","DIPGXIII KO"))
RING1B$sdev=apply(subset(RING1B, select = c(6:11)), 1, sd, na.rm=TRUE)
#RING1B=RING1B[,c(7:18)
names(merge)
H2AK119ub=merge[,c(1:5,89:90,97:102,158:159)]
names(H2AK119ub)=paste(c("Chromosome","Start","End","Width","Strand","BT245","DIPGXIII","G477","pcGBM2","BT245 KO 2","DIPGXIII KO 2","DIPGXIII KO 3","DIPGXIII 2","BT245 KO","DIPGXIII KO"))
H2AK119ub=H2AK119ub[complete.cases(H2AK119ub),]
H2AK119ub=H2AK119ub[is.finite(H2AK119ub$G477),]
H2AK119ub=H2AK119ub[is.finite(H2AK119ub$pcGBM2),]
#DIPG13=data[c(1:4,55,57,59,61,63,65,67,69,71)]
#BT245_KO=data[c(1:4,56,59,62,65,68,71,78)]
#combine overlapping regions for BT245
RING1B=RING1B[order(-RING1B$sdev),]
library(GenomicRanges)
library(rtracklayer)
names(BT245)
BT245_ring1b_ranges=(BT245[1:1000,c(1:14)])
write.table(BT245_ring1b_ranges,"test.bed",sep="\t")
BT245_ring1b_ranges_merged=import.bed('BT245_ring1b_ranges_merged.bed')

BT245_ring1b_ranges_merged=as.data.frame(BT245_ring1b_ranges_merged)
BT245_ring1b_ranges_merged$width=BT245_ring1b_ranges_merged$end-BT245_ring1b_ranges_merged$start
a1=as.data.frame(table(BT245_ring1b_ranges_merged$width))


#merge bin counts
#BT245_ring1b_ranges_merged=makeGRangesFromDataFrame(BT245_ring1b_ranges_merged)
#BT245_ring1b_ranges=makeGRangesFromDataFrame(BT245_ring1b_ranges)
#hits=as.list(findOverlaps(BT245_ring1b_ranges_merged,BT245_ring1b_ranges))
#BT245_ring1b_ranges_merged=as.data.frame(BT245_ring1b_ranges_merged)
#BT245_ring1b_ranges_merged$H3K27me3_BT245 <- sapply(extractList(BT245$H3K27me3_BT245, hits),mean)
#BT245_ring1b_ranges_merged$H3K27ac_BT245 <- sapply(extractList(BT245$H3K27ac_BT245, hits),mean)
#BT245_ring1b_ranges_merged$SUZ12_BT245 <- sapply(extractList(BT245$SUZ12_BT245, hits),mean)
#BT245_ring1b_ranges_merged$H3K27me2_BT245 <- sapply(extractList(BT245$H3K27me2_BT245, hits),mean)
#BT245_ring1b_ranges_merged$H3K4me1_BT245 <- sapply(extractList(BT245$H3K4me1_BT245, hits),mean)
#BT245_ring1b_ranges_merged$H3K4me3_BT245 <- sapply(extractList(BT245$H3K4me3_BT245, hits),mean)
#BT245_ring1b_ranges_merged$H3K36me2_BT245 <- sapply(extractList(BT245$H3K36me2_BT245, hits),mean)
#BT245_ring1b_ranges_merged$H2AK119ub_BT245 <- sapply(extractList(BT245$H2AK119ub_BT245, hits),mean)
#BT245_ring1b_ranges_merged$RING1b_BT245 <- sapply(extractList(BT245$RING1b_BT245, hits),mean)
#BT245_ring1b_ranges_merged$Probe=paste("chr",BT245_ring1b_ranges_merged$seqnames,BT245_ring1b_ranges_merged$start,"-",BT245_ring1b_ranges_merged$end)

#write.table(BT245_ring1b_ranges_merged,"BT245_ring1b_ranges_merged.txt",sep="\t")
BT245_ring1b_ranges_merged=read.table('BT245_ring1b_ranges_merged.txt',header=T)
BT245_ring1b_ranges_merged$Probe=paste(BT245_ring1b_ranges_merged$chromosome,BT245_ring1b_ranges_merged$start,"-",BT245_ring1b_ranges_merged$end)

test=makeGRangesFromDataFrame(BT245_ring1b_ranges_merged)
export.bed(test,"test.bed")
#arrange BT245 KO data set for BT245 top 15K peaks
BT245_KO
BT245_ring1b_ranges_merged=makeGRangesFromDataFrame(BT245_ring1b_ranges_merged)
BT245_ring1b_ranges=makeGRangesFromDataFrame(BT245_KO)
hits=as.list(findOverlaps(BT245_ring1b_ranges_merged,BT245_ring1b_ranges))
BT245ko_ring1b_ranges_merged=as.data.frame(BT245_ring1b_ranges_merged)
BT245ko_ring1b_ranges_merged$H3K27me3_BT245_KO <- sapply(extractList(BT245_KO$H3K27me3_BT245_KO, hits),mean)
BT245ko_ring1b_ranges_merged$H3K27ac_BT245_KO <- sapply(extractList(BT245_KO$H3K27ac_BT245_KO, hits),mean)
BT245ko_ring1b_ranges_merged$SUZ12_BT245_KO<- sapply(extractList(BT245_KO$SUZ12_BT245_KO, hits),mean)
BT245ko_ring1b_ranges_merged$H3K27me2_BT245_KO <- sapply(extractList(BT245_KO$H3K27me2_BT245_KO, hits),mean)
BT245ko_ring1b_ranges_merged$H3K36me2_BT245_KO <- sapply(extractList(BT245_KO$H3K36me2_BT245_KO, hits),mean)
BT245ko_ring1b_ranges_merged$H2AK119ub_BT245_KO <- sapply(extractList(BT245_KO$H2AK119ub_BT245_KO, hits),mean)
BT245ko_ring1b_ranges_merged$RING1b_BT245_KO <- sapply(extractList(BT245_KO$RING1b_BT245_KO, hits),mean)
BT245ko_ring1b_ranges_merged$Probe=paste("chr",BT245ko_ring1b_ranges_merged$seqnames,BT245ko_ring1b_ranges_merged$start,"-",BT245ko_ring1b_ranges_merged$end)
#add genes
a=makeGRangesFromDataFrame(BT245_ring1b_ranges_merged)
b=makeGRangesFromDataFrame(BT245_top15kbins_annotation)
BT245_top15kbins_annotation=read.table("BT245_top15kbins_annotation.txt",header=T)
hits=as.list(findOverlaps(a,b))
BT245_ring1b_ranges_merged$gene=sapply(extractList(BT245_top15kbins_annotation$Gene_Name,hits),paste)
write.table(BT245_ring1b_ranges_merged,"BT245_ring1b_ranges_merged.txt",sep="\t")
#add RPKMs
merge1=BT245_RNA[c(1,12)]
names(merge1)=c("gene","RPKMs")
library(dplyr)
BT245_ring1b_ranges_merged=left_join(BT245_ring1b_ranges_merged,merge1)
BT245_ring1b_ranges_merged=BT245_ring1b_ranges_merged[complete.cases(BT245_ring1b_ranges_merged),]
#plot heatmap (djihad's script)

### Djihad script for Heatmap in R ###
suppressPackageStartupMessages("Need help? Try Stackoverflow: https://stackoverflow.com/tags/ggplot2.")
suppressPackageStartupMessages("Welcome! Related Books: `Practical Guide To Cluster Analysis in R` at https://goo.gl/13EFCZ")

library(RColorBrewer)
library(ggplot2)
library(factoextra)

#my.colours <- c("#313695","#74ADD1","#E0F3F8","#FEE090","#F46D43", "#D73027","#A50026")
my.colours <- rev(brewer.pal(11,"RdYlBu"))
library(dplyr)
data <- RING1B
names(data)
data$sdev=apply(subset(data, select = c(6:11)), 1, sd, na.rm=TRUE)
data$rank=rank(data$sdev)
summary(data$sdev)
data=data[order(-data$sdev),]
dim(data)

names(data)
#data=data[1:10000,]
#random=data[sample(nrow(data),1000), ]
#data=random
cl <- data[1:1000,c(6:11)]
rownames(cl)=cl$Probe
names(cl)


dim(cl)

clt <- t(na.omit(cl[,1:6]))

cltn <- data.matrix(clt)

mode(cltn)<-'numeric'

dim(cltn)

#add clusters
for ( i in seq(1000, 1000, 1000)){
     nbs = i 
     nbr = nbs  
}

distancem <- dist(1-cor(cltn))
hclust_completem <- hclust(distancem, method = "average")
groupe = cutree(hclust_completem , k = 6) #k correspond aux nombres de groupes que nous souhaitons
groupe = groupe[order(groupe)]
groupe1 = names(groupe[groupe[1:length(groupe)] == 1])
groupe2 = names(groupe[groupe[1:length(groupe)] == 2])
groupe3 = names(groupe[groupe[1:length(groupe)] == 3])
groupe4 = names(groupe[groupe[1:length(groupe)] == 4])
groupe5 = names(groupe[groupe[1:length(groupe)] == 5])
groupe6 = names(groupe[groupe[1:length(groupe)] == 6])
filename = paste ("groupe_", nbs,"groupe1.txt")
groupe_file <- file(filename, open = "w")
cat("Groupe 1 :\n", groupe1, "\n", file = groupe_file, sep = "\t")
cat("Groupe 2 :\n", groupe2, "\n", file = groupe_file, sep = "\t")
cat("Groupe 3 :\n", groupe3, "\n", file = groupe_file, sep = "\t")
cat("Groupe 4 :\n", groupe4, "\n", file = groupe_file, sep = "\t")
cat("Groupe 5 :\n", groupe5, "\n", file = groupe_file, sep = "\t")
cat("Groupe 6 :\n", groupe6, "\n", file = groupe_file, sep = "\t")
#cat("Groupe 7 :\n", groupe7, "\n", file = groupe_file, sep = "\t")
close(groupe_file)


names(data)
#range is 1 and all your histone marks
blast <- names(data[1,c(6:11)])

blastv <- as.character(blast)

blastc <- as.character(data[1:1000,1])

names(data)
gr <- data[1,c(6:11)]
grn <- as.numeric(gr)
f.gr <- factor(grn)
c.gr <- rainbow(nlevels(f.gr), start = 0.1, end = 0.9)
gr.color <- rep(length(f.gr), 0)
for (i in 1:length(f.gr))
  gr.color[i] <- c.gr[f.gr[i]==levels(f.gr)]

grc <- data[1:500,1]
grc <- as.numeric (grc)
f.grc <- factor (grc)
c.grc <- c("#313695","#E0F3F8")
c.grc <- rainbow(nlevels(f.grc), start = 0.1, end = 0.9)
gr.colorc <- rep(length(f.grc), 0)
for (i in 1:length(f.grc))
gr.colorc[i] <- c.grc[f.grc[i]==levels(f.grc)]

jpeg("Heatmap_top1Kpeaks_clusterALLRING1b_2.jpeg", width = 12, height =9, units = "cm", res = 1500)

par(mar=c(3, 3, 1, 1))

hcf <- function(cltn)
{
  hclust(cltn, method = "ward.D2")       
}

df <- function(cltn)
{
  dist(cltn, method = "euclidean")
}

heatmap(cltn, distfun = df, hclustfun = hcf, labCol = blastc, labRow = blastv, col = my.colours, RowSideColors = gr.color, cexRow=0.35, cexCol=0.3, scale = "row",revC = T)

title(main="Top 1000 variant RING1b bins", cex.main=0.6, adj = 0)

dev.off()
#extract clusters
hr<-hclust(as.dist(1-cor(cltn,method="pearson")), method="ward.D2")
x=cutree(hr,k=6)
#sort
ggplot(data = data, mapping = aes(x = rank, y = RING1b_DIPG13)) +
  geom_line(color="firebrick1")+theme_classic()+scale_y_continuous(name = "RING1b Enrichment Fold")+
  scale_x_continuous(name="Bin Rank")+ggtitle("DIPGXIII Parental")
ggsave("DIPG13_RING1b_BinsFold.tiff", units="in", width=5, height=4, dpi=600, compression = 'lzw')


#width
library(tidyverse)
library(hrbrthemes)
hist(BT245_ring1b_ranges_merged$width, 
     main="RING1b bin size", 
     xlab="Width", 
     border="black", 
     col="firebrick2",
     xlim=c(0,36000),
     las=1, 
     breaks=10,ylim=c(0,12000))
#export bed based on width
test=BT245_ring1b_ranges_merged[BT245_ring1b_ranges_merged$width>5000,]
test=makeGRangesFromDataFrame(test)
export.bed(test,"test.bed")
#boxplot
if BT245_ring1b_ranges_merged$width>5000 { 
  BT245_ring1b_ranges_merged$type=paste("Domain") 
  } else {
  BT245_ring1b_ranges_merged$type="Peaks"
  }
BT245_ring1b_ranges_merged[BT245_ring1b_ranges_merged$width>5000,]$type="Domain"
library(ggplot2)
ggplot(BT245_ring1b_ranges_merged, aes(x=type, y=width)) + 
  geom_boxplot()+scale_fill_manual(values=c("#999999", "#56B4E9"))+theme_classic()

#prepare ChIP_RNAseq dataset
setwd("C:/Users/egj19/OneDrive - McGill University/Desktop/Analysis/Peakcalls")
BT245_RNA=read.table("seqmonk17.txt",header=T)
BT245_RNA$BT245=(BT245_RNA$BT245.NKO__BT245.NKO.P17.sorted.bam+BT245_RNA$BT245.NKO__BT245.NKO.P18.sorted.bam+BT245_RNA$BT245.NKO__BT245.NKO.P30.sorted.bam)/3 
BT245_RNA$BT245_KO=(BT245_RNA[5]+BT245_RNA[6]+BT245_RNA[7]+BT245_RNA[8])/4
clusters=read.table("bt245_clusters.txt",header=T)
library(GenomicRanges)
cluster1.ranges=makeGRangesFromDataFrame(clusters[clusters$Cluster==1,])
BT245_RNA.ranges=makeGRangesFromDataFrame(BT245_RNA)
cluster1.rna.range=subsetByOverlaps(BT245_RNA.ranges,cluster1.ranges)
cluster1.rna.range=as.data.frame(cluster1.rna.range)
cluster1.rna.range$link=paste(cluster1.rna.range$seqnames,cluster1.rna.range$start,cluster1.rna.range$end)
BT245_RNA$link=paste(BT245_RNA$Chromosome,BT245_RNA$Start,BT245_RNA$End)
names(BT245_RNA)
cluster1.rna.range=left_join(cluster1.rna.range,BT245_RNA[c(14,12)])
mean(cluster1.rna.range$BT245)

cluster2.ranges=makeGRangesFromDataFrame(clusters[clusters$Cluster==2,])
BT245_RNA.ranges=makeGRangesFromDataFrame(BT245_RNA)
cluster2.rna.range=subsetByOverlaps(BT245_RNA.ranges,cluster2.ranges)
cluster2.rna.range=as.data.frame(cluster2.rna.range)
cluster2.rna.range$link=paste(cluster2.rna.range$seqnames,cluster2.rna.range$start,cluster2.rna.range$end)
BT245_RNA$link=paste(BT245_RNA$Chromosome,BT245_RNA$Start,BT245_RNA$End)
names(BT245_RNA)
cluster2.rna.range=left_join(cluster2.rna.range,BT245_RNA[c(14,12)])
mean(cluster2.rna.range$BT245)

cluster3.ranges=makeGRangesFromDataFrame(clusters[clusters$Cluster==3,])
BT245_RNA.ranges=makeGRangesFromDataFrame(BT245_RNA)
cluster3.rna.range=subsetByOverlaps(BT245_RNA.ranges,cluster3.ranges)
cluster3.rna.range=as.data.frame(cluster3.rna.range)
cluster3.rna.range$link=paste(cluster3.rna.range$seqnames,cluster3.rna.range$start,cluster3.rna.range$end)
BT245_RNA$link=paste(BT245_RNA$Chromosome,BT245_RNA$Start,BT245_RNA$End)
names(BT245_RNA)
cluster3.rna.range=left_join(cluster3.rna.range,BT245_RNA[c(14,12)])
mean(cluster3.rna.range$BT245)

##import loops
BT245_loops=read.table("BT245_loops.txt",header=F)
BT245_c2_loops=read.table("BT245_c2_loops.txt",header=F)
BT245_loops_a=makeGRangesFromDataFrame(BT245_loops,seqnames.field = "V1",start.field = "V2",end.field = "V3")

BT245_loops_b=makeGRangesFromDataFrame(BT245_loops,seqnames.field = "V4",start.field = "V5",end.field = "V6")

BT245_gainedloops_a=subsetByOverlaps(BT245_gainedpeaks,BT245_loops_a)
BT245_gainedloops_b=subsetByOverlaps(BT245_gainedpeaks,BT245_loops_b)

#PCA of H2AK119ub data
# Load data
names(H2AK119ub)
data=H2AK119ub[1:100,c(7:12)]
rownames(data)=H2AK119ub$Probe
# Compute distances and hierarchical clustering
library("factoextra")
res.pca <- prcomp(data, scale = TRUE)
fviz_eig(res.pca)

fviz_pca_ind(res.pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = FALSE     # Avoid text overlapping
)
