#library loads
library(GenomicRanges)
library(rtracklayer)
library(ggthemes)
library(hrbrthemes)
library("ggsci")
library("ggplot2")
library("gridExtra")
#load data
setwd("C:/Users/egj19/OneDrive - McGill University/Desktop/Analysis/peakcalls")
cgi_10kb=import.bed("CGI_centred10kb_hg19.bed")
cgi=import.bed("cpgislands_hg19.bed")
BT245_CBX2=import.bed("BT245_CBX2_peaks_seqmonk.bed")
DIPG13_CBX2=import.bed("DIPG13_CBX2_peaks_seqmonk.bed")
BT245_KO_CBX2=import.bed("BT245_KO_CBX2_peaks_seqmonk.bed")
DIPG13_KO_CBX2=import.bed("DIPG13_KO_CBX2_peaks_seqmonk.bed")
setwd("C:/Users/egj19/OneDrive - McGill University/Desktop/Analysis/BEDs")
BT245_K27ac=import.bed("BT245_K27ac.bed")
BT245_K27MKO_K27ac=import.bed("BT245_K27MKO_K27ac.bed")
DIPG13_K27ac=import.bed("DIPG13_K27ac.bed")
DIPG13_K27MKO_K27ac=import.bed("DIPG13_K27MKO_K27ac.bed")
BT245=import.bed("BT245_RING1b_peaks_optimized.bed")
DIPG13=import.bed("DIPG13_RING1b_peaks_optimized.bed")
BT245.ko=import.bed("BT245_C2P8_RING1b_peaks_optimized.bed")
DIPG13.ko=import.bed("DIPG13_C5P6_RING1b_peaks_optimized.bed")
BT245_K27me3=import.bed("BT245_H3K27me3_10kb.bed")
BT245_K27MKO_K27me3=import.bed("BT245_K27MKO_H3K27me3_10kb.bed")
DIPG13_K27me3=import.bed("DIPG13_H3K27me3_10kb.bed")
DIPG13_K27MKO_K27me3=import.bed("DIPG13_K27MKO_H3K27me3_10kb.bed")
BT245_SUZ12=read.table("BT245_SUZ12.narrowPeak.bed")
BT245_SUZ12=makeGRangesFromDataFrame(BT245_SUZ12,seqnames.field= "V1",start.field = "V2",end.field = "V3")
DIPG13_SUZ12=read.table("DIPGXIII_SUZ12.narrowPeak.bed")
DIPG13_SUZ12=makeGRangesFromDataFrame(DIPG13_SUZ12,seqnames.field= "V1",start.field = "V2",end.field = "V3")
setwd("C:/Users/egj19/OneDrive - McGill University/Desktop/Analysis/peakcalls")
HSJ51=import.bed("HSJ51_RING1b_peaks_seqmonk.bed")
HSJ51_K27me3=import.bed("HSJ51_H3K27me3_10Kb.bed")
HSJ51_K27ac=import.bed("HSJ51_H3K27ac_peaks.bed")

pcGBM2=import.bed("pcGBM2_RING1b_peaks_seqmonk.bed")
pcGBM2_K27me3=import.bed("pcGBM2_H3K27me3_10kb.bed")
pcGBM2_K27ac=import.bed("pcGBM2_H3K27ac_peaks.bed")

#ashot's paper for ring1b
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
data$CBX2_DIPG13=data$DIPG13_p73_XxChIP_CBX2_c2.sorted.dup.bam/data$DIPG13.P73_XxChIP_Input.sorted.dup.bam
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
data$H3K27me3_BT245_2=data$BT245.nko.c1p5_H3K27me3_1.sorted.dup.bam/data$BT245.nko.c1p5.Rx_cells_ChIP1_Input_1.sorted.dup.bam
data$H3K27me3_DIPG13_2=data$DIPG13.C12_H3K27me3.sorted.dup.bam/data$DIPG13.C12_Input.sorted.dup.bam
data$H3K27me3_BT245_KO=data$BT245.ko.c2p3.Rx.C_cells_ChIP1_H3K27me3_1.sorted.dup.bam/data$BT245.ko.c2p3.Rx_cells_ChIP1_Input_1.sorted.dup.bam
data$H3K27me3_DIPG13_KO=data$DIPGXIII.ko.c5p5.Rx.C_cells_ChIP1_H3K27me3_1.sorted.dup.bam/data$DIPG13.ko.c5p5_Input.sorted.dup.bam
#hits
library(dplyr)
names(data)
data.range=makeGRangesFromDataFrame(data,seqnames.field = "seqnames",start.field = "Start",end.field = "End")
#region
region=subsetByOverlaps(promoter,gainedpeaks.k27m)
#bt245
hits=subsetByOverlaps(data.range,region)
hits=as.data.frame(hits)
hits$probe=paste(hits$seqnames,hits$start,hits$end)
m1=left_join(hits,data,by="probe")
results1=sum(m1$H2AK119ub_BT245)/sum(data$H2AK119ub_BT245)
#dipg13
results2=sum(m1$H2AK119ub_DIPG13)/sum(data$H2AK119ub_DIPG13)
#dip13 rep2
results10=sum(m1$H2AK119ub_BT245)/sum(data$H2AK119ub_BT245)
#pcgbm2
results5=sum(m1$H2AK119ub_DIPG13)/sum(data$H2AK119ub_DIPG13)
#g477
results6=sum(m1$H3K27me3_DIPG13_KO)/sum(data$H3K27me3_DIPG13_KO)
#bt245 ko 2
results7=sum(m1$H3K27me3_BT245_2)/sum(data$H3K27me3_BT245_2)
#dipg13 ko 2
#results8=sum(m1$H2AK119ub_DIPG13_KO_2)/sum(data$H2AK119ub_DIPG13_KO_2)
#dipg13 ko 3
#results9=sum(m1$H2AK119ub_DIPG13_KO_3)/sum(data$H2AK119ub_DIPG13_KO_3)
#load KO data
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
#bt245 ko1
data.range1=makeGRangesFromDataFrame(data2,seqnames.field = "seqnames",start.field = "Start",end.field = "End")
hits2=subsetByOverlaps(data.range1,region)
hits2=as.data.frame(hits2)
hits2$probe=paste(hits2$seqnames,hits2$start,hits2$end)
m2=left_join(hits2,data2,by="probe")
results3=sum(m2$H2AK119ub_BT245_KO)/sum(data2$H2AK119ub_BT245_KO)
#dipg13 ko1
results4=sum(m2$H2AK119ub_DIPG13_KO)/sum(data2$H2AK119ub_DIPG13_KO)

#results
View(data.frame(results=c(results1,results2,results3,results4,results5,results6,results7,results10)))


df2 <- data.frame(Lines=rep(c("BT245", "DIPGXIII"), each=2),
                  Legend=rep(c("Parental","K27M KO"),2),
                  len=c(results1*100,results3*100,results2*100,results4*100))
library(ggpubr)
library(ggthemes)
library(hrbrthemes)
#plot
ggplot(data=df2, aes(x=Lines, y=len, fill=Legend)) +
  geom_bar(stat="identity", position=position_dodge(),, color="black")+
  geom_text(aes(label=round(len,2)), vjust=1.6, color="white",
            position = position_dodge(0.9), size=3.5)+
  scale_fill_solarized()+
  xlab("")+ylab("SUZ12 reads, %")+ggtitle("Promoters containing CGIs depleted of H3K27me3 ")+
  theme_minimal()

ggsave("BarPlots_Quantification_SUZ12_promotersCpG_ncprc1.tiff", units="in", width=7, height=4, dpi=600, compression = 'lzw')

#RPKM of canonical vs noncanocial 
setwd("C:/Users/egj19/OneDrive - McGill University/Desktop/Analysis/BEDs")
library(readxl)
library(dplyr)
#load annoation
bt245_cprc1_anno=read_excel("bt245_cprc1.annotated.xlsx")
bt245_ncprc1_anno=read_excel("bt245_ncprc1.annotated.xlsx")
bt245_ko_cprc1_anno=read_excel("bt245_ko_cprc1.annotated.xlsx")
bt245_ko_ncprc1_anno=read_excel("bt245_ko_ncprc1.annotated.xlsx")
dipg13_cprc1_anno=read_excel("dipg13_cprc1.annotation.xlsx")
dipg13_ncprc1_anno=read_excel("dipg13_ncprc1.annotation.xlsx")
dipg13_ko_cprc1_anno=read_excel("dipg13_ko_cprc1.annotation.xlsx")
dipg13_ko_ncprc1_anno=read_excel("dipg13_ko_ncprc1.annotation.xlsx")
#load matrices RPKMs
BT245_RPKM=read.table("BT245_RPKM.txt",header=T)
BT245_RPKM$Parental=(BT245_RPKM$BT245_Parental1+BT245_RPKM$BT245_Parental2+BT245_RPKM$BT245_Parental3+BT245_RPKM$BT245_Parental4+BT245_RPKM$BT245_Parental5)/5
BT245_RPKM$KO=(BT245_RPKM$BT245_K27MKO_1+BT245_RPKM$BT245_K27MKO_2+BT245_RPKM$BT245_K27MKO_3+BT245_RPKM$BT245_K27MKO_4)/4
DIPG13_RPKM=read.table("DIPG13_RPKM.txt",header=T)
DIPG13_RPKM$Parental=(DIPG13_RPKM$BK.D13.DMSO.R1.bam+DIPG13_RPKM$BK.D13.DMSO.R2.bam+DIPG13_RPKM$BK.D13.DMSO.R3.bam+DIPG13_RPKM$C12.P2.bam+DIPG13_RPKM$C8.P2.bam)/5
DIPG13_RPKM$KO=(DIPG13_RPKM$BK.D13c5.DMSO.R1.bam+DIPG13_RPKM$BK.D13c5.DMSO.R2.bam+DIPG13_RPKM$C5.P2.bam+DIPG13_RPKM$C10.P2.bam)/4

#merge rpkms
bt245_cprc1_anno=left_join(bt245_cprc1_anno,BT245_RPKM,by="Probe")
bt245_cprc1_anno=bt245_cprc1_anno[complete.cases(bt245_cprc1_anno),]
bt245_ncprc1_anno=left_join(bt245_ncprc1_anno,BT245_RPKM,by="Probe")
bt245_ncprc1_anno=bt245_ncprc1_anno[complete.cases(bt245_ncprc1_anno),]
bt245_ko_cprc1_anno=left_join(bt245_ko_cprc1_anno,BT245_RPKM,by="Probe")
bt245_ko_cprc1_anno=bt245_ko_cprc1_anno[complete.cases(bt245_ko_cprc1_anno),]
bt245_ko_ncprc1_anno=left_join(bt245_ko_ncprc1_anno,BT245_RPKM,by="Probe")
bt245_ko_ncprc1_anno=bt245_ko_ncprc1_anno[complete.cases(bt245_ko_ncprc1_anno),]

dipg13_cprc1_anno=left_join(dipg13_cprc1_anno,DIPG13_RPKM,by="Probe")
dipg13_cprc1_anno=dipg13_cprc1_anno[complete.cases(dipg13_cprc1_anno),]
dipg13_ncprc1_anno=left_join(dipg13_ncprc1_anno,DIPG13_RPKM,by="Probe")
dipg13_ncprc1_anno=dipg13_ncprc1_anno[complete.cases(dipg13_ncprc1_anno),]
dipg13_ko_cprc1_anno=left_join(dipg13_ko_cprc1_anno,DIPG13_RPKM,by="Probe")
dipg13_ko_cprc1_anno=dipg13_ko_cprc1_anno[complete.cases(dipg13_ko_cprc1_anno),]
dipg13_ko_ncprc1_anno=left_join(dipg13_ko_ncprc1_anno,DIPG13_RPKM,by="Probe")
dipg13_ko_ncprc1_anno=dipg13_ko_ncprc1_anno[complete.cases(dipg13_ko_ncprc1_anno),]
#differtial peaks rpkms
DIPG13_gained_cprc1=as.data.frame(DIPG13_gained_cprc1)
DIPG13_gained_cprc1$link=paste(DIPG13_gained_cprc1$seqnames,DIPG13_gained_cprc1$start,DIPG13_gained_cprc1$end)
dipg13_cprc1_anno$link=paste(dipg13_cprc1_anno$Chromosome,dipg13_cprc1_anno$Start,dipg13_cprc1_anno$End)
DIPG13_gained_cprc1=left_join(DIPG13_gained_cprc1,dipg13_cprc1_anno,by="link")
DIPG13_gained_cprc1=DIPG13_gained_cprc1[complete.cases(DIPG13_gained_cprc1),]


DIPG13_maintained_cprc1=as.data.frame(DIPG13_maintained_cprc1)
DIPG13_maintained_cprc1$link=paste(DIPG13_maintained_cprc1$seqnames,DIPG13_maintained_cprc1$start,DIPG13_maintained_cprc1$end)
DIPG13_maintained_cprc1=left_join(DIPG13_maintained_cprc1,dipg13_cprc1_anno,by="link")
DIPG13_maintained_cprc1=DIPG13_maintained_cprc1[,-6]
DIPG13_maintained_cprc1=DIPG13_maintained_cprc1[complete.cases(DIPG13_maintained_cprc1),]

BT245_lost_cprc1=as.data.frame(BT245_lost_cprc1)
BT245_lost_cprc1$link=paste(BT245_lost_cprc1$seqnames,BT245_lost_cprc1$start,BT245_lost_cprc1$end)
bt245_ko_cprc1_anno$link=paste(bt245_ko_cprc1_anno$Chromosome.x,bt245_ko_cprc1_anno$Start.x,bt245_ko_cprc1_anno$End.x)
BT245_lost_cprc1=left_join(BT245_lost_cprc1,bt245_ko_cprc1_anno,by="link")
BT245_lost_cprc1=BT245_lost_cprc1[complete.cases(BT245_lost_cprc1),]

BT245_maintained_cprc1=as.data.frame(BT245_maintained_cprc1)
BT245_maintained_cprc1$link=paste(BT245_maintained_cprc1$seqnames,BT245_maintained_cprc1$start,BT245_maintained_cprc1$end)
BT245_maintained_cprc1=left_join(BT245_maintained_cprc1,bt245_ko_cprc1_anno,by="link")
BT245_maintained_cprc1=BT245_maintained_cprc1[,-6]
BT245_maintained_cprc1=BT245_maintained_cprc1[complete.cases(BT245_maintained_cprc1),]



#create df
mylist <-list(dipg13_cprc1_anno$Parental,dipg13_ncprc1_anno$Parental)
padding  <- function(a,b) {
  zz <-rep(NA,length(X)-length(Y))
  YY  <- c(Y,zz)
  out <- data.frame(cbind(a,YY))
}

Y=bt245_cprc1_anno$Parental
X=bt245_ncprc1_anno$Parental
dd3 <- padding(X,Y)
library(reshape)
names(dd1)=c("Non-canonical","Canonical")
names(dd2)=c("Non-canonical","Canonical")
names(dd3)=c("H3K27me3-depleted","H3K27me3-enriched")
names(dd4)=c("Maintained cPRC1","K27M KO-specifc cPRC1")

dd1=melt(dd1)
dd2=melt(dd2)
dd3=melt(dd3)
dd4=melt(dd4)

names(dd1)=c("PRC1","RPKM")
names(dd2)=c("PRC1","RPKM")
names(dd3)=c("PRC1","RPKM")
names(dd4)=c("PRC1","RPKM")

#violin plot
library("ggsci")
library("ggplot2")
library("gridExtra")
ggplot(dd3, aes(x=PRC1, y=RPKM,fill=PRC1)) + 
  geom_bixplot(trim=FALSE)+  
  labs(title="DIPG13 K27M KO - Expression of PRC1-bound genes",x="PRC1", y = "RPKMs")+ylim(0,50)+
  scale_fill_d3() + theme_classic()
#ggsave("ViolinPlot_DIPG13_K27MKO_PRC1_CanonicalVsNoncanonical_RPKMs.tiff", units="in", width=7, height=4, dpi=600, compression = 'lzw')
#run t-test
t.test(DIPG13_gained_cprc1$KO,DIPG13_gained_cprc1$Parental,paired=TRUE)

# Change box plot line colors by groups
ggplot(dd3, aes(x=PRC1, y=RPKM, fill=PRC1)) +
  geom_boxplot()+scale_fill_d3()+theme_classic()+ylim(0,20)+labs(title="BT245 Parental - Expression of PRC1-bound genes",x="PRC1", y = "RPKMs")
ggsave("BoxPlot_BT245_PRC1_CanonicalVsNoncanonical_RPKMs.tiff", units="in", width=7, height=4, dpi=600, compression = 'lzw')

#venn diagram for intersection
# Load library
library(VennDiagram)

# Generate 3 sets of 200 words
test1=as.list(bt245_ncprc1_anno$Probe,bt245_cprc1_anno$Probe)

set1 <-c(bt245_ncprc1_anno$Probe,bt245_cprc1_anno$Probe)
set1=set1[!is.na(set1)]
set2 <-c(bt245_ko_ncprc1_anno$Probe,bt245_ko_cprc1_anno$Probe)
set2=set2[!is.na(set2)]
set3 <-c(dipg13_cprc1_anno$Probe,dipg13_ncprc1_anno$Probe)
set3=set3[!is.na(set3)]
set4 <-c(dipg13_ko_ncprc1_anno$Probe,dipg13_ko_cprc1_anno$Probe)
set4=set4[!is.na(set4)]

#color

# Chart
venn.diagram(
  x = list(set1,set3,set2,set4),
  category.names = c("BT245 Parental" ,"DIPGXIII Parental", "BT245 K27M KO ", "DIPGXIII K27M KO" ),
  filename = 'Venndiagramm_RING1B_promoters_ParentalVsKO.png',
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 6000 , 
  width = 6000 , 
  resolution = 600,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = c("firebrick2","firebrick3","steelblue2","steelblue3"),
  
  # Numbers
  cex = 1,
  fontface = "bold",
  fontfamily = "sans",)
  


#how many CGIs fall within 5kb of promoter centers
cgi_promoters=subsetByOverlaps(cgi,subsetByOverlaps(cgi,promoter))
cgi_nonpromoters=setdiff(cgi,cgi_promoters)

length(subsetByOverlaps(cgi_promoters,BT245))/length(cgi_promoters)
length(subsetByOverlaps(cgi_nonpromoters,BT245))/length(cgi_nonpromoters)
View(as.data.frame(cgi))

#distance to cgi
library(GenomicRanges)
BT245_RING1b_distancetocgi=read.table("BT245_RING1b_distancetocgi.bed")
BT245_RING1b_distancetocgi$distance=abs(BT245_RING1b_distancetocgi$V13)
cluster1_distancetocgi=read.table("Cluster1_distancetocgi.bed")
cluster2_distancetocgi=read.table("Cluster2_distancetocgi.bed")
cluster3_distancetocgi=read.table("Cluster3_distancetocgi.bed")
cluster4_distancetocgi=read.table("Cluster4_distancetocgi.bed")
cluster5_distancetocgi=read.table("Cluster5_distancetocgi.bed")
K27M_RING1b_distancetocgi=read.table("K27M_RING1b_distancetocgi.bed")
KO_RING1b_distancetocgi=read.table("KO_RING1b_distancetocgi.bed")

#gained ring1b
allpeaks=makeGRangesFromDataFrame(BT245_RING1b_distancetocgi,seqnames.field = "V1",start.field = "V2",end.field = "V3")
allk27acpeaks=as.data.frame(BT245_K27ac)
allk27acpeaks$k27ac="k27ac"
hits=as.list(findOverlaps(BT245.gained,allpeaks))
table1=as.data.frame(BT245.gained)
table1$distance=sapply(extractList(BT245_RING1b_distancetocgi$distance,hits),sum)
hits2=as.list(findOverlaps(BT245.gained,BT245_K27ac))
table1$k27ac=sapply(extractList(allk27acpeaks$k27ac,hits2),paste)
table1=table1[table1$k27ac=="character(0)",]
#table prep
df1=table(cut(KO_RING1b_distancetocgi$V13,breaks=c(0,5000,50000,10000000000)))
#histogram plots distance to cgi
df <- data.frame(
  Distance=c("0-5Kb","5-50Kb","50Kb+"),
  Count=df1)
df$percentage=(df$Count.Freq/sum(df$Count.Freq))*100
library(ggpubr)
ggplot(df, aes(x=Distance,y=percentage)) +
  geom_bar(stat="identity",width=0.7, fill="dodgerblue2",color="black")+
  scale_x_discrete(limits=df$Distance)+ xlab("Distance from closest CGI")+ ylab("Percentage of Peaks")+
  ggtitle("KO RING1b Peak distance from closest CGIs")+
  theme_minimal()
#ggsave("BarPlot_Distancefromcgi_KO.tiff", units="in", width=7, height=4, dpi=600, compression = 'lzw')


#Chi-square of ring1b cgi enrichment
df=data.frame(line=c("BT245","DIPG13","BT245 KO","DIPG13 KO"),cgi=c(13214,7444,12634,4667),noncgi=c(1789,542,1914,592),canonical=c(1524,950,1772,1336),noncanonical=c(10922,6381,10659,3106))
chisq.test(df$cgi,df$noncgi,correct = TRUE)
chisq.test(df$canonical,df$noncanonical)

#k27m vs ko reads of prc1
df2 <- data.frame(Legend=rep(c("H3K27me3-enriched", "H3K27me3-depleted"), each=2),
                  Histone=rep(c("Wild-type","H3K27M"),2),
                  len=c(0.00755,0.03672,0.01105,0.01465))

df2$percentage=df2$len*100
ggplot(data=df2, aes(x=Legend, y=percentage, fill=Histone)) +
  geom_bar(stat="identity", position=position_dodge(),color="black")+
  geom_text(aes(label=round(percentage,2)), vjust=1.6, color="white",
            position = position_dodge(0.9), size=3.5)+
  xlab("")+ylab("SUZ12 reads, %")+ggtitle("Promoters containing CGIs ")+
  theme_minimal()+scale_fill_manual(values=c("firebrick2","dodgerblue2"))

ggsave("BarPlots_Quantification_SUZ12_promoters_K27MvsKO.tiff", units="in", width=7, height=4, dpi=600, compression = 'lzw')

#increased ring1b and hic at cprc1
setwd("C:/Users/egj19/OneDrive - McGill University/Desktop/Analysis/peakcalls")
bt245_cprc1=import.bed("bt245_cprc1.bed")
bt245_ko_cprc1=import.bed("bt245_ko_cprc1.bed")
BT245_loops=read.table("BT245_loops.txt")
BT245_loops$link=paste(BT245_loops$V1,BT245_loops$V2,BT245_loops$V3)
a=makeGRangesFromDataFrame(BT245_loops,seqnames.field = "V1",start.field = "V2",end.field = "V3")
BT245_ko_loops=read.table("BT245_c2_loops.txt")
b=makeGRangesFromDataFrame(BT245_ko_loops,seqnames.field = "V1",start.field = "V2",end.field = "V3")
gained_loops=setdiff(a,subsetByOverlaps(a,b))
region=subsetByOverlaps(bt245_cprc1,bt245_ko_cprc1)
region=subsetByOverlaps(gained_loops,region)
region=as.data.frame(region)
region$link=paste(region$seqnames,region$start,region$end)
region=left_join(region,BT245_loops,by="link")
region=region[complete.cases(region),]
setwd("C:/Users/egj19/OneDrive - McGill University/Desktop/Analysis/BEDs")
K27M_deg=read.table("K27M_DEG.txt",header=T)
hits=as.list(findOverlaps(makeGRangesFromDataFrame(region),makeGRangesFromDataFrame(bt245_cprc1_anno,seqnames.field = "Chromsome",start.field = "Start",end.field = "End")))
region$probe=paste(extractList(bt245_cprc1_anno$Probe,hits))
 

#size of genes intervals
genes=read.table("genes_hg19.bed")
genes1=makeGRangesFromDataFrame(genes,seqnames.field = "V1",start.field = "V2",end.field = "V3")
genes1=as.data.frame(genes1)
genes$promoter_start=genes$V2-2500
genes$promoter_end=genes$V2+2500


promoter=as.data.frame(promoter)
a=sum(promoter$width)*100/3101788170
b=sum(K27M_deg$length)*100/3101788170
c=100-(a+b)

#annotation adjustment to regions
pro.b.k27m.1=sum(31,1.915,1.1)/a
pro.d.k27m.1=sum(42.8984310,3.0435296,1.3861036)/a
pro.b.ko.1=sum(29.833,2.08,1.157)/a
pro.d.ko.1=sum(25.0356108,3.6003733,1.5766983)/a
inter.b.k27m.1=sum(33.355,1.012)/c
inter.d.k27m.1=sum(28.8545476,0.6016279)/c
inter.b.ko.1=sum(36.51,0.869)/c
inter.d.ko.1=sum(40.3261457,0.9037772)/c
intra.b.k27m.1=(100-(pro.b.k27m+inter.b.k27m))/b
intra.d.k27m.1=(100-(pro.d.k27m+inter.d.k27m))/b
intra.b.ko.1=(100-(pro.b.ko+inter.b.ko))/b
intra.d.ko.1=(100-(pro.d.ko+inter.d.ko))/b
df2 <- data.frame(Condition=rep(c("BT245", "BT245 K27M KO","DIPGXIII","DIPGXIII K27M KO"), each=3),
                  Annotation=rep(c("Promoters (within ±2.5Kb from TSS)", "Intergenic (outside ±2.5Kb from TSS)", "Genic"),4),
                  Percentage=c(0.870279573,0.071643395,0.058077032,0.864865481,0.079650456,0.055484063,0.920871634,0.046698511,0.032429855,0.848140089,0.094305507,0.057554404))
df2$Percentage=df2$Percentage*100
library(viridis)
library(hrbrthemes)
ggplot(df2, aes(fill=Annotation, y=Percentage, x=Condition)) + 
  geom_bar(position="stack", stat="identity") +
  ylab("RING1b peaks, %") +
  scale_fill_brewer(palette="YlGnBu") +
  ggtitle("RING1b Distribution - Normalized to region size") +
  theme_minimal() +
  xlab("Cell Lines")+
  geom_text(aes(label = round(Percentage,0)), size = 3, hjust = 0.5, vjust = 0, position = "stack") 

ggsave("BarPlots_RING1b_Annotation_2_normalized.tiff", units="in", width=7, height=4, dpi=600, compression = 'lzw')

#ashot plot
test=read.table("seqmonk53.txt",header=T)
test=makeGRangesFromDataFrame(test)
test=export.bed(test,"BT245_KO_CBX2_peaks_seqmonk.bed")
