#ChIP_RNAseq
setwd("C:/Users/egj19/OneDrive - McGill University/Desktop/Analysis")

#load data
library(readxl)
data <- read_excel("DEG_R.xlsx")

gained_up=intersect(data$genes_upregulated_k27m,data$RING1b_gainedpeaks_K27M)
gained_down=intersect(data$genes_downregulated_k27m,data$RING1b_gainedpeaks_K27M)
lost_up=intersect(data$RING1b_lostpeaks_K27M,data$genes_upregulated_k27m)
lost_down=intersect(data$RING1b_lostpeaks_K27M,data$genes_downregulated_k27m)
maintained_up=intersect(data$RING1b_maintainedpeaks_K27M,data$genes_upregulated_k27m)
maintained_down=intersect(data$RING1b_maintainedpeaks_K27M,data$genes_downregulated_k27m)



#retrieve coordinates of these DEG
genes=read.table("genes_coordinates2.txt",header=T)
genes=as.data.frame(genes)
gained_up=as.data.frame(gained_up)
gained_down=as.data.frame(gained_down)
lost_up=as.data.frame(lost_up)
lost_down=as.data.frame(lost_down)
names(gained_up)=c("Probe")
names(gained_down)=c("Probe")
names(lost_up)=c("Probe")
names(lost_down)=c("Probe")
gained_up=left_join(gained_up,genes)
gained_down=left_join(gained_down,genes)
lost_up=left_join(lost_up,genes)
lost_down=left_join(lost_down,genes)


#create a stacked barplot

# library
library(ggplot2)
library(viridis)
library(hrbrthemes)

# create a dataset
specie <- c(rep("Upregulated" , 3) , rep("Downregulated" , 3))
Sites <- rep(c("Gained RING1b" , "Lost RING1b" , "Maintained RING1b") , 2)
a1=length(gained_up)/sum(!is.na(data$genes_upregulated_k27m))*100
a2=length(lost_up)/sum(!is.na(data$genes_upregulated_k27m))*100
a3=length(maintained_up)/sum(!is.na(data$genes_upregulated_k27m))*100

a4=length(gained_down)/sum(!is.na(data$genes_downregulated_k27m))*100
a5=length(lost_down)/sum(!is.na(data$genes_downregulated_k27m))*100
a6=length(maintained_down)/sum(!is.na(data$genes_downregulated_k27m))*100


value <- abs(c(a1,a2,a3,a4,a5,a6))
data1 <- data.frame(specie,Sites,value)

# Small multiple
library(ggplot2)
library(viridis)
library(hrbrthemes)
BiocManager::install("paletteer")
library(paletteer)

p=ggplot(data1, aes(fill=Sites, y=value, x=specie,label=value)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_viridis(discrete = T) +
  ggtitle("Differential RING1b recruitment compared with Differentially Expressed Genes in K27M vs KO") +
  theme_ipsum() +
  xlab("") +
  ylab("Percentage of Differentially Expressed Genes") +
  theme(plot.title = element_text(color="black", size=25, face="bold"),axis.title.y = element_text(vjust = 5,hjust=0.4,size = 20),axis.text.x = element_text(color="black", size=20, face="italic"))
  
p + geom_text(aes(label = round(value,2)), size = 10, hjust = 0.5, vjust = 3, position ="stack",colour="white")

#BT245 analysis
setwd("C:/Users/egj19/OneDrive - McGill University/Desktop/Analysis")
BT245_deg=read.table("BT245_DEG.txt",header=T)
BT245_deg_up=makeGRangesFromDataFrame(filter(BT245_deg,BT245_deg$Log2_Fold_Change<c(-1)))
BT245_deg_down=filter(BT245_deg,BT245_deg$Log2_Fold_Change>c(1))
BT245_maintainedpeaks=subsetByOverlaps(BT245.1,BT245.ko.1)       
BT245_gainedpeaks=setdiff(BT245.1,subsetByOverlaps(BT245.1,BT245_maintainedpeaks))
length(subsetByOverlaps(BT245.1,BT245_maintainedpeaks))+length(BT245_gainedpeaks)-length(BT245.1)
BT245_lostpeaks=setdiff(BT245.ko.1,subsetByOverlaps(BT245.ko.1,BT245_maintainedpeaks))
length(subsetByOverlaps(BT245.ko.1,BT245_maintainedpeaks))+length(BT245_lostpeaks)-length(BT245.ko.1)
intersect(BT245_gainedpeaks,BT245_lostpeaks)
#prepare ChIP_RNAseq dataset
setwd("C:/Users/egj19/OneDrive - McGill University/Desktop/Analysis/Peakcalls")
BT245_RNA=read.table("seqmonk17.txt",header=T)
BT245_RNA$BT245=(BT245_RNA$BT245.NKO__BT245.NKO.P17.sorted.bam+BT245_RNA$BT245.NKO__BT245.NKO.P18.sorted.bam+BT245_RNA$BT245.NKO__BT245.NKO.P30.sorted.bam)/3 
BT245_RNA$BT245_KO=(BT245_RNA[5]+BT245_RNA[6]+BT245_RNA[7]+BT245_RNA[8])/4
BT245_exp=BT245_RNA[BT245_RNA$BT245<c(1),]
BT245_KO_exp=BT245_RNA[BT245_RNA$BT245_KO<c(1),]

#order-sensitive
BT245_KO_RING1b_genes=read.table("BT245_C2P8_RING1b_peaks_seqmonk_genes.txt",header=T)
BT245_RING1b_genes=read.table("BT245_RING1b_peaks_seqmonk_genes.txt",header=T)
BT245_exp=makeGRangesFromDataFrame(BT245_exp)
BT245_RING1b_genes=makeGRangesFromDataFrame(BT245_RING1b_genes)
BT245_RNA=makeGRangesFromDataFrame(BT245_RNA)
BT245_sil=setdiff(BT245_RNA,BT245_exp)
a1=length(findOverlaps(BT245_RING1b_genes,BT245_exp))/length(BT245_RING1b_genes)

BT245_KO_exp=makeGRangesFromDataFrame(BT245_KO_exp)
BT245_KO_RING1b_genes=makeGRangesFromDataFrame(BT245_KO_RING1b_genes)
BT245_RNA=makeGRangesFromDataFrame(BT245_RNA)
BT245_KO_sil=setdiff(BT245_RNA,BT245_KO_exp)

a2=length(findOverlaps(BT245_KO_RING1b_genes,BT245_KO_exp))/length(BT245_KO_RING1b_genes)
#combine by nature of activity of RING1b
library(rtracklayer)
DEG_RING1b_repressive=read.table("DEG_RING1b_repressive.txt")
DEG_RING1b_repressive=makeGRangesFromDataFrame(DEG_RING1b_repressive,seqnames.field = "V1",start.field = "V2",end.field = "V3")
export.bed(DEG_RING1b_repressive,"DEG_RING1b_repressive.bed")

DEG_RING1b_activating=read.table("DEG_RING1b_activating.txt")
DEG_RING1b_activating=makeGRangesFromDataFrame(DEG_RING1b_activating,seqnames.field = "V1",start.field = "V2",end.field = "V3")
export.bed(DEG_RING1b_activating,"DEG_RING1b_activating.bed")
#match by names
BT245_RNA=read.table("seqmonk17.txt",header=T)
BT245_RNA$BT245=(BT245_RNA$BT245.NKO__BT245.NKO.P17.sorted.bam+BT245_RNA$BT245.NKO__BT245.NKO.P18.sorted.bam+BT245_RNA$BT245.NKO__BT245.NKO.P30.sorted.bam)/3 
BT245_RNA$BT245_KO=(BT245_RNA[5]+BT245_RNA[6]+BT245_RNA[7]+BT245_RNA[8])/4
BT245_exp=BT245_RNA[BT245_RNA$BT245>1,]
BT245_KO_exp=BT245_RNA[BT245_RNA$BT245_KO>1,]
BT245_KO_RING1b_genes=read.table("BT245_C2P8_RING1b_peaks_seqmonk_genes.txt",header=T)
BT245_RING1b_genes=read.table("BT245_RING1b_peaks_seqmonk_genes.txt",header=T)

a3=left_join(BT245_RING1b_genes,BT245_exp,by="Probe")
a4=sum(complete.cases(a3$Chromosome.y))/length(BT245_RING1b_genes$Probe)
a5=left_join(BT245_KO_RING1b_genes,BT245_KO_exp,by="Probe")
a6=sum(complete.cases(a5$Chromosome.y))/length(BT245_KO_RING1b_genes$Probe)
a7=left_join(BT245_RING1b_genes,BT245_RNA,by="Probe")
a8=sum(complete.cases(a7$Chromosome.y))/length(BT245_RING1b_genes$Probe)
a4=a4/a8
a6=a6/a8


#Use Nick's promoters/enhacers lists based on k27ac and DEG
library("rtracklayer")
library("GenomicRanges")
library("ChIPseeker")
library(dplyr)
library(readxl)
chipseeker <- read_excel("chipseeker.xlsx")

BT245.1 <- readPeakFile(chipseeker[[1]])
BT245=makeGRangesFromDataFrame(BT245.1)

DIPG13.1 <- readPeakFile(chipseeker[[3]])
DIPG13=makeGRangesFromDataFrame(DIPG13.1)

BT245.ko.1 <- readPeakFile(chipseeker[[2]])
BT245.ko=makeGRangesFromDataFrame(BT245.ko.1)

DIPG13.ko.1 <- readPeakFile(chipseeker[[3]])
DIPG13.ko=makeGRangesFromDataFrame(DIPG13.ko.1)

k27m=subsetByOverlaps(BT245,DIPG13)
k27m
ko=subsetByOverlaps(DIPG13.ko.1,BT245.ko.1)
ko
#differential peaks
maintainedring1b=subsetByOverlaps(k27m,ko)
gainedpeaks.k27m=setdiff(k27m,subsetByOverlaps(k27m,maintainedring1b))
lostpeaks.k27m=setdiff(ko,subsetByOverlaps(ko,k27m))


enhancers <- read_excel("nick_enhancers_k27ac_rna.xlsx")
promoters<-read.table("promoters.ucsc",header=F)
enhancers=makeGRangesFromDataFrame(enhancers)
promoters=unique(promoters[,c(1,2,3)])
promoters=makeGRangesFromDataFrame(promoters,start.field = "V2",seqnames.field = "V1",end.field = "V3")
promoters=setdiff(promoters,subsetByOverlaps(promoters,enhancers))
enhancers=setdiff(enhancers,subsetByOverlaps(enhancers,promoters))
intersect(promoters,enhancers)

length(subsetByOverlaps(k27m,enhancers))/length(k27m)
length(subsetByOverlaps(ko,enhancers))/length(ko)

length(subsetByOverlaps(k27m,promoters))/length(k27m)
length(subsetByOverlaps(ko,promoters))/length(ko)


#DEG_RING1b_peaks with other marks
DEG_RING1b_activating=import.bed("DEG_RING1b_activating.bed")
DEG_RING1b_repressive=import.bed("DEG_RING1b_repressive.bed")

subsetByOverlaps(makeGRangesFromDataFrame(lost_up[c(1:82),]),gainedk27ac.k27m)

#gained promoters
library(GenomicRanges)
library(rtracklayer)
setwd("C:/Users/egj19/OneDrive - McGill University/Desktop/Analysis/BEDs")
bt245_cprc1=import.bed("bt245_cprc1.bed")
bt245_ko_cprc1=import.bed("bt245_ko_cprc1.bed")
dipg13_cprc1=import.bed("dipg13_cprc1.bed")
dipg13_ko_cprc1=import.bed("dipg13_ko_cprc1.bed")
setwd("C:/Users/egj19/OneDrive - McGill University/Desktop/Analysis/peakcalls")
BT245=import.bed("BT245_RING1b_peaks_optimized.bed")
DIPG13=import.bed("DIPG13_RING1b_peaks_optimized.bed")
BT245.ko=import.bed("BT245_C2P8_RING1b_peaks_optimized.bed")
DIPG13.ko=import.bed("DIPG13_C5P6_RING1b_peaks_optimized.bed")

K27M=subsetByOverlaps(BT245,DIPG13)
KO=subsetByOverlaps(BT245.ko,DIPG13.ko)

gainedpeaks.bt245=setdiff(BT245,subsetByOverlaps(BT245,BT245.ko))
lostpeaks.bt245=setdiff(BT245.ko,subsetByOverlaps(BT245.ko,BT245))
maintainedpeaks.bt245=subsetByOverlaps(BT245,BT245.ko)

gainedpeaks.dipg13=setdiff(DIPG13,subsetByOverlaps(DIPG13,DIPG13.ko))
lostpeaks.dipg13=setdiff(DIPG13.ko,subsetByOverlaps(DIPG13.ko,DIPG13))
maintainedpeaks.dipg13=subsetByOverlaps(DIPG13,DIPG13.ko)

gainedpeaks.k27m=subsetByOverlaps(gainedpeaks.bt245,gainedpeaks.dipg13)
lostpeaks.k27m=subsetByOverlaps(lostpeaks.bt245,lostpeaks.dipg13)
maintainedpeaks=subsetByOverlaps(maintainedpeaks.bt245,maintainedpeaks.dipg13)

library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <-TxDb.Hsapiens.UCSC.hg19.knownGene
promoter <- getPromoters(TxDb=txdb, upstream=2500, downstream=2500)
#promoters with gained/lost/maintained peaks
library(readxl)
library(dplyr)
setwd("C:/Users/egj19/OneDrive - McGill University/Desktop/Analysis/BEDs")

#gained_promoters.anno=subsetByOverlaps(promoter,gainedpeaks.k27m)
#lost_promoters.anno=subsetByOverlaps(promoter,lostpeaks.k27m)
#maintained_promoters.anno=subsetByOverlaps(promoter,maintainedpeaks)

gained_promoters.anno=read_excel("gained_promoters.annotation.xlsx")
lost_promoters.anno=read_excel("lost_promoters.annotation.xlsx")
maintained_promoters.anno=read_excel("maintained_promoters.annotation.xlsx")
#import K27M data
K27M_DEG=read.table("K27M_DEG.txt",header=T)

gained_promoters.anno=left_join(gained_promoters.anno,K27M_DEG,by="Probe")
gained_promoters.anno=gained_promoters.anno[complete.cases(gained_promoters.anno),]

lost_promoters.anno=left_join(lost_promoters.anno,K27M_DEG,by="Probe")
lost_promoters.anno=lost_promoters.anno[complete.cases(lost_promoters.anno),]

maintained_promoters.anno=left_join(maintained_promoters.anno,K27M_DEG,by="Probe")
maintained_promoters.anno=maintained_promoters.anno[complete.cases(maintained_promoters.anno),]
#plot bargraphs
a=length(makeGRangesFromDataFrame(gained_promoters.anno,seqnames.field = "Chromosome.x",start.field = "Start.x",end.field = "End.x"))
b=length(makeGRangesFromDataFrame(lost_promoters.anno,seqnames.field = "Chromosome.x",start.field = "Start.x",end.field = "End.x"))
c=length(makeGRangesFromDataFrame(maintained_promoters.anno,seqnames.field = "Chromosome.x",start.field = "Start.x",end.field = "End.x"))
test1=gained_promoters.anno[gained_promoters.anno$P.adj<0.1,]
test2=lost_promoters.anno[lost_promoters.anno$P.adj<0.1,]
test3=maintained_promoters.anno[maintained_promoters.anno$P.adj<0.1,]

a1=length(makeGRangesFromDataFrame(test1[test1$Log2_FC<c(-1),],seqnames.field = "Chromosome.x",start.field = "Start.x",end.field = "End.x"))
b1=length(makeGRangesFromDataFrame(test2[test2$Log2_FC<c(-1),],seqnames.field = "Chromosome.x",start.field = "Start.x",end.field = "End.x"))
c1=length(makeGRangesFromDataFrame(test3[test3$Log2_FC<c(-1),],seqnames.field = "Chromosome.x",start.field = "Start.x",end.field = "End.x"))

a2=length(makeGRangesFromDataFrame(test1[test1$Log2_FC>c(1),],seqnames.field = "Chromosome.x",start.field = "Start.x",end.field = "End.x"))
b2=length(makeGRangesFromDataFrame(test2[test2$Log2_FC>c(1),],seqnames.field = "Chromosome.x",start.field = "Start.x",end.field = "End.x"))
c2=length(makeGRangesFromDataFrame(test3[test3$Log2_FC>c(1),],seqnames.field = "Chromosome.x",start.field = "Start.x",end.field = "End.x"))

test=rbind(test1,test2)
test=test[abs(test$Log2_FC)>c(1),]
test=makeGRangesFromDataFrame(test,seqnames.field = "Chromosome.x",start.field = "Start.x",end.field = "End.x")
export.bed(test,"test.bed")
#create df
df <- data.frame(Promoters=c("Gained", "Lost", "Maintained"),
                 total=c(a,b,c),Upregulated=c(a1,b1,c1),Downregulated=c(a2,b2,c2))
df$Maintained=df$total-(df$Upregulated+df$Downregulated)
library(reshape2)
df=melt(df,id.vars=c("Promoters","total"))
df$percentage=round((df$value/df$total)*100,1)
library(dplyr)
library(ggpubr)
library(forcats)
df=df[order(df$value),]
df<- ddply(df, "Promoters",
           transform, label_ypos=cumsum(percentage))

#plot
ggplot(data=df,aes(x=Promoters, y=percentage, fill=factor(df$variable,levels=c("Maintained","Upregulated","Downregulated")),label=value)) +
  geom_bar(stat="identity")+
  ylab("Promoters, %")+guides(fill = guide_legend(title=NULL))+
  scale_fill_manual(values=c("gray88","olivedrab4","darkviolet"), 
                    breaks=c("Maintained","Upregulated","Downregulated"),
                    labels=c("No signficant change","Upregulated in K27M","Downregulated in K27M"))+
  theme_classic()+xlab("RING1b at promoters (numbers)")+scale_x_discrete(labels=c("Gained (238)", "Lost (171)", "Maintained (2535)"))+
    theme(axis.text.x = element_text(color="black",size=9))+geom_text(size=3,position = position_stack(vjust = 1.01))
  
#ggsave("BarPlot_K27M_DEG_RING1b.tiff", units="in", width=7, height=4, dpi=600, compression = 'lzw')

#understand gained prc1 peaks
k27m.k27me3=subsetByOverlaps(BT245_K27me3,DIPG13_K27me3)
ko.k27me3=subsetByOverlaps(BT245_K27MKO_K27me3,DIPG13_K27MKO_K27me3)
k27m.k27ac=suppressWarnings(subsetByOverlaps(BT245_K27ac,DIPG13_K27ac))
ko.k27ac=subsetByOverlaps(BT245_K27MKO_K27ac,DIPG13_K27MKO_K27ac)

makeGRangesFromDataFrame(gained_promoters.anno[gained_promoters.anno$Log2_FC<c(-1),],seqnames.field = "Chromosome.x",start.field = "Start.x",end.field = "End.x")
r1=subsetByOverlaps(makeGRangesFromDataFrame(test2[test2$Log2_FC<c(-1),],seqnames.field = "Chromosome.x",start.field = "Start.x",end.field = "End.x"),ko.k27me3)
r2=subsetByOverlaps(makeGRangesFromDataFrame(test2[test2$Log2_FC<c(-1),],seqnames.field = "Chromosome.x",start.field = "Start.x",end.field = "End.x"),ko.k27ac)
length(r1)+length(r2)

r3=subsetByOverlaps(makeGRangesFromDataFrame(test2[test2$Log2_FC>c(1),],seqnames.field = "Chromosome.x",start.field = "Start.x",end.field = "End.x"),ko.k27me3)
r4=subsetByOverlaps(makeGRangesFromDataFrame(test2[test2$Log2_FC>c(1),],seqnames.field = "Chromosome.x",start.field = "Start.x",end.field = "End.x"),ko.k27ac)
length(r3)+length(r4)


#plot heatmap
library(pheatmap)
library(dplyr)
test=test2[abs(test2$Log2_FC)>c(1),]
test=test[order(test$Log2_FC),]
names(test)
test=test[,c(1:10,23:24,16:20,21:22,11:15)]
test=as.data.frame(test)
rownames(test)=test$Probe
test=test[,c(11:24)]

#names(test)=paste(c("BT245 K27M KO 1","BT245 K27M KO 2","BT245 K27M KO 3","BT245 K27M KO 4","BT245 K27M KO 5","BT245 Parental 1","BT245 Parental 2","BT245 Parental 3","BT245 Parental 4","BT245 Parental 5"))
names(test)=paste(c("DIPGXIII Parental 1","DIPGXIII Parental 2","BT245 Parental 1","BT245 Parental 2","BT245 Parental 3","BT245 Parental 4","BT245 Parental 5","DIPGXIII K27M KO 1","DIPGXIII K27M KO 2","BT245 K27M KO 1","BT245 K27M KO 2","BT245 K27M KO 3","BT245 K27M KO 4","BT245 K27M KO 5"))
test=as.matrix(test)


cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}
data_subset_norm <- t(apply(test, 1, cal_z_score))

jpeg("Heatmap_LostPRC1_DEG.jpeg",res=1200,width = 10, height = 8, units = 'in')

pheatmap(test,scale = "row",fontsize_row = 10,cluster_cols = FALSE,cluster_rows = FALSE,main="Lost PRC1 peaks in K27M and associated change in expression")

dev.off()

#volcano plots of cPRC1 vs ncPRC1
K27M_DEG
jpeg("BT245_gainedpeaks_volcanoplot.jpeg",res=600,width = 10, height = 8, units = 'in')
EnhancedVolcano(res,
                lab = res$Probe,
                x = "Log2_FC",
                y = "FDR",
                pCutoff = 0.1,
                xlab=bquote(~Log[2]~"fold change"),
                ylab=bquote(~-Log[10]~adjusted~italic(P)),
                xlim = c(-12, +12),
                title="BT245 K27M ectopic peaks - K27M vs K27M KO",
                pointSize=1.0,
                labSize=2.0,
                colAlpha = 1)
dev.off()



#cPRC1, ncPRC1 and DEGs
library(pheatmap)
library(readxl)
library(dplyr)
setwd("C:/Users/egj19/OneDrive - McGill University/Desktop/Analysis/BEDs")
cprc1.anno=read_excel("cprc1.annotation.xlsx")
ncprc1.anno=read_excel("ncprc1.annotation.xlsx")

K27M_DEG=read.table("K27M_DEG.txt",header=T)
cprc1.anno=left_join(cprc1.anno,K27M_DEG,by="Probe")
cprc1.anno=cprc1.anno[complete.cases(cprc1.anno),]
test=cprc1.anno[cprc1.anno$P.adj<0.1,]
test=test[abs(test$Log2_FC)>1,]
ncprc1.anno=left_join(ncprc1.anno,K27M_DEG,by="Probe")
ncprc1.anno=ncprc1.anno[complete.cases(ncprc1.anno),]
test1=ncprc1.anno[ncprc1.anno$P.adj>0.1,]
test1=test1[abs(test1$Log2_FC)>1,]
names(test1)
test=test[,c(1:10,23:24,16:20,21:22,11:15)]
test=as.data.frame(test)
test=test[order(test$Log2_FC),]
test=test[!duplicated(test$Probe),]
rownames(test)=test$Probe
test=test[,c(11:24)]
names(test)=paste(c("DIPGXIII Parental 1","DIPGXIII Parental 2","BT245 Parental 1","BT245 Parental 2","BT245 Parental 3","BT245 Parental 4","BT245 Parental 5","DIPGXIII K27M KO 1","DIPGXIII K27M KO 2","BT245 K27M KO 1","BT245 K27M KO 2","BT245 K27M KO 3","BT245 K27M KO 4","BT245 K27M KO 5"))
test=as.matrix(test)

cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}
data_subset_norm <- t(apply(test, 1, cal_z_score))

#jpeg("Heatmap_cPRC1_DEG.jpeg",res=1200,width = 10, height = 8, units = 'in')
(col.pal <- rev(RColorBrewer::brewer.pal(9, "RdBu")))
pheatmap(test,scale = "row",fontsize_row =1,cluster_cols = FALSE,cluster_rows = FALSE,main="cPRC1 genes and expression",color=col.pal,show_rownames     = FALSE,)
#dev.off()

#cprc1 rpkm
library(dplyr)
BT245_RPKM=read.table("BT245_RPKM.txt",header=T)
DIPG13_RPKM=read.table("DIPG13_RPKM.txt",header=T)
cprc1.anno=read_excel("cprc1.annotation.xlsx")
ncprc1.anno=read_excel("ncprc1.annotation.xlsx")
cprc1.anno=left_join(cprc1.anno,BT245_RPKM,by="Probe")
cprc1.anno=left_join(cprc1.anno,DIPG13_RPKM,by="Probe")
ncprc1.anno=left_join(ncprc1.anno,BT245_RPKM,by="Probe")
ncprc1.anno=left_join(ncprc1.anno,DIPG13_RPKM,by="Probe")
ncprc1.anno=ncprc1.anno[complete.cases(ncprc1.anno),]

cprc1.anno$K27M=(cprc1.anno$BT245_Parental1+cprc1.anno$BT245_Parental2+cprc1.anno$BT245_Parental3+cprc1.anno$BT245_Parental4+cprc1.anno$BT245_Parental5+cprc1.anno$BK.D13.DMSO.R2.bam+cprc1.anno$BK.D13.DMSO.R3.bam+cprc1.anno$BK.D13.DMSO.R1.bam)/8
cprc1.anno$KO=(cprc1.anno$C12.P2.bam+cprc1.anno$C8.P2.bam+cprc1.anno$BK.D13c5.DMSO.R1.bam+cprc1.anno$BK.D13c5.DMSO.R2.bam+cprc1.anno$C5.P2.bam+cprc1.anno$C10.P2.bam+cprc1.anno$BT245_Parental1+cprc1.anno$BT245_Parental2+cprc1.anno$BT245_Parental3+cprc1.anno$BT245_Parental4+cprc1.anno$BT245_Parental5)/11
cprc1.anno$BT245=(cprc1.anno$BT245_Parental1+cprc1.anno$BT245_Parental2+cprc1.anno$BT245_Parental3+cprc1.anno$BT245_Parental4+cprc1.anno$BT245_Parental5)/5
cprc1.anno$BT245_KO=(cprc1.anno$BT245_K27MKO_1+cprc1.anno$BT245_K27MKO_2+cprc1.anno$BT245_K27MKO_3+cprc1.anno$BT245_K27MKO_4+cprc1.anno$BT245_K27MKO_5)/5
cprc1.anno$DIPG13=(cprc1.anno$DIPGXIII.NKO__DIPGXIII.NKO.C8.sorted.bam+cprc1.anno$DIPGXIII.NKO__DIPGXIII.NKO.C12.sorted.bam)/2
cprc1.anno$DIPG13_KO=(cprc1.anno$DIPGXIII.KO__DIPGXIII.KO.C5.sorted.bam+cprc1.anno$DIPGXIII.KO__DIPGXIII.KO.C10.sorted.bam)/2
ncprc1.anno$K27M=(ncprc1.anno$BT245_Parental1+ncprc1.anno$BT245_Parental2+ncprc1.anno$BT245_Parental3+ncprc1.anno$BT245_Parental4+ncprc1.anno$BT245_Parental5+ncprc1.anno$BK.D13.DMSO.R2.bam+ncprc1.anno$BK.D13.DMSO.R3.bam+ncprc1.anno$BK.D13.DMSO.R1.bam)/8
ncprc1.anno$KO=(ncprc1.anno$C12.P2.bam+ncprc1.anno$C8.P2.bam+ncprc1.anno$BK.D13c5.DMSO.R1.bam+ncprc1.anno$BK.D13c5.DMSO.R2.bam+ncprc1.anno$C5.P2.bam+ncprc1.anno$C10.P2.bam+ncprc1.anno$BT245_Parental1+ncprc1.anno$BT245_Parental2+ncprc1.anno$BT245_Parental3+ncprc1.anno$BT245_Parental4+ncprc1.anno$BT245_Parental5)/11
ncprc1.anno$BT245=(ncprc1.anno$BT245_Parental1+ncprc1.anno$BT245_Parental2+ncprc1.anno$BT245_Parental3+ncprc1.anno$BT245_Parental4+ncprc1.anno$BT245_Parental5)/5
ncprc1.anno$BT245_KO=(ncprc1.anno$BT245_K27MKO_1+ncprc1.anno$BT245_K27MKO_2+ncprc1.anno$BT245_K27MKO_3+ncprc1.anno$BT245_K27MKO_4+ncprc1.anno$BT245_K27MKO_5)/5
ncprc1.anno$DIPG13=(ncprc1.anno$DIPGXIII.NKO__DIPGXIII.NKO.C8.sorted.bam+ncprc1.anno$DIPGXIII.NKO__DIPGXIII.NKO.C12.sorted.bam)/2
ncprc1.anno$DIPG13_KO=(ncprc1.anno$DIPGXIII.KO__DIPGXIII.KO.C5.sorted.bam+ncprc1.anno$DIPGXIII.KO__DIPGXIII.KO.C10.sorted.bam)/2


dd3=data.frame(Line=c('Parental',"K27M KO"),RPKM=c(ncprc1.anno$BT245,ncprc1.anno$BT245_KO))


#violin plot
library(ggplot2)
ggplot(dd3, aes(x=Line, y=RPKM,fill=Line)) + 
  geom_boxplot()+ 
  labs(title="Expression levels of ncPRC1-bound genes",x="PRC1", y = "RPKMs")+ylim(0,10000)+
  scale_fill_manual(values=rev(c("firebrick3","dodgerblue2"))) + theme_classic()

#ggsave("BoxPlot_cPRC1_K27MVsKO_RPKMs.tiff", units="in", width=7, height=4, dpi=600, compression = 'lzw')
t.test(cprc1.anno$K27M,cprc1.anno$KO,paired=F)
summary(ncprc1.anno$K27M)
summary(ncprc1.anno$KO)

#different cbx expressions
names(cbx2)
cbx2=K27M_DEG[K27M_DEG$Probe=="CBX2",]
cbx4=K27M_DEG[K27M_DEG$Probe=="CBX4",]
cbx6=K27M_DEG[K27M_DEG$Probe=="CBX6",]
cbx7=K27M_DEG[K27M_DEG$Probe=="CBX7",]
cbx8=K27M_DEG[K27M_DEG$Probe=="CBX8",]
cbx=rbind(cbx2,cbx4,cbx6,cbx7,cbx8)
cbx=cbx[,c(1:7,20:21,13:17,18:19,8:12)]
names(cbx)=paste(c("Probe","Chromosome","Start","End","P.adj","FDR","Log2_FC","DIPGXIII Parental 1","DIPGXIII Parental 2","BT245 Parental 1","BT245 Parental 2","BT245 Parental 3","BT245 Parental 4","BT245 Parental 5","DIPGXIII K27M KO 1","DIPGXIII K27M KO 2","BT245 K27M KO 1","BT245 K27M KO 2","BT245 K27M KO 3","BT245 K27M KO 4","BT245 K27M KO 5"))
df.bt245=c(cbx=rep(c("CBX2","CBX4","CBX6","CBX7","CBX8"),5),RPKMs=c(as.numeric(cbx[1,10:14]),as.numeric(cbx[2,10:14]),as.numeric(cbx[3,10:14]),as.numeric(cbx[4,10:14])))


df=read_excel("CBX_RPKM.xlsx")
df=melt(df,id.vars = c("CBX","Tissue"))
#violin plot
ggplot(df, aes(x=CBX, y=RPKM,fill=`Cell Line`)) + 
  geom_boxplot()+  
  labs(title="RNA Expression of CBX proteins",x="CBX", y = "RPKMs")+
  scale_fill_manual(values=c("firebrick1","dodgerblue1","firebrick3","dodgerblue3"))+ theme_classic()
ggsave("BoxPlot_CBX_expression.tiff", units="in", width=8, height=4, dpi=600, compression = 'lzw')

#differenr PRC1 components
#different cbx expressions
a=K27M_DEG[K27M_DEG$Probe=="",]
b=K27M_DEG[K27M_DEG$Probe=="CBX2",]
c=K27M_DEG[K27M_DEG$Probe=="CBX4",]
d=K27M_DEG[K27M_DEG$Probe=="CBX6",]
e=K27M_DEG[K27M_DEG$Probe=="CBX7",]
f=K27M_DEG[K27M_DEG$Probe=="CBX8",]
g=K27M_DEG[K27M_DEG$Probe=="MGA",]
h=K27M_DEG[K27M_DEG$Probe=="HDAC1",]
i=K27M_DEG[K27M_DEG$Probe=="HDAC2",]
all=rbind(a,b,c,d,e,f,g,h,i)
all=all[,c(1:7,20:21,13:17,18:19,8:12)]
names(all)=paste(c("Probe","Chromosome","Start","End","P.adj","FDR","Log2_FC","DIPGXIII Parental 1","DIPGXIII Parental 2","BT245 Parental 1","BT245 Parental 2","BT245 Parental 3","BT245 Parental 4","BT245 Parental 5","DIPGXIII K27M KO 1","DIPGXIII K27M KO 2","BT245 K27M KO 1","BT245 K27M KO 2","BT245 K27M KO 3","BT245 K27M KO 4","BT245 K27M KO 5"))
library(reshape2)
#PRC1_components=read_excel("PRC1_components.xlsx")
#PRC1_components.melted=melt(PRC1_components,id.vars=c("Gene","PRC1"))
#PRC1_components.melted=write.table(PRC1_components.melted,"PRC1_components.melted.txt",sep="\t")

PRC1_components.melted=read_excel("PRC1_components.melted.xlsx")

PRC1_components.melted=as.data.frame(PRC1_components.melted)
df=PRC1_components.melted
[
  PRC1_components.melted$PRC1=="cPRC1",]

ggplot(df, aes(x=Gene, y=RPKM,fill=`Cell Line`)) + 
  geom_boxplot()+  
  labs(title="RNA Expression of Non-canonical PRC1 components",x="Components", y = "RPKMs")+
  scale_fill_manual(values=c("firebrick1","dodgerblue1","firebrick3","dodgerblue3"))+ theme_classic()

ggsave("BoxPlot_PRC1_Noncanonical_components_expression.tiff", units="in", width=8, height=4, dpi=600, compression = 'lzw')


#heatmap of PRC1 components
library(ggplot2)
library(hrbrthemes)

ggplot(df, aes(df$`Cell Line`, Gene, fill= log(df$RPKM))) + 
  geom_tile() +
  scale_fill_gradient(low="blue", high="red") 
  
ggplot(df, aes(df$`Cell Line`, Gene)) + geom_tile(aes(fill = scale(df$RPKM)),colour = "white",size=0) + scale_fill_gradient(low = "steelblue2",high = "firebrick2")+theme_ipsum()


#cbx2 bed file dipg13
setwd("C:/Users/egj19/OneDrive - McGill University/Desktop/Analysis/peakcalls")
test=read.table("seqmonk45.txt",header=T)
test=makeGRangesFromDataFrame(test)
export.bed(test,"HSJ51_RING1b_peaks_seqmonk.bed")

#t.test
a=c(0.29467483,
    0.32665711,
    0.1895651,
    0.26849441
    
    
)

b=c(0.11249421,
    0.18383456,
    0.08959224,
    0.15489142
    
    
)
t.test(a,b)

#BT245 Gained peaks at promoters and association
BT245_gainedpeaks=setdiff(BT245,subsetByOverlaps(BT245, BT245.ko))
DIPG13_gainedpeaks=setdiff(DIPG13,subsetByOverlaps(DIPG13,DIPG13.ko))
#at promoters (1885), (4962)
BT245_gainedpromoters=subsetByOverlaps(promoter,BT245_gainedpeaks)
DIPG13_gainedpromoters=subsetByOverlaps(promoter,DIPG13_gainedpeaks)
#K27me3 association (311),(655)
set1=subsetByOverlaps(BT245_gainedpromoters,BT245_H3K27me3_10kb)
set10=subsetByOverlaps(DIPG13_gainedpromoters,DIPG13_H3K27me3_10kb)
#K27ac association (1145),(4491)
set2=subsetByOverlaps(BT245_gainedpromoters,BT245_K27ac)
set20=subsetByOverlaps(DIPG13_gainedpromoters,DIPG13_K27ac)
#what are these genes
peakAnno.test <- annotatePeak(set20, tssRegion=c(-2500, 2500),
                              TxDb=txdb, annoDb="org.Hs.eg.db")
set20.annotated=as.data.frame(peakAnno.test)
set20.annotated=set20.annotated[-c(5:15,17)]
names(set20.annotated)=paste(c("Chromosome","Start","End","Width","Probe"))
library(dplyr)
DEG=read.table("DIPG13_DEG.txt",header=T)
set20.annotated=left_join(set20.annotated,DEG,by="Probe")
#DEG count in set 1 and 2
set20.annotated=set20.annotated[complete.cases(set20.annotated),]
a1=set10.annotated[set10.annotated$P.adj<0.1,]
a2=a1[a1$Log2_FC>1,]
a3=a1[a1$Log2_FC<c(-1),]

b1=set20.annotated[set20.annotated$P.adj<0.1,]
b2=b1[b1$Log2_FC>1,]
b3=b1[b1$Log2_FC<c(-1),]
#60are up and 13 are down out of 240 are DEG in K27me3 gained - BT245
#32are up are 43 are down  and out of 896 are DEG with K27ac -BT245
#activating vs repressive PRC1
#what are these genes
set3=DIPG13_gainedpromoters
peakAnno.test <- annotatePeak(set3, tssRegion=c(-2500, 2500),
                              TxDb=txdb, annoDb="org.Hs.eg.db")
set3.annotated=as.data.frame(peakAnno.test)
set3.annotated=set3.annotated[-c(5:15,17)]
names(set3.annotated)=paste(c("Chromosome","Start","End","Width","Probe"))
library(dplyr)
DEG=read.table("BT245_DEG.txt",header=T)
set3.annotated=left_join(set3.annotated,DEG,by="Probe")
set3.annotated=set3.annotated[complete.cases(set3.annotated),]

c1=set3.annotated[set3.annotated$P.adj<0.1,]
c2=c1[c1$Log2_FC>0,]
c3=c1[c1$Log2_FC<c(0),]

library(GO)
library(ggplot2)
library(ReactomePA)
test=makeGRangesFromDataFrame(c1,seqnames.field = "Chromosome.x",start.field = "Start.x",end.field = "End.x")
export.bed(test,"test.bed")
#peakanno 
peakAnno <- annotatePeak(test, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")
pathway1 <- enrichKEGG(as.data.frame(peakAnno)$geneId)
pathway1 <- enrichGO(as.data.frame(peakAnno)$geneId,"org.Hs.eg.db")
head(pathway1)
dotplot(pathway1)

#gained with loops
BT245_gainedpromoters
BT245_loops=read.table("BT245_loops.txt")
BT245_loops=makeGRangesFromDataFrame(BT245_loops,start.field = "V2",end.field = "V3",seqnames.field = "V1")
BT245_c2_loops=read.table("BT245_c2_loops.txt")
BT245_c2_loops=makeGRangesFromDataFrame(BT245_c2_loops,start.field = "V2",end.field = "V3",seqnames.field = "V1")

BT245_gainedloops=setdiff(BT245_loops,subsetByOverlaps(BT245_loops,BT245_c2_loops))

set4=subsetByOverlaps(BT245_gainedpromoters,BT245_gainedloops)
peakAnno.test <- annotatePeak(set4, tssRegion=c(-2500, 2500),
                              TxDb=txdb, annoDb="org.Hs.eg.db")
set4.annotated=as.data.frame(peakAnno.test)
set4.annotated=set4.annotated[-c(5:15,17)]
names(set4.annotated)=paste(c("Chromosome","Start","End","Width","Probe"))

set4.annotated=left_join(set4.annotated,DEG,by="Probe")
set4.annotated=set4.annotated[complete.cases(set4.annotated),]

d1=set4.annotated[set4.annotated$P.adj<0.1,]
d2=d1[d1$Log2_FC>0,]
d3=d1[d1$Log2_FC<c(0),]
