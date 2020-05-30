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
cPRC1=import.bed("cprc1.bed")
#gained peaks
gainedpeaks.bt245=setdiff(BT245,subsetByOverlaps(BT245,BT245.ko))
gainedpeaks.dipg13=setdiff(DIPG13,subsetByOverlaps(DIPG13,DIPG13.ko))
gainedpeaks.k27m=subsetByOverlaps(gainedpeaks.bt245,gainedpeaks.dipg13)
length(subsetByOverlaps(gainedpeaks.k27m,cgi))/length(gainedpeaks.k27m)

lostpeaks.bt245=setdiff(BT245.ko,subsetByOverlaps(BT245.ko,BT245))
lostpeaks.dipg13=setdiff(DIPG13.ko,subsetByOverlaps(DIPG13.ko,DIPG13))
lostpeaks.k27m=subsetByOverlaps(lostpeaks.bt245,lostpeaks.dipg13)
length(subsetByOverlaps(lostpeaks.k27m,cgi))/length(lostpeaks.k27m)
head(lostpeaks.k27m)
maintainedpeaks.bt245=subsetByOverlaps(BT245,BT245.ko)
maintainedpeaks.dipg13=subsetByOverlaps(DIPG13,DIPG13.ko)
maintainedpeaks.k27m=subsetByOverlaps(maintainedpeaks.bt245,maintainedpeaks.dipg13)
#association with different histone marks
a=suppressWarnings(length(subsetByOverlaps(lostpeaks.k27m,subsetByOverlaps(BT245_K27MKO_K27ac,DIPG13_K27MKO_K27ac))))
b=suppressWarnings(length(subsetByOverlaps(lostpeaks.k27m,subsetByOverlaps(BT245_K27MKO_K27me3,DIPG13_K27MKO_K27me3))))
c=length(lostpeaks.k27m)-(a+b)

#plot pichart of DEG numbers
library(ggplot2)
pie3 <- data.frame(Mark=c("H3K27me3","H3K27ac","Other"),value=c(a,b,c))
library(dplyr)
# Basic piechart
pie3 <- pie3 %>% 
  arrange(desc(Mark)) %>%
  mutate(prop = value/ 506) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

ggplot(pie3, aes(x="", y=prop, fill=Mark)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + ggtitle("Association of lost RING1b Peaks in K27M")+
  geom_text(aes(y = ypos, label = value), color = "white", size=5) +
  scale_fill_manual(values = c("dodgerblue2","firebrick2","gold2"))+
  theme(legend.background = element_rect(size=2, linetype="solid"))
ggsave("PieChart_lostring1b_K27M_markassociation.tiff", units="in", width=5, height=4, dpi=600, compression = 'lzw')
#cgi association
a=suppressWarnings(length(subsetByOverlaps(lostpeaks.k27m,cgi)))
b=length(lostpeaks.k27m)-a

# Basic piechart
pie3 <- data.frame(Site=c("CGI","non-CGI"),value=c(a,b))
pie3 <- pie3 %>% 
  arrange(desc(Site)) %>%
  mutate(prop = value/ 506) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

ggplot(pie3, aes(x="", y=prop, fill=Site)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + ggtitle("Location of lost RING1b Peaks in K27M")+
  geom_text(aes(y = ypos, label = value), color = "white", size=5) +
  scale_fill_manual(values = c("seagreen3","salmon2"))+
  theme(legend.background = element_rect(size=2, linetype="solid"))
ggsave("PieChart_lostring1b_K27M_location.tiff", units="in", width=5, height=4, dpi=600, compression = 'lzw')

#where are these sites?
# library
library(ggplot2)
library(viridis)
library(hrbrthemes)
library(ggpubr)
library("ggsci")
library("ggplot2")
library("gridExtra")

# create a dataset
RING1b <- c(rep("Gained in K27M" , 3) , rep("Lost in K27M" , 3),rep("Maintained in K27M",3))
Annotation <- rep(c("Promoter" , "Genic" , "Intergenic") , 3)
value <- c(28,42,30,28,32,40,70,11,19)
data <- data.frame(RING1b,Annotation,value)

# Small multiple
ggplot(data, aes(fill=Annotation, y=value, x=RING1b)) + 
  geom_bar(position="stack", stat="identity") +
  ggtitle("Annotation of differential K27M peaks") +
  theme_classic() +ylab("Percentage, %")+
  xlab("")+scale_fill_brewer(palette="YlGnBu")

#ggsave("StackedBarplot_RING1b_K27M_annotation.tiff", units="in", width=5, height=4, dpi=600, compression = 'lzw')

peakAnno.test <- annotatePeak(mtaintainedpeaks.k27m, tssRegion=c(-2500, 2500),
                              TxDb=txdb, annoDb="org.Hs.eg.db")
plotAnnoBar(peakAnno.test)

#Percentage of loops mediated by PRC1
setwd("C:/Users/egj19/OneDrive - McGill University/Desktop/Analysis/peakcalls")
BT245_loops=read.table("BT245_loops.txt")
BT245_loops=makeGRangesFromDataFrame(BT245_loops,seqnames.field = "V1",start.field = "V2",end.field = "V3")
BT245_KO_loops=read.table("BT245_c2_loops.txt")
BT245_KO_loops=makeGRangesFromDataFrame(BT245_KO_loops,seqnames.field = "V1",start.field = "V2",end.field = "V3")
gained_loops=setdiff(BT245_loops,subsetByOverlaps(BT245_loops,BT245_KO_loops))
gained_loops=subsetByOverlaps(gained_loops,BT245)
gained_loops=as.data.frame(gained_loops)
BT245_loops=as.data.frame(BT245_loops)
#gained_loops$probe=paste(gained_loops$seqnames,gained_loops$start,gained_loops$end)
BT245_loops=read.table("BT245_loops.txt")
BT245_loops$probe=paste(BT245_loops$V1,BT245_loops$V2,BT245_loops$V3)
#gained_loops=left_join(gained_loops,BT245_loops[,c(4,5)])
#gained_loops=gained_loops[complete.cases(gained_loops),]

gained_loops=makeGRangesFromDataFrame(gained_loops)
super=subsetByOverlaps(BT245,BT245.ko)
super=subsetByOverlaps(super,BT245_K27me3)
length(subsetByOverlaps(gained_loops,super))

super=subsetByOverlaps(gained_loops,super)
super=as.data.frame(super)
super$probe=paste(super$seqnames,super$start,super$end)
super=left_join(super,BT245_loops[,c(4,5)])
super=super[complete.cases(super),]
super.list=as.data.frame(super$V4)
names(super.list)=paste("V4")
super.list2=left_join(super.list,BT245_loops,by="V4")
super.list2=super.list2[!duplicated(super.list2),]
super.list2=super.list2[complete.cases(super.list2),]

super.list2=makeGRangesFromDataFrame(super.list2,seqnames.field = "V1",start.field = "V2",end.field = "V3")
super.list2=subsetByOverlaps(BT245,super.list2)
export.bed(super.list2,"superrecruitment_loops.bed")

BT245_DEG=read.table("BT245_DEG.txt",header=T)
superrecruitment_loops.anno=read_excel("superrecruitment_loops.anno.xlsx")
superrecruitment_loops.anno=left_join(superrecruitment_loops.anno,BT245_DEG,by="Probe")
superrecruitment_loops.anno=superrecruitment_loops.anno[complete.cases(superrecruitment_loops.anno),]
test=superrecruitment_loops.anno[superrecruitment_loops.anno$P.adj<0.1,]
test=makeGRangesFromDataFrame(test,seqnames.field = "Chromosome.x",start.field = "Start.x",end.field = "End.x")
export.bed(test,"test.bed")


#test2=subsetByOverlaps(BT245,bt245_cprc1)
#test2=read_excel("test.anno.xlsx")
#test2=left_join(test2,BT245_DEG,by="Probe")
#test2=test2[complete.cases(test2),]
#test2=test2[test2$P.adj<0.1,]
#test2=makeGRangesFromDataFrame(test2,seqnames.field = "Chromosome.x",start.field = "Start.x",end.field = "End.x")

length(subsetByOverlaps(test2,gained_loops))

#how many of these upregulated


#beyond promoters: H3K27me3 loss at CGIs
K27M_K27me3=subsetByOverlaps(BT245_K27me3,DIPG13_K27me3)
KO_K27me3=subsetByOverlaps(BT245_K27MKO_K27me3,DIPG13_K27MKO_K27me3)
#find lost H3K27me3 at CGI:2291 cgis lost H3K27me3
lost=setdiff(KO_K27me3,subsetByOverlaps(KO_K27me3,K27M_K27me3))
cgi_lost=subsetByOverlaps(cgi,lost)
#how many still had RING1b
lostpeaks=setdiff(KO,subsetByOverlaps(KO,K27M))
maintainedpeaks=subsetByOverlaps(K27M,KO)
gainedpeaks=setdiff(K27M,subsetByOverlaps(K27M,KO))
length(subsetByOverlaps(cgi_lost,maintainedpeaks))/length(cgi_lost)
length(subsetByOverlaps(cgi_lost,gainedpeaks))/length(cgi_lost)
length(subsetByOverlaps(cgi_lost,lostpeaks))/length(cgi_lost)
#how many cgis fall at promoters
length(subsetByOverlaps(cgi,promoter))/length(cgi)
length(subsetByOverlaps(BT245,subsetByOverlaps(cgi,promoter)))/length(subsetByOverlaps(BT245,cgi))
#gained RING1b peaks


#preference for promoters with cgis
a=c(13214/14996,7444/7987,12634/14544,4667/5260)
b=1-a
BiocManager::install("BSDA")
library(BSDA)
z.test(a,b,sigma.x = sd(a),sigma.y = sd(b))

#gained BT245 cPRC1 ring1b
bt245_cprc1=import.bed("bt245_cprc1.bed")
bt245_ko_cprc1=import.bed("bt245_ko_cprc1.bed")
bt245_ncprc1=import.bed("bt245_ncprc1.bed")
bt245_ko_ncprc1=import.bed("bt245_ko_ncprc1.bed")   
cprc1=read_excel("cprc1.annotation.xlsx")
ncprc1=read_excel("ncprc1.annotation.xlsx")       

gained_cprc1_bt245=setdiff(bt245_cprc1,subsetByOverlaps(bt245_cprc1,bt245_ko_cprc1))

cprc1=left_join(cprc1,K27M_DEG,by="Probe")
cprc1=cprc1[complete.cases(cprc1),]
ncprc1=left_join(ncprc1,K27M_DEG,by="Probe")
ncprc1=ncprc1[complete.cases(ncprc1),]
#percentage of cPRC1 that show change in loops
length(subsetByOverlaps(gained_cprc1_bt245,BT245_loops))/length(gained_cprc1_bt245)
length(subsetByOverlaps(bt245_cprc1,BT245_loops))/length(bt245_cprc1)
length(subsetByOverlaps(bt245_ncprc1,BT245_loops))/length(bt245_ncprc1)
test=cprc1[cprc1$P.adj<0.1,]
test=test[abs(test$Log2_FC)>1,]
test2=ncprc1[ncprc1$P.adj<0.1,]
test2=test2[abs(test2$Log2_FC)>1,]
test=makeGRangesFromDataFrame(test,seqnames.field = "Chromosome.x",start.field = "Start.x",end.field = "End.x")
test2=makeGRangesFromDataFrame(test2,seqnames.field = "Chromosome.x",start.field = "Start.x",end.field = "End.x")
length(subsetByOverlaps(test,gained_loops))
length(subsetByOverlaps(test2,gained_loops))
length(subsetByOverlaps(test,BT245_loops))
length(subsetByOverlaps(test2,BT245_loops))

#t-test for H2aK119ub levels
K27M=c(0.00917,0.010,0.01422,0.02517)
KO=c(0.0062,0.0089,0.00752,0.00724)
K27M=c(0.01422,0.02517)
KO=c(0.00752,0.00724)

t.test(K27M,KO,paired=F)
x=data.frame(deg=c(17,28,380),promoters=c(238,171,2535))
x=as.matrix(x)
chisq.test(x)
