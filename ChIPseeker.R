#use ChIPseeker to annotate peaks
## loading packages
library(ChIPseeker)
library(rtracklayer)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <-TxDb.Hsapiens.UCSC.hg19.knownGene
cgi=read.table("cpgislands_hg19.bed")
cgi=makeGRangesFromDataFrame(cgi,start.field = "V2",end.field = "V3",seqnames.field = "V1")
library(clusterProfiler)

#load files
setwd("C:/Users/egj19/OneDrive - McGill University/Desktop/Analysis/peakcalls")

library(readxl)
chipseeker <- read_excel("chipseeker.xlsx")
BT245=import.bed("BT245_RING1b_peaks_optimized.bed")
DIPG13=import.bed("DIPG13_RING1b_peaks_optimized.bed")
BT245.ko=import.bed("BT245_C2P8_RING1b_peaks_optimized.bed")
DIPG13.ko=import.bed("DIPG13_C5P6_RING1b_peaks_optimized.bed")
setwd("C:/Users/egj19/OneDrive - McGill University/Desktop/Analysis/BEDs")
BT245_K27ac=import.bed("BT245_K27ac.bed")
BT245_K27MKO_K27ac=import.bed("BT245_K27MKO_K27ac.bed")
DIPG13_K27ac=import.bed("DIPG13_K27ac.bed")
DIPG13_K27MKO_K27ac=import.bed("DIPG13_K27MKO_K27ac.bed")
BT245_CBX2=import.bed("BT245_CBX2_peaks_seqmonk.bed")
DIPG13_CBX2=import.bed("DIPG13_CBX2_peaks_seqmonk.bed")
BT245_H3K27me3_10kb=import.bed("BT245_H3K27me3_10kb.bed")
BT245_K27MKO_H3K27me3_10kb=import.bed("BT245_K27MKO_H3K27me3_10kb.bed")
DIPG13_H3K27me3_10kb=import.bed("DIPG13_H3K27me3_10kb.bed")
DIPG13_K27MKO_H3K27me3_10kb=import.bed("DIPG13_K27MKO_H3K27me3_10kb.bed")
promoter_cgi=import.bed("promoter_cgi.bed")
promoter_nocgi=import.bed("promoter_nocgi.bed")
K27M=suppressWarnings(subsetByOverlaps(BT245,DIPG13))
KO=suppressWarnings(subsetByOverlaps(DIPG13.ko,BT245.ko))
setwd("C:/Users/egj19/OneDrive - McGill University/Desktop/Analysis/BEDs")
bt245_ko_ncprc1=import.bed("bt245_ko_ncprc1.bed")
bt245_ko_cprc1=import.bed("bt245_ko_cprc1.bed")
bt245_ncprc1=import.bed("bt245_ncprc1.bed")
bt245_cprc1=import.bed("bt245_cprc1.bed")
dipg13_ncprc1=import.bed("dipg13_ncprc1.bed")
dipg13_cprc1=import.bed("dipg13_cprc1.bed")
dipg13_ko_ncprc1=import.bed("dipg13_ko_ncprc1.bed")
dipg13_ko_cprc1=import.bed("dipg13_ko_cprc1.bed")
setwd("C:/Users/egj19/OneDrive - McGill University/Desktop/Analysis/Peakcalls")
cprc1=import.bed("cprc1.bed")
ncprc1=import.bed("ncprc1.bed")
#cell line differential peaks
BT245_cprc1_gained=setdiff(bt245_cprc1,subsetByOverlaps(bt245_cprc1,bt245_ko_cprc1))
DIGP13_cprc1_gained=setdiff(dipg13_cprc1,subsetByOverlaps(dipg13_cprc1,dipg13_ko_cprc1))
K27M_cprc1_gained=subsetByOverlaps(BT245_cprc1_gained,DIGP13_cprc1_gained)

BT245_superrecruitment=subsetByOverlaps(bt245_cprc1,bt245_ko_cprc1)

BT245_ncprc1_gained=setdiff(bt245_ncprc1,subsetByOverlaps(bt245_ncprc1,bt245_ko_ncprc1))
DIGP13_ncprc1_gained=setdiff(dipg13_ncprc1,subsetByOverlaps(dipg13_ncprc1,dipg13_ko_ncprc1))
K27M_ncprc1_gained=subsetByOverlaps(BT245_ncprc1_gained,DIGP13_ncprc1_gained)

#establish ring1b differential peaks -updated
K27M.lost=setdiff(KO,subsetByOverlaps(KO,K27M))
K27M.maintained=subsetByOverlaps(K27M,KO)
K27M.gained=setdiff(K27M,subsetByOverlaps(K27M,KO))
length(K27M)-length(K27M.gained)-length(K27M.maintained)
length(KO)-length(K27M.lost)-length(subsetByOverlaps(KO,K27M))
names(chipseeker)

#heatmap on TSS
promoter <- getPromoters(TxDb=txdb, upstream=2500, downstream=2500)
tagMatrix1 <- getTagMatrix(K27M, windows=promoter_nocgi)
tagMatrix2 <- getTagMatrix(K27M, windows=promoter_cgi)
tagMatrix3 <- getTagMatrix(K27M, windows=promoter)

#aggregate plot on TSS
plotAvgProf(tagMatrix1, xlim=c(-5000, 5000),
            xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
plotAvgProf(tagMatrix2, xlim=c(-5000, 5000),
            xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
plotAvgProf(tagMatrix3, xlim=c(-5000, 5000),
            xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
export.bed(promoter,"promoter.bed")


#95% CI
plotAvgProf(tagMatrix, xlim=c(-5000, 5000), conf = 0.95, resample = 1000)

#annotate peaks
par(mfrow=c(2,2))
peakAnno.b <- annotatePeak(BT245.1, tssRegion=c(-2500, 2500),
                                TxDb=txdb, annoDb="org.Hs.eg.db")
plotAnnoBar(peakAnno.b)

peakAnno.d <- annotatePeak(DIPG13.1, tssRegion=c(-3000, 3000),
                                TxDb=txdb, annoDb="org.Hs.eg.db")
plotAnnoBar(peakAnno.d)

peakAnno.b.ko <- annotatePeak(BT245.ko.1, tssRegion=c(-3000, 3000),
                                TxDb=txdb, annoDb="org.Hs.eg.db")
plotAnnoBar(peakAnno.b.ko)

peakAnno.d.ko <- annotatePeak(DIPG13.ko.1, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")
plotAnnoBar(peakAnno.d.ko)

peakAnno.k27m <- annotatePeak(K27M, tssRegion=c(-2500, 2500),
                           TxDb=txdb, annoDb="org.Hs.eg.db")
plotAnnoBar(peakAnno.k27m)

peakAnno.ko <- annotatePeak(KO, tssRegion=c(-2500, 2500),
                           TxDb=txdb, annoDb="org.Hs.eg.db")
plotAnnoBar(peakAnno.ko)
peakAnno.test <- annotatePeak(set1, tssRegion=c(-2500, 2500),
                            TxDb=txdb, annoDb="org.Hs.eg.db")
#jpeg("AnnoBar_BT245_RING1b_gainedpeaks_AssociatedWithK27me3.jpeg",res=600,width = 10, height = 5, units = 'in')
plotAnnoBar(peakAnno.test)
#dev.off()
upsetplot(peakAnno.test, vennpie=TRUE)
#annotation relative to TSS
#jpeg("AnnoBar_CGIs_allpeaks_TSS.jpeg",res=600,width = 10, height = 5, units = 'in')

plotDistToTSS(peakAnno.test,
              title="Distribution of all CGIs loci\n relative to TSS")
#dev.off()

export.bed(BT245,"test.bed")

b=a$data
b=b[order(b$sign),]

#histogram plots distance to TSS
df <- data.frame(
  Distance=c("<-100Kb","-100-10Kb","-10-5Kb","-5-3Kb","-3-1Kb","-1-0Kb","0-1Kb","1-3Kb","3-5Kb","5-10Kb","10-100Kb",">100Kb"),
  Percentage=c(4.969291,28.894472,7.984366,3.350084,2.707984,5.248465,5.388051,2.205472,2.707983,4.997208,25.321050,6.225572))
)
library(ggpubr)
ggplot(df, aes(x=Distance,y=Percentage)) +
  geom_bar(stat="identity",width=0.7, fill="olivedrab3",color="black")+
  scale_x_discrete(limits=df$Distance)+ xlab("Distance from TSS")+ ylab("Percentage of Peaks")+
  ggtitle("RING1b superrecruitment to H3K27me3 in BT245")+
  theme_minimal()+ylim(0,30)
#ggsave("BarPlots_RING1b_superrecruitment_distribution.tiff", units="in", width=8, height=4, dpi=600, compression = 'lzw')

#promoters without cgi
promoter <- getPromoters(TxDb=txdb, upstream=2500, downstream=2500)
cgi=read.table("cpgislands_hg19.bed")
cgi=makeGRangesFromDataFrame(cgi,start.field = "V2",end.field = "V3",seqnames.field = "V1")
promoter_cgi=subsetByOverlaps(promoter,subsetByOverlaps(promoter,cgi))
promoter_nocgi=subsetByOverlaps(promoter,setdiff(promoter,promoter_cgi))
length(promoter)-length(promoter_cgi)-length(promoter_nocgi)
#export.bed(promoter_nocgi,"promoter_nocgi.bed")
#out of 22801 promoters, 34.535% had no CGIs, 65.46% had CGIs
#out of 14927 promoters with cgis, RING1b occupies:13102 of them,12534 of them
#out of 7874 promoters without cgis, RING1b occupies 1887 of them, 2006 of them
#out of 15474 promters with cgis, K27me3 occupies:5747 of them, 7088 of them
#out of 15474 promoters with cgis, K27ac occupies 9704 of them
#out of 7874 promoters without cgis, K27ac occupies 731 of them
#out of 22801 promoters,RING1b occupies 14934 of them, 15204 of them.
#out of promoters with RING1b,3209 had K27me3, 9997 had K27ac
#out of the 15666 promoter-specific ring1b peaks, at least 1636 were canonical
setwd("C:/Users/egj19/OneDrive - McGill University/Desktop/Analysis/peakcalls")

#13102 RING1b at promoters having CGIs:4675 were cPRC1, 8247 were ncPRC1
#out of 8247 ncPRC1,7178 had K27ac
a=subsetByOverlaps(promoter,subsetByOverlaps(promoter,pcGBM2))
c=subsetByOverlaps(a,subsetByOverlaps(a,pcGBM2_K27ac))
results1=length(c)/length(a)
d=subsetByOverlaps(a,setdiff(a,c))
e=subsetByOverlaps(d,pcGBM2_K27me3)
f=subsetByOverlaps(e,BT245_CBX2)
#export.bed(c,"pcGBM2_ncPRC1.bed")

a1=subsetByOverlaps(promoter_cgi,subsetByOverlaps(promoter_cgi,BT245.ko))
c1=subsetByOverlaps(a1,subsetByOverlaps(a1,BT245_K27MKO_K27ac))
results2=length(c1)/length(a1)
d1=subsetByOverlaps(a1,setdiff(a1,c1))
e1=subsetByOverlaps(d1,BT245_K27MKO_H3K27me3_10kb)

length(subsetByOverlaps(e,e1))/length(e)
length(subsetByOverlaps(c,c1))/length(c)

View(as.data.frame(setdiff(e,e1)))

a2=subsetByOverlaps(promoter_cgi,subsetByOverlaps(promoter_cgi,DIPG13))
c2=subsetByOverlaps(a2,subsetByOverlaps(a2,DIPG13_K27ac))
results3=length(c2)/length(a2)
d2=subsetByOverlaps(a2,setdiff(a2,c2))

a3=subsetByOverlaps(promoter_cgi,subsetByOverlaps(promoter_cgi,DIPG13.ko))
c3=subsetByOverlaps(a3,subsetByOverlaps(a3,DIPG13_K27MKO_K27ac))
results4=length(c3)/length(a3)
d3=subsetByOverlaps(a3,setdiff(a3,c3))


#Promoter RING1b annotation
# library
library(ggplot2)
library(GenomicRanges)
#dataframe
a1=length(subsetByOverlaps(promoter,subsetByOverlaps(promoter_cgi,BT245.ko)))
a2=length(promoter_cgi)-a1
a3=length(promoter_nocgi)
a4=length(subsetByOverlaps(promoter,subsetByOverlaps(promoter,BT245.ko)))
a5=length(promoter)-a4
a6=length(subsetByOverlaps(promoter,subsetByOverlaps(promoter_nocgi,BT245.ko)))
a7=length(promoter_nocgi)-a6
e1=subsetByOverlaps(promoter_cgi,subsetByOverlaps(promoter_cgi,BT245.ko))
e2=suppressWarnings(subsetByOverlaps(e1,BT245_K27MKO_K27ac))
e3=setdiff(e1,subsetByOverlaps(e1,e2))
e4=subsetByOverlaps(e3,BT245_K27MKO_H3K27me3_10kb)
e5=subsetByOverlaps(e4,BT245_KO_CBX2)
e6=subsetByOverlaps(e2,BT245_KO_CBX2)

res1=length(e1)
res2=length(e2)
res3=length(e3)
res4=length(e4)
res5=length(e5)
res6=length(e6)

setwd("C:/Users/egj19/OneDrive - McGill University/Desktop/Analysis/BEDs")
bt245_ko_ncprc1=import.bed("bt245_ko_ncprc1.bed")
bt245_ko_cprc1=import.bed("bt245_ko_cprc1.bed")
bt245_ncprc1=import.bed("bt245_ncprc1.bed")
bt245_cprc1=import.bed("bt245_cprc1.bed")
dipg13_ncprc1=import.bed("dipg13_ncprc1.bed")
dipg13_cprc1=import.bed("dipg13_cprc1.bed")
dipg13_ko_ncprc1=import.bed("dipg13_ko_ncprc1.bed")
dipg13_ko_cprc1=import.bed("dipg13_ko_cprc1.bed")

#make differential sites:
BT245_gained_cprc1=setdiff(bt245_cprc1,subsetByOverlaps(bt245_cprc1,bt245_ko_cprc1))
BT245_maintained_cprc1=subsetByOverlaps(bt245_ko_cprc1,bt245_cprc1)
BT245_lost_cprc1=setdiff(bt245_ko_cprc1,subsetByOverlaps(bt245_ko_cprc1,bt245_cprc1))

DIPG13_gained_cprc1=setdiff(dipg13_cprc1,subsetByOverlaps(dipg13_cprc1,dipg13_ko_cprc1))
DIPG13_maintained_cprc1=subsetByOverlaps(dipg13_ko_cprc1,dipg13_cprc1)
DIPG13_lost_cprc1=setdiff(dipg13_ko_cprc1,subsetByOverlaps(dipg13_ko_cprc1,dipg13_cprc1))


k27m_ncprc1=intersect(bt245_ncprc1,dipg13_ncprc1)
k27m_cprc1=intersect(bt245_cprc1,dipg13_cprc1)
ko_ncprc1=intersect(bt245_ko_ncprc1,dipg13_ko_ncprc1)
ko_cprc1=intersect(bt245_ko_cprc1,dipg13_ko_cprc1)

View(as.data.frame(setdiff(k27m_cprc1,ko_cprc1)))

length(subsetByOverlaps(bt245_ncprc1,bt245_ko_ncprc1))/length(bt245_ncprc1)
length(subsetByOverlaps(dipg13_ncprc1,dipg13_ko_ncprc1))/length(dipg13_ncprc1)
length(subsetByOverlaps(bt245_cprc1,bt245_ko_cprc1))/length(bt245_cprc1)
length(subsetByOverlaps(dipg13_cprc1,dipg13_ko_cprc1))/length(dipg13_cprc1)
   
df2 <- data.frame(Legend=rep(c("RING1b", "No RING1b"), each=3),
                  dose=rep(c("All Promoters","Promoters with CGIs", "Promoters with no CGIs"),2),
                  len=c(a4,a1,a6,a5,a2,a7))

library(ggthemes)
library(hrbrthemes)
library("ggsci")
library("ggplot2")
library("gridExtra")
#plot
ggplot(data=df2, aes(x=dose, y=len, fill=Legend)) +
  geom_bar(stat="identity", position=position_dodge())+
  geom_text(aes(label=len), vjust=1.6, color="white",
            position = position_dodge(0.9), size=3.5)+
    scale_fill_npg()+
  xlab("")+ylab("Number of Promoters")+
  theme_minimal()+ggtitle("BT245 Parental")
ggsave("BarPlots_Promoters_RING1b_HSJ51.tiff", units="in", width=7, height=4, dpi=600, compression = 'lzw')



#cPRC1 vs ncPRC1
df <- data.frame(Legend=c("H3K27me3-enriched","H3K27me3-depleted"),CBX2=c(res5,res6),RING1b=c(res4,res2),percentage=c(res5/res4*100,res6/res2*100))
 df$percentage2=paste(round(df$percentage),"%")                
#plot
library(reshape2)
#copy
ggplot(df, aes(x=Legend, fill=Legend))+
  geom_bar(aes(y=RING1b),stat="identity", color="black")+ggtitle("DIPGXIII K27M KO - CGI Promoters")+ylab("Number of Promoters")+
  theme_minimal()+scale_fill_excel_new()+geom_text(aes(y=RING1b,label=df$RING1b),vjust=1.5, color="white", size=4.5)
ggsave("BarPlots_PromotersCGIs_RING1b_DIPG13_K27MKO_canonicalvsncanonical.tiff", units="in", width=7, height=4, dpi=600, compression = 'lzw')

#optional: superimpose cbx2:   
ggplot(df, aes(x=Legend, fill=Legend))+
  geom_bar(aes(y=RING1b),stat="identity", color="black",alpha=0.5)+
  geom_bar(aes(y=CBX2),stat="identity", color="black",position ="identity",alpha=1)+
  xlab("")+ylab("Percentage of Promoters with CpGs")+ggtitle("CBX2 Occupancy at CGI Promoters - BT245 K27M KO")+
  theme_minimal()+scale_fill_excel_new()+geom_text(aes(y=CBX2,label=df$percentage2),vjust=-2, color="black", size=5.5)
ggsave("BarPlots_PromotersCGIs_RING1b_BT245_KO_canonicalvsncanonical_cbx2.tiff", units="in", width=7, height=4, dpi=600, compression = 'lzw')

#motif analysis on gained RING1b
export.bed(BT245_gainedpeaks,"test.bed")
#heatmap using genomation
library("genomation")
mat=ScoreMatrix("BT245-DMSO-2-Rx_XChIP_RING1B_condition3.bw", windows=promoter, strand.aware = FALSE, weight.col = NULL,
            is.noCovNA = FALSE, type = "auto", rpm = FALSE, unique = FALSE,
            extend = 0, param = NULL, bam.paired.end = FALSE, library.size = NULL)

heatMatrix(mat, grid = FALSE, col = NULL, xcoords = NULL, group = NULL,
           group.col = NULL, order = FALSE, user.order = FALSE, winsorize = c(0,
                                                                              100), clustfun = NULL, main = "", legend.name = NULL, cex.legend = 1,
           xlab = NULL, cex.main = 1, cex.lab = 1, cex.axis = 1,
           newpage = TRUE)


#annotate peaks
peakAnnoList <- lapply(chipseeker, annotatePeak, TxDb=txdb, 
                       tssRegion=c(-2500, 2500), verbose=FALSE)
#gene annotation
data(TSS.human.GRCh37) 
overlaps.anno <- annotatePeakInBatch(overlaps, AnnotationData=annoData, 
                                   output="overlapping", maxgap=5000L)
overlaps.anno <- addGeneIDs(overlaps.anno, "org.Hs.eg.db", "symbol")

peakAnnoList$DIPG13_RING1b
peakAnnoList$DIPG13_KO_RING1b

pro.b.k27m=sum(31,1.915,1.1)
pro.d.k27m=sum(42.8984310,3.0435296,1.3861036)
pro.b.ko=sum(29.833,2.08,1.157)
pro.d.ko=sum(25.0356108,3.6003733,1.5766983)
inter.b.k27m=sum(33.355,1.012)
inter.d.k27m=sum(28.8545476,0.6016279)
inter.b.ko=sum(36.51,0.869)
inter.d.ko=sum(40.3261457,0.9037772)
intra.b.k27m=100-(pro.b.k27m+inter.b.k27m)
intra.d.k27m=100-(pro.d.k27m+inter.d.k27m)
intra.b.ko=100-(pro.b.ko+inter.b.ko)
intra.d.ko=100-(pro.d.ko+inter.d.ko)
#create dataset
df2 <- data.frame(Condition=rep(c("BT245", "BT245 K27M KO","DIPGXIII","DIPGXIII K27M KO"), each=3),
                  Annotation=rep(c("Promoters (within ±2.5Kb from TSS)", "Intergenic (outside ±2.5Kb from TSS)", "Genic"),4),
                  Percentage=c(pro.b.k27m,inter.b.k27m,intra.b.k27m,pro.b.ko,inter.b.ko,intra.b.ko,pro.d.k27m,inter.d.k27m,intra.d.k27m,pro.d.ko,inter.d.ko,intra.d.ko))
# Change the colors manually
library(ggplot2)
p <- ggplot(data=df2, aes(x=Annotation, y=Percentage, fill=Condition)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_minimal()
p + scale_fill_brewer(palette="YlGnBu")

#plot2
library(viridis)
library(hrbrthemes)
ggplot(df2, aes(fill=Annotation, y=Percentage, x=Condition)) + 
  geom_bar(position="stack", stat="identity") +
  ylab("RING1b peaks, %") +
  scale_fill_brewer(palette="YlGnBu") +
  ggtitle("RING1b Distribution") +
  theme_minimal() +
  xlab("Cell Lines")+
  geom_text(aes(label = round(Percentage,0)), size = 3, hjust = 0.5, vjust = 3, position = "stack") 

ggsave("BarPlots_RING1b_Annotation_2.tiff", units="in", width=7, height=4, dpi=600, compression = 'lzw')

#Annotate K27M, KO, Gained, Lost and Maintained peaks
peakAnnoList$K27M_RING1b
peakAnnoList$KO_RING1b
peakAnnoList$K27M_gainedpeaks
peakAnnoList$K27M_lostpeaks
peakAnnoList$K27M_maintainedpeaks

pro.k27m=sum(59.792,1.481,0.771)
pro.ko=sum(55.4,1.9,0.939)
pro.gained=sum(53.7683824,1.7003676,0.8272059)
pro.lost=sum(42.7764326,2.5020178,0.9281679)
pro.main=sum(62.6715093,3.6521388,1.7352704)
inter.k27m=sum(20.1,0.66)
inter.ko=sum(21.939,0.586)
inter.gained=sum(22.4839154,0.7238051)
inter.lost=sum(31.8401937,0.8474576)
inter.main=sum(18.7247780,0.7263923)
intra.main=100-(pro.main+inter.main)
intra.k27m=100-(pro.k27m+inter.k27m)
intra.ko=100-(pro.ko+inter.ko)
intra.gained=100-(pro.gained+inter.gained)
intra.lost=100-(pro.lost+inter.lost)

df3 <- data.frame(Condition=rep(c("K27M", "K27M KO","Gained Peaks","Lost Peaks"), each=3),
                  Annotation=rep(c("Promoters", "Intergenic", "Intragenic"),4),
                  Percentage=c(pro.k27m,inter.k27m,intra.k27m,pro.ko,inter.ko,intra.ko,pro.gained,inter.gained,intra.gained,pro.lost,inter.lost,intra.lost))

# Convert the cyl variable to a factor
df3$Percentage <- as.numeric(df3$Percentage)
df3$i=as.numeric(df3$Annotation)


library("ggpubr")

ggpubr::ggdotchart(df3, x = "Annotation", y = "Percentage",
           color = "Condition",                                # Color by groups
           palette = "Reds", # Custom color palette
           sorting = "descending",                       # Sort value in descending order
           add = "segments",                             # Add segments from y = 0 to dots
           rotate = TRUE,                                # Rotate vertically
           group = "Condition",                                # Order by groups
           dot.size = 6,                                 # Large dot size
           label = round(df3$Percentage),                        # Add mpg values as dot labels
           font.label = list(color = "white", size = 9, 
                             vjust = 0.5),               # Adjust label parameters
           ggtheme = theme_pubr()                        # ggplot2 theme
)
#create vertical plots
library("ggthemes")
p=ggplot(df3, aes(fill=Annotation, y=Percentage, x=Condition)) + 
  geom_bar(position="stack", stat="identity") +
  ggtitle("RING1b Annotation") +
  coord_flip() +
  xlab("")+
  scale_fill_brewer(palette="YlOrRd")

p + theme_hc()+ scale_colour_hc()


#barplots for K27M, KO and Maintained Peaks annotation
df4 <- data.frame(Condition=rep(c("K27M", "K27M KO","Maintained Peaks"), each=3),
                  Annotation=rep(c("Promoters", "Intergenic", "Genic"),3),
                  Percentage=c(pro.gained,inter.gained,intra.gained,pro.lost,inter.lost,intra.lost,pro.main,intra.main,inter.main))

library("ggthemes")
ggplot(df4, aes(fill=Annotation, y=Percentage, x=Condition)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_brewer(palette="YlOrRd") +
  ggtitle("RING1b Peaks Annotation") +
  theme_minimal() +
  xlab("")
ggsave("BarPlots_RING1b_Annotation_3.tiff", units="in", width=6, height=4, dpi=600, compression = 'lzw')


#Annotation of RING1b- K27ac peaks
peakAnnoList$dipg13_k27me3_clusteronpeaks
peakAnnoList$dipg13_KO_k27me3_clusteronpeaks

b.k27ac.1=sum(36.4688204,1.7831205,1.1169547)
b.k27ac.2=sum(28.3295768)
b.k27ac.3=100-(b.k27ac.1+b.k27ac.2)
b.k27ac.4=sum(52.5160600,2.0788722,0.9635974)
b.k27ac.5=(19.5485368)
b.k27ac.6=100-(b.k27ac.4+b.k27ac.5)

d.k27ac.1=sum(45.3249929,2.1352648,0.9359421)
d.k27ac.2=sum(24.6919387)
d.k27ac.3=100-(d.k27ac.1+d.k27ac.2)
d.k27ac.4=sum(40.5911552,2.3691336,1.3312274)
d.k27ac.5=(26.9291516)
d.k27ac.6=100-(d.k27ac.4+d.k27ac.5)

b.k27me3.1=sum(25.4201681,3.0112045,1.6106443)
b.k27me3.2=sum(45.5999066)
b.k27me3.3=100-(b.k27me3.1+b.k27me3.2)
b.k27me3.4=sum(21.2056304,2.8763770,1.4891881)
b.k27me3.5=(47.4500204)
b.k27me3.6=100-(b.k27me3.4+b.k27me3.5)

d.k27me3.1=sum(23.1850437,3.4252183,1.5147380)
d.k27me3.2=sum(46.5884279)
d.k27me3.3=100-(d.k27me3.1+d.k27me3.2)
d.k27me3.4=sum(21.5857929,3.8769385,1.8092380)
d.k27me3.5=(46.0897115)
d.k27me3.6=100-(d.k27me3.4+d.k27me3.5)


#create donut chart
count.data <- data.frame(
  Annotation = c("Promoter", "Intergenic", "Genic"),
  n = c(b.k27me3.1, b.k27me3.2, b.k27me3.3),
  prop = c(d.k27ac.1, d.k27ac.2, d.k27ac.3)
)
count.data
# Add label position
count.data <- count.data %>%
  arrange(desc(Annotation)) %>%
  mutate(lab.ypos = cumsum(prop) - 0.5*prop)
count.data

mycols <- c("#0073C2FF", "#EFC000FF", "#CD534CFF")

ggplot(count.data, aes(x = 2, y = prop, fill = Annotation)) +
  geom_bar(stat = "identity", color = "white") +
  coord_polar(theta = "y", start = 0)+
  geom_text(aes(y = lab.ypos, label = round(prop,0)), color = "white")+
  scale_fill_manual(values = mycols) +
  theme_void()+
  xlim(0.5, 2.5) + 
  ggtitle("RING1b sites Associated with H3K27ac - DIPGXIII Parental")
ggsave("DonutChart_DIPG13_RING1b_H3K27ac.tiff", units="in", width=8, height=7, dpi=600, compression = 'lzw')


#chipseeker pathway analysis
library(GO)
library(ggplot2)
library(ReactomePA)
#peakanno 
peakAnno <- annotatePeak(test, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")
pathway1 <- enrichPathway(as.data.frame(peakAnno)$geneId)
head(pathway1)
#jpeg("GO_BT245_RING1b_superrecruitment.jpeg",res=600,width = 10, height = 8, units = 'in')
dotplot(pathway1)

#dev.off()

#check RINg1b peaks with gene expression K27ac
library(readxl)
setwd("C:/Users/egj19/OneDrive - McGill University/Desktop/Analysis/peakcalls")
BT245_ncPRC1_promoters <- read_excel("BT245_ncPRC1_promoters_annotation.xlsx") 
setwd("C:/Users/egj19/OneDrive - McGill University/Desktop/Analysis/BEDs")
BT245_rna=read.table("BT245_RPKM.txt",header=T)
merge1=BT245_ncPRC1_promoters
names(merge1)=paste(c("Chromosome","Start","End","Probe"))
BT245_ncPRC1_promoters=left_join(merge1,BT245_rna,by="Probe")
BT245_ncPRC1_promoters=BT245_ncPRC1_promoters[complete.cases(BT245_ncPRC1_promoters),]
names(BT245_ncPRC1_promoters)

BT245_ncPRC1_promoters$Parental <- rowMeans(BT245_ncPRC1_promoters[,c(13:17)], na.rm=TRUE)
BT245_ncPRC1_promoters$KO <- rowMeans(BT245_ncPRC1_promoters[,c(8:12)], na.rm=TRUE)
test=BT245_cPRC1_promoters
test=BT245_cPRC1_promoters[BT245_cPRC1_promoters$Parental>1,]

#out of 8147 genes with RING1b and K27ac at promoters, 8011 were expressed with RPKM>1
#out of 2529 genes with RING1b and K27me3 at promoters,1119 were expressed with RPKM>1
#check RING1b peaks with gene expression K27me3
#out of 11386 RING1b peaks with K27ac, 11047 were associated with expressed genes with RPKM>1
#out of 11076 RING1b peaks with H3K27me3, 7138 were associated with expressed genes with RPKM>1
#out of 6908 ncPRC1 promoters,6674 genes were expressed with RPKM>1
#out of 3820 cPRC1 promoters,2964 genes were expressed with RPKM>1

cluster5=import.bed("cluster5.bed")
length(subsetByOverlaps(cluster5,cgi))/length(cluster5)


#how many cPRC1 promoters had CBX2 peaks
CBX2_K27M=subsetByOverlaps(BT245_CBX2,DIPG13_CBX2)
length(subsetByOverlaps(cprc1,CBX2_K27M))/length(cprc1)
length(subsetByOverlaps(cprc1,CBX2_K27M))/length(cprc1)
length(subsetByOverlaps(ncprc1,CBX2_K27M))/length(ncprc1)
length(subsetByOverlaps(ncprc1,CBX2_K27M))/length(ncprc1)

#check for switches
a=subsetByOverlaps(cprc1,bt245_ko_cprc1)
b=subsetByOverlaps(dipg13_ncprc1,dipg13_ko_cprc1)
c=subsetByOverlaps(a,b)

a=subsetByOverlaps(bt245_cprc1,bt245_ko_ncprc1)
b=subsetByOverlaps(dipg13_cprc1,dipg13_ko_ncprc1)
c=subsetByOverlaps(a,b)
