#Update data
setwd("C:/Users/egj19/OneDrive - McGill University/Desktop/Analysis/peakcalls")
data2=read.table("seqmonk42.txt",header=T)
data=read.table("Parental_allmarks_chip_input.txt",header=T)
#Remove Blacklisted regions
blacklist=read.table("blacklist_hg19.bed",header=F)
blacklist=makeGRangesFromDataFrame(blacklist,seqnames.field = "V1",start.field = "V2",end.field = "V3")
data.range=data[c(1:3)]
data.range=makeGRangesFromDataFrame(data.range,start.field = "Start",end.field = "End",seqnames.field = "Seqnames")
data.range.clean=subsetByOverlaps(data.range,setdiff(data.range,subsetByOverlaps(data.range,blacklist)))
data.range.clean=as.data.frame(data.range.clean)
data.range.clean$probe=paste(data.range.clean$seqnames,data.range.clean$start,data.range.clean$end)
data$probe=paste(data$seqnames,data$start,data$end)
library(dplyr)
data=left_join(data.range.clean,data,by="probe")
names(data)


data2$probe=paste(data2$Chromosome,data2$Start,data2$End)
names(data2)
data2=data2[,-c(1:4)]
library(dplyr)
data=left_join(data,data2,by="probe")
data=data[complete.cases(data),]
data=as.data.frame(data)
names(data)
data=data[-c(83)]
write.table(data,"Parental_allmarks_chip_input.txt", sep="\t")


#KO
data=read.table("KO_allmarks_chip_input.txt",header=T)
data2=read.table("seqmonk43.txt",header=T)

blacklist=read.table("blacklist_hg19.bed",header=F)
blacklist=makeGRangesFromDataFrame(blacklist,seqnames.field = "V1",start.field = "V2",end.field = "V3")
data.range=data[c(2,3,4)]
data.range=makeGRangesFromDataFrame(data.range,seqnames.field = "Chromosome.x",start.field = "Start.x",end.field = "End.x")
data.range.clean=subsetByOverlaps(data.range,setdiff(data.range,subsetByOverlaps(data.range,blacklist)))
data.range.clean=as.data.frame(data.range.clean)
data.range.clean$probe=paste(data.range.clean$seqnames,data.range.clean$start,data.range.clean$end)
data$probe=paste(data$Chromosome,data$Start,data$End)
data=left_join(data.range.clean,data)
names(data)

data2$probe=paste(data2$Chromosome,data2$Start,data2$End)
library(dplyr)
data2=data2[,-c(1:4)]
data=left_join(data,data2,by="probe")


data=data[complete.cases(data),]
names(data)
data=data[-c(31)]
write.table(data, "KO_allmarks_chip_input.txt", sep="\t")

#matrix of rpkm from Shakour

matrixrpkm=read.table("matrixRPKM.tsv")
library('biomaRt')
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- rownames(matrixrpkm)
matrixrpkm$ensembl_gene_id=genes
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)
matrixrpkm=merge(matrixrpkm,G_list,by="ensembl_gene_id")
names(matrixrpkm)
BT245_RPKM=matrixrpkm[c(2,3,4,5,6,7,8,9,10,11)]
write.table(BT245_RPKM,"BT245_RPKM.txt",sep="\t")


#make rpkms:
setwd("C:/Users/egj19/OneDrive - McGill University/Desktop/Analysis/peakcalls")

BT245_rna=read.table("DESeq2_BT245.txt",header=T)
BT245_rna$genelength=BT245_rna$End-BT245_rna$Start
total1=sum(BT245_rna$BT245.NKO__BT245.NKO.P31.sorted.bam)
total2=sum(BT245_rna$BT245.NKO__BT245.NKO.P30.sorted.bam)
total3=sum(BT245_rna$BT245.NKO__BT245.NKO.P18.sorted.bam)
total4=sum(BT245_rna$BT245.NKO__BT245.NKO.P17.sorted.bam)
total5=sum(BT245_rna$BT245.NKO__BT245.NKO.P16.sorted.bam)
total6=sum(BT245_rna$BT245.KO.LATE__BT245.KO.C4P6.sorted.bam)
total7=sum(BT245_rna$BT245.KO.LATE__BT245.KO.C2P6.sorted.bam)
total8=sum(BT245_rna$BT245.KO.EARLY__BT245.KO.C5P4.sorted.bam)
total9=sum(BT245_rna$BT245.KO.EARLY__BT245.KO.C4.sorted.bam)
total10=sum(BT245_rna$BT245.KO.EARLY__BT245.KO.C2.sorted.bam)


fac1=total1/1000000
fac2=total2/1000000
fac3=total3/1000000
fac4=total4/1000000
fac5=total5/1000000
fac6=total6/1000000
fac7=total7/1000000
fac8=total8/1000000
fac9=total9/1000000
fac10=total10/1000000

BT245_rna$BT245_Parental1=BT245_rna$BT245.NKO__BT245.NKO.P31.sorted.bam/(fac1*BT245_rna$genelength)
BT245_rna$BT245_Parental2=BT245_rna$BT245.NKO__BT245.NKO.P30.sorted.bam/(fac2*BT245_rna$genelength)
BT245_rna$BT245_Parental3=BT245_rna$BT245.NKO__BT245.NKO.P18.sorted.bam/(fac3*BT245_rna$genelength)
BT245_rna$BT245_Parental4=BT245_rna$BT245.NKO__BT245.NKO.P17.sorted.bam/(fac4*BT245_rna$genelength)
BT245_rna$BT245_Parental5=BT245_rna$BT245.NKO__BT245.NKO.P16.sorted.bam/(fac5*BT245_rna$genelength)
BT245_rna$BT245_K27MKO_1=BT245_rna$BT245.KO.LATE__BT245.KO.C4P6.sorted.bam/(fac6*BT245_rna$genelength)
BT245_rna$BT245_K27MKO_2=BT245_rna$BT245.KO.LATE__BT245.KO.C2P6.sorted.bam/(fac7*BT245_rna$genelength)
BT245_rna$BT245_K27MKO_3=BT245_rna$BT245.KO.EARLY__BT245.KO.C5P4.sorted.bam/(fac8*BT245_rna$genelength)
BT245_rna$BT245_K27MKO_4=BT245_rna$BT245.KO.EARLY__BT245.KO.C4.sorted.bam/(fac9*BT245_rna$genelength)
BT245_rna$BT245_K27MKO_5=BT245_rna$BT245.KO.EARLY__BT245.KO.C2.sorted.bam/(fac10*BT245_rna$genelength)


names(BT245_rna)
BT245_rna=BT245_rna[c(1:7,19:28)]

write.table(BT245_rna,"BT245_RPKM.txt")
