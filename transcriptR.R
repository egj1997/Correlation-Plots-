#ring1b differential peaks
library("rtracklayer")
library("GenomicRanges")
library("ChIPseeker")
library(rtracklayer)
library(dplyr)
#import Ring1b  peak set
setwd("C:/Users/egj19/OneDrive - McGill University/Desktop/Analysis/peakcalls")
library(readxl)
chipseeker <- read_excel("chipseeker.xlsx")
BT245=import.bed("BT245_RING1b_peaks_optimized.bed")
DIPG13=import.bed("DIPG13_RING1b_peaks_optimized.bed")
BT245.ko=import.bed("BT245_C2P8_RING1b_peaks_optimized.bed")
DIPG13.ko=import.bed("DIPG13_C5P6_RING1b_peaks_optimized.bed")

k27m=suppressWarnings(subsetByOverlaps(BT245,DIPG13))
k27m
ko=suppressWarnings(subsetByOverlaps(DIPG13.ko,BT245.ko))
ko
export.bed(ko,"KO_RING1b.bed")
#differential peaks
maintained.ring1b=intersect(k27m,ko)
gainedpeaks.k27m=setdiff(k27m,subsetByOverlaps(k27m,maintained.ring1b))
lostpeaks.k27m=setdiff(ko,subsetByOverlaps(ko,maintained.ring1b))

test1=setdiff(gainedpeaks.k27m,subsetByOverlaps(gainedpeaks.k27m,BT245_K27ac))
test2=setdiff(test1,subsetByOverlaps(test1,BT245_K27MKO_K27ac))
test3=setdiff(test2,subsetByOverlaps(test2,DIPG13_K27ac))

#import SUZ12 peaks
BT245_suz12=read.table("BT245_SUZ12.narrowPeak.bed",header=F)
BT245_suz12=makeGRangesFromDataFrame(BT245_suz12,seqnames.field = "V1",start.field = "V2",end.field = "V3")

DIPG13_suz12=read.table("DIPGXIII_SUZ12.narrowPeak.bed",header=F)
DIPG13_suz12=makeGRangesFromDataFrame(DIPG13_suz12,seqnames.field = "V1",start.field = "V2",end.field = "V3")

#Load K27me3 data
b.k27me3.k27m=import.bed("BT245-C24_H3K27me3_Inp_peaks.broadPeak.bed")
d.k27me3.k27m=import.bed("DIPGXIII-C12-Rx_cells_ChIP1_H3K27me3_1_peaks.broadPeak.bed")
b.k27me3.ko=import.bed("BT245-ko-c2p3-Rx-C_cells_ChIP1_H3K27me3_1_peaks.broadPeak.bed")
d.k27me3.ko=import.bed("DIPGXIII-ko-c5p5-Rx-C_cells_ChIP1_H3K27me3_1_peaks.broadPeak.bed")

k27me3.k27m=suppressWarnings(subsetByOverlaps(b.k27me3.k27m,d.k27me3.k27m))
k27me3.ko=suppressWarnings(subsetByOverlaps(b.k27me3.ko,d.k27me3.ko))

#Load K27ac data
b.k27ac.k27m=import.bed("BT245-nko-c1p5-Rx_cells_ChIP1_H3K27ac_1_peaks.narrowPeak.bed")
b.k27ac.ko=import.bed("BT245-ko-c4p6-Rx_cells_ChIP1_H3K27ac_1_peaks.narrowPeak.bed")
d.k27ac.k27m=import.bed("DIPGXIIIp11_cells_ChIP1_H3K27ac_1_peaks.narrowPeak.bed")
d.k27ac.ko=import.bed("DIPGXIII-ko-c5p5-Rx_cells_ChIP1_H3K27ac_1_peaks.narrowPeak.bed")
k27ac.k27m=suppressWarnings(subsetByOverlaps(b.k27ac.k27m,d.k27ac.k27m))
k27ac.ko=suppressWarnings(subsetByOverlaps(b.k27ac.ko,d.k27ac.ko))

#clean intersects between k27ac and k27me3
k27me3.k27m=setdiff(k27me3.k27m,subsetByOverlaps(k27me3.k27m,intersect(k27me3.k27m,k27ac.k27m)))
k27ac.k27m=setdiff(k27ac.k27m,subsetByOverlaps(k27ac.k27m,intersect(k27ac.k27m,k27me3.k27m)))
k27me3.ko=setdiff(k27me3.ko,subsetByOverlaps(k27me3.ko,intersect(k27me3.ko,k27ac.ko)))
k27ac.ko=setdiff(k27ac.ko,subsetByOverlaps(k27ac.ko,intersect(k27ac.ko,k27me3.ko)))

#loading H2AK119ub data
b.h2ak.k27m=import.bed("BT245-P47-H2AK119ub_peaks.broadPeak.bed")
d.h2ak.k27m=import.bed("DIPG13_P66_Rx_XChIP_H2AK119ub_peaks.broadPeak.bed")
h2ak.k27m=suppressWarnings(subsetByOverlaps(b.h2ak.k27m,d.h2ak.k27m))
h2ak.ko=suppressWarnings(subsetByOverlaps(b.h2ak.ko,d.h2ak.ko))


#differential regions for k27me3
gainedk27me3.k27m=suppressWarnings(setdiff(k27me3.k27m,subsetByOverlaps(k27me3.k27m,k27me3.ko)))
lostk27me3.k27m=suppressWarnings(setdiff(k27me3.ko,subsetByOverlaps(k27me3.ko,k27me3.k27m)))
maintaink27me3.k27m=subsetByOverlaps(k27me3.k27m,k27me3.ko)
#checks
intersect(gainedk27me3.k27m,maintaink27me3.k27m)
intersect(lostk27me3.k27m,maintaink27me3.k27m)
length(k27me3.ko)-length(lostk27me3.k27m)-length(subsetByOverlaps(k27me3.ko,k27me3.k27m))

#differential regions for k27ac
gainedk27ac.k27m=suppressWarnings(setdiff(k27ac.k27m,subsetByOverlaps(k27ac.k27m,k27ac.ko)))
lostk27ac.k27m=suppressWarnings(setdiff(k27ac.ko,subsetByOverlaps(k27ac.ko,k27ac.k27m)))
maintaink27ac.k27m=subsetByOverlaps(k27ac.k27m,k27ac.ko)
#checks
intersect(gainedk27ac.k27m,maintaink27ac.k27m)
intersect(lostk27ac.k27m,maintaink27ac.k27m)
length(k27ac.ko)-length(lostk27ac.k27m)-length(subsetByOverlaps(k27ac.ko,k27ac.k27m))

#check for clean intersects of k27me3 and k27ac
intersect(maintaink27me3.k27m,maintaink27ac.k27m)
intersect(gainedk27ac.k27m,gainedk27me3.k27m)
intersect(lostk27ac.k27m,lostk27me3.k27m)

#dataset for gained ring1b sites
gaink27me3.gainring1b=suppressWarnings(subsetByOverlaps(gainedpeaks.k27m,(subsetByOverlaps(gainedpeaks.k27m,gainedk27me3.k27m))))
gaink27ac.gainring1b=suppressWarnings(subsetByOverlaps(gainedpeaks.k27m,subsetByOverlaps(gainedpeaks.k27m,gainedk27ac.k27m)))
maintaink27me3.gainring1b=suppressWarnings(subsetByOverlaps(gainedpeaks.k27m,subsetByOverlaps(gainedpeaks.k27m,maintaink27me3.k27m)))
maintaink27ac.gainring1b=suppressWarnings(subsetByOverlaps(gainedpeaks.k27m,subsetByOverlaps(gainedpeaks.k27m,maintaink27ac.k27m)))

gainedring1b_k27me3=suppressWarnings(subsetByOverlaps(gainedpeaks.k27m,subsetByOverlaps(gainedpeaks.k27m,k27me3.k27m)))
gainedring1b_k27ac=suppressWarnings(subsetByOverlaps(gainedpeaks.k27m,subsetByOverlaps(gainedpeaks.k27m,k27ac.k27m)))
a6=length(gainedring1b_k27me3)/length(gainedpeaks.k27m)
a7=length(gainedring1b_k27ac)/length(gainedpeaks.k27m)

a1=length(gaink27me3.gainring1b)/length(gainedpeaks.k27m)
a2=length(gaink27ac.gainring1b)/length(gainedpeaks.k27m)
a3=length(maintaink27me3.gainring1b)/length(gainedpeaks.k27m)
a4=length(maintaink27ac.gainring1b)/length(gainedpeaks.k27m)
a5=1-(a1+a2+a3+a4)

#create piechart of ring1b gained sites
# Load ggplot2
library(ggplot2)
library(dplyr)
# Create Data
pie1 <- data.frame(group=c("Gained H3K27me3","Gained H3K27ac","Maintained H3K27me3","Maintained H3K27ac","Other"),value=c(a1,a2,a3,a4,a5))

# Basic piechart

pie1 <- pie1 %>% 
  arrange(desc(group)) %>%
  mutate(prop = value / sum(pie1$value) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )


ggplot(pie1, aes(x="", y=prop, fill=group)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  ggtitle("Gained RING1b peaks") +
  
  geom_text(aes(y = ypos, label = round(value*100,2)), color = "white", size=6) +
  scale_fill_brewer(palette="Set1")+
  theme(legend.background = element_rect(size=4, linetype="solid"))


#dataset for lost ring1b sites
lostk27me3.lostring1b=suppressWarnings(subsetByOverlaps(lostpeaks.k27m,subsetByOverlaps(lostpeaks.k27m,lostk27me3.k27m)))
lostk27ac.lostring1b=suppressWarnings(subsetByOverlaps(lostpeaks.k27m,subsetByOverlaps(lostpeaks.k27m,lostk27ac.k27m)))
maintaink27me3.lostring1b=suppressWarnings(subsetByOverlaps(lostpeaks.k27m,subsetByOverlaps(lostpeaks.k27m,maintaink27me3.k27m)))
maintaink27ac.lostring1b=suppressWarnings(subsetByOverlaps(lostpeaks.k27m,subsetByOverlaps(lostpeaks.k27m,maintaink27ac.k27m)))
gaink27me3.lostring1b=suppressWarnings(subsetByOverlaps(lostpeaks.k27m,gainedk27me3.k27m))

b1=length(lostk27me3.lostring1b)/length(lostpeaks.k27m)
b2=length(lostk27ac.lostring1b)/length(lostpeaks.k27m)
b3=length(maintaink27me3.lostring1b)/length(lostpeaks.k27m)
b4=length(maintaink27ac.lostring1b)/length(lostpeaks.k27m)
b5=1-(b1+b2+b3+b4)

#create piechart for ring1b lost sites
# Create Data
pie2 <- data.frame(group=c("Lost H3K27me3","Lost H3K27ac","Maintained H3K27me3","Maintained H3K27ac","Other"),value=c(b1,b2,b3,b4,b5))

# Basic piechart
pie2 <- pie2 %>% 
  arrange(desc(group)) %>%
  mutate(prop = value / sum(pie2$value) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

ggplot(pie2, aes(x="", y=prop, fill=group)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  ggtitle("Lost RING1b sites") +
  geom_text(aes(y = ypos, label = round(value*100,2)), color = "white", size=5) +
  scale_fill_brewer(palette="Set1")+
  theme(legend.background = element_rect(size=2, linetype="solid"))

#RING1b.peaks association in k27m
RING1b_k27me3.k27m=subsetByOverlaps(k27m,subsetByOverlaps(k27m,k27me3.k27m))
RING1b_k27ac.k27m=subsetByOverlaps(k27m,subsetByOverlaps(k27m,k27ac.k27m))

d1=length(RING1b_k27me3.k27m)/length(k27m)
d2=length(RING1b_k27ac.k27m)/length(k27m)
d3=1-(d2+d1)
# Create Data
pie3 <- data.frame(group=c("H3K27me3","H3K27ac","Other"),value=c(d1,d2,d3))

#create pie chart
pie3 <- pie3 %>% 
  arrange(desc(group)) %>%
  mutate(prop = value / sum(pie3$value) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

ggplot(pie3, aes(x="", y=prop, fill=group)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  ggtitle("RING1b sites in K27M")+
  geom_text(aes(y = ypos, label = round(value*100,2)), color = "white", size=5) +
  scale_fill_brewer(palette="Set1")+
  theme(legend.background = element_rect(size=2, linetype="solid"))

#RING1b.peaks association in KO
RING1b_k27me3.ko=subsetByOverlaps(ko,subsetByOverlaps(ko,k27me3.ko))
RING1b_k27ac.ko=subsetByOverlaps(ko,subsetByOverlaps(ko,k27ac.ko))

e1=length(RING1b_k27me3.ko)/length(ko)
e2=length(RING1b_k27ac.ko)/length(ko)
e3=1-(e2+e1)
# Create Data
pie4 <- data.frame(group=c("H3K27me3","H3K27ac","Other"),value=c(e1,e2,e3))

#create pie chart
pie4 <- pie4 %>% 
  arrange(desc(group)) %>%
  mutate(prop = value / sum(pie3$value) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

ggplot(pie4, aes(x="", y=prop, fill=group)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  ggtitle("RING1b sites in K27M KO")+
  geom_text(aes(y = ypos, label = round(value*100,2)), color = "white", size=5) +
  scale_fill_brewer(palette="Set1")+
  theme(legend.background = element_rect(size=2, linetype="solid"))


#SUZ12
b.suz12.k27m.1=import.bed("BT245p11_SUZ12_Inp_peaks.narrowPeak.bed")
b.suz12.ko.1=import.bed("BT245-ko-c4p9-Rx_cells_ChIP1_SUZ12_1_peaks.narrowPeak.bed")
d.suz12.ko.1=import.bed("DIPGXIII-ko-c5p5-Rx_cells_ChIP1_SUZ12_1_peaks.narrowPeak.bed")
d.suz12.k27m.1=import.bed("DIPG13p14_SUZ12_Inp_peaks.narrowPeak.bed")
b.suz12.k27m.2 <- readPeakFile(chipseeker[[5]])
b.suz12.ko.2 <- readPeakFile(chipseeker[[6]])
d.suz12.ko.2 <- readPeakFile(chipseeker[[8]])
d.suz12.k27m.2 <- readPeakFile(chipseeker[[7]])

b.suz12.k27m=subsetByOverlaps(b.suz12.k27m.1,b.suz12.k27m.2)
d.suz12.k27m=subsetByOverlaps(d.suz12.k27m.1,d.suz12.k27m.2)
b.suz12.ko=b.suz12.ko.2
d.suz12.ko=d.suz12.ko.2
suz12.k27m=subsetByOverlaps(b.suz12.k27m.2,d.suz12.k27m.2)

4\suz12.ko=subsetByOverlaps(b.suz12.ko,d.suz12.ko)

#RING1b, SUZ12 distribution in k27m
RING1b_SUZ12=subsetByOverlaps(k27m,suz12.k27m)
RING1b_only=setdiff(k27m,subsetByOverlaps(k27m,RING1b_SUZ12))
SUZ12_only=setdiff(suz12.k27m,subsetByOverlaps(suz12.k27m,RING1b_SUZ12))

total2=length(SUZ12_only)+length(RING1b_only)+length(RING1b_SUZ12)

f1=length(RING1b_SUZ12)*100/total2
f2=length(RING1b_only)*100/total2
f3=length(SUZ12_only)*100/total2
f1+f2+f3

#RING1b, SUZ12 distribution in ko
RING1b_SUZ12_ko=subsetByOverlaps(ko,suz12.ko)
RING1b_only_ko=setdiff(k27m,subsetByOverlaps(ko,RING1b_SUZ12_ko))
SUZ12_only_ko=setdiff(suz12.k27m,subsetByOverlaps(suz12.ko,RING1b_SUZ12_ko))

total3=length(SUZ12_only_ko)+length(RING1b_only_ko)+length(RING1b_SUZ12_ko)

f4=length(RING1b_SUZ12_ko)*100/total3
f5=length(RING1b_only_ko)*100/total3
f6=length(SUZ12_only_ko)*100/total3
f4+f5+f6
#create barplot
# library
library(ggplot2)
library(viridis)
library(hrbrthemes)
# create a dataset
specie <- c(rep("K27M" , 3) , rep("K27M KO" , 3))
condition <- rep(c("Both RING1b and SUZ12" , "RING1b Only" , "SUZ12 Only") , 2)
value <- abs(c(f1,f2,f3,f4,f5,f6))
data <- data.frame(specie,condition,value)

# Grouped
ggplot(data, aes(fill=condition, y=value, x=specie)) + 
  geom_bar(position="dodge", stat="identity") +
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme_bw()+
  xlab("Condition") +
  ylab("Percentage of sites")+
  theme(axis.text.x = element_text(face = "plain", color = "black", size = 10),axis.title.y = element_text(size=18),legend.title=element_blank(),legend.text = element_text(face="italic",size=14))

#are these overlapping sites rich in k27me3?
answer2=length(subsetByOverlaps(RING1b_SUZ12,gainedk27me3.k27m))/length(RING1b_SUZ12)
answer3=length(subsetByOverlaps(RING1b_SUZ12,maintaink27me3.k27m))/length(RING1b_SUZ12)
#how many of these sites are new? How many have gained k27me3
subsetByOverlaps(RING1b_SUZ12,setdiff(RING1b_SUZ12,RING1b_SUZ12_ko))
answer4=length(subsetByOverlaps(subsetByOverlaps(RING1b_SUZ12,setdiff(RING1b_SUZ12,RING1b_SUZ12_ko)),gainedk27me3.k27m))/length(subsetByOverlaps(RING1b_SUZ12,setdiff(RING1b_SUZ12,RING1b_SUZ12_ko)))
answer5=length(subsetByOverlaps(subsetByOverlaps(RING1b_SUZ12,setdiff(RING1b_SUZ12,RING1b_SUZ12_ko)),maintaink27me3.k27m))/length(subsetByOverlaps(RING1b_SUZ12,setdiff(RING1b_SUZ12,RING1b_SUZ12_ko)))
answer6=length(subsetByOverlaps(subsetByOverlaps(RING1b_SUZ12,setdiff(RING1b_SUZ12,RING1b_SUZ12_ko)),lostk27ac.k27m))/length(subsetByOverlaps(RING1b_SUZ12,setdiff(RING1b_SUZ12,RING1b_SUZ12_ko)))
answer7=length(subsetByOverlaps(subsetByOverlaps(RING1b_SUZ12,setdiff(RING1b_SUZ12,RING1b_SUZ12_ko)),maintaink27ac.k27m))/length(subsetByOverlaps(RING1b_SUZ12,setdiff(RING1b_SUZ12,RING1b_SUZ12_ko)))

#Differential SUZ12 sites
maintainsuz12=suppressWarnings(subsetByOverlaps(suz12.k27m,suz12.ko))
gainedsuz12=setdiff(suz12.k27m,maintainsuz12)
lostsuz12=setdiff(suz12.ko,subsetByOverlaps(suz12.ko,suz12.k27m))
#check
intersect(maintainsuz12,gainedsuz12)
intersect(maintainsuz12,lostsuz12)
length(suz12.ko)-length(subsetByOverlaps(suz12.ko,suz12.k27m))-length(lostsuz12)

gainedsuz12_gainedring1b=suppressWarnings(subsetByOverlaps(gainedsuz12,gainedpeaks.k27m))
maintainsuz12_gainedring1b=suppressWarnings(subsetByOverlaps(maintainsuz12,gainedpeaks.k27m))
lostsuz12_lostring1b=suppressWarnings(subsetByOverlaps(lostsuz12,lostpeaks.k27m))
maintainsuz12_lostring1b=suppressWarnings(subsetByOverlaps(maintainsuz12,lostpeaks.k27m))
gainedsuz12_lostring1b=suppressWarnings(subsetByOverlaps(gainedsuz12,lostpeaks.k27m))
lostsuz12_gainedring1b=suppressWarnings(subsetByOverlaps(lostsuz12,gainedpeaks.k27m))
maintainsuz12_maintainring1b=suppressWarnings(subsetByOverlaps(maintainsuz12,maintainedring1b))
gainedsuz12_maintainring1b=suppressWarnings(subsetByOverlaps(gainedsuz12,maintainedring1b))
lostsuz12_maintainedring1b=suppressWarnings(subsetByOverlaps(lostsuz12,maintainedring1b))

total=length(gainedsuz12_gainedring1b)+length(maintainsuz12_gainedring1b)+length(lostsuz12_lostring1b)+length(maintainsuz12_lostring1b)+length(gainedsuz12_lostring1b)+length(gainedsuz12_lostring1b)+length(lostsuz12_gainedring1b)+length(maintainsuz12_maintainring1b)+length(gainedsuz12_maintainring1b)+length(lostsuz12_maintainedring1b)

c1=length(gainedsuz12_gainedring1b)/total
c2=length(maintainsuz12_gainedring1b)/total
c3=length(lostsuz12_lostring1b)/total
c4=length(maintainsuz12_lostring1b)/total
c5=length(gainedsuz12_lostring1b)/total
c6=length(lostsuz12_gainedring1b)/total
c7=length(maintainsuz12_maintainring1b)/total
c8=length(gainedsuz12_maintainring1b)/total
c9=length(lostsuz12_maintainedring1b)/total
c10=sum(c1,c2,c3,c4,c5,c6,c7,c8,c9)

#create piechart for ring1b lost sites
# Create Data
pie3 <- data.frame(group=c("Gained both SUZ12 and RING1b","Maintained SUZ12 and gained RING1b","Lost both SUZ12 and RING1b","Maintain SUZ12 but lost RING1b","Gained SUZ12 but lost RING1b","Lost SUZ12 but gained RING1b","Maintained both SUZ12 and RING1b","Gained SUZ12 and maintained RING1b","Lost SUZ12 and maintained RING1b"),value=c(c1,c2,c3,c4,c5,c6,c7,c8,c9))

# Basic piechart
pie3 <- pie3 %>% 
  arrange(desc(group)) %>%
  mutate(prop = value / sum(pie1$value) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

ggplot(pie3, aes(x="", y=prop, fill=group)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  geom_text(aes(y = ypos, label = round(value*100,2)), color = "white", size=5) +
  scale_fill_brewer(palette="Set1")+
  theme(legend.background = element_rect(size=2, linetype="solid"))

#sites that maintain both suz12 and ring1b
k27me3_maintainsuz12_maintainring1b=subsetByOverlaps(maintainsuz12_maintainring1b,maintaink27me3.k27m)
k27ac_maintainsuz12_maintainring1b=subsetByOverlaps(maintainsuz12_maintainring1b,maintaink27ac.k27m)
gainedk27me3_maintainsuz12_maintainring1b=subsetByOverlaps(maintainsuz12_maintainring1b,gainedk27me3.k27m)
gainedk27ac_maintainsuz12_maintainring1b=subsetByOverlaps(maintainsuz12_maintainring1b,gainedk27ac.k27m)
lostk27me3_maintainsuz12_maintainring1b=subsetByOverlaps(maintainsuz12_maintainring1b,lostk27me3.k27m)
lostk27ac_maintainsuz12_maintainring1b=subsetByOverlaps(maintainsuz12_maintainring1b,lostk27ac.k27m)

length(k27me3_maintainsuz12_maintainring1b)/length(maintainsuz12_maintainring1b)
length(k27ac_maintainsuz12_maintainring1b)/length(maintainsuz12_maintainring1b)

length(gainedk27me3_maintainsuz12_maintainring1b)/length(maintainsuz12_maintainring1b)
length(gainedk27ac_maintainsuz12_maintainring1b)/length(maintainsuz12_maintainring1b)


#create a treemap of suz12 ring1b gained sites
library(ggplot2)
library(treemapify)

treemap1=data.frame(group=c("Gained H3K27me3","Gained SUZ12","Maintained H3K27me3","Maintained H3K27ac"),values=c(length(gainedk27me3_gainedsuz12_gainedring1b),length(gainedk27ac_gainedsuz12_gainedring1b),length(k27me3_gainedsuz12_gainedring1b),length(k27ac_gainedsuz12_gainedring1b)))
treemap1$percentage=treemap1$values/sum(treemap1$values)*100

ggplot(treemap1, aes(area = values, fill = group,label=round(percentage,2))) +
  geom_treemap()+
geom_treemap_text(fontface = "bold", colour = "white", place = "centre",size=10,
                  grow = TRUE)

#sites with k27me3 and ring1b
gainedk27me3_gainedring1b=intersect(gainedpeaks.k27m,subsetByOverlaps(gainedk27me3.k27m,gainedpeaks.k27m))

export.bed(gainedpeaks.k27m,"gainedk27me3_gainedring1b.bed")

#run go
library(CompGO)
annotateBedFromDb(gRanges = gainedk27me3_gainedring1b, db = TxDb.Hsapiens.UCSC.hg19.knownGene,
                  window = 5000)


#ring1b and cbx2
b.cbx2.k27m=import.bed("BT245_CBX2_peaks_seqmonk.bed")
ring1b_cbx2=subsetByOverlaps(BT245,b.cbx2.k27m)
#overlap with k27me3
length(subsetByOverlaps(ring1b_cbx2,b.k27me3.k27m))/length(ring1b_cbx2)
#overlap with k27ac
length(suppressWarnings(subsetByOverlaps(ring1b_cbx2,b.k27ac.k27m)))/length(ring1b_cbx2)
#overlap with H2AK119ub
b.h2ak119ub.k27m=import.bed("BT245-P47-H2AK119ub_peaks.broadPeak.bed")
length(suppressWarnings(subsetByOverlaps(ring1b_cbx2,b.h2ak119ub.k27m)))/length(ring1b_cbx2)
#overlap with suz12
length(suppressWarnings(subsetByOverlaps(ring1b_cbx2,b.suz12.k27m)))/length(ring1b_cbx2)






#overlap of gained ring1b peaks with k27me3
test1=read.table("seqmonk18.txt",header=T)
test1[4]=test1[3]-test1[2]
hist(table(test1[4]))

export.bed(test1,"GBM002_RING1b_peaks_seqmonk.bed")

test1=import.bed("BT245_RING1b_peaks_seqmonk.bed")
test2=import.bed("DIPG13_RING1b_peaks_seqmonk.bed")
test1=makeGRangesFromDataFrame(test1)
test2=makeGRangesFromDataFrame(test2)

k27m.test=subsetByOverlaps(test1,test2)

test3=makeGRangesFromDataFrame(import.bed("BT245_C2P8_RING1b_peaks_seqmonk.bed"))                               
test4=makeGRangesFromDataFrame(import.bed("DIPG13_C5P6_RING1b_peaks_seqmonk.bed"))                               

ko.test=subsetByOverlaps(test3,test4)

gainedpeaks.test=setdiff(k27m.test,ko.test)
gainedpeaks.test

lostpeaks.test=setdiff(ko.test,k27m.test)
lostpeaks.test

#create enhancer regions
d.k4me1.k27m=import.bed("DIPG13_K4me1_peaks_seqmonk.bed")
d.k4me1.k27m
d.k27ac.k27m=import.bed("DIPG13_K27ac_peaks_seqmonk.bed")
d.k27ac.k27m
d.enhancer.k27m=subsetByOverlaps(d.k27ac.k27m,intersect(d.k4me1.k27m,d.k27ac.k27m))
enhancer.k27m=read.table("enhancers.txt",header=T)
enhancer.k27m=makeGRangesFromDataFrame(enhancer.k27m)
promoter.k27m=read.table("promoters.txt",header=T)
promoter.k27m=makeGRangesFromDataFrame(promoter)
#check for enhancer changes in DIPG13
enhancer_RING1b_k27m=subsetByOverlaps(K27M,subsetByOverlaps(K27M,enhancer.k27m))
enhancer_RING1b_ko=subsetByOverlaps(KO,subsetByOverlaps(KO,enhancer.k27m))
enhancer_gainedpeaks=subsetByOverlaps(gainedpeaks.k27m,subsetByOverlaps(gainedpeaks.k27m,enhancer.k27m))
enhancer_maintainedpeaks=subsetByOverlaps(maintainedring1b,subsetByOverlaps(maintainedring1b,enhancer.k27m))

promoter_RING1b_k27m=subsetByOverlaps(K27M,subsetByOverlaps(K27M,promoter.k27m))
promoter_RING1b_ko=subsetByOverlaps(KO,subsetByOverlaps(KO,promoter.k27m))
promoter_gainedpeaks=subsetByOverlaps(gainedpeaks.k27m,subsetByOverlaps(gainedpeaks.k27m,promoter.k27m))
promoter_maintainedpeaks=subsetByOverlaps(maintainedring1b,subsetByOverlaps(maintainedring1b,promoter.k27m))


length(enhancer_RING1b_k27m)/length(K27M)
length(enhancer_RING1b_ko)/length(KO)
length(enhancer_gainedpeaks)/length(gainedpeaks.k27m)
length(enhancer_maintainedpeaks)/length(maintainedring1b)
length(enhancer.k27m)/length(k27ac.k27m)

length(promoter_RING1b_k27m)/length(K27M)
length(promoter_RING1b_ko)/length(KO)
length(promoter_gainedpeaks)/length(gainedpeaks.k27m)
length(promoter_maintainedpeaks)/length(maintainedring1b)


export.bed(enhancer_gainedpeaks,"enhancer_gainedpeaks.bed")

#RING1b and H2AK119ub association
length(subsetByOverlaps(k27m,h2ak.k27m))/length(k27m)
length(subsetByOverlaps(h2ak.k27m,k27me3.k27m))/length(h2ak.k27m)
length(subsetByOverlaps(h2ak.k27m,k27ac.k27m))/length(h2ak.k27m)

length(subsetByOverlaps(k27m,subsetByOverlaps(k27me3.k27m,h2ak.k27m)))/length(k27m)
length(subsetByOverlaps(ko,h2ak.ko))/length(ko)

#GBM002
g.ring1b=import.bed("GBM002_RING1b_peaks_seqmonk.bed")
k27me3.g34r=import.bed("GBM002p49_H3K27me3_Input_peaks.broadPeak.bed")
length(intersect(g.ring1b,k27me3.g34r))/length(g.ring1b)

#superenhancers
se=read_excel("superenhancers.xlsx")
se=se[se$p.value<c(0.05),]
se=se[se$`Fold (K27M/WT)`>0,]
se=makeGRangesFromDataFrame(se)
length(subsetByOverlaps(se,k27m))/length(se)
length(subsetByOverlaps(se,gainedpeaks.k27m))/length(se)
length(subsetByOverlaps(se,ko))/length(se)

k27m.table=as.data.frame(k27m)
sum(k27m.table$width)/(3*10^9)

#vien diagrams per cell line
library(VennDiagram)

grid.newpage()
draw.pairwise.venn(length(BT245), length(BT245.ko), length(intersect(BT245,BT245.ko)), category = c("Parental", "K27M KO"), lty = rep("blank", 
                                                                                   2), fill = c("red", "light blue"), alpha = rep(0.5, 2), cat.pos = c(0, 
                                                                                                                                                        0), cat.dist = rep(0.025, 2))

grid.newpage()
draw.pairwise.venn(length(DIPG13), length(DIPG13.ko), length(intersect(DIPG13,DIPG13.ko)), category = c("Parental", "K27M KO"), lty = rep("blank", 
                                                                                                                                          2), fill = c("firebrick1", "dodgerblue2"), alpha = rep(0.5, 2), cat.pos = c(0, 
                                                                                                                                                                                                              0), cat.dist = rep(0.025, 2))
#RING1b peaks that are repressive
#Bt245
a1=subsetByOverlaps(BT245,b245_k27me3_clusteronpeaks)
a2=subsetByOverlaps(BT245,bt245_k27ac_clusteronpeaks)
a3=subsetByOverlaps(BT245.ko,b245_KO_k27me3_clusteronpeaks)
a4=subsetByOverlaps(BT245.ko,b245_KO_k27ac_clusteronpeaks)
a5=intersect(a1,a3)
a6=intersect(a2,a4)
#plot repressive
tiff("Venn_BT245_RING1b_k27me3.tiff",res=600,width = 4, height = 4, units = 'in')

grid.newpage()
draw.pairwise.venn(length(a1), length(a3), length((a5)), category = c("Parental", "K27M KO"), lty = rep("blank", 
                                                                                                                                          2), fill = c("red", "light blue"), alpha = rep(0.5, 2), cat.pos = c(0, 
                                                                                                                                                                                                              0), cat.dist = rep(0.025, 2))
dev.off()
#plot active
tiff("Venn_BT245_RING1b_k27ac.tiff",res=600,width = 4, height = 4, units = 'in')

grid.newpage()
draw.pairwise.venn(length(a2), length(a4), length((a6)), category = c("Parental", "K27M KO"), lty = rep("blank", 
                                                                                                        2), fill = c("red", "light blue"), alpha = rep(0.5, 2), cat.pos = c(0, 
                                                                                                                                                                            0), cat.dist = rep(0.025, 2))
dev.off()
#dipg13
a7=subsetByOverlaps(DIPG13,d.k27me3.k27m)
a8=subsetByOverlaps(DIPG13,d.k27ac.k27m)
a9=subsetByOverlaps(DIPG13.ko,d.k27me3.ko)
a10=subsetByOverlaps(DIPG13.ko,d.k27ac.ko)
a11=intersect(a7,a9)
a12=intersect(a8,a10)
#plot repressive
tiff("Venn_DIPG13_RING1b_k27me3.tiff",res=600,width = 4, height = 4, units = 'in')

grid.newpage()
draw.pairwise.venn(length(a7), length(a9), length((a11)), category = c("Parental", "K27M KO"), lty = rep("blank", 
                                                                                                        2), fill = c("red", "light blue"), alpha = rep(0.5, 2), cat.pos = c(0,  
                                                                                                                                                                            0), cat.dist = rep(0.025, 2))
dev.off()
#plot active
tiff("Venn_DIPG13_RING1b_k27ac.tiff",res=600,width = 4, height = 4, units = 'in')

grid.newpage()
draw.pairwise.venn(length(a8), length(a10), length((a12)), category = c("Parental", "K27M KO"), lty = rep("blank", 
                                                                                                         2), fill = c("red", "light blue"), alpha = rep(0.5, 2), cat.pos = c(0,  
                                                                                                                                                                             0), cat.dist = rep(0.025, 2))  
dev.off()
#using chase enrichment on peaks:
#load data
setwd("C:/Users/egj19/OneDrive - McGill University/Desktop/Analysis")
b245_k27me3_clusteronpeaks=import.gff("b245_k27me3_clusteronpeaks.gff")
bt245_k27ac_clusteronpeaks=import.gff("bt245_k27ac_clusteronpeaks.gff")
b245_KO_k27me3_clusteronpeaks=import.gff("b245_KO_k27me3_clusteronpeaks.gff")
b245_KO_k27ac_clusteronpeaks=import.gff("b245_KO_k27ac_clusteronpeaks.gff")
dipg13_k27ac_clusteronpeaks=import.gff("dipg13_k27ac_clusteronpeaks.gff")
dipg13_k27me3_clusteronpeaks=import.gff("dipg13_k27me3_clusteronpeaks.gff")
dipg13_KO_k27ac_clusteronpeaks=import.gff("dipg13_KO_k27ac_clusteronpeaks.gff")
dipg13_KO_k27me3_clusteronpeaks=import.gff("dipg13_KO_k27me3_clusteronpeaks.gff")
b245_k27me3_clusteronpeaks=makeGRangesFromDataFrame(b245_k27me3_clusteronpeaks)


#plot bt245
#k27me3
tiff("Venn_BT245_RING1b_k27me3.tiff",res=600,width = 4, height = 4, units = 'in')

grid.newpage()
draw.pairwise.venn(length(b245_k27me3_clusteronpeaks), length(b245_KO_k27me3_clusteronpeaks), length(intersect(b245_k27me3_clusteronpeaks,b245_KO_k27me3_clusteronpeaks)), category = c("Parental", "K27M KO"), lty = rep("blank", 
                                                                                                          2), fill = c("red", "light blue"), alpha = rep(0.5, 2), cat.pos = c(0,  
                                                                                                                                                                              0), cat.dist = rep(0.025, 2))
dev.off()
#k27ac                                                                                                                                                                              
tiff("Venn_BT245_RING1b_k27ac.tiff",res=600,width = 4, height = 4, units = 'in')

grid.newpage()
draw.pairwise.venn(length(bt245_k27ac_clusteronpeaks), length(b245_KO_k27ac_clusteronpeaks), length(intersect(bt245_k27ac_clusteronpeaks,b245_KO_k27ac_clusteronpeaks)), category = c("Parental", "K27M KO"), lty = rep("blank", 
                                                                                                                                                                                                                          2), fill = c("red", "light blue"), alpha = rep(0.5, 2), cat.pos = c(0,  
                                                                                                                                                                                                                                                                                              0), cat.dist = rep(0.025, 2))
dev.off()
#dipg13
#k27me3
tiff("Venn_DIPG13_RING1b_k27me3.tiff",res=600,width = 4, height = 4, units = 'in')

grid.newpage()
draw.pairwise.venn(length(dipg13_k27me3_clusteronpeaks), length(dipg13_KO_k27me3_clusteronpeaks), length(intersect(dipg13_k27me3_clusteronpeaks,dipg13_KO_k27me3_clusteronpeaks)), category = c("Parental", "K27M KO"), lty = rep("blank", 
                                                                                                                                                                                                                          2), fill = c("red", "light blue"), alpha = rep(0.5, 2), cat.pos = c(0,  
                                                                                                                                                                                                                                                                                              0), cat.dist = rep(0.025, 2))
dev.off()
#k27ac 
tiff("Venn_DIPG13_RING1b_k27ac.tiff",res=600,width = 4, height = 4, units = 'in')

grid.newpage()
draw.pairwise.venn(length(dipg13_k27ac_clusteronpeaks), length(dipg13_KO_k27ac_clusteronpeaks), length(intersect(dipg13_k27ac_clusteronpeaks,dipg13_KO_k27ac_clusteronpeaks)), category = c("Parental", "K27M KO"), lty = rep("blank", 
                                                                                                                                                                                                                        2), fill = c("red", "light blue"), alpha = rep(0.5, 2), cat.pos = c(0, 
                                                                                                                                                                                                                                                                                            0), cat.dist = rep(0.025, 2)) 
dev.off()
#percentage of RING1b as repressive/activating complex
a1=length(b245_k27me3_clusteronpeaks)*100/length(BT245)
a2=length(bt245_k27ac_clusteronpeaks)*100/length(BT245)
a3=length(b245_KO_k27me3_clusteronpeaks)*100/length(BT245.ko)
a4=length(b245_KO_k27ac_clusteronpeaks)*100/length(BT245.ko)

a5=length(dipg13_k27me3_clusteronpeaks)*100/length(DIPG13)
a6=length(dipg13_k27ac_clusteronpeaks)*100/length(DIPG13)
a7=length(dipg13_KO_k27me3_clusteronpeaks)*100/length(DIPG13.ko)
a8=length(dipg13_KO_k27ac_clusteronpeaks)*100/length(DIPG13.ko)
ggsave("BarPlot_RING1b_association_BT245.tiff", units="in", width=5, height=4, dpi=600, compression = 'lzw')

#create df 
df <- data.frame(Condition=c("BT245 Parental","BT245 Parental", "BT245 K27M KO","BT245 K27M KO"),Complex=c("Repressive", "Activating"),
                 Percentage=c(a1,a2,a3,a4))
head(df)

df1 <- data.frame(Condition=c("DIPGXIII Parental","DIPGXIII Parental", "DIPGXIII K27M KO","DIPGXIII K27M KO"),Complex=c("Repressive", "Activating"),
                 Percentage=c(a5,a6,a7,a8))

library(ggplot2)
p <- ggplot(data=df1, aes(x=Complex, y=Percentage, fill=Condition)) +
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  theme_minimal()+ xlab("PRC1 Complex")+
  geom_text(aes(label=paste(round(Percentage,0))), vjust=1.6, color="white",
            position = position_dodge(0.9), size=5)

p + scale_fill_manual(values=c('dodgerblue3','brown1'))
ggsave("BarPlot_RING1b_association_DIPG13.tiff", units="in", width=5, height=4, dpi=600, compression = 'lzw')


#New peaks association?
b.gainedpeaks=setdiff(BT245,subsetByOverlaps(BT245,b.maintained))
b.maintained=intersect(BT245,BT245.ko)
b.lostpeaks=setdiff(BT245.ko,subsetByOverlaps(BT245.ko,b.maintained))
b1=length(subsetByOverlaps(b.gainedpeaks,bt245_k27ac_clusteronpeaks))/length(b.gainedpeaks)
b2=length(subsetByOverlaps(b.gainedpeaks,b245_k27me3_clusteronpeaks))/length(b.gainedpeaks)

#ashot's paper for ring1b
data.range=makeGRangesFromDataFrame(data)
hits=subsetByOverlaps(data.range,BT245)
hits=as.data.frame(hits)
hits$probe=paste(hits$seqnames,hits$start,hits$end)
data$probe=paste(data$Chromosome,data$Start,data$End)
c1=left_join(hits,data)
sum(c1$H3K27me3_BT245)/sum(data$H3K27me3_BT245)

hits1=subsetByOverlaps(data.range,DIPG13)
hits1=as.data.frame(hits1)
hits1$probe=paste(hits1$seqnames,hits1$start,hits1$end)
data$probe=paste(data$Chromosome,data$Start,data$End)
c2=left_join(hits1,data)
sum(c2$H3K27me3_DIPG13)/sum(data$H3K27me3_DIPG13)

data.range1=makeGRangesFromDataFrame(data2)
hits3=subsetByOverlaps(data.range1,BT245.ko)
hits3=as.data.frame(hits3)
hits3$probe=paste(hits3$seqnames,hits3$start,hits3$end)
data2$probe=paste(data2$Chromosome,data2$Start,data2$End)
c3=left_join(hits3,data2)
sum(c3$H3K27me3_BT245)/sum(data2$H3K27me3_BT245)

hits4=subsetByOverlaps(data.range1,DIPG13.ko)
hits4=as.data.frame(hits4)
hits4$probe=paste(hits4$seqnames,hits4$start,hits4$end)
data2$probe=paste(data2$Chromosome,data2$Start,data2$End)
c4=left_join(hits4,data2)
sum(c3$H3K27me3_DIPG13)/sum(data2$H3K27me3_DIPG13)

#add replicates
data3=read.table("seqmonk21.txt",header=T)
data3$H3K27me3_BT245_1=data3$BT245.nko.c1p5_H3K27me3_1.sorted.dup.bam/data3$BT245.nko.c1p5.Rx_cells_ChIP1_Input_1.sorted.dup.bam
data3$H3K27me3_DIPG13_1=data3$DIPG13.C12_H3K27me3.sorted.dup.bam/data3$DIPG13.C12_Input.sorted.dup.bam
data3$H3K27me3_BT245_KO_1=data3$BT245.ko.c2p3.Rx.C_cells_ChIP1_H3K27me3_1.sorted.dup.bam/data3$BT245.ko.c2p3.Rx_cells_ChIP1_Input_1.sorted.dup.bam
data3$H3K27me3_DIPG13_KO_1=data3$DIPGXIII.ko.c5p5.Rx.C_cells_ChIP1_H3K27me3_1.sorted.dup.bam/data3$DIPGXIII.ko.c5p5.Rx_cells_ChIP1_Input_1.sorted.dup.bam
data3$probe=paste(data3$Chromosome,data3$Start,data3$End)
data3.range=data3[c(2,3,4)]
data3.range=makeGRangesFromDataFrame(data3.range)
hits5=subsetByOverlaps(data3.range,BT245)
hits6=subsetByOverlaps(data3.range,DIPG13)
hits7=subsetByOverlaps(data3.range,BT245.ko)
hits8=subsetByOverlaps(data3.range,DIPG13.ko)
hits5=as.data.frame(hits5)
hits5$probe=paste(hits5$seqnames,hits5$start,hits5$end)
c5=left_join(hits5,data3)
sum(c5$H3K27me3_BT245_1)/sum(data3$H3K27me3_BT245_1)

hits6=as.data.frame(hits6)
hits6$probe=paste(hits6$seqnames,hits6$start,hits6$end)
c6=left_join(hits6,data3)
sum(c6$H3K27me3_DIPG13_1)/sum(data3$H3K27me3_DIPG13_1)

hits7=as.data.frame(hits7)
hits7$probe=paste(hits7$seqnames,hits7$start,hits7$end)
c7=left_join(hits7,data3)
sum(c7$H3K27me3_BT245_KO)/sum(data3$H3K27me3_BT245_KO)

hits8=as.data.frame(hits8)
hits8$probe=paste(hits8$seqnames,hits8$start,hits8$end)
c8=left_join(hits8,data3)
sum(c8$H3K27me3_DIPG13_KO)/sum(data3$H3K27me3_DIPG13_KO)

test=import.gff('testk27ac.gff')
test=makeGRangesFromDataFrame(test)
export.bed(test,"testk27ac.bed")

#import BT245 Haifen Data
b.k27me3.haifen.k27m=import.bed('BT245-nko-c1p5-Rx_H3K27me3-50kb.filtered.haifen.bed')
b.k27me3.haifen.ko=import.bed("BT245-ko-c4p6-Rx_H3K27me3-50kb.filtered.haifen.bed")

maintaink27me3=intersect(b.k27me3.haifen.k27m,b.k27me3.haifen.ko)
gaink27me3=setdiff(b.k27me3.haifen.k27m,subsetByOverlaps(b.k27me3.haifen.k27m,maintaink27me3))
losek27me3=setdiff(b.k27me3.haifen.ko,subsetByOverlaps(b.k27me3.haifen.ko,maintaink27me3))
b.gainedpeaks=setdiff(BT245,subsetByOverlaps(BT245,intersect(BT245,BT245.ko)))
b.lostpeaks=setdiff(BT245.ko,subsetByOverlaps(BT245.ko,intersect(BT245,BT245.ko)))

a1=subsetByOverlaps(b.gainedpeaks,b.k27me3.haifen.k27m)
length(a1)/length(b.gainedpeaks)
a2=subsetByOverlaps(b.gainedpeaks,b.k27ac.k27m)
length(a2)/length(b.gainedpeaks)

a3=subsetByOverlaps(b.lostpeaks,b.k27me3.haifen.ko)
length(a3)/length(b.lostpeaks)
a4=subsetByOverlaps(b.lostpeaks,b.k27ac.ko)
length(a4)/length(b.lostpeaks)

#size of peak calls
setwd("C:/Users/egj19/OneDrive - McGill University/Desktop/Analysis/peakcalls")
library(rtracklayer)
library(readxl)
library('ChIPseeker')
library(GenomicRanges)
#broad peaks
BT245 <- import.bed("BT245_RING1b_peaks_optimized.bed")
DIPG13=import.bed("DIPG13_RING1b_peaks_seqmonk.bed")
BT245.ko=import.bed("BT245_C2P8_RING1b_peaks_optimized.bed")
DIPG13.ko=import.bed("DIPG13_C5P6_RING1b_peaks_seqmonk.bed")

#ectopic distribution of RING1b
BT245_maintainedpeaks=subsetByOverlaps(BT245, BT245.ko)
BT245_maintainedpeaks
BT245_gainedpeaks=setdiff(BT245,subsetByOverlaps(BT245,BT245_maintainedpeaks))
BT245_gainedpeaks

BT245_lostpeaks=setdiff(BT245.ko,subsetByOverlaps(BT245.ko,BT245_maintainedpeaks))
BT245_lostpeaks


#import annotated tables
BT245_gainedpeaks.annotation=read.table('BT245_gainedpeaks.annotation.txt',header=T)
merge1=BT245_gainedpeaks.annotation
names(merge1)=paste(c("Chromosome","Start","End","Annotation","Probe"))
BT245_lostpeaks.annotation=read.table('BT245_lostpeaks.annotation.txt',header=T)
merge2=BT245_lostpeaks.annotation
names(merge2)=paste(c("Chromosome","Start","End","Annotation","Probe"))

#load DEG info
setwd("C:/Users/egj19/OneDrive - McGill University/Desktop/Analysis/peakcalls")
BT245_deg=read.table("DESeq2_BT245.txt",header=T)
names(BT245_deg)
#How many of these new sites change in expression?
library(dplyr)
BT245_gainedpeaks.annotation.rna=left_join(merge1,BT245_deg,by="Probe")
BT245_gainedpeaks.annotation.rna=BT245_gainedpeaks.annotation.rna[complete.cases(BT245_gainedpeaks.annotation.rna),]

#How many of these lost sites change in expression
BT245_lostpeaks.annotation.rna=left_join(merge2,BT245_deg,by="Probe")
BT245_lostpeaks.annotation.rna=BT245_lostpeaks.annotation.rna[complete.cases(BT245_lostpeaks.annotation.rna),]

#volcano plot
library(ggplot2)
library(EnhancedVolcano)
names(BT245_gainedpeaks.annotation.rna)
res=BT245_gainedpeaks.annotation.rna
names(res)
res$Log2_FC=-res$Log2_FC
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


#How many of gained K27M ring1b peaks have k27ac?
setwd("C:/Users/egj19/OneDrive - McGill University/Desktop/Analysis/peakcalls")
BT245_K27ac_annotation=read.table("BT245_K27ac_annotation.txt",header=T)
BT245_K27ac_annotation$K27ac=paste("k27ac")
BT245_K27ac_range=makeGRangesFromDataFrame(BT245_K27ac_annotation)
BT245_gainedpeaks.range=makeGRangesFromDataFrame(BT245_gainedpeaks.annotation.rna,seqnames.field = "Chromosome.x",start.field = "Start.x",end.field = "End.x")
hits=as.list(findOverlaps(BT245_gainedpeaks.range,BT245_K27ac_range))
BT245_gainedpeaks.annotation.rna$K27ac=sapply(extractList(BT245_K27ac_annotation$K27ac,hits),paste)
#how many of these sites have k27me3 retained k27me3?
b.k27me3=import.bed("BT245-nko-c1p5-Rx_H3K27me3-50kb.filtered.haifen.bed")
b.ko.k27me3=import.bed('BT245-ko-c4p6-Rx_H3K27me3-50kb.filtered.haifen.bed')

bt245_k27me3_maintained=subsetByOverlaps(b.k27me3,b.ko.k27me3)
bt245_k27me3_maintained$K27me3="K27me3"

BT245_K27me3_range=makeGRangesFromDataFrame(bt245_k27me3_maintained)
BT245_gainedpeaks.range=makeGRangesFromDataFrame(BT245_gainedpeaks.annotation.rna,seqnames.field = "Chromosome.x",start.field = "Start.x",end.field = "End.x")
hits=as.list(findOverlaps(BT245_gainedpeaks.range,BT245_K27me3_range))
BT245_gainedpeaks.annotation.rna$K27me3=sapply(extractList(bt245_k27me3_maintained$K27me3,hits),paste)

#import loop tables
BT245_loops=read.table("BT245_FDR80_merged_loops.bed",header=F)
BT245_c2_loops=read.table("BT245_c2a_FDR80_merged_loops.bed",header=F)

BT245_loops_all=makeGRangesFromDataFrame(BT245_loops,seqnames.field = "V1",start.field = "V2",end.field = "V3")
BT245_loops_c2_all=makeGRangesFromDataFrame(BT245_c2_loops,seqnames.field = "V1",start.field = "V2",end.field = "V3")

#add loop dimension
hits=as.list(findOverlaps(makeGRangesFromDataFrame(BT245_gainedpeaks.annotation.rna,seqnames.field = "Chromosome.x",start.field = "Start.x",end.field = "End.x"),BT245_loops_all))
BT245_gainedpeaks.annotation.rna$Parental_loop=sapply(extractList(BT245_loops$V4,hits),paste)

hits2=as.list(findOverlaps(makeGRangesFromDataFrame(BT245_gainedpeaks.annotation.rna,seqnames.field = "Chromosome.x",start.field = "Start.x",end.field = "End.x"),BT245_loops_c2_all))
BT245_gainedpeaks.annotation.rna$KO_loop=sapply(extractList(BT245_c2_loops$V4,hits2),paste)
#add cbx2
b.cbx2=import.bed("BT245_CBX2_peaks_seqmonk.bed")
hits=as.list(findOverlaps(makeGRangesFromDataFrame(BT245_gainedpeaks.annotation.rna,seqnames.field = "Chromosome.x",start.field = "Start.x",end.field = "End.x"),b.cbx2))
b.cbx2=as.data.frame(b.cbx2)
b.cbx2$cbx2="CBX2"
BT245_gainedpeaks.annotation.rna$CBX2=sapply(extractList(b.cbx2$cbx2,hits),paste)

#differential loops
BT245_maintainedloops=intersect(BT245_loops_all,BT245_loops_c2_all)
BT245_gainedloops=setdiff(BT245_loops_all,subsetByOverlaps(BT245_loops_all,BT245_maintainedloops))

hits3=as.list(findOverlaps(makeGRangesFromDataFrame(BT245_gainedpeaks.annotation.rna,seqnames.field = "Chromosome.x",start.field = "Start.x",end.field = "End.x"),BT245_gainedloops))
BT245_gainedloops=as.data.frame(BT245_gainedloops)
BT245_gainedloops$probe=rownames(BT245_gainedloops)
BT245_gainedpeaks.annotation.rna$gained_loop=sapply(extractList(BT245_gainedloops$probe,hits),paste)

#plot
names(BT245_gainedpeaks.annotation.rna)
res=BT245_gainedpeaks.annotation.rna
res$Log2_FC=-res$Log2_FC
select.lab=res[res$K27ac=="character(0)",]
select.lab=select.lab[!select.lab$K27me3=="character(0)",]
select.lab=select.lab[abs(select.lab$Log2_FC)>2,]
select.lab=unique(select.lab$Probe[select.lab$P.adj<10^-5])

jpeg("Volcanoplot_superrecruitment_theory_BT245.jpeg",res=600,width = 10, height = 8, units = 'in')
EnhancedVolcano(test,
                lab = test$Probe,
                x = 'Log2_FC',
                y = 'P.adj',
                selectLab = select.lab,
                xlim = c(-12,12),
                title="BT245 K27M K27me3 superrecruitment sites - K27M vs K27M KO",
                xlab = bquote(~Log[2]~ 'fold change'),
                ylab=bquote(~-Log[10]~adjusted~italic(P)),
                pointSize = 1.0,
                labSize = 1.5,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = FALSE,
                colAlpha = 1)
dev.off()



p1 <- EnhancedVolcano(test,
                      lab = test$Probe,
                      x = 'Log2_FC',
                      y = 'P.adj',
                      selectLab = select.lab,
                      xlim = c(-12,12),
                      xlab = bquote(~Log[2]~ 'fold change'),
                      title = 'K27me3 superrecruitment',
                      labSize = 3.0,
                      labCol = 'purple',
                      labFace = 'bold',
                      boxedLabels = TRUE,
                      shape = 42,
                      colCustom = c("blue"),
                      colAlpha = 1,
                      legendPosition = 'top',
                      legendLabSize = 15,
                      legendIconSize = 5.0,
                      shade = select.lab,
                      shadeLabel = 'Cell-type I',
                      shadeAlpha = 1/2,
                      shadeFill = 'purple',
                      shadeSize = 1,
                      shadeBins = 5,
                      drawConnectors = TRUE,
                      widthConnectors = 1.0,
                      colConnectors = 'grey30',
                      gridlines.major = TRUE,
                      gridlines.minor = FALSE,
                      border = 'partial',
                      borderWidth = 1.5,
                      borderColour = 'black')

p1


#how many genes were gained ring1b associated with?810 genes DEG out of (logfold1, pval of 0.1) 3091 genes.
#how many of these genes were upregulated/activated: 585 genes were activated .310 had k27me3 reduction
#how many of these genes were repressed: 225 genes
#how many were associated with k27me3 superrecruitment: 407 genes out of 810 genes 
#how many were associated with loops: 226 genes (100 with K27me3 superrecruitment,34 with acetylation,104 with neither)
#207 genes with loops: 42 were downregulated, 165 were upregulated
#how many were associated with k27ac increase: 120 genes
#out of all ectopic genes,1288 were associated with K27me3, 725 were associated with K27ac

test=BT245_gainedpeaks.annotation.rna[BT245_gainedpeaks.annotation.rna$K27ac=="character(0)",]
test=test[!test$K27me3=="character(0)",]
test=test[!duplicated(test$Probe),]
test=BT245_gainedpeaks.annotation.rna[BT245_gainedpeaks.annotation.rna$P.adj<c(0.1),]
test=test[abs(test$Log2_FC)>1,]
test=test[test$K27me3=="character(0)",]
test=test[!test$K27ac=="character(0)",]
test=test[test$CBX2=="character(0)",]
test=test[!test$Parental_loop=="character(0)",]
test=test[!duplicated(test$Probe),]
test=makeGRangesFromDataFrame(test,seqnames.field = "Chromosome.x",start.field = "Start.x",end.field = "End.x")
export.bed(test,"test.bed")

#create piechart of gained ring1b associations
# Load ggplot2
library(ggplot2)
library("ggsci")
library("ggplot2")
library("gridExtra")
# Create Data
data <- data.frame(
  a=c("H3K27me3","H3K27ac","Others"),
  value=c(407,120,283)
)

# Compute the position of labels
data <- data %>% 
  arrange(desc(a)) %>%
  mutate(prop = value / sum(data$value) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

# Basic piechart
ggplot(data, aes(x="", y=prop, fill=a)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  theme(legend.position="none") +
  ggtitle("K27M-specific RING1b peaks in BT245")+
  
  geom_text(aes(y = ypos, label = a), color = "white", size=6) +
  scale_fill_startrek()
ggsave("PieChart_BT245_gainedRING1b_association.tiff", units="in", width=5, height=4, dpi=600, compression = 'lzw')

export.bed(BT245_gainedpeaks,"test.bed")
#Make a Grange of these RING1b superrecruitment sites
supertable1=BT245_gainedpeaks.annotation.rna[BT245_gainedpeaks.annotation.rna$K27ac=="character(0)",]
supertable1=supertable1[!supertable1$K27me3=="character(0)",]
supertable1=supertable1[supertable1$P.adj<0.1,]
supertable1=supertable1[abs(supertable1$Log2_FC)>1,]

#how many of these loops had k27ac
hits=as.list(findOverlaps(BT245_loops_all.1,makeGRangesFromDataFrame(BT245_K27ac_annotation)))
BT245_loops_all.1$k27ac=sapply(extractList(BT245_K27ac_annotation$K27ac,hits),paste)
BT245_loops_all.1=as.data.frame(BT245_loops_all.1)
BT245_loops_all.2=BT245_loops_all.1[!BT245_loops_all.1$k27ac=="character(0)",]
BT245_loops_all.2=makeGRangesFromDataFrame(BT245_loops_all.2)

hits=as.list(findOverlaps(makeGRangesFromDataFrame(supertable1,seqnames.field = "Chromosome.x",start.field = "Start.x",end.field = "End.x"),BT245_loops_all))
BT245_loops_all=as.data.frame(BT245_loops_all)
BT245_loops_all$loop="loop"
supertable1$loop=sapply(extractList(BT245_loops_all$loop,hits),paste)
supertable2=supertable1[!supertable1$loop=="character(0)",]

#plot heatmap
library(pheatmap)
library(dplyr)
setwd("C:/Users/egj19/OneDrive - McGill University/Desktop/Analysis/peakcalls")
BT245_deg=read.table("DESeq2_BT245.txt",header=T)
test=DIPG13_gained_cprc1
test=test[,c(1:10,16:19,11:15,20:21)]
#test=left_join(test,BT245_deg,by="Probe")
#test=test[order(-abs(test$Log2_FC)),]
test$diff=test$Parental-test$KO
test=test[order(-test$diff),]
test=test[1:15,]
rownames(test)=test$Probe
names(test)

test=test[,c(11:19)]

#names(test)=paste(c("BT245 K27M KO 1","BT245 K27M KO 2","BT245 K27M KO 3","BT245 K27M KO 4","BT245 K27M KO 5","BT245 Parental 1","BT245 Parental 2","BT245 Parental 3","BT245 Parental 4","BT245 Parental 5"))
names(test)=paste(c("DIPGXIII K27M KO 1","DIPGXIII K27M KO 2","DIPGXIII K27M KO 3","DIPGXIII K27M KO 4","DIPGXIII Parental 1","DIPGXIII Parental 2","DIPGXIII Parental 3","DIPGXIII Parental 4","DIPGXIII Parental 5"))
test=as.matrix(test)


cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}
data_subset_norm <- t(apply(test, 1, cal_z_score))

jpeg("Heatmap_TOP15genes_DIPG13_RING1b_cPRC1.jpeg",res=1200,width = 10, height = 8, units = 'in')

pheatmap(test,scale = "row",fontsize_row = 10,cluster_cols = FALSE,cluster_rows = FALSE,main="Top 15 DEGs associated with RING1b super-recruitment to CGI promoters")

dev.off()
#how many are up/down regulated
#251 are downregulated 443 are upregulated
#how many had loops associated with them? 155 had loops. 53 were downregulated 96 were upregulated
#how many ring1b sites have loops associated
a=BT245_gainedpeaks.annotation.rna[!BT245_gainedpeaks.annotation.rna$loops=="character(0)",]
length(subsetByOverlaps(BT245,makeGRangesFromDataFrame(BT245_loops,seqnames.field= "V1",start.field = "V2",end.field = "V3")))/length(BT245)


#Chipseeker for these k27me3 superrecruitment sites with deg
b=makeGRangesFromDataFrame(a,seqnames.field = "Chromosome.x",start.field = "Start.x",end.field = "End.x")
library(ChIPseeker)
library(org.Hs.eg.db)
txdb=TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
peakAnno <- annotatePeak(b, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")
gene <- seq2gene(b, tssRegion = c(-3000, 3000), flankDistance = 3000, TxDb=txdb)
pathway2 <- enrichPathway(gene)
head(pathway2, 2)
dotplot(pathway1)
library(ReactomePA)
pathway1 <- enrichPathway(as.data.frame(peakAnno)$geneId)
dotplot(pathway1)
head(pathway1, 2)

#optimize gained peaks call
library(rtracklayer)
BT245_gainedpeaks=import.bed("BT245_gainedpeaks.bed")
setwd("C:/Users/egj19/OneDrive - McGill University/Desktop/Analysis/peakcalls")
data=read.table("seqmonk25.txt",header=T)
data$RING1b_BT245=data$BT245.DMSO.2.Rx_XChIP_RING1B_condition3.sorted.dup.bam/data$BT245.DMSO.3.Rx_cells_ChIP1_Input_1.sorted.dup.bam
data1=read.table("seqmonk27.txt",header=T)
data1$RING1b_BT245_KO=data1$BT245.ko.C2P8.Rx_XChIP_RING1B.sorted.dup.bam/data1$BT245.ko.c2p3.Rx_cells_ChIP1_Input_1.sorted.dup.bam
data1=data1[!is.na(data1$RING1b_BT245_KO),]
data1=data1[!is.infinite(data1$RING1b_BT245_KO),]

library(dplyr)
test3=left_join(data,data1,by="Probe")
test3=test3[complete.cases(test3),]
test3=test3[!test3$RING1b_BT245_KO==0,]
test3$diff=log2(test3$RING1b_BT245/test3$RING1b_BT245_KO)
summary(test3$diff)

test3=test3[test3$diff>0.9,]

test4=makeGRangesFromDataFrame(test3,seqnames.field = "Chromosome.x",start.field = "Start.x",end.field = "End.x")
length(test4)

test5=subsetByOverlaps(BT245_gainedpeaks,test4)
test5
test7=makeGRangesFromDataFrame(test3[test3$diff>3,],seqnames.field = "Chromosome.x",start.field = "Start.x",end.field = "End.x")

test6=union(test5,test7)
test6

export.bed(test6,"test.bed")

#import ring1b regions on promoters with k27me3 from chase
setwd("C:/Users/egj19/OneDrive - McGill University/Desktop")
test=import.gff("test.gff")
test=makeGRangesFromDataFrame(test)
hits=as.list(findOverlaps(test,test4))
test=as.data.frame(test)
test$diff=sapply(extractList(as.numeric(test3$diff),hits),mean)
test=test[!test$diff=="numeric(0)",]
summary(test$diff)


#optimize peakcalls
setwd("C:/Users/egj19/OneDrive - McGill University/Desktop/Analysis/peakcalls")
DIPG13=import.bed("DIPG13_C5P6_RING1b_peaks_seqmonk.bed")
data=read.table("seqmonk31.txt",header=T)
data$RING1b_DIPG13=data$DIPGXIII.ko.c5p6.7.Rx_XChIP_RING1B.sorted.dup.bam/data$DIPGXIII.ko.c5p6.7.Rx_XChIP_Input.sorted.dup.bam
data=data[!is.na(data$RING1b_DIPG13),]
data=data[!is.infinite(data$RING1b_DIPG13),]
test1=data[data$RING1b_DIPG13>2.3,]

test1=makeGRangesFromDataFrame(test1)
length(test1)

test2=subsetByOverlaps(DIPG13,test1)
length(test2)

test3=makeGRangesFromDataFrame(data[data$RING1b_DIPG13>2.9,])
test4=makeGRangesFromDataFrame(data[data$DIPGXIII.ko.c5p6.7.Rx_XChIP_RING1B.sorted.dup.bam>2.9,])
summary(data$DIPGXIII.ko.c5p6.7.Rx_XChIP_RING1B.sorted.dup.bam)
test5=subsetByOverlaps(test4,test3)
length(test5)
test6=setdiff(test5,subsetByOverlaps(test5,test2))
length(test6)

test7=c(test6,test2)
test7=makeGRangesFromDataFrame(test7)
length(test7)
export.bed(test7,"DIPG13_C5P6_RING1b_peaks_optimized.bed")



