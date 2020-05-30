library(DiffBind)
#load data
setwd("C:/Users/egj19/OneDrive - McGill University/Desktop/Analysis/peakcalls")
samples <- read.csv('diffbind_ring1b.csv')
dbObj <- dba(sampleSheet=samples)
dbObj
dbObj <- dba.count(dbObj, bUseSummarizeOverlaps=FALSE)
dbObj
#plot PCA
dba.plotPCA(dbObj,  attributes=DBA_CONDITION, label=DBA_ID)
dev.off()
options(device = "RStudioGD")
par(mar=c(1,1,1,1))
plot(dbObj)

#establish contrast
dbObj <- dba.contrast(dbObj, categories=c(DBA_CONDITION), minMembers = 2)
dbObj
dbObj <- dba.analyze(dbObj, method=DBA_ALL_METHODS)
results <- dba.report(dbObj, contrast=1, th=1)
View(as.data.frame(results))
dba.show(dbObj, bContrasts=T)	
dba.plotMA(dbObj, method=DBA_EDGER)

#ven diagram
olap.rate <- dba.overlap(dbObj,mode=DBA_OLAP_ALL)

names(dbObj$masks)
dba.plotVenn(dbObj,dbObj$masks$Parental)
dbobj_consensus <- dba.peakset(dbObj, consensus=c(DBA_CONDITION),minOverlap=0.9)
dbobj_consensus
#jpeg("DiffBind_VennDiagram_RING1b.jpeg", res = 600,width = 10, height = 5, units = 'in')
dba.plotVenn(dbObj,dbObj$masks$H3K27me3 & dbObj$masks$RING1b)
#dev.off()
