tempData <- read.table('/Volumes/Hard Drive/ChIA-PET/scRNA-seq/count.csv', sep = "")
LIF.2i.rep1 = 1:82
LIF.2i.rep2 = 83:141
LIF.2i.rep3 = 142:213
LIF.2i.rep4 = 214:295
LIF.2i = 1:295
LIF.2ai.rep1 = 296:388
LIF.2ai.rep2 = 389:454
serum.rep1 = 455:535
serum.rep2 = 536:626
serum.rep3 = 627:704
serum = 455:704
tempData = as.matrix(tempData)
idx = 1720
Data = tempData[,c(LIF.2i,serum)]
colnames(Data) = c(rep('2i',length(LIF.2i)),rep('serum',length(serum)))
# RNA-Seq
library(emdbook)

mESC.2i <- Data[,which(colnames(Data)=='2i')]
mESC.serum <- Data[,which(colnames(Data)=='serum')]

mESC.2i.SUM <- apply(as.matrix(mESC.2i),sum,MARGIN = 2)
mESC.serum.SUM <- apply(as.matrix(mESC.serum),sum,MARGIN = 2)

CPM.2i <- t(t(mESC.2i) / mESC.2i.SUM) *1e6
CPM.serum <- t(t(mESC.serum) / mESC.serum.SUM) *1e6

mESC.2i.MEAN <- apply(CPM.2i,MARGIN = 1,mean)
mESC.2i.SD <- apply(CPM.2i,MARGIN = 1,sd)
mESC.2i.CV <- mESC.2i.SD/mESC.2i.MEAN

mESC.serum.MEAN <- apply(CPM.serum,MARGIN = 1,mean)
mESC.serum.SD <- apply(CPM.serum,MARGIN = 1,sd)
mESC.serum.CV <- mESC.serum.SD/mESC.serum.MEAN


plot(mESC.2i.MEAN,mESC.2i.CV,log = 'x')
plot(mESC.serum.MEAN,mESC.serum.CV,log = 'x')

bins = lseq(10,1e4,100)

Fun <- function(x, argv){ 
  idx = data.frame(matrix(ncol = length(argv)-1,nrow = 20000))
  for (i in seq(58)){
    temp = which(x>=argv[i] & x<argv[i+1]) 
    idx[1:length(temp),i] = temp
  }
  return(idx)}

mESC.2i.idx <- Fun(mESC.2i.MEAN,bins)
mESC.2i.CV.MAT = apply(mESC.2i.idx,c(1,2),FUN = function(x, V) V[x] , V = mESC.2i.CV)
V1 = apply(mESC.2i.CV.MAT,MARGIN = 2, sd,na.rm = TRUE)
U1 = apply(mESC.2i.CV.MAT,MARGIN = 2, median,na.rm = TRUE)
mESC.2i.selected.CV = t((t(mESC.2i.CV.MAT)-U1) > 1.2*V1)
mESC.2i.selected.idx = subset.matrix(mESC.2i.idx,subset = mESC.2i.selected.CV)



mESC.serum.idx <- Fun(mESC.serum.MEAN,bins)
mESC.serum.CV.MAT = apply(mESC.serum.idx,c(1,2),FUN = function(x, V) V[x] , V = mESC.serum.CV)
V2 = apply(mESC.serum.CV.MAT,MARGIN = 2, sd,na.rm = TRUE)
U2 = apply(mESC.serum.CV.MAT,MARGIN = 2, median,na.rm = TRUE)
mESC.serum.selected.CV = t((t(mESC.serum.CV.MAT)-U2) > 1.2*V2)
mESC.serum.selected.idx = subset.matrix(mESC.serum.idx,subset = mESC.serum.selected.CV)


idx.filtered = sort(x = union(as.matrix(mESC.2i.selected.idx),as.matrix(mESC.serum.selected.idx)))
tempV2 = as.numeric(as.character(mESC.2i.CV[idx.filtered]))
tempV4 = as.numeric(as.character(mESC.serum.CV[idx.filtered]))
tempV1 = as.numeric(as.character(mESC.2i.MEAN[idx.filtered]))
tempV3 = as.numeric(as.character(mESC.serum.MEAN[idx.filtered]))
myData <- data.frame(tempV1,tempV2,tempV3,tempV4)
row.names(myData) = row.names(Data)[idx.filtered]

colnames(myData) = c('mESC.2i.MEAN','mESC.2i.CV','mESC.serum.MEAN','mESC.serum.CV')

temp = read.csv('/Volumes/Hard Drive/ChIA-PET/scRNA-seq/mESC_Den_CV.csv')
temp = temp[-which(duplicated(temp$Ensembl.ID)),]
#myGENE.IDs = temp$GeneID.1[which(as.character(temp$Ensembl.ID)%in%row.names(Data))]
lost.ids = which(!row.names(myData)%in%as.character(temp$Ensembl.ID))
#FC = mESC.2i.MEAN/mESC.serum.MEAN
myData = myData[-lost.ids,]
idx = which(as.character(temp$Ensembl.ID)%in%row.names(myData))
myData = myData[as.character(temp$Ensembl.ID[idx]),]
myData$GeneSymbol = temp$GeneID.1[idx]

myData$Chr = temp$Chr[idx]
myData$start = temp$start[idx]
myData$end = temp$end[idx]
myData$strand = temp$strand[idx]
myData$mESC.2i.Density = NA
myData$mESC.serum.Density = NA

beta.cat = read.csv('/Volumes/Hard Drive/ChIA-PET/ChIP-seq/betacat_bindingSite.csv')
beta.cat = toupper(beta.cat$X.4[-c(1,2)])

FC = abs(myData$mESC.2i.MEAN-myData$mESC.serum.MEAN)/pmin(myData$mESC.2i.MEAN,myData$mESC.serum.MEAN)
myData$FC = FC
select.idx = which(myData$FC<1 & as.character(myData$GeneSymbol)%in%beta.cat)
myData.filtered = myData[select.idx,]

boxplot(myData.filtered$mESC.2i.CV,myData.filtered$mESC.serum.CV,outline = FALSE)
boxplot(myData$mESC.2i.CV,myData$mESC.serum.CV,outline = FALSE)

plot(myData.filtered$mESC.2i.CV,myData.filtered$mESC.serum.CV, log = 'xy')



idx2 = which((myData.filtered$DM.2i>1) | (myData.filtered$DM.serum >1))
boxplot(myData.filtered$DM.2i[idx2],myData.filtered$DM.serum[idx2],outline = FALSE,names = c('2i','serum'),ylim = c(-1,3.5))
title('DM')
text(x = c(0,3.2),c('','***'),cex = 2)
t.test(myData.filtered$DM.2i[idx2],myData.filtered$DM.serum[idx2])
wilcox.test(myData.filtered$DM.2i[idx2],myData.filtered$DM.serum[idx2])









### Read DEs from file
non.DE.genes = read.csv('/Volumes/Hard Drive/ChIA-PET/RNAseq.csv')
non.DE.genes = non.DE.genes$Gene[which(non.DE.genes$Q.value>0.01)]

require(genefilter)
require(qvalue)
M = data.frame(CPM.2i,CPM.serum)
M = M[idx.filtered,]
factor = as.factor(c(rep('mESC.2i',dim(CPM.2i)[2]),rep('mESC.serum',dim(CPM.serum)[2])))
p.values<-rowttests(as.matrix(M),factor)$p.value
q.values<-qvalue(p.values)$qvalues
sum(q.values<0.01)
idx.filtered.ByDE = idx.filtered[q.values>0.01]
myGENE.IDs2 = temp$GeneID.1[which(as.character(temp$Ensembl.ID)%in%row.names(Data[idx.filtered.ByDE,]))]
lost.ids2 = which(!row.names(Data[idx.filtered.ByDE,])%in%as.character(temp$Ensembl.ID))
myData2 <- data.frame(myGENE.IDs2,tempV1[idx.filtered.ByDE],tempV2[idx.filtered.ByDE],tempV3[idx.filtered.ByDE],tempV4[idx.filtered.ByDE])
colnames(myData) = c('GeneID','mESC.2i.MEAN','mESC.2i.CV','mESC.serum.MEAN','mESC.serum.CV')


### filter DEs
myData.ByDE = myData[which(myData$GeneID%in%non.DE.genes),]

###### Some Plots
dev.off()
svg()
plot(myData.ByDE$mESC.2i.CV,myData.ByDE$mESC.serum.CV, col = 'blue', xlab = '2i', ylab = 'serum', main = 'CV During naive to primed')
abline(c(0,1))
Fit1 <- lm(myData.ByDE$hNPC.CV~myData.ByDE$hESC.CV)
abline(Fit1$coefficients,col = 'blue', lty = 2)
dev.off()
file.size("Rplot001.svg")/1000

plot(myData$hESC.CV,myData$hNPC.CV, col = 'red', xlab = 'hESC', ylab = 'NPC', main = 'CV During Differentiation')
abline(c(0,1))
Fit1 <- lm(myData.ByDE$hNPC.CV~myData.ByDE$hESC.CV)
abline(Fit1$coefficients,col = 'blue', lty = 2)
Fit2 <- lm(as.numeric(as.character(hNPC.CV[-idx.filtered.ByDE]))~as.numeric(as.character(H1.hESC.CV[-idx.filtered.ByDE])))
abline(Fit2$coefficients,col = 'red', lty = 2)
lines(myData.ByDE$hESC.CV,myData.ByDE$hNPC.CV,xlim = c(0,5),col = 'blue',type = 'p')
legend(0.5,5,legend = c('Differentially Expressed', 'House Keeping'), col = c('red','blue'), pch = 1,text.width = 2.3,cex = 0.8)

require(affy)
MAdat = cbind(myData$hESC.CV,myData$hNPC.CV)
colnames(MAdat) = c('hESC','NPS')
ma.plot( rowMeans(log2(MAdat)), log2(MAdat[, 2])-log2(MAdat[, 1]), cex=1)

require(gplots)
heatmap.2(MAdat[,c(2,1)],cexCol = 1.75)
#### TSS coordinates ####

coordinate.info = read.csv('Desktop/DATA.csv')
coordinate.info = coordinate.info[-which(coordinate.info$strand == 'n/a'),]
coordinate.info = coordinate.info[-which(duplicated(coordinate.info$GeneID)),]
idx = lapply(myData$GeneID,FUN = function(x) which(x == coordinate.info$GeneID))
idx[grep('integer', idx)] = NA
idx = unlist(idx)
myData$Chr = sapply(idx, FUN = function(x) coordinate.info$Chr[x])
myData$start = sapply(idx, FUN = function(x) coordinate.info$start[x])
myData$end = sapply(idx, FUN = function(x) coordinate.info$end[x])
myData$strand = sapply(idx, FUN = function(x) coordinate.info$strand[x])

#### DE Filtered TSS
coordinate.info = read.csv('Desktop/DATA.csv')
coordinate.info = coordinate.info[-which(coordinate.info$strand == 'n/a'),]
coordinate.info = coordinate.info[-which(duplicated(coordinate.info$GeneID)),]
idx = lapply(myData.ByDE$GeneID,FUN = function(x) which(x == coordinate.info$GeneID))
idx[grep('integer', idx)] = NA
idx = unlist(idx)
myData.ByDE$Chr = sapply(idx, FUN = function(x) coordinate.info$Chr[x])
myData.ByDE$start = sapply(idx, FUN = function(x) coordinate.info$start[x])
myData.ByDE$end = sapply(idx, FUN = function(x) coordinate.info$end[x])
myData.ByDE$strand = sapply(idx, FUN = function(x) coordinate.info$strand[x])




