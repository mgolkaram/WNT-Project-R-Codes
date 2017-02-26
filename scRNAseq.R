library(scater)
sce <- read.table('/Volumes/Hard Drive/ChIA-PET/scRNA-seq/count.csv', sep = "")
dim(sce)
checkForSpike <- grepl("^ERCC", rownames(sce))
checkForMito <- grepl("^mt-", rownames(sce))
sce <- calculateQCMetrics(sce, feature_controls=list(ERCC=checkForSpike, Mt=checkForMito))
head(colnames(pData(sce)))
library(scran)
isSpike(sce) <- "ERCC"

libsize.drop <- isOutlier(sce$total_counts, nmads=3, type="lower", log=TRUE)
feature.drop <- isOutlier(sce$total_features, nmads=3, type="lower", log=TRUE)

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

mito.drop <- isOutlier(sce$pct_counts_feature_controls_Mt, nmads=3, type="higher")
spike.drop <- isOutlier(sce$pct_counts_feature_controls_ERCC, nmads=3, type="higher")

sce <- sce[,!(libsize.drop | feature.drop | mito.drop | spike.drop)]
data.frame(ByLibSize=sum(libsize.drop), ByFeature=sum(feature.drop),
           ByMito=sum(mito.drop), BySpike=sum(spike.drop), Remaining=ncol(sce))

fontsize <- theme(axis.text=element_text(size=12), axis.title=element_text(size=16))
plotPCA(sce, pca_data_input="pdata") + fontsize

mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))
library(org.Mm.eg.db)
anno <- select(org.Mm.eg.db, column="ENSEMBL", keytype="SYMBOL", keys=rownames(sce))
ensembl <- anno$ENSEMBL[match(rownames(sce), anno$SYMBOL)]
assignments <- cyclone(sce, mm.pairs, gene.names=ensembl)
plot(assignments$score$G1, assignments$score$G2M, xlab="G1 score", ylab="G2/M score", pch=16)

sce <- sce[,assignments$phases=="G1"]

ave.counts <- rowMeans(counts(sce))
keep <- ave.counts >= 1
sum(keep)

dev.off()
svg()
require(sm)
GENENAME = myData$GeneID[order(myData$hESC.CV,decreasing = TRUE)[7:15]]
#par(new = TRUE)
par(mar=c(1,1,1,1))
par(mfrow = c(3,3))
lapply(X = GENENAME,FUN = function(arg){
  vioplot::vioplot(CPM.H1.hESC[Data$X==arg,],CPM.hNPC[Data$X==arg,], names = c('ESC','NPC'))
  title(arg, ylab = 'Expression Level')
})
par(new = FALSE)
dev.off()
file.size('Rplot001')/1000

numcells <- nexprs(sce, byrow=TRUE)
alt.keep <- numcells >= 10
sum(alt.keep)

sce <- computeSumFactors(sce, sizes=c(20, 40, 60, 80))
summary(sizeFactors(sce))

plot(sizeFactors(sce), sce$total_counts/1e6, log="xy",
     ylab="Library size (millions)", xlab="Size factor")

sce <- computeSpikeFactors(sce, type="ERCC", general.use=FALSE)

sce <- normalize(sce)


plotExplanatoryVariables(sce, variables=c("counts_feature_controls_ERCC", 
                                          "log10_counts_feature_controls_ERCC")) + fontsize

var.fit <- trendVar(sce, trend="loess", use.spikes=FALSE, span=0.2)
var.out <- decomposeVar(sce, var.fit)

plot(var.out$mean, var.out$total, pch=16, cex=0.6, xlab="Mean log-expression", 
     ylab="Variance of log-expression")

naive = read.csv('/Volumes/Hard Drive/ChIA-PET/ChIA-PET_Naive.csv')
primed = read.csv('/Volumes/Hard Drive/ChIA-PET/ChIA-PET_Primed.csv')
naive$sample = 1
primed$sample = 2

hg19 = read.csv('/Volumes/Hard Drive/ChIA-PET/GO.csv')

cross <- function(ID, data.set,SHIFT = 500000){
  chr = hg19$Chr[ID]
  x = ifelse(hg19$strand[ID] == '+',hg19$start[ID], hg19$end[ID])
  data.set = data.set[which(as.character(data.set$Chr_1) == chr),]
  idx1 = which((data.set$Start_1 - SHIFT)<x & (data.set$End_1 + SHIFT)>x)
  idx2 = which((data.set$Start_2 - SHIFT)<x & (data.set$End_2 + SHIFT)>x)
  idx = union(idx1, idx2)
  return(sum(data.set$Score[idx]))
  
}

naive.Density = unlist(lapply(X = seq(dim(hg19)[1]),FUN =  cross, data.set =  naive))
primed.Density = unlist(lapply(X = seq(dim(hg19)[1]), FUN = cross, data.set = primed))

boxplot(as.matrix(cbind(naive.Density,primed.Density)),outline = FALSE,names=c('Naive','Primed'), ylim = c(0,480),main = 'Density obtained by ChIA-PET Scores')
text(c(0,470),c('','***'),cex = 2)

o <- order(var.out$mean)
lines(var.out$mean[o], var.out$tech[o], col="green")
cur.spike <- isSpike(sce)
points(var.out$mean[cur.spike], var.out$total[cur.spike], col="red", pch=16)

