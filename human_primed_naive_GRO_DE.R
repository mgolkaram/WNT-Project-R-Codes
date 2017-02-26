require('Sushi')
require('rtracklayer')
require('ChIPseeker')
require('GenomicRanges')

#GRO.naive.df = read.csv('/Volumes/Hard Drive/human_naive/GRO(RNAseq)/GSM1579369_Wnt3a_1.ucsc.bigWig')
#GRO.PET.primed.df = read.csv('/Volumes/Hard Drive/human_naive/GRO(RNAseq)/GSM1579367_untreated_1.ucsc.bigWig')
GRO.naive.gr = import('/Volumes/Hard Drive/human_naive/GRO(RNAseq)/GSM1579369_Wnt3a_1.ucsc.bigWig',format = 'bigWig')
GRO.primed.gr = import('/Volumes/Hard Drive/human_naive/GRO(RNAseq)/GSM1579367_untreated_1.ucsc.bigWig',format = 'bigWig')

gr <- ChIP.seq.naive.gr
ChIP.seq.naive.df <- data.frame(chr=seqnames(gr),
                                starts=start(gr)-1,
                                ends=end(gr),
                                scores=gr$score,
                                strands=strand(gr))

gr <- ChIP.seq.primed.gr
ChIP.seq.primed.df <- data.frame(chr=seqnames(gr),
                                 starts=start(gr)-1,
                                 ends=end(gr),
                                 scores=gr$score,
                                 strands=strand(gr))


ChIA.PET.naive.df1 = ChIA.PET.naive.df[,c(1,2,3,8)]
colnames(ChIA.PET.naive.df1) = c('chr','start','end','score')

ChIA.PET.naive.df2 = ChIA.PET.naive.df[,c(4,5,6,8)]
colnames(ChIA.PET.naive.df2) = c('chr','start','end','score')

ChIA.PET.naive.df.total = rbind(ChIA.PET.naive.df1,ChIA.PET.naive.df2)

ChIA.PET.primed.df1 = ChIA.PET.primed.df[,c(1,2,3,8)]
colnames(ChIA.PET.primed.df1) = c('chr','start','end','score')

ChIA.PET.primed.df2 = ChIA.PET.primed.df[,c(4,5,6,8)]
colnames(ChIA.PET.primed.df2) = c('chr','start','end','score')

ChIA.PET.primed.df.total = rbind(ChIA.PET.primed.df1,ChIA.PET.primed.df2)

ChIA.PET.naive.gr = makeGRangesFromDataFrame(ChIA.PET.naive.df.total)
ChIA.PET.primed.gr = makeGRangesFromDataFrame(ChIA.PET.primed.df.total)


#### KLF5
chrom            = "chr13"
chromstart       = 73503142
chromend         = 73751676

plot(new = TRUE)
plotBedpe(ChIA.PET.naive.df,chrom,chromstart,chromend,color = 'red',heights = ChIA.PET.naive.df$Score,offset=0,flip=FALSE,bty='n',
          lwd=1,plottype="ribbons",border="red",transparency=.50)
plotBedgraph(ChIP.seq.naive.df,chrom,chromstart,chromend,transparency=.50,flip=FALSE,color="blue",linecol="blue",overlay=TRUE,rescaleoverlay=TRUE)
labelgenome(chrom,chromstart,chromend,side=1,scipen=20,n=3,line=.18,chromline=.5,scaleline=0.5,scale="Mb")
#axis(side=2,las=2,tcl=.2)
#mtext("Z-score",side=2,line=1.75,cex=.75,font=2)
title(main = 'KLF5')
arrows(x0 = 73633142,x1 = 73651676, y0 = 3, y1= 3,length = 0.1,angle = 10,lwd = 2)
lines(c(73633142,73651676),c(3,3),lwd = 2,pch = '|',type = 'p')
legend("topleft",inset =0.01,legend=c(expression(paste(beta,'-cat (Wnt3a)')),'', "CTCF-CTCF loops"),col=c('blue','white','red'),pch=19,bty='n',text.font=2)

plot(new = TRUE)
plotBedpe(ChIA.PET.primed.df,chrom,chromstart,chromend,color = 'red',heights = ChIA.PET.primed.df$Score,offset=0,flip=FALSE,bty='n',
          lwd=1,plottype="ribbons",border="red",transparency=.50)
plotBedgraph(ChIP.seq.primed.df,chrom,chromstart,chromend,transparency=.50,flip=FALSE,color="blue",linecol="blue",overlay=FALSE,rescaleoverlay=FALSE,range = c(0,200))
labelgenome(chrom,chromstart,chromend,side=1,scipen=20,n=3,line=.18,chromline=.5,scaleline=0.5,scale="Mb")
#axis(side=2,las=2,tcl=.2)
#mtext("Z-score",side=2,line=1.75,cex=.75,font=2)
title(main = 'KLF5')
arrows(x0 = 73633142,x1 = 73651676, y0 = 100, y1= 100,length = 0.1,angle = 10,lwd = 2)
lines(c(73633142,73651676),c(100,100),lwd = 2,pch = '|',type = 'p')
legend("topleft",inset =0.01,legend=c(expression(paste(beta,'-cat (Untreated)')),'', "CTCF-CTCF loops"),col=c('blue','white','red'),pch=19,bty='n',text.font=2)


#### LEFTY1
chrom            = "chr1"
chromstart       = 226053982
chromend         = 226096846

plot(new = TRUE)
plotBedpe(ChIA.PET.naive.df,chrom,chromstart,chromend,color = 'red',heights = ChIA.PET.naive.df$Score,offset=0,flip=FALSE,bty='n',
          lwd=1,plottype="ribbons",border="red",transparency=.50)
plotBedgraph(ChIP.seq.naive.df,chrom,chromstart,chromend,transparency=.50,flip=FALSE,color="blue",linecol="blue",overlay=TRUE,rescaleoverlay=TRUE)
labelgenome(chrom,chromstart,chromend,side=1,scipen=20,n=3,line=.18,chromline=.5,scaleline=0.5,scale="Mb")
#axis(side=2,las=2,tcl=.2)
#mtext("Z-score",side=2,line=1.75,cex=.75,font=2)
title(main = 'LEFTY1')
arrows(x0 = 226076846,x1 = 226073982, y0 = 3, y1= 3,length = 0.1,angle = 10,lwd = 2)
lines(c(226076846,226073982),c(3,3),lwd = 2,pch = '|',type = 'p')
legend("topleft",inset =0.01,legend=c(expression(paste(beta,'-cat (Wnt3a)')),'', "CTCF-CTCF loops"),col=c('blue','white','red'),pch=19,bty='n',text.font=2)

plot(new = TRUE)
plotBedpe(ChIA.PET.primed.df,chrom,chromstart,chromend,color = 'red',heights = ChIA.PET.primed.df$Score,offset=0,flip=FALSE,bty='n',
          lwd=1,plottype="ribbons",border="red",transparency=.50)
plotBedgraph(ChIP.seq.primed.df,chrom,chromstart,chromend,transparency=.50,flip=FALSE,color="blue",linecol="blue",overlay=FALSE,rescaleoverlay=FALSE,range = c(0,200))
labelgenome(chrom,chromstart,chromend,side=1,scipen=20,n=3,line=.18,chromline=.5,scaleline=0.5,scale="Mb")
#axis(side=2,las=2,tcl=.2)
#mtext("Z-score",side=2,line=1.75,cex=.75,font=2)
title(main = 'LEFTY1')
arrows(x0 = 226076846,x1 = 226073982, y0 = 100, y1= 100,length = 0.1,angle = 10,lwd = 2)
lines(c(226076846,226073982),c(100,100),lwd = 2,pch = '|',type = 'p')
legend("topleft",inset =0.01,legend=c(expression(paste(beta,'-cat (Untreated)')),'', "CTCF-CTCF loops"),col=c('blue','white','red'),pch=19,bty='n',text.font=2)

#### PRDM14
chrom            = "chr8"
chromstart       = 70953886
chromend         = 71023562

plot(new = TRUE)
plotBedpe(ChIA.PET.naive.df,chrom,chromstart,chromend,color = 'red',heights = ChIA.PET.naive.df$Score,offset=0,flip=FALSE,bty='n',
          lwd=1,plottype="ribbons",border="red",transparency=.50)
plotBedgraph(ChIP.seq.naive.df,chrom,chromstart,chromend,transparency=.50,flip=FALSE,color="blue",linecol="blue",overlay=TRUE,rescaleoverlay=TRUE)
labelgenome(chrom,chromstart,chromend,side=1,scipen=20,n=3,line=.18,chromline=.5,scaleline=0.5,scale="Mb")
#axis(side=2,las=2,tcl=.2)
#mtext("Z-score",side=2,line=1.75,cex=.75,font=2)
title(main = 'PRDM14')
arrows(x0 = 70983562,x1 = 70963886, y0 = 2.3, y1= 2.3,length = 0.1,angle = 10,lwd = 2)
lines(c(70983562,70963886),c(2.3,2.3),lwd = 2,pch = '|',type = 'p')
legend("topleft",inset =0.01,legend=c(expression(paste(beta,'-cat (Wnt3a)')),'', "CTCF-CTCF loops"),col=c('blue','white','red'),pch=19,bty='n',text.font=2)

plot(new = TRUE)
plotBedpe(ChIA.PET.primed.df,chrom,chromstart,chromend,color = 'red',heights = ChIA.PET.primed.df$Score,offset=0,flip=FALSE,bty='n',
          lwd=1,plottype="ribbons",border="red",transparency=.50)
plotBedgraph(ChIP.seq.primed.df,chrom,chromstart,chromend,transparency=.50,flip=FALSE,color="blue",linecol="blue",overlay=FALSE,rescaleoverlay=FALSE,range = c(0,200))
labelgenome(chrom,chromstart,chromend,side=1,scipen=20,n=3,line=.18,chromline=.5,scaleline=0.5,scale="Mb")
#axis(side=2,las=2,tcl=.2)
#mtext("Z-score",side=2,line=1.75,cex=.75,font=2)
title(main = 'PRDM14')
arrows(x0 = 70983562,x1 = 70963886, y0 = 100, y1= 100,length = 0.1,angle = 10,lwd = 2)
lines(c(70983562,70963886),c(100,100),lwd = 2,pch = '|',type = 'p')
legend("topleft",inset =0.01,legend=c(expression(paste(beta,'-cat (Untreated)')),'', "CTCF-CTCF loops"),col=c('blue','white','red'),pch=19,bty='n',text.font=2)

#### Overlap test
require(ChIPseeker)
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
require(clusterProfiler)

write.table(ChIP.seq.naive.df, file="/Volumes/Hard Drive/ChIA-PET/ChIP_seq_naive.bed", quote=F, sep="\t", row.names=F, col.names=F)
write.table(ChIA.PET.naive.df.total, file="/Volumes/Hard Drive/ChIA-PET/ChIA_PET_naive.bed", quote=F, sep="\t", row.names=F, col.names=F)
write.table(ChIP.seq.primed.df, file="/Volumes/Hard Drive/ChIA-PET/ChIP_seq_primed.bed", quote=F, sep="\t", row.names=F, col.names=F)
write.table(ChIA.PET.primed.df.total, file="/Volumes/Hard Drive/ChIA-PET/ChIA_PET_primed.bed", quote=F, sep="\t", row.names=F, col.names=F)

myFiles.naive = list('/Volumes/Hard Drive/ChIA-PET/ChIP_seq_naive.bed.gz','/Volumes/Hard Drive/ChIA-PET/ChIA_PET_naive.bed.gz')
names(myFiles.naive) = c('ChIP.seq','ChIA.PET')
enrichPeakOverlap(queryPeak     = myFiles.naive[[1]], 
                  targetPeak    = myFiles.naive[[2]], 
                  TxDb          = txdb, 
                  pAdjustMethod = "BH", 
                  nShuffle      = 50, 
                  chainFile     = NULL,
                  verbose       = TRUE)

myFiles.primed = list('/Volumes/Hard Drive/ChIA-PET/ChIP_seq_primed.bed.gz','/Volumes/Hard Drive/ChIA-PET/ChIA_PET_primed.bed.gz')
names(myFiles.primed) = c('ChIP.seq','ChIA.PET')
enrichPeakOverlap(queryPeak     = myFiles.primed[[1]], 
                  targetPeak    = myFiles.primed[[2]], 
                  TxDb          = txdb, 
                  pAdjustMethod = "BH", 
                  nShuffle      = 50, 
                  chainFile     = NULL,
                  verbose       = TRUE)

### Only significant peaks 

ChIP.seq.primed.df.sig = ChIP.seq.primed.df[which(ChIP.seq.primed.df$scores>100),]
ChIA.PET.primed.df.total.sig = ChIA.PET.primed.df.total[which(ChIA.PET.primed.df.total$score>100),]

write.table(ChIP.seq.primed.df.sig, file="/Volumes/Hard Drive/ChIA-PET/ChIP_seq_primed_sig.bed", quote=F, sep="\t", row.names=F, col.names=F)
write.table(ChIA.PET.primed.df.total.sig, file="/Volumes/Hard Drive/ChIA-PET/ChIA_PET_primed_sig.bed", quote=F, sep="\t", row.names=F, col.names=F)

myFiles.primed.sig = list('/Volumes/Hard Drive/ChIA-PET/ChIP_seq_primed_sig.bed.gz','/Volumes/Hard Drive/ChIA-PET/ChIA_PET_primed_sig.bed.gz')
names(myFiles.primed.sig) = c('ChIP.seq','ChIA.PET')
enrichPeakOverlap(queryPeak     = myFiles.primed.sig[[1]], 
                  targetPeak    = myFiles.primed[[2]], 
                  TxDb          = txdb, 
                  pAdjustMethod = "BH", 
                  nShuffle      = 1000, 
                  chainFile     = NULL,
                  verbose       = TRUE)


ChIP.seq.naive.df.sig = ChIP.seq.naive.df[which(ChIP.seq.naive.df$scores>100),]
ChIA.PET.naive.df.total.sig = ChIA.PET.naive.df.total[which(ChIA.PET.naive.df.total$score>100),]

write.table(ChIP.seq.naive.df.sig, file="/Volumes/Hard Drive/ChIA-PET/ChIP_seq_naive_sig.bed", quote=F, sep="\t", row.names=F, col.names=F)
write.table(ChIA.PET.naive.df.total.sig, file="/Volumes/Hard Drive/ChIA-PET/ChIA_PET_naive_sig.bed", quote=F, sep="\t", row.names=F, col.names=F)

myFiles.naive.sig = list('/Volumes/Hard Drive/ChIA-PET/ChIP_seq_naive_sig.bed.gz','/Volumes/Hard Drive/ChIA-PET/ChIA_PET_naive_sig.bed.gz')
names(myFiles.naive.sig) = c('ChIP.seq','ChIA.PET')
enrichPeakOverlap(queryPeak     = myFiles.naive.sig[[1]], 
                  targetPeak    = myFiles.naive[[2]], 
                  TxDb          = txdb, 
                  pAdjustMethod = "BH", 
                  nShuffle      = 1000, 
                  chainFile     = NULL,
                  verbose       = TRUE)



#####################################################################
gr1<-ChIP.seq.naive.gr[which(ChIP.seq.naive.gr$score>10)]
gr2<-ChIP.seq.primed.gr[which(ChIP.seq.primed.gr$score>10)]
gr3<-subsetByOverlaps(gr1,gr2)

new.gr = gr1[!gr1%in%gr3]

gr4<-subsetByOverlaps(ChIA.PET.primed.gr,ChIA.PET.naive.gr)
gr5<-ChIA.PET.primed.gr[!ChIA.PET.primed.gr%in%gr4]



Beta.Cat<-data.frame(chr=seqnames(new.gr),
                     starts=start(new.gr)-1,
                     ends=end(new.gr),
                     scores=new.gr$score,
                     strands=strand(new.gr))

loop<-data.frame(chr=seqnames(gr5),
                 starts=start(gr5)-1,
                 ends=end(gr5),
                 strands=strand(gr5))

write.table(Beta.Cat, file="/Volumes/Hard Drive/ChIA-PET/Beta_Cat.bed", quote=F, sep="\t", row.names=F, col.names=F)
write.table(loop, file="/Volumes/Hard Drive/ChIA-PET/loop.bed", quote=F, sep="\t", row.names=F, col.names=F)

myFiles = list('/Volumes/Hard Drive/ChIA-PET/Beta_Cat.bed.gz','/Volumes/Hard Drive/ChIA-PET/loop.bed.gz')
names(myFiles) = c('Beta.Cat','loop')
enrichPeakOverlap(queryPeak     = myFiles[[1]], 
                  targetPeak    = myFiles[[2]], 
                  TxDb          = txdb, 
                  pAdjustMethod = "BH", 
                  nShuffle      = 1000, 
                  chainFile     = NULL,
                  verbose       = TRUE)


