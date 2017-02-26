## loading packages
require('ChIPseeker')
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
require('clusterProfiler')
require('rtracklayer')
files <- list('/Volumes/Hard Drive/human_naive/GSE64758_RAW/GSM1579346_Betacatenin_Wnt3a.ucsc.bigWig')
names(files) = 'WNT'
peakAnno <- annotatePeak(files[[1]], tssRegion=c(-3000, 3000), 
                         TxDb=txdb, annoDb="org.Hs.eg.db",verbose = TRUE)




library(BiocParallel)

eByg <- exonsBy(txdb, by="gene")
??bfl <- BamFileList(outpaths(args), yieldSize=50000, index=character())
multicoreParam <- MulticoreParam(workers=4); register(multicoreParam); registered()
??counteByg <- bplapply(bfl, function(x) summarizeOverlaps(eByg, x, mode="Union", ignore.strand=TRUE, inter.feature=TRUE, singleEnd=FALSE)) # Note: for strand-specific RNA-Seq set 'ignore.strand=FALSE' and for PE data set 'singleEnd=FALSE'

countDFeByg <- sapply(seq(along=counteByg), function(x) assays(counteByg[[x]])$counts)

rownames(countDFeByg) <- names(rowRanges(counteByg[[1]])); colnames(countDFeByg) <- names(bfl)
rpkmDFeByg <- apply(countDFeByg, 2, function(x) returnRPKM(counts=x, ranges=eByg))
??write.table(countDFeByg, "results/countDFeByg.xls", col.names=NA, quote=FALSE, sep="\t")
??write.table(rpkmDFeByg, "results/rpkmDFeByg.xls", col.names=NA, quote=FALSE, sep="\t")