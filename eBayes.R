require(limma)
Jiwon  = read.csv('Desktop/Jiwon_Cell_RNAseq.csv')
fac1 = factor(c('hESC','hESC','hESC','NPC','NPC','NPC'))
fac2 = factor(c('hESC','hESC','hESC','MES','MES','MES'))
fit1 = lmFit(Jiwon[-c(1,2,9:11)],design = model.matrix(~fac1))
fit1 <- eBayes(fit1)
fit2 = lmFit(Jiwon[-c(1,2,6:8)],design = model.matrix(~fac2))
fit2 <- eBayes(fit2)
#tt <- topTable(fit, coef=2, number=Inf, sort.by="none")
#Jiwon.RNAseq.geneNames = Jiwon[which(tt$adj.P.Val > 0.05),]
res1 <- topTable(fit1,coef=2, confint=0.95,number = Inf,sort.by = 'none') # 95% CIs for the log-fold changes
res2 <- topTable(fit2,coef=2, confint=0.95,number = Inf,sort.by = 'none')
keep1 <- which(res1$CI.L > -1 & res1$CI.R < 1)
keep2 <- which(res2$CI.L > -1 & res2$CI.R < 1)
keep = intersect(keep1,keep2)
res1[keep,]
Jiwon.RNAseq.geneNames = Jiwon[keep,1]
temp.idx = which(DATA$GeneID %in% Jiwon.RNAseq.geneNames)
m = data.matrix(DATA[temp.idx,c(6:8)])
svg()
hm = heatmap.2(m,cexCol = 0.8)
dev.off()
file.size('Rplot001')/1000
hc<- hclust(dist(m))
GO.idx = which(cutree(hc, k=2) [hc$order]==2)
GO.genes = DATA[as.array(GO.idx),1]
write.csv('Desktop/Jiwon_Cell_GO1.csv',x = GO.genes)
