library(NMF)
args = commandArgs(trailingOnly=TRUE)

if (dir.exists(dirname(args[2]))== FALSE) {
  dir.create(dirname(args[2]))
}

#load training data
NMFtab = readRDS(args[1])
nmfSeed('nndsvd')
meth <- nmfAlgorithm(version='R')
meth <- c(names(meth), meth)
NMFtabSlim=NMFtab[,1:15]
estim.r <- nmf(NMFtabSlim, 2:10, nrun=10, seed=123456, .opt='vp3')
V.random <- randomize(NMFtabSlim)

# estimate quality measures from the shuffled data (use default NMF algorithm)
estim.r.random <- nmf(V.random, 2:10, nrun=10, seed=123456, .opt='vp3')
DeltaSil=estim.r$measures$silhouette.consensus-estim.r.random$measures$silhouette.consensus
DeltaCoph=estim.r$measures$cophenetic-estim.r.random$measures$cophenetic
ChoseRank=min(which(DeltaSil==max(DeltaSil))+1,which(DeltaCoph==max(DeltaCoph))+1)
val_matrix = matrix(, nrow = length(DeltaCoph), ncol = 2)
val_matrix[,1] = DeltaSil 
val_matrix[,2] = DeltaCoph
pdf(args[4])
barplot(t(val_matrix), names.arg = 2:(length(DeltaCoph)+1), beside = TRUE, col = c("#009E73", "#D55E00"), legend.text = c("Silhouette", "Cophenetic"),ylim=range(pretty(c(0, DeltaSil))), xlab="Rank", 
        ylab="Delta", main = "Difference between original and randomized data",cex.names=2,cex.axis=2)
abline(h=c(max(DeltaSil) , max(DeltaCoph)) , col=c("#009E73", "#D55E00"))
dev.off()
pdf(args[3])
plot(estim.r,estim.r.random)
dev.off()
#get patterns
res <- nmf(NMFtabSlim, ChoseRank, nrun=10, seed=123456, .opt='vp3')
saveRDS(res,file=args[2])

