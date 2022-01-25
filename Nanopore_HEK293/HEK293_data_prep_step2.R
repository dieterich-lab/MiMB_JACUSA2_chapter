library(NMF)

NMFtab = readRDS(snakemake@input[[1]])

if (dir.exists(dirname(snakemake@output[[1]]))== FALSE) {
  dir.create(dirname(snakemake@output[[1]]))
}

nmfSeed('nndsvd')
meth <- nmfAlgorithm(version='R')
meth <- c(names(meth), meth)
NMFtabSlim=NMFtab[,1:15]
print(dim(NMFtabSlim))
estim.r <- nmf(NMFtabSlim, 2:10, nrun=10, seed=123456, .opt='vp3')

V.random <- randomize(NMFtabSlim)
# estimate quality measures from the shuffled data (use default NMF algorithm)
estim.r.random <- nmf(V.random, 2:10, nrun=10, seed=123456, .opt='vp3')

DeltaSil=estim.r$measures$silhouette.consensus-estim.r.random$measures$silhouette.consensus
DeltaCoph=estim.r$measures$cophenetic-estim.r.random$measures$cophenetic

ChoseRank=min(which(DeltaSil==max(DeltaSil))+1,which(DeltaCoph==max(DeltaCoph))+1)

pdf(snakemake@output[[2]])
plot(estim.r,estim.r.random)
dev.off()
#exit(0);
#how can we select factorization rank ?
print(ChoseRank)
res <- nmf(NMFtabSlim, ChoseRank, nrun=10, seed=123456, .opt='vp3')
saveRDS(res,file=snakemake@output[[1]])

