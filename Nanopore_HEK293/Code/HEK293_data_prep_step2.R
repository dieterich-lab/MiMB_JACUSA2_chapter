#!/usr/bin/env Rscript
#read in base call error profiles
args = commandArgs(trailingOnly=TRUE)

#BigTable=readRDS(paste0(args[1],"/BigTable.rds"))
#read in miCLIP data
#Noverlap=read.delim(args[2],header=F,as.is=T)

BigTable=readRDS("BigTable.rds")

#read in miCLIP data
Noverlap=read.delim("/prj/JACUSA2_TestField/Nanopore_HEK293/miCLIP2/miCLIP_union_flat_exclude_Y_chromosome.bed",header=F,as.is=T)

#TODO: setup parameters: number of cores

Noverlap.df=data.frame(ID=paste0(Noverlap[,1],":",Noverlap[,2]-2,"_",Noverlap[,3]+2,":",Noverlap[,6]),Type=Noverlap[,4])
rownames(Noverlap.df)=Noverlap.df[,1]

ins=intersect(rownames(BigTable),subset(Noverlap.df$ID,Noverlap.df$Type=="Boulias,Koertel,Koh")); #
#ins=intersect(rownames(BigTable),Noverlap.df$ID);

NMFtab=BigTable[ins,]
                                        #train NMF based on core miCLIP data WT_vs_KO shared across all 3 experiment

library(NMF)
nmfSeed('nndsvd')
meth <- nmfAlgorithm(version='R')
meth <- c(names(meth), meth)
NMFtabSlim=NMFtab[,args[1]]

estim.r <- nmf(NMFtabSlim, 2:10, nrun=10, seed=123456, .opt='vp3')

V.random <- randomize(NMFtabSlim)
# estimate quality measures from the shuffled data (use default NMF algorithm)
estim.r.random <- nmf(V.random, 2:10, nrun=10, seed=123456, .opt='vp3')

DeltaSil=estim.r$measures$silhouette.consensus-estim.r.random$measures$silhouette.consensus
DeltaCoph=estim.r$measures$cophenetic-estim.r.random$measures$cophenetic

ChoseRank=min(which(DeltaSil==max(DeltaSil))+1,which(DeltaCoph==max(DeltaCoph))+1)
pdf("NMF_assess.pdf")
plot(estim.r,estim.r.random)
dev.off()
#exit(0);
#how can we select factorization rank ?
res <- nmf(NMFtabSlim, ChoseRank, nrun=10, seed=123456, .opt='vp3')

saveRDS(res,file="NMF.rds")

