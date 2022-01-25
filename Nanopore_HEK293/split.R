BigTable=readRDS(snakemake@input[[1]])
Noverlap=read.delim(snakemake@params[[1]],header=F,as.is=T)
Noverlap.df=data.frame(ID=paste0(Noverlap[,1],":",Noverlap[,2]-2,"_",Noverlap[,3]+2,":",Noverlap[,6]),Type=Noverlap[,4])
rownames(Noverlap.df)=Noverlap.df[,1]
ins=intersect(rownames(BigTable),subset(Noverlap.df$ID,Noverlap.df$Type=="Boulias,Koertel,Koh")); #
saveRDS(BigTable[ins,],file=snakemake@output[[1]])
saveRDS(BigTable[!(rownames(BigTable)%in%ins),],file=snakemake@output[[2]])

