args = commandArgs(trailingOnly=TRUE)
BigTable=readRDS(args[1])
Noverlap=read.delim(args[2],header=F,as.is=T)
Noverlap.df=data.frame(ID=paste0(Noverlap[,1],":",Noverlap[,2]-2,"_",Noverlap[,3]+2,":",Noverlap[,6]),Type=Noverlap[,4])
rownames(Noverlap.df)=Noverlap.df[,1]
ins=intersect(rownames(BigTable),subset(Noverlap.df$ID,Noverlap.df$Type==args[5])); #
saveRDS(BigTable[ins,],file=args[3])
saveRDS(BigTable[!(rownames(BigTable)%in%ins),],file=args[4])

