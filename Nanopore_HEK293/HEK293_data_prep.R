#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
                                        #Read preprocessed JACUSA2 output
print("Rebuild matrix of features ...")
print("Read")
prefix<-""; #/Volumes
Exp1=read.delim(paste0(args[1],"/data_reformat.txt"),as.is=T,header=F)
#
#gzip#
                                        #Formatting
print("Format")
colnames(Exp1)=c("ID","contig","position","call2.score","deletion.score","insertion.score","Base","Anchor","Strand")
if (dir.exists(dirname(args[2]))== FALSE) {
  dir.create(dirname(args[2]))
}
#Split by experiment#
# print(eins$RT)
# print(args[2])
# Exp1=subset(eins,eins$RT==args[2])
# Exp2=subset(eins,eins$RT==args[3])
# Exp3=subset(eins,eins$RT==args[4])

print("Exp1")
print(dim(Exp1))
Call2=tapply(Exp1$call2.score,list(Exp1$ID,Exp1$Anchor),sum)
Call2[is.na(Call2)]<-0
colnames(Call2)<-paste0("Exp1Call2Score_",colnames(Call2))
print('ste')

Deletion=tapply(Exp1$deletion.score,list(Exp1$ID,Exp1$Anchor),sum)
Deletion[is.na(Deletion)]<-0
colnames(Deletion)<-paste0("Exp1DeletionScore_",colnames(Deletion))
print('ste')

Insertion=tapply(Exp1$insertion.score,list(Exp1$ID,Exp1$Anchor),sum)
Insertion[is.na(Insertion)]<-0
colnames(Insertion)<-paste0("Exp1InsertionScore_",colnames(Insertion))
print('ste')
BigTable=merge(data.frame(ID=rownames(Call2),Call2),data.frame(ID=rownames(Deletion),Deletion),by.x=1,by.y=1)
BigTable=merge(BigTable,data.frame(ID=rownames(Insertion),Insertion),by.x=1,by.y=1)
# print("Exp2")#

# Call2=tapply(Exp2$call2.score,list(Exp2$ID,Exp2$Anchor),sum)
# Call2[is.na(Call2)]<-0
# colnames(Call2)<-paste0("Exp2Call2Score_",colnames(Call2))

# Deletion=tapply(Exp2$deletion.score,list(Exp2$ID,Exp2$Anchor),sum)
# Deletion[is.na(Deletion)]<-0
# colnames(Deletion)<-paste0("Exp2DeletionScore_",colnames(Deletion))

# Insertion=tapply(Exp2$insertion.score,list(Exp2$ID,Exp2$Anchor),sum)
# Insertion[is.na(Insertion)]<-0
# colnames(Insertion)<-paste0("Exp2InsertionScore_",colnames(Insertion))

#                                         #That is why, we get different results .. force coverage on both exp.

# BigTable=merge(BigTable,data.frame(ID=rownames(Call2),Call2),by.x=1,by.y=1)
# BigTable=merge(BigTable,data.frame(ID=rownames(Deletion),Deletion),by.x=1,by.y=1)
# BigTable=merge(BigTable,data.frame(ID=rownames(Insertion),Insertion),by.x=1,by.y=1)

# print("Exp3")#

# Call2=tapply(Exp3$call2.score,list(Exp3$ID,Exp3$Anchor),sum)
# Call2[is.na(Call2)]<-0
# colnames(Call2)<-paste0("Exp3Call2Score_",colnames(Call2))

# Deletion=tapply(Exp3$deletion.score,list(Exp3$ID,Exp3$Anchor),sum)
# Deletion[is.na(Deletion)]<-0
# colnames(Deletion)<-paste0("Exp3DeletionScore_",colnames(Deletion))

# Insertion=tapply(Exp3$insertion.score,list(Exp3$ID,Exp3$Anchor),sum)
# Insertion[is.na(Insertion)]<-0
# colnames(Insertion)<-paste0("Exp3InsertionScore_",colnames(Insertion))

                                        #That is why, we get different results .. force coverage on both exp.

# BigTable=merge(BigTable,data.frame(ID=rownames(Call2),Call2),by.x=1,by.y=1)
# BigTable=merge(BigTable,data.frame(ID=rownames(Deletion),Deletion),by.x=1,by.y=1)
# BigTable=merge(BigTable,data.frame(ID=rownames(Insertion),Insertion),by.x=1,by.y=1)


motif=read.table(paste0(args[1],"/checkMotif_reformat.txt"),as.is=T,header=F)
# motif[,1]=gsub("-",":-",motif[,1])
# motif[,1]=gsub("\\+",":\\+",motif[,1])

#*************************************BEGIN. Modified part *********************************************
motif[,1]=gsub("-","_",motif[,1])
motif[,1]=gsub("\\(\\_\\)",":-",motif[,1])
motif[,1]=gsub("\\(\\+\\)",":\\+",motif[,1])
#*************************************END. Modified part *********************************************

BigTable=merge(BigTable,motif,by.x=1,by.y=1);

#BigTable=merge(BigTable,data.frame(ID=names(OverLapWithStudies),miCLIP=OverLapWithStudies),by.x=1,by.y=1,all.x=T)
##
                                        #miCLIP is not positioned on the centre of 5mer / recompute overlap !!
#exit(0);
colnames(BigTable)[ncol(BigTable)]="Motif"
rownames(BigTable)=BigTable[,1]
BigTable=BigTable[,-1]

BigTable$DRACH<-rep(0,nrow(BigTable))
BigTable$DRACH[grep("[AGT][AG]AC[ACT]",BigTable$Motif)]<-1


saveRDS(BigTable,file=args[2])


