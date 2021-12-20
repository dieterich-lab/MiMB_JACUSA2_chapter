library(NMF)
BigTable=readRDS("BigTable.rds")
res=readRDS("NMF.rds");
##
pdf("NMF_6cluster_new.pdf")
layout(cbind(1,2))
# basis components
basismap(res)
# mixture coefficients
aga<-coefmap(res)
dev.off()

w<-basis(res)
h <- coef(res)

# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

                                        #revcovered pattern
for(k in 1:nrow(h))
    {
pdf(paste0("Pattern_",k,"_barplot_NMF.pdf"),width=14,height=7)
#split.screen(c(1, 2))
#screen(1)
rep1=matrix(h[k,],ncol=5,byrow=T)[1:3,]
rownames(rep1)=c("BASE","DEL","INS")
colnames(rep1)=c("Pos1","Pos2","Pos3","Pos4","Pos5")
barplot(rep1,legend=F, main="Experiment 1",col=cbbPalette,,cex.names=2,cex.axis=2)
#screen(2)
#rep2=matrix(h[k,],ncol=5,byrow=T)[4:6,]
#rownames(rep2)=c("BASE","DEL","INS")
#colnames(rep2)=c("Pos1","Pos2","Pos3","Pos4","Pos5")
#barplot(rep2,legend=T, main="Experiment 2",col=cbbPalette,,cex.names=2,cex.axis=2)
dev.off()
}
#barplot(hplot[1,sort(colnames(hplot))],beside=T,las=2,col=cbbPalette[1],legend=c("Pattern1","Pattern2","Pattern3","Pattern4","Pattern5")[1])

#order by position and plot each profile separately #colors !

#ids=names(which(wScale[,3]>0.8))
# ok, we have 5 factors by now... we do not see a strong interaction across positions.
#bla, bla-.
#Story line 5mer->NMF->Factor2->ROC->MostDiscrFeature for both data sets.


wExtended=data.frame(w,BigTable[rownames(w),"Motif"])
table(apply(wExtended[,1:nrow(h)],1,function(x){which(x==max(x))}))
tt=apply(wExtended[,1:nrow(h)],1,function(x){which(x==max(x))})

                                        #BigTable[names(which(tt==k)),"Motif"]
                                        #BigTable[names(which(tt!=k)),"Motif"]

CountMatrix<-function(inp,pseudo=0.01){
factor3=t(as.matrix((sapply(sapply(inp,strsplit,""),unlist))))
factor3=apply(factor3,2,table)
#add function
factor3Mat=matrix(0,nrow=4,ncol=5)
rownames(factor3Mat)<-c("A","C","G","T")
colnames(factor3Mat)<-1:5

for(k in 1:5)
{
    for(l in names(factor3[[k]]))
    {
        factor3Mat[l,k]<-as.numeric(factor3[[k]][l])
    }
}
return(factor3Mat)
}

                                        #motif search TODO
library(Logolas)
for(k in 1:nrow(h))
    {
        pdf(paste0("Logo_",k,"_barplot_NMF.pdf"),width=7,height=7)
        logomaker(BigTable[names(which(tt==k)),"Motif"], type = "Logo")
        logomaker(apply(CountMatrix(BigTable[names(which(tt==k)),"Motif"]),2,function(x){x/sum(x)}), type = "EDLogo", bg=apply(CountMatrix(BigTable[names(which(tt!=k)),"Motif"]),2,function(x){x/sum(x)}))
        #apply(CountMatrix(BigTable[names(which(tt!=k)),"Motif"]),2,function(x){x/sum(x)})
        dev.off()
        }

exit(0);
#überprüfen von links und rechts
ScoresExp1=as.matrix(BigTable[,1:15])
NMFScoreExp1=ScoresExp1%*%t(h[,1:15])

#überprüfen von links und rechts IVT => Simulation#
ScoresIVT=as.matrix(BigTable[,16:30])
NMFScoreIVT=ScoresIVT%*%t(h[,1:15])


Noverlap=read.delim("/Volumes//prj/JACUSA2_TestField/Nanopore_HEK293/miCLIP2/miCLIP_union_flat_exclude_Y_chromosome.bed",header=F,as.is=T)
Noverlap.df=data.frame(ID=paste0(Noverlap[,1],":",Noverlap[,2]-2,"_",Noverlap[,3]+2,":",Noverlap[,6]),Type=Noverlap[,4])
rownames(Noverlap.df)=Noverlap.df[,1]

ins=intersect(rownames(BigTable),Noverlap.df$ID);

dataGG=data.frame(NMFScoreIVT[ins,],Noverlap.df[ins,2])
colnames(dataGG)=c(paste("NMF",1:nrow(h),sep=""),"CLIP")
library(ggplot2)
xdensity <- ggplot(dataGG, aes(NMF3, color=CLIP)) + 
stat_ecdf() +
scale_color_manual(values = cbbPalette) + 
theme_bw()

#critical need to check 
ggsave("NMF4_ecdf.pdf",device="pdf")
saveRDS(dataGG,file="ScoreProfile_NMFall.rds")
saveRDS(bla3,file="ScoreProfile_NMFall_plusNonCLIP.rds")

#19:3976465-3976470:- ##
