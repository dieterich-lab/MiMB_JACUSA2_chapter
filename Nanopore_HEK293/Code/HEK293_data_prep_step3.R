library(NMF)
args = commandArgs(trailingOnly=TRUE)

BigTable=readRDS(paste0(args[1],"/BigTable.rds"))
res=readRDS(paste0(args[1],"/NMF.rds"));
##
pdf(paste0(args[1],"/NMF_cluster_new.pdf"))
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
pdf(paste0(args[1],"/Pattern_",k,"_barplot_NMF.pdf"),width=14,height=7)
rep1=matrix(h[k,],ncol=5,byrow=T)[1:3,]
rownames(rep1)=c("BASE","DEL","INS")
colnames(rep1)=c("Pos1","Pos2","Pos3","Pos4","Pos5")
barplot(rep1,legend=F, main=paste0("NMF_Pattern_",k),col=cbbPalette,cex.names=2,cex.axis=2)
dev.off()
}

                                        #Choice of best pattern for subsequent analysis

wExtended=data.frame(w,BigTable[rownames(w),"Motif"])
print("Table of number of best hits for each pattern")
TabPatternInstances=table(apply(wExtended[,1:nrow(h)],1,function(x){which(x==max(x))}))
print(TabPatternInstances)
#chosenPattern=which(TabPatternInstances==max(TabPatternInstances))
#print(paste0("We suggest to use pattern:",chosenPattern))
tt=apply(wExtended[,1:nrow(h)],1,function(x){which(x==max(x))})

#Compute CountMatrix                                  
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

library(ggplot2)                                        #replace with ggseqlogo
library(ggseqlogo)
for(k in 1:nrow(h))
    {
         ggseqlogo(BigTable[names(which(tt==k)),"Motif"])
#        logomaker(apply(CountMatrix(BigTable[names(which(tt==k)),"Motif"]),2,function(x){x/sum(x)}), type = "EDLogo", bg=apply(CountMatrix(BigTable[names(which(tt!=k)),"Motif"]),2,function(x){x/sum(x)}))
        #apply(CountMatrix(BigTable[names(which(tt!=k)),"Motif"]),2,function(x){x/sum(x)})
        ggsave(paste0(args[1],"/SeqLogo_",k,"_NMF.pdf"),width=7,height=7,device="pdf")
       
        }
AllSites=as.matrix(BigTable[,1:15])
AllSitesScores=AllSites%*%t(h)

Noverlap=read.delim(args[2],header=F,as.is=T)
Noverlap.df=data.frame(ID=paste0(Noverlap[,1],":",Noverlap[,2]-2,"_",Noverlap[,3]+2,":",Noverlap[,6]),Type=Noverlap[,4])
rownames(Noverlap.df)=Noverlap.df[,1]

ins=intersect(rownames(BigTable),Noverlap.df$ID);

dataGG=data.frame(AllSitesScores[ins,],Noverlap.df[ins,2])
colnames(dataGG)=c(paste("NMF",1:nrow(h),sep=""),"CLIP")

#*************************************BEGIN. Modified part *********************************************
dataNG = data.frame(AllSitesScores,rep('NO.CLIP',dim(AllSitesScores)[1]))
colnames(dataNG)=c(paste("NMF",1:nrow(h),sep=""),"CLIP")
dataNG[ins,"CLIP"]='CLIP'
#*************************************END. Modified part *********************************************

colnames(AllSitesScores)=c(paste("NMF",1:nrow(h),sep=""))
AllSitesScores=data.frame(AllSitesScores[rownames(BigTable),],TotalScore=apply(AllSitesScores[rownames(BigTable),],1,sum),BigTable[,c("Motif","DRACH")])

for(dude in colnames(dataGG)[1:nrow(h)])
    {
xdensity <- ggplot(dataGG, aes_string(x=dude, color="CLIP")) + 
stat_ecdf() +
scale_color_manual(values = cbbPalette) + 
theme_bw()

#critical need to check 
ggsave(paste0(args[1],'/',dude,"_ecdf.pdf"),device="pdf")
}

#*************************************BEGIN. Modified part *********************************************
for(dude in colnames(dataNG)[1:nrow(h)])
    {
xdensity <- ggplot(dataNG, aes_string(x=dude, color="CLIP")) + 
stat_ecdf() +
scale_color_manual(values = cbbPalette) + 
theme_bw()

#critical need to check 
ggsave(paste0(args[1],'/',dude,"_ecdf_2.pdf"),device="pdf")
}
#*************************************END. Modified part *********************************************

saveRDS(AllSitesScores,file=paste0(args[1],"/ScoreProfile_NMFall_plusNonCLIP.rds"))

