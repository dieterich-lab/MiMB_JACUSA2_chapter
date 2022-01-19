library(NMF)
library(pROC)
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
pdf(paste0(args[1],"/Pattern_",k,"_barplot_NMF.pdf"),width=7,height=7)
rep1=matrix(h[k,],ncol=5,byrow=T)[1:3,]
rownames(rep1)=c("BASE","DEL","INS")
colnames(rep1)=c("Pos1","Pos2","Pos3","Pos4","Pos5")
if (k ==1) rep = rep1 else rep = rep1 + rep
barplot(rep1,legend=F, main=paste0("NMF_Pattern_",k),col=cbbPalette,cex.names=2,cex.axis=2)
dev.off()
}
pdf(paste0(args[1],"/Pattern_Sum_barplot_NMF.pdf"),width=7,height=7)

barplot(rep,legend=F, main=paste0("NMF_Sum_Pattern"),col=cbbPalette,cex.names=2,cex.axis=2)
dev.off()


#                                         #Choice of best pattern for subsequent analysis

wExtended=data.frame(w,BigTable[rownames(w),"Motif"])
print("Table of number of best hits for each pattern")
TabPatternInstances=table(apply(wExtended[,1:nrow(h)],1,function(x){which(x==max(x))}))
print(TabPatternInstances)
#chosenPattern=which(TabPatternInstances==max(TabPatternInstances))
#print(paste0("We suggest to use pattern:",chosenPattern))
tt=apply(wExtended[,1:nrow(h)],1,function(x){which(x==max(x))})

# # Compute CountMatrix                                  
# CountMatrix<-function(inp,pseudo=0.01){
# factor3=t(as.matrix((sapply(sapply(inp,strsplit,""),unlist))))
# factor3=apply(factor3,2,table)
# #add function
# factor3Mat=matrix(0,nrow=4,ncol=5)
# rownames(factor3Mat)<-c("A","C","G","T")
# colnames(factor3Mat)<-1:5

# for(k in 1:5)
# {
#     for(l in names(factor3[[k]]))
#     {
#         factor3Mat[l,k]<-as.numeric(factor3[[k]][l])
#     }
# }
# return(factor3Mat)
# }

library(ggplot2)                                        #replace with ggseqlogo
library(ggseqlogo)
# for(k in 1:nrow(h))
#     {
#          ggseqlogo(BigTable[names(which(tt==k)),"Motif"])
# #        logomaker(apply(CountMatrix(BigTable[names(which(tt==k)),"Motif"]),2,function(x){x/sum(x)}), type = "EDLogo", bg=apply(CountMatrix(BigTable[names(which(tt!=k)),"Motif"]),2,function(x){x/sum(x)}))
#         #apply(CountMatrix(BigTable[names(which(tt!=k)),"Motif"]),2,function(x){x/sum(x)})
#         ggsave(paste0(args[1],"/SeqLogo_",k,"_NMF.pdf"),width=7,height=7,device="pdf")
       
#         }
AllSites=as.matrix(BigTable[,16:30])
AllSitesScores=AllSites%*%t(h)
AllSitesScores_sum=AllSites%*% colSums(h)

Noverlap=read.delim(args[2],header=F,as.is=T)
Noverlap.df=data.frame(ID=paste0(Noverlap[,1],":",Noverlap[,2]-2,"_",Noverlap[,3]+2,":",Noverlap[,6]),Type=Noverlap[,4])
rownames(Noverlap.df)=Noverlap.df[,1]

ins=intersect(rownames(BigTable),Noverlap.df$ID); #
ins2= intersect(rownames(BigTable),subset(Noverlap.df$ID,Noverlap.df$Type=="Boulias,Koertel,Koh"))
ins3= intersect(rownames(BigTable),subset(Noverlap.df$ID,Noverlap.df$Type!="Boulias,Koertel,Koh"))

dataGG=data.frame(AllSitesScores[ins,],Noverlap.df[ins,2])
print(dim(dataGG))
colnames(dataGG)=c(paste("NMF",1:nrow(h),sep=""),"CLIP")

#*************************************BEGIN. Modified part *********************************************
TestSitesScores = AllSitesScores[!(row.names(AllSitesScores) %in% ins2),]

dataNG = data.frame(TestSitesScores,rep("NO",dim(TestSitesScores)[1]))
colnames(dataNG)=c(paste("NMF",1:nrow(h),sep=""),"CLIP")
dataNG[ins,"CLIP"]=Noverlap.df[ins,2]
write.csv(dataNG,paste0(args[1],"/wt_ivt.csv"),row.names = FALSE) 
for(dude in colnames(dataNG)[1:nrow(h)])
    {
xdensity <- ggplot(dataNG, aes_string(x=dude, color="CLIP")) + 
stat_ecdf() +
scale_color_manual(values = cbbPalette) + 
theme_bw()

# critical need to check 
ggsave(paste0(args[1],'/',dude,"_ecdf_3.pdf"),device="pdf")
}

# for(dude in colnames(dataNG)[1:nrow(h)])
#     {
# rocobj = roc(dataNG[,"CLIP"], dataNG[,dude], plot= TRUE)
# auc <- round(auc(dataNG[,"CLIP"], dataNG[,dude]),4)
# g = ggroc(rocobj, colour = 'steelblue', size = 2) 
# g+ ggtitle(paste0(dude,'_ROC Curve ', '(AUC = ', auc, ')'))+    theme_minimal()

# # #critical need to check 
# ggsave(paste0(args[1],"/", dude,"_roc.pdf"),device="pdf")
# }
dataNG = data.frame(AllSitesScores_sum,rep("NO",dim(AllSitesScores_sum)[1]))
colnames(dataNG)=c("Sum_NMF","CLIP")
dataNG[ins,"CLIP"]=Noverlap.df[ins,2]
write.csv(dataNG,paste0(args[1],"/wt_ivt_sum.csv"),row.names = FALSE) 

xdensity <- ggplot(dataNG, aes_string(x='Sum_NMF', color="CLIP")) + 
stat_ecdf() +
scale_color_manual(values = cbbPalette) + 
theme_bw()

#critical need to check 
ggsave(paste0(args[1],"/NMF_SUM_ecdf_2.pdf"),device="pdf")
print(dim(TestSitesScores))

# dataGG = data.frame(rowSums(AllSitesScores[ ins3,]),Noverlap.df[ins3,2])
# colnames(dataGG)=c("Sum_NMF","CLIP")
# xdensity <- ggplot(dataGG, aes_string(x='Sum_NMF', color="CLIP")) + 
# stat_ecdf() +
# scale_color_manual(values = cbbPalette) + 
# theme_bw()

# #critical need to check 
# ggsave(paste0(args[1],"/NMF_SUM_ecdf.pdf"),device="pdf")

#*************************************END. Modified part *********************************************

# colnames(AllSitesScores)=c(paste("NMF",1:nrow(h),sep=""))
# AllSitesScores=data.frame(AllSitesScores[rownames(BigTable),],TotalScore=apply(AllSitesScores[rownames(BigTable),],1,sum),BigTable[,c("Motif","DRACH")])


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
# for(dude in colnames(dataNG)[1:nrow(h)])
#     {
# xdensity <- ggplot(dataNG, aes_string(x=dude, color="CLIP")) + 
# stat_ecdf() +
# scale_color_manual(values = cbbPalette) + 
# theme_bw()

# #critical need to check 
# ggsave(paste0(args[1],'/',dude,"_ecdf_2.pdf"),device="pdf")
# }


#*************************************END. Modified part *********************************************

# saveRDS(AllSitesScores,file=paste0(args[1],"/ScoreProfile_NMFall_plusNonCLIP.rds"))

#***************************************ROC curves ***************************************************
# print(colnames(dataNG))

rocobj = roc(dataNG$CLIP, dataNG$Sum_NMF, plot= TRUE)
auc <- round(auc(dataNG$CLIP, dataNG$Sum_NMF),4)
g = ggroc(rocobj, colour = 'steelblue', size = 2) 
g+ ggtitle(paste0('SUM_NMF ROC Curve ', '(AUC = ', auc, ')'))+    theme_minimal()

#critical need to check 
ggsave(paste0(args[1],"/SUM_NMF_roc2.pdf"),device="pdf")
