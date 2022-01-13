#library(NMF)
#BigTable=readRDS("BigTable.rds")
#res=readRDS("NMF.rds");

AllSitesScores=readRDS("ScoreProfile_NMFall_plusNonCLIP.rds")

# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

Noverlap=read.delim("/prj/JACUSA2_TestField/Nanopore_HEK293/miCLIP2/miCLIP_union_flat_exclude_Y_chromosome.bed",header=F,as.is=T)
Noverlap.df=data.frame(ID=paste0(Noverlap[,1],":",Noverlap[,2]-2,"_",Noverlap[,3]+2,":",Noverlap[,6]),Type=Noverlap[,4])
rownames(Noverlap.df)=Noverlap.df[,1]

ins=intersect(rownames(AllSitesScores),Noverlap.df$ID);
NoCLIPsites=setdiff(rownames(AllSitesScores),ins)

print(length(ins))
print(length(NoCLIPsites))

dataGG=data.frame(AllSitesScores[c(ins,NoCLIPsites),-c(ncol(AllSitesScores)-1,ncol(AllSitesScores))],CLIP=c(Noverlap.df[ins,2],rep("NoCLIP",length(NoCLIPsites))))
NumberOfFactors=ncol(AllSitesScores)-2
colnames(dataGG)=c(paste("NMF",1:NumberOfFactors,sep=""),"CLIP")

library(ggplot2)
require(plotROC)
p<-ggplot(dataGG, aes(m = NMF2, d = CLIP)) + geom_roc(labels=F)
#color = Base_1
p<-p+geom_abline(slope=1,intercept=0)
p
ggsave("roc.pdf")

calc_auc(p)

for(dude in colnames(dataGG)[1:NumberOfFactors])
    {
xdensity <- ggplot(dataGG, aes_string(dude, color="CLIP")) + 
stat_ecdf() +
scale_color_manual(values = cbbPalette) + 
theme_bw()

#critical need to check 
ggsave(paste0(dude,"_ecdf.pdf"),device="pdf")
}


#19:3976465-3976470:- ##
