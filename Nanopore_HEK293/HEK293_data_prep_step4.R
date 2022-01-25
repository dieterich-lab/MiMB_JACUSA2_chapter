library(NMF)
library(ggplot2)   

print("Predict modification ...")
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

BigTable=readRDS(snakemake@input[[1]])
res=readRDS(snakemake@input[[2]]);
AllSites=as.matrix(BigTable[,1:15])
h = coef(res)
h = rbind(h, t(colSums(h)))
h = rbind(h, t(colSums(h[c(2,3,6,7),])))
h = rbind(h, t(colSums(h[c(2,3,6,7,4),])))

AllSitesScores=AllSites%*%t(h)
print(dim(AllSitesScores))
colnames(AllSitesScores)=c(paste("Pattern",1:nrow(h),sep=""))
AllSitesScores=data.frame(AllSitesScores[rownames(BigTable),],TotalScore=apply(AllSitesScores[rownames(BigTable),],1,sum),BigTable[,c("Motif","DRACH")])

if (dir.exists(snakemake@output[[1]])== FALSE) {
  dir.create(snakemake@output[[1]])
}
saveRDS(AllSitesScores,file=paste0(snakemake@output[[1]],"/ScoreProfile_NMFall_plusNonCLIP.rds"))


Noverlap=read.delim(snakemake@params[[1]],header=F,as.is=T)
Noverlap.df=data.frame(ID=paste0(Noverlap[,1],":",Noverlap[,2]-2,"_",Noverlap[,3]+2,":",Noverlap[,6]),Type=Noverlap[,4])
rownames(Noverlap.df)=Noverlap.df[,1]

ins2= intersect(rownames(AllSitesScores),subset(Noverlap.df$ID,Noverlap.df$Type!="Boulias,Koertel,Koh"))
ins=intersect(rownames(AllSitesScores),Noverlap.df$ID);
NoCLIPsites=setdiff(rownames(AllSitesScores),ins)
print(length(ins2))
print(length(ins))
print(length(NoCLIPsites))

dataGG=data.frame(AllSitesScores[c(ins,NoCLIPsites),-c(ncol(AllSitesScores)-1,ncol(AllSitesScores))],CLIP=c(as.character(Noverlap.df[ins,2]),rep("NoCLIP",length(NoCLIPsites))))
NumberOfFactors=ncol(AllSitesScores)-2
colnames(dataGG)=c(paste("Pattern",1:NumberOfFactors,sep=""),"CLIP")


for(dude in colnames(dataGG)[9:NumberOfFactors])
    {
xdensity <- ggplot(dataGG, aes_string(dude, color="CLIP")) + 
stat_ecdf() +
scale_color_manual(values = cbbPalette) + theme_bw()+  scale_x_continuous(name = "score")+theme(
    axis.title.x = element_text(size=rel(2)),      
    axis.title.y = element_text(size=rel(2)),
    axis.text = element_text(size = rel(1.5)))

#critical need to check 
ggsave(paste0(snakemake@output[[1]],"/",dude,"_ecdf.pdf"),device="pdf")
}
######################################### eCDF with two classws###################################
dataGG_<-data.frame(dataGG[,1:NumberOfFactors],CLIP=c(rep("CLIP",length(ins)),rep("noCLIP",length(NoCLIPsites))))

for(dude in colnames(dataGG_)[8:NumberOfFactors])
    {
xdensity <- ggplot(dataGG_, aes_string(dude, color="CLIP")) + 
stat_ecdf() +
scale_color_manual(values = cbbPalette) + 
theme_bw()

#critical need to check 
ggsave(paste0(snakemake@output[[1]],"/",dude,"_ecdf_2.pdf"),device="pdf")
}


##############################################PPV PLot ############################################
library(pROC)

# A few constants
temperatureColor <- "#009E73"
priceColor <- "#D55E00"
dataGG=data.frame(AllSitesScores[c(ins2,NoCLIPsites),-c(ncol(AllSitesScores)-1,ncol(AllSitesScores))],CLIP=c(as.character(Noverlap.df[ins2,2]),rep("NoCLIP",length(NoCLIPsites))))

x<-data.frame(dataGG[,-ncol(dataGG)],CLIP=c(rep(1,length(ins2)),rep(0,length(NoCLIPsites))))
for(dude in colnames(x)[9:NumberOfFactors])
    {
print(dude)
# outliers <- boxplot(dataGG[dude], plot=FALSE)$out
# x<- x[-which(x[,dude] %in% outliers),]
r <- roc(x[,"CLIP"], x[,dude])
coordinates <- coords(r, x = "all", input = "threshold", ret = c("threshold","ppv","tpr","tp", "tn","fp","fn"))
coordinates$sum <- log2(coordinates$tp + coordinates$fp)
coordinates$ppv <- coordinates$ppv*20 
# coordinates$thre <- coordinates$threshold
ggplot(coordinates, aes(x=threshold)) +  
  geom_line( aes(y=sum), size=2, color="#009E73") + 
  geom_line( aes(y=ppv), size =2,color="#D55E00") +
  scale_y_continuous(    
    # Features of the first axis
    name = "Log2(#predictions)",
    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*5,name="Positive Predictive Value")
  ) +  scale_x_continuous(name = "Cutoff")+  theme_bw() + theme(
    axis.title.x = element_text(size=rel(2)),      
    axis.title.y = element_text(color = "#009E73", size=rel(2)),
    axis.title.y.right = element_text(color = "#D55E00", size=rel(0.9)), 
    axis.text = element_text(size = rel(1.5))

  ) 
ggsave(paste0(snakemake@output[[1]],"/",dude,"_ppv.pdf"),device="pdf")

    
}
