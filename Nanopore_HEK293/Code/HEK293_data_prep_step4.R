#library(NMF)
#BigTable=readRDS("BigTable.rds")
#res=readRDS("NMF.rds");
args = commandArgs(trailingOnly=TRUE)

AllSitesScores=readRDS(paste0(args[1],"/ScoreProfile_NMFall_plusNonCLIP.rds"))

# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

Noverlap=read.delim(args[2],header=F,as.is=T)
Noverlap.df=data.frame(ID=paste0(Noverlap[,1],":",Noverlap[,2]-2,"_",Noverlap[,3]+2,":",Noverlap[,6]),Type=Noverlap[,4])
rownames(Noverlap.df)=Noverlap.df[,1]

ins=intersect(rownames(AllSitesScores),Noverlap.df$ID);
NoCLIPsites=setdiff(rownames(AllSitesScores),ins)

print(length(ins))
print(length(NoCLIPsites))

dataGG=data.frame(AllSitesScores[c(ins,NoCLIPsites),-c(ncol(AllSitesScores)-1,ncol(AllSitesScores))],CLIP=c(as.character(Noverlap.df[ins,2]),rep("NoCLIP",length(NoCLIPsites))))
NumberOfFactors=ncol(AllSitesScores)-2
colnames(dataGG)=c(paste("Pattern",1:NumberOfFactors,sep=""),"CLIP")

rocca<-data.frame(NMF2=dataGG[,"NMF2"],CLIP=c(rep(1,length(ins)),rep(0,length(NoCLIPsites))))
library(ggplot2)
# require(plotROC)
# p<-ggplot(rocca, aes(m = NMF2, d = CLIP)) + geom_roc(labels=F)
# #color = Base_1
# p<-p+geom_abline(slope=1,intercept=0)
# p
# ggsave(paste0(args[1],"roc.pdf"))
# #ok, kinda works
# calc_auc(p)

# for(dude in colnames(dataGG)[1:NumberOfFactors])
#     {
# xdensity <- ggplot(dataGG, aes_string(dude, color="CLIP")) + 
# stat_ecdf() +
# scale_color_manual(values = cbbPalette) + 
# theme_bw()

# #critical need to check 
# ggsave(paste0(args[1],"/",dude,"_ecdf.pdf"),device="pdf")
# }
########################################## eCDF with two classws###################################
# dataGG_<-data.frame(dataGG[,1:NumberOfFactors],CLIP=c(rep("CLIP",length(ins)),rep("noCLIP",length(NoCLIPsites))))

# for(dude in colnames(dataGG_)[1:NumberOfFactors])
#     {
# xdensity <- ggplot(dataGG_, aes_string(dude, color="CLIP")) + 
# stat_ecdf() +
# scale_color_manual(values = cbbPalette) + 
# theme_bw()

# #critical need to check 
# ggsave(paste0(args[1],"/",dude,"_ecdf_2.pdf"),device="pdf")
# }


##############################################PPV PLot ############################################
library(pROC)

# A few constants
temperatureColor <- "#009E73"
priceColor <- "#D55E00"
x<-data.frame(dataGG[,-ncol(dataGG)],CLIP=c(rep(1,length(ins)),rep(0,length(NoCLIPsites))))
for(dude in colnames(x)[1:NumberOfFactors])
    {
print(dude)
# outliers <- boxplot(dataGG[dude], plot=FALSE)$out
# x<- x[-which(x[,dude] %in% outliers),]
r <- roc(x[,"CLIP"], x[,dude])
coordinates <- coords(r, x = "all", input = "threshold", ret = c("threshold","ppv","tp", "tn","fp","fn"))
coordinates$sum <- (coordinates$tp + coordinates$fp)/(coordinates$tp + coordinates$fp+ coordinates$tn+ coordinates$fn)*100
coordinates$ppv <- coordinates$ppv*100 
# coordinates$thre <- coordinates$threshold
ggplot(coordinates, aes(x=threshold)) +  
  geom_line( aes(y=sum), size=2, color="#009E73") + 
  geom_line( aes(y=ppv), size =2,color="#D55E00") +
  scale_y_continuous(    
    # Features of the first axis
    name = "Predicted Targets Rate",
    # Add a second axis and specify its features
    sec.axis = sec_axis(~.,name="Positive Predictive Value")
  ) +  scale_x_continuous(name = paste0('Score (',dude,')') )+  theme_bw() + theme(
    axis.title.x = element_text(size=rel(1.5)),      
    axis.title.y = element_text(color = "#009E73", size=rel(2)),
    axis.title.y.right = element_text(color = "#D55E00", size=rel(1)), 
    axis.text = element_text(size = rel(1.5))

  ) 
  ggtitle(dude)
ggsave(paste0(args[1],"/",dude,"_ppv.pdf"),device="pdf")}
#19:3976465-3976470:- ##