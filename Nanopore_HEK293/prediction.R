library(NMF)
library(ggplot2)   
library(stringr)
library(pROC)

print("Predict modification ...")

if (dir.exists(snakemake@output[[1]])== FALSE) {
      dir.create(snakemake@output[[1]])
    }

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
########################################### get scores #############################################
#load data and patterns
BigTable=readRDS(snakemake@input[[1]])
res=readRDS(snakemake@input[[2]]);
AllSites=as.matrix(BigTable[,1:15])
h = coef(res)
h = rbind(h, t(colSums(h)))
#add patterns combination
for (pt in snakemake@params[[2]]){
    h = rbind(h, t(colSums(h[pt,])))
    }

#generate scores for all sites
AllSitesScores=AllSites%*%t(h)
colnames(AllSitesScores)=c(paste("Pattern",1:nrow(h),sep=""))
AllSitesScores=data.frame(AllSitesScores[rownames(BigTable),],BigTable[,c("Motif","DRACH")])


# saveRDS(AllSitesScores,file=paste0(snakemake@output[[1]],"/ScoreProfile_NMFall_plusNonCLIP.rds"))
#get modified sites
Noverlap=read.delim(snakemake@params[[1]],header=F,as.is=T)
Noverlap.df=data.frame(ID=paste0(Noverlap[,1],":",Noverlap[,2]-2,"_",Noverlap[,3]+2,":",Noverlap[,6]),Type=Noverlap[,4])
rownames(Noverlap.df)=Noverlap.df[,1]

ins2= intersect(rownames(AllSitesScores),subset(Noverlap.df$ID,Noverlap.df$Type!=snakemake@params[[3]]))
ins=intersect(rownames(AllSitesScores),Noverlap.df$ID);
NoCLIPsites=setdiff(rownames(AllSitesScores),ins)

#get scores for the test set
TestSitesScores = AllSitesScores[c(ins2,NoCLIPsites),]
TestSitesScores=data.frame(str_split_fixed(rownames(TestSitesScores), "[:_]", 4),TestSitesScores)
colnames(TestSitesScores)[1:4]=c('chr','start','end','strand')

#add CLIP label to the scores
dataGG=data.frame(AllSitesScores[c(ins,NoCLIPsites),-c(ncol(AllSitesScores)-1,ncol(AllSitesScores))],CLIP=c(as.character(Noverlap.df[ins,2]),rep("NoCLIP",length(NoCLIPsites))))
NumberOfFactors=ncol(AllSitesScores)-2
# colnames(dataGG)=c(paste("Pattern",1:NumberOfFactors,sep=""),"CLIP")

######################################### plot eCDF ###################################

for(dude in colnames(dataGG)[8:NumberOfFactors])
    {
    xdensity <- ggplot(dataGG, aes_string(dude, color="CLIP")) + 
    stat_ecdf() +
    scale_color_manual(values = cbbPalette) + theme_bw()+  scale_x_continuous(name = "score")+theme(
    axis.title.x = element_text(size=rel(2)),      
    axis.title.y = element_text(size=rel(2)),
    axis.text = element_text(size = rel(1.5)))
    ggsave(paste0(snakemake@output[[1]],"/",dude,"_ecdf.pdf"),device="pdf")
    }
# ######################################### plot eCDF with two classes###################################
# dataGG_<-data.frame(dataGG[,1:NumberOfFactors],CLIP=c(rep("CLIP",length(ins)),rep("noCLIP",length(NoCLIPsites))))
# print(colnames(dataGG_))

# for(dude in colnames(dataGG_)[9:NumberOfFactors])
#     {
#     xdensity <- ggplot(dataGG_, aes_string(dude, color="CLIP")) + 
#     stat_ecdf() +
#     scale_color_manual(values = cbbPalette) + 
#     theme_bw()
#     ggsave(paste0(snakemake@output[[1]],"/",dude,"_ecdf_2.pdf"),device="pdf")
#     }


# ##############################################plot PPV ############################################
x=data.frame(TestSitesScores[,-c(1,2,3,4,ncol(TestSitesScores)-1,ncol(TestSitesScores))],CLIP=c(rep(1,length(ins2)),rep(0,length(NoCLIPsites))))
# x<-data.frame(dataGG[,-ncol(dataGG)],CLIP=c(rep(1,length(ins2)),rep(0,length(NoCLIPsites))))
for(dude in colnames(x)[8:NumberOfFactors])
    {
    print(dude)
    r <- roc(x[,"CLIP"], x[,dude])
    coordinates <- coords(r, x = "all", input = "threshold", ret = c("threshold","ppv","tpr","tp", "tn","fp","fn"))
    coordinates$sum <- log2(coordinates$tp + coordinates$fp)
    coordinates$ppv <- coordinates$ppv*20 
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
    scores= TestSitesScores[, c('chr','start','end','Motif',dude,'strand')]
    colnames(scores)[4:5] = c('name','score')
    write.table(scores,paste0(snakemake@output[[1]],"/",dude,"_prediction.bed"), row.names= FALSE)
    }
