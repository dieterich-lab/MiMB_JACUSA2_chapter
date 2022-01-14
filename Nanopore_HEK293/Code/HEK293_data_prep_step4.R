#library(NMF)
#BigTable=readRDS("BigTable.rds")
#res=readRDS("NMF.rds");

AllSitesScores=readRDS("ScoreProfile_NMFall_plusNonCLIP.rds")

#Just focus on DRACH site
AllSitesScores=subset(AllSitesScores,AllSitesScores$DRACH==1)

# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

Noverlap=read.delim("/Volumes/prj/JACUSA2_TestField/Nanopore_HEK293/miCLIP2/miCLIP_union_flat_exclude_Y_chromosome.bed",header=F,as.is=T)
Noverlap.df=data.frame(ID=paste0(Noverlap[,1],":",Noverlap[,2]-2,"_",Noverlap[,3]+2,":",Noverlap[,6]),Type=Noverlap[,4])
rownames(Noverlap.df)=Noverlap.df[,1]

ins=intersect(rownames(AllSitesScores),Noverlap.df$ID);
NoCLIPsites=setdiff(rownames(AllSitesScores),ins)

print(length(ins))
print(length(NoCLIPsites))

dataGG=data.frame(AllSitesScores[c(ins,NoCLIPsites),-c(ncol(AllSitesScores)-1,ncol(AllSitesScores))],CLIP=c(as.character(Noverlap.df[ins,2]),rep("NoCLIP",length(NoCLIPsites))))
NumberOfFactors=ncol(AllSitesScores)-2
colnames(dataGG)=c(paste("NMF",1:NumberOfFactors,sep=""),"CLIP")

library(ggplot2)

                                        #All experimental below

#require(plotROC)
#ROC plots do not make sense.
#for(dude in colnames(dataGG)[1:NumberOfFactors])
#    {
#        rocca<-data.frame(NMF2=dataGG[,dude],CLIP=c(rep(1,length(ins)),rep(0,length(NoCLIPsites))))
#        p<-ggplot(rocca, aes(m = NMF2, d = CLIP)) + geom_roc(labels=F)
#        p<-p+geom_abline(slope=1,intercept=0)
#        ggsave(paste0(dude,"_roc.pdf"),device="pdf")
#        calc_auc(p)
#    }


p <- ggplot(data.frame(PerfbgSub[which(PerfbgSub[,2]>10),]), aes(x = Score))
p <- p + ylab("% predicted targets")

p <- p + ggtitle("AATF (55315, 55317) vs eCLIP et al. (GFPcdAdar subtracted)")
#p <- p + ggtitle("AATF (52653, 52655) vs eCLIP et al. (GFPcdAdar subtracted)")

  p <- p + geom_line(aes(y = PercentOverlap, colour = "% overlap with eCLIP"))
  
  # adding the relative humidity data, transformed to match roughly the range of the temperature
  p <- p + geom_line(aes(y = PercentOriginal, colour = "% predicted targets"))

  # now adding the secondary axis, following the example in the help file ?scale_y_continuous
  # and, very important, reverting the above transformation
  p <- p + scale_y_continuous(sec.axis = sec_axis(~.*1, name = "% overlap with eCLIP"))
  
  # modifying colours and theme options
  p <- p + scale_colour_manual(values = c("blue", "red"))

for(dude in colnames(dataGG)[1:NumberOfFactors])
    {
        xdensity <- ggplot(dataGG, aes_string(dude, color="CLIP")) + 
            stat_ecdf() +
            scale_color_manual(values = cbbPalette) + 
            theme_bw() + xlim(0, 10)
#critical need to check 
        ggsave(paste0(dude,"_ecdf.pdf"),device="pdf")
    }


#19:3976465-3976470:- ##
