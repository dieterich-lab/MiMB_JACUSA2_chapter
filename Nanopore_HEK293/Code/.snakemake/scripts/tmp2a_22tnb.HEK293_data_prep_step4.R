
######## snakemake preamble start (automatically inserted, do not edit) ########
library(methods)
Snakemake <- setClass(
    "Snakemake",
    slots = c(
        input = "list",
        output = "list",
        params = "list",
        wildcards = "list",
        threads = "numeric",
        log = "list",
        resources = "list",
        config = "list",
        rule = "character",
        bench_iteration = "numeric",
        scriptdir = "character",
        source = "function"
    )
)
snakemake <- Snakemake(
    input = list('output/analysis/bivariate_HEK293_WT_KO/with_Cond1/test_features.rds', 'output/analysis/bivariate_HEK293_WT_KO/with_Cond1/NMF.rds', "bigtable" = 'output/analysis/bivariate_HEK293_WT_KO/with_Cond1/test_features.rds', "nmf" = 'output/analysis/bivariate_HEK293_WT_KO/with_Cond1/NMF.rds'),
    output = list('output/analysis/bivariate_HEK293_WT_KO/with_Cond1/prediction', "data" = 'output/analysis/bivariate_HEK293_WT_KO/with_Cond1/prediction'),
    params = list('/prj/MiMB_book_chapter_Amina_Isabel/Nanopore/Amina/data/miCLIP_union_flat_exclude_Y_chromosome.bed', "mod" = '/prj/MiMB_book_chapter_Amina_Isabel/Nanopore/Amina/data/miCLIP_union_flat_exclude_Y_chromosome.bed'),
    wildcards = list('HEK293_WT_KO', 'Cond1', "var1" = 'HEK293_WT_KO', "var2" = 'Cond1'),
    threads = 1,
    log = list(),
    resources = list('tmpdir', "tmpdir" = '/tmp'),
    config = list("label" = 'HEK293_WT_KO', "path_jar" = 'JACUSA_v2.0.2-RC.jar', "picard_jar" = '/biosw/picard-tools/2.5.0/picard.jar', "jacusa_params" = list("p" = 16, "D" = '', "I" = '', "P1" = 'FR-SECONDSTRAND', "P2" = 'FR-SECONDSTRAND', "m" = 1, "q" = 1, "c" = 4, "a" = 'D,Y'), "java_params" = list("Xmx20g" = '', "XX:ParallelGCThreads=10" = ''), "nanosim_params" = list("threads" = 10, "fast_file" = '', "ref_genome" = '', "annotation" = '', "num_reads" = '', "seed" = '', "basecaller" = ''), "path_out" = 'output', "path_inp" = '../../data', "path_ref" = '/biodb/genomes/homo_sapiens/GRCh38_96/GRCh38_96.fa', "mod_file" = '/prj/MiMB_book_chapter_Amina_Isabel/Nanopore/Amina/data/miCLIP_union_flat_exclude_Y_chromosome.bed', "hg38" = '/prj/MiMB_book_chapter_Amina_Isabel/Nanopore/HEK293/JACUSA2/hg38.genome', "id" = c('1', '20', '2', '5'), "filter" = list("base" = 'A', "min_mis" = 1), "target" = c(), "data" = list("cond1" = c('HEK293T-WT-rep2', 'HEK293T-WT-rep3'), "cond2" = c('HEK293T-KO-rep2', 'HEK293T-KO-rep3')), "NMF_params" = list("NMF_permuation" = 30, "nmf_rank" = 3, "nrun" = 10, "rank_range" = c(2, 5), "thread" = 10), "traning" = NULL, "combination" = list("pattern" = c(1, 2, 3), "external_pattern" = '', "internal_pattern" = c('1', '2', '5')), "sampling_params" = list("seed" = c('0', '2'), "coverage" = c('300', '1000')), "mixing_params" = list("sampling_cov" = '1000', "sampling_seed" = '42', "mixing_thre" = c('0.1', '0.5'), "mixing_seed" = c('0', '2')), "analysis_params" = list("method" = 'LOF', "LOF_neighbors" = c(20, 60), "LOF_contamination" = 0.05)),
    rule = 'predict_mod',
    bench_iteration = as.numeric(NA),
    scriptdir = '',
    source = function(...){
        wd <- getwd()
        setwd(snakemake@scriptdir)
        source(...)
        setwd(wd)
    }
)


######## snakemake preamble end #########
library(NMF)
#BigTable=readRDS("BigTable.rds")
#res=readRDS("NMF.rds");
args = commandArgs(trailingOnly=TRUE)

library(ggplot2)                                        #replace with ggseqlogo

BigTable=readRDS(snakemake@input[[1]])
res=readRDS(snakemake@input[[2]]);
AllSites=as.matrix(BigTable[,1:15])
h = coef(res)
h = rbind(h, t(colSums(h)))
h = rbind(h, t(colSums(h[c(2,3,6,7),])))
h = rbind(h, t(colSums(h[c(2,3,6,7,4),])))

# print(coef(res))
AllSitesScores=AllSites%*%t(h)

Noverlap=read.delim(snakemake@params[[1]],header=F,as.is=T)
Noverlap.df=data.frame(ID=paste0(Noverlap[,1],":",Noverlap[,2]-2,"_",Noverlap[,3]+2,":",Noverlap[,6]),Type=Noverlap[,4])
rownames(Noverlap.df)=Noverlap.df[,1]

print(dim(AllSitesScores))
colnames(AllSitesScores)=c(paste("Pattern",1:nrow(h),sep=""))
AllSitesScores=data.frame(AllSitesScores[rownames(BigTable),],TotalScore=apply(AllSitesScores[rownames(BigTable),],1,sum),BigTable[,c("Motif","DRACH")])

if (dir.exists(snakemake@output[[1]])== FALSE) {
  dir.create(snakemake@output[[1]])
}
saveRDS(AllSitesScores,file=paste0(snakemake@output[[1]],"/ScoreProfile_NMFall_plusNonCLIP.rds"))
#AllSitesScores
# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

ins2= intersect(rownames(AllSitesScores),subset(Noverlap.df$ID,Noverlap.df$Type!="Boulias,Koertel,Koh"))
ins=intersect(rownames(AllSitesScores),Noverlap.df$ID);
NoCLIPsites=setdiff(rownames(AllSitesScores),ins)
print(length(ins2))
print(length(ins))
print(length(NoCLIPsites))

dataGG=data.frame(AllSitesScores[c(ins,NoCLIPsites),-c(ncol(AllSitesScores)-1,ncol(AllSitesScores))],CLIP=c(as.character(Noverlap.df[ins,2]),rep("NoCLIP",length(NoCLIPsites))))
NumberOfFactors=ncol(AllSitesScores)-2
colnames(dataGG)=c(paste("Pattern",1:NumberOfFactors,sep=""),"CLIP")

# rocca<-data.frame(NMF2=dataGG[,"NMF2"],CLIP=c(rep(1,length(ins)),rep(0,length(NoCLIPsites))))
library(ggplot2)
# require(plotROC)
# p<-ggplot(rocca, aes(m = NMF2, d = CLIP)) + geom_roc(labels=F)
# #color = Base_1
# p<-p+geom_abline(slope=1,intercept=0)
# p
# ggsave(paste0(args[1],"roc.pdf"))
# #ok, kinda works
# calc_auc(p)

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
#   ggtitle(dude)
ggsave(paste0(snakemake@output[[1]],"/",dude,"_ppv.pdf"),device="pdf")
# # coordinates <- coords(r, x = "all", input = "threshold", ret = c("threshold","ppv","tpr","tp", "tn","fp","fn"))
# # coordinates$ppv <- coordinates$ppv*100 
# # coordinates$tpr <- coordinates$tpr*100 
# # ggplot(coordinates, aes(x=tpr)) +  
# # #   geom_line( aes(y=npv), size=2, color="#009E73") + 
# #   geom_line( aes(y=ppv), size =2,color="#D55E00") +
# #   scale_y_continuous(    
# #     # Features of the first axis
# #     name = "Predicted Targets Rate",
# #     # Add a second axis and specify its features
# # #     sec.axis = sec_axis(~.,name="Positive Predictive Value")
# #   ) +  scale_x_continuous(name = "NPV" )+  theme_bw() + theme(
# #     axis.title.x = element_text(size=rel(1.5)),      
# #     axis.title.y = element_text(color = "#009E73", size=rel(2)),
# # #     axis.title.y.right = element_text(color = "#D55E00", size=rel(1)), 
# #     axis.text = element_text(size = rel(1.5))

# #   ) 
# #   ggtitle(dude)
# # ggsave(paste0(args[1],"/",dude,"_npv.pdf"),device="pdf")
    
}
# #19:3976465-3976470:- ##