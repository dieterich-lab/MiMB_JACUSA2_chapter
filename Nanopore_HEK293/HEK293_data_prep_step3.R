library(NMF)
print("Vizualize patterns ...")
if (dir.exists(snakemake@output[[1]])== FALSE) {
  dir.create(snakemake@output[[1]])
}
BigTable=readRDS(snakemake@input[[1]])
res=readRDS(snakemake@input[[2]]);
##

colnames(res)= c("Mismatch_Pos1","Mismatch_Pos2","Mismatch_Pos3","Mismatch_Pos4","Mismatch_Pos5","Deletion_Pos1","Deletion_Pos2","Deletion_Pos3","Deletion_Pos4","Deletion_Pos5","Insertion_Pos1","Insertion_Pos2","Insertion_Pos3","Insertion_Pos4","Insertion_Pos5")
pdf(paste0(snakemake@output[[1]],"/NMF_matrices.pdf"))
layout(cbind(1,2))
# basis components
basismap(res)
# mixture coefficients
aga<-coefmap(res)
dev.off()

w<-basis(res)
h <- coef(res)
h = rbind(h, t(colSums(h)))
for (pt in snakemake@params[[2]]){
h = rbind(h, t(colSums(h[pt,])))
}
# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

for(k in 1:nrow(h))
    {
pdf(paste0(snakemake@output[[1]],"/Pattern_",k,"_barplot_NMF.pdf"),width=7,height=7)
rep1=matrix(h[k,],ncol=5,byrow=T)[1:3,]
rownames(rep1)=c("Mismatch","Deletion","Insertion")
colnames(rep1)=c("Pos1","Pos2","Pos3","Pos4","Pos5")
barplot(rep1,main=paste0("NMF_Pattern_",k),col=cbbPalette,cex.names=2,cex.axis=2, legend = TRUE)
dev.off()
}

                                        #Choice of best pattern for subsequent analysis

wExtended=data.frame(w,BigTable[rownames(w),"Motif"])
print("Table of number of best hits for each pattern")
TabPatternInstances=table(apply(wExtended[,1:(dim(w)[2])],1,function(x){which(x==max(x))}))
print(dim(TabPatternInstances))

pdf(paste0(snakemake@output[[1]],"/pattern_scores.pdf"),width=7,height=7)
barplot(height=TabPatternInstances, names=1:(dim(w)[2]), main="NMF_Patterns_Scoring",xlab="Patterns", 
        ylab="Membership Indicator (basis matrix)",ylim=range(pretty(c(0, TabPatternInstances))))
dev.off()

