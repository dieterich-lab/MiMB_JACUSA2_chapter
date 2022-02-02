library(NMF)
args <- commandArgs(trailingOnly=TRUE)
res<-readRDS(args[1])
colnames(res)<- c("Mismatch_Pos1","Mismatch_Pos2","Mismatch_Pos3","Mismatch_Pos4","Mismatch_Pos5","Deletion_Pos1","Deletion_Pos2","Deletion_Pos3","Deletion_Pos4","Deletion_Pos5","Insertion_Pos1","Insertion_Pos2","Insertion_Pos3","Insertion_Pos4","Insertion_Pos5")
pdf(args[2])
layout(cbind(1,2))
# basis components
basismap(res)
# mixture coefficients
aga<-coefmap(res)
dev.off()
#Choice of best pattern for subsequent analysis
w<-basis(res)
print("Table of number of best hits for each pattern")
TabPatternInstances<-table(apply(w,1,function(x){which(x==max(x))}))
print(TabPatternInstances)
pdf(args[3],width=7,height=7)
barplot(height=TabPatternInstances, names=1:(dim(w)[2]), main="NMF_Patterns_Scoring",xlab="Patterns", 
        ylab="Membership Score (basis matrix)",ylim=range(pretty(c(0, TabPatternInstances))))
dev.off()

