library(NMF)
BigTable=readRDS("BigTable.rds")
res=readRDS("NMF.rds");
##
pdf("NMF_6cluster_new.pdf")
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
pdf(paste0("Pattern_",k,"_barplot_NMF.pdf"),width=14,height=7)
rep1=matrix(h[k,],ncol=5,byrow=T)[1:3,]
rownames(rep1)=c("BASE","DEL","INS")
colnames(rep1)=c("Pos1","Pos2","Pos3","Pos4","Pos5")
barplot(rep1,legend=F, main="Experiment 1",col=cbbPalette,,cex.names=2,cex.axis=2)
dev.off()
}


wExtended=data.frame(w,BigTable[rownames(w),"Motif"])
print("Table of number of best hits for each pattern")  
print(table(apply(wExtended[,1:nrow(h)],1,function(x){which(x==max(x))})))

tt=apply(wExtended[,1:nrow(h)],1,function(x){which(x==max(x))})

                                        #BigTable[names(which(tt==k)),"Motif"]
                                        #BigTable[names(which(tt!=k)),"Motif"]

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

                                        #motif search TODO
library(Logolas)
for(k in 1:nrow(h))
    {
        pdf(paste0("Logo_",k,"_barplot_NMF.pdf"),width=7,height=7)
        logomaker(BigTable[names(which(tt==k)),"Motif"], type = "Logo")
        logomaker(apply(CountMatrix(BigTable[names(which(tt==k)),"Motif"]),2,function(x){x/sum(x)}), type = "EDLogo", bg=apply(CountMatrix(BigTable[names(which(tt!=k)),"Motif"]),2,function(x){x/sum(x)}))
        #apply(CountMatrix(BigTable[names(which(tt!=k)),"Motif"]),2,function(x){x/sum(x)})
        dev.off()
        }

