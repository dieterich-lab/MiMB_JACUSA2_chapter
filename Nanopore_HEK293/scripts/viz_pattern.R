library(NMF)
print("Vizualize patterns ...")

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


if (dir.exists(snakemake@output[[1]])== FALSE) {
      dir.create(snakemake@output[[1]])
    }

#load data and patterns
res<-readRDS(snakemake@input[[1]]);
w<-basis(res)
print(dim(w))
h <- coef(res)
h <- rbind(h, t(colSums(h)))
#add the combined patterns
for (pt in snakemake@params[[1]]){
    h <- rbind(h, t(colSums(h[pt,])))
    }
print(dim(h))

# barplots of patterns
for(k in 1:nrow(h))
    {
    pdf(paste0(snakemake@output[[1]],"/pattern_",k,"_barplot_NMF.pdf"),width=7,height=7)
    rep1<-matrix(h[k,],ncol=5,byrow=T)[1:3,]
    rownames(rep1)<-c("Mismatch","Deletion","Insertion")
    colnames(rep1)<-c("Pos1","Pos2","Pos3","Pos4","Pos5")
    barplot(rep1,main=paste0("NMF_Pattern_",k),col=cbbPalette,cex.names=2,cex.axis=2, legend = TRUE)
    dev.off()
    }
