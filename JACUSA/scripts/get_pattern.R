library(NMF)

plotPattern <- function(coeff,path_out)
{ colpalette <- c( "#F58634","#FFCC29", "#00AF91")
    print(dim(coeff))
    pdf(path_out,width=5,height=20)
    split.screen(c(dim(coeff)[1],1))        
    for (i in 1:dim(coeff)[1]){
    pat=matrix(coeff[i,],ncol=5,byrow=F)
       if (dim(coeff)[2]>5){rownames(pat)=c("Mis","Ins","Del")
       }else{rownames(pat)=c("Mis")}
        
    colnames(pat)=c("Pos1","Pos2","Pos3","Pos4","Pos5")
    screen(i)
    barplot(pat,legend=T, main=paste("Pattern ",i),col=colpalette)

}
    dev.off()
}
path_out=snakemake@output[[1]]
data_matrix = read.csv(snakemake@input[[1]])
rank = snakemake@params[[1]]
for (lab in unique(data_matrix$label))
{
data_matrix_tmp = data_matrix[data_matrix$label == lab,]
data_matrix_tmp=data_matrix_tmp[,which(colnames(data_matrix_tmp)=="Mis_x_1"):dim(data_matrix_tmp)[2]]
data_matrix_tmp=data_matrix_tmp[, colSums(abs(data_matrix_tmp)) != 0]
data_matrix_tmp = data_matrix_tmp +  abs(min(data_matrix_tmp))
fit_NMF= nmf(data_matrix_tmp,rank, nrun=snakemake@params[[2]],  seed=42,.options=paste('vp',snakemake@params[[3]],sep = ""))
plotPattern(coef(fit_NMF),gsub('Cond1_Cond2', lab,path_out))
write.csv(coef(fit_NMF),gsub('pdf','csv',gsub('Cond1_Cond2', lab,path_out)), row.names = FALSE)
}
