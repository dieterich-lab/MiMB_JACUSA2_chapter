library(NMF)
estimateRank <- function(data, nb_permutation = 40, interval = seq(2,9), method='snmf/r', nrun=30, label= '',path ='') {
  cc= matrix(0,nb_permutation,length(interval))
    for (j in 1:nb_permutation)
    { nmf_random= nmf(randomize(data), interval, nrun=nrun,.options=paste('vp',np,sep = ""), seed=42)
       cc[j,]= nmf_random$measures$cophenetic
    }

    val = nmf(data,interval, nrun=nrun)$measures$cophenetic
    m = colMeans(cc)
    sdd= apply(cc, 2, sd)
    error=1.96* sdd/sqrt(nb_permutation)
    inf = m-error
    sup = m+error

    a= 0.5
    b=1.5 
    pdf(paste(path_output,"/",label,"_NMFestimateRank_CophCorr.pdf",sep = ""))
    boxplot(cc,col='#F58634', names=c(interval), xlab='Rank',ylab="Cophenetic Correlation", ylim = c(min(cc)-0.05, 1) , frame.plot = FALSE)
    len = length(interval)
    for (i in 1:len)
    {
      lines(c(a,b),c(inf[i],inf[i]),col="#00AF91")
      lines(c(a,b),c(m[i],m[i]),col="#FFCC29",lwd=2)
      lines(c(a,b),c(sup[i],sup[i]),col="#00AF91")
      lines(c(a,b),c(val[i],val[i]),col='#AC0D0D',lwd=2)
      a= a+1
      b= b+1
    }
    legend("bottomleft", c("True_CC","95% CI", "Mean"), lty=1,col = c(2,4,3),bty ="n")
   dev.off()
    
  return(cc)
  }


data_matrix = read.csv(snakemake@input[[1]])

path_output = snakemake@output[[1]]

if (dir.exists(path_output)== FALSE) {
  dir.create(path_output)
}
nr  = snakemake@params[[2]] 
np= snakemake@params[[3]]
min_rank = snakemake@params[[4]]
max_rank = snakemake@params[[5]]

for (lab in unique(data_matrix$label))
{
data_matrix_tmp = data_matrix[data_matrix$label == lab,]
data_matrix_tmp=data_matrix_tmp[,which(colnames(data_matrix_tmp)=="Mis_x_1"):dim(data_matrix_tmp)[2]]
data_matrix_tmp=data_matrix_tmp[, colSums(abs(data_matrix_tmp)) != 0]
data_matrix_tmp = data_matrix_tmp +  abs(min(data_matrix_tmp))
if(snakemake@params[[1]] >0)
{estRank = estimateRank(data_matrix_tmp, nb_permutation = snakemake@params[[1]], interval = seq(min_rank,max_rank), nrun = nr,label=lab, path=path_output )}
nmf_true = nmf(data_matrix_tmp,min_rank:max_rank, nrun=nr,.options=paste('vp',np,sep = "") , seed=42)
nmf_random= nmf(randomize(data_matrix_tmp), min_rank:max_rank, nrun=nr,.options=paste('vp',np,sep = ""), seed=42)
pdf(paste(path_output,"/",lab,"_NMFestimateRank_Statistics.pdf",sep = ""))
plot(nmf_true,nmf_random)
dev.off()
}
