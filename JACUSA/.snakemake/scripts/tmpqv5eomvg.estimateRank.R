
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
    input = list('output/snakemake/analysis/tetR_metvstetR/Features.out'),
    output = list('output/snakemake/analysis/tetR_metvstetR/NMFestimateRank_Statistics.pdf'),
    params = list(30, 30, 6, c(2, 4), c(2, 4)),
    wildcards = list(),
    threads = 1,
    log = list(),
    resources = list('tmpdir', "tmpdir" = '/tmp'),
    config = list("label" = 'tetR_metvstetR', "path_jar" = '/prj/Isabel_ONT_rRNA/software/JACUSA_v2.0.0-RC22.jar', "picard_jar" = '/biosw/picard-tools/2.5.0/picard.jar', "jacusa_params" = list("p" = 16, "D" = '', "I" = '', "P1" = 'FR-SECONDSTRAND', "P2" = 'FR-SECONDSTRAND'), "path_out" = 'output/snakemake', "path_inp" = '/prj/Claudia_Hoebartner_WUE/tetR_and_blaR/Amina/data/tetR_blaR', "path_ref" = '/prj/Claudia_Hoebartner_WUE/tetR_and_blaR/Amina/data/tetR_blaR/Wuerzburg_number1.fa', "data" = list("cond1" = c('Ribozyme_teR-met_blaR'), "cond2" = c('tetR-I_blaR-met')), "NMF_params" = list("NMF_permuation" = 30, "nmf_rank" = 5, "nrun" = 30, "rank_range" = c(2, 4), "thread" = 6)),
    rule = 'estimate_rank',
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
estimateRank <- function(data, nb_permutation = 10, interval = seq(2,9), method='snmf/r', nrun=30) {
  cc= matrix(0,nb_permutation,length(interval))
    for (j in 1:nb_permutation)
    { nmf_random= nmf(randomize(data), interval, nrun=nrun,.options=paste('vp',np), seed=42)
       cc[j,]= nmf_random$measures$cophenetic
    }
    nmf_true = nmf((data),interval, nrun=nrun,.options='vp30' , seed=42)
    pdf("NMFestimateRank_Statistics.pdf")
    plot(nmf_true,nmf_random)
    dev.off()
    
    val = nmf((data),interval, nrun=nrun)$measures$cophenetic
    m = colMeans(cc)
    sdd= apply(cc, 2, sd)
    error=1.96* sdd/sqrt(nb_permutation)
    inf = m-error
    sup = m+error

    a= 0.5
    b=1.5 
    pdf("NMFestimateRank_CophCorr.pdf")
    boxplot(cc,col='#F58634', names=c(interval), xlab='Rank',ylab="Cophenetic Correlation", ylim = c(min(cc)-0.05, 1) , frame.plot = FALSE)
    len = len(interval)
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
nr  = snakemake@params[[2]] 
np= snakemake@params[[3]]
data_matrix = read.csv(snakemake@input[[1]])
data_matrix = data_matrix +  abs(min(data_matrix))
if (snakemake@params[[1]] >29)
{estRank = estimateRank(data_matrix, nb_permutation = snakemake@params[[1]], interval = seq(2,9), nrun = nr )}
else 
    {nmf_true = nmf((data),interval, nrun=nr,.options=paste('vp',np) , seed=42)
     nmf_random= nmf(randomize(data), interval, nrun=nr,.options=paste('vp',np), seed=42)
     pdf("NMFestimateRank_Statistics.pdf")
     plot(nmf_true,nmf_random)
     dev.off()}

