import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
from lib import *


downsampling_thres = snakemake.params[0]     # List of downsampling Read  Coverage 
target = snakemake.params[1]    # Target site
seeds = snakemake.params[2]       # List of samtools seeds for downsampling and mixing
ref = snakemake.params[3]           # Sequence ID
mixing_thres = snakemake.params[4]                 # List of modification rates
mixing_downsampling_seed = snakemake.params[5]    # Downsampling seed for the mixing samples
mixing_downsampling_thre = snakemake.params[6]   # Downsampling covcerage for the mixing samples
output =  snakemake.output[0]   # output folder
if not os.path.exists(output):
    os.makedirs(output, exist_ok=True)
output = output + "/"
label = snakemake.params[7]     # Analysis label
method = snakemake.params[8]      # Outlier detection method
contamination = str(snakemake.params[9])  # LOF contamination value
path_inp = snakemake.params[10]     #Path to inputs
features = snakemake.params.features
baseline= 'median'
mix_seeds = snakemake.params.mix_seed       # List of samtools seeds for downsampling and mixing


#Generate barplots showing the diffrence between th target position and the median across different seeds of downsampling in terms of LOF scores.
for feat in features:
    for sl in target: 
        ar =[]
        df =[]
        state = True
        for thre in downsampling_thres: 
                tab = []
                for ij in seeds:
                    inp1 = path_inp+ '.sampled'+thre+'/DowS'+ij+'/with_Cond1/'+method+contamination+'/'+ feat+'_lof_scores.csv'
                    t = pd.read_csv(inp1, sep = ',',header = 0)
                    maxi = t.sort_values('scores').reset_index(drop= True).loc[0,'scores']
                    if baseline == 'median':
                        median = t['scores'].median()
                        val = t.loc[(t.pos - sl)==0,'scores'].values[0]
                        tab.append((median - val)/abs(maxi))
                    else:
                        if state == True: 
                            arr = t.sort_values('scores').reset_index(drop= True)
                            ref = arr[abs(arr.pos - sl) > 4].reset_index(drop = True).loc[0,'scores']
                            pos = arr[abs(arr.pos - sl) > 4].reset_index(drop = True).loc[0,'pos']
                            state = False
                        val = t.loc[(t.pos - sl)==0,'scores'].values[0]
                        tab.append((ref - val)/abs(maxi))
                df.append([int(thre),np.mean(tab),np.std(tab)] )  
        ar= pd.DataFrame(df, columns=['cov', 'scores', 'std'])
        title = "LOF_"+feat+"_Site:"+str(sl)
        plot(ar = ar, x_values= 'cov',y_values ="scores", xerr_values = 'std',xlabel = 'Read Coverage',  ylabel ='Distance', title = title, path = output)
        
#Generate barplots showing the diffrence between th target position and the median in terms of combined features from JACUSA2 estimated scores.
for feat in features:
    for sl in target: 
        ar =[]
        df =[]
        for thre in downsampling_thres: 
                tab = []
                for ij in seeds:
                    inp1 = path_inp+ '.sampled'+thre+'/DowS'+ij+'/with_Cond1/'+method+contamination+'/prediction.out'         
                    t = pd.read_csv(inp1, sep = ',',header = 0)
                    t =t.loc[t.label == 'Cond1_Cond2', :]
#                     t =t.loc[t.label == label, :]
                    median = t[feat].median()
                    val = t.loc[(t.Ref_Pos == ref+ "_"+str(sl)),feat].values[0]
                    maxi = t.sort_values(feat,ascending = False ).reset_index(drop= True).loc[0,feat]
                    mini = t.sort_values(feat,ascending = True ).reset_index(drop= True).loc[0,feat]
                    tab.append((val - median)/(maxi-mini))
                
                df.append([int(thre),np.mean(tab),np.std(tab)] )  
        ar= pd.DataFrame(df, columns=['cov', feat, 'std'])
        title = "JACUSA_"+feat+"_Site:"+str(sl)
        plot(ar = ar, x_values= 'cov',y_values = feat, xerr_values = 'std', xlabel = 'Read Coverage', ylabel = 'Distance', title = title, path= output)
#Generate barplots showing the fraction of detection of target position and its neighbors as outleirs across different seeds of downsamplings for different amount of reads.

for feat in features:
    for sl in target: 
        ar =[]
        df =[]
        for thre in downsampling_thres: 
                tab = pd.DataFrame([])
                for ij in seeds:
                    inp1 = path_inp+ '.sampled'+thre+'/DowS'+ij+'/with_Cond1/'+method+contamination+'/'+ feat+'_lof_scores.csv'
                    t = pd.read_csv(inp1, sep = ',',header = 0)
                    t =t.loc[t.outlier == True, :]
                    zz = t.loc[(t.pos - sl)== 0,:]
                    qq = t.loc[abs(t.pos - sl) < 3,:]
                    tab =tab.append([[qq.shape[0]/t.shape[0],zz.shape[0],qq.shape[0]!=0]])
                df.append([thre,np.mean(tab[0]),np.mean(tab[1]),np.mean(tab[2])] )  
        ar= pd.DataFrame(df, columns=['thre','5mer outliers fraction', 'site '+str(sl), '5mer'])
        title = "Downsapling Detection Frequency_"+feat+"_Site:"+str(sl)
        plot(ar, x_values='thre',y_values =['5mer', 'site '+str(sl)], xlabel = 'Read Coverage', ylabel = "Fraction detected as outliers", title = title, path = output, legend = True)

#Generate barplots showing the fraction of detection of target position and its neighbors as outleirs across different seeds of mixing for different mixing rates.
mix_seeds = mix_seeds[mix_seeds!=mixing_downsampling_seed]
for feat in features:
    for sl in target: 
        ar =[]
        df =[]
        for thre in mixing_thres: 
                tab = pd.DataFrame([])
                for ij in seeds:
                    inp1 = path_inp+'.sampled'+mixing_downsampling_thre+'/DowS'+mixing_downsampling_seed+'/with_MixS'+ij+'/'+str(thre)+'_Cond1/'+method+contamination+'/'+feat+'_lof_scores.csv'         
                    t = pd.read_csv(inp1, sep = ',',header = 0)
                    t =t.loc[t.outlier == True, :]
                    zz = t.loc[(t.pos - sl)== 0,:]
                    qq = t.loc[abs(t.pos - sl) < 3,:]
                    tab =tab.append([[qq.shape[0]/t.shape[0],zz.shape[0],qq.shape[0]!=0]])
                df.append([thre,np.mean(tab[0]),np.mean(tab[1]),np.mean(tab[2])] )  
        ar= pd.DataFrame(df, columns=['thre','5mer outliers fraction', 'site '+str(sl), '5mer'])
        title = "Mixing Detection Frequency_"+feat+"_Site:"+str(sl)
        plot(ar, x_values='thre',y_values =['5mer', 'site '+str(sl)], xlabel = 'Modification Rate', ylabel = "Fraction detected as outliers", title = title, path = output, legend = True)


