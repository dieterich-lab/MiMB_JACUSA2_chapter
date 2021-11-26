import pandas as pd
import numpy as np
from lib import *
if not os.path.exists(snakemake.output[0]):
    os.makedirs(snakemake.output[0], exist_ok=True)
out = snakemake.output[0] + "/"
ref = snakemake.params.ref
title = snakemake.params.lab+"_"+ref
describe = pd.DataFrame([])
features_df = pd.read_csv(snakemake.input[0], sep = ',')
table = features_df.loc[:,['label','Ref_Pos','Ref','Pos', 'ModStatus','Status']]
try: 
           pattern = pd.read_csv(combination['external_pattern'], sep = ',')                 
except Exception:
           pattern = pd.read_csv(snakemake.input[1], sep = ',')        
for f in  snakemake.params.combination['pattern']:
            p = pattern.iloc[f,:]
            tmp=np.zeros((features_df.shape[0]))
            
            for colname in pattern.columns:
                tmp = tmp+p[colname]*features_df[colname]
            table['Pattern_'+str(f)]=tmp
patterns_names = list(map('Pattern_'.__add__,map(str, snakemake.params.combination['pattern'])))
if snakemake.params.analysis_type == "bivariate_":
    for feature in patterns_names:
                df_= table[table["Ref"] == ref ].reset_index(drop=True)
                res = BarPlot(df_[['Pos',feature, 'ModStatus','Status']],feature,title,snakemake.params.target,outlier=snakemake.params.lof_thre,neigh=snakemake.params.lof_neigh,path = out, method = snakemake.params.meth)
                describe=pd.concat([describe,res], axis =1)
else:

    for feature in patterns_names:
        label1 = table['label'].unique()[0]
        label2 = table['label'].unique()[1]
        label3 = table['label'].unique()[2]

        x= table[(table["label"] == label1) & (table["Ref"] == ref)]
        y= table[(table["label"] == label2) & (table["Ref"] == ref)]
        z= table[(table["label"] == label3) & (table["Ref"] == ref)]

        df_ = pd.merge(x,y, on = 'Ref_Pos')
        df_ = pd.merge(df_,z, on = 'Ref_Pos').reset_index(drop=True)

        res = ScatterPlot(df_[['Pos',feature+'_x',feature+'_y',feature]],feature,title,snakemake.params.target,outlier=snakemake.params.lof_thre,neigh=snakemake.params.lof_neigh,path = out)
        describe=pd.concat([describe,res], axis =1)
describe.columns = patterns_names
describe.to_csv(out+'summary.out')
table.to_csv(out+'prediction.out')

