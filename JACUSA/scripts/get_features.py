import pandas as pd
import numpy as np
from lib import * 
df3=pd.DataFrame()
try:
    mods = pd.read_csv(snakemake.params[0], sep = ',')
    mods['Ref_Pos'] = mods[mods.columns[0]]+ "_" + mods[mods.columns[1]].astype(str)
except Exception:
    pass 
labels = ['Cond1_Cond2','Cond2_Cond3','Cond1_Cond3']
for i in range(len(snakemake.input)):
            inp1= snakemake.input[i]
            df0 = pd.read_csv(inp1, sep = '\t',skiprows=1)
            df1 = ExtractFeatures(df0).reset_index(drop = True)
            df1['Ref_Pos'] = df1["Ref"]+ "_" + df1["Pos"].astype(str) 
            try:
                df1 = pd.merge(df1,mods[['Ref_Pos','ModStatus','Status']], on='Ref_Pos')  
            except Exception:
                df1['ModStatus'] = 'Unm'
                df1['Status'] = 'Unm'
            df2 = KmerFeatures(df1)
            df2.insert(0, "label", labels[i])
            df3 = df3.append(df2)
df3.to_csv(snakemake.output[0], index=False)