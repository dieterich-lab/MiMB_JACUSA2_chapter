import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.neighbors import LocalOutlierFactor
import os
from hyperopt import hp
from hyperopt import fmin, tpe
from scipy import stats
import hyperopt
"""
Main Functions

* ExtractFeatures(...)
* KmerFeatures(...)
* BarPlot(...)
* ScatterPlot(...)
"""
def ExtractFeatures(df) :      
    """
    ExtractFeatures function extracts features from JACUSA2 CALL2 output as a dataframe columns: 
    'Ref','Pos','Base','strand','Cov1','Cov2','Min_cov','Mis','Ins','Del','dA','dC',
    'dG','dT','A1','C1','G1','T1','A2','C2','G2','T2','Ins11','Ins21','Del11','Del21'

    :param df: output of JACUSA2 Call2
    :return: a dataframe of columns representing features: Reference, Base, coverage from both conditions,
    min coverage, Base frequencies, Mismatch, deletion and insertion scores and counts
    """ 
    # select features from the JACUSA2 CALL 2 output
    df=df[df['strand']=='+']
    df2= df.rename(columns={"#contig": "Ref", "end": "Pos", "score":"Mis", 'ref':'Base'})
    df2 =  df2[df2['Mis']!='*']
    df2['Mis'] =  df2['Mis'].astype(float)
    tmp=np.zeros((df2.shape[0],6))
    try: 
        df2[['A2','C2','G2','T2']] = df['bases21'].str.split(',',expand=True).astype(float)
        df2["Cov2"]= df2["A2"] +df2["C2"]+df2["G2"]+df2["T2"]
        df2[['A1','C1','G1','T1']] = df['bases11'].str.split(',',expand=True).astype(float)
        df2["Cov1"]= df2["A1"] +df2["C1"]+df2["G1"]+df2["T1"]
        df2["Min_cov"]=pd.DataFrame([df2['Cov1'], df2['Cov2']]).min()
    except Exception:
        df2[['A1','C1','G1','T1']] = df['bases11'].str.split(',',expand=True).astype(float)
        df2["Cov1"]= df2["A1"] +df2["C1"]+df2["G1"]+df2["T1"]
#     df2[["A1","C1","G1","T1"]]=df2[["A1","C1","G1","T1"]].div(df2.Cov1, axis=0)
#     df2[["A2","C2","G2","T2"]]=df2[["A2","C2","G2","T2"]].div(df2.Cov2, axis=0)
#     df2[["dA","dC","dG","dT"]] = abs(df2[["A1","C1","G1","T1"]].values - df2[["A2","C2","G2","T2"]].values)

    # extract features from info field
    for idx in range(df2.shape[0]):
        strg = df2['info'].iloc[idx].replace(';resetP=0','')
        if strg.find("deletion_score") != -1:
                interm = strg[strg.find("deletion_score=")+len("deletion_score="):]
                if interm.find(";") != -1:
                    tmp[idx][0] = float(interm[0 : interm.find(";")])
                else:
                    tmp[idx][0] = float(interm)
        if strg.find("insertion_score") != -1:
                interm = strg[strg.find("insertion_score=")+len("insertion_score="):]
                if interm.find(";") != -1:
                    tmp[idx][1] = float(interm[0 : interm.find(";")])
                else:
                    tmp[idx][1] = float(interm)
        if strg.find("ins11") != -1:
                interm = strg[strg.find("ins11=")+len("ins11="):]
                tmp[idx][2] = float(interm[0 : interm.find(",")])

        if strg.find("ins21") != -1:
                interm = strg[strg.find("ins21=")+len("ins21="):]
                tmp[idx][3] = float(interm[0 : interm.find(",")])

        if strg.find("del11") != -1:
                interm = strg[strg.find("del11=")+len("del11="):]
                tmp[idx][4] = float(interm[0 : interm.find(",")])

        if strg.find("del21") != -1:
                interm = strg[strg.find("del21=")+len("del21="):]
                tmp[idx][5] = float(interm[0 : interm.find(",")])

    df2.loc[:,['Del','Ins','Ins11','Ins21','Del11','Del21']]=np.nan_to_num(tmp)
    try:
        df = df2[['Ref', 'Pos', 'Base', 'strand','Cov1','Cov2','Min_cov',"Mis","Ins","Del","A1","C1","G1","T1","A2","C2","G2","T2",'Ins11','Ins21','Del11','Del21']]
    except Exception:
        df = df2[['Ref', 'Pos', 'Base', 'strand','Cov1',"Mis","Ins","Del","A1","C1","G1","T1",'Ins11','Ins21','Del11','Del21']]
    return df

def KmerFeatures(df,K=5):
    """
    KmerFeatures function extracts features from JACUSA2 CALL2 output in Kmer context as a dataframe columns: 
    'Mis_x_$k','Ins_x_$k','Del_x_$k','5mer'

    :param df: JACUSA2 Call2 extracted features 
    :param K: mer size
    :return: a dataframe of features in Kmer context
    """     
    size = df.shape[0]
    for k in range(K):
            df[['Mis_x_'+str(k+1),'Ins_x_'+str(k+1),'Del_x_'+str(k+1)]] =np.zeros((size,3))

    # Add Mis, Ins and Del in 5mer context
    for i in range(size):
        dff = df[df['Ref'] == df.loc[i,'Ref']]
        dff2 = dff[dff['Pos'].isin(range(dff.loc[i,'Pos']-4,dff.loc[i,'Pos']+5))]
        for jj in range(K):
             dff3 = dff2[dff2['Pos'].isin(range(dff2.loc[i,'Pos']-K+jj+1,dff2.loc[i,'Pos']+jj+1))]
             df.loc[i,'Mis_x_'+str(K-jj)]= dff3['Mis'].sum()   
             df.loc[i,'Ins_x_'+str(K-jj)]= dff3['Ins'].sum()   
             df.loc[i,'Del_x_'+str(K-jj)]= dff3['Del'].sum()   
    return df

def LOF(outlier,  rangeOut, val,filterPos):
    def objective(n):
        lof = LocalOutlierFactor(contamination =outlier, n_neighbors = n)
        ohat = lof.fit_predict(val)
        return lof.negative_outlier_factor_[filterPos].min()

    space = hp.choice('n', range(rangeOut[0],rangeOut[1]))
    best = fmin(objective, space, algo=tpe.suggest, max_evals=50)
    return hyperopt.space_eval(space, best)


def BarPlot(df_,feat,title,target,  indx=[],outlier = 0.001,neigh=[10,50], path ='',method = 'IQR'):
    """
    BarPlot function perform bar plot of one features(score) from JACUSA2 CALL2 output

    :param df_: a dataframe of features including the score to plot  
    :param feat: feature's name
    :param title: title of the barplot
    :param indx: outliers indices, if provided, no need to do the outlier detection
    :param outlier: outlier contamination for the LocalOutlierFactor method
    :param neigh: range of neighborhood size for the LocalOutlierFactor method    
    :param zoom: provided if we want to zoom in a specific rRNA squence region 
    :param path: path to save plots as eps/pdf format
    :param method: method used to detect outliers : IQR or LOF
    
    """     
    val = df_[feat]
    pos =df_['Pos']
    filterPos = val > val.median()
    idx_bool=np.full((len(val),), False, dtype=bool)
    status = np.full((len(val),), False, dtype=bool)
    max_lof_score = 0
    result_lof = []
    # outliers indices
    if len(indx)>0:
        # outliers are given
        indices = indx
    else: 
        if method == 'LOF':
           # perform outlier detection with LocalOutlierFactor method
            n_best = LOF(outlier,neigh,val.values.reshape(-1, 1),filterPos)
            print(n_best)
            lof = LocalOutlierFactor(contamination =0.5, n_neighbors = n_best)
            ohat = lof.fit_predict(val.values.reshape(-1, 1))
            scores = lof.negative_outlier_factor_
            result_lof=pd.DataFrame([scores]).T.describe()
            result_lof.loc['median',:]=[np.median(scores)]            
            result_lof.loc['skew',:]=[stats.skew(scores)]
            result_lof.loc['kurtosis',:]=[stats.kurtosis(scores)]
            sortedarray = np.argsort(scores)
            nb = round(outlier*len(val))
            elemts = np.where(val < val.median())[0]
            i=0
            j=0
            while i < nb:
                if sortedarray[j] in elemts:
                    idx_bool[sortedarray[j]]=True
                    status[sortedarray[j]]=False
                else:
                    idx_bool[sortedarray[j]]=True
                    status[sortedarray[j]]=True 
                    if i==0:
                        result_lof.loc['min_modified',:]=[scores[sortedarray[j]]]                        
                    i= i+1
                j=j+1
            outl = df_.iloc[status,:]
#             idx_bool = ohat == -1
            result_lof.loc['context_size',:]=[n_best]
            lofscores = pd.DataFrame({'pos': pos,'scores':lof.negative_outlier_factor_, 'outlier':idx_bool, 'modified':status})
            lofscores.to_csv(path+feat+'_lof_scores.csv', index = False)
        else: 
        # Using IQR method, extract the upper and lower quantiles
            pokemon_HP_lq = df_[feat].quantile(0.25)
            pokemon_HP_uq = df_[feat].quantile(0.75)
            #extract the inter quartile range
            pokemon_HP_iqr = pokemon_HP_uq - pokemon_HP_lq#get the upper and lower bounds
            lower_bound = pokemon_HP_lq - 1.5*pokemon_HP_iqr
            upper_bound = pokemon_HP_uq + 1.5*pokemon_HP_iqr#extract values outside these bounds 
            outl = df_[(df_[feat] <= lower_bound) | (df_[feat] >= upper_bound)]
            
                
    # barplots
    _, ax = plt.subplots(figsize=(20,10), dpi = 300)
    ax.stem(pos,val, 'silver',markerfmt=" ", label = 'Inliers')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # Highlight each classe of outliers with a diffrent color : modfied, neighbors and unmodified positions
    ind = outl[outl['ModStatus' ] != 'Unm'].index
    if len(ind) != 0:
                m, s, b = ax.stem(pos[ind],val[ind], markerfmt='.', label = 'Outlier & Modified')
                plt.setp([m, s], color='#2CBDFE')
    ind = outl[outl['Status' ].isna()].index
    if len(ind) != 0:
                m, s, b = ax.stem(pos[ind],val[ind], markerfmt=".", label = 'Outlier & Neighbor')
                plt.setp([m, s], color='#47DBCD')
    ind = outl[outl['Status' ] == 'Unm'].index
    if len(ind) != 0:
                m, s, b = ax.stem(pos[ind],val[ind], markerfmt=".", label = 'Outlier & Not modified')             
                plt.setp([m, s], color='#F5B10C')
    for t in target:
                ax.annotate(str(t), (t, val[pos == t].values[0]), xytext=(0, 1), textcoords="offset points", fontweight='bold',fontstyle='italic', size = 10)          
          
    ax.legend()
    ax.set_xlabel("Position (nt)",fontsize = 16)
    ax.set_ylabel(feat,fontsize = 16)
    ax.set_title(title, fontsize = 16)
    for label in (ax.get_xticklabels() + ax.get_yticklabels()):
        label.set_fontsize(14)    
    ax.margins(x=0.01, y = 0.05)
    plt.savefig(path+ title+'_'+feat+'_'+ str(outlier)+'_'+ str(n_best) + '.eps') 
    plt.savefig(path+ title+'_'+feat+'_'+ str(outlier)+'_'+ str(n_best)  + '.pdf') 
    return result_lof
    
def ScatterPlot(df,feat,lab,target, outlier=0.001,neigh=[10,50], path =''):
    """
    ScatterPlot_3Cond function perform scatter plot with three features (scores) from three JACUSA2 CALL2 outputs

    :param df: a dataframe of two features (scores) and position columns
    :param lab: analysis label 
    :param feat: feature's name
    :param ref: rRNA type 18S/28S
    :param target: target position  
    :param outlier: outlier contamination for the LocalOutlierFactor method   
    :param neigh: range of neighborhood size for the LocalOutlierFactor method        
    :param path: path to save plots in eps/pdf format
    
    """     
    
    # select condition features
    X = df.loc[:,feat+'_x'].values.reshape(-1, 1)  
    Y = df.loc[:,feat+'_y'].values.reshape(-1, 1)  
    Z = df.loc[:,feat].values.reshape(-1, 1)     
    
    # perform linear regression
    linear_regressor = LinearRegression()  
    linear_regressor.fit(X, Y)  
    Y_pred = linear_regressor.predict(X)

    # identify outliers
#     print(Z >Z.median())
    filterPos = (df[feat+'_x'] >  df[feat+'_x'].median()) | (df[feat] >df[feat].median())
    idx_bool=np.full((len(X),), False, dtype=bool)
    status = np.full((len(X),), False, dtype=bool)
    max_lof_score = 0
    result_lof = []  
    val = df.loc[:,[feat+'_x',feat+'_y', feat]]
    n_best = LOF(outlier,neigh,val,filterPos)
    lof = LocalOutlierFactor(contamination =0.5, n_neighbors = n_best)
    ohat = lof.fit_predict(val)
    scores = lof.negative_outlier_factor_
    result_lof=pd.DataFrame([scores]).T.describe()
    result_lof.loc['median',:]=[np.median(scores)]            
    result_lof.loc['skew',:]=[stats.skew(scores)]
    result_lof.loc['kurtosis',:]=[stats.kurtosis(scores)]
    sortedarray = np.argsort(scores)
    nb = round(outlier*val.shape[0])
    elemts = np.where(filterPos == False)[0]
    i=0
    j=0
    while i < nb:
            if sortedarray[j] in elemts:
                    idx_bool[sortedarray[j]]=True
                    status[sortedarray[j]]=False
            else:
                    idx_bool[sortedarray[j]]=True
                    status[sortedarray[j]]=True 
                    if i==0:
                        result_lof.loc['min_modified',:]=[scores[sortedarray[j]]]                        
                    i= i+1
            j=j+1
    result_lof.loc['context_size',:]=[n_best]
    lofscores = pd.DataFrame({'pos': df['Pos'],'scores':lof.negative_outlier_factor_, 'outlier':idx_bool, 'modified':status})
    indices = np.where(status == True)    
    zz = np.zeros(len(X)).astype(int)
    zz[indices[0]]=int(1)
    lofscores = pd.DataFrame({'pos': df['Pos'],'scores':lof.negative_outlier_factor_, 'outlier':idx_bool})
    lofscores.to_csv(path+feat+'_lof_scores.csv', index = False)
    
    # scatter plot
    color_map = plt.cm.get_cmap('Spectral')
    reversed_color_map = color_map.reversed()    
    fig, ax = plt.subplots(dpi = 1200)

    zz = np.zeros(len(X) + len(target)).astype(int)
    for sz in range(len(target)):
        dff = df[df['Pos'] == target[sz]].index.values
        X= np.append(X,X[dff[0]])
        Y=np.append(Y,Y[dff[0]])
        Z=np.append(Z,Z[dff[0]])
        Y_pred=np.append(Y_pred,Y_pred[dff[0]])
        zz[-sz-1]=int(1)
    
    colors = np.array(['None', 'lime'])
    sp = ax.scatter(X, Y, c= Z, cmap=reversed_color_map, edgecolor=colors[zz])   
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    fig.colorbar(sp, label = ' Cond1vsCond3', shrink=0.5, fraction=0.046, pad = 0.04)
    ax.plot(X, Y_pred, color='red')

    
    # add labels to outliers
    dff = df[df['Pos'].isin(target)].index.values
    for ii in indices[0]:
        if ii in dff :
            ax.annotate(str(df.loc[ii,'Pos']), (X[ii], Y[ii]), xytext=(0, 1), textcoords="offset points", fontweight='bold',fontstyle='italic', size = 6)          
        else:
            if abs(ii - dff[0]) < 3:
                ax.annotate(str(df.loc[ii,'Pos']), (X[ii], Y[ii]), xytext=(0, 1), textcoords="offset points", fontweight='bold', size = 6)
            else:
                ax.annotate(str(df.loc[ii,'Pos']), (X[ii], Y[ii]), xytext=(0, 1), textcoords="offset points", size = 6)
            
    for position in  dff:
        if (position in indices[0]) == False :
            ax.annotate(str(df.loc[position,'Pos']), (X[position], Y[position]),  xytext=(0, 1), textcoords="offset points",fontstyle='italic',size = 6)     
        
    ax.set_title(lab +'_' + feat)
    ax.set_xlabel(' Cond1vsCond2')
    ax.set_ylabel(' Cond2vsCond3')
    plt.savefig(path+ lab+'_'+feat+'_'+ str(outlier) + '.eps',  bbox_inches = "tight") 
    plt.savefig(path+ lab+'_'+feat+'_'+ str(outlier) + '.pdf',  bbox_inches = "tight") 
    return result_lof

def plot (ar, x_values, y_values, xerr_values=None,xlabel='', ylabel='',title ='', color ='#3288bd', path = '', legend = False):
    """
    plot function plot horizontal barplots from a dataframe columns
    """ 
    fig, axi= plt.subplots(figsize=(3,3), dpi = 300)
    ax = ar.plot(kind = "barh",x=x_values, y = y_values, legend = False, title = title, xerr = xerr_values , xlabel = xlabel,ax = axi)#     barp(arr) 
    ax.set_xlabel(ylabel)
    for s in ['top', 'right']:
        ax.spines[s].set_visible(False)
    if legend == True:
        plt.legend()
#     plt.show()
    plt.savefig(path+title + '.pdf' ,  bbox_inches = "tight") 
    plt.savefig(path+title+ '.eps'  ,  bbox_inches = "tight")  