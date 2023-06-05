import pandas as pd
import matplotlib.pylab as plt
import seaborn as sns
import numpy as np
'''  
def volcano(names, fold_change, pval):
    """
    Parameters
    ----------
    names : str[]
    fold_change : int[]
    pval : int[]
    list of p-values
    Return
    ----------
    Volcano Plot
    """
    plt.scatter(x=fold_change,y=pval.apply(lambda x:-np.log10(x)),s=1)
    # highlight down- or up- regulated genes
    size = fold_change.len()
    down_names[size]
    down_fc[size]
    down_pv[size]
    up_fc[size]
    up_pv[size]
    up_names[size]
    for i in fold_change.len():
        if fold_change[i] <= 0:
            down_fc.append(fold_change[i])
            down_pv.append(pval[i]) 
            down_names.append(names[i])
        if fold_change[i] >= 0:
            up_fc.append(fold_change[i])
            up_pv.append(pval[i])
            up_names.append(names[i])

            plt.scatter(x=down_fc,y=down_pv.apply(lambda x:-np.log10(x)),s=3,label="Down-regulated",color="blue")
            plt.scatter(x=up_fc,y=up_pv.apply(lambda x:-np.log10(x)),s=3,label="Up-regulated",color="red")
  
    for i in down_fc.len():
        if -np.log10(down_pv[i]) > 20:
            plt.text(x=down_fc[i],y=-np.log10(down_pv[i]),s=down_names[i],fontdict=dict(color='blue',size=6))
      
    for i in up_fc.len():
        if -np.log10(up_pv[i]) > 20:
            plt.text(x=up_fc[i],y=-np.log10(up_pv[i]),s=up_names[i],fontdict=dict(color='blue',size=6))

    plt.xlabel("log2 fold change")
    plt.ylabel("-log10 pvalue")
    plt.axvline(0,color="grey",linestyle="--")
    plt.axhline(2,color="grey",linestyle="--")

    plt.legend()
      
'''   
def convertParameters(mean, variance):
    '''
    Parameters
    ----------
    mean : int
        the mean of the control dataset
        
    dispersion : int
        the dispersion of the control dataset
        
    Returns
    ----------
    r: int
        the converted binomial distribution parameter(number of successes)
        
    p: int
        the converted binomial distribution parameter(probability of success)
    '''

    
    p = (variance - mean)/variance
    r = (mean * mean)/(variance - mean)
    return [p,r]
        
    
