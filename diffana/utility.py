import pandas as pd
import matplotlib.pylab as plt
import seaborn as sns
import numpy as np
import math

def volcano(names, fold_change, pval, file):
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
	
	# highlight down- or up- regulated genes
	size = len(fold_change)
	down_names = []
	down_fc = []
	down_pv = []
	up_fc = []
	up_pv = []
	up_names = []
	
	for i in range(size):
		if pval[i] != 'NA' and pval[i] != 0 and 1 - pval[i] != 0:
			if fold_change[i] <= 0:
				down_fc.append(fold_change[i])
				down_pv.append(-math.log(1-pval[i],10)) 
				down_names.append(names[i])
			if fold_change[i] >= 0:
				up_fc.append(fold_change[i])
				up_pv.append(-math.log(pval[i],10))
				up_names.append(names[i])
		else:
			continue;
	plt.scatter(x=down_fc,y=down_pv,s=3,label="Down-regulated",color="blue")
	plt.scatter(x=up_fc,y=up_pv,s=3,label="Up-regulated",color="red")
  
	#Label names of differentially expressed genes up- or down regulated
	for i in range(len(down_fc)):
		if down_pv[i] > 20:
			plt.text(x=down_fc[i],y=down_pv[i],s=down_names[i],fontdict=dict(color='blue',size=6))
      
	for i in range(len(up_fc)):
		if up_pv[i] > 20:
			plt.text(x=up_fc[i],y=up_pv[i],s=up_names[i],fontdict=dict(color='red',size=6))

	plt.xlabel("log2 fold change")
	plt.ylabel("-log10 pvalue")
	plt.axvline(0,color="grey",linestyle="--")
	plt.axhline(2,color="grey",linestyle="--")

	plt.legend()
	plt.savefig(file)

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
    n:int
    	the converted binomial distribution parameter(total number of trials)
	'''
	print(variance, mean)
	if variance==0 and mean ==0:
		return [0,0,0]
	p = mean/variance
	r = mean * p /(1-p)
	n = r/p
	print([p,int(r),int(n)])
	return [p,int(r),int(n)]
        
    
