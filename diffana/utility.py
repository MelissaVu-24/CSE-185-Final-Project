import pandas as pd
import matplotlib.pylab as plt
import seaborn as sns
import numpy as np
import math
import statistics

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
	if variance==0 and mean ==0:
		return [0,0,0]
	p = mean/variance
	r = (mean **2)/(variance-mean)
	n = r/p
	return [p,int(r),int(n)]
        
def sizeFactor(counts):
	'''
	Parameters
	------------
	counts : 2D-array
		the 2D-array containing arrays of counts for a gene in a condition
		
	Returns
	-----------
	sizefactor : array	
		an array containing the size factor of each sample
	'''
	sizefactor = []
	total = []
	for gene in counts:
		product = 1
		for sample in gene:
			product = product * sample
		total.append(product ** (1/len(gene)))
	for j in range(len(counts[0])):	
		rooc = []
		for i in range(len(counts)):
			if(total[i] == 0):
				rooc.append(1)
			else:
				rooc.append(counts[i][j]/total[i])
		sizefactor.append(statistics.median(rooc))
	return sizefactor
			
def meanCond(sizefactor, counts):
	'''
	Parameters
	------------
	counts : 2D-array
		the 2D-array containing arrays of counts for a gene in a condition
		
	sizefactor : array	
		an array containing the size factor of each sample
		
	Returns
	-----------
	mean : array	
		an array containing the adjusted means of the genes in one condition
	'''
	mean = []
	for gene in counts:
		total = 0
		for i in range(len(gene)):
			total = total + gene[i]/sizefactor[i]
		mean.append(total/len(gene))
	return mean
			
def varCond(sizefactor, counts, mean):
	'''
	Parameters
	------------
	counts : 2D-array
		the 2D-array containing arrays of counts for a gene in a condition
		
	sizefactor : array	
		an array containing the size factor of each sample
		
	mean : array	
		an array containing the adjusted means of the genes in one condition
		
	Returns
	-----------
	var : array
		array containing the raw variance of the genes in one condition
	'''
	w = []
	z = []
	for i in range(len(counts)):
		wtotal = 0
		ztotal = 0
		for j in range(len(counts[i])):
			wtotal = wtotal + (counts[i][j]/sizefactor[j] - mean[i]) ** 2
			ztotal = ztotal + (1/sizefactor[j])
		w.append(wtotal/(len(counts[i])-1))
		z.append(mean[i] * ztotal/(len(counts[i])))
	localvar = []
	for i in range(len(w)):
		localvar.append(w[i]-z[i])
	var = []
	for i in range(len(w)):
		var.append(w[i]-z[i])
	return var		

def sampleVar(sizefactor, mean, variance):
	'''
	Parameters
	------------
	counts : 2D-array
		the 2D-array containing arrays of counts for a gene in a condition
		
	sizefactor : array	
		an array containing the size factor of each sample
		
	variance : array
		an array containing the raw variance of the genes in one condition
		
	Returns
	-----------
	[sampleVar, sampleMean] : 2D-array
		an array containg sampleVar, the sum of the variance of the samples for one condition and sampleMean, the sum of the mean for one condition.
	'''
	sampleVar = []
	sampleMean = []
	for i in range(len(mean)):
		sampleMean.append(mean[i] * len(sizefactor))
		var = mean[i] * len(sizefactor)
		for j in range(len(sizefactor)):
			var = var + (sizefactor[j]**2) * variance[i]
		sampleVar.append(var)
	return [sampleVar, sampleMean]
	
			
