"""
Command-line script to perform differential expression analysis on RSEM files

Similar to DESeq2
"""

import argparse
import os
import sys
import math
import numpy
from . import utility
from scipy.stats import poisson
from scipy.stats import nbinom

def main():
	parser = argparse.ArgumentParser(
		prog="diffana",
	description="Command-line script to perform differential expression analysis on RSEM files"
	)
      
	# Input controls
	parser.add_argument("--RSEMcon", action = 'append', help="control RSEM files", type=str)

	# Input experiments 
	parser.add_argument("--RSEMexp", action = 'append', help="experiment RSEM files", type=str)

	# Output
	parser.add_argument("-o", "--out", help="Write output to file." "Default: stdout", metavar="file", type=str, required=False)
	parser.add_argument("-vp", "--vp", help="Plot volcano plot." "Default: stout", metavar="file", type=str, required=False)
	
	args = parser.parse_args()

	#output file
	if args.out is None:
  		outf = sys.stdout
	else: 
		outf = open(args.out, "w")
		
	#input control file1
	inputc = open(args.RSEMcon[0], "r")
	#input control file2
	inpute = open(args.RSEMexp[0], "r")

	linec = inputc.readlines()
	linee = inpute.readlines()

	countsc = []
	countse = []
	for i in range(len(linec)-1):
		countsc.append([])
		countse.append([])
	name = []
	#filling in the control counts for genes
	track = 0
	linenum = 0
	for line in linec:
		if track == 0:
			track = 1
			continue
		countsc[linenum].append(float(line.split('	')[4]))
		name.append(line.split('	')[0])
		linenum+=1
      
	inputc.close()
	for i in range(1, len(args.RSEMcon)):
		track = 0
		linenum = 0
		inputc = open(args.RSEMcon[i], "r")
		linec = inputc.readlines()
		for line in linec:
			if track == 0:
				track = 1
				continue
			countsc[linenum].append(float(line.split('	')[4]))
			linenum+=1
			inputc.close()
        
	#filling in the experiment counts for genes
	track = 0
	linenum = 0
    
	for line in linee:
		if track == 0:
			track = 1
			continue
		countse[linenum].append(float(line.split('	')[4]))
		linenum+=1
      
	inpute.close()
	for i in range(1, len(args.RSEMexp)):
		track = 0
		linenum = 0
		inpute = open(args.RSEMexp[i], "r")
		linee = inpute.readlines()
		for line in linee:
			if track == 0:
				track = 1
			continue
		countse[linenum].append(float(line.split('	')[4]))
		linenum+=1
		inpute.close()
    
	# The mean counts for the two data sets
	meanc = []
	meane = []
    
	for i in countsc:
		total = 0
		for j in i:
			total = total + j
		mean = total/len(i)
		meanc.append(mean)

	for i in countse:
		total = 0
		for j in i:
			total = total + j
		mean = total/len(i)
		meane.append(mean)
		
	# The variance of the data sets	
	variancec = []
	variancee = []
	for i in countsc:
		variancee.append(numpy.var(countsc[i]))

	for i in countse:
		variancee.append(numpy.var(countsc[i]))
		
	#Filling in the logfc
	logfc =[]
	
	for i in range(len(meanc)):
		if(meanc[i] == 0):
			fc = 0
		else:
			fc = math.log2(meane[i]/meanc[i])
		logfc.append(fc)
	prob = []
	#Negative Binomial Distribution probability of all genes
	parametersc = []
	parameterse = []
	for i in range(len(meanc)):
		parametersc.append(utility.convertParameters(meanc[i], variancec[i]))
	for i in range(len(meane)):
		parameterse.append(utility.convertParameters(meane[i], variancee[i]))
	for i in range(len(meanc)):
		pr = nbinom.pmf(1-parameterse[i][1],parameterse[i][1],parameterse[i][0])
		prob.append(pr)
'''				    
	#Poisson probability of all genes
	for i in range(len(meanc)):
		pr = poisson.pmf(meane[i],meanc[i])
		prob.append(pr)
'''   	
	
	# outputting the files
	outf.write("\"\",\"log2FoldChange\",\"pvalue\"" + "\n")
	for i in range(len(prob)):
		outf.write("\"" + name[i] + "\"," + str(logfc[i])+ "," + str(prob[i]) + "\n")
    
	outf.close()
	
	# volcano plot, still in progress
	if args.vp is not None:
		utility.volcano(name, logfc, prob)
		
	
	sys.exit(0)	
		   
if __name__ == "__main__":
	main()
