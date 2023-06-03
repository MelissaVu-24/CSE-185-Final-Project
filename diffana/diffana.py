"""
Command-line script to perform differential expression analysis on RSEM files

Similar to DESeq2
"""

import argparse
import os
import sys
import math
#from . import plot as volcano
from scipy.stats import poisson

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
#	parser.add_argument("-vp", "--vp", help="Plot volcano plot." "Default: stout", metavar="file", type=str, required=False)
	
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
		
	#Poisson probability of all genes
	prob = []
	for i in range(len(meanc)):
		pr = 1-poisson.cdf(meane[i],meanc[i])
		prob.append(pr)
	logfc =[]
	
	for i in range(len(meanc)):
		if(meanc[i] == 0):
			fc = 0
		else:
			fc = math.log2(meane[i]/meanc[i])
		logfc.append(fc)
    	
	# outputting the files
	outf.write("\"\",\"log2FoldChange\",\"pvalue\"" + "\n")
	for i in range(len(prob)):
		outf.write("\"" + name[i] + "\"," + str(logfc[i])+ "," + str(prob[i]) + "\n")
    
	outf.close()
	
	# volcano plot, still in progress
#	if args.vp is None:
#		continue
#	else:
#		volcano(name, logfc, prob)
	
	sys.exit(0)	
		   
if __name__ == "__main__":
	main()
