"""
Command-line script to perform differential expression analysis on RSEM files

Similar to DESeq2
"""

import argparse
import os
import sys
import math
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
  
    args = parser.parse_args()
    
    #output file
    if args.out is None:
      outf = sys.stdout
    else: outf = open(args.out, "w")
    

    #input control file1
    inputc = open(args.RSEMcon[0], "r")
    #input control file2
    inpute = open(args.RSEMexp[0], "r")
    
    countsc = [len(inputc.readlines())-1][len(args.RSEMcon)]
    countse = [len(inpute.readlines())-1][len(args.RSEMexp)]
    name = []
    #filling in the control counts for genes
    track = 0
    linenum = 0
    for line in inputc:
      if track == 0:
        track = 1
        continue
      countsc[linenum].append(line.split('	')[4])
      name.append(line.split('	')[0])
      linenum+=1
      
    inputc.close()
    for i in range(1, len(args.RSEMcon)):
      track = 0
      linenum = 0

      for line in inputc:
        if track == 0:
          track = 1
          continue
        countsc[linenum].append(line.split('	')[4])
        linenum+=1
      inputc.close()
        
    #filling in the experiment counts for genes
    track = 0
    linenum = 0
    for line in inpute:
      if track == 0:
        track = 1
        continue
      countse[linenum].append(float(line.split('	')[4]))
      linenum+=1
      
    inpute.close()
    for i in range(1, len(args.RSEMexp)):
      track = 0
      linenum = 0

      for line in inpute:
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
      meanc.append(mean)

    #Poisson probability of all genes
    prob = []
    for i in range(len(meanc)):
      pr = (meanc[i]**meane[i]) * (math.exp(-meanc[i])) / math.factorial(int(meane[i]))
      
    outf.write("\"\",\"baseMean\",\"log2FoldChange\",\"pvalue\"" + "\n")
    for i in range(len(meanc)):
      outf.write("\"" + name[i] + "\"," + math.log2(meane[i]/meanc[i]) + "," + prob[i] + "," + "\n")
