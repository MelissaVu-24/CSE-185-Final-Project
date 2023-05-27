# Diffana CSE 185 Final Project (in progress)
This project tool is aimed to conduct a differential expression analysis on an RSEM dataset that is in the format *.genes.results of genes for two different conditions. It is supposed to be similar to the Deseq2 tool, where it will compare the datasets of the different conditions to find which genes are differentially expressed. The tool will output a file in the format of .csv which includes a log2FoldChange and a P-value column. The tool also outputs a volcano plot visualizing the differentially expressed genes that are upstream and downstream.

## Installation Instructions:
The tool requires STAR and RSEM. You can install these tools with the command:

`pip install STAR RSEM`

After installing the necessary libraries, install the tool using:

  `python setup.py install`
  
To test that the tool runs use:

  `diffana --help`
  
## Basic Usage Instructions:
The basic usage is:

  `python diffana --RSEMcon in.genes.results`
  
To run diffana on test examples from files in this repo:

  `python diffana --RSEMcon ~/1/test1.txt --RSEMcon ~/1/test3.txt --RSEMexp ~/1/test2.txt --RSEMexp ~/1/test4.txt` 
  
## Full Usage Instructions:

