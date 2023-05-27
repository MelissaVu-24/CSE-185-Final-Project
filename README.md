# Diffana CSE 185 Final Project (in progress)
This project tool is aimed to conduct a differential expression analysis on an RSEM dataset that is in the format *.genes.results of genes for two different conditions. It is supposed to be similar to the Deseq2 tool, where it will compare the datasets of the different conditions to find which genes are differentially expressed. The tool will output a file in the format of .csv which includes a log2FoldChange and a p-value column. The p-value is predicted using Poisson distribution.

## Installation Instructions:

You may install the tool using this way:
```
git clone https://github.com/MelissaVu-24/CSE-185-Final-Project
cd CSE-185-Final-Project
python setup.py install`
```  
If you do not have root acces, instead use:
  `python setup.py install --user`
 
access the tool by:
  `cd .local/bin`
  
To test that the tool runs use:

  `diffana --help`
  
## Basic Usage Instructions:
The basic usage is:

  `python diffana [--RSEMcon control.genes.results] [--RSEMexp experiment.genes.results] [other options]`
  
  
To run diffana on test examples from files in this repo:

  `python diffana --RSEMcon ~/1/test1.txt --RSEMcon ~/1/test3.txt --RSEMexp ~/1/test2.txt --RSEMexp ~/1/test4.txt` 
  
It should produce output like this:
```
"","baseMean","log2FoldChange","pvalue"
"ENSMUSG00000000001",-1.218388824925326,1.0
"ENSMUSG00000000003",0,0.0
```
  
## diffana options:

    *`-o File`, `-out File`:Write the output to a file in the csv form. Default will print to stdout.
    
# Contributors

This repository was generated by Owin Gong and Melissa Vu, with inspiration from the [Deseq2]([https://github.com/gymreklab/trtools].


