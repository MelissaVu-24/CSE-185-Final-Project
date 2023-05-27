Description:
This project tool is aimed to conduct a differential expression analysis on an RSEM dataset that is in the format BLANK. It is supposed to be similar 
to the Deseq2 tool. The tool will output a file in the format BLANK which includes a log2FoldChange and a P-value column. THe tool also outputs a volcano plot visualizing the differentially expressed genes that are upstream and downstream.

Installation Instructions:
The tool requires STAR and RSEM. You can install these tools with the command:
  pip install STAR RSEM
After installing the necessary libraries, install the tool using:
  python setup.py install --prefix=$HOME
To test that the tool runs use:
  diffana --help
  
Basic Usage Instructions:
1. 

