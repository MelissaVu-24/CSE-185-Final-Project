import pandas as pd
import matplotlib.pylab as plt
import seaborn as sns
import numpy as np

def main():
  
def volcano(csv):
  """
    Parameters
    ----------
    csv : str
      Name of csv file diffana.py produces
    Return
    ----------
    Volcano Plot
  """
  df = pd.read_csv(csv,sep=",",index_col=0)
  df.head()

  plt.scatter(x=df['log2FoldChange'],y=df['pvalue'].apply(lambda x:-np.log10(x)),s=1)

  # highlight down- or up- regulated genes
  down = df[(df['log2FoldChange']<=0)]
  up = df[(df['log2FoldChange']>=0)]

  plt.scatter(x=down['log2FoldChange'],y=down['pvalue'].apply(lambda x:-np.log10(x)),s=3,label="Down-regulated",color="blue")
  plt.scatter(x=up['log2FoldChange'],y=up['pvalue'].apply(lambda x:-np.log10(x)),s=3,label="Up-regulated",color="red")
  
  for i,r in down.iterrows():
      if -np.log10(r['pvalue']) > 20:
          plt.text(x=r['log2FoldChange'],y=-np.log10(r['pvalue']),s=i,fontdict=dict(color='blue',size=6))

  for i,r in up.iterrows():
      if -np.log10(r['pvalue']) > 20:
          plt.text(x=r['log2FoldChange'],y=-np.log10(r['pvalue']),s=i,fontdict=dict(color='red',size=6))

  plt.xlabel("log2 fold change")
  plt.ylabel("-log10 pvalue")
  plt.axvline(0,color="grey",linestyle="--")
  plt.axhline(2,color="grey",linestyle="--")

  plt.legend()
  
if __name__ == "__main__":
  main()
  
  
