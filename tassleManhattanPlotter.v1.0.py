#create manhattan plots with both FDR & bonferroni corrections displayed
#output to pdfResWriter.py
import pandas as pd, numpy as np
import matplotlib.pyplot as plt
import math

#set params
DATA = 'MLM_tassle_BLB_all.txt' # a Tassel results output file
outdir = '.' #output dir
dpi = 300 #figure dpi
colors = ['darkred', 'red'] #plot cols
figsize=(15, 5) #figure dimensions (cm)
fontsize = 18
ab = 1 #bootstrap & FDR ablines (0 = no abline)

#read data
data = pd.read_csv(DATA, header = 0, sep='\t')
traits = data.Trait.unique() #column name of phenotypic data tested against genotypic

def makePlot(df, trait, fthr):
    df_grouped = df.groupby(('Chr')) #plot grouping
    bthr = math.log10(0.05 / df.shape[0]) * -1 #calc bootstrap threshold
    fig = plt.figure()
    ax = fig.add_subplot(111)
    x_labels = []
    x_labels_pos = []
    for num, (name, group) in enumerate(df_grouped):
        group.plot(kind='scatter', x='ind', y='p', figsize=figsize, color=colors[num % len(colors)], ax=ax, s=12)
        x_labels.append(name)
        x_labels_pos.append((group['ind'].iloc[-1] - (group['ind'].iloc[-1] - group['ind'].iloc[0])/2))
    plt.rcParams['xtick.labelsize'] = fontsize 
    plt.rcParams['ytick.labelsize'] = fontsize 
    ax.set_xticks(x_labels_pos)
    ax.set_xticklabels(x_labels)
    ax.set_xlim([0, len(df)])
    ax.set_ylim([0, math.ceil(pd.concat([pd.DataFrame([bthr]), df['p']], axis=0).max() + 1)]) #ylim = val above largest score or threshold + 1
    ax.set_xlabel('Chromosome', fontsize = fontsize)
    ax.set_ylabel('-log10(p)', fontsize = fontsize)
    if ab != 0 and fthr != 0:
        plt.axhline(y=bthr, color='k', linestyle='dashed') 
        plt.axhline(y=fthr, color='k', linestyle='dashdot')
    fig.savefig(f"{outdir}/{trait}.mnhtn_{dpi}dpi.png", dpi = dpi)
    plt.show()

#evaluate False Dicovery Rate
def fdr(p):
    p_values = np.sort(p)
    N = len(p_values)
    i = np.arange(1, N+1) # the 1-based i index of the p values, as in p(i)
    q = 0.05

    below = p_values < (q * i / N) # True where p(i)<qi/N
    
    if len(np.where(below)[0]) == 0:
        fdr = 0
    else:
        max_below = np.max(np.where(below)[0]) # Max Python array index where p(i)<qi/N
        fdr = math.log10(p_values[max_below]) * -1

    return fdr

nSnps = 0
for trait in traits: #make plot for each strain
    df = data.loc[data.Trait == trait] #make seperate df for each strain
    df = df.dropna(subset=['p'])
    fthr = fdr(df.p) #calc FDR thresh
    df.p = np.log10(df.p) * -1 #convert p-vals to neg logs
    df.Chr = df.Chr.astype(int) #make sure are int as will become x-axis labels
    df['ind'] = range(len(df)) #renumber indices
    makePlot(df, trait, fthr)
