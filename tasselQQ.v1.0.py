import pandas as pd, numpy as np
import matplotlib.pyplot as plt

#set params
DATA = 'MLM_tassle_BLB_all.txt' # a Tassel results output file
outdir = '.' #output dir
dpi = 300 #figure dpi
figsize = (6, 5) #figure dimensions (cm)
fontsize = 18
linetype = "--"

data = pd.read_csv(DATA, header = 0, sep='\t')
traits = data.Trait.unique()

def ppoints(n, a): #x-axis in QQ plot
    """ numpy analogue of `R`'s `ppoints` function
        see https://stackoverflow.com/questions/20292216/imitating-ppoints-r-function-in-python
    """
    try:
        n = np.float(len(n))
    except TypeError:
        n = np.float(n)
    return np.log10((np.arange(n) + 1 - a)/(n + 1 - 2*a)) * -1

def makePlot(unifp, log_vals, trait):
    a, b = int(unifp[0][0]) + 1, int(log_vals[-1]) + 1
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111)
    plt.rcParams['xtick.labelsize'] = fontsize
    plt.rcParams['ytick.labelsize'] = fontsize
    plt.scatter(unifp, log_vals[::-1], c="0", s=20)
    ax.set_xlim([0, a])
    ax.set_ylim([0, b]) #ylim = val above largest score or threshold + 1
    ax.plot([0, a], [0, b], ls=linetype, c=".3")
    ax.set_xlabel('Expected -log10(p)', fontsize = fontsize)
    ax.set_ylabel('Observed -log10(p)', fontsize = fontsize)
    fig.savefig(f"{outdir}/{trait}.QQ_{dpi}dpi.png", dpi = dpi)
    plt.show()

for trait in traits:
    df = data.loc[data.Trait == trait]
    df = df.dropna(subset=['p'])
    df.p = np.log10(df.p) * -1 #convert p-vals to neg logs
    log_vals = np.sort(df.p)
    n = len(log_vals)
    unifp = pd.DataFrame(ppoints(n, 0.5)) #2nd value can change according to n (see R ppoints function)
    makePlot(unifp, log_vals, trait)
