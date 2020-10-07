#Code developed to identify MSU genomic regions that correlate with
#significantly segregating SNP loci identified with TASSEL analysis
#across Xoo variants

import pandas as pd, numpy as np
import xlsxwriter, math, sys, getopt
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.font_manager import FontProperties

#SET PARAMS & INPUT FILES:
THRESH = ['BOOT', 'FDR'][0] #set threshold method: bootstrap vs. false discovery rate
reduce_table = 1 # choose 1 for table reduction (sign results nly)
PREFIX = 'outFile' #output files prefix
file = 'tassel.txt' #tassel output file e.g. 'MLM_tassle_BLB_all.txt'
GENES_DF = 'genes.csv' #ref genome gene location info

#read any new params
if '-h' in sys.argv:
    print("\n --Flag: Variable Name - Explanation - Default Value\n\n '-t': THRESH - choose a threshold method 'bootstrap' = 0; 'false discovery [FDR]' rate = 1 - default = '{}'\n '-r': reduce_table - select reduced output (sign results only) '1' or full output '0'' - default = '{}'\n '-f': file - tassel output file name - default = '{}'\n '-p': prefix - output file prefix - default = '{}'\n '-g': GENES_DF - list of MSU gene regions - default = '{}'\n"
          .format(THRESH, reduce_table, file, PREFIX, GENES_DF))
    sys.exit(2)
try: #input sys args
    myopts, args = getopt.getopt(sys.argv[1:],"t:r:f:p:g:")
except:
    print('I/O ERROR!')
    sys.exit(2)

#Sort non-default params and check for input errors
print('\n### Check inputted params ###')
for a, b in myopts:
    print(a, b)
    if a == '-t': THRESH = ['BOOT', 'FDR'][int(b)]
    if a == '-r': reduce_table = int(b)
    if a == '-f': file = b
    if a == '-p': PREFIX = b
    if a == '-g': GENES_DF = b

swch = 0
for sysarg in sys.argv[1:]:
    if sysarg[0] == '-': swch = 1
    if sysarg[0] != '-' and swch == 1:
        swch = 0
        continue
    if sysarg[0] != '-' and swch == 0:
        print('No flag added for: ', sysarg, sep = '')
        swch = 0
        sys.exit(2)

if reduce_table == 1:
    print('Table reduction selected!\nOutputting columns for traits with significant QTLs only.\n')
    txt = ''
else: txt = '_LongFormat'
if PREFIX[-1] != '_' or PREFIX[-1] != '-' or PREFIX[-1] != '.': PREFIX += '_'
print(f"You have selected '{THRESH}' significance thresholding\n")
print(f"\n{PREFIX}{file}\n")

#Load data
tmp_data = pd.read_csv(file, header = 0, sep='\t')
tbr = ['Trait', 'Marker', 'Chr', 'Pos', 'p']
data = pd.DataFrame()

for t in tbr: data = pd.concat([data, tmp_data[t]], axis = 1)
traits = data.Trait.unique()

genes_df = pd.read_csv(GENES_DF, header=0)

def main():
    #find the best marker at ea chrom for ea trait (blb strain) and test its significance
    print(f"Traits: {traits}")
    mat, locs_mat, all_sign_locs, thrs = [], [], [], []
    for trait in traits:
        print(f"Trait no.{np.where(traits==trait)[0][0] + 1}: {trait}")
        df = data.loc[data.Trait == trait] #subset of file
        df = df.dropna(subset=['p']) #remve NANs
        df = df.reset_index(drop=True) #reset index
        newcol = df.iloc[:, 4:5].applymap(lambda x: math.log10(x) * -1) #new col is -log10
        newcol.columns = ['logP']
        df = pd.concat([df, newcol], axis=1)
        if THRESH == 'FDR':
            thr, txt2 = fdr(df.p), '_fdr' #false discovery rate threshold
            thrs.append(thr)
        if THRESH == 'BOOT': thr, txt2 = np.log10(0.05 / len(df)) * -1, '_boot' #bootstrap threshold
        genes, locs = [], []
        for gene in genes_df.locus.unique(): #go thru all loci ca. 55000 - this is [unavoidably!] slow
            gene_df = genes_df.loc[genes_df.locus == gene] #make DF of just the gene
            end_locs = gene_df.iloc[0,3:5].sort_values() #id 5'-3' posns
            chrom = float(gene_df.iloc[0, 0][3:]) #convert chsome format btwn DFs
            indxs = df.index[(df.Chr == chrom) & (df.logP > thr) & (df.Pos > end_locs[0]) & (df.Pos < end_locs[1])] #id sign SNPs with -log10 > threshold
            if len(indxs) > 0: genes.append(gene)
            else: genes.append('-') #if no sign SNPs add a blank
            locs.append(list(df.Marker[indxs])), all_sign_locs.extend(list(df.Marker[indxs])) #add to list of loci
        mat.append(genes), locs_mat.append(locs) #record genes and loci

    tab = pd.DataFrame(mat).T #tab for final output
    tab.columns = traits
    locs_tab = pd.DataFrame(locs_mat).T #tab for setting colours
    locs_tab.columns = traits

    if reduce_table == 1: tab, locs_tab = reduceTab(tab, locs_tab)

    to_drop = []
    for i in tab.index:
        if '-' in tab.iloc[i, :].unique() and len(tab.iloc[i, :].unique()) == 1: to_drop.append(i)

    tab = tab.drop(to_drop)
    locs_tab = locs_tab.drop(to_drop)
    gene_funcs = tab.replace(genes_df.set_index('locus')['annotation'])

    #change col 'rs#' to 'rs'
    #hapmap = hapmap.rename(columns = {'rs#':'rs'})
    colours = genCols(tab)
    writeTable(tab, colours, txt, txt2)
    writeXL(tab, colours, gene_funcs, txt, txt2)
    chroms = data.Chr.unique()[~np.isnan(data.Chr.unique())]
    if THRESH == 'BOOT': thrs = [thr] * chroms.shape[0]
    restLoci(data, all_sign_locs, txt, txt2, thrs, chroms)

#make colours DF to match final DF
def genCols(tab):
    #create list of lists of mapped colours
    cmap = ['red', 'yellow', 'b', 'green', 'cyan', 'm', 'crimson', 'pink', 'limegreen', 'royalblue', 'deeppink', 
        'firebrick', 'coral', 'olive', 'mediumpurple', 'khaki', 'greenyellow', 'gray', 'forestgreen', 'gold', 
        'darkturquoise', 'dodgerblue', 'darkviolet', 'orangered', 'paleturquoise'] #this needs to be long enough; can also use: cmap = plt.get_cmap('gnuplot') // colors = [cmap(i) for i in np.linspace(0, 1, number)] to automate properly but need to work it out

    colours = pd.DataFrame().reindex_like(tab)
    colours.iloc[:,:] = 'silver'

    for chrom in data.Chr.unique()[~np.isnan(data.Chr.unique())]:
        for col in tab.columns:
            colours[col][tab[col].str.contains(f"Os{str(int(chrom)).zfill(2)}")] = cmap[int(chrom) - 1]
    colours = colours.values.tolist()

    return colours

#make colour coded table & write as PDF
def writeTable(tab, colours, txt, txt2):
    fig, ax = plt.subplots(figsize=(tab.shape[1], 4)) #try (8, 4) & alter font size!!!
    ax.axis('tight')
    ax.axis('off')
    the_table = ax.table(cellText=tab.values, colColours=['lime'] * tab.shape[1],
                     cellColours=colours, colLabels=tab.columns, loc='center',
                     cellLoc='center')
    the_table.auto_set_font_size(False)
    [cell.set_text_props(fontproperties=FontProperties(size='5', weight='bold')) for key, cell in the_table.get_celld().items()] #make bold font 
    pp = PdfPages(f"{PREFIX}MLM_BLB_tassel_LDgenes{txt}{txt2}_msu7.pdf")
    pp.savefig(fig, bbox_inches='tight')
    pp.close()
    plt.close(fig)

    pd.DataFrame(colours).to_csv(f"{PREFIX}MLM_BLB_tassel_LDgenes_colours{txt}{txt2}_msu7.csv", index = False)

#calculate false discovery rate
def fdr(p):
    p_values = np.sort(p)
    N = len(p_values)
    i = np.arange(1, N+1) # the 1-based i index of the p values, as in p(i)
    q = 0.05

    below = p_values < (q * i / N) # True where p(i)<qi/N
    
    if len(np.where(below)[0]) == 0: fdr = 100000
    else:
        max_below = np.max(np.where(below)[0]) # Max Python array index where p(i)<qi/N
        fdr = math.log10(p_values[max_below]) * -1

    return fdr

def reduceTab(tab, locs_tab):
    for col in tab.columns:
        if '-' in list(tab[col].unique()) and len(list(tab[col].unique())) == 1:
            del tab[col]
            del locs_tab[col]

    return tab, locs_tab

#write colour coded XL sheet matching PDF
def writeXL(tab, colours, gene_funcs, txt, txt2):
    colours = pd.DataFrame(colours)
    writer = pd.ExcelWriter(f"{PREFIX}MLM_BLB_tassel_LDgenes{txt}{txt2}_msu7.xlsx", engine='xlsxwriter')
    writer2 = pd.ExcelWriter(f"{PREFIX}MLM_BLB_tassel_LDfuncGenes{txt}{txt2}_msu7.xlsx", engine='xlsxwriter')

    tab.to_excel(writer, sheet_name='Sheet1', index=False)
    gene_funcs.to_excel(writer2, sheet_name='Sheet1', index=False)

    workbook = writer.book
    workbook2 = writer2.book
    worksheet = writer.sheets['Sheet1']
    worksheet2 = writer2.sheets['Sheet1']

    dik = {}
    for i in range(tab.shape[1]):
        for j in range(tab.shape[0]):
            if tab.iloc[j, i] not in dik: dik[f'"{tab.iloc[j, i]}"'] = colours.iloc[j, i]

    dik2 = {}
    for i in range(gene_funcs.shape[1]):
        for j in range(gene_funcs.shape[0]):
            if gene_funcs.iloc[j, i] not in dik2: dik2[f'"{gene_funcs.iloc[j, i]}"'] = colours.iloc[j, i]

    txt3 = f"A2:{chr(tab.shape[1] + 64)}{tab.shape[0] + 1}"
    for val, color in dik.items():
        fmt = workbook.add_format({'bg_color': color})
        worksheet.conditional_format(txt3, {'type': 'cell',
                                           'criteria': '=',
                                           'value': val,
                                           'format': fmt})
        for val, color in dik2.items():
            fmt2 = workbook2.add_format({'bg_color': color})
            worksheet2.conditional_format(txt3, {'type': 'cell',
                                           'criteria': '=',
                                           'value': val,
                                           'format': fmt2})
    writer.save()
    writer2.save()

#work out the rest of the sign. SNPs not found in MSU regions
def restLoci(data, all_sign_locs, txt, txt2, thrs, chroms):
    other_sgnLocs = pd.DataFrame()
    rest = data.loc[~data.Marker.isin(all_sign_locs)]
    new = rest.iloc[:, 4:5].applymap(lambda x: math.log10(x) * -1)
    new.columns = ['logP']
    rest = pd.concat([rest, new], axis=1)
    for trait in traits:
        other_inds = rest.index[(rest.logP > thrs[traits.tolist().index(trait)]) & (rest.Trait == trait)]
        other_sgnLocs = pd.concat([other_sgnLocs, rest.loc[other_inds, :].sort_values('Marker')])

    other_sgnLocs.to_csv(f"{PREFIX}MLM_BLB_tassel_LDgenes_other_sign_locs{txt}{txt2}_msu7.csv", index = False)

    all_sign_locs = pd.DataFrame(all_sign_locs)
    all_sign_locs.columns = ['markers']
    all_sign_locs.to_csv(f"{PREFIX}{THRESH}_identifed_loci_msu7.csv", index=False)

if __name__ == '__main__': main()
