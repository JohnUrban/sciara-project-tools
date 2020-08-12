import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns

def plotEdgeR_XvA_old(f):
    classes = f['chr']
    X = f['logCPM']
    y = f['logFC']

    df = pd.DataFrame(map(list, zip(*[X.T, y.ravel().T])))
    df = df.reset_index()
    df['index'] = classes[:]

    colors = {'X':'red', 'A':'blue'}
    lw = {'X':2, 'A':4}
    g = sns.jointplot(X, y, kind='reg', scatter = False, fit_reg=False, marginal_kws=dict(hist=False, kde=False))
    plt.xlim(-2,14)
    plt.ylim(-15.5,15.5)
    for i, subdata in df.groupby("index"):
        sns.kdeplot(subdata.iloc[:,1], ax=g.ax_marg_x, legend=False, color=colors[i], label=i, lw=lw[i])
        sns.kdeplot(subdata.iloc[:,2], ax=g.ax_marg_y, vertical=True, legend=False, color=colors[i], label=i, lw=lw[i])
        g.ax_joint.plot(subdata.iloc[:,1], subdata.iloc[:,2], "o", ms=2, color=colors[i], label=i)
    plt.tight_layout()


def plotEdgeR_XvA(f):
    df = f.copy()

    colors = {'X':'red', 'A':'blue'}
    lw = {'X':2, 'A':4}

    g = sns.jointplot(df['logCPM'], df['logFC'], kind='reg', scatter = False, fit_reg=False, marginal_kws=dict(hist=False, kde=False))
    plt.xlim(-2,14)
    plt.ylim(-15.5,15.5)

    for i, subdata in df.groupby("chr"):
        #print(subdata)
        sns.kdeplot(subdata['logCPM'], ax=g.ax_marg_x, legend=False, color=colors[i], label=i, lw=lw[i])
        sns.kdeplot(subdata['logFC'], ax=g.ax_marg_y, vertical=True, legend=False, color=colors[i], label=i, lw=lw[i])
        g.ax_joint.plot(subdata['logCPM'], subdata['logFC'], "o", ms=2, color=colors[i], label=i)
    plt.tight_layout()


def plotDEseq_XvA(f, logFCtype='default'):
    ''' logFCtype in (default, shrink, ashr, apeglm)'''
    df = f.copy()
    lfcdict = {'default':'log2FoldChange', 'shrink':'log2FoldChange.shrink', 'ashr':'log2FoldChange.ashr', 'apeglm':'log2FoldChange.apeglm'}
    logFC = lfcdict[logFCtype]
    
    colors = {'X':'red', 'A':'blue'}
    lw = {'X':2, 'A':4}

    g = sns.jointplot(np.log2(df['baseMean']), df[logFC], kind='reg', scatter = False, fit_reg=False, marginal_kws=dict(hist=False, kde=False))
    plt.xlim(-1,20)
    plt.ylim(-10.5,10.5)

    for i, subdata in df.groupby("chr"):
        #print(subdata)
        sns.kdeplot(np.log2(subdata['baseMean']), ax=g.ax_marg_x, legend=False, color=colors[i], label=i, lw=lw[i])
        sns.kdeplot(subdata[logFC], ax=g.ax_marg_y, vertical=True, legend=False, color=colors[i], label=i, lw=lw[i])
        g.ax_joint.plot(np.log2(subdata['baseMean']), subdata[logFC], "o", ms=2, color=colors[i], label=i)
    plt.tight_layout()



def violin(d,stages = None, rmleg=False, save=False, ax=False, logFC=None):
    if stages is None:
        stages = d.keys()
    ALL = pd.concat([d[e] for e in stages], sort=False)

    sns.set_style('darkgrid')
    if ax:
        g = sns.violinplot(x='stage', y=logFC, hue='chr', data=ALL, 
                       notch=True, orient='v', split=True, inner='box', palette={'X':'red','A':'blue'},
                          ax=ax)
    else:
        g = sns.violinplot(x='stage', y=logFC, hue='chr', data=ALL, 
                       notch=True, orient='v', split=True, inner='box', palette={'X':'red','A':'blue'})
    if rmleg:
        g.legend_.remove()
    plt.tight_layout()
    if save:
        plt.savefig(save) #'pyFigures/violin-dosage-compensation.pdf'

def violin_edgeR(d,stages = None, rmleg=False, save=False, ax=False):
    violin(d,stages, rmleg, save, ax, logFC='logFC')

def violin_DEseq(d,stages = None, rmleg=False, save=False, ax=False):
    violin(d,stages, rmleg, save, ax, logFC='log2FoldChange')


def volcano(df, ax, FC=0.1,FDR=2, stage="", method="",xcol='logFC', ycol='log10fdr', y2log=False, x2log=False):
    x = df[xcol]
    y = df[ycol]
    if y2log:
        y = -1*np.log10(y)
    if x2log:
        x = -1*np.log2(x)
    gate = (x.abs()>=FC) & (y>=FDR)
    ax.scatter(x=x, y=y, c='blue')
    ax.scatter(x=x[gate], y=y[gate], c='red')
    ax.set_title(method + ':\n' + stage + '\nAbs Log2FC >= '+str(FC) + '\nFDR <= ' + str(100*(10**(-FDR))) + '%')

def pvalDistro(df, ax, stage='', method='', pval='PValue'):
    sns.distplot(df.dropna()[pval], ax=ax)
    ax.set_title(method + ':\n' + stage)


def setanal(A,B,C):
    print("A", len(A))
    print("B", len(B))
    print("C", len(C))
    print("Union AUB", len(A.union(B)))
    print("Union AUC", len(A.union(C)))
    print("Union BUC", len(B.union(C)))
    print("Union AUBUC", len(A.union(B).union(C)))
    
    print("intersection A^B", len(A.intersection(B)))
    print("intersection A^C", len(A.intersection(C)))
    print("intersection B^C", len(B.intersection(C)))
    print("intersection A^B^C", len(A.intersection(B).intersection(C)))
    
    print("Jaccard A,B", len(A.intersection(B)) / len(A.union(B)) )
    print("Jaccard A,C", len(A.intersection(C)) / len(A.union(C)) )
    print("Jaccard B,C", len(B.intersection(C)) / len(B.union(C)) )
    print("Jaccard A,B,C", len(A.intersection(B).intersection(C)) / len(A.union(B).union(C)) )
    
    print("Percent of A in B", len(A.intersection(B)) / len(A))
    print("Percent of A in C", len(A.intersection(C)) / len(A))
    
    print("Percent of B in A", len(A.intersection(B)) / len(B))
    print("Percent of B in C", len(B.intersection(C)) / len(B))
    
    print("Percent of C in A", len(A.intersection(C)) / len(C))
    print("Percent of C in B", len(B.intersection(C)) / len(C))
    
    print("Num in 2+", len(A.intersection(B).union(A.intersection(C)).union(B.intersection(C))) )
    print("Pct in 2+", len(A.intersection(B).union(A.intersection(C)).union(B.intersection(C))) / len(A.union(B).union(C)) ) 
