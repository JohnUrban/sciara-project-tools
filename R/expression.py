import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.cluster.hierarchy import *


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


def identity(x):
    return x

def generic_joint_plot(f, x, y, group, ax, xlim=None, ylim=None, lw=(1,1), colors=None, defaultcolors=["red", "black", "blue", "green"], fxnX=identity, fxnY=identity,):
    ''' f = pandas dataframe
        x = column name
        y = column name
        group = column name
        ax = ax from something like "fig, ax = plt.subplots()"
        xlim = tuple
        ylim = tuple
        lw = tuple
        colors = dict with unique values in group as keys and colors as values.
            if left as None, tolerates up to 4 values -- should only be 2 (e.g. True/False, 0/1)
    '''
    
    df = f.copy()

    if colors is None:
        print("CHOOSING COLORS\n")
        uni = df[group].unique()
        colors = {}
        for i in range(len(uni)):
            colors[uni[i]] = defaultcolors[i]
            

    g = sns.jointplot(fxnX(df[x]), fxnY(df[y]), kind='reg', scatter = False, fit_reg=False, marginal_kws=dict(hist=False, kde=False))

    if xlim is not None:
        plt.xlim(xlim[0], xlim[1])
    if ylim is not None:
        plt.ylim(ylim[0], ylim[1])

    for i, subdata in df.groupby(group):
        print(i)
        sns.kdeplot(fxnX(subdata[x]), ax=g.ax_marg_x, legend=False, color=colors[i], label=i, lw=lw[0])
        sns.kdeplot(fxnY(subdata[y]), ax=g.ax_marg_y, vertical=True, legend=False, color=colors[i], label=i, lw=lw[1])
        g.ax_joint.plot(fxnX(subdata[x]), fxnY(subdata[y]), "o", ms=2, color=colors[i], label=i)
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




def volcano(df, ax, FC=0.1, FDR=2, stage="", method="", xcol='logFC', ycol='log10fdr', y2log=False, x2log=False, cinsig='blue', csig='red', alpha=0.1):
    x = df[xcol]
    y = df[ycol]
    if y2log:
        y = -1*np.log10(y)
    if x2log:
        x = -1*np.log2(x)
    gateup = (x.abs()>=FC) & (y>=FDR)
    gatedown = (x.abs()<FC) | (y<FDR)
    ax.scatter(x=x[gatedown], y=y[gatedown], c=cinsig, alpha=alpha) 
    ax.scatter(x=x[gateup], y=y[gateup], c=csig, alpha=alpha)
    ax.set_title(method + ':\n' + stage + '\nAbs Log2FC >= '+str(FC) + '\nFDR <= ' + str(100*(10**(-FDR))) + '%')


def pvalDistro(df, ax, stage='', method='', pval='PValue'):
    sns.distplot(df.dropna()[pval], ax=ax)
    ax.set_title(method + ':\n' + stage)




def marginalized_volcano(df, ax, FC=0.1, FDR=2, stage="", method="", xcol='logFC', ycol='log10fdr', y2log=False, x2log=False, alpha=0.1):
    x = df[xcol]
    y = df[ycol]
    if y2log:
        y = -1*np.log10(y)
    if x2log:
        x = -1*np.log2(x)
    gateup = (x.abs()>=FC) & (y>=FDR)
    gatedown = (x.abs()<FC) | (y<FDR)
    ax.scatter(x=x[gatedown], y=y[gatedown], c='blue', alpha=alpha) 
    ax.scatter(x=x[gateup], y=y[gateup], c='red', alpha=alpha)
    ax.set_title(method + ':\n' + stage + '\nAbs Log2FC >= '+str(FC) + '\nFDR <= ' + str(100*(10**(-FDR))) + '%')



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









def create_margin_plot_parameters(left=0.1, width=0.65, bottom=0.1, height=0.65, spacing=0.005, marginToMainProp=0.1):# definitions for the axes
    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom + height + spacing, width, marginToMainProp*height]
    rect_histy = [left + width + spacing, bottom, marginToMainProp*height, height]
    return rect_scatter, rect_histx, rect_histy
    
    
def create_margin_plot_template(rect_scatter, rect_histx, rect_histy, figsize=(8,8)):
    # start with a square Figure
    fig = plt.figure(figsize=figsize)

    ax = fig.add_axes(rect_scatter)
    ax_histx = fig.add_axes(rect_histx, sharex=ax)
    ax_histy = fig.add_axes(rect_histy, sharey=ax)

    return (ax, ax_histx, ax_histy)

    
def scatter_hist(x, y, ax, ax_histx, ax_histy, color, alpha, binwidth = 0.25):
    # no labels
    ax_histx.tick_params(axis="x", labelbottom=False)
    ax_histy.tick_params(axis="y", labelleft=False)

    # the scatter plot:
    ax.scatter(x, y, c=color, alpha=alpha)

    # now determine nice limits by hand:
    xymax = max(np.max(np.abs(x)), np.max(np.abs(y)))
    lim = (int(xymax/binwidth) + 1) * binwidth

    bins = np.arange(-lim, lim + binwidth, binwidth)
    ax_histx.hist(x, bins=bins)
    ax_histy.hist(y, bins=bins, orientation='horizontal')




def single_marginalized_volcano(df, FC=1, FDR=2, stage="", method="", xcol='logFC', ycol='log10fdr',
                                y2log=False, x2log=False, alpha=0.02, csig='red', cinsig='blue', figsize=(8,8),
                                left=0.1, width=0.65, bottom=0.1, height=0.65,
                                spacing=0.005, xnbins=50, ynbins=50, xbinwidth = None, ybinwidth = None,
                                xhistcolor="black", yhistcolor="black", marginToMainProp=0.1,
                                xColForGates=None, yColForGates=None):
    '''
    df      =   pandas dataframe
    FC      =   cutoff for abs log2 fold change. default = 1
    FDR     =   cutoff for BH FDR. default = 2 (1\%).
    stage   =   optional words for title
    method  =   optional words for title
    xcol    =   df key for x-axis. default = logFC
    ycol    =   def key for y-axis. default = log10fdr
    y2log   =   bool. whether or not to take log10 of ycol. Default: False (normally already done).
    x2log   =   bool. whether or not to take log2 of ycol. Default: False (normally already done)
    alpha   =   transparency level of dots. default = 0.02
    csig    =   color for dots above cutoffs. default = red.
    cinsig  =   color for dots below 1 or both cutoffs. default = blue.
    figsize =   tuple figsize for matplotlib. default = (8,8)
    left    =   float for constructing plot. default = 0.1
    width   =   float for constructing plot. default = 0.65
    bottom  =   float for constructing plot. default = 0.1
    height  =   float for constructing plot. default = 0.65
    spacing =   float for constructing plot. default = 0.005
    xnbins   =   number of bins to use on x hist. default = 50.
    ynbins   =   number of bins to use in y hist. default = 50.
    xbinwidth =  float for width of bins. default - None. Uses nbins instead. Suggestion to try = 1. When xlims and ylims very differnt using the same binwidth for each  can look funny.
    ybinwidth =  float for width of bins. default - None. Uses nbins instead. Suggestion to try = 0.25. When xlims and ylims very differnt using the same binwidth for each  can look funny.
    xhistcolor = color of x hist bars. default = grey.
    yhistcolor = color of y hist bars. default = grey.    
    marginToMainProp = float between 0 and 1. proportional heights of margin plots relative to main plot. default = 0.1 (10%).
    xColForGates    =   Default is to use xcol.
    yColForGates    =   Default is to use ycol.
    Use seaborn sns.set_style to get dark grid.
    '''
    rect_scatter, rect_histx, rect_histy = create_margin_plot_parameters(left, width, bottom, height, spacing, marginToMainProp)
    ax, ax_histx, ax_histy = create_margin_plot_template(rect_scatter, rect_histx, rect_histy, figsize=figsize)

    x = df[xcol]
    y = df[ycol]

        
    if y2log:
        y = -1*np.log10(y)
    if x2log:
        x = -1*np.log2(x)


    # Create cutoff gates
    if xColForGates is None:
        xColForGates = xcol
    if yColForGates is None:
        yColForGates = ycol
    gateup = (df[xColForGates].abs() >= FC) & (df[yColForGates] >= FDR)
    gatedown = (df[xColForGates].abs() < FC) | (df[yColForGates] < FDR)


    # no labels
    ax_histx.tick_params(axis="x", labelbottom=False)
    ax_histy.tick_params(axis="y", labelleft=False)

    # the scatter plot:
    ax.scatter(x=x[gatedown], y=y[gatedown], c=cinsig, alpha=alpha) 
    ax.scatter(x=x[gateup], y=y[gateup], c=csig, alpha=alpha)
    #ax.set_title(method + ':\n' + stage + '\nAbs Log2FC >= '+str(FC) + '\nFDR <= ' + str(100*(10**(-FDR))) + '%')

    # now determine nice limits by hand:
    if xbinwidth is not None:
        xbins = np.arange(x.min(), x.max(), xbinwidth)
    else:
        xbins = np.linspace(x.min(), x.max(), xnbins)

    if ybinwidth is not None:
        ybins = np.arange(y.min(), y.max(), ybinwidth)
    else:
        ybins = np.linspace(y.min(), y.max(), ynbins)


    ## The histograms
    ax_histx.hist(x, bins=xbins, color=xhistcolor)
    ax_histy.hist(y, bins=ybins, color=yhistcolor, orientation='horizontal')






def single_marginalized_scatter_any_ax(df, ax, FC=1, FDR=2, stage="", method="", xcol='logFC', ycol='log10fdr',
                                    y2log=False, x2log=False, alpha=0.02, csig='red', cinsig='blue', 
                                    xnbins=50, ynbins=50, xbinwidth = None, ybinwidth = None,
                                    xhistcolor="black", yhistcolor="black", marginToMainProp=0.35,
                                       xColForGates=None, yColForGates=None):
    '''
    df      =   pandas dataframe
    FC      =   cutoff for abs log2 fold change. default = 1
    FDR     =   cutoff for BH FDR. default = 2 (1\%).
    stage   =   optional words for title
    method  =   optional words for title
    xcol    =   df key for x-axis. default = logFC
    ycol    =   def key for y-axis. default = log10fdr
    y2log   =   bool. whether or not to take log10 of ycol. Default: False (normally already done).
    x2log   =   bool. whether or not to take log2 of ycol. Default: False (normally already done)
    alpha   =   transparency level of dots. default = 0.02
    csig    =   color for dots above cutoffs. default = red.
    cinsig  =   color for dots below 1 or both cutoffs. default = blue.
    xnbins   =   number of bins to use on x hist. default = 50.
    ynbins   =   number of bins to use in y hist. default = 50.
    xbinwidth =  float for width of bins. default - None. Uses nbins instead. Suggestion to try = 1. When xlims and ylims very differnt using the same binwidth for each  can look funny.
    ybinwidth =  float for width of bins. default - None. Uses nbins instead. Suggestion to try = 0.25. When xlims and ylims very differnt using the same binwidth for each  can look funny.
    xhistcolor = color of x hist bars. default = grey.
    yhistcolor = color of y hist bars. default = grey.    
    marginToMainProp = float between 0 and 1. proportional heights of margin plots relative to main plot. default = 0.35 (35%).
    xColForGates    =   Default is to use xcol.
    yColForGates    =   Default is to use ycol.
    Use seaborn sns.set_style() to get "darkgrid', "white", "whitegrid", etc.


    This function requires you to first define a (fig, ax) -- for example, using:
        fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(16,16))

    Thus any axis can be provided to this one, and it may ultimately be the more useful function for making figures.
    
    '''

    # Define X and Y
    x = df[xcol]
    y = df[ycol]

        
    # Optionally log X and Y
    if y2log:
        y = -1*np.log10(y)
    if x2log:
        x = -1*np.log2(x)

    # Create cutoff gates
    if xColForGates is None:
        xColForGates = xcol
    if yColForGates is None:
        yColForGates = ycol
    gateup = (df[xColForGates].abs() >= FC) & (df[yColForGates] >= FDR)
    gatedown = (df[xColForGates].abs() < FC) | (df[yColForGates] < FDR)


    # the scatter plot:
    ax.scatter(x=x[gatedown], y=y[gatedown], c=cinsig, alpha=alpha) 
    ax.scatter(x=x[gateup], y=y[gateup], c=csig, alpha=alpha)
    
    # create new axes on the right and on the top of the current axes
    # The first argument of the new_vertical(new_horizontal) method is
    # the height (width) of the axes to be created in inches.
    divider = make_axes_locatable(ax)
    axHistx = divider.append_axes("top", marginToMainProp, pad=0.1, sharex=ax)
    axHisty = divider.append_axes("right", marginToMainProp, pad=0.1, sharey=ax)

    # make some labels invisible
    axHistx.xaxis.set_tick_params(labelbottom=False)
    axHisty.yaxis.set_tick_params(labelleft=False)

    # now determine nice limits by hand:
    if xbinwidth is not None:
        xbins = np.arange(x.min(), x.max(), xbinwidth)
    else:
        xbins = np.linspace(x.min(), x.max(), xnbins)

    if ybinwidth is not None:
        ybins = np.arange(y.min(), y.max(), ybinwidth)
    else:
        ybins = np.linspace(y.min(), y.max(), ynbins)

    # Create Hist
    axHistx.hist(x, bins=xbins, color=xhistcolor)
    axHisty.hist(y, bins=ybins, color=yhistcolor, orientation='horizontal')

    # the xaxis of axHistx and yaxis of axHisty are shared with axScatter,
    # thus there is no need to manually adjust the xlim and ylim of these
    # axis.

    axHistx.set_yticks([0, 50, 100])

    axHisty.set_xticks([0, 50, 100])



def single_marginalized_volcano_any_ax(df, ax, FC=1, FDR=2, stage="", method="", xcol='logFC', ycol='log10fdr',
                                    y2log=False, x2log=False, alpha=0.02, csig='red', cinsig='blue', 
                                    xnbins=50, ynbins=50, xbinwidth = None, ybinwidth = None,
                                    xhistcolor="black", yhistcolor="black", marginToMainProp=0.35,
                                       xColForGates=None, yColForGates=None):
    '''
    df      =   pandas dataframe
    FC      =   cutoff for abs log2 fold change. default = 1
    FDR     =   cutoff for BH FDR. default = 2 (1\%).
    stage   =   optional words for title
    method  =   optional words for title
    y2log   =   bool. whether or not to take log10 of ycol. Default: False (normally already done).
    x2log   =   bool. whether or not to take log2 of ycol. Default: False (normally already done)
    alpha   =   transparency level of dots. default = 0.02
    csig    =   color for dots above cutoffs. default = red.
    cinsig  =   color for dots below 1 or both cutoffs. default = blue.
    xnbins   =   number of bins to use on x hist. default = 50.
    ynbins   =   number of bins to use in y hist. default = 50.
    xbinwidth =  float for width of bins. default - None. Uses nbins instead. Suggestion to try = 1. When xlims and ylims very differnt using the same binwidth for each  can look funny.
    ybinwidth =  float for width of bins. default - None. Uses nbins instead. Suggestion to try = 0.25. When xlims and ylims very differnt using the same binwidth for each  can look funny.
    xhistcolor = color of x hist bars. default = grey.
    yhistcolor = color of y hist bars. default = grey.    
    marginToMainProp = float between 0 and 1. proportional heights of margin plots relative to main plot. default = 0.35 (35%).

    Use seaborn sns.set_style() to get "darkgrid', "white", "whitegrid", etc.

    This function differs from single_marginalized_volcano() in a few surface-level ways.
    Mainly, it requires you to first define a (fig, ax) -- for example, using:
        fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(16,16))

    Thus any axis can be provided to this one, and it may ultimately be the more useful function for making figures.
    
    '''
    single_marginalized_scatter_any_ax(df=df, ax=ax, FC=FC, FDR=FDR, stage=stage, method=method,
                                       xcol='logFC', ycol='log10fdr', xColForGates='logFC', yColForGates='log10fdr',
                                       y2log=y2log, x2log=x2log, alpha=0.02, csig=csig, cinsig=cinsig,
                                       xnbins=xnbins, ynbins=ynbins, xbinwidth = xbinwidth, ybinwidth = ybinwidth,
                                       xhistcolor=xhistcolor, yhistcolor=yhistcolor,
                                       marginToMainProp=marginToMainProp)
    


def single_marginalized_MA_volcano_any_ax(df, ax, FC=1, FDR=2, stage="", method="",
                                    y2log=False, x2log=False, alpha=0.02, csig='red', cinsig='blue', 
                                    xnbins=50, ynbins=50, xbinwidth = None, ybinwidth = None,
                                    xhistcolor="black", yhistcolor="black", marginToMainProp=0.35):
    '''
    This returns an MA plot of logCPM vs logFC that is colored by log10fdr and logFC cutoffs that
        are often used in Volcano plots.

    For more flexibility on what is plotted on the X and Y axes, and how dots are colored,
        use single_marginalized_scatter_any_ax

    df      =   pandas dataframe
    FC      =   cutoff for abs log2 fold change. default = 1
    FDR     =   cutoff for BH FDR. default = 2 (1\%).
    stage   =   optional words for title
    method  =   optional words for title
    y2log   =   bool. whether or not to take log10 of ycol. Default: False (normally already done).
    x2log   =   bool. whether or not to take log2 of ycol. Default: False (normally already done)
    alpha   =   transparency level of dots. default = 0.02
    csig    =   color for dots above cutoffs. default = red.
    cinsig  =   color for dots below 1 or both cutoffs. default = blue.
    xnbins   =   number of bins to use on x hist. default = 50.
    ynbins   =   number of bins to use in y hist. default = 50.
    xbinwidth =  float for width of bins. default - None. Uses nbins instead. Suggestion to try = 1. When xlims and ylims very differnt using the same binwidth for each  can look funny.
    ybinwidth =  float for width of bins. default - None. Uses nbins instead. Suggestion to try = 0.25. When xlims and ylims very differnt using the same binwidth for each  can look funny.
    xhistcolor = color of x hist bars. default = grey.
    yhistcolor = color of y hist bars. default = grey.    
    marginToMainProp = float between 0 and 1. proportional heights of margin plots relative to main plot. default = 0.35 (35%).


    Use seaborn sns.set_style() to get "darkgrid', "white", "whitegrid", etc.

    This function differs from single_marginalized_volcano() in a few surface-level ways.
    Mainly, it requires you to first define a (fig, ax) -- for example, using:
        fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(16,16))

    Thus any axis can be provided to this one, and it may ultimately be the more useful function for making figures.
    
    '''
    single_marginalized_scatter_any_ax(df=df, ax=ax, FC=FC, FDR=FDR, stage=stage, method=method,
                                       xcol='logCPM', ycol='logFC', xColForGates='logFC', yColForGates='log10fdr',
                                       y2log=y2log, x2log=x2log, alpha=0.02, csig=csig, cinsig=cinsig,
                                       xnbins=xnbins, ynbins=ynbins, xbinwidth = xbinwidth, ybinwidth = ybinwidth,
                                       xhistcolor=xhistcolor, yhistcolor=yhistcolor,
                                       marginToMainProp=marginToMainProp)






colorlist = ['red','orange','yellow','green','blue', 'indigo', 'violet', 'purple', 'pink', 'black', 'magenta', 'cyan', 'grey']

def cor_cluster_map(df, columns=None, distancemethod='ward', nrowclusts=2, ncolclusts=2, nrows='all',
                    colors = colorlist,
                cmap='coolwarm', center=None, vmin=None, vmax=None, colSameAsRow=False):
    '''
    df              =   pandas dataframe with gene name index and logFC columns from diff conditions.
    columns         =   list of columns to use/subset from df. Default is to use all as is.
    distancemethod  =   Default: ward. Other: centroid, single, average, median, complete, weighted.
                        Ward's method has worked best in my hands in some limited testing.
                        It is also the default method for hierarchical clustering in scikitlearn.
                        I believe it is also the default clustering for clustermap.......?
    nrowclusts      =   number of clusters to identify and color in rows of df. default = 2.
    ncolclusts      =   number of clusters to identify and color in cols of df. default = 2.
    nrows            =   integer number of genes to sample. making the plot takes a long time. so exploratory stuff should be small samples. Default = "all".
    colors          =   color list to sample from for cluster colors. default list has 10 colors, so can support 10 clusters by default.
                        ['red','orange','yellow','green','blue', 'indigo', 'violet', 'purple', 'pink', 'black']
    cmap            =   color scheme for heatmap. Default="coolwarm". Many others - see matplotlib and seaborn doc.
                        Recommended tries: vlag, magma.
    center          =   As in sns.heatmap().
    vmin            =   As in sns.heatmap().
    vmax            =   As in sns.heatmap().

    '''
    ## Create separate object
    df = df.copy()
    
    ## Optional subsetting
    if columns is not None:
        df = df[columns]
        

    ## Make gene-by-gene correlation matrix
    if nrows == "all":
        cordata = df.transpose().corr()
    elif nrows < df.shape[0]:
        cordata = df.sample(nrows).transpose().corr()
    else:
        return 'nrows needs to either be less than the number of rows in df, or "all"'

    return cluster_map(cordata, columns=None, distancemethod=distancemethod,
                nrowclusts=nrowclusts, ncolclusts=ncolclusts, nrows='all',
                colors = colors, cmap=cmap, center=center, vmin=vmin, vmax=vmax,
                colSameAsRow=True)
    


def cluster_map(df, columns=None, distancemethod='ward', nrowclusts=2, ncolclusts=2, nrows='all',
                    colors = colorlist,
                cmap='coolwarm', center=None, vmin=None, vmax=None, colSameAsRow=False):
    '''
    df              =   pandas dataframe with gene name index and logFC columns from diff conditions.
    columns         =   list of columns to use/subset from df. Default is to use all as is.
    distancemethod  =   Default: ward. Other: centroid, single, average, median, complete, weighted.
                        Ward's method has worked best in my hands in some limited testing.
                        It is also the default method for hierarchical clustering in scikitlearn.
                        I believe it is also the default clustering for clustermap.......?
    nrowclusts      =   number of clusters to identify and color in rows of df. default = 2.
    ncolclusts      =   number of clusters to identify and color in cols of df. default = 2.
    nrows            =   integer number of genes to sample. making the plot takes a long time. so exploratory stuff should be small samples. Default = "all".
    colors          =   color list to sample from for cluster colors. default list has 10 colors, so can support 10 clusters by default.
                        ['red','orange','yellow','green','blue', 'indigo', 'violet', 'purple', 'pink', 'black']
    cmap            =   color scheme for heatmap. Default="coolwarm". Many others - see matplotlib and seaborn doc.
                        Recommended tries: vlag, magma.
    center          =   As in sns.heatmap().
    vmin            =   As in sns.heatmap().
    vmax            =   As in sns.heatmap().
    colSameAsRow    =   When given a matrix that has the same rows and columns and identical across the diagonal, use this to avoid computing the same linkages twice.
                        An example is a correlation matrix.
                        Use/see cor_cluster_map() that wraps over this function to easily produce correlation clusters.
    '''

    ## Create separate object
    df = df.copy()
    
    ## Optional subsetting
    if columns is not None:
        df = df[columns]
        

    ## Make gene-by-gene correlation matrix
    if nrows == "all":
        cordata = df.transpose().corr()
    elif nrows < df.shape[0]:
        cordata = df.sample(nrows).transpose().corr()
    else:
        return 'nrows needs to either be less than the number of rows in df, or "all"'


    ## Learn clustering of correlation matrix, which will be the same for rows and cols
    row_linkage = linkage(df, method=distancemethod, metric='euclidean')
    col_linkage = row_linkage if colSameAsRow else linkage(df.transpose(), method=distancemethod, metric='euclidean')
    
    ## Cutting tree
    row_clusters = cut_tree(row_linkage, n_clusters=nrowclusts)
    col_clusters = row_clusters if colSameAsRow else cut_tree(col_linkage, n_clusters=ncolclusts)

    ## Create color vector based on cluster assignments
    row_colors = list([colors[e] for e in row_clusters.flatten()])
    col_colors = row_colors if colSameAsRow else list([colors[e] for e in col_clusters.flatten()])

    ## Generate cluster map
    sns.clustermap(df, cmap=cmap, center=center, vmin=vmin, vmax=vmax, row_colors=row_colors, 
                   col_colors=col_colors, row_linkage=row_linkage, col_linkage=col_linkage)

    ## Return df with cluster numbers
    df['cluster'] = row_clusters.flatten()
    df['color'] = row_colors
    return df































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
