#!/usr/bin/env python

import os
import numpy as np
import tools
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as shc
import pandas as pd
import seaborn as sns
import argparse

parser = argparse.ArgumentParser(description="Calculate the enrichment and ΔΔG/RT of differing cores and their flanking bp, then generate figures based on the results.")
parser.add_argument("max_core_length", type=int, help="Length of the longest core used for alignment")
parser.add_argument("R0_input")
parser.add_argument("RN_input_and_cycle", nargs="+", help="Pairs of aligned reads along with their corresponsing cycle number.")
parser.add_argument("-A", "--autoscale", action="store_true", default=False, help = "Predict how to rescale the enrichment of later inputs to match the scale of the first input given. Alternatively, the given cycle numbers will be used for rescaling.")
parser.add_argument("-B", "--bounds", default=None, type=float, required=False, help='Specifies the color bar bounds ± the mean for plotting the ΔΔG/RT of flanking positions. Helpful for consistency accross comparisons.')
args = parser.parse_args()

corelen = args.max_core_length
r0file = args.R0_input
rnfiles = args.RN_input_and_cycle[::2]
cycles = np.array(args.RN_input_and_cycle[1::2]).astype(float)
autoscale = args.autoscale
flankBound = args.bounds

labelspace = 0.9
dendsize = 0.9
cbarheight = 2.4

imgext = '.png'
# imgext = '.svg'

baseLookup = {"A":0, "C":1, "G":2, "T":3}
def encode1mers(seqs):
    X = np.zeros((len(seqs), len(seqs[0])*4), dtype=bool)
    for i, seq in enumerate(seqs):
        for j,c in enumerate(seq):
            if(not c == "N"): 
                X[i,j * 4 + baseLookup[c]] = True
    return X

def sumFlank(df):
    shift = df['shift'].iloc[0]
    left = (encode1mers(df['seq'].str[flankLen-shift:flankLen].values) * df['count'].values.reshape(-1,1)).sum(axis=0)
    right = (encode1mers(df['seq'].str[flankLen+corelen:flankLen+corelen+flankLen-shift].values) * df['count'].values.reshape(-1,1)).sum(axis=0)
    return np.append(left, right)

def anyZero(df):
    return np.any(np.vstack(df) == 0)

def center(matrix, bound):
    # bound = np.nanmax(np.abs(matrix))
    # print("abs max: " + str(bound))
    # bound = 0.75
    # print("Bound set to +/- " + str(bound) + " for matplotlib's ""RdBu"" cmap")
    return (matrix + bound)/(2*bound)

print("Reading input . . .")
R0 = pd.read_csv(r0file, sep='\t')
RNs = [pd.read_csv(f, sep='\t') for f in rnfiles]
for i, cycle in enumerate(cycles):
    RNs[i]['round'] = cycle
RN = pd.concat(RNs)

varLen = len(R0.iloc[0,0].replace("N", ""))
flankLen = varLen-corelen

rndir = rnfiles[0].split("_")[0]
tools.mkdir(rndir + "_analysis")
os.chdir(rndir + "_analysis")

R0['core'] = R0.iloc[:,0].str.slice(flankLen, flankLen + corelen)
RN['core'] = RN.iloc[:,0].str.slice(flankLen, flankLen + corelen)

print("Grouping shifts . . .")
R0groups = R0.groupby(['core', 'shift', 'strand'])
RNgroups = RN.groupby(['core', 'shift', 'strand', 'round'])

print("Calculating core enrichment . . .")
R0cores = R0groups['count'].sum()
RNcores = RNgroups['count'].sum()
coreE = RNcores.divide(R0cores)
maxCore = coreE.groupby(level=('core')).mean().idxmax()
coreE = coreE.divide(coreE.loc[maxCore])
coreE = coreE.reorder_levels(['core', 'shift', 'strand', 'round'])

scalingFactors = None
if(autoscale):
    logcoreE = np.log(coreE)
    scalingFactors = logcoreE.groupby("round").mean()
    scalingFactors = scalingFactors / scalingFactors.loc[cycles[0]]
    scalingFactors.name = 'Scaling Factors'
else:
    scalingFactors = pd.Series(cycles, index=pd.Index(cycles, name='round'), name='Scaling Factors')

print(scalingFactors)

coreE = -np.log(coreE)
coreE = coreE/scalingFactors

print("Calculating flanking enrichment . . .")
R0flanks = R0groups.apply(sumFlank)
RNflanks = RNgroups.apply(sumFlank)

# Remove a core if no information at a flanking position (could be from pseudopalindomic core)
argfilter = RNflanks.groupby(level=('core', 'shift', 'strand')).agg(anyZero)
coreE = coreE[~argfilter]
argfilter = coreE.groupby(level='core').count() < flankLen
print("Removed Cores:")
print('\n'.join(argfilter[argfilter==True].index))
coreE = coreE[argfilter[argfilter==False].index]

# Remove outliers (Tukey's method)
Q1 = coreE.groupby(level='core').quantile(0.25)
Q3 = coreE.groupby(level='core').quantile(0.75)
IQR = Q3 - Q1
coreE = coreE[(coreE - (Q3 + IQR * 1.5)) <= 0]
coreE = coreE[(coreE - (Q1 - IQR * 1.5)) >= 0]

# Calculate core means and 95% confidence intervals
coreE_mean = coreE.groupby(level='core').mean().sort_values()
coreE.name = "ΔΔG/RT"
ucores = coreE_mean.index
coreE_mean = coreE_mean.loc[ucores]
coreE = coreE.loc[ucores]
coreE_CI = coreE.groupby(level='core').std() * 1.96 / coreE.groupby(level='core').count().pow(1/2)

# Align flanking positions
R0flanks[:] = [np.concatenate((np.full((flankLen-shift)*4, np.nan), x, np.full(shift*4, np.nan))) for x, shift in zip(R0flanks, R0flanks.index.get_level_values('shift'))]
RNflanks[:] = [np.concatenate((np.full((flankLen-shift)*4, np.nan), x, np.full(shift*4, np.nan))) for x, shift in zip(RNflanks, RNflanks.index.get_level_values('shift'))]

R0flanks = pd.DataFrame(R0flanks.tolist(), index=R0flanks.index)
RNflanks = pd.DataFrame(RNflanks.tolist(), index=RNflanks.index)

cols = pd.MultiIndex.from_tuples([(int(c/4), c%4) for c in R0flanks.columns], names = ['pos', 'base'])
R0flanks.columns = cols
RNflanks.columns = cols

flankE = RNflanks.divide(R0flanks)

flankE = flankE.loc[ucores]
flankE = flankE.replace(0, np.nan)

# Nullify positions affected by nan
flankE[flankE.isna().groupby(axis=1, level='pos').any()] = np.nan

# Calculate enrichment and center
flankE = -np.log(flankE)
flankE = flankE.divide(scalingFactors, level='round', axis='index')
flankE = flankE.subtract(flankE.groupby(axis=1, level=0).mean(), axis=0, level=0)

flankE.name = "ΔΔG/RT"
flankE_mean = flankE.groupby(level='core').mean()
flankE_CI = flankE.groupby(level='core').std() * 1.96 / np.sqrt(np.isfinite(flankE).groupby(level='core').sum())

# Provide outputs based on the core enrichments
fig, ax = plt.subplots(figsize=(5, 15))
ax = sns.violinplot(x=coreE.name, y='core', data=coreE[ucores[1:]].reset_index(), linewidth = 0, color="C0", cut=0, inner='stick')
plt.setp(ax.collections, alpha=0.4)
ax.set_yticklabels(list(ucores[1:]), rotation="horizontal", fontname="Courier New", fontsize=12)
ax.errorbar(y=np.arange(len(coreE_mean[ucores[1:]])), x = coreE_mean[ucores[1:]], xerr=coreE_CI[ucores[1:]], capsize=3, linewidth=0, elinewidth=1.0, ecolor="C0")
ax.set_ylim(len(ucores)-1, -1)
plt.savefig(rndir + "_core_ddG" + imgext, bbox_inches="tight", dpi=600)
plt.close()

coreE_mean = coreE_mean.to_frame()
coreE_mean['CI'] = coreE_CI
coreE_mean.to_csv(rndir + "_core_ddG.tsv", sep='\t', header='ΔΔG/RT\t95%_CI')
np.exp(-coreE_mean.iloc[:,0]).to_csv(rndir + "_core_RelE.tsv", sep='\t', header='RelE')

# Provide outputs based on the flank enrichments

# flankBound = 0.75
if(not flankBound):
    flankBound = np.ceil((np.nanmax(np.abs(flankE_mean.loc[ucores].values))*100))/100

posLabels = ["-" + str(s) + "\n" + b if b == "C" else "\n" + b for b, s in zip(tools.N*flankLen, np.repeat(np.arange(flankLen)[::-1] + 1, 4))]
posLabels += ["+" + str(s) + "\n" + b if b == "C" else "\n" + b for b, s in zip(tools.N*flankLen, np.repeat(np.arange(flankLen) + 1, 4))]

plotmat = center(flankE_mean.loc[ucores].values, flankBound)
figsize = np.array(plotmat.transpose().shape)/4
fig = plt.figure(figsize=figsize)
ax1 = fig.add_axes((0, 0, 1, 1))
ax2 = fig.add_axes((1 + 1.4/figsize[0], 0.6/figsize[1], 0.3/figsize[0], cbarheight/figsize[1]))
tools.plotGrid(None, plotmat[::-1], ax=ax1, xticks=posLabels, yticks=ucores[::-1], gridstridex=4, gridstridey=None, vmin=0, vmax=1, vline=flankLen*4-0.5, cmap="RdBu_r")
plt.colorbar(mpl.cm.ScalarMappable(cmap='RdBu'), cax=ax2, pad=0, ticks=[0, 1])
ax2.set_yticklabels([-flankBound, flankBound], fontsize=18)
ax2.set_ylabel("$\it{-ΔΔG/RT}$", fontsize=18, rotation=270, labelpad=4, fontweight='bold')
plt.savefig(rndir + "_flank_ddG" + imgext, bbox_inches="tight", dpi=600)
plt.close()


ref = ucores[0]
colmean = center(flankE.loc[ref], flankBound).mean().values.reshape(1, -1)
flankExample = np.concatenate((center(flankE.loc[ref].values, flankBound), colmean), axis=0)
yticks = [f'R{int(t[2])}{t[1]}' for t in flankE.loc[ref].index.tolist()]
yticks += ["Mean"]
plotmat = flankExample[:,:36]
figsize = np.array(plotmat.transpose().shape)/4
fig = plt.figure(figsize=figsize)
ax1 = fig.add_axes((0, 0, 1, 1))
ax2 = fig.add_axes((1+0.7/figsize[0], 0, 1, 1))
ax3 = fig.add_axes((2+1.0/figsize[0], 0.6/figsize[1], 0.3/figsize[0], cbarheight/figsize[1]))
tools.plotGrid(None, plotmat[::-1], ax=ax1, xticks=posLabels[:36], yticks = yticks[::-1], gridstridex=4,  vmin=0, vmax=1, vline=(flankLen)*4-0.5, cmap="RdBu_r")
ax1.set_yticks([-0.5, 0.5, len(plotmat)-0.5], minor=True)
plotmat = flankExample[:,36:]
tools.plotGrid(None, plotmat[::-1], ax=ax2, xticks=posLabels[36:], gridstridex=4,  vmin=-0.25, vmax=1.25, vline=-0.5, cmap="RdBu_r")
ax2.set_yticks([-0.5, 0.5, len(plotmat)-0.5], minor=True)
plt.colorbar(mpl.cm.ScalarMappable(cmap='RdBu'), cax=ax3, pad=0, ticks=[0, 1])
ax3.set_yticklabels([-flankBound, flankBound], fontsize=18)
ax3.set_ylabel("$\it{-ΔΔG/RT}$", fontsize=18, rotation=270, labelpad=4, fontweight='bold')
plt.savefig(rndir + "_flank_ddG_" + ref + imgext, bbox_inches="tight", dpi=600)
plt.close()

flankPos = [b + " [-" + str(s) + "]" for b, s in zip(tools.N*flankLen, np.repeat(np.arange(flankLen)[::-1] + 1, 4))]
flankPos += [b + " [+" + str(s) + "]" for b, s in zip(tools.N*flankLen, np.repeat(np.arange(flankLen) + 1, 4))]
flankE_mean_relmax = flankE_mean - flankE_mean.groupby(axis=1, level='pos').min()
flankE_mean_relmax.to_csv(rndir + "_flank_ddG.tsv", sep='\t', header=flankPos)
flankE_CI.to_csv(rndir + "_flank_ddG_CI.tsv", sep='\t', header=flankPos)

maxFlanks = flankE.groupby(axis=1, level='pos').max() - flankE.groupby(axis=1, level='pos').min()
maxFlanks_mean = maxFlanks.mean()
maxFlanks_var = maxFlanks.std() ** 2
fig, ax = plt.subplots(figsize=(8,5))
width=0.7
x = np.arange(len(maxFlanks_mean))
ax.bar(x, maxFlanks_mean, width, label=rndir)
ax.errorbar(x, maxFlanks_mean, maxFlanks_var, capsize=2, linewidth=0, elinewidth=1, ecolor="black")
xlabels = ["-" + str(x) for x in np.arange(1,flankLen+1)[::-1]] + ["+" + str(x) for x in np.arange(1,flankLen+1)]
ax.set_xticks(x)
ax.set_xticklabels(xlabels, fontsize=10)
outdata = np.array([xlabels, maxFlanks_mean, maxFlanks_var]).transpose()
tools.saveTable(rndir + "_flank_range.tsv", outdata, header="position\trange\tvariance")
plt.savefig(rndir + "_flank_range" + imgext, bbox_inches="tight", dpi=600)
plt.close()


ltrim = 3
rtrim = 3
method = 'average'
argsort = None

left = center(flankE_mean.loc[ucores].values[:,ltrim*4:flankLen*4], flankBound)
linkage = shc.linkage(left, method=method, metric="cityblock", optimal_ordering=True)
figsize = np.array(left.transpose().shape)/4
figsize[0] /= 0.92
fig = plt.figure(figsize=figsize)
ax1 = fig.add_axes((0, 0, 1, 1))
ax2 = fig.add_axes((1 + labelspace/figsize[0], 0, dendsize/figsize[0], 1))
with plt.rc_context({'lines.linewidth': 1.0}):
    argsort = shc.dendrogram(linkage, orientation="right", ax=ax2, no_labels=True, color_threshold=0, above_threshold_color="black",)["leaves"]
tools.plotGrid(None, left[argsort], xticks=posLabels[ltrim*4:flankLen*4], yticks=ucores[argsort], gridstridex=4, gridstridey=None, vmin=-0, vmax=1, vline=(flankLen)*4-0.5-ltrim*4, cmap="RdBu_r", ax=ax1)
ax2.axis("off")

right = center(flankE_mean.loc[ucores].values[:,flankLen*4:-rtrim*4], flankBound)
linkage = shc.linkage(right, method=method, metric="cityblock", optimal_ordering=True)
ax3 = fig.add_axes((1+labelspace/figsize[0]+dendsize/figsize[0] + 0.3/figsize[0], 0, dendsize/figsize[0], 1))
ax4 = fig.add_axes((1 + 2*labelspace/figsize[0] + 2*dendsize/figsize[0] + 0.3/figsize[0], 0, 1, 1))
with plt.rc_context({'lines.linewidth': 1.0}):
    argsort = shc.dendrogram(linkage, orientation="left", ax=ax3, no_labels=True, color_threshold=0, above_threshold_color="black",)["leaves"]
tools.plotGrid(None, right[argsort], xticks=posLabels[flankLen*4:-rtrim*4], yticks=ucores[argsort], gridstridex=4, gridstridey=None, vmin=0, vmax=1, vline=-0.5, cmap="RdBu_r", ax=ax4)
ax4.yaxis.tick_left()
ax3.axis("off")
ax5 = fig.add_axes((2 + 2*labelspace/figsize[0] + 2*dendsize/figsize[0] + 0.5/figsize[0], 0.6/figsize[1], 0.3/figsize[0], cbarheight/figsize[1]))
plt.colorbar(mpl.cm.ScalarMappable(cmap='RdBu'), cax=ax5, pad=0, ticks=[0, 1])
ax5.set_yticklabels([-flankBound, flankBound], fontsize=18)
ax5.set_ylabel("$\it{-ΔΔG/RT}$", fontsize=18, rotation=270, labelpad=4, fontweight='bold')
plt.savefig(rndir + "_flank_clustered_" + method + imgext, bbox_inches="tight", dpi=600)
plt.close()