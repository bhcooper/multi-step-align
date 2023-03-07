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

args = parser.parse_args()

corelen = args.max_core_length
r0file = args.R0_input
rnfiles = args.RN_input_and_cycle[::2]
cycles = np.array(args.RN_input_and_cycle[1::2]).astype(float)

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

def center(matrix):
    bound = np.nanmax(np.abs(matrix))
    print("abs max: " + str(bound))
    bound = 0.75
    print("Bound set to +/- " + str(bound) + " for matplotlib's ""RdBu"" cmap")
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

R2filter = coreE.index.get_level_values('round') == 2
scalingFactor = np.log(coreE[R2filter]).mean()/np.log(coreE[~R2filter]).mean()
# scalingFactor = 2
print(f'Scaling Factor: {"%.2f" % scalingFactor}')
coreE = -np.log(coreE)
coreE[R2filter] = coreE[R2filter]/scalingFactor

print("Calculating flanking enrichment . . .")
R0flanks = R0groups.apply(sumFlank)
RNflanks = RNgroups.apply(sumFlank)

# Remove a core if no information at a flanking position (could be from pseudopalindomic core)
argfilter = RNflanks.groupby(level=('core', 'shift', 'strand')).agg(anyZero)
# print("Removed Shifts:")
# print(argfilter[argfilter==True])
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

# Print outliers
# print(coreE[(coreE - (Q3 + IQR * 1.5)) > 0])
# print(coreE[(coreE - (Q1 - IQR * 1.5)) < 0])

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

# R0flanks = R0flanks.sum(axis=0, level="strand")
flankE = RNflanks.divide(R0flanks)

flankE = flankE.loc[ucores]
R2filter = flankE.index.get_level_values('round') == 2
flankE = flankE.replace(0, np.nan)
# Nullify positions affected by nan
flankE[flankE.isna().groupby(axis=1, level='pos').any()] = np.nan

# Normalize each position against max per position
# maxcols = flankE[flankE.mean(axis=0).groupby(level=0).idxmax()]
# maxcols.columns = maxcols.columns.droplevel('base')
# flankE = flankE.divide(maxcols, axis=0, level=0)

flankE = -np.log(flankE)
flankE[R2filter] = flankE[R2filter]/scalingFactor
flankE = flankE.subtract(flankE.groupby(axis=1, level=0).mean(), axis=0, level=0)

flankE.name = "ΔΔG/RT"
flankE_mean = flankE.groupby(level='core').mean()
flankE_CI = flankE.groupby(level='core').std() * 1.96 / np.sqrt(np.isfinite(flankE).groupby(level='core').sum())

# Output core info
fig, ax = plt.subplots(figsize=(5, 15))
ax = sns.violinplot(x=coreE.name, y='core', data=coreE[ucores[1:]].reset_index(), linewidth = 0, color="C0", cut=0, inner='stick')
# ax = sns.violinplot(x=coreE.name, y='core', data=coreE.xs(1, level='round').reset_index(), linewidth = 0, color="lightskyblue", cut=0)
# ax = sns.violinplot(x=coreE.name, y='core', data=coreE.xs(2, level='round').reset_index(), linewidth = 0, color="lightcoral", cut=0)
plt.setp(ax.collections, alpha=0.4)
ax.set_yticklabels(list(ucores[1:]), rotation="horizontal", fontname="Courier New", fontsize=12)
ax.errorbar(y=np.arange(len(coreE_mean[ucores[1:]])), x = coreE_mean[ucores[1:]], xerr=coreE_CI[ucores[1:]], capsize=3, linewidth=0, elinewidth=1.0, ecolor="C0")
ax.set_ylim(len(ucores)-1, -1)
# with open(rndir + "_core_ddG.pkl", 'wb') as f:
    # pickle.dump([coreE, coreE_mean, coreE_CI, ucores, rndir], f)

plt.savefig(rndir + "_core_ddG" + imgext, bbox_inches="tight", dpi=600)
plt.close()

coreE_mean = coreE_mean.to_frame()
coreE_mean['CI'] = coreE_CI
coreE_mean.to_csv(rndir + "_core_ddG.tsv", sep='\t', header='ΔΔG/RT\t95%_CI')
np.exp(-coreE_mean.iloc[:,0]).to_csv(rndir + "_core_RelE.tsv", sep='\t', header='RelE')

# Output flank info
posLabels = ["-" + str(s) + "\n" + b if b == "C" else "\n" + b for b, s in zip(tools.N*flankLen, np.repeat(np.arange(flankLen)[::-1] + 1, 4))]
posLabels += ["+" + str(s) + "\n" + b if b == "C" else "\n" + b for b, s in zip(tools.N*flankLen, np.repeat(np.arange(flankLen) + 1, 4))]
print("Max edge accross all cores")
tools.plotGrid(rndir + "_edge_ddG" + imgext, center(flankE_mean.loc[ucores].values), xticks=posLabels, yticks=ucores, gridstridex=4, gridstridey=None, vmin=-0.25, vmax=1.25, vline=flankLen*4-0.5, cmap="RdBu_r")
print("Max edge for best cores")
ref = ucores[0]
# ref = "CGAAACA"
colmean = center(flankE.loc[ref]).mean().values.reshape(1, -1)
flankExample = np.concatenate((center(flankE.loc[ref].values), colmean), axis=0)
yticks = [f'R{t[2]}{t[1]}' for t in flankE.loc[ref].index.tolist()]
yticks += ["Mean"]
tools.plotGrid(rndir + "_edge_ddG_" +   ref + "_left" + imgext, flankExample[:,12:36], xticks=posLabels[12:36], yticks = yticks, gridstridex=4, gridstridey=len(flankExample)-1,  vmin=-0.25, vmax=1.25, vline=(flankLen-3)*4-0.5, cmap="RdBu_r")
tools.plotGrid(rndir + "_edge_ddG_" + ref + "_right" + imgext, flankExample[:,36:-12], xticks=posLabels[36:-12], gridstridex=4, gridstridey=len(flankExample)-1,  vmin=-0.25, vmax=1.25, vline=-0.5, cmap="RdBu_r")
# tools.plotGrid(rndir + "_edge_ddG_allShifts.png", center(flankE.loc[ucores[0]].xs("+", level="strand").values[:,4:-4]), xticks=posLabels[4:-4], gridstridex=4, gridstridey=None,  vmin=-0.25, vmax=1.25, vline=flankLen*4-0.5-4, cmap="RdBu_r")

ltrim = 3
rtrim = 3

argsort = np.argsort(ucores)
left = center(flankE_mean.loc[ucores[argsort]].values[:,ltrim*4:flankLen*4])
xticks = posLabels[ltrim*4:flankLen*4]
figsize = np.array(left.transpose().shape)/4
figsize[0] /= 0.92
tools.plotGrid(rndir + "_edges_left_lexi" + imgext, left, xticks=xticks, yticks=[u.replace('-', '') for u in ucores[argsort]], gridstridex=4, gridstridey=len(ucores), vmin=-0.25, vmax=1.25, vline=(flankLen-ltrim)*4-0.5, cmap="RdBu_r")

argsort = np.argsort([u[::-1].replace('-', '') for u in ucores])
right = center(flankE_mean.loc[ucores[argsort]].values[:,flankLen*4:-rtrim*4])
xticks = posLabels[flankLen*4:-rtrim*4]
figsize = np.array(right.transpose().shape)/4
figsize[0] /= 0.92
tools.plotGrid(rndir + "_edges_right_lexi" + imgext, right, xticks=xticks, yticks=[u.replace('-', '') for u in ucores[argsort]], gridstridex=4, gridstridey=len(ucores), vmin=-0.25, vmax=1.25, vline=-0.5, cmap="RdBu_r", tick_left=True)

# Output flank contributions

flankPos = [b + " [-" + str(s) + "]" for b, s in zip(tools.N*flankLen, np.repeat(np.arange(flankLen)[::-1] + 1, 4))]
flankPos += [b + " [+" + str(s) + "]" for b, s in zip(tools.N*flankLen, np.repeat(np.arange(flankLen) + 1, 4))]
flankE_mean_relmax = flankE_mean - flankE_mean.groupby(axis=1, level='pos').min()
flankE_mean_relmax.to_csv(rndir + "_edge_ddG.tsv", sep='\t', header=flankPos)
flankE_CI.to_csv(rndir + "_edge_ddG_CI.tsv", sep='\t', header=flankPos)

# Calculate maximum contributions (A little noisy, maybe better to try range of average only)
# flankRange = flankE_mean.mean()
# flankRange = flankRange.groupby(level='pos').max() - flankRange.groupby(level='pos').min()
# fig, ax = plt.subplots(figsize=(8,5))
# width=0.7
# x = np.arange(len(flankRange))
# ax.bar(x, flankRange, width, label=rndir)
# # ax.errorbar(x, flankRange_mean, flankRange_var, capsize=2, linewidth=0, elinewidth=1, ecolor="black")
# xlabels = ["-" + str(x) for x in np.arange(1,flankLen+1)[::-1]] + ["+" + str(x) for x in np.arange(1,flankLen+1)]
# ax.set_xticks(x)
# ax.set_xticklabels(xlabels, fontsize=10)
# outdata = np.array([xlabels, flankRange]).transpose()
# tools.saveTable(rndir + "_edge_range.tsv", outdata, header="position\trange")
# plt.savefig(rndir + "_edge_range" + imgext, bbox_inches="tight", dpi=600)
# plt.close()

# # Calculate maximum contributions (A little noisy, maybe better to try range of average only)
# argmax = flankE_mean.mean(axis=0).groupby(level='pos').idxmax()
# argmax = pd.MultiIndex.from_tuples(argmax, names = ['pos', 'base'])
# argmin = flankE_mean.mean(axis=0).groupby(level='pos').idxmin()
# argmin = pd.MultiIndex.from_tuples(argmin, names = ['pos', 'base'])
# flankRange = flankE[argmax].droplevel(level='base', axis=1) - flankE[argmin].droplevel(level='base', axis=1)
# flankRange_mean = flankRange.mean()
# flankRange_var = flankRange.std() ** 2
# fig, ax = plt.subplots(figsize=(8,5))
# width=0.7
# x = np.arange(len(flankRange_mean))
# ax.bar(x, flankRange_mean, width, label=rndir)
# ax.errorbar(x, flankRange_mean, flankRange_var, capsize=2, linewidth=0, elinewidth=1, ecolor="black")
# xlabels = ["-" + str(x) for x in np.arange(1,flankLen+1)[::-1]] + ["+" + str(x) for x in np.arange(1,flankLen+1)]
# ax.set_xticks(x)
# ax.set_xticklabels(xlabels, fontsize=10)
# outdata = np.array([xlabels, flankRange_mean, flankRange_var]).transpose()
# tools.saveTable(rndir + "_edge_range.tsv", outdata, header="position\trange\tvariance")
# plt.savefig(rndir + "_edge_range" + imgext, bbox_inches="tight", dpi=600)
# plt.close()

# # Calculate maximum contributions (Low noise, shuffled bases not counted, only span)
maxFlanks = flankE.groupby(axis=1, level='pos').max() - flankE.groupby(axis=1, level='pos').min()
maxFlanks_mean = maxFlanks.mean()
maxFlanks_var = maxFlanks.std() ** 2
# maxFlanks_CI = maxFlanks.std() * 1.96 / np.sqrt(np.sum(np.isfinite(maxFlanks), axis=0))
fig, ax = plt.subplots(figsize=(8,5))
width=0.7
x = np.arange(len(maxFlanks_mean))
ax.bar(x, maxFlanks_mean, width, label=rndir)
ax.errorbar(x, maxFlanks_mean, maxFlanks_var, capsize=2, linewidth=0, elinewidth=1, ecolor="black")
xlabels = ["-" + str(x) for x in np.arange(1,flankLen+1)[::-1]] + ["+" + str(x) for x in np.arange(1,flankLen+1)]
ax.set_xticks(x)
ax.set_xticklabels(xlabels, fontsize=10)
outdata = np.array([xlabels, maxFlanks_mean, maxFlanks_var]).transpose()
tools.saveTable(rndir + "_edge_range.tsv", outdata, header="position\trange\tvariance")
plt.savefig(rndir + "_edge_range" + imgext, bbox_inches="tight", dpi=600)
plt.close()

labelspace = 1.2
dendsize = 0.9

methods = ["single", "complete", "average", "weighted"]
# methods = ["average"]
for m in methods:
    left = center(flankE_mean.loc[ucores].values[:,ltrim*4:flankLen*4])
    linkage = shc.linkage(left, method=m, metric="cityblock", optimal_ordering=True)
    figsize = np.array(left.transpose().shape)/4
    figsize[0] /= 0.92
    fig = plt.figure(figsize=figsize)
    ax1 = fig.add_axes((0, 0, 1, 1))
    ax2 = fig.add_axes((1 + labelspace/figsize[0], 0, dendsize/figsize[0], 1))
    with plt.rc_context({'lines.linewidth': 1.0}):
        argsort = shc.dendrogram(linkage, orientation="right", ax=ax2, no_labels=True, color_threshold=0, above_threshold_color="black",)["leaves"][::-1]
    tools.plotGrid(None, left[argsort], xticks=posLabels[ltrim*4:flankLen*4], yticks=ucores[argsort], gridstridex=4, gridstridey=None, vmin=-0.25, vmax=1.25, vline=(flankLen-ltrim)*4-0.5, cmap="RdBu_r", ax=ax1)
    ax2.axis("off")
    plt.savefig(rndir + "_edges_left_clustered_" + m + imgext, bbox_inches="tight", dpi=600)
    plt.close()

    right = center(flankE_mean.loc[ucores].values[:,flankLen*4:-rtrim*4])
    linkage = shc.linkage(right, method=m, metric="cityblock", optimal_ordering=True)
    figsize = np.array(right.transpose().shape)/4
    figsize[0] /= 0.92
    fig = plt.figure(figsize=figsize)
    ax1 = fig.add_axes((0, 0, 1, 1))
    ax2 = fig.add_axes((-dendsize/figsize[0]-labelspace/figsize[0], 0, dendsize/figsize[0], 1))
    with plt.rc_context({'lines.linewidth': 1.0}):
        argsort = shc.dendrogram(linkage, orientation="left", ax=ax2, no_labels=True, color_threshold=0, above_threshold_color="black",)["leaves"][::-1]
    tools.plotGrid(None, right[argsort], xticks=posLabels[flankLen*4:-rtrim*4], yticks=ucores[argsort], gridstridex=4, gridstridey=None, vmin=-0.25, vmax=1.25, vline=-0.5, cmap="RdBu_r", ax=ax1)
    ax1.yaxis.tick_left()
    ax2.axis("off")
    plt.savefig(rndir + "_edges_right_clustered_" + m + imgext, bbox_inches="tight", dpi=600)
    plt.close()