#!/usr/bin/env python

import sys
import os
import numpy as np
import tools
import pickle
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as shc

config = tools.loadConfig(sys.argv[1])
ladapter = config['ladapter']
radapter = config['radapter']
ladapter_rc = tools.rc([radapter])[0]
radapter_rc = tools.rc([ladapter])[0]
markovOrder = config['markovOrder']
corelen = config['corelen']

rndir = sys.argv[2]
r0kmerDist = tools.loadVariables(config['r0kmerDist'])
r0dist = r0kmerDist[1].values.transpose().flatten()

Xlist = []
Xcorelist = []
shifts = []
shiftnames = []
ucoreslist = []
seqlen = None

picklefile = rndir + "_savedScan.pkl"
if(os.path.exists(picklefile)):
    with open(picklefile, 'rb') as f:
        Xlist, Xcorelist, shifts, shiftnames, ucoreslist, seqlen = pickle.load(f)
    print("Loaded " + picklefile + "\n")
else:
    os.chdir(rndir)
    for i, rnf in enumerate([f for f in sorted(os.listdir(".")) if f.endswith(".tsv")]):

        shiftnames += [rnf[:-4].replace("_", " ")]

        print("Processing " + rnf + " . . .")
        rnseqs, rny = tools.readScores(rnf)
        seqlen = len(rnseqs[0])
        nwindows = seqlen - corelen + 1

        shift = int(shiftnames[-1].split()[-3][5:])
        shifts += [shift]

        # Encode cores
        rncores = np.array([seq[shift:shift+corelen] for seq in rnseqs])
        ucores = np.unique(rncores)
        corelookup = dict(zip(ucores, np.arange(len(ucores))))
        rncores = [corelookup[c] for c in rncores]
        ucores = [c.replace("-", "") for c in ucores]
        ucoreslist += [ucores]
        corelens = [len(c) for c in ucores]

        # Count cores
        rncorecounts = np.zeros(len(ucores))
        for c,y in zip(rncores, rny):
            rncorecounts[c] += y

        # Calculate core enrichments
        onRC = False
        if("-" in rnf):
            onRC = True
        Edist = tools.markov(ucores, r0kmerDist, markovOrder, shift, onRC)

        # Core enrichment
        Xcore = rncorecounts/Edist
        Xcorelist += [Xcore]

        r0d = np.array(r0dist)
        if("-" in shiftnames[-1]):
            r0d = r0d[::-1]
        r0X = np.zeros((len(ucores), (nwindows - 1) * 4))

        for i,l in enumerate(corelens):
            gap = corelen - l
            if (gap == 0):
                r0X[i] = np.append(r0d[:shift*4],r0d[(shift+l)*4:])
            if (gap > 0):
                r0X[i] = np.append(r0d[:shift*4],r0d[(shift+l)*4:-gap*4])

        rnedges = np.array([seq[:shift] + seq[shift+corelen:] for seq in rnseqs])

        filts = [rncores == i for i in np.arange(len(ucores))]
        rnedgelist = [rnedges[filt] for filt in filts]
        rnylist = [rny[filt] for filt in filts]
        rnedgelist = [tools.encode1mers(rnedges) for rnedges in rnedgelist]

        rnX = np.array([np.dot(rnedges.transpose(), rny) for rnedges, rny in zip(rnedgelist, rnylist)]).astype(float)
        # Nullify all bases if one is nan
        filt = np.argwhere(rnX == 0)
        filt[:,1] = (filt[:,1]/4).astype(int)
        for f in filt:
            rnX[f[0],f[1]*4:f[1]*4 + 4] = np.nan
        rnX /= r0X
        Xlist += [rnX]
    Xlist = np.array(Xlist, dtype=object)
    os.chdir("..")
    with open(picklefile, 'wb') as f:
            pickle.dump([Xlist, Xcorelist, shifts, shiftnames, ucoreslist, seqlen], f)

tools.mkdir(rndir + "_analysis")
os.chdir(rndir + "_analysis")

shiftnames = np.array(shiftnames)
shifts = np.array(shifts)
Xcorelist = np.array(Xcorelist, dtype=object)
ucoreslist = np.array(ucoreslist, dtype=object)

allcores = np.unique(np.concatenate([np.unique(c) for c in ucoreslist]))
ucores = tools.getOverlapSet(ucoreslist)

argsort = np.argsort(shifts, kind="stable")
shifts = shifts[argsort]
shiftnames = shiftnames[argsort]
Xlist = Xlist[argsort]
Xcorelist = Xcorelist[argsort]
ucoreslist = ucoreslist[argsort]


R2filter = np.array(["R2" in s for s in shiftnames])
plusfilter = np.array(["+" in s for s in shiftnames])

# Nullify shift of core if affected by nans in flank
argfilters = [np.array([np.any(np.isnan(x.astype(float))) for x in X]) for X in Xlist]

for i in range(len(Xcorelist)):
    Xcorelist[i][argfilters[i]] = np.nan

# Align all shifts
numshifts = np.max(shifts)
Xlist = [np.concatenate((np.full((len(X), (numshifts-shift)*4), np.nan), X, np.full((len(X), shift*4), np.nan)), axis=1) for shift,X in zip(shifts, Xlist)]
Xlist = np.array(Xlist, dtype=object) 

# Restrict to cores occuring in all shifts
ucores = tools.getOverlapSet(list(ucoreslist))
argfilters = [np.intersect1d(X, ucores, return_indices=True)[1] for X in ucoreslist]
Xlist = np.array([X[filt] for X, filt in zip(Xlist, argfilters)], dtype=float)
Xcorelist = np.array([X[filt] for X, filt in zip(Xcorelist, argfilters)], dtype=float)

# Remove a core if no information at a flanking position (could be from pseudopalindomic core)
argfilter = np.array([np.any(np.all(np.isnan(Xlist[:,i,:]), axis=0)) for i in range(len(ucores))])
print("Removed Cores:")
print("\n".join(ucores[argfilter]))
print()

ucores = ucores[~argfilter]
Xlist = Xlist[:,~argfilter,:]
Xcorelist = Xcorelist[:,~argfilter]

# Calculate core ddG/RT
Xcoremeans = np.nanmean(Xcorelist, axis=0)
argmaxcore = np.argmax(Xcoremeans)
Xcorelist /= Xcorelist[:,argmaxcore].reshape(-1,1)
Xcorelist = np.log(Xcorelist)
Xcorelist[R2filter] /= 2

# Get average core ddG/RT and savec
core_means = np.nanmean(Xcorelist, axis=0)
core_CI = np.nanstd(Xcorelist, axis=0) * 1.96 / np.sqrt(len(Xcorelist))
outdata = np.array([ucores, core_means, core_CI]).transpose()
outdata = outdata[np.argsort(-core_means)]
tools.saveTable(rndir + "_core_ddG.tsv", outdata, header="Core\tddG/RT\t0.95_CI")

# # # Plot ddG/RT violin plot
fig, ax = plt.subplots(figsize=(5, 15))
argsort = np.argsort(core_means)[:-1]
nonan = Xcorelist[:,argsort].transpose()
nonan = [-X[np.isfinite(X)] for X in nonan]

ax.violinplot(nonan, showextrema=False, vert=False)
ax.errorbar(y=np.arange(len(core_CI)-1)+1, x=-core_means[argsort], xerr=-core_CI[argsort], capsize=3, linewidth=0, elinewidth=1.0, ecolor="C0")
ax.set_yticks(np.arange(len(ucores)-1)+1)
ax.set_yticklabels(list(ucores[argsort]), rotation="horizontal", fontname="Courier New", fontsize=18)
plt.xticks(fontsize=16)
plt.ylim((0,len(ucores)))
plt.title(rndir, fontsize=18)
plt.savefig(rndir + "_core_ddG.png", bbox_inches="tight", dpi=600)

# Find ddG/RT per postion
for i in np.arange(0, Xlist.shape[2], 4):
    Xlist[:,:,i:i+4] /= np.mean(Xlist[:,:,i:i+4], axis=2)[:,:,np.newaxis]
Xlist = np.log(Xlist)
Xlist[R2filter] /= 2

# Find average ddG/RT per core, and save
edge_means = np.nanmean(Xlist, axis=0)
for i in np.arange(0, edge_means.shape[1], 4):
    edge_means[:,i:i+4] -= np.nanmean(edge_means[:,i:i+4], axis=1).reshape(-1,1)
edge_CI = np.nanstd(Xlist, axis=0) * 1.96 / np.sqrt(np.sum(np.isfinite(Xlist), axis=0))
outdata = np.concatenate((ucores.reshape(-1,1), edge_means), axis=1)
posLabels = [b + " [-" + str(s) + "]" for b, s in zip(tools.N*numshifts, np.repeat(np.arange(numshifts)[::-1] + 1, 4))]
posLabels += [b + " [+" + str(s) + "]" for b, s in zip(tools.N*numshifts, np.repeat(np.arange(numshifts) + 1, 4))]
tools.saveTable(rndir + "_edge_ddG.tsv", outdata, header="\t".join(["core"] + posLabels))
outdata = np.concatenate((ucores.reshape(-1,1), edge_CI), axis=1)
tools.saveTable(rndir + "_edge_ddG_CI.tsv", outdata, header="\t".join(["core"] + posLabels))

# Plot all shifts for one core
centered = np.array(Xlist[:,argmaxcore,:])
for i in np.arange(0, centered.shape[1], 4):
    centered[:,i:i+4] -= np.mean(centered[:,i:i+4], axis=1).reshape(-1,1)
bound = np.nanmax(np.abs(centered))
print(rndir + " abs max: " + str(bound))
# bound = 0.75
print("Bound set to +/- " + str(bound) + " for matplotlib's ""RdBu"" cmap")
centered = (centered + bound)/(2*bound)
centered = np.concatenate((centered, np.nanmean(centered, axis=0).reshape(1, -1)), axis=0)
posLabels = [" -" + str(s) + "\n" + b if b == "C" else "\n" + b for b, s in zip(tools.N*numshifts, np.repeat(np.arange(numshifts)[::-1] + 1, 4))]
posLabels += [" +" + str(s) + "\n" + b if b == "C" else "\n" + b for b, s in zip(tools.N*numshifts, np.repeat(np.arange(numshifts) + 1, 4))]
yticks = shiftnames
yticks = ["".join(np.take(t.split(), [3,2])) for t in yticks] + ["Mean"]
mid = numshifts*4
tools.plotGrid(rndir + "_" + ucores[argmaxcore] + "_all_shifts.png", centered, ylabel="", xticks=posLabels, yticks=yticks, gridstridex = 4, gridstridey = len(shiftnames), vmin=-0.25, vmax=1.25, vline=numshifts*4-0.5, cmap="RdBu")
tools.plotGrid(rndir + "_" + ucores[argmaxcore] + "_all_shifts_left.png", centered[:,:int(mid)], ylabel="", xticks=posLabels[:int(mid)], yticks=yticks, gridstridex = 4, gridstridey = len(shiftnames), vmin=-0.25, vmax=1.25, vline=mid-0.5, cmap="RdBu")
tools.plotGrid(rndir + "_" + ucores[argmaxcore] + "_all_shifts_right.png", centered[:,int(mid):], ylabel="", xticks=posLabels[int(mid):], gridstridex = 4, gridstridey = len(shiftnames), vmin=-0.25, vmax=1.25, vline=-0.5, cmap="RdBu")

# Subtract averages per position per core and plot
centered = np.array(edge_means)
bound = np.nanmax(np.abs(centered))
print(rndir + " abs max: " + str(bound))
# bound = 0.75
print("Bound set to +/- " + str(bound) + " for matplotlib's ""RdBu"" cmap")
centered = (centered + bound)/(2*bound)
tools.plotGrid(rndir + "_edge_ddG.png", centered, xticks=posLabels, yticks=ucores, gridstridex=4, gridstridey=None, vmin=-0.25, vmax=1.25, vline=numshifts*4-0.5, cmap="RdBu")
# Plot average only over all cores
tools.plotGrid(rndir + "_edge_ddG_average.png", np.mean(centered,axis=0).reshape(1, -1), xticks=posLabels, yticks=["7-mer Average"], gridstridex=4, gridstridey=None, vmin=-0.25, vmax=1.25, vline=numshifts*4-0.5, cmap="RdBu")

# Calculate maximum contributions
colave = np.zeros((len(edge_means), (numshifts)*2))
for i in np.arange((numshifts)*2):
    colave[:,i] = np.max(edge_means[:,i*4:(i+1)*4], axis=1) - np.min(edge_means[:,i*4:(i+1)*4], axis=1)

var = np.nanstd(colave, axis=0) ** 2
CI = np.nanstd(colave, axis=0) * 1.96 / np.sqrt(np.sum(np.isfinite(colave), axis=0))
colave_mean = np.nanmean(colave, axis=0)

fig, ax = plt.subplots(figsize=(8,5))
width=0.7
x = np.arange(len(colave_mean))
ax.bar(x, colave_mean, width, label=rndir)
# ax.violinplot(colave1, showextrema=False, vert=True)
ax.errorbar(x, colave_mean, var, capsize=2, linewidth=0, elinewidth=1, ecolor="black")
xlabels = ["-" + str(x) for x in np.arange(1,numshifts+1)[::-1]] + ["+" + str(x) for x in np.arange(1,numshifts+1)]

ax.set_xticks(x)
ax.set_xticklabels(xlabels, fontsize=10)

print(np.sum(colave_mean[5:-7]))
outdata = np.array([xlabels, colave_mean, var]).transpose()
tools.saveTable(rndir + "_edge_range.tsv", outdata, header="position\trange\tvariance")
plt.savefig(rndir + "_edge_range.png", bbox_inches="tight", dpi=600)
plt.close()

# Hierarchical clustering of plots / trimmed region
labelspace = 1.2
dendsize = 0.9
endtrim = numshifts
labels = posLabels[:endtrim*4]
temp = centered[:,:endtrim*4]

figsize = np.array(temp.transpose().shape)/4
figsize[0] /= 0.92
argsort = np.argsort(ucores)
fig = plt.figure(figsize=figsize)
ax1 = fig.add_axes((0, 0, 1, 1))
tools.plotGrid(None, temp[argsort], xticks=labels, yticks=ucores[argsort], gridstridex=4, gridstridey=None, vmin=-0.25, vmax=1.25, vline=numshifts*4-0.5, cmap="RdBu", ax=ax1)
plt.savefig(rndir + "_edges_left_lexi_sort.png", bbox_inches="tight", dpi=600)

# methods = ["single", "complete", "average", "weighted"]
methods = ["average"]
for m in methods:
    linkage = shc.linkage(temp[:,-4:], method=m, metric="cityblock", optimal_ordering=True)
    fig = plt.figure(figsize=figsize)
    ax1 = fig.add_axes((0, 0, 1, 1))
    ax2 = fig.add_axes((1 + labelspace/figsize[0], 0, dendsize/figsize[0], 1))
    with plt.rc_context({'lines.linewidth': 1.0}):
        argsort = shc.dendrogram(linkage, orientation="right", ax=ax2, no_labels=True, color_threshold=0, above_threshold_color="black",)["leaves"][::-1]
    tools.plotGrid(None, temp[argsort], xticks=labels, yticks=ucores[argsort], gridstridex=4, gridstridey=None, vmin=-0.25, vmax=1.25, vline=numshifts*4-0.5, cmap="RdBu", ax=ax1)
    ax2.axis("off")
    plt.savefig(rndir + "_edges_left_clustered_" + m + ".png", bbox_inches="tight", dpi=600)

starttrim = numshifts
labels = posLabels[starttrim*4:]
temp = centered[:,starttrim*4:]

figsize = np.array(temp.transpose().shape)/4
figsize[0] /= 0.92
argsort = np.argsort([u[::-1] for u in ucores])
fig = plt.figure(figsize=figsize)
ax1 = fig.add_axes((0, 0, 1, 1))
tools.plotGrid(None, temp[argsort], xticks=labels, yticks=ucores[argsort], gridstridex=4, gridstridey=None, vmin=-0.25, vmax=1.25, vline=(numshifts-starttrim)*4-0.5, cmap="RdBu", ax=ax1)
plt.savefig(rndir + "_edges_right_lexi_r_sort.png", bbox_inches="tight", dpi=600)

for m in methods:
    linkage = shc.linkage(temp, method=m, metric="cityblock", optimal_ordering=True)
    fig = plt.figure(figsize=figsize)
    ax1 = fig.add_axes((0, 0, 1, 1))
    ax2 = fig.add_axes((-dendsize/figsize[0]-labelspace/figsize[0], 0, dendsize/figsize[0], 1))
    with plt.rc_context({'lines.linewidth': 1.0}):
        argsort = shc.dendrogram(linkage, orientation="left", ax=ax2, no_labels=True, color_threshold=0, above_threshold_color="black",)["leaves"][::-1]
    tools.plotGrid(None, temp[argsort], xticks=labels, yticks=ucores[argsort], gridstridex=4, gridstridey=None, vmin=-0.25, vmax=1.25, vline=(numshifts-starttrim)*4-0.5, cmap="RdBu", ax=ax1)
    ax1.yaxis.tick_left()
    ax2.axis("off")
    plt.savefig(rndir + "_edges_right_" + m + ".png", bbox_inches="tight", dpi=600)