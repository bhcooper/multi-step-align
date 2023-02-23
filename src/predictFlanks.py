#!/usr/bin/env python

import sys
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.use('Agg')

from matplotlib import rcParams
rcParams['font.family'] = 'Arial'
rcParams['xtick.major.size'] = 0
rcParams['xtick.minor.size'] = 0
rcParams['ytick.major.size'] = 0
rcParams['ytick.minor.size'] = 0
Nlookup = {"A":0, "C":1, "G":2, "T":3}

def MSE(x, y):
    return np.sum((x - y)**2)

N = ['A', 'C', 'G', 'T']

obs = pd.read_csv(sys.argv[1], sep='\t', index_col=0).transpose()
obs /= obs.sum()

bg = pd.read_csv(sys.argv[2], sep='\t', index_col=0).transpose()
bg /= bg.sum()

core = sys.argv[4]
flankddG = pd.read_csv(sys.argv[3], sep='\t', index_col=0).loc[core]
flankddG = pd.DataFrame(flankddG.values.reshape(-1, 4), columns=N)
flankddG.index = list(range(int(-len(flankddG)/2), 0)) + list(range(1, int(len(flankddG)/2+1)))
flankddG = flankddG.loc[obs.columns].transpose()
selexE = np.exp(-flankddG)
selexE /= selexE.sum()

selexPred = bg * selexE
selexPred /= selexPred.sum()

beesemE = pd.read_csv(sys.argv[5], sep='\t', skiprows=1).transpose()
beesemE.columns = obs.columns
beesemE /= beesemE.sum()

beesemPred = bg * beesemE
beesemPred /= beesemPred.sum()

obs = np.log(obs)
bd = np.log(bg)
selexPred = np.log(selexPred)
beesemPred = np.log(beesemPred)

x = np.corrcoef(obs.values.flatten(), selexPred.values.flatten())[0,1]
print(f'SELEX r =  {x}')
x = np.corrcoef(obs.values.flatten(), beesemPred.values.flatten())[0,1]
print(f'BEESEM r = {x}')
x = spearmanr(obs.values.flatten(), selexPred.values.flatten())[0]
print(f'SELEX ρ =  {x}')
x = spearmanr(obs.values.flatten(), beesemPred.values.flatten())[0]
print(f'BEESEM ρ = {x}')
x = MSE(obs.values.flatten(), selexPred.values.flatten())
print(f'SELEX MSE =  {x}')
x = MSE(obs.values.flatten(), beesemPred.values.flatten())
print(f'BEESEM MSE = {x}')

matrix = np.array([obs.values.transpose().flatten(), selexPred.values.transpose().flatten(), beesemPred.values.transpose().flatten()])
vmin = -np.ceil(abs(np.min(matrix)*10))/10
vmax = -np.floor(abs(np.max(matrix)*10))/10
print(f'vmin: {vmin}')
print(f'vmax: {vmax}')

lFlankLen = np.sum(obs.columns < 0)
rFlankLen = np.sum(obs.columns > 0)

figsize = np.array(matrix.transpose().shape)/4
fig, ax = plt.subplots(figsize=figsize)
mat = ax.matshow(matrix, cmap='RdBu', aspect="auto", vmin = vmin, vmax = vmax)
ax.axvline(lFlankLen * 4 - 0.5, color='red', linewidth=2)

posLabels = ["−" + str(s) + "\n" + b if b == "C" else "\n" + b for b, s in zip(N*lFlankLen, np.repeat(np.arange(lFlankLen)[::-1] + 1, 4))]
posLabels += ["+" + str(s) + "\n" + b if b == "C" else "\n" + b for b, s in zip(N*rFlankLen, np.repeat(np.arange(rFlankLen) + 1, 4))]

ax.set_xticks(range(matrix.shape[1]))
ax.set_yticks(range(matrix.shape[0]))
ax.set_xticks(np.arange(matrix.shape[1]+1)[::4] - 0.5, minor=True)
ax.set_yticks(np.arange(matrix.shape[0]+1)[::len(matrix)] - 0.5, minor=True)
ax.set_xticklabels(posLabels, fontsize=14, rotation=0)
ax.set_yticklabels(['observed', 'alignment-based', 'BEESEM-based'], fontsize=14, rotation=0, fontname="Arial")
for edge, spine in ax.spines.items():
    spine.set_visible(False)
ax.grid(which="minor", color="black", linewidth=0.5)
ax.yaxis.set_label_position('right')
ax.set_ylabel('[ln($p$)]', fontsize=12, labelpad = 7)
ax.set_title(f"ChIP-exo {core} Flanking Preferences", fontsize=16, y=2)
cbar = plt.colorbar(mat, aspect=5, ticks = [vmin, vmax], pad=0.07)
cbar.outline.set_linewidth(0.3)
plt.savefig('predFlanks_' + core + ".png", bbox_inches="tight", pad_inches=0, dpi=600)
plt.close()
