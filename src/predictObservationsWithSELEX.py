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

imgext = '.png'
# imgext = '.svg'

def MSE(x, y):
    return np.sum((x - y)**2)/len(x)

Nlookup = {"A":0, "C":1, "G":2, "T":3}

comp = {"A":"T", "C":"G", "G":"C", "T":"A", "N":"N"}
def rc(seqs):
    return np.array(["".join([comp[x] for x in seq][::-1]) for seq in seqs])

obs = pd.read_csv(sys.argv[1], sep='\t', index_col=0, usecols=(0,1), names=['', 'p'], skiprows=1)
obs = obs/obs.sum()

bg = pd.read_csv(sys.argv[2], sep='\t', index_col=0, usecols=(0,1), names=['', 'p'], skiprows=1)
bg /= bg.sum()

relE = pd.read_csv(sys.argv[3], sep='\t', index_col=0, usecols=(0,1), names=['', 'p'], skiprows=1)
relE /= relE.sum()

pred = bg * relE
pred /= pred.sum()
pred = pred.sort_values(pred.columns[0], ascending=False)

df = obs.merge(pred, left_index=True, right_index=True)
df = np.log(df)

fig, ax = plt.subplots(figsize=(3,3))
ax.scatter(df.iloc[:,0], df.iloc[:,1], s=3)
low = df.min().min()
high = df.max().max()
length = high-low
low -= 0.05 * length
high += 0.05 * length
r = np.corrcoef(df.iloc[:,0], df.iloc[:,1])[0,1]
print("r: " + str(r))
mse = MSE(df.iloc[:,0], df.iloc[:,1])
print("MSE: " + str(mse))
ax.text(low+(0.3*length),high-(0.15*length), '${r}$ = %.2f' % r, fontsize=12, ha='center')
ax.text(low+(0.3*length),high-(0.23*length), '${MSE}$ = %.2f' % mse, fontsize=12, ha='center')
plt.ylim((low,high))
ax.set_xticks(ax.get_yticks())
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.xlim((low,high))
plt.title('ChIP-exo Core Preferences', fontsize=14)
plt.xlabel('observed [ln($p$)]', fontsize=14)
plt.ylabel('alignment-based [ln($p$)]', fontsize=14)
# plt.plot([low, high], [low, high], color = 'firebrick', alpha=0.5, linewidth = 1)
plt.savefig("predSELEX" + imgext, bbox_inches="tight", dpi=1200)
plt.close()