#!/usr/bin/env python

import sys
import numpy as np
import pandas as pd
import tools

from sklearn.linear_model import LassoCV
from sklearn.model_selection import cross_validate
from sklearn.preprocessing import StandardScaler
from sklearn.utils import shuffle
from sklearn.model_selection import train_test_split

import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.use('Agg')

from matplotlib import rcParams
rcParams['font.family'] = 'Arial'
rcParams['xtick.minor.size'] = 0
rcParams['ytick.minor.size'] = 0

N = ["A", "C", "G", "T"]

infile = sys.argv[1]
lflank = int(sys.argv[2])
rflank = int(sys.argv[3])

seqs, y = tools.readScores(infile)
argsort = np.argsort(-y)
seqs = seqs[argsort]
y = y[argsort]
seqs, y = shuffle(seqs, y, random_state=0)
y = np.log(y)

# Encode 1-hot cores
cores = None
if(rflank > 0):
    cores = np.array([seq[lflank:-rflank] for seq in seqs])
else:
    cores = np.array([seq[lflank:] for seq in seqs])
Xcore, ucores = tools.oneHotEncode(cores)

# Encode 1-mers of flanks
edges = None
if(rflank > 0):
    edges = np.array([seq[:lflank] + seq[-rflank:] for seq in seqs])
else:
    edges = np.array([seq[:lflank] for seq in seqs])
Xedges = tools.encode1mers(edges)

# 1-mer + 2-mer edges
Xedges_2mer = np.array([np.matmul(X.reshape(-1,1), X.reshape(1,-1)) for X in Xedges])
triu = np.triu_indices(Xedges_2mer.shape[1], k=1) 
Xedges_2mer = np.array([X[triu] for X in Xedges_2mer])

# Core dependent edges
Xedges_core = np.array([np.matmul(xc.reshape(-1,1), xe.reshape(1,-1)).flatten() for xc, xe in zip(Xcore, Xedges)])

# Full 1-mer encoding
X_1mer = tools.encode1mers(seqs)

# Feature concatenation
Xs = [
    np.concatenate((Xcore, Xedges), axis=1),
    np.concatenate((Xcore, Xedges_2mer), axis=1),
    np.concatenate((Xcore, Xedges_core), axis=1),
    np.concatenate((Xcore, Xedges_core, Xedges_2mer), axis=1)]

featTable = pd.DataFrame({
    "Core Identity": [True, True, True, True],
    "Flank Mononucleotides": [True, True, True, True],
    "Flank Dinucleotides": [False, True, False, True],
    "Core-Dependent Mononucleotides": [False, False, True, True]})

scores = []
for i,X in enumerate(Xs):
    print("\nTraining model . . .")
    print(featTable.iloc[i,:])
    X = StandardScaler().fit_transform(X)
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3)
    model = LassoCV(max_iter = int(1e6), cv=5, n_jobs=-1)
    model = cross_validate(model, X, y, cv=5, scoring="r2", return_train_score=True)
    print('Train R²: %.3f' % np.mean(model["train_score"]))
    print('Test  R²: %.3f' % np.mean(model["test_score"]))
    scores += [np.mean(model["test_score"])]

fig, ax = plt.subplots(2, 1, figsize=(3,4), gridspec_kw={'height_ratios': [3, 1]}, sharex=True)
plt.subplots_adjust(hspace=0.2)

ax[0].bar(range(len(featTable)), scores)
ax[0].set_xticks([])
ax[0].set_ylim([0.0,1.0])
ax[0].set_ylabel('Model Performance [${R^2}$]', fontsize=10)
ax[1].matshow(featTable.values.transpose(), cmap="Blues", aspect='auto', alpha=0.5)
ax[1].set_xticks([])
ax[1].set_yticks(range(featTable.shape[1]))
ax[1].set_yticklabels([s + ":" for s in featTable.columns])
ax[1].spines[:].set_visible(False)
ax[1].set_xticks(np.arange(len(featTable)+1)-0.5, minor=True)
ax[1].set_yticks(np.arange(featTable.shape[1]+1)-0.5, minor=True)
ax[1].grid(which="minor", axis='x', color="white", linewidth=8)
ax[1].grid(which="minor", axis='y', color="white", linewidth=4)
ax[1].set_title("Model Features", fontsize=10, y=0.95)
plt.tick_params(left=False)

fig.suptitle(f"MLR Performance for Core-aligned {len(seqs[0])}-mers", y=0.95)

plt.savefig(infile[:-4] + "_MLR.png", bbox_inches='tight', dpi=600)
# plt.savefig(sys.argv[1][:-4] + "_MLR.svg", bbox_inches='tight', dpi=600)
plt.close()

featTable['Mean Test R^2'] = scores
featTable.to_csv(infile[:-4] + "_MLR.tsv", sep='\t', index=False)
