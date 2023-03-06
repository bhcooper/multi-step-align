#!python

#install keras from https://github.com/kundajelab/keras/tree/keras_1
import keras
import keras_genomics
import numpy as np
import pandas as pd
import sys
import pickle
np.random.seed(1)
from scipy.stats import spearmanr
from sklearn.linear_model import LassoCV

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import logomaker as lm

from matplotlib import rcParams
rcParams['font.family'] = 'Arial'
rcParams['xtick.major.size'] = 0
rcParams['xtick.minor.size'] = 0
rcParams['ytick.major.size'] = 0
rcParams['ytick.minor.size'] = 0

encode = {"A":[1,0,0,0], "C":[0,1,0,0], "G":[0,0,1,0], "T":[0,0,0,1], "N":[0.25,0.25,0.25,0.25]}

N = ['A', 'C', 'G', 'T']
NX = np.eye(4)
def permute(plen):
    pX = NX.copy()
    for i in range(1, plen):
        pX = np.concatenate([np.concatenate((pX, np.broadcast_to(row, (len(pX), 4))), axis=1) for row in NX], axis=0)
    return pX

def plotGrid(filename, matrix, xticks=[], yticks=[], xtickrotation="horizontal", ytickrotation="horizontal", xlabel="", ylabel="", title="", cmap="RdYlBu", vmin=None, vmax=None, gridstridex=None, gridstridey = None, figsize=None, vline=None, ax=None, vlines = None, tick_left = None):
    axGiven = True
    if(figsize == None):
        if(len(matrix.shape) == 1):
            matrix = matrix.reshape(1, -1)
        figsize = np.array(matrix.transpose().shape)/4
    if(ax == None):
        fig, ax = plt.subplots(figsize=figsize)
        axGiven = False
    ax.pcolormesh(matrix, cmap=cmap, vmin=vmin, vmax=vmax)
    ax.invert_yaxis()
    # ax.matshow(matrix, cmap=cmap, vmin=vmin, vmax=vmax, aspect="auto")
    if(vline):
        # ax.axvline(vline+0.5, color='red', linewidth=2)
        ax.axvline(vline, color='red', linewidth=2)
    if(vlines):
        [ax.axvline(line+0.5, color='firebrick', linewidth=4) for line in vlines]
        # [ax.axvline(line, color='firebrick', linewidth=4) for line in vlines]
    ax.set_xticks(np.arange(matrix.shape[1])+0.5)
    # ax.set_xticks(np.arange(matrix.shape[1]))
    ax.set_yticks(np.arange(matrix.shape[0])+0.5)
    # ax.set_yticks(np.arange(matrix.shape[0]))
    if(gridstridex == None):
        gridstridex = matrix.shape[1]
    if(gridstridey == None):
        gridstridey = matrix.shape[0]
    ax.set_xticks(np.append(np.arange(matrix.shape[1]+1)[::gridstridex], matrix.shape[1]), minor=True)
    # ax.set_xticks(np.arange(matrix.shape[1]+1)[::gridstridex] - 0.5, minor=True)
    ax.set_yticks(np.append(np.arange(matrix.shape[0]+1)[::gridstridey], matrix.shape[0]), minor=True)
    # ax.set_yticks(np.arange(matrix.shape[0]+1)[::gridstridey] - 0.5, minor=True)
    ax.set_xticklabels(xticks, fontsize=16, rotation=xtickrotation)
    ax.xaxis.set_ticks_position('top') 
    ax.set_yticklabels(yticks, fontsize=18, rotation=ytickrotation, fontname="Courier New")
    ax.yaxis.tick_right()
    for edge, spine in ax.spines.items():
        spine.set_visible(False)
    ax.grid(which="minor", color="black", linewidth=0.5)
    ax.set_xlabel(xlabel, fontsize=18, fontweight="bold")
    ax.xaxis.set_label_position('top') 
    ax.set_ylabel(ylabel, fontsize=18, fontweight="bold")
    ax.set_title(title)
    if(tick_left):
        ax.yaxis.tick_left()

    if(not axGiven):
        plt.savefig(filename, bbox_inches="tight", pad_inches=0, dpi=600)
        plt.close()


lFlankLen = 4
rFlankLen = 2
modelPadLen = 4

lFlank = permute(lFlankLen)
rFlank = permute(rFlankLen)
nperm = len(lFlank) * len(rFlank)
modelPad = np.full(modelPadLen*4, 0.25)
lFlank = np.concatenate((np.broadcast_to(modelPad, (len(lFlank), len(modelPad))), lFlank), axis=1)
rFlank = np.concatenate((rFlank, np.broadcast_to(modelPad, (len(rFlank), len(modelPad)))), axis=1)

cores = pd.read_csv(sys.argv[1], sep='\t', usecols=(0,)).values.ravel()
flanks = -pd.read_csv(sys.argv[2], sep='\t', index_col=0)
corelen = len(cores[0])
matrix = np.zeros((len(cores), (lFlankLen+rFlankLen)*4))
edges = np.zeros((len(cores), (lFlankLen+rFlankLen)*4))

model = keras.models.load_model("trained_cnn.h5", custom_objects={
        "RevCompConv1D":keras_genomics.layers.convolutional.RevCompConv1D,
        "RevCompConv1DBatchNorm":keras_genomics.layers.normalization.RevCompConv1DBatchNorm,
        "DenseAfterRevcompConv1D":keras_genomics.layers.core.DenseAfterRevcompConv1D})

# modelMin = modelMax = None
# with open("trained_minMax.pkl", 'rb') as f:
    # modelMin, modelMax = pickle.load(f)
    
for i,core in enumerate(cores):
    # No reason to encode code, can remove
    Xcore = np.array([encode[c] for c in core]).flatten()
    X = np.concatenate((np.tile(lFlank, (len(rFlank), 1)), np.broadcast_to(Xcore, (nperm, len(Xcore))), np.repeat(rFlank, len(lFlank), axis=0)), axis=1)
    X = X.reshape(len(X), -1, 4)
    # pred_test = np.concatenate([model.predict(X_test[:,i:i+X_train.shape[1],:]) for i in range(X_test.shape[1] - X_train.shape[1] + 1)], axis=1)
    # pred_test = np.max(pred_test, axis=1)

    pred = model.predict(X).ravel()
    # pred = pred * (modelMax - modelMin) + modelMin
    # df = pd.DataFrame({"seq":seqs, "pred":pred}).sort_values("pred", ascending=False)
    # df.to_csv(sys.argv[1][:-4] + "_DeepBind.tsv", sep='\t')

    X = X[:,modelPadLen:-modelPadLen,:]
    X = np.delete(X, np.arange(X.shape[1])[lFlankLen:-rFlankLen], axis=1)
    X = X.reshape((len(X), -1))
    MLR = LassoCV().fit(X, pred)
    print(core)
    print("R² = " + str(MLR.score(X, pred)))
    coef = MLR.coef_.reshape(-1,4)
    coef -= np.mean(coef, axis=1).reshape(-1,1)
    flank = flanks.loc[core].values.reshape(-1,4)[5:-7]
    flank -= np.mean(flank, axis=1).reshape(-1,1)
    matrix[i] = coef.flatten()        
    edges[i] = flank.flatten()        
    
    coef = pd.DataFrame(coef, columns=['A','C','G','T'])
    flank = pd.DataFrame(flank, columns=['A','C','G','T'])

    fig, ax = plt.subplots(2, 1, figsize=(0.55 * len(coef), 4.0))
    logo1 = lm.Logo(coef, font_name="Roboto", ax=ax[0], fade_below=0.5, 
        color_scheme={
            'A': (16/255, 150/255, 72/255),
            'C':(37/255, 92/255, 153/255),
            'G':(247/255, 179/255, 43/255),
            'T':(214/255, 40/255, 57/255)})
    logo1.ax.set_ylabel("MLR-based\n$−ΔΔG/RT$", fontsize=15)
    logo1.ax.set_xticks(np.arange(6))
    logo1.ax.set_xticklabels(['−4', '−3', '−2', '−1', '+1', '+2'])
    logo1.ax.axvline(3.5, color='red', linewidth=2)
    logo1.style_glyphs_below(flip=False)
    logo1.ax.set_ylim([-0.75, 0.75])
    logo2 = lm.Logo(flank, font_name="Roboto", ax=ax[1], fade_below=0.5, 
        color_scheme={
            'A': (16/255, 150/255, 72/255),
            'C':(37/255, 92/255, 153/255),
            'G':(247/255, 179/255, 43/255),
            'T':(214/255, 40/255, 57/255)})
    logo2.ax.set_ylabel("alignment-based\n$−ΔΔG/RT$", fontsize=15)
    logo2.ax.set_xticks(np.arange(6))
    logo2.ax.set_xticklabels(['−4', '−3', '−2', '−1', '+1', '+2'])
    logo2.ax.axvline(3.5, color='red', linewidth=2)
    logo2.style_glyphs_below(flip=False)
    logo2.ax.set_ylim([-0.75, 0.75])
    plt.suptitle(core)
    plt.savefig(f"DeepBind_" + core + ".png", bbox_inches="tight", dpi=600)
    plt.savefig(f"DeepBind_" + core + ".svg", bbox_inches="tight", dpi=600)
    plt.close()

vmax = np.max(np.abs(matrix))
print(vmax)
vmax = np.max(np.abs(edges))
print(vmax)
vmax = 0.75
print(vmax)

posLabels = ["−" + str(s) + "\n" + b if b == "C" else "\n" + b for b, s in zip(N*4, np.repeat(np.arange(4)[::-1] + 1, 4))]
posLabels += ["+" + str(s) + "\n" + b if b == "C" else "\n" + b for b, s in zip(N*2, np.repeat(np.arange(2) + 1, 4))]

plotGrid("DeepBind_matrix.png", matrix, vmin=-vmax, vmax=vmax, cmap='RdBu', gridstridex = 4, xticks = posLabels, yticks=cores, vline=4*4)
plotGrid("DeepBind_matrix.svg", matrix, vmin=-vmax, vmax=vmax, cmap='RdBu', gridstridex = 4, xticks = posLabels, yticks=cores, vline=4*4)

plotGrid("Edge_matrix.png", edges, vmin=-vmax, vmax=vmax, cmap='RdBu', gridstridex = 4, xticks = posLabels, yticks=cores, vline=4*4)
plotGrid("Edge_matrix.svg", edges, vmin=-vmax, vmax=vmax, cmap='RdBu', gridstridex = 4, xticks = posLabels, yticks=cores, vline=4*4)

