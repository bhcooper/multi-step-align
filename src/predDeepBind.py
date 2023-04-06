#! /usr/bin/env python

#install keras from https://github.com/kundajelab/keras/tree/keras_1
import keras
import keras_genomics
import numpy as np
import pandas as pd
import sys
import tools
np.random.seed(1)
from sklearn.linear_model import LassoCV

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import logomaker as lm
import argparse

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

def plotGrid(filename, matrix, xticks=[], yticks=[], xtickrotation="horizontal", ytickrotation="horizontal", xlabel="", ylabel="", title="", cmap="RdYlBu", bounds=None, gridstridex=None, gridstridey = None, figsize=None, vline=None, ax=None, tick_left = None):

    figsize = np.array(matrix.transpose().shape)/4
    fig = plt.figure(figsize=figsize)
    ax1 = fig.add_axes((0, 0, 1, 1))
    ax2 = fig.add_axes((1 + 1.4/figsize[0], 0.6/figsize[1], 0.3/figsize[0], 2.3/figsize[1]))
    
    # ax1.pcolormesh(matrix, cmap=cmap, vmin=vmin, vmax=vmax)
    ax1.invert_yaxis()
    ax1.matshow(matrix, cmap=cmap, vmin=-bounds, vmax=bounds, aspect="auto")
    ax1.axvline(vline+0.5, color='red', linewidth=2)
    # ax.set_xticks(np.arange(matrix.shape[1])+0.5)
    ax1.set_xticks(np.arange(matrix.shape[1]))
    # ax.set_yticks(np.arange(matrix.shape[0])+0.5)
    ax1.set_yticks(np.arange(matrix.shape[0]))
    if(gridstridex == None):
        gridstridex = matrix.shape[1]
    if(gridstridey == None):
        gridstridey = matrix.shape[0]
    # ax.set_xticks(np.append(np.arange(matrix.shape[1]+1)[::gridstridex], matrix.shape[1]), minor=True)
    ax1.set_xticks(np.arange(matrix.shape[1]+1)[::gridstridex] - 0.5, minor=True)
    # ax.set_yticks(np.append(np.arange(matrix.shape[0]+1)[::gridstridey], matrix.shape[0]), minor=True)
    ax1.set_yticks(np.arange(matrix.shape[0]+1)[::gridstridey] - 0.5, minor=True)
    ax1.set_xticklabels(xticks, fontsize=16, rotation=xtickrotation)
    ax1.xaxis.set_ticks_position('top') 
    ax1.set_yticklabels(yticks, fontsize=18, rotation=ytickrotation, fontname="Courier New")
    ax1.yaxis.tick_right()
    for edge, spine in ax1.spines.items():
        spine.set_visible(False)
    ax1.grid(which="minor", color="black", linewidth=0.5)
    ax1.set_xlabel(xlabel, fontsize=18, fontweight="bold")
    ax1.xaxis.set_label_position('top') 
    ax1.set_ylabel(ylabel, fontsize=18, fontweight="bold")
    ax1.set_title(title)
    if(tick_left):
        ax1.yaxis.tick_left()

    plt.colorbar(mpl.cm.ScalarMappable(cmap=cmap), cax=ax2, pad=0, ticks=[0, 1])
    ax2.set_yticklabels([-bounds, bounds], fontsize=18)
    ax2.set_ylabel("$\it{-ΔΔG/RT}$", fontsize=18, rotation=270, labelpad=2, fontweight='bold')

    plt.savefig(filename, bbox_inches="tight", pad_inches=0, dpi=600)
    plt.close()

parser = argparse.ArgumentParser(description="Use MLR to interpret flanking preferences according to a trained DeepBind model.")
parser.add_argument("model_file", help="Trained model file in .h5 format")
parser.add_argument("cores_table", help="Tab-delimited file containing a list of cores in the first column")
parser.add_argument("flanks_table", help="Tab-delimited file containing cores and their corresponding flanking preferences.")
parser.add_argument("left_flank_len", type=int, help="Number of bp to consider 5' of the core.")
parser.add_argument("right_flank_len", type=int, help="Number of bp to consider 3' of the core.")
parser.add_argument("-B", "--bounds", default=None, type=float, required=False, help='Specifies the color bar bounds ± the mean for plotting the ΔΔG/RT of flanking positions. Helpful for consistency accross comparisons.')
args = parser.parse_args()


modelfile = args.model_file
corefile = args.cores_table
flankfile = args.flanks_table
lFlankLen = args.left_flank_len
rFlankLen = args.right_flank_len
bounds = args.bounds

cores = pd.read_csv(corefile, sep='\t', usecols=(0,)).values.ravel()
flanks = -pd.read_csv(flankfile, sep='\t', index_col=0)

lFlank = permute(lFlankLen)
rFlank = permute(rFlankLen)
nperm = len(lFlank) * len(rFlank)

corelen = len(cores[0])
matrix = np.zeros((len(cores), (lFlankLen+rFlankLen)*4))
edges = np.zeros((len(cores), (lFlankLen+rFlankLen)*4))

model = keras.models.load_model(modelfile, custom_objects={
        "RevCompConv1D":keras_genomics.layers.convolutional.RevCompConv1D,
        "RevCompConv1DBatchNorm":keras_genomics.layers.normalization.RevCompConv1DBatchNorm,
        "DenseAfterRevcompConv1D":keras_genomics.layers.core.DenseAfterRevcompConv1D})


# pip install pydot graphviz
from keras.utils import plot_model
plot_model(model, to_file='model.svg', show_shapes=True, show_layer_names=False, expand_nested=False, dpi=600)

model.summary()
tools.mkdir("DeepBind_logos")

for i,core in enumerate(cores):
    Xcore = np.array([encode[c] for c in core]).flatten()
    X = np.concatenate((np.tile(lFlank, (len(rFlank), 1)), np.broadcast_to(Xcore, (nperm, len(Xcore))), np.repeat(rFlank, len(lFlank), axis=0)), axis=1)
    X = X.reshape(len(X), -1, 4)

    pred = model.predict(X).ravel()

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
    plt.savefig(f"DeepBind_logos/DeepBind_" + core + ".png", bbox_inches="tight", dpi=600)
    plt.close()

if(not bounds):
    bounds = np.ceil((max(np.nanmax(np.abs(matrix)), np.nanmax(np.abs(edges)))*100))/100

posLabels = ["−" + str(s) + "\n" + b if b == "C" else "\n" + b for b, s in zip(N*4, np.repeat(np.arange(4)[::-1] + 1, 4))]
posLabels += ["+" + str(s) + "\n" + b if b == "C" else "\n" + b for b, s in zip(N*2, np.repeat(np.arange(2) + 1, 4))]

plotGrid("DeepBind_flank_matrix.png", matrix, bounds=bounds, cmap='RdBu', gridstridex = 4, xticks = posLabels, yticks=cores, vline=4*4)
plotGrid("alignment_flank_matrix.png", edges, bounds=bounds, cmap='RdBu', gridstridex = 4, xticks = posLabels, yticks=cores, vline=4*4)

