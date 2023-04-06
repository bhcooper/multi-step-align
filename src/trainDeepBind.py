#! /usr/bin/env python

import tools
import keras
import keras_genomics

kernel_size = 8

import numpy as np
import pandas as pd
import sys
np.random.seed(1)
from sklearn import preprocessing
from sklearn import model_selection
from sklearn.metrics import r2_score

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

# Model parameters
kernel_size = 8

seqs_train, y_train = pd.read_csv(sys.argv[1], sep='\t', usecols=(0,1)).transpose().values
y_train = np.log(y_train.astype(float)).reshape(-1,1)
X_train = tools.encode1mers(seqs_train).reshape((len(seqs_train), -1, 4))
X_train = tools.encode1mers(seqs_train).reshape((len(seqs_train), -1, 4))

if(len(sys.argv) == 2):
    X_train, X_test, y_train, y_test = model_selection.train_test_split(X_train, y_train, test_size=0.3, random_state = 0)
else:
    seqs_test, y_test = pd.read_csv(sys.argv[2], sep='\t', usecols=(0,1)).transpose().values
    y_test = np.log(y_test.astype(float)).reshape(-1,1)
    X_test = tools.encode1mers(seqs_test).reshape((len(seqs_test), -1, 4))

def getDensityColors(x, y):
    bins = [1000, 1000]
    hh, locx, locy = np.histogram2d(x, y, bins=bins)
    c = np.array([hh[np.argmax(a<=locx[1:]),np.argmax(b<=locy[1:])] for a,b in zip(x,y)])
    c = preprocessing.scale(c)
    argsort = np.argsort(c)
    x = x[argsort]
    y = y[argsort]
    c = c[argsort]
    return c, x, y

def plotR2(filename, y1, y2, xlabel="", ylabel="", title="", c=[], cmap="RdYlBu", vmin=None, vmax=None, high = None, low = None, s = 3):
    fig, ax = plt.subplots(figsize=(3,3))
    if(len(c) == 0):
        c, y1, y2 = getDensityColors(y1, y2)
        cmap = "jet"
    ax.scatter(y1, y2, s=s, c = c, cmap=cmap, vmin=vmin, vmax=vmax)
    if(low == None):
        low = np.min(np.append(y1,y2))
    if(high == None):
        high = np.max(np.append(y1,y2))
    length = high-low
    low -= 0.05 * length
    high += 0.05 * length
    r2 = r2_score(y1, y2)
    print("R2: " + str(r2))
    ax.text(high-0.45*length,low+0.05*length, '${R^2}$ = %.2f' % r2, fontsize=12)
    plt.ylim((low,high))
    ax.set_xticks(ax.get_yticks())
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.xlim((low,high))
    plt.title(title, fontsize=18)
    plt.xlabel(xlabel, fontsize=14)
    plt.ylabel(ylabel, fontsize=14)
    plt.plot([low, high], [low, high], color = 'firebrick', alpha=0.5, linewidth = 1)
    plt.savefig(filename, bbox_inches="tight", dpi=1200)
    plt.close()
    return r2


model = keras.models.Sequential()
model.add(keras_genomics.layers.convolutional.RevCompConv1D(input_shape=(X_train.shape[1],4),
    filters=10,
    kernel_size=kernel_size,
    activation='relu',
    padding='same'))
model.add(keras_genomics.layers.normalization.RevCompConv1DBatchNorm())
model.add(keras_genomics.layers.convolutional.RevCompConv1D(filters=10,
    kernel_size=kernel_size,
    activation='relu',
    padding='same'))
model.add(keras_genomics.layers.normalization.RevCompConv1DBatchNorm())
model.add(keras_genomics.layers.convolutional.RevCompConv1D(filters=10,
    kernel_size=kernel_size,
    activation='relu',
    padding='same'))
model.add(keras_genomics.layers.normalization.RevCompConv1DBatchNorm())
model.add(keras.layers.MaxPooling1D(pool_size=10, padding='same'))
model.add(keras_genomics.layers.core.DenseAfterRevcompConv1D(units=10, activation='relu'))
model.add(keras.layers.Dense(units=1, activation='linear'))

model.compile(optimizer="adam", loss="MSE")

history = model.fit(X_train, y_train, epochs=10, verbose=1)

model.save("trained_cnn.h5")

# model = keras.models.load_model("trained_cnn.h5", custom_objects={
#         "RevCompConv1D":keras_genomics.layers.convolutional.RevCompConv1D,
#         "RevCompConv1DBatchNorm":keras_genomics.layers.normalization.RevCompConv1DBatchNorm,
#         "DenseAfterRevcompConv1D":keras_genomics.layers.core.DenseAfterRevcompConv1D})

pred_test = np.concatenate([model.predict(X_test[:,i:i+X_train.shape[1],:]) for i in range(X_test.shape[1] - X_train.shape[1] + 1)], axis=1)
pred_test = np.max(pred_test, axis=1)

plotR2("DeepBind_R2_train.png", y_train.ravel(), model.predict(X_train).ravel(), 'observed', 'predicted')
# plotR2("R2_train.svg", y_train.ravel(), model.predict(X_train).ravel(), 'observed', 'predicted')
plotR2("DeepBind_R2_test.png", y_test.ravel(), pred_test, 'observed', 'predicted')
# plotR2("R2_test.svg", y_test.ravel(), pred_test, 'observed', 'predicted')
