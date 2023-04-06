#!/usr/bin/env python

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.use('Agg')
import multiprocessing as mp
import sharedmem as sm
import yaml
import pickle
import regex as re

from sklearn import metrics
from sklearn import preprocessing
from collections import defaultdict
from pandas import read_csv
from pandas import DataFrame as df

from matplotlib import rcParams
rcParams['font.family'] = 'Arial'
rcParams['xtick.major.size'] = 0
rcParams['xtick.minor.size'] = 0
rcParams['ytick.major.size'] = 0
rcParams['ytick.minor.size'] = 0
rcParams['ytick.direction'] = 'in'
rcParams['xtick.direction'] = 'in'
rcParams['svg.fonttype'] = 'none'

N = ["A", "C", "G", "T"]
IUPAC = {"A":["A"], "C":["C"], "G":["G"], "T":["T"], "R":["A", "G"],"Y":["C", "T"],"K":["G", "T"],"M":["A", "C"],"S":["C", "G"],
    "W":["A", "T"],"B":["C", "G", "T"],"D":["A", "G", "T"],"H":["A", "C", "T"],
    "V":["A", "C", "G"],"N":["A", "C", "G", "T"]}

def loadConfig(fname):
    with open(fname, 'r') as cfile:
        return yaml.safe_load(cfile)

def saveVariables(pkl, tup):
    with open(pkl, 'wb') as f:
        pickle.dump(tup, f)

def loadVariables(pkl):
    tup = None
    with open(pkl, 'rb') as f:
        tup = pickle.load(f)
    return tup

def saveTable(filename, data, header="Seqs\ty"):
    np.savetxt(filename, data, fmt="%s", delimiter="\t", header=header, comments="")

base = {"A":0, "C":1, "G":2, "T":3}
def calcOffset(kmer):
    order = len(kmer)
    offset = 0
    for i, c in enumerate(kmer):
        offset += base[c] * 4 ** (order-1)
        order -= 1
    return offset

comp = {"A":"T", "C":"G", "G":"C", "T":"A", "N":"N"}
def rc(seqs):
    return np.array(["".join([comp[x] for x in seq][::-1]) for seq in seqs])

def getPermutations(seq):
    seqs = seq.upper()
    permutations = IUPAC[seqs[0]]
    for c in seqs[1:]:
        permutations = [seq + p for p in IUPAC[c] for seq in permutations]
    return np.array(permutations)

def mkdir(path):
    if not os.path.exists(path):
        os.mkdir(path)

def readScores(file, nrows=None):
    data = None
    if(nrows == None):
        data = read_csv(file, sep="\t", usecols=(0,1))
    else:
        data = read_csv(file, sep="\t", usecols=(0,1), nrows=nrows)
    seqs = data.iloc[:,0].values
    y = data.iloc[:,1].values   
    return seqs, y

def getOverlapSet(seqslist):
    numlists = len(seqslist)
    seqs, counts = np.unique(np.concatenate(seqslist), return_counts=True)
    return seqs[counts == numlists]

def plotGrid(filename, matrix, xticks=[], yticks=[], xtickrotation="horizontal", ytickrotation="horizontal", xlabel="", ylabel="", title="", cmap="RdYlBu", vmin=None, vmax=None, gridstridex=None, gridstridey = None, figsize=None, vline=None, ax=None, vlines = None, tick_left = None):
    axGiven = True
    if(figsize == None):
        if(len(matrix.shape) == 1):
            matrix = matrix.reshape(1, -1)
        figsize = np.array(matrix.transpose().shape)/4
    if(ax == None):
        fig, ax = plt.subplots(figsize=figsize)
        axGiven = False
    ax.matshow(matrix, cmap=cmap, vmin=vmin, vmax=vmax, aspect='equal')
    ax.invert_yaxis()
    if(vline):
        ax.axvline(vline, color='red', linewidth=2)
    if(vlines):
        [ax.axvline(line, color='firebrick', linewidth=4) for line in vlines]
    ax.set_xticks(np.arange(matrix.shape[1]))
    # ax.set_xticks(np.arange(matrix.shape[1]))
    ax.set_yticks(np.arange(matrix.shape[0]))
    # ax.set_yticks(np.arange(matrix.shape[0]))
    if(gridstridex == None):
        gridstridex = matrix.shape[1]
    if(gridstridey == None):
        gridstridey = matrix.shape[0]
    # ax.set_xticks(np.append(np.arange(matrix.shape[1]+1)[::gridstridex], matrix.shape[1]), minor=True)
    ax.set_xticks(np.arange(matrix.shape[1]+1)[::gridstridex] - 0.5, minor=True)
    # ax.set_yticks(np.append(np.arange(matrix.shape[0]+1)[::gridstridey], matrix.shape[0]), minor=True)
    ax.set_yticks((np.arange(matrix.shape[0]+1)[::gridstridey] - 0.5), minor=True)
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

def consolidate(seqs, enrichment):
    lookup = defaultdict(float, zip(seqs, enrichment))
    for key in list(lookup.keys()):
        score = lookup[key]
        rcScore = lookup[rc([key])[0]]
        if(not rcScore == 0.0 and not rcScore == score):
            ave = (score + rcScore) / 2
            lookup[key] = ave
            lookup[rc([key])[0]] = ave
    return np.array([lookup[seq] for seq in seqs])

def oneHotEncode(Xlist):
    if(not type(Xlist) == list):
        Xlist = [Xlist]
    Xlist = [X.reshape(-1,1) if len(X.shape) == 1 else X for X in Xlist]
    encoder = preprocessing.OneHotEncoder(sparse=False).fit(np.concatenate(Xlist, axis=0))
    Xlist = tuple(encoder.transform(X) for X in Xlist)
    return Xlist + (encoder.categories_[0],)

def consolidatedf(df):
    loc = df.index
    df_rc = df.copy()
    df_rc.index = rc(df.index)
    df = pd.concat((df, df_rc))
    df = df.groupby(df.index).mean()
    return df.loc[loc]

def parallelAsync(target, ncpu, *args):
    if(ncpu == 0):
        ncpu = mp.cpu_count()
    pList = [mp.Process(target=target, args=list(args) + [indices]) for indices in np.array_split(np.arange(len(args[1])), ncpu)]
    [p.start() for p in pList]
    [p.join() for p in pList]
    return args[0]

def parallelSync(target, ncpu, shared, *args):
    if(ncpu == 0):
        ncpu = mp.cpu_count()
    temps = []
    for s in shared:
        temps += [sm.empty((ncpu,) + s.shape, dtype=s.dtype)]
    pList = [mp.Process(target=target, args=temps + list(args) + [indices,p]) for p, indices in enumerate(np.array_split(np.arange(len(args[0])), ncpu))]
    [p.start() for p in pList]
    [p.join() for p in pList]
    for s,t in zip(shared, temps):
        s[:] = np.sum(t, axis=0)
    return shared

def encode1mers_task(X, seqs, indices):    
    for i in indices:
        for j,c in enumerate(seqs[i]):
            if(not c == "N"): 
                X[i,j * 4 + calcOffset(c)] = True
    
def encode1mers(seqs, ncpu=0):
    if(ncpu == 0):
        ncpu = mp.cpu_count()
    X = sm.empty((len(seqs), len(seqs[0])*4), dtype=bool)
    parallelAsync(encode1mers_task, ncpu, X, seqs)
    return X

def getKmerDist_task(counts, seqs, seq_counts, seqlookup, k, nwindows, indices, pid):
    for i in range(nwindows):
        kmers = np.array([seqlookup[seq[i:i+k]] for seq in seqs[indices]])   
        for t,count in zip(kmers, seq_counts[indices]):
            counts[pid,t,i] += count

def getKmerDist(seqs, seq_counts, k, ncpu=0):
    if(ncpu == 0):
        ncpu = mp.cpu_count()
    useqs = getPermutations("N"*k)
    seqlookup = dict(zip(useqs, np.arange(len(useqs))))
    seqlen = len(seqs[0])
    nwindows = seqlen - k + 1
    counts = sm.empty((len(useqs), nwindows))
    parallelSync(getKmerDist_task, ncpu, (counts,), seqs, seq_counts, seqlookup, k, nwindows)
    minCount = np.min(counts)
    counts = counts / np.sum(seq_counts)
    counts = df(counts, index=useqs)
    return counts, minCount

def countKmers_task(seqs, seq_counts, k, nwindows, indices, q):
    counts = defaultdict(float)
    for i in range(nwindows):
        for seq,count in zip(seqs[indices], seq_counts[indices]):
            counts[seq[i:i+k]] += count
    q.put(counts)

def countKmers(seqs, seq_counts, k, ncpu=0):
    if(ncpu == 0):
        ncpu = mp.cpu_count()
    seqlen = len(seqs[0])
    nwindows = seqlen - k + 1
    q = mp.Queue()
    pList = [mp.Process(target=countKmers_task, args=[seqs, seq_counts, k, nwindows, indices, q]) for indices in np.array_split(np.arange(len(seqs)), ncpu)]
    [p.start() for p in pList]
    counts = defaultdict(float)
    for p in pList:
        for key, value in q.get().items():
            counts[key] += value
    [p.join() for p in pList]
    return counts

# Update core ID and start position if unique core (including gap) fully in variable region
def getHit(cores, starts, seqs, ladapter, radapter, query, ncores, corelookup, corelens, maxcorelen, indices):
    ladapterlen = len(ladapter)
    seqlen = len(seqs[0])
    for i in indices:
        matches = list(re.finditer(query, ladapter + seqs[i] + radapter, overlapped=True))
        if(len(matches) == 1):
            c = corelookup[matches[0][0]]
            ngap = maxcorelen - corelens[c]
            if ngap:
                if (c < ncores):
                    if(matches[0].start() >= ladapterlen and matches[0].end() <= ladapterlen + seqlen - ngap):
                        cores[i] = c
                        starts[i] = matches[0].start() - ladapterlen
                # onRC
                else: 
                    if(matches[0].start() >= ladapterlen + ngap and matches[0].end() <= ladapterlen + seqlen):
                        cores[i] = c
                        starts[i] = matches[0].start() - ladapterlen
            else:
                if (matches[0].start() >= ladapterlen and matches[0].end() <= ladapterlen + seqlen):
                    cores[i] = c
                    starts[i] = matches[0].start() - ladapterlen

def countHits(counts, seqs, ladapter, radapter, query, indices):
    for i in indices:
        matches = list(re.finditer(query, ladapter + seqs[i] + radapter, overlapped=True))
        counts[i] = len(matches)

def encodeCores(coresX, seqs, corelookup, ladapter, llen, indices):
    for j in range(coresX.shape[1]):
        coresX[indices,j] = [corelookup[tuple(seq[ladapter+j:llen+j])] for seq in seqs[indices]]

def getExpectedDist(seqs, order, trainKCounts):
    counts1 = trainKCounts[order]
    counts2 = trainKCounts[order-1]
    counts1 = dict(zip(counts1.keys(), np.array(list(counts1.values())) / np.sum(list(counts1.values()))))
    counts2 = dict(zip(counts2.keys(), np.array(list(counts2.values())) / np.sum(list(counts2.values()))))

    pred = np.array([np.prod([counts1[seq[j:j+order+1]] for j in range(len(seq) - order)]) for seq in seqs])
    argfilter = pred > 0
    if(order > 0):
        pred[argfilter] /= np.array([np.prod([counts2[seq[j:j+order]] for j in range(1,len(seq) - order)]) for seq in seqs])[argfilter]
    return pred

def getShiftExpectedDist(seqs, trainKCounts, order, shift, onRC):
    counts1 = trainKCounts[order]
    counts2 = trainKCounts[order-1]
    if(onRC):
        counts1 = counts1.iloc[:,::-1]
        counts1.index = rc(counts1.index.values)
        if(order > 0):
            counts2 = counts2.iloc[:,::-1]
            counts2.index = rc(counts2.index.values)
    counts1 = counts1.iloc[:,shift:]
    if(order > 0):
        counts2 = counts2.iloc[:,shift:]

    pred = np.array([[seq[j:j+order+1] for j in range(len(seq) - order)] for seq in seqs])
    pred = np.array([counts1.loc[pred[:,j]].iloc[:,j].values for j in range(len(seqs[0]) - order)])
    pred = np.product(pred, axis=0)
    argfilter = pred > 0
    if(order > 0):
        div = np.array([[seq[j:j+order] for j in range(1,len(seq) - order)] for seq in seqs[argfilter]])
        div = np.array([counts2.loc[div[:,j]].iloc[:,j].values for j in range(len(seqs[0]) - order-1)])
        div = np.product(div, axis=0)
        pred[argfilter] /= div
    return pred

def trainMarkov(r0trainfile, r0testfile, k):
    r0train = []
    r0test = []
    trainKCounts = []
    testKCounts = []

    for k_ in range(1, k+1):
        if(len(r0test) == 0):
            print("Reading " + r0testfile + " . . .")
            r0test, r0testcounts = readScores(r0testfile)
        print("Counting k = " + str(k_))
        testCounts = countKmers(r0test, r0testcounts, k_, ncpu=0)
        if(np.min(list(testCounts.values())) < 100):
            break
        testKCounts += [testCounts]
        if(len(r0train) == 0):
            print("Reading " + r0trainfile + " . . .")
            r0train, r0traincounts = readScores(r0trainfile)
        trainCounts = countKmers(r0train, r0traincounts, k_, ncpu=0)
        trainKCounts += [trainCounts]

    mmkMax = 0
    r2max = np.NINF
    for k_ in range(1,len(testKCounts)+1):
        print("Markov Model Order: " + str(k_ - 1))
        test = testKCounts[-1]
        EDist = getExpectedDist(list(test.keys()), k_-1, trainKCounts)
        testDist = np.array(list(test.values()))
        testDist /= np.sum(testDist)
        r2 = metrics.r2_score(testDist, EDist)
        
        print("R² = " + str(r2))
        if(r2 > r2max):
            mmkMax = k_
            r2max = r2
        else: 
            break
    
    orderMax = mmkMax - 1
    print("Best Markov Order: " + str(orderMax))
    print("R² = " + str(r2max) + "\n")

    # Use all counts in final model
    for i in range(len(trainKCounts)):
        for key, value in testKCounts[i].items():
            trainKCounts[i][key] += value

    return orderMax, trainKCounts

