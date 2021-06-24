#!/usr/bin/env python

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.use('Agg')
import multiprocessing as mp
import sharedmem as sm
import yaml
import pickle
import regex as re

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

def plotGrid(filename, matrix, xticks=[], yticks=[], xtickrotation="horizontal", ytickrotation="horizontal", xlabel="", ylabel="", title="", cmap="RdYlBu", vmin=None, vmax=None, gridstridex=None, gridstridey = None, figsize=None, vline=None, ax=None, vlines = None):
    axGiven = True
    if(figsize == None):
        figsize = np.array(matrix.transpose().shape)/4
    if(ax == None):
        fig, ax = plt.subplots(figsize=figsize)
        axGiven = False
    ax.matshow(matrix, cmap=cmap, vmin=vmin, vmax=vmax, aspect="auto")
    if(vline):
        ax.axvline(vline, color='red', linewidth=2)
    if(vlines):
        [ax.axvline(line, color='firebrick', linewidth=4) for line in vlines]
    ax.set_xticks(range(matrix.shape[1]))
    ax.set_yticks(range(matrix.shape[0]))
    if(gridstridex != None):
        ax.set_xticks(np.arange(matrix.shape[1]+1)[::gridstridex] - 0.5, minor=True)
    if(gridstridey != None):
        ax.set_yticks(np.arange(matrix.shape[0]+1)[::gridstridey] - 0.5, minor=True)
    ax.set_xticklabels(xticks, fontsize=16, rotation=xtickrotation)
    ax.set_yticklabels(yticks, fontsize=18, rotation=ytickrotation, fontname="Courier New")
    ax.yaxis.tick_right()
    for edge, spine in ax.spines.items():
        spine.set_visible(False)
    ax.grid(which="minor", color="black", linewidth=0.5)
    ax.set_xlabel(xlabel, fontsize=18, fontweight="bold")
    ax.set_ylabel(ylabel, fontsize=18, fontweight="bold")
    ax.set_title(title)
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

def markov(seqs, kmerCounts, order, shift, onRC):
    counts1 = kmerCounts[order+1]
    counts2 = kmerCounts[order]
    if(onRC):
        counts1 = counts1.iloc[:,::-1]
        counts2 = counts2.iloc[:,::-1]
        counts1.index = rc(counts1.index.values)
        counts2.index = rc(counts2.index.values)
    counts1 = counts1.iloc[:,shift:]
    counts2 = counts2.iloc[:,shift:]

    pred = np.array([np.prod([counts1.loc[seq[j:j+order+1]].iloc[j] for j in range(len(seq) - order)]) for seq in seqs])
    argfilter = pred > 0
    if(order > 0):
        pred[argfilter] /= np.array([np.prod([counts2.loc[seq[j:j+order]].iloc[j] for j in range(1,len(seq) - order)]) for seq in seqs])[argfilter]
    return pred


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
    X = sm.empty((len(seqs), len(seqs[0])*4))
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
    for i in indices:
        for j in range(coresX.shape[1]):
            x = corelookup[tuple(seqs[i][ladapter+j:llen+j])]
            if(x):
                coresX[i,j] = x