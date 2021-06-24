#!/usr/bin/env python

import sys
import numpy as np
import tools
from pandas import read_csv

seqs = read_csv(sys.argv[1], sep='\t')
y = seqs.iloc[:,1]
onRC = seqs.iloc[:,2] == '-'
seqs = seqs.iloc[:,0]
seqLen = len(seqs[0])
varLen = len(seqs[0].replace("N", ""))
numShifts = seqLen - varLen + 1

for i in range(numShifts):
    argfilter = np.logical_not(np.array([seq[i] == "N" for seq in seqs]))
    dirname = sys.argv[1].split("_")[0]
    round = sys.argv[1].split("_")[1]
    tools.mkdir(dirname)
    filt = np.logical_and(argfilter, np.logical_not(onRC))
    tools.saveTable(dirname + "/" + dirname + "_shift" + str(numShifts-1-i) + "_+_" + round + ".tsv", np.array([[seq.replace("N", "") for seq in seqs[filt]], y[filt]]).transpose())
    filt = np.logical_and(argfilter, onRC)
    tools.saveTable(dirname + "/" + dirname + "_shift" + str(numShifts-1-i) + "_-_" + round + ".tsv", np.array([[seq.replace("N", "") for seq in seqs[filt]], y[filt]]).transpose())
    seqs = seqs[~argfilter]
    y = y[~argfilter]
    onRC = onRC[~argfilter]