#!/usr/bin/env python

import sys
import numpy as np
import sharedmem as sm
import regex as re
import tools

from pandas import read_csv

config = tools.loadConfig(sys.argv[1])
ncpu = config['ncpu']
ladapter = config['ladapter']
radapter = config['radapter']

# Using top results from SELEXisolate.py
ucores = read_csv(sys.argv[3], sep="\t")
ucores = ucores.iloc[:,0]
ncores = len(ucores)
ucores = np.append(ucores, tools.rc(ucores))
corelens = np.array([len(c) for c in ucores])
maxcorelen = max(corelens)
query = "|".join(ucores)

# Check adapters for cores
if(not re.search(query, ladapter) == None):
    print("Error: core " + re.search(query, ladapter)[0] + " found in 5' adapter")
    exit()
if(not re.search(query, radapter) == None):
    print("Error: core" + re.search(query, radapter)[0] + " found in 3' adapter")
    exit()

# Only include adapter if it involves part of variable region
ladapter = ladapter[-(maxcorelen-1):]
radapter = radapter[:(maxcorelen-1)]
ladapterlen = len(ladapter)
radapterlen = len(radapter)

print("Reading input . . .")

seqs = read_csv(sys.argv[2], sep="\t")
y = sm.copy(seqs.iloc[:,1].values)
seqs = sm.copy(seqs.iloc[:,0].values)
seqlen = len(seqs[0])
nwindows = seqlen - maxcorelen + 1

output = np.zeros((len(np.arange(5, ncores+1, 5)), 4))
for i,ncores_temp in enumerate(np.arange(5, ncores+1, 5)):
    print("Counting hits using " + str(ncores_temp) + " cores . . .")
    ucores_temp = ucores[:ncores_temp]
    ucores_temp = np.append(ucores_temp, tools.rc(ucores_temp))
    corelens_temp = np.array([len(c) for c in ucores_temp])
    query_temp = "|".join(ucores_temp)
    counts = sm.empty(len(seqs), dtype=np.uint8)
    tools.parallelAsync(tools.countHits, ncpu, counts, seqs, ladapter, radapter, query_temp)
    counts[counts > 1] = 2
    counts = np.unique(counts, return_counts=True)
    output[i,0] = ncores_temp
    output[i,counts[0] + 1] = counts[1]

tools.saveTable(sys.argv[2][:-4] + "_hitCounts.tsv", output.astype(str), header = "# cores\t0 hits\t1 hit\t>1 hit")
