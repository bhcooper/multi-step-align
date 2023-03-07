#!/usr/bin/env python

import numpy as np
import pandas as pd
import sharedmem as sm
import tools
import argparse
from collections import defaultdict

parser = argparse.ArgumentParser(description="Sorts candidate cores based on an interative framework in which reads containing identified cores are iteratively removed so that the enrichment of one core is independent of another.")
parser.add_argument('R0_input', help='Counts for unique reads from the initial library')
parser.add_argument("selected_input", help="Tab-delimited file containing counts for unique reads")
parser.add_argument('SELEX_round', type=int, help='Number of rounds of selection between input and R0')
parser.add_argument("candidate_cores", help="List of candidate cores to be prioritized")
parser.add_argument("left_adapter", help="Fixed adapter 5' of the variable region")
parser.add_argument("right_adapter", help="Fixed adapter 3' of the variable region")
parser.add_argument("stopping_parameter", type=float, help="Should be a float between 0 and 1, definining when the algorithm should stop adding to the list of prioritized cores, based on the enrichment of the most recently prioritized cores and the average enrichment of the  previously prioritized cores")

args = parser.parse_args()
r0file = args.R0_input
rnfile = args.selected_input 
cycle = args.SELEX_round
candidates = args.candidate_cores
ladapter = args.left_adapter
radapter = args.right_adapter
stop = args.stopping_parameter

ucores = tools.read_csv(candidates, sep="\t", usecols=(0,)).iloc[:,0]
rccores = tools.rc(ucores)
ucores = np.append(ucores, rccores[np.logical_not(np.isin(rccores, ucores))])

corelen = len(ucores[0])
ucores = np.append(["N" * corelen],ucores)
corelookup = defaultdict(np.int32, zip([tuple(u) for u in ucores], np.arange(len(ucores))))

ladapter = ladapter[-(corelen-1):]
radapter = radapter[:(corelen-1)]
ladapterlen = len(ladapter)
radapterlen = len(radapter)

ncpu = 0

print("Reading input . . .")

r0seqs, r0y = tools.readScores(r0file)
all7, counts7 = np.unique([s[:7] for s in r0seqs], return_counts=True)
all7 = {k:v for k, v in zip(all7, counts7)}
r0seqs = sm.copy(np.array([ladapter + s + radapter for s in r0seqs]))

rnseqs, rny = tools.readScores(rnfile)
rnseqs = sm.copy(np.array([ladapter + s + radapter for s in rnseqs]))
seqlen = len(rnseqs[0])
nwindows = seqlen - corelen + 1
varcols = np.arange(nwindows)[ladapterlen:-radapterlen]
nonvar = np.delete(np.arange(nwindows), varcols)

total = np.sum(rny)
total0 = np.sum(r0y)

print("Encoding cores . . .")

coresX = sm.empty((len(rnseqs), nwindows), dtype=np.int32)
tools.parallelAsync(tools.encodeCores, ncpu, coresX, rnseqs, corelookup, 0, corelen)
cores0 = sm.empty((len(r0seqs), nwindows), dtype=np.int32)
tools.parallelAsync(tools.encodeCores, ncpu, cores0, r0seqs, corelookup, 0, corelen)

print("Counting cores . . .")

corecounts = np.apply_along_axis(np.bincount, axis=0, arr=coresX[:,varcols], weights = rny, minlength=len(ucores))
corecounts = corecounts[1:]
corecounts = pd.DataFrame(corecounts, index=ucores[1:])

corecounts0 = np.apply_along_axis(np.bincount, axis=0, arr=cores0[:,varcols], weights = r0y, minlength=len(ucores))
corecounts0 = corecounts0[1:]
corecounts0 = pd.DataFrame(corecounts0, index=ucores[1:])

tops = []
nums = []
maxE = None
while(len(coresX) > 0):
    coreE = (corecounts/total)/(corecounts0/total0)
    coreE = coreE.mean(axis=1, skipna=True).to_frame()
    coreE = tools.consolidatedf(coreE)
    coretop = coreE.idxmax(skipna=True)[0]
    if(not maxE):
        maxE = coreE.loc[coretop][0]
    coreE = coreE / maxE
    coreE = coreE.pow(1/cycle)
    argtop = corelookup[tuple(coretop)]
    coretop_rc = tools.rc([coretop])[0]
    argrc = corelookup[tuple(coretop_rc)]
    hits = np.logical_or(np.any(coresX == argtop, axis=1), np.any(coresX == argrc, axis=1))
    if(np.sum(rny[hits]) == 0):
        break 
    hits0 = np.logical_or(np.any(cores0 == argtop, axis=1), np.any(cores0 == argrc, axis=1))
    tops += [coretop]
    nums += [coreE.loc[coretop][0]]

    print(str(len(tops)) + ": " + tops[-1])
    print('E = %.2f' % nums[-1])
    if(len(nums) >= 5):
        r = (nums[-1]/np.mean(nums[-5:-1]))
        print('R = %.2f' % r)
        if(r > stop):
            break

    hitcounts = np.apply_along_axis(np.bincount, axis=0, arr=coresX[hits][:,varcols], weights = rny[hits], minlength=len(ucores))
    hitcounts = hitcounts[1:]
    hitcounts = pd.DataFrame(hitcounts, index=ucores[1:])

    # hitcounts0 = np.apply_along_axis(np.bincount, axis=0, arr=cores0[hits0][:,varcols], weights = r0y[hits0], minlength=len(ucores))
    # hitcounts0 = hitcounts0[1:]
    # hitcounts0 = pd.DataFrame(hitcounts0, index=ucores[1:])

    corecounts -= hitcounts
    # corecounts0 -= hitcounts0
    coresX = coresX[~hits]
    # cores0 = cores0[~hits0]
    rny = rny[~hits]
    # r0y = r0y[~hits0]

    corecounts = corecounts.drop(index = [coretop, coretop_rc])
    corecounts0 = corecounts0.drop(index = [coretop, coretop_rc])

print("Coverage: %.2f" % ((total-sum(rny))/total))
np.savetxt(rnfile[:-4] + "_topcores.tsv",np.array(tops)[:,np.newaxis], fmt="%s", delimiter="\t", header="core", comments="")