#!python

import sys
import numpy as np
import sharedmem as sm
import tools
import argparse
from collections import defaultdict

parser = argparse.ArgumentParser(description="Sorts candidate cores based on an interative framework in which reads containing identified cores are iteratively removed so that the enrichment of one core is independent of another.")
parser.add_argument('R0_train', help='Counts for unique reads from the initial library, split for Markov model training/testing')
parser.add_argument('R0_test', help='Counts for unique reads from the initial library, split for Markov model training/testing')
parser.add_argument("selected_input", help="Tab-delimited file containing counts for unique reads")
parser.add_argument('SELEX_round', type=int, help='Number of rounds of selection between input and R0')
parser.add_argument("candidate_cores", help="List of candidate cores to be prioritized")
parser.add_argument("left_adapter", help="List of candidate cores to be prioritized")
parser.add_argument("right_adapter", help="List of candidate cores to be prioritized")
parser.add_argument("stopping_parameter", type=float, help="Should be a float between 0 and 1, definining when the algorithm should stop adding to the list of prioritized cores, based on the enrichment of the most recently prioritized cores and the average enrichment of the  previously prioritized cores")

args = parser.parse_args()
r0trainfile = args.R0_train
r0testfile = args.R0_test
rnfile = args.selected_input 
cycle = args.SELEX_round
candidates = args.candidate_cores
ladapter = args.left_adapter
radapter = args.right_adapter
stop = args.stopping_parameter

ucores = tools.read_csv(candidates, sep="\t", usecols=(0,)).iloc[:,0]
ucores = np.append(ucores, tools.rc(ucores))
ncores = len(ucores)
corelen = len(ucores[0])
ucores = np.append(["N" * corelen],ucores)
corelookup = defaultdict(np.int32, zip([tuple(u) for u in ucores], np.arange(len(ucores))))

ladapter = ladapter[-(corelen-1):]
radapter = radapter[:(corelen-1)]
ladapterlen = len(ladapter)
radapterlen = len(radapter)

ncpu = 0

print("Training Markov Model . . .")

orderMax, trainKCounts = tools.trainMarkov(r0trainfile, r0testfile, corelen)

print("Reading input . . .")

r0seqs1, r0y1 = tools.readScores(r0trainfile)
r0seqs2, r0y2 = tools.readScores(r0trainfile)
r0seqs = np.concatenate((r0seqs1, r0seqs2), axis=0)
r0y = np.concatenate((r0y1, r0y2), axis=0)
del r0seqs1, r0seqs2, r0y1, r0y2

rnseqs, rny = tools.readScores(rnfile)
rnseqs = sm.copy(np.array([ladapter + s + radapter for s in rnseqs]))
seqlen = len(rnseqs[0])
nwindows = seqlen - corelen + 1
varcols = np.arange(nwindows)[ladapterlen:-radapterlen]
nonvar = np.delete(np.arange(nwindows), varcols)
total = np.sum(rny)

print("Encoding cores . . .")

coresX = sm.empty((len(rnseqs), nwindows), dtype=np.int32)
tools.parallelAsync(tools.encodeCores, ncpu, coresX, rnseqs, corelookup, 0, corelen)
# cores0 = sm.empty((len(r0seqs), nwindows), dtype=np.int32)
# tools.parallelAsync(tools.encodeCores, ncpu, cores0, r0seqs, corelookup, 0, corelen)

print("Counting cores . . .")

corecounts = sm.empty((ncores, len(varcols)))
for i,col in enumerate(varcols):
    for c,y in zip(coresX[:,col], rny):
        if(c):
            corecounts[c-1,i] += y

# Ecounts = sm.empty(corecounts.shape)
# for i in range(Ecounts.shape[1]):
    # Ecounts[:,i] = tools.markov(ucores[1:], r0kmerDist, markovOrder, i, False)

Edist = tools.getExpectedDist(ucores[1:], orderMax, trainKCounts)

tops = []
nums = []
while(len(coresX) > 0):
    coreE = (corecounts/total)
    coreE = np.mean(coreE, axis=1)
    coreE = coreE/Edist
    coreE = tools.consolidate(ucores[1:], coreE)
    argtop = np.argmax(coreE)+1
    # if(coreE[argtop-1] == 0):
        # break

    argrc = corelookup[tuple(tools.rc([ucores[argtop]])[0])]
    hits = np.logical_or(np.any(coresX == argtop, axis=1), np.any(coresX == argrc, axis=1))
    if(np.sum(rny[hits]) == 0):
        break 
    tops += [ucores[argtop]]
    nums += [coreE[argtop-1]]
    # nums += [int(np.sum(rny[hits]))]
    
    print(str(len(tops)) + ": " + tops[-1])
    print('E = %.2f' % nums[-1])
    if(len(nums) >= 5):
        r = (nums[-1]/np.mean(nums[-5:-1]))
        print('R = %.2f' % r)
        if(r > stop):
            break

    hitcounts = sm.empty(corecounts.shape)
    for i,col in enumerate(varcols):
        for c,y in zip(coresX[hits][:,col], rny[hits]):
            if(c):
                hitcounts[c-1,i] += y

    # Update Markov counts
    

    corecounts -= hitcounts
    coresX = coresX[~hits]
    rny = rny[~hits]

print("Coverage: %.2f" % ((total-sum(rny))/total))
np.savetxt(sys.argv[2][:-4] + "_topcores.tsv",np.array(tops)[:,np.newaxis], fmt="%s", delimiter="\t", header="core", comments="")