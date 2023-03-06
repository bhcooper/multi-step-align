#!python

import sys
import numpy as np
import multiprocessing as mp
import sharedmem as sm
import tools

from collections import defaultdict

config = tools.loadConfig(sys.argv[1])
ncpu = config['ncpu']
ladapter = config['ladapter']
radapter = config['radapter']
markovOrder = config['markovOrder']
numisolates = config['numisolates']
r0kmerDist = tools.loadVariables(config['r0kmerDist'])

ucores = tools.read_csv(sys.argv[3], sep="\t", usecols=(0,)).iloc[:,0]
ucores = [seq.strip("N") for seq in ucores]

# ucores = np.unique([u[1:] for u in ucores])

ucores = np.append(ucores, tools.rc(ucores))
ncores = len(ucores)
numisolates = min(ncores, numisolates)
corelen = len(ucores[0])
ucores = np.append(["N" * corelen],ucores)
corelookup = defaultdict(np.int32, zip([tuple(u) for u in ucores], np.arange(len(ucores))))

ladapter = ladapter[-(corelen-1):]
radapter = radapter[:(corelen-1)]
ladapterlen = len(ladapter)
radapterlen = len(radapter)

print("Reading input . . .")

rnseqs, rny = tools.readScores(sys.argv[2])
rnseqs = sm.copy(np.array([ladapter + s + radapter for s in rnseqs]))
seqlen = len(rnseqs[0])
nwindows = seqlen - corelen + 1
varcols = np.arange(nwindows)[ladapterlen:-radapterlen]
nonvar = np.delete(np.arange(nwindows), varcols)
total = np.sum(rny)

print("Encoding cores . . .")

coresX = sm.empty((len(rnseqs), nwindows), dtype=np.int32)
tools.parallelAsync(tools.encodeCores, ncpu, coresX, rnseqs, corelookup, 0, corelen)
del rnseqs

print("Counting cores . . .")

corecounts = sm.empty((ncores, len(varcols)))
for i,col in enumerate(varcols):
    for c,y in zip(coresX[:,col], rny):
        if(c):
            corecounts[c-1,i] += y

Ecounts = sm.empty(corecounts.shape)
for i in range(Ecounts.shape[1]):
    Ecounts[:,i] = tools.markov(ucores[1:], r0kmerDist, markovOrder, i, False)

tops = []
nums = []
while(len(tops) < numisolates and len(coresX) > 0):
# while(len(coresX) > 0):
    coreE = (corecounts/total)/Ecounts
    coreE = np.mean(coreE, axis=1)
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
        # if(r > 0.95):
            # break


    hitcounts = sm.empty(corecounts.shape)
    for i,col in enumerate(varcols):
        for c,y in zip(coresX[hits][:,col], rny[hits]):
            if(c):
                hitcounts[c-1,i] += y
    corecounts -= hitcounts
    coresX = coresX[~hits]
    rny = rny[~hits]

print("Unmatched: %.2f" % (sum(rny)/total))
np.savetxt(sys.argv[2][:-4] + "_k" + str(corelen) + "_isolates.tsv",np.array([tops, nums]).transpose(), fmt="%s", delimiter="\t", header="core\tE", comments="")
# np.savetxt(sys.argv[2][:-4] + "_k" + str(corelen) + "_isolates.tsv",np.array(tops)[:,np.newaxis], fmt="%s", delimiter="\t", header="core", comments="")
