#!python

import sys
import numpy as np
import multiprocessing as mp
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
# ucores = ucores.iloc[:40,0]
ucores = ucores.iloc[:,0]

# Prepend N
# ucores = ["A" + u[:-1] for u in ucores]
# ucores = [c + u for u in ucores for c in tools.N]
# ucores = np.unique([u[1:] for u in ucores])
# ucores = [c + u for u in ucores for c in tools.N]

ncores = len(ucores)
ucores = np.append(ucores, tools.rc(ucores))
corelens = np.array([len(c) for c in ucores])
maxcorelen = max(corelens)
corelookup = dict(zip(ucores, np.arange(len(ucores))))
query = "|".join(ucores)

# Check adapters for cores
# if(not re.search(query, ladapter) == None):
#     print("Error: core " + re.search(query, ladapter)[0] + " found in 5' adapter")
#     exit()
# if(not re.search(query, radapter) == None):
#     print("Error: core" + re.search(query, radapter)[0] + " found in 3' adapter")
#     exit()

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

print("Finding hits . . .")

cores = sm.empty(len(seqs), dtype=np.uint16)
starts = sm.full(len(seqs), -1, dtype=np.int8)
tools.parallelAsync(tools.getHit, ncpu, cores, starts, seqs, ladapter, radapter, query, ncores, corelookup, corelens, maxcorelen)

print("Filtering unalignable reads . . .")
filt = np.logical_not(starts == -1)
seqs = seqs[filt]
cores = cores[filt]
starts = starts[filt]
print(("%.2f" % (100 * np.sum(y[filt])/np.sum(y))) + "% of reads can be validly aligned (excluding cores that span the fixed adapters)")
y = y[filt]

print("Orienting reverse complements . . .")
onRC = cores >= ncores
seqs[onRC] = tools.rc(seqs[onRC])
cores[onRC] -= ncores
starts[onRC] = seqlen - starts[onRC] - corelens[cores[onRC]]

print("Gapping short cores . . .")
ngaps = maxcorelen - corelens[cores]
shortCore = np.logical_not(ngaps == 0)
seqs[shortCore] = [seq[:start+clen] + "-" * ngap + seq[start+clen:][:-ngap] for seq,start,clen,ngap in zip(seqs[shortCore], starts[shortCore], corelens[cores[shortCore]], ngaps[shortCore])]

print("Aligning all seqs . . .")
seqs = np.array(["N" * (nwindows - start - 1) + seq + "N" * (start) for seq,start in zip(seqs, starts)])

print("Saving output . . .")
queryname = sys.argv[3][:-4]
strand = np.full(len(seqs), "+")
strand[onRC] = "-"
np.savetxt(sys.argv[2][:-4] + "_" + queryname + ".tsv", np.array([seqs, y, strand, starts]).transpose(), delimiter="\t", header="seq\tcount\tstrand\tshift", comments="", fmt="%s")