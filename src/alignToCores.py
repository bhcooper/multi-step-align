#!/usr/bin/env python

import numpy as np
import sharedmem as sm
import tools
import argparse
import pandas as pd

parser = argparse.ArgumentParser(description="aligns full-length reads to a given list of cores, discarding reads with multiple hits and hits which span the adapter")
parser.add_argument("selected_input", help="Tab-delimited file containing counts for unique reads")
parser.add_argument("cores", help="Newline separated list of cores with a header row")
parser.add_argument("left_adapter", help="Fixed adapter 5' of the variable region")
parser.add_argument("right_adapter", help="Fixed adapter 3' of the variable region")

args = parser.parse_args()

readfile = args.selected_input
corefile = args.cores
ladapter = args.left_adapter
radapter = args.right_adapter

ncpu = 0

# Using top results from SELEXisolate.py
ucores = pd.read_csv(corefile, sep="\t")
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

seqs = pd.read_csv(readfile, sep="\t")
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
queryname = corefile[:-4]
strand = np.full(len(seqs), "+")
strand[onRC] = "-"
np.savetxt(readfile[:-4] + "_" + queryname + ".tsv", np.array([seqs, y, strand, starts]).transpose(), delimiter="\t", header="seq\tcount\tstrand\tshift", comments="", fmt="%s")