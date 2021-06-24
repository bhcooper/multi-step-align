#!/usr/bin/env python

import sys
import numpy as np
import pandas as pd
import tools
from collections import defaultdict

expectedSeq = "GTAAACA"

alts = {"A":["C", "G", "T"],"C":["A", "G", "T"],"G":["A", "C", "T"],"T":["A", "C", "G"]}

def countOccur(seqs, target):
    return sum(1 for seq in seqs if target in seq)

def fill(i, shift, change, iteration, crawl):
    crawl["Shift"][i] = shift
    crawl["Change"][i] = change
    crawl["Iteration"][i] = iteration
    crawl["Filled"][i] = True
    return crawl

def crawlSNPs(i, crawl, lookup, iteration):
    seq = crawl["Seqs"][i]
    for j,c in enumerate(seq):
        for d in alts[c]:
            temp = seq[0:j] + d + seq[j+1:]
            desti = lookup.pop(temp, -1)
            rci = lookup.pop(tools.rc([temp])[0], -1)
            if(not desti == -1):
                crawl["Checked"][rci] = True
                crawl = fill(desti, crawl["Shift"][i], abs(crawl["Score"][i]-crawl["Score"][desti]), iteration, crawl)
    return crawl

def crawlLeft(i, crawl, lookup, iteration):
    for c in tools.N:
        temp = c + crawl["Seqs"][i][:-1]
        desti = lookup.pop(temp, -1)
        rci = lookup.pop(tools.rc([temp])[0], -1)
        if(not desti == -1):
            crawl["Checked"][rci] = True
            crawl = fill(desti, crawl["Shift"][i]-1, abs(crawl["Score"][i]-crawl["Score"][desti]), iteration, crawl)
    return crawl
    
def crawlRight(i, crawl, lookup, iteration):
    for c in tools.N:
        temp = crawl["Seqs"][i][1:] + c
        desti = lookup.pop(temp, -1)
        rci = lookup.pop(tools.rc([temp])[0], -1)
        if(not desti == -1):
            crawl["Checked"][rci] = True
            crawl = fill(desti, crawl["Shift"][i]+1, abs(crawl["Score"][i]-crawl["Score"][desti]), iteration, crawl)
    return crawl

def consolidate(lookup):
    for key in list(lookup.keys()):
        if(lookup[tools.rc([key])[0]] == 0.0):
            lookup[tools.rc([key])[0]] = lookup[key]
        else:
            ave = (lookup[key] + lookup[tools.rc([key])[0]]) / 2
            lookup[key] = ave
            lookup[tools.rc([key])[0]] = ave
    return lookup
        
seqs, scores = tools.readScores(sys.argv[1])
seqs = np.array([seq.replace("N", "") for seq in seqs])
lookup = defaultdict(float, zip(seqs, scores))
lookup = consolidate(lookup)
seqs = np.array(list(lookup.keys()))
scores = np.array(list(lookup.values()))

argsort = np.argsort(-scores)
seqs = seqs[argsort]
scores = scores[argsort]
crawl = {"Seqs":seqs, "Score":scores, "Shift":np.zeros(len(lookup), dtype=int), "Change":np.zeros(len(lookup)),
    "Iteration":np.zeros(len(lookup), dtype=int), "Filled":np.full(len(lookup), False), "Checked":np.full(len(lookup), False)}

lookup = defaultdict(lambda: -1, zip(crawl["Seqs"], range(len(crawl["Seqs"]))))
iteration = 0
crawl = fill(0,0,0,0,crawl)
print("Aligning", len(lookup), "sequences . . .")
lookup.pop(crawl["Seqs"][0])
if(not crawl["Seqs"][0] == tools.rc([crawl["Seqs"][0]])[0]):
    crawl["Checked"][lookup.pop(tools.rc([crawl["Seqs"][0]])[0])] = True
nextMax = np.argmax(~crawl["Checked"] * crawl["Filled"])
while(crawl["Checked"][nextMax] == False and len(lookup) > 0):
    crawl["Checked"][nextMax] = True
    iteration += 1
    crawl = crawlSNPs(nextMax, crawl, lookup, iteration)
    crawl = crawlLeft(nextMax, crawl, lookup, iteration)
    crawl = crawlRight(nextMax, crawl, lookup, iteration)
    nextMax = np.argmax(~crawl["Checked"] * crawl["Filled"])
print("Alignment complete")
print("Unable to align",len(lookup),"sequences")


# Remove unfilled
filter = np.where(crawl["Filled"] == True)
crawl.pop("Filled")
crawl.pop("Checked")
crawl["Seqs"] = crawl["Seqs"][filter]
crawl["Score"] = crawl["Score"][filter]
crawl["Shift"] = crawl["Shift"][filter]
crawl["Change"] = crawl["Change"][filter]
crawl["Iteration"] = crawl["Iteration"][filter]

rcs = tools.rc(crawl["Seqs"])
if(countOccur(crawl["Seqs"], expectedSeq) < countOccur(rcs, expectedSeq)):
    crawl["Seqs"] = rcs
    crawl["Shift"] = -crawl["Shift"]

crawl = pd.DataFrame(crawl)
minShift = -np.min(crawl["Shift"])
maxShift = np.max(crawl["Shift"])
for i, item in crawl.iterrows():
    shift = item[2]
    crawl.at[i, "Seqs"] = "N" * (minShift+shift) + item[0] + "N" * (maxShift-shift)

crawl.to_csv(sys.argv[1][:-4] + "_crawled.tsv", index=False, header=crawl.columns, sep="\t")