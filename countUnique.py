#!/usr/bin/env python

import sys
import gzip
import numpy as np

# Read sequence lines from gzip file
seqs = np.array([line.strip() for i, line in enumerate(gzip.open(sys.argv[1], "rt", encoding='utf-8')) if i % 4 == 1])

unique,counts = np.unique(seqs, return_counts=True)
argfilter = np.array([not "N" in u for u in unique])
unique = unique[argfilter]
counts = counts[argfilter]

# Combine seqs and counts into a 2D array
temp = np.array([unique, counts]).transpose()
# Sort by counts, descending
temp = temp[np.argsort(-counts)]

np.savetxt(sys.argv[1][:-9] + ".tsv", temp, delimiter="\t", header="Seq\tCount", comments="", fmt="%s")