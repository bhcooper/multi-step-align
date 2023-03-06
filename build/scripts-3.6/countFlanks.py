#!python

import sys
import numpy as np
import pandas as pd
import regex as re
from Bio import SeqIO

N = ["A", "C", "G", "T"]
comp = {"A":"T", "C":"G", "G":"C", "T":"A", "N":"N"}
def rc(seq):
    return "".join([comp[x] for x in seq][::-1])

seqs = [str(seq_record.seq) for seq_record in SeqIO.parse(sys.argv[1], "fasta") if not seq_record.id == "chrM"]
core = sys.argv[2]
lFlankLen = int(sys.argv[3])
rFlankLen = int(sys.argv[4])
query = "[ACGT]" * lFlankLen + core + "[ACGT]" * rFlankLen
rcQuery = "[ACGT]" * rFlankLen + rc(core) + "[ACGT]" * lFlankLen

matches = np.array([[m[0] for m in re.finditer(query, seq, overlapped=True)] for seq in seqs], dtype=object)
rcMatches = np.array([[rc(m[0]) for m in re.finditer(rcQuery, seq, overlapped=True)] for seq in seqs], dtype=object)

matches = np.concatenate(matches)
rcMatches = np.concatenate(rcMatches)
matches = np.concatenate((matches, rcMatches))
matches = np.array([list(m) for m in matches])
X = np.array([np.sum(matches == c, axis = 0) for c in N]).transpose()
X = np.concatenate((X[:lFlankLen], X[-rFlankLen:]), axis = 0)
index = [str(i) for i in range(-lFlankLen, 0)] + [f"+{i}" for i in range(1, rFlankLen + 1)]
X = pd.DataFrame(X, index = index, columns = N)
X.to_csv(sys.argv[1].split(".")[0] + f"_{core}_-{lFlankLen}_+{rFlankLen}.tsv", sep="\t")