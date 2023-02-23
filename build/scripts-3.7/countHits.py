#!python

import sys
import numpy as np
import pandas as pd
import regex as re
from Bio import SeqIO

N = ["A", "C", "G", "T"]
comp = {"A":"T", "C":"G", "G":"C", "T":"A", "N":"N"}
def rc(seqs):
    return np.array(["".join([comp[x] for x in seq][::-1]) for seq in seqs])

seqs = [str(seq_record.seq) for seq_record in SeqIO.parse(sys.argv[1], "fasta") if not seq_record.id == "chrM"]

cores = pd.read_csv(sys.argv[2], sep='\t').iloc[:,0]
cores = [seq.replace("-", "") for seq in cores]
query = ("|").join(cores) + ("|") + ("|").join(rc(cores))
rcset = set(rc(cores))

matches = np.array([[m for m in re.finditer(query, seq, overlapped=True)] for seq in seqs], dtype=object)
hitCounts = np.array([len(m) for m in matches])

matches = np.concatenate(matches)
matches = [m[0] for m in matches]
matches = [rc([seq])[0] if seq in rcset else seq for seq in matches]
matches, counts = np.unique(matches, return_counts=True)
matches = pd.DataFrame({"seq":matches, "count":counts})
matches = matches.sort_values("count", ascending=False)
matches.to_csv(sys.argv[1].split(".")[0] + "_hits.tsv", sep='\t', index=False)