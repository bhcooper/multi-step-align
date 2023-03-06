#!/usr/bin/env python

import numpy as np
import pandas as pd
import argparse
import re
import tools

parser = argparse.ArgumentParser(description='Trims the given alignment to the region covering the consensus based on the location of its occurence within the most enriched sequence, reverse-complementing all sequences if necessary.')
parser.add_argument('input_alignment', help='tab-delimited alignment from TopDownCrawl with sequences in the first column and enrichment in the second')
parser.add_argument('consensus', help='DNA sequence to be trimmed to')

args = parser.parse_args()
fname = args.input_alignment
consensus = args.consensus
rcConsensus = tools.rc([consensus])[0]

df = pd.read_csv(fname, sep="\t", usecols=[0,1])
df = df.sort_values(df.columns[1], ascending=False)
pattern = f'{consensus}|{rcConsensus}'

start = -1
onRc = False
for seq in df.iloc[:,0]:
    match = re.search(pattern, seq)
    if match:
        start = match.start()
        if(seq[start:start+len(consensus)] == rcConsensus):
            onRc = True
        break

if(start == -1):
    print("Consensus was not found. Please check your input.")
else:
   df.iloc[:,0] = df.iloc[:,0].str[start:start+len(consensus)]

df = df.drop(columns=[df.columns[1]])
df = df.drop_duplicates()
df = df[np.logical_not(df.iloc[:,0].str.contains("_"))]
if(onRc):
    df.iloc[:,0] = tools.rc(df.iloc[:,0].values)
df.to_csv(f'{fname[:-4]}_{consensus}.tsv', sep='\t', index=False)
