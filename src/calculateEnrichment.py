#!/usr/bin/env python

import numpy as np
import pandas as pd
import tools
import argparse

parser = argparse.ArgumentParser(description="Calculates the round-corrected relative enrichment of k-mers in the provided input, relative to the expected frequency based on a Markov model trained on the initial library. The optimal Markov model order is determined to be that which performs best on the test set, which is not seen during training. ")
parser.add_argument('R0_train', help='Counts for unique reads from the initial library, split for Markov model training/testing')
parser.add_argument('R0_test', help='Counts for unique reads from the initial library, split for Markov model training/testing')
parser.add_argument('selected_input', help='Counts for unique reads from the selected sample')
parser.add_argument('SELEX_round', type=int, help='Number of rounds of selection between input and R0')
parser.add_argument('klen', type=int, help='Desired k-mer length')
args = parser.parse_args()

r0trainfile = args.R0_train
r0testfile = args.R0_test
rnfile = args.selected_input
cycles = args.SELEX_round
k = args.klen

rnname = rnfile.split("_")[0]
rntype = "_".join(rnfile[:-4].split("_")[1:])

orderMax, trainKCounts = tools.trainMarkov(r0trainfile, r0testfile, k)

print("Reading " + rnfile + " . . .")
rn, rny = tools.readScores(rnfile)
print("Counting k = " + str(k))
rnCounts = tools.countKmers(rn, rny, k, ncpu=0)
df = pd.DataFrame({"seqs":rnCounts.keys(), "counts":rnCounts.values()})
df['dist'] = df['counts']/df['counts'].sum()
df = df[df['counts'] >= 100]
df['dist_expected'] = tools.getExpectedDist(df['seqs'].values, orderMax, trainKCounts)
df['enrichment'] = df['dist']/df['dist_expected']
df['enrichment'] = np.power(df['enrichment'], 1/cycles)
df['enrichment_rc_ave'] = tools.consolidate(df['seqs'].values, df['enrichment'].values)
df['enrichment_relative'] = df['enrichment_rc_ave']/df['enrichment_rc_ave'].max()
df = df.sort_values('enrichment_rc_ave', ascending=False)
df = df[['seqs', 'enrichment_relative', 'enrichment_rc_ave', 'enrichment']]
df.to_csv(rnfile[:-4] + "_k" + str(k) + ".tsv", index=False, sep='\t')