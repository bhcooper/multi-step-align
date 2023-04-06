#!/usr/bin/env python

import sys
import tools
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description="Calculated the enrichment of aligned sequences spanning the core and several flanking positions, using R0 to train a Markov model.")
parser.add_argument('R0_input', help='Unaligned counts for unique reads from the initial library')
parser.add_argument("aligned_input", help="Selected reads aligned with alignToCores.py")
parser.add_argument('SELEX_round', type=int, help='Number of rounds of selection between aligned input and R0')
parser.add_argument("max_core_len", type=int, help="Max core length provided to alignToCores.py, for this analysis only sequences with cores of this length will be returned.")
parser.add_argument("left_flank_len", type=int, help="Number of bp to include 5' of the core")
parser.add_argument("right_flank_len", type=int, help="Number of bp to include 3' of the core")

args = parser.parse_args()
r0file = args.R0_input
rnfile = args.aligned_input
round = args.SELEX_round
coreLen = args.max_core_len
lflank = args.left_flank_len
rflank = args.right_flank_len

print("Reading input . . .")
r0seqs, r0y= tools.readScores(r0file)
r0kmerCounts = []
for i in range(1, len(r0seqs[0])):
    kmerCount, minCount = tools.getKmerDist(r0seqs, r0y, i)
    if minCount < 100:
        break
    r0kmerCounts += [kmerCount]
markovOrder = len(r0kmerCounts) - 1

RN =  pd.read_csv(rnfile, sep='\t')

print("Removed gapped cores")
RN = RN[~RN.iloc[:,0].str.contains('-')]

varLen = len(RN.iloc[0,0].replace("N", ""))
flankLen = varLen-coreLen

RN = RN[RN['shift'] >= lflank]
RN = RN[RN['shift'] <= flankLen-rflank]

RN.iloc[:,0] = RN.iloc[:,0].str.slice(flankLen-lflank, flankLen + coreLen + rflank)
RN.iloc[:,0] = RN.iloc[:,0].str.replace("-", "")

RN = RN.groupby(['seq', 'strand', 'shift']).sum()

# Only keep sequences that occur in at least 2 shifts
filt = RN.index.get_level_values(level=0).value_counts()
filt = filt[filt > 1].index
RN = RN.loc[filt]

RN = RN.reset_index().set_index(['shift', 'strand']).sort_index()
RN['E0'] = 0.0
for group in RN.index.unique():
    # print(f'Calculating background for group: {group}')
    shift, strand = group
    onRC = strand == "-"
    RN.loc[group, 'E0'] = tools.getShiftExpectedDist(RN.loc[group, 'seq'].values, r0kmerCounts, markovOrder, shift-lflank, onRC)


RN['E'] = RN.iloc[:,1]/RN['E0']
RN['E'] = RN['E'].pow(1/round)
RN = RN.dropna()
RN = RN.set_index('seq')
RN = RN.groupby(level='seq').mean()
RN['E'] /= RN['E'].max()
RN = RN.sort_values('E', ascending=False)
RN[['E']].to_csv(f'{sys.argv[2][:-4]}_-{lflank}_+{rflank}.tsv', sep='\t')