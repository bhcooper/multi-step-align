#!python

import sys
import tools
import pandas as pd

config = tools.loadConfig(sys.argv[1])
ladapter = config['ladapter']
radapter = config['radapter']
ladapter_rc = tools.rc([radapter])[0]
radapter_rc = tools.rc([ladapter])[0]
markovOrder = config['markovOrder']
coreLen = config['corelen']
r0kmerDist = tools.loadVariables(config['r0kmerDist'])

print("Reading input . . .")
# R0 = pd.read_csv(sys.argv[2], sep='\t')
RN =  pd.read_csv(sys.argv[2], sep='\t')
round = int(sys.argv[3])

print("Removed gapped cores")
# R0 = R0[~R0.iloc[:,0].str.contains('-')]
RN = RN[~RN.iloc[:,0].str.contains('-')]

# print("Remove non-gapped cores")
# RN = RN[RN.iloc[:,0].str.contains('-')]

lflank = int(sys.argv[4])
rflank = int(sys.argv[5])

varLen = len(RN.iloc[0,0].replace("N", ""))
flankLen = varLen-coreLen

# R0 = R0[R0['shift'] >= lflank]
RN = RN[RN['shift'] >= lflank]
# R0 = R0[R0['shift'] <= flankLen-rflank]
RN = RN[RN['shift'] <= flankLen-rflank]

# R0.iloc[:,0] = R0.iloc[:,0].str.slice(flankLen-lflank, flankLen + coreLen + rflank)
RN.iloc[:,0] = RN.iloc[:,0].str.slice(flankLen-lflank, flankLen + coreLen + rflank)
RN.iloc[:,0] = RN.iloc[:,0].str.replace("-", "")

# R0 = R0.groupby(['seq', 'strand', 'shift']).sum()
RN = RN.groupby(['seq', 'strand', 'shift']).sum()

# Only keep sequences that occur in all shifts
filt = RN.index.get_level_values(level=0).value_counts()
filt = filt[filt > 1].index
# filt = filt[filt == filt.max()].index
RN = RN.loc[filt]

RN = RN.reset_index().set_index(['shift', 'strand']).sort_index()
RN['E0'] = 0.0
for group in RN.index.unique():
    print(f'Calculating background for group: {group}')
    shift, strand = group
    onRC = strand == "-"
    RN.loc[group, 'E0'] = tools.markov(RN.loc[group, 'seq'].values, r0kmerDist, markovOrder, shift-lflank, onRC)


RN['E'] = RN.iloc[:,1]/RN['E0']
RN['E'] = RN['E'].pow(1/round)
RN = RN.dropna()
RN = RN.set_index('seq')
RN = RN.groupby(level='seq').mean()
RN['E'] /= RN['E'].max()
RN = RN.sort_values('E', ascending=False)
RN[['E']].to_csv(f'{sys.argv[2][:-4]}_-{lflank}_+{rflank}.tsv', sep='\t')

# Use R0 counts for background

# R0 = R0.iloc[:,:2].groupby(R0.columns[0]).sum()
# RN = RN.iloc[:,:2].groupby(RN.columns[0]).sum()

# E = RN/R0
# # E = RN
# E = E.pow(1/round)
# E /= E.max()
# E = E.sort_values(E.columns[0], ascending=False)
# E = E.dropna()
# E.to_csv(f'{sys.argv[3][:-4]}_-{lflank}_+{rflank}.tsv', sep='\t', header="seq\tE")