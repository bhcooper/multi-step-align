#!python

import sys
import os
import numpy as np
import tools

from sklearn import metrics

config = tools.loadConfig(sys.argv[1])
ncpu = config['ncpu']
r0trainfile = config['r0trainfile']
r0testfile = config['r0testfile']
rnfile = sys.argv[2]
cycles = int(sys.argv[3])

rnname = rnfile.split("_")[0]
rntype = "_".join(rnfile[:-4].split("_")[1:])

tools.mkdir("savedCounts")
tools.mkdir("enrichment")
tools.mkdir("enrichment/" + rnname)

def getExpectedDist(seqs, mmk, trainKCounts):
    counts1 = trainKCounts[mmk-1]
    counts2 = trainKCounts[mmk-2]
    counts1 = dict(zip(counts1.keys(), np.array(list(counts1.values())) / np.sum(list(counts1.values()))))
    counts2 = dict(zip(counts2.keys(), np.array(list(counts2.values())) / np.sum(list(counts2.values()))))

    pred = np.array([np.prod([counts1[seq[j:j+mmk]] for j in range(len(seq) - mmk + 1)]) for seq in seqs])
    argfilter = pred > 0
    if(mmk > 1):
        pred[argfilter] /= np.array([np.prod([counts2[seq[j:j+mmk-1]] for j in range(1,len(seq) - mmk + 1)]) for seq in seqs])[argfilter]

    return pred

varlen = len(tools.read_csv(r0trainfile, nrows=1, sep='\t', usecols=[0]).values[0][0])
kmax = 0
r0train = []
r0test = []
trainKCounts = []
testKCounts = []

for k in range(1, varlen):
    kmax = k - 1
    testCounts = None
    trainCounts = None
    minCountSummed = 0
    print("Counting k = " + str(k))
    picklefile = "savedCounts/" + r0testfile + ".k" + str(k) + ".pkl"
    if(os.path.exists(picklefile)):
        testCounts = tools.loadVariables(picklefile)
        print("Saved counts loaded")
    else:
        if(len(r0test) == 0):
            print("Reading " + r0testfile + " . . .")
            r0test, r0testcounts = tools.readScores(r0testfile)
        testCounts = tools.countKmers(r0test, r0testcounts, k, ncpu)
        tools.saveVariables(picklefile, testCounts)
    if(np.min(list(testCounts.values())) < 100):
        break
    testKCounts += [testCounts]
    picklefile = "savedCounts/" + r0trainfile + ".k" + str(k) + ".pkl"
    if(os.path.exists(picklefile)):
        trainCounts = tools.loadVariables(picklefile)
        print("Saved counts loaded")
    else:
        if(len(r0train) == 0):
            print("Reading " + r0trainfile + " . . .")
            r0train, r0traincounts = tools.readScores(r0trainfile)
        trainCounts = tools.countKmers(r0train, r0traincounts, k, ncpu)
        tools.saveVariables(picklefile, trainCounts)
    trainKCounts += [trainCounts]

print("\nk-max: " + str(kmax) + "\n")

# import pandas as pd

# print(type(trainKCounts[5]))
# x = pd.DataFrame.from_dict(testKCounts[4], orient='index').sort_values(0, ascending=False)
# x.iloc[:,0] /= x.iloc[:,0].sum()

# kmax = 6
# print("\nk/-max set to: " + str(kmax))

print("\n")
mmkMax = 0
r2max = np.NINF
for k in range(1,kmax+1):
    print("Markov Model Order: " + str(k - 1))
    test = testKCounts[-1]
    EDist = getExpectedDist(list(test.keys()), k, trainKCounts)
    # x = pd.DataFrame({"seqs": list(test.keys()), "EDist":EDist}).sort_values("EDist", ascending=False)
    testDist = np.array(list(test.values()))
    testDist /= np.sum(testDist)
    r2 = np.sqrt(metrics.r2_score(testDist, EDist))       
    
    print("R^2 = " + str(r2))
    if(r2 > r2max):
        mmkMax = k
        r2max = r2
    else: 
        break
    

print("Best Markov Order: " + str(mmkMax-1) + "\n")
print("R^2 = " + str(r2max) + "\n")

# Use all counts in final model
# for i in range(len(trainKCounts)):
    # for key, value in testKCounts[i].items():
        # trainKCounts[i][key] += value

rn = []
infokMax = 0
infoMax = 0
maxCounts = None
maxEDist = None
# for k in range(mmkMax, varlen+1):
print("k-mer max restricted to 13 bp")
for k in range(mmkMax, 14):
    picklefile = "savedCounts/" + rnfile + ".k" + str(k) + ".pkl"
    if(os.path.exists(picklefile)):
        rnCounts = tools.loadVariables(picklefile)
        print("Saved counts loaded")
    else:
        if(len(rn) == 0):
            print("Reading " + rnfile + " . . .")
            rn, rny = tools.readScores(rnfile)
        print("Counting k = " + str(k))
        rnCounts = tools.countKmers(rn, rny, k, ncpu)
        tools.saveVariables(picklefile, rnCounts)
    rnDist = np.array(list(rnCounts.values()))
    # print("100 count restriction removed")
    mask100 = rnDist >= 100
    if(np.sum(mask100) == 0):
        break
    total = np.sum(list(rnCounts.values()))
    rnDist /= np.sum(rnDist)
    rnDist = rnDist[mask100]
    # seqs = np.array(list(rnCounts.keys()))
    seqs = np.array(list(rnCounts.keys()))[mask100]
    # counts = np.array(list(rnCounts.values()))
    counts = np.array(list(rnCounts.values()))[mask100]
    Er0Dist = getExpectedDist(seqs, mmkMax, trainKCounts)
    Er0Counts = Er0Dist * total
    # x = pd.DataFrame({"seqs": seqs, "count": counts, "EDist":Er0Dist, "Ecounts":Er0Counts}).sort_values("count", ascending=False)
    info = np.sum(rnDist * np.log2(rnDist/Er0Dist))
    rnDistSum = np.sum(rnDist)
    Er0DistSum = np.sum(Er0Dist)
    if(not rnDistSum == 1 and not Er0DistSum == 1):
        info += (1 - rnDistSum) * np.log2((1 - rnDistSum)/(1 - Er0DistSum))
    print("Info Gain: k = " + str(k))
    print(info)
    Er0Counts = Er0Dist * total
    enrichment = rnDist/Er0Dist
    
    enrichment = np.power(enrichment, 1/cycles)

    se = enrichment * np.sqrt(2/counts)
    rcEnrichment = tools.consolidate(seqs, enrichment)
    argSort = np.argsort(-rcEnrichment)
    seqs = seqs[argSort]
    rnDist = rnDist[argSort]
    counts = counts[argSort]
    Er0Counts = Er0Counts[argSort]
    Er0Dist = Er0Dist[argSort]
    enrichment = enrichment[argSort]
    rcEnrichment = rcEnrichment[argSort]
    relEnrichment = rcEnrichment/np.max(rcEnrichment)
    se = se[argSort]
    path = "enrichment/" + rnname + "/k" + str(k)
    tools.mkdir(path)
    outfile = open(path + "/info_" + rntype + "=" + "%.2f" % info, 'w')
    outfile.close()
    outfile = open(path + "/affinities_" + rntype + ".tsv", 'w')
    outfile.write("Seqs\tRelative\tEnrichment_rc_ave\tEnrichment\trnCounts\tE(r0Counts)\trnDist\tE(r0Dist)\tSE\n")
    for line in zip(seqs.astype(str), relEnrichment.astype(str), rcEnrichment.astype(str), enrichment.astype(str), counts.astype(str), Er0Counts.astype(str), rnDist.astype(str), Er0Dist.astype(str), se.astype(str)):
        outfile.write("\t".join(line) + "\n")
    outfile.close()