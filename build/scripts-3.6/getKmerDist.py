#!python

import sys
import tools

seqs, y = tools.readScores(sys.argv[1])

kmerCounts = [[]]
for i in range(1, len(seqs[0])):
    print("Counting " + str(i) + " . . .")
    kmerCount, minCount = tools.getKmerDist(seqs, y, i)
    if minCount < 100:
        break
    kmerCounts += [kmerCount]
tools.saveVariables(sys.argv[1][:-4] + "_kmerDist.pkl", kmerCounts)