#!python

import sys
import numpy as np
from pandas import read_csv

data = read_csv(sys.argv[1], sep="\t", usecols=[0,1])
subset1, subset2 = np.array_split(data.sample(frac=1), 2)
np.savetxt(sys.argv[1][:-4] + "_split1.tsv", subset1.iloc[np.argsort(-subset1.iloc[:,1])], delimiter="\t", header="seq\tcount", comments="", fmt="%s")
np.savetxt(sys.argv[1][:-4] + "_split2.tsv", subset2.iloc[np.argsort(-subset2.iloc[:,1])], delimiter="\t", header="seq\tcount", comments="", fmt="%s")