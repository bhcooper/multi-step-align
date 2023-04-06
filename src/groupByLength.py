#!/usr/bin/env python

import sys
import numpy as np
import pandas as pd

data = pd.read_csv(sys.argv[1], sep="\t")
data.iloc[:,0] = data.iloc[:,0].str.replace("-", "")
data['klen'] = data.iloc[:,0].str.len()
for index, group in data.groupby('klen'):
    group.iloc[:,:-1].to_csv(sys.argv[1][:-4] + "_k" + str(index) + ".tsv", sep='\t', index=False)