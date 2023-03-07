#!/usr/bin/env python

import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='Combines provided lists, removing duplicates')
parser.add_argument('core_list', nargs='+', help='Newline separated list of cores each with a header row')

args = parser.parse_args()
files = args.core_list
for f in files:
    print(f)

df = pd.concat([pd.read_csv(f, sep='\t') for f in files]).drop_duplicates()
df.to_csv("allcores.tsv", sep='\t', index=False, header="core")
