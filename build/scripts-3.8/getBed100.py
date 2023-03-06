#!python

import sys
import numpy as np
import pandas as pd
from Bio import SeqIO

bed = [str(seq_record.id).split("__") for seq_record in SeqIO.parse(sys.argv[1], "fasta")]

# Get center of peaks and add span to both sides
span = 50
bed = [[seq[0], int(np.mean([int(seq[1]),int(seq[2])]))] for seq in bed]
bed = [[seq[0], seq[1]-span, seq[1]+span] for seq in bed]
bed = pd.DataFrame(bed)

# Sort bed output
bed = bed.sort_values([bed.columns[0], bed.columns[1], bed.columns[2]])

# Output bed file
bed.to_csv(sys.argv[1].split(".")[0] + "_100.bed", header=None, index=False, sep='\t')