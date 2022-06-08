
# Count unique reads and save as tsv
# Input is expected to end with .fastq.gz extension
# For next steps, adapters should be trimmed before this step
# Reads including "N" (ambiguous bases) are removed
countUnique.py R0.fastq.gz
countUnique.py Fkh1_R1.fastq.gz
countUnique.py Fkh1_R2.fastq.gz

# Divide unique reads into two non-overlapping sets of equally sized datasets (subsets 1 & 2)
# To be used for training and testing various order markov models
split.py R0.tsv

# Get the distribution of k-mers along the variable region and save for use in position specific Markov modeling of the core
# Markov model order set in config.yml
getKmerDist.py R0.tsv

# Set R0 train and test files and use best Markov model to calculate enrichment for various length k-mers
# R-th root of the enrichment is provided for later rounds of sequencing (round number provided as final argument)
analyzeSELEXCounts.py config.yml Fkh1_R1.tsv 1
analyzeSELEXCounts.py config.yml Fkh1_R2.tsv 2

# Aligns sequences based on enrichment scores using top-down crawl method described in publication (Cooper et al, 2022)
# R1 used to include more moderate-low affinity sequences
pip install TopDownCrawl
TopDownCrawl enrichment/Fkh1/k9/affinities_R1.tsv

# Larger k-mer length was used for alignment, but is trimmed to 7-bp region centered on known consensus for list of candidate cores
# Trimmed shift 0 from affinities_R1_aligned.tsv to 7-bp region covering GTAAACA (or 6-bp region covering GACGCA for Fhl1)
# Reprioritize the list of aligned sequences by iterative filtering as described in publication (Cooper et al, 2022)
# Only top 100 candidates are returned by default
SELEXisolate.py config.yml Fkh1_R1.tsv enrichment/Fkh1/k9/affinities_R1_aligned_trimmed.tsv

# Count how many sequences can be aligned given a specified list of cores
countAlignToCores.py config.yml Fkh1_R1.tsv Fkh1_R1_k7_isolates.tsv
countAlignToCores.py config.yml Fkh1_R2.tsv Fkh1_R1_k7_isolates.tsv

# The final list of candidate cores includes isolates from Fkh1, Fkh2, Hcm1, and Fhl1 datasets based on the criteria described in publication (Cooper et al, 2022)
# This full list is used to align all reads from each sample
alignToCores.py config.yml R0.tsv allcores.tsv
alignToCores.py config.yml Fkh1_R1.tsv allcores.tsv
alignToCores.py config.yml Fkh1_R2.tsv allcores.tsv

# Calculate the core and edge ddG/RT for each protein as described in the publication (Cooper et al, 2022)
# Also generates several figures to summarize the calculated values
# Removes any cores that are missing information for any edge position (ie: pseudopalindromes, insufficient read depth)
# The square root of the enrichment is used in all calculations for R2
analyzeAlignedSELEX.py config.yml R0_allcores.tsv Fkh1_R1_allcores.tsv Fkh1_R2_allcores.tsv