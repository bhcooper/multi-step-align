
# Count unique reads and save as tsv
# Input is expected to end with .fastq.gz extension
# For next steps, adapters should be trimmed before this step
# Reads including "N" (ambiguous bases) are removed
countUnique.py R0.fastq.gz
countUnique.py Fkh1_R1.fastq.gz
countUnique.py Fkh1_R2.fastq.gz
countUnique.py Fkh2_R1.fastq.gz
countUnique.py Fkh2_R2.fastq.gz

# Divide unique reads into two non-overlapping sets of equally sized datasets (subsets 1 & 2)
# To be used for training and testing various order markov models
split.py R0.tsv

# Set R0 train and test files and use best Markov model to calculate enrichment for various length k-mers
# If "R2" is in file name, the square root of the enrichment is provided in the output to better represent relative affinities
analyzeSELEXCounts.py config.yml Fkh1_R1.tsv
analyzeSELEXCounts.py config.yml Fkh1_R2.tsv
analyzeSELEXCounts.py config.yml Fkh2_R1.tsv
analyzeSELEXCounts.py config.yml Fkh2_R2.tsv

# Aligns sequences based on enrichment scores using top-down crawl method described in publication (Cooper et al, 2021)
# R1 used to include more moderate-low affinity sequences
topDownCrawl.py enrichment/Fkh1/k7/affinities_R1.tsv
topDownCrawl.py enrichment/Fkh2/k7/affinities_R1.tsv

# Get the distribution of k-mers along the variable region and save for use in position specific Markov modeling of the core
# Markov model order set in config.yml
getKmerDist.py R0.tsv

# Reprioritize the list of aligned sequences by iterative filtering as described in publication (Cooper et al, 2021)
# Only top 100 candidates are returned by default
SELEXisolate.py config.yml Fkh1_R1.tsv enrichment/Fkh1/k7/affinities_R1_crawled.tsv
SELEXisolate.py config.yml Fkh2_R1.tsv enrichment/Fkh2/k7/affinities_R1_crawled.tsv

# Count how many sequences can be aligned given a specified list of cores
countAlignToCores.py config.yml Fkh1_R1.tsv Fkh1_R1_7mer_isolates.tsv
countAlignToCores.py config.yml Fkh1_R2.tsv Fkh1_R1_7mer_isolates.tsv
countAlignToCores.py config.yml Fkh2_R1.tsv Fkh2_R1_7mer_isolates.tsv
countAlignToCores.py config.yml Fkh2_R2.tsv Fkh2_R1_7mer_isolates.tsv

# Since additional hits taper off around 40 cores, created new file called topisolates.tsv containing the union of top 40 from Fkh1 and Fkh2
# Based on previous knowledge, one core, TGACGCA, was trimmed to the 6-bp Fhl1 motif, GACGCA
# The final list of candidate cores includes 46 unique 7-mers and one 6-mer
# This full list is used to align all reads from each sample
alignToCores.py config.yml Fkh1_R1.tsv topisolates.tsv
alignToCores.py config.yml Fkh1_R2.tsv topisolates.tsv
alignToCores.py config.yml Fkh2_R1.tsv topisolates.tsv
alignToCores.py config.yml Fkh2_R2.tsv topisolates.tsv

# Split aligned reads into independent windows, based on location and orientation of the core
# Creates folder for the output named using the prefix before the first underscore (ie: Fkh1 and Fkh2)
splitShifts.py Fkh1_R1_topisolates.tsv
splitShifts.py Fkh1_R2_topisolates.tsv
splitShifts.py Fkh2_R1_topisolates.tsv
splitShifts.py Fkh2_R2_topisolates.tsv

# Calculate the core and edge ddG/RT for each protein as described in the publication (Cooper et al, 2021)
# Also generates several figures to summarize the calculated values
# Removes any cores that are missing information for any edge position (ie: pseudopalindromes, insufficient read depth)
# If "R2" is in file name, the square root of the enrichment is used in all calculations
# Saved processed data to pkl file ending in _savedScan.pkl depending on input name
SELEXscanAll.py config.yml Fkh1
SELEXscanAll.py config.yml Fkh2