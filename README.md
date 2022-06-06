# multi-step-align

Following the workflow.sh script, you should be able to repeat all analysis steps described in Cooper et al (to be published). Data files can be downloaded from GEO with accession GSE178811.\
\
All scripts were tested using Python 3.8.8 on linux. Scripts must be made executable (chmod +x <filename>) and can be added to your PATH variable for easy execution. See requirements.txt for required packages and version numbers used during testing.\
 \
Note: If you would like to apply this framework to your dataset, please consider the following caveats:\
• This framework was designed to analyze round 1 (R1) or round 2 (R2) SELEX-seq data since these rounds will have more information about moderate to low affinity binding sites.\
• Since this framework requires the alignment of every read to one core (including the reverse complement strand), this framework is not suitable for palindromic or pseudopalindromic cores\
  • A palindrome-specific version may be released at a later time that would treat flanking nucleotide contributions as symmetric\
• This framework does not work well for poorly defined cores, or cores that do not faithfully indicate a binding site (ie: very short cores <5 bp)
