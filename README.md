# Analysis of SELEX-seq Data for Understanding the Binding of FOX TFs
Here we go over several analytical frameworks for the analysis of SELEX-seq data which we obtained for all *S. verecisiae* Forkhead box (FOX) transcription factors (TFs). Steps of the framework can be readily modified for application to other SELEX-seq datasets. 

*Cooper, B. H., Dantas Machado, A. C., Gan, Y., Aparicio, O. M., & Rohs, R. (in revision). DNA Binding Specificity of all four Saccharomyces cerevisiae Forkhead Transcription Factors.*

## Data availability
Processed fastq files can be found on GEO using accession GSE178811, where all trimming and filtering steps are described. 

## Multi-step Alignment
This section details the steps used to align full-length SELEX-seq reads and calculate the relative enrichment and *ΔΔG/RT* of identified cores and flanking positions. This framework depends on the identification of "core" seqeunces, which are *k*-mers which can faithfully indicate true binding sites and which are in alignment with each other to cover a fixed region of the binding site. We utilize cores for alignment rather than a PWM because they do not assume independence between positions of the core and can better account for interdependencies. Analyzing cores separately can also better reveal when "false" cores, which do not faithfully indicate true binding sites, or unaligned cores were included in alignment.

### Identification of candidate cores

We expect cores to be highly enriched within the SELEX-seq data, but not all enriched *k*-mers will cover the same region of the binding site. For this reason, we utilize a previously published tool called Top-Down Crawl (TDC), which was developed for the alignment of quantitative binding data from experiments such as SELEX-seq. 

We start with the alignment of 9-mers, since this is the longest *k*-mer in which a majority of unique *k*-mers occur at least 100 times. Long *k*-mers provide more positions to inform the alignment process, but increasing *k*-mers length decreases the signal-to-noise ratio of enrichment measurements and decreases coverage as described in the text.

*Cooper, B. H., Chiu, T. P., & Rohs, R. (2022). Top-Down Crawl: a method for the ultra-rapid and motif-free alignment of sequences with associated binding metrics. Bioinformatics, 38(22), 5121-5123.*

First, we calculate the relative enrichment of 9-mers using a modified R script from the SELEX package available on bioconductor. The script also depends on the R package, 'stringr'. 

https://bioconductor.org/packages/release/bioc/html/SELEX.html

```
./calculateEnrichment.R <R0 input> <R# input> <round #> <k length>

./calculateEnrichment.R Fkh1_Fkh2_R0.fastq.gz Fkh1_R2.fastq.gz 2 9
# Output: Fkh1_R2_k9.tsv
```

The outputs are then aligned using TDC, which is available for download through pip. 

```
pip install TopDownCrawl
TopDownCrawl Fkh1_R2_k9.tsv
# Output: Fkh1_R2_k9_aligned.tsv
```

Although our alignment includes 9 bp sequences, our final list of cores only need to cover enough bp to faithfully indicate true binding sites. Including additional bp limits the analysis of flanking positions which can no longer be varied across all 4 bases. Furthermore, additional bp increase the noise in flanking *ΔΔG/RT* measurements since they come from longer *k*-mers. For Fkh1, Fkh2, and Hcm1, we trim the 9-mer alignment to the 7 bp region covering the canonical binding site, GTAAACA. For Fhl1, we trim the alignment to the 6 bp region covering the consensus GACGCA. Unique sequences from this process make up the list of candidate cores for each TF.

```
trimToConsensus.py Fkh1_R2_k9_aligned.tsv GTAAACA
trimToConsensus.py Fhl1_R2_k9_aligned.tsv GACGCA
# Output: Fkh1_R2_k9_aligned_GTAAACA.tsv Fhl1_R2_k9_aligned_GACGCA.tsv
```

### Reprioritization of candidate cores
To avoid complications resulting from combinatorial effects between multiple binding sites, we restrict our analysis to reads that only align to one core. This is also key to ensuring that observed flanking preferences are acting to modulate the given core rather than creating additional cores. However, this creates a trade-off between the number of cores we choose to analyze and the number of reads we can align. Therefore, rather than including the entire list of candidate cores in the alignment process, we must prioritize a subset of these sequences. This is done using a framework we call iterative prioritization as described in the manuscript. The included script incorporates the 95% stopping rule as desribed in the text to only return cores which appear to be enriched above background. 

```
prioritize.py config.yml Fkh1_R1.tsv Fkh1_R2_k9_aligned_GTAAACA.tsv
# Output: Fkh1_R1_topcores.tsv
```

In our case, we wanted to use a consistent list of cores accross all TFs tested, so we concatenated the lists of identified cores for each TF and removed duplicates. This can be done with the included helper script as shown below. 

```
mergeCores.py Fkh1_R1_topcores.tsv Fkh2_R1_topcores.tsv Hcm1_R1_topcores.tsv Fhl1_R1_topcores.tsv
# Output: allcores.tsv
```

### Alignment to selected cores and calculation of *ΔΔG/RT*
At this point, full-length reads are aligned to our set of cores, discarding reads which do not align to exactly one core, considering adapter sequences up to 6 bp away upstream or downstream the variable region. Aligned sequences can then be used to calculate the the relative enrichment and *ΔΔG/RT* of differing cores and flanking bp. Several graphics are generated as well as described in the manuscript. 

```
alignToCores.py config.yml Fkh1_Fkh2_R0.tsv allcores.tsv
alignToCores.py config.yml Fkh1_R1.tsv allcores.tsv
alignToCores.py config.yml Fkh1_R2.tsv allcores.tsv

analyzeAlignedSELEX.py config.yml Fkh1_Fkh2_R0_allcores.tsv Fkh1_R1_allcores.tsv 1 Fkh1_R2_allcores.tsv 2
# Output folder: Fkh1_analysis
```

## Analysis of ChIP-exo data



## Multiple Linear Regression (MLR)

## Predicting flanking preferences using DeepBind and MLR

## Caveats
If you would like to apply this framework to your dataset, please consider the following caveats:\
• This framework does not work well for poorly defined cores which do not faithfully indicate true binding sites
• Since this framework requires the alignment of every read to one core (including the reverse complement strand), this framework is not suitable for palindromic or pseudopalindromic cores\
&nbsp;&nbsp;• A palindrome-specific version may be released at a later time that would treat flanking nucleotide contributions as symmetric\