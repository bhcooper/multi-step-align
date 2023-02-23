# Analysis of SELEX-seq Data for Understanding the Binding of FOX TFs
Here we go over several analytical frameworks for the analysis of SELEX-seq data which we obtained for all *S. verecisiae* Forkhead box (FOX) transcription factors (TFs), as described in our manuscript<sup>1</sup>. The framework can be readily modified for application to other SELEX-seq datasets. 

<sup>1</sup> *Cooper, B. H., Dantas Machado, A. C., Gan, Y., Aparicio, O. M., & Rohs, R. (in revision). DNA Binding Specificity of all four Saccharomyces cerevisiae Forkhead Transcription Factors.*

## Setup / Installation
It is best to install all required packaged into an isolated environment in order to avoid dependency conflicts. I recommend doing this using conda (https://anaconda.org/). 

```
conda create -c conda-forge -c bioconda -c defaults -n multi-step-align python=3.7 bedtools bioconductor-selex r-stringr r-rjava
conda activate multi-step-align

pip install .
```

## Data availability
Processed fastq files can be found on GEO using accession GSE178811, where all trimming and filtering steps are described. 

## Multi-step Alignment
This section details the steps used to align full-length SELEX-seq reads and calculate the relative enrichment and *ΔΔG/RT* of identified cores and flanking positions. This framework depends on the identification of "core" seqeunces, which are *k*-mers which can faithfully indicate true binding sites and which are in alignment with each other to cover a fixed region of the binding site. We utilize cores for alignment rather than a PWM because they do not assume independence between positions of the core and can better account for interdependencies. Analyzing cores separately can also better reveal when "false" cores, which do not faithfully indicate true binding sites, or unaligned cores were included in alignment.

### Identification of candidate cores

We expect cores to be highly enriched within the SELEX-seq data, but not all enriched *k*-mers will cover the same region of the binding site. For this reason, we utilize a previously published tool called Top-Down Crawl<sup>2</sup> (TDC), which was developed for the alignment of quantitative binding data from experiments such as SELEX-seq. 


We start with the alignment of 9-mers, since this is the longest *k*-mer in which a majority of unique *k*-mers occur at least 100 times. Long *k*-mers provide more positions to inform the alignment process, but increasing *k*-mers length decreases the signal-to-noise ratio of enrichment measurements and decreases coverage as described in the text<sup>1</sup>.

<sup>2</sup> *Cooper, B. H., Chiu, T. P., & Rohs, R. (2022). Top-Down Crawl: a method for the ultra-rapid and motif-free alignment of sequences with associated binding metrics. Bioinformatics, 38(22), 5121-5123.*

First, we calculate the relative enrichment of 9-mers using a modified R script from the SELEX package available on bioconductor. The script also depends on the R package, 'stringr'. 

https://bioconductor.org/packages/release/bioc/html/SELEX.html

```
calculateEnrichment.R <R0 input> <R# input> <round #> <k length>

calculateEnrichment.R Fkh1_Fkh2_R0.fastq.gz Fkh1_R2.fastq.gz 2 9
# Output: Fkh1_R2_k9.tsv
```

The outputs are then aligned using TDC<sup>2</sup>, which is available for download through pip. 

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
To avoid complications resulting from combinatorial effects between multiple binding sites, we restrict our analysis to reads that only align to one core. This is also key to ensuring that observed flanking preferences are acting to modulate the given core rather than creating additional cores. However, this creates a trade-off between the number of cores we choose to analyze and the number of reads we can align. Therefore, rather than including the entire list of candidate cores in the alignment process, we must prioritize a subset of these sequences. This is done using a framework we call iterative prioritization as described in the manuscript<sup>1</sup>. The included script incorporates the 95% stopping rule as desribed in the text to only return cores which appear to be enriched above background. 

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
At this point, full-length reads are aligned to our set of cores, discarding reads which do not align to exactly one core, considering adapter sequences up to 6 bp away upstream or downstream the variable region. Aligned sequences can then be used to calculate the the relative enrichment and *ΔΔG/RT* of differing cores and flanking bp. Several graphics are generated as well as described in the manuscript<sup>1</sup>. 

```
alignToCores.py config.yml Fkh1_Fkh2_R0.tsv allcores.tsv
alignToCores.py config.yml Fkh1_R1.tsv allcores.tsv
alignToCores.py config.yml Fkh1_R2.tsv allcores.tsv

analyzeAlignedSELEX.py config.yml Fkh1_Fkh2_R0_allcores.tsv Fkh1_R1_allcores.tsv 1 Fkh1_R2_allcores.tsv 2
# Output folder: Fkh1_analysis
```

## Analysis of ChIP-exo data

The goal of this analysis is to determine if our SELEX-seq based measurements of enrichment can be used to predict the enrichment of sequences bound *in vivo*. In a previously published work, researchers performed ChIP-exo experiments targeting Fkh1 in S. cerevisiae<sup>3</sup>. In the corresponding GitHub page, a fasta file is provided by the authors containing genomic sequences spanning ±250 bp around each peak identified in the study. In our case, we were interested in counting the frequency of bound sites rather than identifying de novo motifs. For this reason, we aimed to utilize a narrower window spanning ±50 around each peak. Furthermore, to avoid double counting bound sites, we then merged overlapping peaks.

<sup>3</sup> *Mondeel, T. D. G. A., Holland, P., Nielsen, J., & Barberis, M. (2019). ChIP-exo analysis highlights Fkh1 and Fkh2 transcription factors as hubs that integrate multi-scale networks in budding yeast. Nucleic acids research, 47(15), 7825-7841.*

```
wget https://raw.githubusercontent.com/barberislab/ChIP-exo_Fkh1_Fkh2/master/Data/binding_motif_analysis/fkh1_combined/Fkh1_combined.fasta


getBed100.py Fkh1_combined.fasta
Output: Fkh1_combined_100.bed

wget https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.fa.gz
gunzip sacCer3.fa.gz
bedtools merge -i Fkh1_combined_100.bed > Fkh1_combined_100_merged.bed 
bedtools getfasta -fi sacCer3.fa -bed Fkh1_combined_100_merged.bed > Fkh1_combined_100_merged.fa
```

From this output, we can calculate the relative frequency of cores which appear to have been bound by Fkh1 *in vivo*. We can then compare these values to those that would be predicted based on the relative frequency of each core accross the entire genome and the enrichment measurements obtained from the SELEX-seq experiment. 

```
# Fkh1_core_RelE.tsv in Fkh1_analysis folder
countHits.py Fkh1_combined_100_merged.fa Fkh1_core_RelE.tsv
countHits.py sacCer3.fa Fkh1_core_RelE.tsv
# Output: Fkh1_combined_100_merged_hits.tsv, sacCer3_hits.tsv

predictObservationsWithSELEX.py Fkh1_combined_100_merged_hits.tsv sacCer3_hits.tsv Fkh1_core_RelE.tsv
# Output: predSELEX.png
```

Alternatively, we can predict the expected enrichment using predicted frequencies according to a BEESEM-derived PFM. BEESEM can be downloaded from the GitHub page linked below. 

https://github.com/sx-ruan/BEESEM

``` 
beesem.py -s <most enriched 7-mer as seed> <output name> <processed prior> <processed R# input>
beesem.py -s GTAAACA beesem_k7 Fkh1_Fkh2_R0.tsv Fkh1_R1.tsv
# Output: beesem_k7_rep=1_phs=10/results/pfm_r0_p9_w7.txt

predictObservationsWithBEESEM.py Fkh1_combined_100_merged_hits.tsv sacCer3_hits.tsv beesem_k7_rep=1_phs=10/results/pfm_r0_p9_w7.txt
# Output: predBEESEM.png
```

We also investigated whether the *in vivo* enrichment of bp flanking the core could be estimated using our SELEX-seq data. For this purpose, we focused on the positions flanking the most enriched core GTAAACA.

```
countFlanks.py Fkh1_combined_100_merged.fa GTAAACA 4 2
countFlanks.py sacCer3.fa GTAAACA 4 2
# Output: Fkh1_combined_100_merged_GTAAACA_-4_+2.tsv, sacCer3_GTAAACA_-4_+2.tsv

# Fkh1_edge_ddG.tsv in Fkh1_analysis folder
predictFlanks.py Fkh1_combined_100_merged_GTAAACA_-4_+2.tsv sacCer3_GTAAACA_-4_+2.tsv Fkh1_edge_ddG.tsv GTAAACA BEESEM_flanks.txt
# Output: predFlanks_GTAAACA.png, metrics printed to STDOUT
```

In our case, SELEX-based predictions were able to better predict the observed enrichment of core sequences and flanking bp compared to BEESEM-based predictions. Additional details are discussed in the text<sup>1</sup>. 

## Multiple Linear Regression (MLR)

The goal of this section is to use MLR to probe the importance of interdependencies accross differing regions of the binding site. By using this simple model, we can easily add or withold higher order features, as described in the text<sup>1</sup>, to determine what extent these features impact performance. 

```
getMaskedEnrichment.py <input alignment> <round #> <5' flank length> <3' flank length>
getMaskedEnrichment.py Fkh1_R2_allcores.tsv 2 4 2
# Output: Fkh1_R2_allcores_-4_+2.tsv

getMaskedEnrichment.py <input> <5' flank length> <3' flank length>
scoreMLR.py Fkh1_R2_allcores_-4_+2.tsv 4 2
# Output: Fkh1_R2_allcores_-4_+2_MLR.png
```

We can also evaluate the performance of models predicting *ΔΔG/RT* of shorter *k*-mers, which cannot capture as many interependencies, but are generally less noisy.

```
getMaskedEnrichment.py <input alignment> <round #> <5' flank length> <3' flank length>
getMaskedEnrichment.py Fkh1_R2_allcores.tsv 2 4 0
# Output: Fkh1_R2_allcores_-4_+0.tsv

getMaskedEnrichment.py <input> <5' flank length> <3' flank length>
scoreMLR.py Fkh1_R2_allcores_-4_+0.tsv 4 0
# Output: Fkh1_R2_allcores_-4_+0_MLR.png
```

## Predicting flanking preferences using DeepBind and MLR

To compare our alignment-based framework with a deep learning approach, we trained a model based on DeepBind's reverse-complement weight sharing framework<sup>4</sup> to predict high-resolution estimates of the *ΔΔG/RT* for any given 13-mer. To interpret the model in a matter that is comparable to our approach, we first predicted the *ΔΔG/RT* of all 13-mers covering four bp 5' and two bp 3' each of our selected seven bp cores. We then used these predictions to train an MLR model for each core, using 1-mer sequence features of the flanks as input. For each flanking position, model weights are centered and plotted for comparison with our flanking *ΔΔG/RT* measurements.

```
# I used conda to ensure compatible versioning with DeepBind
conda create -n deepbind tensorflow=1 h5py=2.10 pandas matplotlib logomaker scikit-learn
conda activate deepbind

git clone --branch keras_1 https://github.com/kundajelab/keras.git
cd keras
pip install .
cd ..

# Removing the 100 count minimum seemed to greatly improve predictive power in my case
calculateEnrichmentNoMin.R Fkh1_Fkh2_R0.fastq.gz Fkh1_R2.fastq.gz 2 13
# Output: Fkh1_R2_k13_noMin.tsv

trainDeepBind.py Fkh1_R2_k13_noMin.tsv
# Output: trained_cnn.h5, trained_minMax.pkl

predDeepBind.py trained_cnn.h5 traind_minMax.pkl Fkh1_core_ddG.tsv Fkh1_edge_ddG.tsv
# Output: DeepBind_MLR.png
```

We found that our framework was more sensitive to flanking preferences, especially for lower affinity cores. Furthermore, our model was better able to distinguish preferences for flanking bp which modulate the affinity of the selected core vs flanking bp which create additional cores. 

## Caveats
If you would like to apply the multi-step alignment framework to your dataset, please consider the following caveats:\
• This framework does not work well for poorly defined cores which do not faithfully indicate true binding sites\
• Since this framework requires the alignment of every read to one core (including the reverse complement strand), this framework is not suitable for the analysis of palindromic or pseudopalindromic cores\
&nbsp;&nbsp;&nbsp;&nbsp;• A palindrome-specific version may be released at a later time that would treat flanking nucleotide contributions as symmetric\