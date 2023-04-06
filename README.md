# Analysis of SELEX-seq Data for Understanding the Binding of FOX TFs
Here we go over several analytical frameworks for the analysis of SELEX-seq data which we obtained for all *S. verecisiae* Forkhead box (FOX) transcription factors (TFs), as described in our manuscript<sup>1</sup>. The framework can be readily modified for application to other SELEX-seq datasets. 

<sup>1</sup> *Cooper, B. H., Dantas Machado, A. C., Gan, Y., Aparicio, O. M., & Rohs, R. (in revision). DNA Binding Specificity of all four Saccharomyces cerevisiae Forkhead Transcription Factors.*

## Setup / Installation
It is best to install all required packaged into an isolated environment in order to avoid dependency conflicts. I recommend doing this using conda/miniconda (https://anaconda.org/). 

```
conda create -c conda-forge -c bioconda -c defaults -n multi-step-align python=3.8 bedtools bioconductor-selex r-stringr r-rjava
conda activate multi-step-align

pip install .

cd keras_genomics
pip install .
cd ..
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
countUnique.py Fkh1_R1.fastq.gz
countUnique.py Fkh1_R2.fastq.gz
countUnique.py Fkh1_Fkh2_R0.fastq.gz
# Output: Fkh1_R1.tsv, Fkh1_R2.tsv, Fkh1_Fkh2_R0.tsv

split.py Fkh1_Fkh2_R0.tsv
# Output: Fkh1_Fkh2_R0_split1.tsv, Fkh1_Fkh2_R0_split2.tsv

# Help: calculateEnrichment.py -h
calculateEnrichment.py Fkh1_Fkh2_R0_split1.tsv Fkh1_Fkh2_R0_split2.tsv Fkh1_R2.tsv 2 9
# Output: Fkh1_R2_k9.tsv
# Can be repeated for Fkh2, Hcm1, and Fhl1
```

The outputs are then aligned using TDC<sup>2</sup>, which is available for download through pip. 

```
TopDownCrawl Fkh1_R2_k9.tsv
# Output: Fkh1_R2_k9_aligned.tsv
# Can be repeated for Fkh2, Hcm1, and Fhl1
```
 
Although our alignment includes 9 bp sequences, our final list of cores only need to cover enough bp to faithfully indicate true binding sites. Including additional bp limits the analysis of flanking positions which can no longer be varied across all 4 bases. Furthermore, additional bp increase the noise in flanking *ΔΔG/RT* measurements since they come from longer *k*-mers. For Fkh1, Fkh2, and Hcm1, we trim the 9-mer alignment to the 7 bp region covering the most enriched sequence containing the canonical binding site, GTAAACA, or GACGCA for Fhl1. Unique sequences from this process make up the list of candidate cores for each TF.

```
# Help: trimToConsensus.py -h
trimToConsensus.py Fkh1_R2_k9_aligned.tsv GTAAACA
trimToConsensus.py Fhl1_R2_k9_aligned.tsv GACGCA
# Output: Fkh1_R2_k9_aligned_GTAAACA.tsv Fhl1_R2_k9_aligned_GACGCA.tsv
```

### Reprioritization of candidate cores
To avoid complications resulting from combinatorial effects between multiple binding sites, we restrict our analysis to reads that only align to one core. This is also key to ensuring that observed flanking preferences are acting to modulate the given core rather than creating additional cores. However, this creates a trade-off between the number of cores we choose to analyze and the number of reads we can align. Therefore, rather than including the entire list of candidate cores in the alignment process, we must prioritize a subset of these sequences. This is done using a framework we call iterative prioritization as described in the manuscript<sup>1</sup>. The included script incorporates the 95% stopping rule as desribed in the text to only return cores which appear to be enriched above background. 

```
# Help: prioritize.py -h
prioritize.py Fkh1_Fkh2_R0.tsv Fkh1_R1.tsv 1 Fkh1_R2_k9_aligned_GTAAACA.tsv GAGTTCTACAGTCCGACGATCCAG TCCGTATCGCTCCTCCAATG 0.95
prioritize.py Fkh1_Fkh2_R0.tsv Fkh2_R1.tsv 1 Fkh2_R2_k9_aligned_GTAAACA.tsv GAGTTCTACAGTCCGACGATCCAG TCCGTATCGCTCCTCCAATG 0.95
prioritize.py Hcm1_Fhl1_R0.tsv Hcm1_R1.tsv 1 Hcm1_R2_k9_aligned_GTAAACA.tsv GAGTTCTACAGTCCGACGATCCAG TCCGTATCGCTCCTCCAATG 0.95
prioritize.py Hcm1_Fhl1_R0.tsv Fhl1_R1.tsv 1 Fhl1_R2_k9_aligned_GACGCA.tsv GAGTTCTACAGTCCGACGATCCAG TCCGTATCGCTCCTCCAATG 0.95
# Output: Fkh1_R1_topcores.tsv
```

In our case, we wanted to use a consistent list of cores accross all TFs tested, so we concatenated the lists of identified cores for each TF and removed duplicates. This can be done with the included helper script as shown below. 

```
# Help: mergeCores.py -h
mergeCores.py *_topcores.tsv
# Output: allcores.tsv
```

### Alignment to selected cores and calculation of *ΔΔG/RT*
At this point, full-length reads are aligned to our set of cores, discarding reads which do not align to exactly one core, considering adapter sequences up to 6 bp away upstream or downstream the variable region. Aligned sequences can then be used to calculate the the relative enrichment and *ΔΔG/RT* of differing cores and flanking bp. Several graphics are generated as well as described in the manuscript<sup>1</sup>. 

```
# Help: alignToCores.py -h
alignToCores.py Fkh1_Fkh2_R0.tsv allcores.tsv GAGTTCTACAGTCCGACGATCCAG TCCGTATCGCTCCTCCAATG
alignToCores.py Fkh1_R1.tsv allcores.tsv GAGTTCTACAGTCCGACGATCCAG TCCGTATCGCTCCTCCAATG
`alignToCores.py Fkh1_R2.tsv allcores.tsv GAGTTCTACAGTCCGACGATCCAG TCCGTATCGCTCCTCCAATG

# Help: analyzeAlignedSELEX.py -h
# Bound was set to ensure a consistent color scale accross datasets
analyzeAlignedSELEX.py 7 Fkh1_Fkh2_R0_allcores.tsv Fkh1_R1_allcores.tsv 1 Fkh1_R2_allcores.tsv 2 --autoscale --bound 0.75
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
# Fkh1_core_RelE.tsv in Fkh1_analysis folder, we focus only on 7-mers for this analysis
groupByLength.py Fkh1_core_RelE.tsv
# Output: Fkh1_core_RelE_k7.tsv
countHits.py Fkh1_combined_100_merged.fa Fkh1_core_RelE_k7.tsv
countHits.py sacCer3.fa Fkh1_core_RelE_k7.tsv
# Output: Fkh1_combined_100_merged_hits.tsv, sacCer3_hits.tsv

predictObservationsWithSELEX.py Fkh1_combined_100_merged_hits.tsv sacCer3_hits.tsv Fkh1_core_RelE_k7.tsv
# Output: predSELEX.png
```

Alternatively, we can predict the expected enrichment using predicted frequencies according to a BEESEM-derived PFM. The orignal BEESEM<sup>4</sup> GitHub page is linked below, however we provide a forked repository with minimal changes to allow execution of beesem.py and formatter.py after installation with setup.py, and to allow generation of motifs up to 13 bp long.

https://github.com/sx-ruan/BEESEM
https://github.com/bhcooper/BEESEM

<sup>4</sup> *Ruan, S., Swamidass, S. J., & Stormo, G. D. (2017). BEESEM: estimation of binding energy models using HT-SELEX data. Bioinformatics, 33(15), 2288-2295.*

Additionally, we create a new environment for BEESEM since it requires Python 2.7. 

``` 
conda create -c conda-forge -c bioconda -c defaults -n BEESEM python=2.7 seqtk
conda activate BEESEM
pip install numpy scipy matplotlib
pip install openopt funcdesigner
git clone https://github.com/bhcooper/BEESEM.git
cd BEESEM
pip install .
cd ..

# We were not able to run BEESEM with our full set of reads due to the extreme computational demands, so we used a random subset of  10% of our reads
seqtk sample Fkh1_Fkh2_R0.fastq.gz 0.1 | gzip > Fkh1_Fkh2_R0_0.1.fastq.gz
seqtk sample Fkh1_R1.fastq.gz 0.1 | gzip > Fkh1_R1_0.1.fastq.gz

conda activate multi-step-align
countUnique.py Fkh1_Fkh2_R0_0.1.fastq.gz
countUnique.py Fkh1_R1_0.1.fastq.gz

# Remove headers for BEESEM
tail -n +2 Fkh1_Fkh2_R0_0.1.tsv > Fkh1_Fkh2_R0_0.1_BEESEM.tsv
tail -n +2 Fkh1_R1_0.1.tsv > Fkh1_R1_0.1_BEESEM.tsv

conda activate BEESEM
# beesem.py -s <most enriched 7-mer as seed> -f <ladapter> <radapter> <output name> <processed prior> <processed R# input>
beesem.py -s GTAAACA -f GAGTTCTACAGTCCGACGATCCAG TCCGTATCGCTCCTCCAATG beesem_GTAAACA Fkh1_Fkh2_R0_0.1_BEESEM.tsv Fkh1_R1_0.1_BEESEM.tsv
# Output: beesem_k7_rep=1_phs=10/results/pfm_r0_p9_w7.txt

conda activate multi-step-align
predictObservationsWithBEESEM.py Fkh1_combined_100_merged_hits.tsv sacCer3_hits.tsv beesem_GTAAACA_rep=1_phs=10/results/pfm_r0_p9_w7.txt
# Output: predBEESEM.png
```

We also investigated whether the *in vivo* enrichment of bp flanking the core could be estimated using our SELEX-seq data. For this purpose, we focused on the positions flanking the most enriched core GTAAACA. This script performs both the alignment-based and BEESEM-based comparisons. With limited memory, BEESEM can be used to to calculate 5' and 3' preferences separately, then preferenced from each can be used. 

```
conda activate BEESEM
beesem.py -s AAAAGTAAACA -f GAGTTCTACAGTCCGACGATCCAG TCCGTATCGCTCCTCCAATG beesem_AAAAGTAAACA Fkh1_Fkh2_R0_0.1_BEESEM.tsv Fkh1_R1_0.1_BEESEM.tsv
beesem.py -s GTAAACAAA -f GAGTTCTACAGTCCGACGATCCAG TCCGTATCGCTCCTCCAATG beesem_GTAAACAAA Fkh1_Fkh2_R0_0.1_BEESEM.tsv Fkh1_R1_0.1_BEESEM.tsv

conda activate multi-step-align
countFlanks.py Fkh1_combined_100_merged.fa GTAAACA 4 2
countFlanks.py sacCer3.fa GTAAACA 4 2
# Output: Fkh1_combined_100_merged_GTAAACA_-4_+2.tsv, sacCer3_GTAAACA_-4_+2.tsv

# Fkh1_flank_ddG.tsv from Fkh1_analysis folder
# Help: predictFlanks.py -h
predictFlanks.py Fkh1_combined_100_merged_GTAAACA_-4_+2.tsv sacCer3_GTAAACA_-4_+2.tsv Fkh1_flank_ddG.tsv GTAAACA \
    beesem_AAAAGTAAACA_rep=1_phs=10/results/pfm_r0_p9_w11.txt 4 beesem_GTAAACAAA_rep=1_phs=10/results/pfm_r0_p9_w9.txt 2
# Output: predFlanks_GTAAACA.png, performance metrics are printed to STDOUT
```

If you have a computer with a large amount of memory (≈128 GB), you can use BEESEM to generate a 13-bp motif covering all flanking positions of interest. 

```
conda activate BEESEM
beesem.py -s AAAAGTAAACAAA -f GAGTTCTACAGTCCGACGATCCAG TCCGTATCGCTCCTCCAATG beesem_AAAAGTAAACAAA Fkh1_Fkh2_R0_0.1_BEESEM.tsv Fkh1_R1_0.1_BEESEM.tsv

conda activate multi-step-align
predictFlanks.py Fkh1_combined_100_merged_GTAAACA_-4_+2.tsv sacCer3_GTAAACA_-4_+2.tsv Fkh1_flank_ddG.tsv GTAAACA \
    beesem_AAAAGTAAACAAA_rep=1_phs=10/results/pfm_r0_p9_w13.txt 4 beesem_AAAAGTAAACAAA_rep=1_phs=10/results/pfm_r0_p9_w13.txt 2
```

In our case, SELEX-based predictions were able to better predict the observed enrichment of core sequences and flanking bp compared to BEESEM-based predictions. Additional details are discussed in the text<sup>1</sup>. 

## Multiple Linear Regression (MLR)

The goal of this section is to use MLR to probe the importance of interdependencies accross differing regions of the binding site. By using this simple model, we can easily add or withold higher order features, as described in the text<sup>1</sup>, to determine what extent these features impact performance. 

```
# Help: getMaskedEnrichment.py -h
# getMaskedEnrichment.py <R0 input> <input alignment> <round #> <max core length> <5' flank length> <3' flank length>
getMaskedEnrichment.py Fkh1_Fkh2_R0.tsv Fkh1_R2_allcores.tsv 2 7 4 2
# Output: Fkh1_R2_allcores_-4_+2.tsv

# scoreMLR.py <enrichment table> <5' flank length> <3' flank length>
scoreMLR.py Fkh1_R2_allcores_-4_+2.tsv 4 2
# Output: Fkh1_R2_allcores_-4_+2_MLR.png, Fkh1_R2_allcores_-4_+2_MLR.tsv
```

We can also evaluate the performance of models predicting *ΔΔG/RT* of shorter *k*-mers, which cannot capture as many interependencies, but are generally less noisy.

```
getMaskedEnrichment.py Fkh1_Fkh2_R0.tsv Fkh1_R2_allcores.tsv 2 7 4 0
# Output: Fkh1_R2_allcores_-4_+0.tsv

scoreMLR.py Fkh1_R2_allcores_-4_+0.tsv 4 0
# Output: Fkh1_R2_allcores_-4_+0_MLR.png, Fkh1_R2_allcores_-4_+0_MLR.tsv
```

## Predicting flanking preferences using DeepBind and MLR

To compare our alignment-based framework with a deep learning approach, we trained a model based on DeepBind's reverse-complement parameter sharing framework<sup>5,6</sup> to predict high-resolution estimates of the *ΔΔG/RT* for any given 13-mer. To interpret the model in a matter that is comparable to our approach, we first predicted the *ΔΔG/RT* of all 13-mers covering four bp 5' and two bp 3' each of our selected seven bp cores. We then used these predictions to train an MLR model for each core, using 1-mer sequence features of the flanks as input. For each flanking position, model weights are centered and plotted for comparison with our flanking *ΔΔG/RT* measurements.

<sup>5</sup> *Alipanahi, B., Delong, A., Weirauch, M. T., & Frey, B. J. (2015). Predicting the sequence specificities of DNA-and RNA-binding proteins by deep learning. Nature biotechnology, 33(8), 831-838.*
<sup>6</sup> *Shrikumar, A., Greenside, P., & Kundaje, A. (2017). Reverse-complement parameter sharing improves deep learning models for genomics. BioRxiv, 103663.*

```
# Removing the 100 count minimum greatly improves predictive power, feel free to evaluate the filtered dataset by using calculateEnrichment.R instead
calculateEnrichmentNoMin.py Fkh1_Fkh2_R0_split1.tsv Fkh1_Fkh2_R0_split2.tsv Fkh1_R2.tsv 2 13
# Output: Fkh1_R2_k13_noMin.tsv

trainDeepBind.py Fkh1_R2_k13_noMin.tsv
# Output: trained_cnn.h5

# From Fkh1_analysis folder
groupByLength.py Fkh1_core_ddG.tsv
# Output: Fkh1_core_ddG_k7.tsv

# Help: predDeepBind.py -h
predDeepBind.py trained_cnn.h5 Fkh1_core_ddG_k7.tsv Fkh1_flank_ddG.tsv 4 2 --bound 0.75
# Output: DeepBind_logos, DeepBind_flank_matrix.png, alignment_flank_matrix.png
```

We found that our framework was more sensitive to flanking preferences, especially for lower affinity cores. Furthermore, our model was better able to distinguish preferences for flanking bp which modulate the affinity of the selected core vs flanking bp which create additional cores. 

## Caveats
If you would like to apply the multi-step alignment framework to your dataset, please consider the following caveats:\
• This framework does not work well for poorly defined cores which do not faithfully indicate true binding sites\
• Since this framework requires the alignment of every read to one core (including the reverse complement strand), this framework is not suitable for the analysis of palindromic or pseudopalindromic cores\
&nbsp;&nbsp;&nbsp;&nbsp;• A palindrome-specific version may be released at a later time that would treat flanking nucleotide contributions as symmetric\