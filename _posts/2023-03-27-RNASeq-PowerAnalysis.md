---
layout: post
title: Power analysis for RNAseq
date: '2023-03-27'
categories: Analysis
tags: [RNASeq, Bioinformatics]
projects: Sediment stress
---

## Conducting a power analysis for RNASeq data 

My sediment stress paper is in review at PeerJ and they need info about power analysis calculation and info on biological & technical replicates. They recommended using the R package [RNASeqPower](https://bioconductor.org/packages/release/bioc/html/RNASeqPower.html) for power analysis/sample size calculation. 

From PeerJ:

Please edit your manuscript to include BOTH:

1. a power analysis calculation (software tools are available, e.g. https://doi.org/doi:10.18129/B9.bioc.RNASeqPower (see usage instructions)) AND
2. information on biological and technical replicates used to achieve the claimed statistical power

 "A suitable wording for the power analysis would be something like: 'The statistical power of this experimental design, calculated in RNASeqPower is 0.84' (replace the program name if you used a different tool and provide the actual calculated power)."

### Power analysis 

What is an RNASeq power analysis? Why do a power analysis?

- An RNASeq power analysis is a statistical method used to estimate the number of samples required to detect differential gene expression with a certain degree of statistical pwoer 
	- Selecting optimal number of biological replicates to achieve desired statistical power (sample size estimation)
	- Estimating the likelihood of successfully finding statistical significance in dataset (power estimation)

From the RNASeqPower vignette: 


	- Depth of sequencing refers to number of reads generated for each sample 
	- Expected count is the average number of reads expected to align to a given transcript in a single sample 
	- cv is a measure of the dispersion of the data
	- Reflects the degree of biological and technical variability within each group
	- Refers to the fold change in expression level between two groups 
	- False positive rate (α) is the probability of incorrectly rejecting the null hypothesis when it is actually true 
	- False negative rate (β) is the probability of incorrectly failing to reject the null hypothesis when its actually false 
	- Power is the probability of correctly rejecting the null when its actually false  

*Note: there are options for n and n2, as well as cv and cv2 if there are two groups. If only one group is specified, the program will assume that the n and cv values are equal for both groups.*

#### Examples 

Let's do some examples first:

```
rnapower(depth = 20, cv = 0.4, effect = c(1.25, 1.5, 1.75, 2), alpha = 0.05, power = c(0.8, 0.9))
```

We see here to detect a 25% increase with high power (0.9), we need ~89 samples. But if we want to detect a 100% change with high power, we only need ~10 samples. As the true biological impact grows, the sample size decreases. Let's increase the depth: 

```
rnapower(depth = 100, cv = 0.4, effect = c(1.25, 1.5, 1.75, 2), alpha = 0.05, power = c(0.8, 0.9))

0.8      0.9
```

Increasing the depth decreases the sample sizes needed. What if we changed the cv?

```
rnapower(depth = 20, cv = 0.3, effect = c(1.25, 1.5, 1.75, 2), alpha = 0.05, power = c(0.8, 0.9))
```

Decreasing the cv also decreases the sample sizes needed. 


#### My data 

In my study, I have 6 species from 2 locations. Does this mean I should do a power analysis for each species? I'll start w/ Acerv from FL. 

First, I'm going to calculate the coefficient of variation of each gene in the filtered gene count matrix.

```
##filtered gene count matrix 

Now I have the CVs for each group for Acerv. It looks like w/ RNASeqPower, there can only be two groups. I'm going to pool the treatment samples. 



Now lets try the code! First, I'm including all the things: 

```


Try running w/o specifying effect size:


This means that, with the given parameters, I will have a 2.6 effect size at 0.8 power. But is this for each gene? Or the mean of the genes for each group?

Let's try running w/o specifying n sizes.

```

This is saying that, for instance, I need a sample size of 9 to achieve an effect size of 2 and a power of 0.8? Does this mean a sample size of 10 for all groups? 

Now let's try running w/o specifying power  

```
rnapower(depth = 5, n = 3, n2 = 10, cv = 0.2569478, cv2 = 0.3105552, effect = c(1.25, 1.5, 1.75, 2), alpha = 0.05)
1.25        1.5       1.75          2 
```

To detect 100% effect (effect of 2), the power will be 0.52. That seems pretty low, given that power is usually set to 0.8 or 0.9.


##### How is depth of coverage calculated? 

- From [Sims et al. 2014](https://www.nature.com/articles/nrg3642): "The average depth of sequencing coverage can be defined theoretically as LN/G, where L is the read length, N is the number of reads and G is the haploid genome length."
	- For Acerv: (50 bp read length * 16,499,909 avg reads) / 430,000,000 bp genome length = 1.92 coverage. What does this mean? Each gene has 1.92 average coverage? That doesn't make any sense, but need to check actually genome bp length 
- From [Hart et al. 2013](https://www.liebertpub.com/doi/10.1089/cmb.2012.0283): "We refer to depth and coverage interchangeably—both meaning how many reads are assigned to a particular gene"
	- This article corresponds to the RNASeqPower R package 
	- "Overall, 85%–96% of all annotated genes were sequenced at a rate of 0.1 reads per million or greater, regardless of the biospecimen or depth. In other words, a sequencing depth of 10 million reads will ensure that approximately 90% of all genes will be covered by at least 10 reads."
- This is the advice that chatGPT gave me (based on [Conesa et al. 2016](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0881-8)): 
	- 1) Determine total number of reads generated by sequencing platform for a single sample 
		- Avg number of reads for Acerv = 16,499,909 
	- 2) Estimate size of transcriptome being analyzed
		- Number of protein coding genes in Acerv = 33,322
	- 3) Divide the total number of reads by the estimated size of the transcriptome to obtain depth of sequencing 
		- 16,499,909 avg reads / 33,322 genes = 495.17 reads per gene? 

Given this calculation, lets calculate n size using depth=495. 

```
rnapower(depth = 495, cv = 0.2569478, cv2 = 0.3105552, effect = c(1.25, 1.5, 1.75, 2), alpha = 0.05, power = c(0.8, 0.9))
		0.8       0.9
```

To achieve an effect size of 2 and a power of 0.8, I would need ~2-3 samples per group?

This [tutorial](https://scienceparkstudygroup.github.io/rna-seq-lesson/02-experimental-design-considerations/index.html#17-power-analysis-for-rna-seq) says that "Replicates are almost always preferred to greater sequencing depth for bulk RNA-Seq" and provides these recommendations for general gene-level differential expression:

- ENCODE guidelines suggest 30 million SE reads per sample (stranded).
- 15 million reads per sample is often sufficient, if there are a good number of replicates (>3).
- Spend money on more biological replicates, if possible.
- Generally recommended to have read length >= 50 bp

### Issues / questions / notes

- Samples in sediment stress study have different read lengths - FL samples are all 50 bp long, but HI samples are either 50 or 125 bp long 
- I have more than 2 groups (in HI, 3 groups; in FL, 5 groups), but RNASeqPower can only handle 2 groups. Do I pool the samples from the treatments?
- What does PeerJ mean by the following statement?
	- "information on biological and technical replicates used to achieve the claimed statistical power" - what information do they mean here?
- Should all genes be held to the same effect size? What subset of the genes (i.e., all genes, filtered genes, differentially expressed genes) should we be considering here?
- If the analysis says I need more replicates, what can I do??
	- From my reading, it seems like biological replicates are more important than sequencing depth

#### Other papers to look at

- [Ching et al. 2014](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4201821/) - Power analysis and sample size estimation for RNA-Seq differential expression
- [Zhao et al. 2018](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2191-5) - RnaSeqSampleSize: real data based sample size estimation for RNA sequencing
	- This is another R package, but uses different parameters (see [vignette](http://www.bioconductor.org/packages/release/bioc/vignettes/RnaSeqSampleSize/inst/doc/RnaSeqSampleSize.pdf))



