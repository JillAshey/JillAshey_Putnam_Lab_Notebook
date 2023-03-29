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

From PeerJ: "A suitable wording for the power analysis would be something like: 'The statistical power of this experimental design, calculated in RNASeqPower is 0.84' (replace the program name if you used a different tool and provide the actual calculated power)."

### Power analysis 

From the vignette: 
A formal sample size calculation for comparison of two groups will involve five factors, each of which is an argument to the rnapower function.
- The depth of sequencing and consequent expected count μ for a given transcript, argument depth.- The coefficient of variation of counts within each of the two groups, argument cv.- The relative expression that we wish to detect ∆, argument effect.- The target false positive rate α and false negative rate β desired (or power = 1 − β), arguments alpha and power.- The number of samples n in each group, argument n

First, I'm going to calculate the coefficient of variation of each gene in the filtered gene count matrix.

```
##filtered gene count matrix acerv_filt <- read.csv("Output/DESeq2/acerv/acerv_counts_sub_filt_20210327.csv")colnames(acerv_filt)[1] <- "gene_id"for ( col in 1:ncol(acerv_filt)){  colnames(acerv_filt)[col] <-  gsub("X", "", colnames(acerv_filt)[col]) # remvoing X from col names}# Separate into Control, T1, T2, T3 & T4 dfsacerv_filt_control <- acerv_filt[,c("gene_id","25_ctl1_Ac_GF_1", "27_ctl2_Ac_YG_1", "41_ctl3_Ac_RN_1")]acerv_filt_T1 <- acerv_filt[,c("gene_id", "37_T13_Ac_ML", "52_T11_Ac_II")] # removed sample 24 from dataset during filtering acerv_filt_T2 <- acerv_filt[,c("gene_id","31_T22_Ac_UV", "38_T23_Ac_IN", "53_T21_Ac_NH")]acerv_filt_T3 <- acerv_filt[,c("gene_id","19_T33_Ac_WK", "47_T31_Ac_JB", "57_T32_Ac_NM")]acerv_filt_T4 <- acerv_filt[,c("gene_id","35_T43_Ac_MT", "54_T42_Ac_JQ")] # removed sample 45 from dataset during filtering # Calculate row mean, std dev, and coefficient of variationacerv_filt_control$mean <- rowMeans(acerv_filt_control[, c(2:4)], na.rm = T)acerv_filt_control$sd <- rowSds(as.matrix(acerv_filt_control[, c(2:4)], na.rm = T))acerv_filt_control$cv <- rowCVs(acerv_filt_control[,2:4])mean(acerv_filt_control$cv, na.rm = T) # 0.2569478median(acerv_filt_control$cv, na.rm = T) # 0.2067829acerv_filt_T1$mean <- rowMeans(acerv_filt_T1[, c(2:3)], na.rm = T)acerv_filt_T1$sd <- rowSds(as.matrix(acerv_filt_T1[, c(2:3)], na.rm = T))acerv_filt_T1$cv <- rowCVs(acerv_filt_T1[,2:3])mean(acerv_filt_T1$cv, na.rm = T) # 0.2148519median(acerv_filt_T1$cv, na.rm = T) # 0.1741535acerv_filt_T2$mean <- rowMeans(acerv_filt_T2[, c(2:4)], na.rm = T)acerv_filt_T2$sd <- rowSds(as.matrix(acerv_filt_T2[, c(2:4)], na.rm = T))acerv_filt_T2$cv <- rowCVs(acerv_filt_T2[,2:4])mean(acerv_filt_T2$cv, na.rm = T) # 0.3396199median(acerv_filt_T2$cv, na.rm = T) # 0.2994043acerv_filt_T3$mean <- rowMeans(acerv_filt_T3[, c(2:4)], na.rm = T)acerv_filt_T3$sd <- rowSds(as.matrix(acerv_filt_T3[, c(2:4)], na.rm = T))acerv_filt_T3$cv <- rowCVs(acerv_filt_T3[,2:4])mean(acerv_filt_T3$cv, na.rm = T) # 0.5463428median(acerv_filt_T3$cv, na.rm = T) # 0.5210751acerv_unfilt_T4$mean <- rowMeans(acerv_unfilt_T4[, c(2:3)], na.rm = T)acerv_unfilt_T4$sd <- rowSds(as.matrix(acerv_unfilt_T4[, c(2:3)], na.rm = T))acerv_unfilt_T4$cv <- rowCVs(acerv_unfilt_T4[,2:3])mean(acerv_unfilt_T4$cv, na.rm = T) # 0.5909885median(acerv_filt_T4$cv, na.rm = T) # Null ? # Average cv for T1, T2, T3, T4df <- c(0.2148519, 0.3396199, 0.5463428, 0.5909885)mean(df) # 0.4229508 ```

Now I have the CVs for each group for Acerv. It looks like w/ RNASeqPower, there can only be two groups. I'm going to pool the treatment samples. 
Parameters: 
- Depth = 5 (?)- n = 3 (control)- n2 = 10 (T1, T2, T3, T4 - not including outliers)- cv = 0.2569478- cv2 = 0.4229508- effect = 1.25, 1.5, 1.75, 2- alpha = 0.05- power = 0.8, 0.9

Now lets try the code!
```rnapower(depth = 5, n = 3, n2 = 10, cv = 0.2569478, cv2 = 0.4229508, effect = c(1.25, 1.5, 1.75, 2), alpha = 0.05, power = c(0.8, 0.9))## ^^ When I try to run the code above with all the info, I get this error: Exactly one of n, cv, effect, alpha, or power should be unspecified. So it looks like I need to not include one of the arguments so the function can calculate it for me. But what am I trying to calculate if I have all the information already?### try to # w/o effect size rnapower(depth = 5, n = 3, n2 = 10, cv = 0.2569478, cv2 = 0.4229508, alpha = 0.05, power = c(0.8, 0.9))#     0.8      0.9 # 2.709279 3.168286 # What does this mean?# w/o nrnapower(depth = 5, cv = 0.2569478, cv2 = 0.4229508, effect = c(1.25, 1.5, 1.75, 2), alpha = 0.05, power = c(0.8, 0.9))#            0.8       0.9# 1.25 101.65712 136.09004# 1.5   30.78928  41.21811# 1.75  16.16317  21.63790# 2     10.53551  14.10406# this is saying that, for instance, I need a sample size of 10 to achieve an effect size of 2 and a power of 0.8?```

How is depth of coverage calculated? 

- From [Sims et al. 2014](https://www.nature.com/articles/nrg3642): "The average depth of sequencing coverage can be defined theoretically as LN/G, where L is the read length, N is the number of reads and G is the haploid genome length."
	- For Acerv: (50 bp read length * 15,000,000 reads) / 430,000,000 bp genome length = 1.74 coverage. What does this mean? Each gene has 1.74 average coverage. 
- From [Hart et al. 2013](https://www.liebertpub.com/doi/10.1089/cmb.2012.0283): "We refer to depth and coverage interchangeably—both meaning how many reads are assigned to a particular gene"
	- This article corresponds to the RNASeqPower R package 
	- "Overall, 85%–96% of all annotated genes were sequenced at a rate of 0.1 reads per million or greater, regardless of the biospecimen or depth. In other words, a sequencing depth of 10 million reads will ensure that approximately 90% of all genes will be covered by at least 10 reads."

This [tutorial](https://scienceparkstudygroup.github.io/rna-seq-lesson/02-experimental-design-considerations/index.html#17-power-analysis-for-rna-seq) says that "Replicates are almost always preferred to greater sequencing depth for bulk RNA-Seq" and provides these recommendations for general gene-level differential expression:

- ENCODE guidelines suggest 30 million SE reads per sample (stranded).
- 15 million reads per sample is often sufficient, if there are a good number of replicates (>3).
- Spend money on more biological replicates, if possible.
- Generally recommended to have read length >= 50 bp



### Issues 

- Samples in sediment stress study have different read lengths - FL samples are all 50 bp long, but HI samples are either 50 or 125 bp long 
- Should all genes be held to the same effect size? What subset of the genes (i.e., all genes, filtered genes, differentially expressed genes) should we be considering here?
- If the analysis says I need more replicates, what can I do? 
	- From my reading, it seems like biological replicates are more important than sequencing depth

#### Other papers to look at

- [Ching et al. 2014](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4201821/) - Power analysis and sample size estimation for RNA-Seq differential expression
- [Zhao et al. 2018](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2191-5) - RnaSeqSampleSize: real data based sample size estimation for RNA sequencing
	- This is another R package, but uses different parameters (see [vignette](http://www.bioconductor.org/packages/release/bioc/vignettes/RnaSeqSampleSize/inst/doc/RnaSeqSampleSize.pdf))




