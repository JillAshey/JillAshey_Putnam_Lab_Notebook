---
layout: post
title: Mcap Ribofree expression test
date: '2024-01-31'
categories: Analysis
tags: [Bioinformatics, Mcap]
projects: Mcap Developmental Timeseries 2023
---

I am deciding what kind of depth/coverage I will need for sequencing if I proceed with the Zymo Ribofree library prep kit. I am using the gene count matrices generated through the [Express Compare](https://github.com/hputnam/Express_Compare/tree/main) project. In this post, I looked at all of the Mcap Ribofree gene count matrices that I could find to assess how many genes were retained before and after filtering. I did a two-step filtering process - first, I removed all rows whose sum was zero (ie those genes were not expressed at all). Then I applied the pOverA(0.75,5).

All of the code was run on my personal computer under R version 4.3.1





















The count matrices had different sample names in them. Here's the summary info for matrices that included Sample4, Sample5, and Sample6.

| Matrix   | Library Prep | Genes present | Genes retained after removing rows == 0 | Genes retained after pOverA | Bioinformatics processing | Link to csv                                                                                                                                                                                                                                          |
| -------- | ------------ | ------------- | --------------------------------------- | --------------------------- | ------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Matrix 1 | RiboFree     | 55086         | 53498                                   | 26476                       | Unknown                   | [https://github.com/hputnam/Express_Compare/tree/main/RAnalysis/Data/Mcap](https://github.com/hputnam/Express_Compare/tree/main/RAnalysis/Data/Mcap)                                                                                                 |
| Matrix 2 | RiboFree     | 55086         | 53498                                   | 25708                       | Hisat2                    | [https://github.com/hputnam/Express_Compare/blob/main/hisat2_workflow/counts_matrices/MCap_gene_matrix.csv](https://github.com/hputnam/Express_Compare/blob/main/hisat2_workflow/counts_matrices/MCap_gene_matrix.csv)                               |
| Matrix 5 | RiboFree     | 55086         | 53723                                   | 28835                       | Stringtie                 | [https://github.com/hputnam/Express_Compare/blob/main/stringtie_assembly/counts_matrices/Mcap/redo/gene_count_matrix.csv](https://github.com/hputnam/Express_Compare/blob/main/stringtie_assembly/counts_matrices/Mcap/redo/gene_count_matrix.csv)   |
| Matrix 6 | RiboFree     | 54993         | 53638                                   | 28825                       | Stringtie                 | [https://github.com/hputnam/Express_Compare/blob/main/stringtie_assembly/counts_matrices/Mcap/redo/gene_count_matrix2.csv](https://github.com/hputnam/Express_Compare/blob/main/stringtie_assembly/counts_matrices/Mcap/redo/gene_count_matrix2.csv) |

Here's the summary info for matrices that included X54, X57, and X65.

| -------- | ------------ | ------------- | --------------------------------------- | --------------------------- | ---------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| Matrix 3 | RiboFree     | 63227         | 62673                                   | 46806                       | Fastp, STAR, stringtie       | [https://github.com/hputnam/Express_Compare/blob/main/star_kallisto_workflow/Output/Fastp_STAR_Stringtie_Matrices/Mcap_RiboFree/gene_count_matrix.csv](https://github.com/hputnam/Express_Compare/blob/main/star_kallisto_workflow/Output/Fastp_STAR_Stringtie_Matrices/Mcap_RiboFree/gene_count_matrix.csv)             |
| Matrix 4 | RiboFree     | 63227         | 62707                                   | 48366                       | Trimmomatic, STAR, stringtie | [https://github.com/hputnam/Express_Compare/blob/main/star_kallisto_workflow/Output/Trimmomatic_STAR_Stringtie_Matrices/Mcap_RiboFree/gene_count_matrix.csv](https://github.com/hputnam/Express_Compare/blob/main/star_kallisto_workflow/Output/Trimmomatic_STAR_Stringtie_Matrices/Mcap_RiboFree/gene_count_matrix.csv) |