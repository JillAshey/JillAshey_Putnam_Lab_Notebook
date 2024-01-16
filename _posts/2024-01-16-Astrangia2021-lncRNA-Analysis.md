---
layout: post
title: Astrangia 2021 lncRNA analysis 
date: '2024-01-16'
categories: Analysis
tags: [Bioinformatics, Astrangia, lncRNA]
projects: Astrangia 2021
---

## Astrangia 2021 lncRNA analysis

These data came from my Astrangia 2021 experiment, during which adult Astrangia colonies were exposed to ambient and high temperatures for ~9 months. Sample files were trimmed with fastp, aligned to the reference genome with bowtie2, and assembled with stringtie and gffcompare (code for those steps is [here](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/_posts/2023-06-15-Astrangia2021-mRNA-Analysis.md)). Using the merged GTF file made during the gffcompare step, I will now extract lncRNAs from the data following the [e5 lncRNA discovery workflow](https://github.com/urol-e5/deep-dive/blob/main/D-Apul/code/05.32-lncRNA-discovery-overview.md). 

### 20240116

Make a new lncRNA folder in the Astrangia2021 folder on Andromeda. 

```
cd /data/putnamlab/jillashey/Astrangia2021
mkdir lncRNA
cd lncRNA
mkdir scripts data output
```

Copy the merged GTF file (created during the gffcompare step) in the lncRNA data folder. 

```
cd data

cp /data/putnamlab/jillashey/Astrangia2021/mRNA/stringtie/bowtie/merged.annotated.gtf .
```

Filters the combined GTF output from GFFcompare to select only the lines representing "transcripts" and excluding lines starting with "#" (these are lines in the output format from GFFcompare that don't contain transcript information). This step further filters to keep only those with lengths greater than 199 bases. The size filter of +200nt is a common filtering step for isolating lncRNAs.

```
awk '$3 == "transcript" && $1 !~ /^#/' merged.annotated.gtf | awk '($5 - $4 > 199) || ($4 - $5 > 199)' > Apoc_lncRNA_candidates.gtf

head Apoc_lncRNA_candidates.gtf 
chromosome_1	StringTie	transcript	34636	40489	.	+	.	transcript_id "evm.model.chromosome_1.2"; gene_id "MSTRG.2"; gene_name "evm.TU.chromosome_1.2"; xloc "XLOC_000001"; ref_gene_id "evm.TU.chromosome_1.2"; cmp_ref "evm.model.chromosome_1.2"; class_code "="; tss_id "TSS1";
chromosome_1	StringTie	transcript	43758	53463	.	+	.	transcript_id "evm.model.chromosome_1.3"; gene_id "MSTRG.3"; gene_name "evm.TU.chromosome_1.3"; xloc "XLOC_000002"; ref_gene_id "evm.TU.chromosome_1.3"; cmp_ref "evm.model.chromosome_1.3"; class_code "="; tss_id "TSS2";
chromosome_1	StringTie	transcript	62282	74713	.	+	.	transcript_id "evm.model.chromosome_1.5"; gene_id "MSTRG.5"; gene_name "evm.TU.chromosome_1.5"; xloc "XLOC_000003"; ref_gene_id "evm.TU.chromosome_1.5"; cmp_ref "evm.model.chromosome_1.5"; class_code "="; tss_id "TSS3";
chromosome_1	StringTie	transcript	78661	89242	.	+	.	transcript_id "evm.model.chromosome_1.7.1.5f15db9d"; gene_id "MSTRG.8"; gene_name "evm.TU.chromosome_1.7"; xloc "XLOC_000004"; ref_gene_id "evm.TU.chromosome_1.7"; cmp_ref "evm.model.chromosome_1.7.1.5f15db9d"; class_code "="; tss_id "TSS4";
chromosome_1	StringTie	transcript	78661	92518	.	+	.	transcript_id "evm.model.chromosome_1.7"; gene_id "MSTRG.8"; gene_name "evm.TU.chromosome_1.7"; xloc "XLOC_000004"; ref_gene_id "evm.TU.chromosome_1.7"; cmp_ref "evm.model.chromosome_1.7"; class_code "="; tss_id "TSS4";
chromosome_1	EVM	transcript	158727	159275	.	+	.	transcript_id "evm.model.chromosome_1.12"; gene_id "evm.TU.chromosome_1.12"; gene_name "evm.TU.chromosome_1.12"; xloc "XLOC_000005"; ref_gene_id "evm.TU.chromosome_1.12"; cmp_ref "evm.model.chromosome_1.12"; class_code "="; tss_id "TSS5";
chromosome_1	StringTie	transcript	276514	277413	.	+	.	transcript_id "evm.model.chromosome_1.16"; gene_id "MSTRG.14"; gene_name "evm.TU.chromosome_1.16"; xloc "XLOC_000006"; ref_gene_id "evm.TU.chromosome_1.16"; cmp_ref "evm.model.chromosome_1.16"; class_code "="; tss_id "TSS6";
chromosome_1	StringTie	transcript	295649	299431	.	+	.	transcript_id "evm.model.chromosome_1.18"; gene_id "MSTRG.18"; gene_name "evm.TU.chromosome_1.18"; xloc "XLOC_000007"; ref_gene_id "evm.TU.chromosome_1.18"; cmp_ref "evm.model.chromosome_1.18"; class_code "="; tss_id "TSS7";
chromosome_1	StringTie	transcript	304272	305150	.	+	.	transcript_id "evm.model.chromosome_1.19"; gene_id "MSTRG.17"; gene_name "evm.TU.chromosome_1.19"; xloc "XLOC_000008"; ref_gene_id "evm.TU.chromosome_1.19"; cmp_ref "evm.model.chromosome_1.19"; class_code "="; tss_id "TSS8";
chromosome_1	StringTie	transcript	350318	367541	.	+	.	transcript_id "evm.model.chromosome_1.25"; gene_id "MSTRG.23"; gene_name "evm.TU.chromosome_1.25"; xloc "XLOC_000009"; ref_gene_id "evm.TU.chromosome_1.25"; cmp_ref "evm.model.chromosome_1.25"; class_code "="; tss_id "TSS9";

wc -l Apoc_lncRNA_candidates.gtf 
45484 Apoc_lncRNA_candidates.gtf
```

Use [bedtools](https://bedtools.readthedocs.io/en/latest/content/tools/getfasta.html) to extract sequence data from the reference genome based on the coordinates from the GTF. The resulting seqs represent potential lncRNA candidates. In scripts, `nano fastaFromBed.sh`

```
#!/bin/bash -i
#SBATCH -t 120:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=500GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Astrangia2021/lncRNA/scripts              
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

module load BEDTools/2.30.0-GCC-11.3.0 

echo "Obtaining potential lncRNA seqs based on coordinates" $(date)

fastaFromBed -fi /data/putnamlab/jillashey/Astrangia_Genome/apoculata.assembly.scaffolds_chromosome_level.fasta -bed /data/putnamlab/jillashey/Astrangia2021/lncRNA/data/Apoc_lncRNA_candidates.gtf -fo /data/putnamlab/jillashey/Astrangia2021/lncRNA/data/Apoc_lncRNA_candidates.fasta -name -split

echo "Potential lncRNA fasta created" $(date)
```

Submitted batch job 292598. Took about a min. 

```
zgrep -c ">" Apoc_lncRNA_candidates.fasta 
45484
```

Using python, run CPC2.py, which predicts whether a transcript is coding or non-coding. CPC2 uses ORF (Open Reading Frame) Analysis, Isometric Feature Mapping (Isomap), Sequence Homology, RNA Sequence Features, and Quality of Translation to assess coding potential and flag any transcripts we would want to exclude using the FASTA generated in the previous step. The CPC2 python script is [here](https://github.com/biocoder/CPC2/blob/master/bin/CPC2.py). I copy and pasted it into a .py script in the `scripts` folder on the server. CPC2 also has a [website version](http://cpc2.gao-lab.org/index.php). 

```
module purge
module load Python/2.7.18-GCCcore-9.3.0

python CPC2.py -i /data/putnamlab/jillashey/Astrangia2021/lncRNA/data/Apoc_lncRNA_candidates.fasta -o /data/putnamlab/jillashey/Astrangia2021/lncRNA/output/CPC2
```

Getting this error: 

```
Traceback (most recent call last):
  File "CPC2.py", line 10, in <module>
    import numpy as np
ImportError: No module named numpy
```

Python should have numpy preloaded...maybe need a different version of python? Tried with python 3 didn't work. Gave me a similar error but commands instead of numpy. Let's try this version `Python/2.7.15-foss-2018b`. Okay now getting this error: `Illegal instruction (core dumped)`. Let's try this version `Python/2.7.18-GCCcore-10.2.0`. Now getting the numpy module error again. This [issue](https://github.com/gao-lab/CPC2_standalone/issues/3) on the CPC2 github page was running into a similar error with the commands module. The author said "'commands' is a default module in python 2. If you are using python 3, please try this version https://github.com/gao-lab/CPC2_standalone/releases/tag/v1.0.1". 


These are all of the python 2 versions on Andromeda. Try them all. 

```
module load Python/2.7.15-foss-2018b # illegal instruction error
module load Python/2.7.15-GCCcore-7.3.0-bare # illegal instruction error 
module load Python/2.7.16-GCCcore-8.3.0 # numpy error 
module load Python/2.7.18-GCCcore-10.2.0 # numpy error
module load Python/2.7.18-GCCcore-10.3.0-bare # numpy error 
module load Python/2.7.18-GCCcore-11.2.0-bare # illegal instruction error
module load Python/2.7.18-GCCcore-9.3.0 # numpy error 
```

Okay none of the python 2 versions worked. There is a [CPC2 conda package](https://anaconda.org/cbp44/cpc2/files), but one of the github [issues](https://github.com/gao-lab/CPC2_standalone/issues/5) said that this package is not released or maintained by the Gao Lab so I am hesitant to use it. Going to try to install it anyway. 

```
module load Miniconda3/4.9.2

conda create --prefix /data/putnamlab/cpc2
conda activate /data/putnamlab/cpc2

conda install cbp44::cpc2 
Collecting package metadata (current_repodata.json): done
Solving environment: failed with initial frozen solve. Retrying with flexible solve.
Solving environment: failed with repodata from current_repodata.json, will retry with next repodata source.
Collecting package metadata (repodata.json): done
Solving environment: failed with initial frozen solve. Retrying with flexible solve.
Solving environment: - 
Found conflicts! Looking for incompatible packages.
This can take several minutes.  Press CTRL-C to abort.
failed

ResolvePackageNotFound: 
  - libsvm=3.16
```

Install failed. Going to email Kevin Bryan to see if he can install. 

Once we have the lncRNA fasta, what do we do with it? How do we identify lncRNAs in the samples themselves?

Possible lncRNA-mRNA interaction software

- [RNAplex](https://academic.oup.com/bioinformatics/article/24/22/2657/184477)
- [lncTar](http://www.cuilab.cn/lnctar)
- RNAhybrid