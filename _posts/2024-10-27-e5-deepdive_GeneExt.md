---
layout: post
title: e5 closest
date: '2024-08-19'
categories: Analysis
tags: [Bioinformatics]
projects: e5 deep dive
---

## e5 deep dive - extending genes to identify 3' UTRs 

For the e5 deep dive project, we want to identify the genes that the miRNAs bind to. Mature miRNAs usually bind to the 3' UTR portion of an mRNA but our coral gtfs/gffs do not have 3' UTRs annotated. To annotate them, I'm going to use [GeneExt](https://github.com/sebepedroslab/GeneExt), which is a program that extends the genes to obtain the 3' UTR information. [Zoe](https://github.com/zdellaert/LaserCoral/blob/main/code/RNA-seq-bioinf.md) the gene extending queen, has gotten gene ext on unity (I was having a lot of trouble installing it on Andromeda), so I am going to run gene ext on [Unity](https://ood.unity.rc.umass.edu/pun/sys/dashboard/batch_connect/sys/bc_desktop/session_contexts/new). 


### Apulchra 

I am going to test gene ext on Apulchra first. I downloaded the [Apul bam files](https://gannet.fish.washington.edu/seashell/bu-github/deep-dive-expression/D-Apul/output/07.2-Apul-Hisat/) from gannet -- SR alined the deep dive reads to the Apul genome. 

Make folders on unity 

```
cd /project/pi_hputnam_uri_edu/
mkdir jillashey
cd jillashey 
mkdir e5_deepdive
cd e5_deepdive
mkdir data scripts output 
```

Download the bam files from gannet to unity. 

```
cd data 
wget https://gannet.fish.washington.edu/seashell/bu-github/deep-dive-expression/D-Apul/output/07.2-Apul-Hisat/RNA-ACR-140.sorted.bam
https://gannet.fish.washington.edu/seashell/bu-github/deep-dive-expression/D-Apul/output/07.2-Apul-Hisat/RNA-ACR-145.sorted.bam
https://gannet.fish.washington.edu/seashell/bu-github/deep-dive-expression/D-Apul/output/07.2-Apul-Hisat/RNA-ACR-150.sorted.bam
https://gannet.fish.washington.edu/seashell/bu-github/deep-dive-expression/D-Apul/output/07.2-Apul-Hisat/RNA-ACR-173.sorted.bam
https://gannet.fish.washington.edu/seashell/bu-github/deep-dive-expression/D-Apul/output/07.2-Apul-Hisat/RNA-ACR-178.sorted.bam
```

In the scripts folder: `nano merge_bam.sh`

```
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --mem=100GB
#SBATCH --export=NONE
#SBATCH --error="%x_error.%j" #write out slurm error reports
#SBATCH --output="%x_output.%j" #write out any program outpus
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH -D /project/pi_hputnam_uri_edu/jillashey/e5_deepdive/scripts #set working directory
#SBATCH --constraint=avx512

#load modules
module load uri/main
module load SAMtools/1.16.1-GCC-11.3.0

cd /project/pi_hputnam_uri_edu/jillashey/e5_deepdive/data

#use samtools merge to merge all the files
samtools merge merge.bam *.bam
```

Submitted batch job 25761547. The resulting file, `merge.bam`, should be as large or larger than the sum of the individual bam files. Once this is done running, run gene ext! 

In the scripts folder: `nano GeneExt.sh`

```
#!/bin/bash
#SBATCH -t 18:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=500GB
#SBATCH --export=NONE
#SBATCH --error=../scripts/outs_errs/“%x_error.%j” #write out slurm error reports
#SBATCH --output=../scripts/outs_errs/“%x_output.%j” #write out any program outpus
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH -D /project/pi_hputnam_uri_edu/jillashey/e5_deepdive/scripts #set working directory

# activate environment
module load conda/latest #load miniconda
conda activate /work/pi_hputnam_uri_edu/conda/envs/GeneExt/geneext

echo "Environment activated, run Gene ext" $(date)

# use --clip_strand both to not allow GeneExt to create overlaps on the same strand
python /work/pi_hputnam_uri_edu/conda/envs/GeneExt/geneext.py --verbose=3 \
    -g /project/pi_hputnam_uri_edu/jillashey/e5_deepdive/data/Apulchra-genome.gtf \
    -b /project/pi_hputnam_uri_edu/jillashey/e5_deepdive/data/merge.bam \
    -o /project/pi_hputnam_uri_edu/jillashey/e5_deepdive/output/Apul_GeneExt.gtf \
    -j $SLURM_CPUS_PER_TASK --clip_strand both --force
    
echo "Gene Ext complete, deactivate env" $(date)
conda deactivate
```

Submitted batch job 25762636. I got output but not what I thought...There is not gtf file in the output folder, only a log file (`Apul_GeneExt.gtf.GeneExt.log`. The log file has this info at the end: 

```
FUN_044369-T1
FUN_044370-T1
1 transcripts found.
FUN_044370-T1
FUN_044371-T1
1 transcripts found.
FUN_044371-T1
Fix done, annotation with gene features: tmp/genome.fixed.gtf
Added gene feature names. New file name: tmp/genome.fixed.gtf
```

Maybe it needed the gene features added before it could extend the genes? In the tmp folder:

```
wc -l genome.fixed.gtf 
298279 genome.fixed.gtf

head genome.fixed.gtf 
ntLink_0        funannotate     gene    1105    7056    .       +       .       gene_id "FUN_000001"
ntLink_0        funannotate     transcript      1105    7056    .       +       .       transcript_id "FUN_000001-T1"; gene_id "FUN_000001";
ntLink_0        funannotate     exon    1105    1188    .       +       .       transcript_id "FUN_000001-T1"; gene_id "FUN_000001";
ntLink_0        funannotate     exon    1861    1941    .       +       .       transcript_id "FUN_000001-T1"; gene_id "FUN_000001";
ntLink_0        funannotate     exon    2762    2839    .       +       .       transcript_id "FUN_000001-T1"; gene_id "FUN_000001";
ntLink_0        funannotate     exon    5044    7056    .       +       .       transcript_id "FUN_000001-T1"; gene_id "FUN_000001";
ntLink_0        funannotate     gene    10215   15286   .       +       .       gene_id "FUN_000002"
ntLink_0        funannotate     transcript      10215   15286   .       +       .       transcript_id "FUN_000002-T1"; gene_id "FUN_000002";
ntLink_0        funannotate     exon    13074   14383   .       +       .       transcript_id "FUN_000002-T1"; gene_id "FUN_000002";
ntLink_0        funannotate     exon    14722   14900   .       +       .       transcript_id "FUN_000002-T1"; gene_id "FUN_000002";
```

Compared to the actual gtf: 

```
wc -l Apulchra-genome.gtf 
455521 Apulchra-genome.gtf

head Apulchra-genome.gtf 
ntLink_0        funannotate     transcript      1105    7056    .       +       .       transcript_id "FUN_000001-T1"; gene_id "FUN_000001"
ntLink_0        funannotate     exon    1105    1188    .       +       .       transcript_id "FUN_000001-T1"; gene_id "FUN_000001";
ntLink_0        funannotate     exon    1861    1941    .       +       .       transcript_id "FUN_000001-T1"; gene_id "FUN_000001";
ntLink_0        funannotate     exon    2762    2839    .       +       .       transcript_id "FUN_000001-T1"; gene_id "FUN_000001";
ntLink_0        funannotate     exon    5044    7056    .       +       .       transcript_id "FUN_000001-T1"; gene_id "FUN_000001";
ntLink_0        funannotate     CDS     1105    1188    .       +       0       transcript_id "FUN_000001-T1"; gene_id "FUN_000001";
ntLink_0        funannotate     CDS     1861    1941    .       +       0       transcript_id "FUN_000001-T1"; gene_id "FUN_000001";
ntLink_0        funannotate     CDS     2762    2839    .       +       0       transcript_id "FUN_000001-T1"; gene_id "FUN_000001";
ntLink_0        funannotate     CDS     5044    7056    .       +       0       transcript_id "FUN_000001-T1"; gene_id "FUN_000001";
ntLink_0        funannotate     transcript      10215   15286   .       +       .       transcript_id "FUN_000002-T1"; gene_id "FUN_000002"
```

It looks like gene ext added gene rows and removed CDSs from the file? Am I supposed to use this gtf to run gene ext? Maybe. Going to move it from the tmp directory into the data directory in case gene ext deletes the tmp directory. 

```
mv genome.fixed.gtf ../../data/
```

Edit the `Gene_Ext.sh` script so that the `-g` option is directed to the fixed gtf. Submitted batch job 25767632


