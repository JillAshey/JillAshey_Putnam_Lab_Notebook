---
layout: post
title: e5 closest
date: '2024-11-01'
categories: Analysis
tags: [Bioinformatics]
projects: e5 deep dive
---

## e5 deep dive expression - miRNA target prediction 

I am redoing the miRNA target prediction analysis for the deep dive expression project. Kathleen reran the miRNA short stack with a different version of the software and the results changed slightly (see [issue](https://github.com/urol-e5/deep-dive-expression/issues/3)). 

For miRNA target prediction, the 3'UTR sequences of the genes are needed. These typically aren't annotated in coral genomes but I manually added them by adding 1000bp to the right flank of the gene and subtracting any overlap with nearby genes (see [post](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/_posts/2024-06-15-e5-deepdive-miRNA-TargetPrediction.md)). I already have the 3'UTR sequences for our species of interest woohoo. Download updated miRNA sequences onto Andromeda for target prediction with miranda software. 

Note - miRNA sequences are not named yet as of 11/1/24. 

### Apulchra 

Download new sequences and rename fasta files. Copied fasta sequences from [here](https://github.com/urol-e5/deep-dive-expression/blob/main/D-Apul/output/11-Apul-sRNA-ShortStack_4.1.0-pulchra_genome/ShortStack_out/mir.fasta). 

```
cd /data/putnamlab/jillashey/e5/mirna_seqs
nano Apul_updated_results.fasta # paste fasta sequences 
grep -A 1 'mature' Apul_updated_results.fasta | grep -v '^--$' > Apul_updated_results_mature.fasta

zgrep -c ">" Apul_updated_results_mature.fasta
39
```

In original run (using old short stack version and Amil genome), 38 mature miRNAs were identified. With updated short stack and new Apul genome, 39 mature miRNAs were identified. 

I already have the 3'UTR fasta file here: `/data/putnamlab/jillashey/e5/refs/Apul/Apul_3UTR_1kb.fasta`. In the scripts folder: `nano miranda_strict_all_1kb_apul_updated.sh`

```
#!/bin/bash -i
#SBATCH -t 48:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --export=NONE
#SBATCH --mem=500GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/e5/scripts
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

echo "Apul starting miranda run with all genes and miRNAs with energy cutoff <-20 and strict binding invoked"$(date)
echo "Using updated miRNAs from newer short stack version"

module load Miniconda3/4.9.2
conda activate /data/putnamlab/conda/miranda 

miranda /data/putnamlab/jillashey/e5/mirna_seqs/Apul_updated_results_mature.fasta /data/putnamlab/jillashey/e5/refs/Apul/Apul_3UTR_1kb.fasta -en -20 -strict -out /data/putnamlab/jillashey/e5/output/miranda/miranda_strict_all_1kb_apul_updated.tab

conda deactivate

echo "miranda run finished!"$(date)
echo "counting number of interactions attempted" $(date)

zgrep -c "Performing Scan" /data/putnamlab/jillashey/e5/output/miranda/miranda_strict_all_1kb_apul_updated.tab

echo "Parsing output" $(date)
grep -A 1 "Scores for this hit:" /data/putnamlab/jillashey/e5/output/miranda/miranda_strict_all_1kb_apul.tab | sort | grep '>' > /data/putnamlab/jillashey/e5/output/miranda/miranda_strict_all_1kb_parsed_apul_updated.txt

echo "counting number of putative interactions predicted" $(date)
wc -l /data/putnamlab/jillashey/e5/output/miranda/miranda_strict_all_1kb_parsed_apul_updated.txt

echo "Apul miranda script complete" $(date)
```

Submitted batch job 346351. Results: 

```
counting number of interactions attempted Fri Nov 1 14:02:45 EDT 2024
2062788
Parsing output Fri Nov 1 14:02:46 EDT 2024
counting number of putative interactions predicted Fri Nov 1 14:02:46 EDT 2024
6109 /data/putnamlab/jillashey/e5/output/miranda/miranda_strict_all_1kb_parsed_apul_updated.txt
```

With previous miRNAs and Amil genome, there were 4144 putative interactions predicted. With updated miRNAs and Apul genome, 6109 putative interactions were predicted. 

It is possible that miRNAs also bind to coding regions as well. Perform miranda run using coding sequences as targets. First, extract sequences from genome using bedtools.

```
interactive
cd /data/putnamlab/jillashey/e5/refs/Apul
module load BEDTools/2.30.0-GCC-11.3.0

awk '$3 == "transcript"' genome.fixed.gtf > Apul_transcripts.gtf

bedtools getfasta -fi /data/putnamlab/REFS/Apul/apul.hifiasm.s55_pa.p_ctg.fa.k32.w100.z1000.ntLink.5rounds.masked.fa -bed Apul_transcripts.gtf -fo apul_transcripts.fasta

grep -c ">" apul_transcripts.fasta
44371

wc -l Apul_transcripts.gtf
44371
```

Run miranda. In the scripts folder: `nano miranda_strict_all_mrna_apul_updated.sh`

```
#!/bin/bash -i
#SBATCH -t 48:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --export=NONE
#SBATCH --mem=500GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/e5/scripts
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

echo "Apul starting miranda run with mRNAs as targets and miRNAs with energy cutoff <-20 and strict binding invoked"$(date)
echo "Using updated miRNAs from newer short stack version"

module load Miniconda3/4.9.2
conda activate /data/putnamlab/conda/miranda 

miranda /data/putnamlab/jillashey/e5/mirna_seqs/Apul_updated_results_mature.fasta /data/putnamlab/jillashey/e5/refs/Apul/apul_transcripts.fasta -en -20 -strict -out /data/putnamlab/jillashey/e5/output/miranda/miranda_strict_all_mrna_apul_updated.tab

conda deactivate

echo "miranda run finished!"$(date)
echo "counting number of interactions attempted" $(date)

zgrep -c "Performing Scan" /data/putnamlab/jillashey/e5/output/miranda/miranda_strict_all_mrna_apul_updated.tab

echo "Parsing output" $(date)
grep -A 1 "Scores for this hit:" /data/putnamlab/jillashey/e5/output/miranda/miranda_strict_all_mrna_apul_updated.tab | sort | grep '>' > /data/putnamlab/jillashey/e5/output/miranda/miranda_strict_all_mrna_parsed_apul_updated.txt

echo "counting number of putative interactions predicted" $(date)
wc -l /data/putnamlab/jillashey/e5/output/miranda/miranda_strict_all_mrna_parsed_apul_updated.txt

echo "Apul miranda script complete" $(date)
```

Submitted batch job 361539

### Pevermanni 

Download new sequences and rename fasta files. Copied fasta sequences from [here](https://github.com/urol-e5/deep-dive-expression/blob/main/E-Peve/output/05-Peve-sRNA-ShortStack_4.1.0/ShortStack_out/mir.fasta). 

```
cd /data/putnamlab/jillashey/e5/mirna_seqs
nano Peve_updated_results.fasta # paste fasta sequences 
grep -A 1 'mature' Peve_updated_results.fasta | grep -v '^--$' > Peve_updated_results_mature.fasta

zgrep -c ">" Peve_updated_results_mature.fasta
45
```

In original run (using old short stack version), 46 mature miRNAs were identified. With updated short stack, 45 mature miRNAs were identified. 

I already have the 3'UTR fasta file here: `/data/putnamlab/jillashey/genome/Peve/peve_3UTR_1kb.fasta`. In the scripts folder: `nano miranda_strict_all_1kb_peve_updated.sh`

```
#!/bin/bash -i
#SBATCH -t 48:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --export=NONE
#SBATCH --mem=500GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/e5/scripts
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

echo "Peve starting miranda run with all genes and miRNAs with energy cutoff <-20 and strict binding invoked"$(date)
echo "Using updated miRNAs from newer short stack version"

module load Miniconda3/4.9.2
conda activate /data/putnamlab/conda/miranda 

miranda /data/putnamlab/jillashey/e5/mirna_seqs/Peve_updated_results_mature.fasta /data/putnamlab/jillashey/genome/Peve/peve_3UTR_1kb.fasta -en -20 -strict -out /data/putnamlab/jillashey/e5/output/miranda/miranda_strict_all_1kb_peve_updated.tab

conda deactivate

echo "miranda run finished!"$(date)
echo "counting number of interactions attempted" $(date)

zgrep -c "Performing Scan" /data/putnamlab/jillashey/e5/output/miranda/miranda_strict_all_1kb_peve_updated.tab

echo "Parsing output" $(date)
grep -A 1 "Scores for this hit:" /data/putnamlab/jillashey/e5/output/miranda/miranda_strict_all_1kb_peve_updated.tab | sort | grep '>' > /data/putnamlab/jillashey/e5/output/miranda/miranda_strict_all_1kb_parsed_peve_updated.txt

echo "counting number of putative interactions predicted" $(date)
wc -l /data/putnamlab/jillashey/e5/output/miranda/miranda_strict_all_1kb_parsed_peve_updated.txt

echo "Peve miranda script complete" $(date)
```

Submitted batch job 346352. Results: 

```
counting number of interactions attempted Fri Nov 1 14:05:10 EDT 2024
1750500
Parsing output Fri Nov 1 14:05:17 EDT 2024
counting number of putative interactions predicted Fri Nov 1 14:05:17 EDT 2024
5067 /data/putnamlab/jillashey/e5/output/miranda/miranda_strict_all_1kb_parsed_peve_updated.txt
```

With previous miRNAs, there were 7187 putative interactions predicted. With updated miRNAs, 5067 putative interactions were predicted. 

It is possible that miRNAs also bind to coding regions as well. Perform miranda run using coding sequences as targets. First, extract sequences from genome using bedtools. 

```
interactive
cd /data/putnamlab/jillashey/genome/Peve
module load BEDTools/2.30.0-GCC-11.3.0

bedtools getfasta -fi Porites_evermanni_v1.fa -bed peve_GFFannotation.mRNA_sorted.gff -fo peve_mRNA.fasta
grep -c ">" peve_mRNA.fasta
40389
wc -l peve_GFFannotation.mRNA_sorted.gff
40389
```

Run miranda. In the scripts folder: `nano miranda_strict_all_mrna_peve_updated.sh `

```
#!/bin/bash -i
#SBATCH -t 48:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --export=NONE
#SBATCH --mem=500GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/e5/scripts
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

echo "Peve starting miranda run with mRNAs as targets and miRNAs with energy cutoff <-20 and strict binding invoked"$(date)
echo "Using updated miRNAs from newer short stack version"

module load Miniconda3/4.9.2
conda activate /data/putnamlab/conda/miranda 

miranda /data/putnamlab/jillashey/e5/mirna_seqs/Peve_updated_results_mature.fasta /data/putnamlab/jillashey/genome/Peve/peve_mRNA.fasta -en -20 -strict -out /data/putnamlab/jillashey/e5/output/miranda/miranda_strict_all_mrna_peve_updated.tab

conda deactivate

echo "miranda run finished!"$(date)
echo "counting number of interactions attempted" $(date)

zgrep -c "Performing Scan" /data/putnamlab/jillashey/e5/output/miranda/miranda_strict_all_mrna_peve_updated.tab

echo "Parsing output" $(date)
grep -A 1 "Scores for this hit:" /data/putnamlab/jillashey/e5/output/miranda/miranda_strict_all_mrna_peve_updated.tab | sort | grep '>' > /data/putnamlab/jillashey/e5/output/miranda/miranda_strict_all_mrna_parsed_peve_updated.txt

echo "counting number of putative interactions predicted" $(date)
wc -l /data/putnamlab/jillashey/e5/output/miranda/miranda_strict_all_mrna_parsed_peve_updated.txt

echo "Peve miranda script complete" $(date)
```

Submitted batch job 361535


### Ptuahiniensis

Download new sequences and rename fasta files. Copied fasta sequences from [here](https://github.com/urol-e5/deep-dive-expression/blob/main/F-Ptuh/output/05-Ptuh-sRNA-ShortStack_4.1.0/ShortStack_out/mir.fasta). 

```
cd /data/putnamlab/jillashey/e5/mirna_seqs
nano Ptuh_updated_results.fasta # paste fasta sequences 
grep -A 1 'mature' Ptuh_updated_results.fasta | grep -v '^--$' > Ptuh_updated_results_mature.fasta

zgrep -c ">" Ptuh_updated_results_mature.fasta
37
```

In original run (using old short stack version), 37 mature miRNAs were identified. With updated short stack, 37 mature miRNAs were identified. 

I already have the 3'UTR fasta file here: `/data/putnamlab/jillashey/genome/Pmea/Pmea_3UTR_1kb.fasta`. In the scripts folder: `nano miranda_strict_all_1kb_ptuh_updated.sh`

```
#!/bin/bash -i
#SBATCH -t 48:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --export=NONE
#SBATCH --mem=500GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/e5/scripts
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

echo "Ptuh starting miranda run with all genes and miRNAs with energy cutoff <-20 and strict binding invoked"$(date)
echo "Using updated miRNAs from newer short stack version"

module load Miniconda3/4.9.2
conda activate /data/putnamlab/conda/miranda 

miranda /data/putnamlab/jillashey/e5/mirna_seqs/Ptuh_updated_results_mature.fasta /data/putnamlab/jillashey/genome/Pmea/Pmea_3UTR_1kb.fasta -en -20 -strict -out /data/putnamlab/jillashey/e5/output/miranda/miranda_strict_all_1kb_ptuh_updated.tab

conda deactivate

echo "miranda run finished!"$(date)
echo "counting number of interactions attempted" $(date)

zgrep -c "Performing Scan" /data/putnamlab/jillashey/e5/output/miranda/miranda_strict_all_1kb_ptuh_updated.tab

echo "Parsing output" $(date)
grep -A 1 "Scores for this hit:" /data/putnamlab/jillashey/e5/output/miranda/miranda_strict_all_1kb_ptuh_updated.tab | sort | grep '>' > /data/putnamlab/jillashey/e5/output/miranda/miranda_strict_all_1kb_parsed_ptuh_updated.txt

echo "counting number of putative interactions predicted" $(date)
wc -l /data/putnamlab/jillashey/e5/output/miranda/miranda_strict_all_1kb_parsed_ptuh_updated.txt

echo "Peve miranda script complete" $(date)
```

Submitted batch job 346353. Results: 

```
counting number of interactions attempted Fri Nov 1 14:04:16 EDT 2024
1204979
Parsing output Fri Nov 1 14:04:19 EDT 2024
counting number of putative interactions predicted Fri Nov 1 14:04:19 EDT 2024
4105 /data/putnamlab/jillashey/e5/output/miranda/miranda_strict_all_1kb_parsed_ptuh_updated.txt
```

With previous miRNAs, there were 3863 putative interactions predicted. With updated miRNAs, 4105 putative interactions were predicted. 

It is possible that miRNAs also bind to coding regions as well. Perform miranda run using coding sequences as targets. First, extract sequences from genome using bedtools. 

```
interactive
cd /data/putnamlab/jillashey/genome/Pmea
module load BEDTools/2.30.0-GCC-11.3.0

bedtools getfasta -fi Pocillopora_meandrina_HIv1.assembly.fasta -bed Pmea_transcript.gtf -fo pmea_mrna.fasta
grep -c ">" pmea_mrna.fasta
31840
wc -l Pmea_transcript.gtf
31840
```

Run miranda. In the scripts folder: `nano miranda_strict_all_mrna_ptuh_updated.sh`

```
#!/bin/bash -i
#SBATCH -t 48:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --export=NONE
#SBATCH --mem=500GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/e5/scripts
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

echo "Pmea starting miranda run with mRNAs as targets and miRNAs with energy cutoff <-20 and strict binding invoked"$(date)
echo "Using updated miRNAs from newer short stack version"

module load Miniconda3/4.9.2
conda activate /data/putnamlab/conda/miranda 

miranda /data/putnamlab/jillashey/e5/mirna_seqs/Ptuh_updated_results_mature.fasta /data/putnamlab/jillashey/genome/Pmea/pmea_mrna.fasta -en -20 -strict -out /data/putnamlab/jillashey/e5/output/miranda/miranda_strict_all_mrna_ptuh_updated.tab

conda deactivate

echo "miranda run finished!"$(date)
echo "counting number of interactions attempted" $(date)

zgrep -c "Performing Scan" /data/putnamlab/jillashey/e5/output/miranda/miranda_strict_all_mrna_ptuh_updated.tab

echo "Parsing output" $(date)
grep -A 1 "Scores for this hit:" /data/putnamlab/jillashey/e5/output/miranda/miranda_strict_all_mrna_ptuh_updated.tab | sort | grep '>' > /data/putnamlab/jillashey/e5/output/miranda/miranda_strict_all_mrna_parsed_ptuh_updated.txt

echo "counting number of putative interactions predicted" $(date)
wc -l /data/putnamlab/jillashey/e5/output/miranda/miranda_strict_all_mrna_parsed_ptuh_updated.txt

echo "Pmea miranda script complete" $(date)
```

Submitted batch job 361538. Ran in a little under an hour. Results: 

```
counting number of interactions attempted Thu Feb 13 11:34:11 EST 2025
1178080
Parsing output Thu Feb 13 11:34:14 EST 2025
counting number of putative interactions predicted Thu Feb 13 11:34:15 EST 2025
26311 /data/putnamlab/jillashey/e5/output/miranda/miranda_strict_all_mrna_parsed_ptuh_updated.txt
```

WAY may putative interactions predicted using mRNAs than 3'UTR. 