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

### 20240321

It's been a while. Kevin Bryan did install CPC2 on the server (`CPC2/1.0.1-foss-2022a`) so I can run it there. I already copied the [CPC2](https://github.com/gao-lab/CPC2_standalone/blob/master/bin/CPC2.py) python script to the server. In the scripts folder: `nano cpc2.sh`

```
#!/bin/bash 
#SBATCH -t 100:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --export=NONE
#SBATCH --mem=250GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Astrangia2021/lncRNA/scripts              
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

module avail Python/3.7.0-foss-2018b
module avail CPC2/1.0.1-foss-2022a

echo "Evaluating coding potential of lncRNA candidates" $(date)

python CPC2.py -i /data/putnamlab/jillashey/Astrangia2021/lncRNA/data/Apoc_lncRNA_candidates.fasta -o /data/putnamlab/jillashey/Astrangia2021/lncRNA/output/CPC2

echo "Evaluation complete!" $(date)
```

Submitted batch job 309769. Ran for about 20 seconds and got this error: 

```
------------------------------- /opt/modules/all -------------------------------
Python/3.7.0-foss-2018b

------------------------------- /opt/modules/all -------------------------------
CPC2/1.0.1-foss-2022a
Traceback (most recent call last):
  File "CPC2.py", line 6, in <module>
    import commands
ModuleNotFoundError: No module named 'commands'
```

This looks like a python issue...maybe try a different version. On the [CPC2 github](https://github.com/gao-lab/CPC2_standalone), it says that "This is a python 2 verison of CPC2" so I'll try a python 2 version: `Python/2.7.15-foss-2018b`. I am also now realizing the I put 'module avail' instead of 'module load' in the code so fixing that. Submitted batch job 309770. Failed, got a lot of module conflicting with other modules. 

```
foss/2022a(24):ERROR:150: Module 'foss/2022a' conflicts with the currently loaded module(s) 'foss/2018b'
foss/2022a(24):ERROR:102: Tcl command execution failed: conflict foss

Python/3.10.4-GCCcore-11.3.0(62):ERROR:150: Module 'Python/3.10.4-GCCcore-11.3.0' conflicts with the currently loaded module(s) 'Python/2.7.15-foss
-2018b'
Python/3.10.4-GCCcore-11.3.0(62):ERROR:102: Tcl command execution failed: conflict Python

GCCcore/11.3.0(24):ERROR:150: Module 'GCCcore/11.3.0' conflicts with the currently loaded module(s) 'GCCcore/7.3.0'
GCCcore/11.3.0(24):ERROR:102: Tcl command execution failed: conflict GCCcore

binutils/.2.38-GCCcore-11.3.0(22):ERROR:150: Module 'binutils/.2.38-GCCcore-11.3.0' conflicts with the currently loaded module(s) 'binutils/2.30-GCCcore-7.3.0'
binutils/.2.38-GCCcore-11.3.0(22):ERROR:102: Tcl command execution failed: conflict binutils

foss/2022a(24):ERROR:150: Module 'foss/2022a' conflicts with the currently loaded module(s) 'foss/2018b'
foss/2022a(24):ERROR:102: Tcl command execution failed: conflict foss
```

Adding `module purge` above the module load. Nope didn't like that. I may have to load a specific version of GCCcore. I'm going to do: 

```
module load GCCcore/7.3.0 
module load Python/2.7.15-foss-2018b 
module load CPC2/1.0.1-foss-2022a
```

Submitted batch job 309773. Now getting this error: 

```
------------------------------- /opt/modules/all -------------------------------
CPC2/1.0.1-foss-2022a
Traceback (most recent call last):
  File "CPC2.py", line 11, in <module>
    from Bio.Seq import Seq
ImportError: No module named Bio.Seq
```

Try with this version: `Python/2.7.15-GCCcore-7.3.0-bare`. Submitted batch job 309774. Now getting this error: 

```
------------------------------- /opt/modules/all -------------------------------
CPC2/1.0.1-foss-2022a
Traceback (most recent call last):
  File "CPC2.py", line 10, in <module>
    import numpy as np
ImportError: No module named numpy
```

Let's try this: 

```
module load GCCcore/9.3.0 
module load Python/3.8.2-GCCcore-9.3.0
```

Submitted batch job 309775. Back to the commands error...ChatGPT recommended that I try this: 

```
module purge  # Unload conflicting modules
module load GCCcore/7.3.0  # Load compatible GCCcore version
module load Python/2.7.15-foss-2018b  # Load Python 2 environment
pip install numpy  # Install numpy package
pip install biopython  # Install biopython package (which includes Bio.Seq)
```

Added it to the script. Submitted batch job 309776. Failed with this error: 

```
You are using pip version 10.0.1, however version 20.3.4 is available.
You should consider upgrading via the 'pip install --upgrade pip' command.
Command "python setup.py egg_info" failed with error code 1 in /tmp/pip-install-PyLt8e/biopython/
You are using pip version 10.0.1, however version 20.3.4 is available.
You should consider upgrading via the 'pip install --upgrade pip' command.
Traceback (most recent call last):
  File "CPC2.py", line 11, in <module>
    from Bio.Seq import Seq
ImportError: No module named Bio.Seq
```

ChatGPT recommended this: `pip install --upgrade pip`. Adding to code. Submitted batch job 309777. Basically the same error as above along with: 

```
Could not install packages due to an EnvironmentError: [Errno 13] Permission denied: '/opt/software/Python/2.7.15-foss-2018b/bin/pip'
Consider using the `--user` option or check the permissions.
```

Bleh. I could try on the online version. If submitting in fasta format, it must be 50 Mb max. Going to try to break it up: 

```
cd ../data
awk '/^>/{s=++d".fasta"} {print > s}' RS= d=0 Apoc_lncRNA_candidates.fasta
```

Didn't work. I'm seeing that when I load the CPC2 module, it can autofill `CPC2.py`, so maybe I don't even need to load python. Submitted batch job 309778 (also removed `python` from beginning of code). Got this error: 

```
Traceback (most recent call last):
  File "/opt/software/CPC2/1.0.1-foss-2022a/bin/CPC2.py", line 374, in <module>
    sys.exit(__main())
  File "/opt/software/CPC2/1.0.1-foss-2022a/bin/CPC2.py", line 47, in __main
    if calculate_potential(options.fasta,strand,output_orf,options.outfile):
  File "/opt/software/CPC2/1.0.1-foss-2022a/bin/CPC2.py", line 262, in calculate_potential
    ftmp_result = open(outfile,"w")
IsADirectoryError: [Errno 21] Is a directory: '/data/putnamlab/jillashey/Astrangia2021/lncRNA/output/CPC2'
```

Going to remove the `-o` flag. Submitted batch job 309780. This appears to have worked!!!!! Ran very fast and output a file that includes each transcript and its coding/noncoding potential. Move this file to the output folder. 

```
mv cpc2output.txt.txt ../output/CPC2
cd ../output/CPC2
```

Remove transcripts with the label 'coding'. We only keeping noncoding!

```
awk '$8 != "coding"' cpc2output.txt.txt > apoc_noncoding_transcripts_info.txt

awk '$8 == "noncoding" {print $1}' cpc2output.txt.txt > apoc_noncoding_transcripts_ids.txt

grep -Fwf apoc_noncoding_transcripts_ids.txt /data/putnamlab/jillashey/Astrangia2021/lncRNA/data/Apoc_lncRNA_candidates.fasta > apoc_merged_final_lncRNAs.gtf

wc -l *
   28981 apoc_merged_final_lncRNAs.gtf
   28981 apoc_noncoding_transcripts_ids.txt
   28982 apoc_noncoding_transcripts_info.txt
   45485 cpc2output.txt.txt
```

Hooray!! Remove duplicate lines from gtf

```
awk '!seen[$0]++' apoc_merged_final_lncRNAs.gtf > apoc_deduplicated_final_lncRNAs.gtf

wc -l apoc_deduplicated_final_lncRNAs.gtf 
28810 apoc_deduplicated_final_lncRNAs.gtf
``` 

Started with 45485 potential lncRNA candidates and finished with 28810 putative lncRNAs. Format gtf to bed file. 

```
awk -F":|-" '{print $3 "\t" $4 "\t" $5}' apoc_deduplicated_final_lncRNAs.gtf > apoc_deduplicated_final_lncRNAs.bed
```

The bed file should have the chromosome name, start position and stop position. Use bedtools `getfasta` to extract lncRNA sequences from genome. 

```
module load BEDTools/2.30.0-GCC-11.3.0 

bedtools getfasta -fi /data/putnamlab/jillashey/Astrangia_Genome/apoculata.assembly.scaffolds_chromosome_level.fasta -bed apoc_deduplicated_final_lncRNAs.bed -fo apoc_bedtools_lncRNAs.fasta -name

zgrep -c ">" apoc_bedtools_lncRNAs.fasta 
28810
```

Yay! We now have a fasta of our lncRNAs. Now we want to quantify them using [kallisto](https://pachterlab.github.io/kallisto/manual) (following Zach's [workflow](https://github.com/zbengt/oyster-lnc/blob/main/code/10-count-matrices-DESeq2-final.Rmd) for quantifying lncRNAs). First, make a kallisto folder in the output folder. 

```
cd /data/putnamlab/jillashey/Astrangia2021/lncRNA/output
mkdir kallisto
cd kallisto
```

The lncRNA fasta that was just created (`apoc_bedtools_lncRNAs.fasta`) will be used to generate the kallisto index. Then, the RNA-seq reads (located `/data/putnamlab/jillashey/Astrangia2021/mRNA/data/trim`) will be aligned to the lncRNA index using kallisto. In the scripts folder: `nano kallisto.sh`

```
#!/bin/bash 
#SBATCH -t 100:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --export=NONE
#SBATCH --mem=250GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Astrangia2021/lncRNA/scripts              
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

module load kallisto/0.48.0-gompi-2022a 

echo "Creating kallisto index from lncRNA fasta" $(date)

kallisto index -i /data/putnamlab/jillashey/Astrangia2021/lncRNA/output/kallisto/apoc_lncRNA_index /data/putnamlab/jillashey/Astrangia2021/lncRNA/output/CPC2/apoc_bedtools_lncRNAs.fasta 

echo "lncRNA index creation complete, starting read alignment" $(date)

array=($(ls /data/putnamlab/jillashey/Astrangia2021/mRNA/data/trim/*R1_001.fastq.gz)) # call the clean sequences - make an array to align

for i in ${array[@]}; do
kallisto quant -i /data/putnamlab/jillashey/Astrangia2021/lncRNA/output/kallisto/apoc_lncRNA_index -o /data/putnamlab/jillashey/Astrangia2021/lncRNA/output/kallisto/kallisto.${i} ${i} $(echo ${i}|sed s/_R1/_R2/) 
done 

echo "lncRNA alignment complete!" $(date)
```

Submitted batch job 309785. Index was create successfully, but quant failed. I got this error for all of my samples: 

```
[quant] fragment length distribution will be estimated from the data
Error: could not create directory /data/putnamlab/jillashey/Astrangia2021/lncRNA/output/kallisto/kallisto./data/putnamlab/jillashey/Astrangia2021/mRNA/data/trim/test.AST-1065_R1_001.fastq.gz
```

I modified the for loop in the code (and commented out the indexing step, as that was already done)

```
for i in ${array[@]}; do
    # Extract just the filename from the input FASTQ file path
    filename=$(basename ${i})
    # Construct the output directory path
    output_dir="/data/putnamlab/jillashey/Astrangia2021/lncRNA/output/kallisto/kallisto.${filename}"
    # Run kallisto quant
    kallisto quant -i /data/putnamlab/jillashey/Astrangia2021/lncRNA/output/kallisto/apoc_lncRNA_index -o ${output_dir} ${i} $(echo ${i} | sed 's/_R1/_R2/')
done
```

Submitted batch job 309786. finished running in about 5 hours. Now I need to convert the gene abundances that kallisto generated into a count matrix for DESeq2. It appears that Zach did this two different ways: one way using [Trinity](https://github.com/zbengt/oyster-lnc/blob/main/code/10-count-matrices-DESeq2-final.Rmd) and one way using an [R script](https://zbengt.github.io/2023-03-09-Mergin-Kallisto_Abundance/). I'll come back to this. 

```
Trinity/2.12.0-foss-2019b-Python-3.7.4
abundance_estimates_to_matrix
```

### 20240322 

Now that I have identified the lncRNAs, I can run blast with lncRNAs as the query and protein, genome, mRNAs and smRNAs in input dbs (similar to what has been done in the e5 [code](https://github.com/urol-e5/deep-dive/blob/main/E-Peve/code/02-Peve-lncRNA-align.md)). 

```
cd /data/putnamlab/jillashey/Astrangia2021/lncRNA/
mkdir blast 
```

All of my lncRNA blast output will go into the blast folder above. For now, I will not set any specific e-value or bitscore cutoffs. I will write 4 different scripts for blasting protein, genome, mRNAs and smRNAs. In the scripts folder: `nano blastx_lncRNA.sh`

```
#!/bin/bash 
#SBATCH -t 100:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --export=NONE
#SBATCH --mem=250GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Astrangia2021/lncRNA/scripts              
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

module load BLAST+/2.13.0-gompi-2022a

echo "Use the prot db already created in small RNA folder, do not have to make new db. Blasting lncRNAs against prot db" $(date)

blastx -query /data/putnamlab/jillashey/Astrangia2021/lncRNA/output/CPC2/apoc_bedtools_lncRNAs.fasta -db /data/putnamlab/jillashey/Astrangia2021/smRNA/data/apoc_prot_db -out /data/putnamlab/jillashey/Astrangia2021/lncRNA/blast/blastx_prot_lncRNA_query.tab -outfmt 6

echo "Number of hits?"
wc -l /data/putnamlab/jillashey/Astrangia2021/lncRNA/blast/blastx_prot_lncRNA_query.tab

echo "File header"
head /data/putnamlab/jillashey/Astrangia2021/lncRNA/blast/blastx_prot_lncRNA_query.tab

echo "Blast complete!" $(date)
```

Submitted batch job 309831. Next, blast lncRNAs to genome. In the scripts folder: `nano blastn_genome_lncRNA.sh`

```
#!/bin/bash 
#SBATCH -t 100:00:00
#SBATCH --nodes=1 --ntasks-per-node=5
#SBATCH --export=NONE
#SBATCH --mem=250GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Astrangia2021/lncRNA/scripts              
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

module load BLAST+/2.13.0-gompi-2022a

echo "Making blast db from genome" $(date)

makeblastdb -in /data/putnamlab/jillashey/Astrangia_Genome/apoculata.assembly.scaffolds_chromosome_level.fasta -dbtype nucl -out /data/putnamlab/jillashey/Astrangia_Genome/apoc_genome_db

echo "Blast db creation complete, blasting lncRNAs against db" $(date)

blastn -query /data/putnamlab/jillashey/Astrangia2021/lncRNA/output/CPC2/apoc_bedtools_lncRNAs.fasta -db /data/putnamlab/jillashey/Astrangia_Genome/apoc_genome_db -out /data/putnamlab/jillashey/Astrangia2021/lncRNA/blast/blastn_genome_lncRNA_query.tab -outfmt 6

echo "Number of hits?"
wc -l /data/putnamlab/jillashey/Astrangia2021/lncRNA/blast/blastn_genome_lncRNA_query.tab

echo "File header"
head /data/putnamlab/jillashey/Astrangia2021/lncRNA/blast/blastn_genome_lncRNA_query.tab

echo "Blast complete!" $(date)
```

Submitted batch job 309837. Next, blast lncRNAs to mRNAs. This is probably the one that I am most interested in. In the scripts folder: `nano blastn_mRNA_lncRNA.sh`

```
#!/bin/bash 
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --export=NONE
#SBATCH --mem=250GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Astrangia2021/lncRNA/scripts              
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

module load BLAST+/2.13.0-gompi-2022a

echo "Use the mRNA db already created in small RNA folder, do not have to make new db. Blasting lncRNAs against mRNA db" $(date)

blastn -query /data/putnamlab/jillashey/Astrangia2021/lncRNA/output/CPC2/apoc_bedtools_lncRNAs.fasta -db /data/putnamlab/jillashey/Astrangia2021/smRNA/data/apoc_mRNA_db -out /data/putnamlab/jillashey/Astrangia2021/lncRNA/blast/blastn_mRNA_lncRNA_query.tab -outfmt 6

echo "Number of hits?"
wc -l /data/putnamlab/jillashey/Astrangia2021/lncRNA/blast/blastn_mRNA_lncRNA_query.tab

echo "File header"
head /data/putnamlab/jillashey/Astrangia2021/lncRNA/blast/blastn_mRNA_lncRNA_query.tab

echo "Blast complete!" $(date)
```

Submitted batch job 309836. Lastly, blast lncRNAs to smRNAs. In the scripts folder: `nano blastn_smRNA_lncRNA.sh`

```
#!/bin/bash 
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --export=NONE
#SBATCH --mem=250GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Astrangia2021/lncRNA/scripts              
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

module load BLAST+/2.13.0-gompi-2022a

echo "Making blast db from smRNA" $(date)

makeblastdb -in /data/putnamlab/jillashey/Astrangia2021/smRNA/mirdeep2/all/mirna_results_04_02_2024_t_11_15_57/mature_all.fa -dbtype nucl -out /data/putnamlab/jillashey/Astrangia2021/smRNA/mirdeep2/all/mirna_results_04_02_2024_t_11_15_57/apoc_smRNA_db

echo "Blast db creation complete, blasting lncRNAs against db" $(date)

blastn -query /data/putnamlab/jillashey/Astrangia2021/lncRNA/output/CPC2/apoc_bedtools_lncRNAs.fasta -db /data/putnamlab/jillashey/Astrangia2021/smRNA/mirdeep2/all/mirna_results_04_02_2024_t_11_15_57/apoc_smRNA_db -out /data/putnamlab/jillashey/Astrangia2021/lncRNA/blast/blastn_smRNA_lncRNA_query.tab -outfmt 6

echo "Number of hits?"
wc -l /data/putnamlab/jillashey/Astrangia2021/lncRNA/blast/blastn_smRNA_lncRNA_query.tab

echo "File header"
head /data/putnamlab/jillashey/Astrangia2021/lncRNA/blast/blastn_smRNA_lncRNA_query.tab

echo "Blast complete!" $(date)
```

Submitted batch job 309839

### 20240325

The blast scripts all pended for a while and then ran. 

```
cd /data/putnamlab/jillashey/Astrangia2021/lncRNA/blast

ls -othr
total 6.3G
-rw-r--r--. 1 jillashey 102M Mar 22 19:39 blastn_mRNA_lncRNA_query.tab
-rw-r--r--. 1 jillashey 5.3G Mar 22 21:07 blastn_genome_lncRNA_query.tab
-rw-r--r--. 1 jillashey    0 Mar 22 21:18 blastn_smRNA_lncRNA_query.tab
-rw-r--r--. 1 jillashey 870M Mar 22 23:31 blastx_prot_lncRNA_query.tab
```

There were no hits for the lncRNAs against the small RNAs which is strange. I would have assumed there would be some hits. Let's look at the output file for each blast result. 

lncRNA blasting to protein sequences: 

```
Number of hits?
7966677 /data/putnamlab/jillashey/Astrangia2021/lncRNA/blast/blastx_prot_lncRNA_query.tab
File header
::chromosome_1:34635-40489      protein|evm.model.chromosome_1.2        98.113  53      1       0       252     410     32      84      3.13e-27
        111
::chromosome_1:34635-40489      protein|evm.model.chromosome_1.2        95.238  42      2       0       5726    5851    158     199     5.91e-18
        84.7
::chromosome_1:34635-40489      protein|evm.model.chromosome_1.2        98.077  52      1       0       2918    3073    85      136     2.18e-16
        80.1
::chromosome_1:34635-40489      protein|evm.model.chromosome_1.2        72.000  50      12      1       4       147     1       50      2.55e-13
        71.2
::chromosome_1:34635-40489      protein|evm.model.chromosome_1.2        100.000 23      0       0       5175    5243    137     159     2.6     32.7
::chromosome_1:34635-40489      protein|evm.model.chromosome_9.2265     44.615  65      28      2       899     1093    202     258     1.78e-06
        52.0
::chromosome_1:34635-40489      protein|evm.model.chromosome_3.1096     44.615  65      28      2       899     1093    336     392     2.93e-06
        52.4
::chromosome_1:34635-40489      protein|evm.model.chromosome_1.366      60.526  38      15      0       980     1093    1007    1044    4.44e-05
        49.3
::chromosome_1:34635-40489      protein|evm.model.chromosome_11.2373    58.333  36      15      0       983     1090    346     381     7.66e-05
        47.8
::chromosome_1:34635-40489      protein|evm.model.chromosome_7.1885     40.000  50      29      1       947     1093    17      66      8.63e-05
```

lncRNA blasting to genome: 

```
Number of hits?
55347126 /data/putnamlab/jillashey/Astrangia2021/lncRNA/blast/blastn_genome_lncRNA_query.tab
File header
::chromosome_1:34635-40489      chromosome_1    100.000 5854    0       0       1       5854    34636   40489   0.0     10811
::chromosome_1:34635-40489      chromosome_1    81.850  1427    101     75      3229    4529    1503669 1502275 0.0     1055
::chromosome_1:34635-40489      chromosome_1    89.669  242     15      5       495     730     4524407 4524644 1.07e-78        300
::chromosome_1:34635-40489      chromosome_1    87.552  241     19      8       497     729     13968647        13968410        3.03e-69        268
::chromosome_1:34635-40489      chromosome_1    87.554  233     21      4       495     724     17099195        17098968        1.41e-67        263
::chromosome_1:34635-40489      chromosome_1    84.644  267     36      5       1790    2055    17004614        17004352        5.07e-67        261
::chromosome_1:34635-40489      chromosome_1    87.773  229     17      8       505     727     3663790 3664013 6.56e-66        257
::chromosome_1:34635-40489      chromosome_1    87.069  232     22      4       495     723     3648353 3648127 2.36e-65        255
::chromosome_1:34635-40489      chromosome_1    87.215  219     17      3       4195    4412    13516147        13516355        2.38e-60        239
::chromosome_1:34635-40489      chromosome_1    76.860  484     57      19      3707    4161    13515581        13516038        2.39e-55        222
```

lncRNA blasting to mRNA sequences: 

```
Number of hits?
1071243 /data/putnamlab/jillashey/Astrangia2021/lncRNA/blast/blastn_mRNA_lncRNA_query.tab
File header
::chromosome_1:34635-40489      evm.model.chromosome_1.2        100.000 159     0       0       253     411     98      256     6.28e-78        294
::chromosome_1:34635-40489      evm.model.chromosome_1.2        100.000 157     0       0       2919    3075    257     413     8.12e-77        291
::chromosome_1:34635-40489      evm.model.chromosome_1.2        100.000 124     0       0       5731    5854    480     603     1.80e-58        230
::chromosome_1:34635-40489      evm.model.chromosome_1.2        100.000 99      0       0       1       99      1       99      1.42e-44        183
::chromosome_1:34635-40489      evm.model.chromosome_1.2        100.000 68      0       0       5175    5242    412     479     2.42e-27        126
::chromosome_1:62281-74713      evm.model.chromosome_3.1716     99.485  1165    6       0       5287    6451    1329    165     0.0     2119
::chromosome_1:62281-74713      evm.model.chromosome_9.2364     91.733  1137    84      2       5287    6414    1136    1       0.0     1570
::chromosome_1:62281-74713      evm.model.chromosome_14.123     96.715  822     13      9       5630    6451    1857    1050    0.0     1356
::chromosome_1:62281-74713      evm.model.chromosome_5.1437     98.503  735     3       1       5717    6451    2838    2112    0.0     1290
::chromosome_1:62281-74713      evm.model.chromosome_1.5        100.000 204     0       0       11093   11296   1151    1354    1.29e-102       377
```

lncRNA blasting to smRNA sequences: 

```
Number of hits?
0 /data/putnamlab/jillashey/Astrangia2021/lncRNA/blast/blastn_smRNA_lncRNA_query.tab
```

Maybe I should blast to the collapsed smRNA file instead of the mature miRNA file. I'm not going to do that for now but I am going to download the blast lncRNA results to my local computer. I only moved the mRNA results to my github and the others to my local computer, as they were too large to be stored on github. I need to make an OSF repo and move the larger files there.  






Possible lncRNA-mRNA interaction software

- [RNAplex](https://academic.oup.com/bioinformatics/article/24/22/2657/184477)
- [lncTar](http://www.cuilab.cn/lnctar)
- RNAhybrid