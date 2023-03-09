---
layout: post
title: Point Judith Oyster Nutrient DNA methylation 
date: '2023-03-06'
categories: Analysis
tags: [MBD-BS, Bioinformatics, oyster, nutrients, DNA, methylation]
projects: Point Judith Oyster
---

## Point Judith oyster DNA methylation

As a lab, we have been working on an oyster paper that details the responses of Eastern oysters to differing nutrient regimes. [Here](https://github.com/hputnam/Cvir_Nut_Int) is the project github page. 

ES ran the methylation data through the [nf-core methylseq](https://nf-co.re/methylseq/usage) pipeline, but there are some discrepancies in the QC outputs that we need to resolve. Methylseq doesn't appear to be recognizing the adapter bp and so isn't properly trimming them. [Here](https://github.com/emmastrand/EmmaStrand_Notebook/blob/master/_posts/2022-05-09-Point-Judith-Oyster-DNA-Methylation-(MBD-BS).md) is ES notebook post detailing her methylseq runs. To address this, I will run the methylseq pipeline but run each step individually (as opposed to the steps being run in the pipeline wrapper). 

20230306

### Set up data in my directory 

Original data lives here: `/data/putnamlab/KITT/hputnam/20200119_Oyst_Nut/MBDBS`

#### Make new folder in my own directory 

```
cd /data/putnamlab/jillashey

mkdir Oys_Nutrient
cd Oys_Nutrient

mkdir MBDBS RNASeq
cd MBDBS

mkdir data scripts fastqc_results
cd data 
mkdir raw trim

cd ../fastqc_results
mkdir raw trim 
```

#### Copy data over into my directory 

There are 24 total files, but 12 total samples (each sample has a R1 and R2 file).

```
cd /data/putnamlab/KITT/hputnam/20200119_Oyst_Nut/MBDBS
cp *fastq.gz /data/putnamlab/jillashey/Oys_Nutrient/MBDBS/data/raw

cd /data/putnamlab/jillashey/Oys_Nutrient/MBDBS/data/raw
```

#### Count number of reads per file 

Number of reads can be counted by counting the number of '@__' terms there are in the file, as each one corresponds to one read.

```
zgrep -c "@M00" *fastq.gz

HPB10_S44_L001_R1_001.fastq.gz:1036453
HPB10_S44_L001_R2_001.fastq.gz:1036453
HPB11_S45_L001_R1_001.fastq.gz:976844
HPB11_S45_L001_R2_001.fastq.gz:976844
HPB12_S46_L001_R1_001.fastq.gz:1541665
HPB12_S46_L001_R2_001.fastq.gz:1541665
HPB1_S35_L001_R1_001.fastq.gz:965843
HPB1_S35_L001_R2_001.fastq.gz:965843
HPB2_S36_L001_R1_001.fastq.gz:990538
HPB2_S36_L001_R2_001.fastq.gz:990538
HPB3_S37_L001_R1_001.fastq.gz:953246
HPB3_S37_L001_R2_001.fastq.gz:953246
HPB4_S38_L001_R1_001.fastq.gz:864752
HPB4_S38_L001_R2_001.fastq.gz:864752
HPB5_S39_L001_R1_001.fastq.gz:923141
HPB5_S39_L001_R2_001.fastq.gz:923141
HPB6_S40_L001_R1_001.fastq.gz:909119
HPB6_S40_L001_R2_001.fastq.gz:909119
HPB7_S41_L001_R1_001.fastq.gz:1053144
HPB7_S41_L001_R2_001.fastq.gz:1053144
HPB8_S42_L001_R1_001.fastq.gz:923398
HPB8_S42_L001_R2_001.fastq.gz:923398
HPB9_S43_L001_R1_001.fastq.gz:1067829
HPB9_S43_L001_R2_001.fastq.gz:1067829
```

### Run fastqc on raw data 

In scripts folder: `nano fastqc_raw.sh`

```
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Oys_Nutrient/MBDBS/scripts              
#SBATCH --error="fastqc_raw_error" #if your job fails, the error report will be put in this file
#SBATCH --output="fastqc_raw_output" #once your job is completed, any final job report comments will be put in this file

source /usr/share/Modules/init/sh # load the module function

cd /data/putnamlab/jillashey/Oys_Nutrient/MBDBS

module load FastQC/0.11.9-Java-11
module load MultiQC/1.9-intel-2020a-Python-3.8.2

for file in /data/putnamlab/jillashey/Oys_Nutrient/MBDBS/data/raw/*fastq.gz
do 
fastqc $file --outdir /data/putnamlab/jillashey/Oys_Nutrient/MBDBS/fastqc_results/raw
done

multiqc --interactive fastqc_results
```

`sbatch fastqc_raw.sh`; Submitted batch job 239141

#### Copy multiQC report to Cvir github and look at results 

On local computer: 

```
scp jillashey@ssh3.hac.uri.edu:/data/putnamlab/jillashey/Oys_Nutrient/MBDBS/multiqc_report.html /Users/jillashey/Desktop/PutnamLab/Repositories/Cvir_Nut_Int/output/MBDBS/JA
```

- This is ES initial report: https://github.com/hputnam/Cvir_Nut_Int/blob/master/output/MBDBS/initial_multiqc_report.html
- This is my initial report: https://github.com/hputnam/Cvir_Nut_Int/blob/master/output/MBDBS/JA/raw_multiqc_report.html 




Mine & ES report seem to be the same. That's a good sign (I believe)! Now can move onto trimming w/ Trim Galore.

### Trim sequences 

Going to trim w/ [Trim Galore](https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md), as this is what methylseq uses for trimming purposes. 

Two versions of trim galore on URI Andromeda: `all/Trim_Galore/0.6.6-GCC-10.2.0-Python-2.7.18` and `all/Trim_Galore/0.6.7-GCCcore-11.2.0`. I'm going to use the newest one. 

I'm going to do several iterations of trimming so I'm going to make different trimming folders in my data and fastqc_results folder.

```
cd /data/putnamlab/jillashey/Oys_Nutrient/MBDBS/fastqc_results
mkdir trim1 trim2 trim3 trim4

cd /data/putnamlab/jillashey/Oys_Nutrient/MBDBS/data
mkdir trim1 trim2 trim3 trim4
```

##### Attempt 1

I'm not going to add any sort of cutoffs for this initial trimming run. 

In scripts folder: `nano trim_galore1.sh`

```
#!/bin/bash
#SBATCH -t 48:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Oys_Nutrient/MBDBS/scripts              
#SBATCH --error="trim_galore1_error" #if your job fails, the error report will be put in this file
#SBATCH --output="trim_galore1_output" #once your job is completed, any final job report comments will be put in this file

source /usr/share/Modules/init/sh # load the module function

cd /data/putnamlab/jillashey/Oys_Nutrient/MBDBS/data/raw

module load Trim_Galore/0.6.7-GCCcore-11.2.0

for file in "HPB10_S44" "HPB11_S45" "HPB12_S46" "HPB1_S35" "HPB2_S36" "HPB3_S37" "HPB4_S38" "HPB5_S39" "HPB6_S40" "HPB7_S41" "HPB8_S42" "HPB9_S43"
do 
trim_galore --paired ${file}_L001_R1_001.fastq.gz ${file}_L001_R2_001.fastq.gz -o /data/putnamlab/jillashey/Oys_Nutrient/MBDBS/data/trim1/
done
```

`sbatch trim_galore1.sh`; Submitted batch job 240029

For each read, the script generates a new/trimmed .fq.gz file and a trimming report .txt file. Here's an example of what the .txt file looks like for sample `HPB10_S44_L001_R1_001`: 

```
Using Illumina adapter for trimming (count: 804887). Second best hit was Nextera (count: 4)
Adapter sequence: 'AGATCGGAAGAGC' (Illumina TruSeq, Sanger iPCR; auto-detected)
Maximum trimming error rate: 0.1 (default)
Minimum required adapter overlap (stringency): 1 bp
Minimum required sequence length for both reads before a sequence pair gets removed: 20 bp
Output file will be GZIP compressed


This is cutadapt 3.5 with Python 3.9.6
Command line parameters: -j 1 -e 0.1 -q 20 -O 1 -a AGATCGGAAGAGC HPB10_S44_L001_R1_001.fastq.gz
Processing reads on 1 core in single-end mode ...
Finished in 77.59 s (75 µs/read; 0.80 M reads/minute).

=== Summary ===

Total reads processed:               1,036,453
Reads with adapters:                   984,111 (94.9%)
Reads written (passing filters):     1,036,453 (100.0%)

Total basepairs processed:   311,706,539 bp
Quality-trimmed:              51,156,210 bp (16.4%)
Total written (filtered):    173,587,926 bp (55.7%)

=== Adapter 1 ===

Sequence: AGATCGGAAGAGC; Type: regular 3'; Length: 13; Trimmed: 984111 times

Minimum overlap: 1
No. of allowed errors:
1-9 bp: 0; 10-13 bp: 1

Bases preceding removed adapters:
  A: 32.8%
  C: 12.2%
  G: 21.1%
  T: 33.8%
  none/other: 0.0%
```

##### Attempt 2

In ES notebook [post](https://github.com/emmastrand/EmmaStrand_Notebook/blob/master/_posts/2022-05-09-Point-Judith-Oyster-DNA-Methylation-(MBD-BS).md), she ran a couple of different iterations of methylseq where she used different cutoffs for trimming. 

- `PJ_methylseq1.sh: --clip_r1 10 \ --clip_r2 10 \ --three_prime_clip_r1 10 --three_prime_clip_r2 10 \`
	- From ES notebook: Run this first to assess m-bias and then decide if we need more trial runs. See Emma Strand and Kevin Wong's notebook posts for methylation scripts and how they dealt with these issues.
- `PJ_methylseq2.sh: --clip_r1 10 \ --clip_r2 10 \ --three_prime_clip_r1 20 --three_prime_clip_r2 20 \`
	- From ES notebook: Because the run time seemed oddly short, I wanted to run this again to be sure it worked.
- `PJ_methylseq3.sh: --clip_r1 150 \ --clip_r2 150 \ --three_prime_clip_r1 150 --three_prime_clip_r2 150 \`
	- From ES notebook: We think that methylseq isn't properly recognizing the adapters because 1.) the multiqc report says 300 bp length when usually for methylation data we use 2x150 bp. 2.) If the program isn't trimming the adapter fully that would explain the extreme m-bias, adapter content, and decreased quality all after 150 bp length. 3.) Usually Illumina adapters are ~150 bp of the read length and are not included in multiqc reports (if recognized appropriately).
I'm cutting 150 here to see if this makes a difference but unsure if we trust methylseq if it's not recognizing the adapters (i.e. what other problems might we have and don't see yet?)

Parameters 

- `--clip_R1` removes XX bps from the 5' end of read 1. 
- `--clip_R2` removes XX bps from the 5' end of read 2. 
- `--three_prime_clip_R1` removes XX bp after the 3' end of read 1 AFTER adapter/quality trimming has been performed
- `--three_prime_clip_R2` removes XX bp after the 3' end of read 2 AFTER adapter/quality trimming has been performed

Running her first iteration (`--clip_r1 10 \ --clip_r2 10 \ --three_prime_clip_r1 10 --three_prime_clip_r2 10 \`). I'm clipping 10 bp from the 5' end of reads 1 and 2. After adapter/quality trimming has been done, 10 bp will also be clipped from the 3' end of reads 1 and 2. 

In scripts folder: `nano trim_galore2.sh`

```
#!/bin/bash
#SBATCH -t 48:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Oys_Nutrient/MBDBS/scripts              
#SBATCH --error="trim_galore2_error" #if your job fails, the error report will be put in this file
#SBATCH --output="trim_galore2_output" #once your job is completed, any final job report comments will be put in this file

source /usr/share/Modules/init/sh # load the module function

cd /data/putnamlab/jillashey/Oys_Nutrient/MBDBS/data/raw

module load Trim_Galore/0.6.7-GCCcore-11.2.0

for file in "HPB10_S44" "HPB11_S45" "HPB12_S46" "HPB1_S35" "HPB2_S36" "HPB3_S37" "HPB4_S38" "HPB5_S39" "HPB6_S40" "HPB7_S41" "HPB8_S42" "HPB9_S43"
do 
trim_galore --paired ${file}_L001_R1_001.fastq.gz ${file}_L001_R2_001.fastq.gz --clip_r1 10 --clip_r2 10 --three_prime_clip_r1 10 --three_prime_clip_r2 10 -o /data/putnamlab/jillashey/Oys_Nutrient/MBDBS/data/trim2/
done
```

`sbatch trim_galore2.sh`; Submitted batch job 241069

##### Attempt 3

Running her second iteration (`--clip_r1 10 \ --clip_r2 10 \ --three_prime_clip_r1 20 --three_prime_clip_r2 20 \`)

In scripts folder: `nano trim_galore3.sh`

```
#!/bin/bash
#SBATCH -t 48:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Oys_Nutrient/MBDBS/scripts              
#SBATCH --error="trim_galore3_error" #if your job fails, the error report will be put in this file
#SBATCH --output="trim_galore3_output" #once your job is completed, any final job report comments will be put in this file

source /usr/share/Modules/init/sh # load the module function

cd /data/putnamlab/jillashey/Oys_Nutrient/MBDBS/data/raw

module load Trim_Galore/0.6.7-GCCcore-11.2.0

for file in "HPB10_S44" "HPB11_S45" "HPB12_S46" "HPB1_S35" "HPB2_S36" "HPB3_S37" "HPB4_S38" "HPB5_S39" "HPB6_S40" "HPB7_S41" "HPB8_S42" "HPB9_S43"
do 
trim_galore --paired ${file}_L001_R1_001.fastq.gz ${file}_L001_R2_001.fastq.gz --clip_r1 10 --clip_r2 10 --three_prime_clip_r1 20 --three_prime_clip_r2 20 -o /data/putnamlab/jillashey/Oys_Nutrient/MBDBS/data/trim3/
done
```

`sbatch trim_galore3.sh`; Submitted batch job 241074

##### Attempt 4

Running her third iteration: `--clip_r1 150 \ --clip_r2 150 \ --three_prime_clip_r1 150 --three_prime_clip_r2 150 \`

In scripts folder: `nano trim_galore4.sh`

```
#!/bin/bash
#SBATCH -t 48:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Oys_Nutrient/MBDBS/scripts              
#SBATCH --error="trim_galore4_error" #if your job fails, the error report will be put in this file
#SBATCH --output="trim_galore4_output" #once your job is completed, any final job report comments will be put in this file

source /usr/share/Modules/init/sh # load the module function

cd /data/putnamlab/jillashey/Oys_Nutrient/MBDBS/data/raw

module load Trim_Galore/0.6.7-GCCcore-11.2.0

for file in "HPB10_S44" "HPB11_S45" "HPB12_S46" "HPB1_S35" "HPB2_S36" "HPB3_S37" "HPB4_S38" "HPB5_S39" "HPB6_S40" "HPB7_S41" "HPB8_S42" "HPB9_S43"
do 
trim_galore --paired ${file}_L001_R1_001.fastq.gz ${file}_L001_R2_001.fastq.gz --clip_r1 150 --clip_r2 150 --three_prime_clip_r1 150 --three_prime_clip_r2 150 -o /data/putnamlab/jillashey/Oys_Nutrient/MBDBS/data/trim4/
done
```

`sbatch trim_galore4.sh`; Submitted batch job 241075

This failed - error file said that "The 5' clipping value for read 1 should have a sensible value (> 0 and < read length)". So the value should be greater than 0 and less than the read length (300?). But isn't our read length 300 and 150 is less than 300? Changing `clip_r1`, `clip_r2`, `three_prime_clip_r1`, and `three_prime_clip_r2` to 100 instead of 150 and seeing how that works. Submitted batch job 241088. Failed agian, same error as before.

In the trim galore source [code](https://github.com/FelixKrueger/TrimGalore/blob/master/trim_galore), it looks like trimming is only possible if its between 0 and less than 100. See code here: 

```
  if (defined $clip_r1){ # trimming 5' bases of read 1
      unless ($clip_r1 > 0 and $clip_r1 < 100){
      die "The 5' clipping value for read 1 should have a sensible value (> 0 and < read length)\n\n";
      }
  }

  if (defined $clip_r2){ # trimming 5' bases of read 2
      unless ($clip_r2 > 0 and $clip_r2 < 100){
      die "The 5' clipping value for read 2 should have a sensible value (> 0 and < read length)\n\n";
      }
  }

  ### Trimming at the 3' end
  if (defined $three_prime_clip_r1){ # trimming 3' bases of read 1
      unless ($three_prime_clip_r1 > 0 and $three_prime_clip_r1 < 100){
      die "The 3' clipping value for read 1 should have a sensible value (> 0 and < read length)\n\n";
      }
  }

  if (defined $three_prime_clip_r2){ # trimming 3' bases of read 2
      unless ($three_prime_clip_r2 > 0 and $three_prime_clip_r2 < 100){
      die "The 3' clipping value for read 2 should have a sensible value (> 0 and < read length)\n\n";
      }
```

So I don't think I can trim more than 100 bp. Let's try 99 bp in trim_galore4.sh just to see if it will work. Submitted batch job 241090.

##### Attempt 5 

Based on the fastqc results, it looks like the 

### Run Fastqc on trimmed data 

##### Attempt 1 

In scripts folder: `nano fastqc_trim1.sh`

```
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Oys_Nutrient/MBDBS/scripts              
#SBATCH --error="fastqc_trim1_error" #if your job fails, the error report will be put in this file
#SBATCH --output="fastqc_trim1_output" #once your job is completed, any final job report comments will be put in this file

source /usr/share/Modules/init/sh # load the module function

cd /data/putnamlab/jillashey/Oys_Nutrient/MBDBS

module load FastQC/0.11.9-Java-11
module load MultiQC/1.9-intel-2020a-Python-3.8.2

for file in /data/putnamlab/jillashey/Oys_Nutrient/MBDBS/data/trim1/*fq.gz
do 
fastqc $file --outdir /data/putnamlab/jillashey/Oys_Nutrient/MBDBS/fastqc_results/trim1
done

multiqc --interactive fastqc_results/trim1
```

`sbatch fastqc_trim1.sh`; Submitted batch job 241058

##### Attempt 2

In scripts folder: `nano fastqc_trim2.sh`

```
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Oys_Nutrient/MBDBS/scripts              
#SBATCH --error="fastqc_trim2_error" #if your job fails, the error report will be put in this file
#SBATCH --output="fastqc_trim2_output" #once your job is completed, any final job report comments will be put in this file

source /usr/share/Modules/init/sh # load the module function

cd /data/putnamlab/jillashey/Oys_Nutrient/MBDBS

module load FastQC/0.11.9-Java-11
module load MultiQC/1.9-intel-2020a-Python-3.8.2

for file in /data/putnamlab/jillashey/Oys_Nutrient/MBDBS/data/trim2/*fq.gz
do 
fastqc $file --outdir /data/putnamlab/jillashey/Oys_Nutrient/MBDBS/fastqc_results/trim2
done

multiqc --interactive fastqc_results/trim2
```

`sbatch fastqc_trim2.sh`; Submitted batch job 241100

##### Attempt 3

In scripts folder: `nano fastqc_trim3.sh`

```
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Oys_Nutrient/MBDBS/scripts              
#SBATCH --error="fastqc_trim3_error" #if your job fails, the error report will be put in this file
#SBATCH --output="fastqc_trim3_output" #once your job is completed, any final job report comments will be put in this file

source /usr/share/Modules/init/sh # load the module function

cd /data/putnamlab/jillashey/Oys_Nutrient/MBDBS

module load FastQC/0.11.9-Java-11
module load MultiQC/1.9-intel-2020a-Python-3.8.2

for file in /data/putnamlab/jillashey/Oys_Nutrient/MBDBS/data/trim3/*fq.gz
do 
fastqc $file --outdir /data/putnamlab/jillashey/Oys_Nutrient/MBDBS/fastqc_results/trim3
done

multiqc --interactive fastqc_results/trim3
```

`sbatch fastqc_trim3.sh`; Submitted batch job 241125

##### Attempt 4

In scripts folder: `nano fastqc_trim4.sh`

```
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Oys_Nutrient/MBDBS/scripts              
#SBATCH --error="fastqc_trim4_error" #if your job fails, the error report will be put in this file
#SBATCH --output="fastqc_trim4_output" #once your job is completed, any final job report comments will be put in this file

source /usr/share/Modules/init/sh # load the module function

cd /data/putnamlab/jillashey/Oys_Nutrient/MBDBS

module load FastQC/0.11.9-Java-11
module load MultiQC/1.9-intel-2020a-Python-3.8.2

for file in /data/putnamlab/jillashey/Oys_Nutrient/MBDBS/data/trim4/*fq.gz
do 
fastqc $file --outdir /data/putnamlab/jillashey/Oys_Nutrient/MBDBS/fastqc_results/trim4
done

multiqc --interactive fastqc_results/trim4
```

`sbatch fastqc_trim4.sh`; Submitted batch job 241148