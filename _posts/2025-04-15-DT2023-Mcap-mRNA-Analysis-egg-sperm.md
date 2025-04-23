---
layout: post
title: Developmental 2023 Timeseries mRNA analysis - egg and sperm analysis 
date: '2025-04-15'
categories: Analysis
tags: [Bioinformatics, mRNA, Montipora capitata]
projects: Developmental Timseries 2023
---

## Developmental 2023 Timeseries mRNA analysis - egg and sperm

I have sequenced and analyzed samples from the following time points: 1, 4, 9, 14, 22, 28, 48 and 72 hpf. The github for this project is [here](https://github.com/JillAshey/Hawaii_Developmental_TimeSeries). During the closed portion of my defense, we discussed incorporating unfertilized egg and sperm samples into the analyses to understand the contribution of the sperm v. egg to the mRNA complement. In my 2023 experiment, I did not collect sperm or unfertilized egg samples so I am going to use the sperm and egg samples that were used in [Van et Etten et al. 2020](https://peerj.com/articles/9739/#supp-1), which analyzed sperm and unfertilized egg samples in Mcap. In the paper, they state the following on accessing the sample fastq files: "The egg data are publicly available under NCBI BioProject [PRJNA616341](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA616341) (SAMN14486762, SAMN14486763, SAMN14486764) and the sperm data are publicly available under NCBI BioProject [PRJNA339779](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA339779)." They also state 

Here are the samples that I need to download from NCBI: 

- [SAMN05607941](https://www.ncbi.nlm.nih.gov/sra/SRX2039373[accn]) - sperm
- [SAMN14486762](https://www.ncbi.nlm.nih.gov/sra/?term=SAMN14486762) - egg
- [SAMN14486763](https://www.ncbi.nlm.nih.gov/sra/?term=SAMN14486763) - egg 
- [SAMN14486764](https://www.ncbi.nlm.nih.gov/sra/?term=SAMN14486764) - egg

The egg and sperm libraries were generated using different methods, which is confusing/annoying. 

|                  | Sperm                                                           | Egg                                                                 |
| ---------------- | --------------------------------------------------------------- | ------------------------------------------------------------------- |
| Library prep kit | Illumina TruSeq RNA Library Prep Kit v2                         | Standard Illumina strand-specific RNA-seq prep with polyA selection |
| Sequencer        | Illumina MiSeq flowcell using the Illumina MiSeq Reagent Kit v3 | Illumina HiSeq                                                      |
| Configuation     | Single end                                                      | Paired end                                                          |

To download these sequences from NCBI, I need to run SRA toolkit on Unity (see [example](https://github.com/zdellaert/multi-sp-snRNA/blob/main/notes/mitochondrial-reads.md) from ZD notebook). 

`nano prefetch.sh`

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=5
#SBATCH --partition=uri-cpu
#SBATCH --no-requeue
#SBATCH --mem=200GB
#SBATCH -t 24:00:00
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error
#SBATCH -D /work/pi_hputnam_uri_edu/jillashey/Mcap_2023/polyA/scripts

module load uri/main
module load SRA-Toolkit/3.0.3-gompi-2022a

cd /project/pi_hputnam_uri_edu/jillashey/Mcap_2023/egg_sperm/

prefetch --max-size 30GB SRR4048723 # sperm sample
prefetch --max-size 30GB SRR11452263 # egg sample
prefetch --max-size 30GB SRR11452262 # egg sample
prefetch --max-size 30GB SRR11452251 # egg sample
```

Submitted batch job 33030565. Downloaded successfully. I now need to convert the sra file to fastq also with SRA toolkit. 

`nano fasterq.sh`

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=5
#SBATCH --partition=uri-cpu
#SBATCH --no-requeue
#SBATCH --mem=200GB
#SBATCH -t 24:00:00
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error
#SBATCH -D /work/pi_hputnam_uri_edu/jillashey/Mcap_2023/polyA/scripts

module load uri/main
module load SRA-Toolkit/3.0.3-gompi-2022a

cd /project/pi_hputnam_uri_edu/jillashey/Mcap_2023/egg_sperm/

fasterq-dump SRR11452251
fasterq-dump SRR11452262
fasterq-dump SRR11452263
fasterq-dump SRR4048723
```

Submitted batch job 33035417.

Make directories on Unity. Sym link the raw data in the project directory to the work directory.

```
cd /work/pi_hputnam_uri_edu/jillashey/Mcap_2023/egg_sperm/data/raw
ln -s /project/pi_hputnam_uri_edu/jillashey/Mcap_2023/egg_sperm/*fastq .
```

Run fastqc on the raw reads. `nano raw_fastqc.sh`

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=5
#SBATCH --signal=2
#SBATCH --no-requeue
#SBATCH --mem=100GB
#SBATCH -t 12:00:00
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error
#SBATCH -D /work/pi_hputnam_uri_edu/jillashey/Mcap_2023/egg_sperm/scripts 

# load modules needed
module load parallel/20240822
module load fastqc/0.12.1
module load uri/main
module load all/MultiQC/1.12-foss-2021b

cd /work/pi_hputnam_uri_edu/jillashey/Mcap_2023/egg_sperm/data/raw

# Create an array of fastq files to process
files=($('ls' *.fastq)) 

# Run fastqc in parallel
echo "Starting fastqc..." $(date)
parallel -j 20 "fastqc {} -o /work/pi_hputnam_uri_edu/jillashey/Mcap_2023/egg_sperm/output/fastqc/raw && echo 'Processed {}'" ::: "${files[@]}"
echo "fastQC done." $(date)

cd /work/pi_hputnam_uri_edu/jillashey/Mcap_2023/egg_sperm/output/fastqc/raw

echo "Starting multiqc..." $(date)
multiqc *

echo "Initial QC of egg/sperm data complete." $(date)
```

Submitted batch job 33036909. Data is super clean! I need to do some adapter trimming but other than that, its very high quality. 

`nano trim_cutadapt.sh`

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=5
#SBATCH --partition=uri-cpu
#SBATCH --no-requeue
#SBATCH --mem=100GB
#SBATCH -t 12:00:00
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error
#SBATCH -D /work/pi_hputnam_uri_edu/jillashey/Mcap_2023/egg_sperm/scripts 

# load modules needed
module load parallel/20240822
module load fastqc/0.12.1
module load uri/main
module load all/cutadapt/3.5-GCCcore-11.2.0
module load all/MultiQC/1.12-foss-2021b

cd /work/pi_hputnam_uri_edu/jillashey/Mcap_2023/egg_sperm/data/raw

# Adapter sequence (Illumina universal)
ADAPTER=AGATCGGAAGAGC

echo "Trim egg samples (PE)" $(date)
# Paired-end files
for ID in SRR11452251 SRR11452262 SRR11452263; do
  cutadapt -a $ADAPTER -A $ADAPTER  -q 20,20 --minimum-length=20 -o ${ID}_1_AdapterTrimmed.fastq -p ${ID}_2_AdapterTrimmed.fastq ${ID}_1.fastq ${ID}_2.fastq
done

echo "Trim sperm samples (SE)" $(date)
# Single-end file
cutadapt -a $ADAPTER  -q 20,20 --minimum-length=20 -o SRR4048723_AdapterTrimmed.fastq SRR4048723.fastq

mv *AdapterTrimmed.fastq /scratch/workspace/jillashey_uri_edu-ashey_scratch/Mcap2023/egg_sperm_trim
cd /scratch/workspace/jillashey_uri_edu-ashey_scratch/Mcap2023/egg_sperm_trim

# Create an array of fastq files to process
files=($('ls' *.fastq)) 

# Run fastqc in parallel
echo "Starting fastqc..." $(date)
parallel -j 20 "fastqc {} -o /work/pi_hputnam_uri_edu/jillashey/Mcap_2023/egg_sperm/output/fastqc/trim && echo 'Processed {}'" ::: "${files[@]}"
echo "fastQC done." $(date)

cd /work/pi_hputnam_uri_edu/jillashey/Mcap_2023/egg_sperm/output/fastqc/trim

echo "Starting multiqc..." $(date)
multiqc *
echo "MultiQC complete..." $(date)
```

Submitted batch job 33041964





The egg and sperm samples were prepared in different ways so I will have to trim them in different ways. In the Van Etten et al. paper, they used CLC Genomics Workbench to trim. 


