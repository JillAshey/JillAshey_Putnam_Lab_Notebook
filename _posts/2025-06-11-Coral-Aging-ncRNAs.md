---
layout: post
title: Aging ncRNAs 
date: '2025-06-11'
categories: Analysis
tags: [Bioinformatics]
projects: 
---

## ncRNA characterization for coral aging 

Emma and I have been talking about using ncRNA expression as a 'clock' to age a coral. We decided to move forward with an initial analysis of Mcap (since we have so many samples across life stages) RNAseq data. For this first pass, we are going to characterize lncRNAs (and maybe smRNAs) in Mcap samples from Emma's work (adults) and my work (larvae). I'm not sure if we will be able to derive smRNAs from RNAseq data, but we can try! 

All of our sequenced samples are on Unity in various places. Here are the paths for the data we will use: 

```
/project/pi_hputnam_uri_edu/raw_sequencing_data/20220203_BleachedPairs_RNASeq
/project/pi_hputnam_uri_edu/raw_sequencing_data/20240328_Mcap_RNASeq_Devo
/project/pi_hputnam_uri_edu/raw_sequencing_data/20240920_Ashey_Mcap_Devo
/project/pi_hputnam_uri_edu/raw_sequencing_data/20191104_HoloInt/30-274849628
/project/pi_hputnam_uri_edu/raw_sequencing_data/20200404_HoloInt_Batch3
```

First, I am going to sym link all the folders to one central folder. This might take a while so will do as a job. In the scripts folder: `nano symlinks.sh`

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
#SBATCH -D /work/pi_hputnam_uri_edu/jillashey/coral_aging/scripts

# Path to your TSV file
FILE_LIST="/work/pi_hputnam_uri_edu/jillashey/coral_aging/mcap_rna_seq_files.tsv"

# Directory where you want to place all symlinks
LINK_DIR="/work/pi_hputnam_uri_edu/jillashey/coral_aging/data/raw"

# Skip header and read file
tail -n +2 "$FILE_LIST" | while IFS=$'\t' read -r R1 R2 DIR; do
    # Remove any trailing carriage return characters
    R1=$(echo "$R1" | tr -d '\r')
    R2=$(echo "$R2" | tr -d '\r')
    DIR=$(echo "$DIR" | tr -d '\r')

    SRC1="${DIR}/${R1}"
    SRC2="${DIR}/${R2}"
    
    ln -sf "$SRC1" "$LINK_DIR/$(basename "$R1")"
    ln -sf "$SRC2" "$LINK_DIR/$(basename "$R2")"
done
```

Submitted batch job 38064879. Ran super fast lol. QC the data. In the scripts folder: `nano raw_qc.sh`

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
#SBATCH -D /work/pi_hputnam_uri_edu/jillashey/coral_aging/scripts

# load modules needed
module load fastqc/0.12.1
module load uri/main
module load all/MultiQC/1.12-foss-2021b

# Set directory where all fastq.gz files are located
DATA_DIR="/work/pi_hputnam_uri_edu/jillashey/coral_aging/data/raw"

# Output directory for FastQC results
FASTQC_OUT="/work/pi_hputnam_uri_edu/jillashey/coral_aging/output/fastqc/raw"

echo "Running FastQC on all FASTQ files in $DATA_DIR..."

# Run FastQC on all .fastq.gz files
fastqc -o "$FASTQC_OUT" "$DATA_DIR"/*.fastq.gz

echo "FastQC complete. Running MultiQC..."

# Move into the FastQC output directory and run MultiQC
cd "$FASTQC_OUT" || exit
multiqc .

echo "MultiQC report generated in $FASTQC_OUT"
```

Submitted batch job 38065181. Raw QC is [here](https://github.com/JillAshey/Coral_Aging/blob/main/data/mcap/raw_mcap_age_multiqc_report.html). Overall, the data looks quite good. There are some differences in sequencing depth, as well as read length. We will keep an eye on these as we move forward. The primary issue with the raw reads is the adapter content. Emma and I both used fastp to trim so I am going to stick with that. Also going to make a scratch directory to store the trimmed fastq files and bam files. 

```
ws_allocate -n coral_age -G pi_hputnam_uri_edu -d 30 -r 5 -m jillashey@uri.edu
Info: creating workspace.
/scratch3/workspace/jillashey_uri_edu-coral_age
remaining extensions  : 5
remaining time in days: 30

cd /scratch3/workspace/jillashey_uri_edu-coral_age
mkdir trim bam
```

In the scripts folder: `nano trim_and_qc.sh`

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
#SBATCH -D /work/pi_hputnam_uri_edu/jillashey/coral_aging/scripts

# load modules needed
module load fastqc/0.12.1
module load uri/main
module load all/MultiQC/1.12-foss-2021b
module load fastp/0.23.2-GCC-11.2.0

# Set paths
RAW_DIR="/work/pi_hputnam_uri_edu/jillashey/coral_aging/data/raw"
TRIM_DIR="/scratch3/workspace/jillashey_uri_edu-coral_age/trim"
FASTQC_DIR="/data/putnamlab/jillashey/DT_Mcap_2023/mRNA/output/fastqc/trim"
MULTIQC_DIR="/data/putnamlab/jillashey/DT_Mcap_2023/mRNA/output/fastqc/trim"

echo "Starting trimming and QC pipeline..."

# Get all R1 files
for R1 in "$RAW_DIR"/*_R1_001.fastq.gz; do
    BASE=$(basename "$R1" _R1_001.fastq.gz)
    R2="$RAW_DIR/${BASE}_R2_001.fastq.gz"

    OUT1="$TRIM_DIR/trim.${BASE}_R1_001.fastq.gz"
    OUT2="$TRIM_DIR/trim.${BASE}_R2_001.fastq.gz"

    echo "--------------------------------------------"
    echo "Processing sample: $BASE"
    echo "Input R1: $R1"
    echo "Input R2: $R2"
    echo "Output R1: $OUT1"
    echo "Output R2: $OUT2"

    # Run fastp
    echo "Running fastp..."
    fastp \
        --in1 "$R1" \
        --in2 "$R2" \
        --out1 "$OUT1" \
        --out2 "$OUT2" \
        --detect_adapter_for_pe \
        --qualified_quality_phred 30 \
        --unqualified_percent_limit 10 \
        #--length_required 100 \
        --cut_right \
        --cut_right_window_size 5 \
        --cut_right_mean_quality 20 
    echo "Finished fastp for $BASE."

    # Run fastqc
    echo "Running FastQC..."
    fastqc -o "$FASTQC_DIR" "$OUT1" "$OUT2"
    echo "Finished FastQC for $BASE."

done

# Run multiqc
echo "Running MultiQC on all FastQC outputs..."
multiqc "$FASTQC_DIR" -o "$MULTIQC_DIR"
echo "MultiQC complete. All done!"
```

Submitted batch job 38139904
