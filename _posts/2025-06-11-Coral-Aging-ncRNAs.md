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
FASTQC_DIR="/work/pi_hputnam_uri_edu/jillashey/coral_aging/output/fastqc/trim"
MULTIQC_DIR="/work/pi_hputnam_uri_edu/jillashey/coral_aging/output/fastqc/trim"

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

trim batch job 38139904. QC batch job 38154121. Trimmed QC is [here](https://github.com/JillAshey/Coral_Aging/blob/main/data/mcap/trim_mcap_age_multiqc_report.html). Trimming looks great!! 

Align trimmed reads. In the scripts folder: `nano align.sh`

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=2
#SBATCH --partition=uri-cpu
#SBATCH --no-requeue
#SBATCH --mem=200GB
#SBATCH -t 100:00:00
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error
#SBATCH -D /work/pi_hputnam_uri_edu/jillashey/coral_aging/scripts

# Load modules 
module load uri/main
module load HISAT2/2.2.1-gompi-2022a
module load SAMtools/1.18-GCC-12.3.0

# Set directories
TRIM_DIR="/scratch3/workspace/jillashey_uri_edu-coral_age/trim"
BAM_DIR="/scratch3/workspace/jillashey_uri_edu-coral_age/bam"
HISAT2_INDEX="/work/pi_hputnam_uri_edu/HI_Genomes/MCapV3/McapV3_hisat2_ref" 

## hisat2 reference for Mcap V3 already built and in HI_Genomes folder 

# Get list of trimmed R1 files
array=($(ls "$TRIM_DIR"/trim.*_R1_001.fastq.gz))

# Alignment loop
for R1 in "${array[@]}"; do
    # Get basename (e.g., trim.sample -> sample)
    basefile=$(basename "$R1")
    sample_name=$(echo "$basefile" | cut -d '.' -f2)

    R2="${R1/_R1_/_R2_}"

    echo "Aligning sample: $sample_name"
    echo "  Input R1: $R1"
    echo "  Input R2: $R2"

    # Run HISAT2 alignment
    hisat2 -p 8 --rna-strandness RF --dta \
        -x "$HISAT2_INDEX" \
        -1 "$R1" -2 "$R2" \
        -S "${sample_name}.sam"

    # Convert to sorted BAM
    samtools sort -@ 8 -o "$BAM_DIR/${sample_name}.bam" "${sample_name}.sam"
    echo "  -> ${sample_name}.bam created in $BAM_DIR"

    # Remove intermediate SAM
    rm "${sample_name}.sam"
done

echo "Alignment complete!" $(date)

# Mapping summary
SUMMARY_FILE="$BAM_DIR/mapped_reads_counts_Mcap.txt"
echo "Sample Mapping Summary - $(date)" > "$SUMMARY_FILE"

for BAM in "$BAM_DIR"/*.bam; do
    echo "$(basename "$BAM")" >> "$SUMMARY_FILE"
    samtools flagstat "$BAM" | grep "mapped (" >> "$SUMMARY_FILE"
done

echo "Mapping summary saved to $SUMMARY_FILE"
```

Submitted batch job 38207703. High mapping percentages for all samples (except for a few early life stage samples from my work, this is to be expected). See mapping output [here](https://github.com/JillAshey/Coral_Aging/blob/main/data/mcap/mapped_reads_counts_Mcap.txt). 

Run stringtie to assemble reads. In the scripts folder: `nano assemble.sh`

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=2
#SBATCH --partition=uri-cpu
#SBATCH --no-requeue
#SBATCH --mem=200GB
#SBATCH -t 100:00:00
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error
#SBATCH -D /work/pi_hputnam_uri_edu/jillashey/coral_aging/scripts

# Load modules 
module load uri/main
module load StringTie/2.2.1-GCC-11.2.0

# Set directories
BAM_DIR="/scratch3/workspace/jillashey_uri_edu-coral_age/bam"
STRING_DIR="/scratch3/workspace/jillashey_uri_edu-coral_age/stringtie"
GFF="/work/pi_hputnam_uri_edu/HI_Genomes/MCapV3/Montipora_capitata_HIv3.genes_fixed.gff3"

# Get list of bam files 
array=($(ls "$BAM_DIR"/*bam))

echo "Assembling transcripts using stringtie" $(date)

# Loop through and run StringTie
for i in "${array[@]}"; do
    sample_name=$(basename "$i" .bam)

    echo "Running StringTie on $sample_name"

    stringtie "$i" \
        -p 8 \
        -e \
        -B \
        -G "$GFF" \
        -A "$STRING_DIR/${sample_name}.gene_abund.tab" \
        -o "$STRING_DIR/${sample_name}.gtf"

    echo "Done with $sample_name"
done

echo "All transcript assemblies complete!" $(date)
```

Submitted batch job 38291875. 

Once assembly is done, make a list of gtfs. 

```
cd /scratch3/workspace/jillashey_uri_edu-coral_age/stringtie
ls *.gtf > gtf_list.txt
```

Merge gtfs into single gtf to check accuracy

```
module load uri/main
module load StringTie/2.2.1-GCC-11.2.0

stringtie --merge -e -p 8 -G /work/pi_hputnam_uri_edu/HI_Genomes/MCapV3/Montipora_capitata_HIv3.genes_fixed.gff3 -o mcap_age_merged.gtf gtf_list.txt 
```

Use gffcompare to look at accuracy of assembly 

```
module load GffCompare/0.12.6-GCC-11.2.0
gffcompare -r /work/pi_hputnam_uri_edu/HI_Genomes/MCapV3/Montipora_capitata_HIv3.genes_fixed.gff3 -G -o merged mcap_age_merged.gtf 
54384 reference transcripts loaded.
  2 duplicate reference transcripts discarded.
  54384 query transfrags loaded.
```

Look at merge stats 

```
#= Summary for dataset: mcap_age_merged.gtf 
#     Query mRNAs :   54384 in   54185 loci  (36023 multi-exon transcripts)
#            (141 multi-transcript loci, ~1.0 transcripts per locus)
# Reference mRNAs :   54382 in   54185 loci  (36023 multi-exon)
# Super-loci w/ reference transcripts:    54185
#-----------------| Sensitivity | Precision  |
        Base level:   100.0     |   100.0    |
        Exon level:    99.9     |    99.9    |
      Intron level:   100.0     |   100.0    |
Intron chain level:   100.0     |   100.0    |
  Transcript level:    99.9     |    99.9    |
       Locus level:   100.0     |   100.0    |
     Matching intron chains:   36023
       Matching transcripts:   54307
              Matching loci:   54173
          Missed exons:       0/256029  (  0.0%)
           Novel exons:       0/256026  (  0.0%)
        Missed introns:       0/201643  (  0.0%)
         Novel introns:       0/201643  (  0.0%)
           Missed loci:       0/54185   (  0.0%)
            Novel loci:       0/54185   (  0.0%)
 Total union super-loci across all input datasets: 54185 
54384 out of 54384 consensus transcripts written in merged.annotated.gtf (0 discarded as redundant)
```

Assembly looks good! Now on to lncRNA identification (using previous Mcap [workflow](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/_posts/2024-03-28-DT2023-Mcap-mRNA-Analysis.md)). From the merged gtf file, pull out any transcripts that are >199bp (the typical cutoff for lncRNA length). 

```
salloc 
awk '$3 == "transcript" && $1 !~ /^#/' merged.annotated.gtf | awk '($5 - $4 > 199) || ($4 - $5 > 199)' > mcap_age_lncRNA_candidates.gtf

head mcap_age_lncRNA_candidates.gtf
Montipora_capitata_HIv3___Scaffold_1    StringTie       transcript      82397   95409   .       +       .       transcript_id "Montipora_capitata_HIv3___RNAseq.g4584.t1"; gene_id "MSTRG.4"; gene_n
ame "Montipora_capitata_HIv3___RNAseq.g4584.t1"; xloc "XLOC_000001"; ref_gene_id "Montipora_capitata_HIv3___RNAseq.g4584.t1"; cmp_ref "Montipora_capitata_HIv3___RNAseq.g4584.t1"; class_code "="; t
ss_id "TSS1";
Montipora_capitata_HIv3___Scaffold_1    StringTie       transcript      109801  163388  .       +       .       transcript_id "Montipora_capitata_HIv3___RNAseq.g4586.t1"; gene_id "MSTRG.7"; gene_n
ame "Montipora_capitata_HIv3___RNAseq.g4586.t1"; xloc "XLOC_000002"; ref_gene_id "Montipora_capitata_HIv3___RNAseq.g4586.t1"; cmp_ref "Montipora_capitata_HIv3___RNAseq.g4586.t1"; class_code "="; t
ss_id "TSS2";
Montipora_capitata_HIv3___Scaffold_1    StringTie       transcript      162568  163156  .       +       .       transcript_id "Montipora_capitata_HIv3___TS.g26272.t1"; gene_id "MSTRG.7"; gene_name
 "Montipora_capitata_HIv3___TS.g26272.t1"; xloc "XLOC_000003"; ref_gene_id "Montipora_capitata_HIv3___TS.g26272.t1"; cmp_ref "Montipora_capitata_HIv3___TS.g26272.t1"; class_code "="; tss_id "TSS3"
;
Montipora_capitata_HIv3___Scaffold_1    StringTie       transcript      169950  170561  .       +       .       transcript_id "Montipora_capitata_HIv3___TS.g26273.t1"; gene_id "MSTRG.5"; gene_name
 "Montipora_capitata_HIv3___TS.g26273.t1"; xloc "XLOC_000004"; ref_gene_id "Montipora_capitata_HIv3___TS.g26273.t1"; cmp_ref "Montipora_capitata_HIv3___TS.g26273.t1"; class_code "="; tss_id "TSS4"
;
Montipora_capitata_HIv3___Scaffold_1    StringTie       transcript      170982  172082  .       +       .       transcript_id "Montipora_capitata_HIv3___RNAseq.g4588.t1"; gene_id "MSTRG.13"; gene_
name "Montipora_capitata_HIv3___RNAseq.g4588.t1"; xloc "XLOC_000005"; ref_gene_id "Montipora_capitata_HIv3___RNAseq.g4588.t1"; cmp_ref "Montipora_capitata_HIv3___RNAseq.g4588.t1"; class_code "="; 
tss_id "TSS5";
Montipora_capitata_HIv3___Scaffold_1    StringTie       transcript      176400  276376  .       +       .       transcript_id "Montipora_capitata_HIv3___RNAseq.g4589.t1"; gene_id "MSTRG.14"; gene_
name "Montipora_capitata_HIv3___RNAseq.g4589.t1"; xloc "XLOC_000006"; ref_gene_id "Montipora_capitata_HIv3___RNAseq.g4589.t1"; cmp_ref "Montipora_capitata_HIv3___RNAseq.g4589.t1"; class_code "="; 
tss_id "TSS6";
Montipora_capitata_HIv3___Scaffold_1    StringTie       transcript      204191  204451  .       +       .       transcript_id "Montipora_capitata_HIv3___TS.g26276.t1"; gene_id "MSTRG.15"; gene_nam
e "Montipora_capitata_HIv3___TS.g26276.t1"; xloc "XLOC_000007"; ref_gene_id "Montipora_capitata_HIv3___TS.g26276.t1"; cmp_ref "Montipora_capitata_HIv3___TS.g26276.t1"; class_code "="; tss_id "TSS7
";
Montipora_capitata_HIv3___Scaffold_1    StringTie       transcript      223366  223638  .       +       .       transcript_id "Montipora_capitata_HIv3___TS.g26277.t1"; gene_id "MSTRG.16"; gene_nam
e "Montipora_capitata_HIv3___TS.g26277.t1"; xloc "XLOC_000008"; ref_gene_id "Montipora_capitata_HIv3___TS.g26277.t1"; cmp_ref "Montipora_capitata_HIv3___TS.g26277.t1"; class_code "="; tss_id "TSS8
";
Montipora_capitata_HIv3___Scaffold_1    StringTie       transcript      330242  330991  .       +       .       transcript_id "Montipora_capitata_HIv3___RNAseq.g4592.t1"; gene_id "MSTRG.20"; gene_
name "Montipora_capitata_HIv3___RNAseq.g4592.t1"; xloc "XLOC_000009"; ref_gene_id "Montipora_capitata_HIv3___RNAseq.g4592.t1"; cmp_ref "Montipora_capitata_HIv3___RNAseq.g4592.t1"; class_code "="; 
tss_id "TSS9";
Montipora_capitata_HIv3___Scaffold_1    StringTie       transcript      396628  400449  .       +       .       transcript_id "Montipora_capitata_HIv3___TS.g26284.t1"; gene_id "MSTRG.23"; gene_nam
e "Montipora_capitata_HIv3___TS.g26284.t1"; xloc "XLOC_000010"; ref_gene_id "Montipora_capitata_HIv3___TS.g26284.t1"; cmp_ref "Montipora_capitata_HIv3___TS.g26284.t1"; class_code "="; tss_id "TSS1
0";

wc -l mcap_age_lncRNA_candidates.gtf
54314 mcap_age_lncRNA_candidates.gtf
```

Use bedtools to extract sequence data from the reference genome based on the coordinates from the GTF. The resulting seqs represent potential lncRNA candidates.

```
module load bedtools2/2.31.1
bedtools getfasta -fi /work/pi_hputnam_uri_edu/HI_Genomes/MCapV3/Montipora_capitata_HIv3.assembly.fasta -bed mcap_age_lncRNA_candidates.gtf -fo mcap_age_lncRNA_candidates.fasta -fullHeader -split

zgrep -c ">" mcap_age_lncRNA_candidates.fasta
54314
```

I will need to run [CPC2](https://github.com/gao-lab/CPC2_standalone) to identify transcripts that are coding v. noncoding, but it is not installed on Unity. I asked them to install it and they said: 

```
This is essentially a script that can be set up by following the instructions in the link you provided.  Below is an example to show how to set it up. (Everything run from the terminal)

"""
# Get a compute node
salloc -c 2 --mem=10g

# Create and activate a python virtual environment called venv
python -m venv venv
source venv/bin/activate

# Download prerequisites dor CPC2
pip install biopython six

# Download CPC2
wget https://github.com/gao-lab/CPC2_standalone/archive/refs/tags/v1.0.1.tar.gz
gzip -dc v1.0.1.tar.gz | tar xf -
cd CPC2_standalone-1.0.1

# Set up CPC2
export CPC_HOME="$PWD"
cd libs/libsvm
gzip -dc libsvm-3.18.tar.gz | tar xf -
cd libsvm-3.18
make clean && make

# Run example
cd $CPC_HOME
bin/CPC2.py -i data/example.fa -o example_output
```

Ran the above code and installed here: `/work/pi_hputnam_uri_edu/pgrams/CPC2_standalone-1.0.1`. 

Run CPC2 using my data 

```
cd /work/pi_hputnam_uri_edu/pgrams/CPC2_standalone-1.0.1
source ../venv/bin/activate
export CPC_HOME="$PWD"

python bin/CPC2.py \
  -i /scratch3/workspace/jillashey_uri_edu-coral_age/stringtie/mcap_age_lncRNA_candidates.fasta -o /scratch3/workspace/jillashey_uri_edu-coral_age/stringtie/mcap_age_CPC2.txt
```  
  
Look at results! 

```
head mcap_age_CPC2.txt
#ID     transcript_length       peptide_length  Fickett_score   pI      ORF_integrity   coding_probabi
lity    label
Montipora_capitata_HIv3___Scaffold_1:82396-95409        13013   184     0.34138 9.893076515197755    1
0.817981        coding
Montipora_capitata_HIv3___Scaffold_1:109800-163388      53588   1078    0.2921  7.529540824890138    1
1       coding
Montipora_capitata_HIv3___Scaffold_1:162567-163156      589     140     0.3788  8.867577934265139    1
0.688879        coding
Montipora_capitata_HIv3___Scaffold_1:169949-170561      612     204     0.46457000000000004     8.5183
52699279784     1       0.998454        coding
Montipora_capitata_HIv3___Scaffold_1:170981-172082      1101    367     0.45627 9.624243354797361    1
1       coding
Montipora_capitata_HIv3___Scaffold_1:176399-276376      99977   349     0.29367000000000004     9.5263
15879821777     1       0.995895        coding
Montipora_capitata_HIv3___Scaffold_1:204190-204451      261     87      0.35153 10.518226432800294   1
0.105992        noncoding
Montipora_capitata_HIv3___Scaffold_1:223365-223638      273     91      0.45305 9.851301002502442    1
0.43786 noncoding
Montipora_capitata_HIv3___Scaffold_1:330241-330991      750     53      0.44283 4.5459486007690435   1
0.100661        noncoding

wc -l mcap_age_CPC2.txt
54315 mcap_age_CPC2.txt
```

Filter by noncoding transcripts 

```
awk '$8 != "coding"' mcap_age_CPC2.txt > mcap_age_noncoding_transcripts_info.txt
wc -l mcap_age_noncoding_transcripts_info.txt
31444 mcap_age_noncoding_transcripts_info.txt

awk '$8 == "noncoding" {print $1}' mcap_age_CPC2.txt > mcap_age_noncoding_transcripts_ids.txt
wc -l mcap_age_noncoding_transcripts_ids.txt
31443 mcap_age_noncoding_transcripts_ids.txt
```

Extract putative lncRNAs from fasta

```
sed 's/^transcript:://' mcap_age_noncoding_transcripts_ids.txt > mcap_age_noncoding_transcripts_ids_modified.txt

grep -A 1 -Ff mcap_age_noncoding_transcripts_ids_modified.txt mcap_age_lncRNA_candidates.fasta | sed '/^--$/d' > mcap_age_lncRNA_putative.fasta
```

- extract putative lncrna from gtf 
- filter fasta and gtf by seqs identified by blast 
- gtf to bed 





















To remove potential contimation from rRNAs or tRNAs, download rfam database and blast putative lncRNAs against it. In the scripts folder: `nano rfam_blast.sh`

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=2
#SBATCH --partition=uri-cpu
#SBATCH --no-requeue
#SBATCH --mem=200GB
#SBATCH -t 100:00:00
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error
#SBATCH -D /work/pi_hputnam_uri_edu/jillashey/coral_aging/scripts

# Load modules 
module load uri/main
module load BLAST+/2.15.0-gompi-2023a

echo "Download rfam sequences" $(date)

cd /scratch3/workspace/jillashey_uri_edu-coral_age/rfam_ref

wget -r -np -nd -A .fa.gz ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/fasta_files/

echo "Download complete, cat all sequences" $(date)

cat *fa.gz > rfam_seqs.fa.gz
rm RF*
gunzip rfam_seqs.fa.gz

echo "Cat and unzip complete" $(date)
echo "Blasting putative lncRNAs against rfam db" $(date)
echo "Making blast db" $(date)

makeblastdb -in /scratch3/workspace/jillashey_uri_edu-coral_age/rfam_ref/rfam_seqs.fa -out /scratch3/workspace/jillashey_uri_edu-coral_age/rfam_ref/rfam_db -dbtype nucl

echo "Db creation complete, blast putative lncRNAs against rfam db" $(date)

blastn -query /scratch3/workspace/jillashey_uri_edu-coral_age/stringtie/mcap_age_lncRNA_putative.fasta -db /scratch3/workspace/jillashey_uri_edu-coral_age/rfam_ref/rfam_db -out /scratch3/workspace/jillashey_uri_edu-coral_age/blast/mcap_age_lncRNA_rfam_blastn.tab -evalue 1E-40 -num_threads 10 -max_target_seqs 1 -max_hsps 1 -outfmt 6

echo "Blast complete" $(date)
```

Submitted batch job 38305252. Look at blast output 

```
cd /scratch3/workspace/jillashey_uri_edu-coral_age/blast/

head mcap_age_lncRNA_rfam_blastn.tab 
Montipora_capitata_HIv3___Scaffold_10:50987091-51005681 MU825396.1/212773-212653        96.694  121     4       0       10869   10989   121     1       3.87e-48        202
Montipora_capitata_HIv3___Scaffold_1013:6609-13142      MU826834.1/1912895-1913082      95.531  179     6       2       3345    3522    188     11      1.32e-73        285
Montipora_capitata_HIv3___Scaffold_1013:13927-20177     MU827310.1/718302-718465        96.951  164     5       0       2981    3144    164     1       7.61e-71        276
Montipora_capitata_HIv3___Scaffold_11:3287015-3291500   MU827309.1/2350248-2350386      96.063  127     5       0       4035    4161    127     1       2.04e-50        207
Montipora_capitata_HIv3___Scaffold_13:24569367-24573028 MU827310.1/717380-717262        95.798  119     5       0       1580    1698    1       119     4.65e-46        193
Montipora_capitata_HIv3___Scaffold_13:24617531-24626060 LJWW01001071.1/6624-6820        95.690  116     3       2       2761    2876    142     29      1.80e-43        185
Montipora_capitata_HIv3___Scaffold_13:24738989-24746769 MU826834.1/1912895-1913082      95.506  178     5       2       3075    3251    188     13      2.04e-72        281
Montipora_capitata_HIv3___Scaffold_13:24798335-24808343 LSMT01000541.1/6279-6088        97.917  192     4       0       2639    2830    192     1       7.14e-88        333
Montipora_capitata_HIv3___Scaffold_4:36816332-36834523  MU827817.1/523255-523552        90.333  300     25      4       17249   17546   1       298     7.52e-105       390
Montipora_capitata_HIv3___Scaffold_4:31137676-31146429  MU826834.1/1905971-1906134      93.252  163     11      0       4963    5125    163     1       3.89e-60        241
15 mcap_age_lncRNA_rfam_blastn.tab
```

Only 15 sequences got hits as possible contamination. Make a list of these

```
awk '{print $1}' mcap_age_lncRNA_rfam_blastn.tab > sequences_to_remove.txt
```

Remove these sequences from putative fasta file 

```
cd ../stringtie
grep -v -F -f <(sed 's/^/>/' /scratch3/workspace/jillashey_uri_edu-coral_age/blast/sequences_to_remove.txt) mcap_age_lncRNA_putative.fasta | sed '/^$/d' > mcap_age_lncRNA_putative_purge.fasta

zgrep -c ">" mcap_age_lncRNA_putative.fasta
31443
zgrep -c ">" mcap_age_lncRNA_putative_purge.fasta
31428
```

Quantify putative lncRNA counts with kallisto. Steven quantified lncRNAs from [deep dive](https://github.com/urol-e5/deep-dive-expression/blob/main/E-Peve/code/07-Peve-lncRNA-matrix.qmd) with feature counts but I think either will work. In the scripts folder: `nano kallisto.sh`

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=2
#SBATCH --partition=uri-cpu
#SBATCH --no-requeue
#SBATCH --mem=200GB
#SBATCH -t 100:00:00
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error
#SBATCH -D /work/pi_hputnam_uri_edu/jillashey/coral_aging/scripts

# Load modules 
module load uri/main
module avail Trinity/2.15.1-foss-2022a
module avail kallisto/0.48.0-gompi-2022a

echo "Creating kallisto index from lncRNA fasta" $(date)

kallisto index -i /scratch3/workspace/jillashey_uri_edu-coral_age/kallisto/mcap_age_lncRNA_index /scratch3/workspace/jillashey_uri_edu-coral_age/stringtie/mcap_age_lncRNA_putative_purge.fasta

echo "lncRNA index creation complete, starting read alignment" $(date)






array=($(ls /data/putnamlab/jillashey/DT_Mcap_2023/mRNA/data/trim/*R1_001.fastq.gz)) # call the clean sequences - make an array to align

for i in ${array[@]}; do
    # Extract just the filename from the input FASTQ file path
    filename=$(basename ${i})
    # Construct the output directory path
    output_dir="/data/putnamlab/jillashey/DT_Mcap_2023/mRNA/output/kallisto/kallisto.${filename}"
    # Run kallisto quant
    kallisto quant -i /data/putnamlab/jillashey/DT_Mcap_2023/mRNA/output/kallisto/mcap_lncRNA_index -o ${output_dir} ${i} $(echo ${i} | sed 's/_R1/_R2/')
done

echo "lncRNA alignment complete!" $(date)

```





