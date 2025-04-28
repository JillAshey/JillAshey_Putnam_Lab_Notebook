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

Submitted batch job 33041964. Trimming looks great! Trimmed QC can be found [here](https://github.com/JillAshey/DevelopmentalTimeseries/blob/main/data/Molecular/egg_sperm/trim_multiqc_report.html). 

Time to align! `nano align.sh`

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
module load uri/main
module load all/HISAT2/2.2.1-gompi-2022a
module load all/SAMtools/1.18-GCC-12.3.0

echo "Building genome reference" $(date)

# index the reference genome for Pacuta output index to working directory
hisat2-build -f /work/pi_hputnam_uri_edu/HI_Genomes/MCapV3/Montipora_capitata_HIv3.assembly.fasta /work/pi_hputnam_uri_edu/HI_Genomes/MCapV3/McapV3_hisat2_ref
echo "Referece genome indexed. Starting alingment" $(date)

# This script exports alignments as bam files
# sorts the bam file because Stringtie takes a sorted file for input (--dta)
# removes the sam file because it is no longer needed

cd /scratch/workspace/jillashey_uri_edu-ashey_scratch/Mcap2023/egg_sperm_trim

echo "Align egg samples (PE)" $(date)
# Paired-end files
for ID in SRR11452251 SRR11452262 SRR11452263; do
    R1="${ID}_1_AdapterTrimmed.fastq"
    R2="${ID}_2_AdapterTrimmed.fastq"
    sample_name="${ID}"
    
    hisat2 -p 8 --rna-strandness RF --dta \
        -x /work/pi_hputnam_uri_edu/HI_Genomes/MCapV3/McapV3_hisat2_ref \
        -1 ${R1} -2 ${R2} -S ${sample_name}.sam

    samtools sort -@ 8 -o ${sample_name}.bam ${sample_name}.sam
    echo "${sample_name} bam-ified!"
    rm ${sample_name}.sam
done
echo "Egg samples (PE) aligned!" $(date)

echo "Align sperm samples (SE)" $(date)
ID="SRR4048723"
READ="${ID}_AdapterTrimmed.fastq"

hisat2 -p 8 --dta \
    -x /work/pi_hputnam_uri_edu/HI_Genomes/MCapV3/McapV3_hisat2_ref \
    -U ${READ} -S ${ID}.sam

samtools sort -@ 8 -o ${ID}.bam ${ID}.sam
echo "${ID} bam-ified!"
rm ${ID}.sam
echo "Sperm samples (SE) aligned!" $(date)

echo "Alignment complete!" $(date)
```

Submitted batch job 33167366. Important note: For the sperm sample: do not specify `--rna-strandness`, as it is unstranded. For the egg samples, set `--rna-strandness RF`. 

Alingment complete! Took 1.5 hours. 

```
19257744 reads; of these:
  19257744 (100.00%) were paired; of these:
    7874384 (40.89%) aligned concordantly 0 times
    7691727 (39.94%) aligned concordantly exactly 1 time
    3691633 (19.17%) aligned concordantly >1 times
    ----
    7874384 pairs aligned concordantly 0 times; of these:
      146546 (1.86%) aligned discordantly 1 time
    ----
    7727838 pairs aligned 0 times concordantly or discordantly; of these:
      15455676 mates make up the pairs; of these:
        14657075 (94.83%) aligned 0 times
        599112 (3.88%) aligned exactly 1 time
        199489 (1.29%) aligned >1 times
61.94% overall alignment rate
[bam_sort_core] merging from 3 files and 8 in-memory blocks...
22408812 reads; of these:
  22408812 (100.00%) were paired; of these:
    6983155 (31.16%) aligned concordantly 0 times
    8980110 (40.07%) aligned concordantly exactly 1 time
    6445547 (28.76%) aligned concordantly >1 times
    ----
    6983155 pairs aligned concordantly 0 times; of these:
      226933 (3.25%) aligned discordantly 1 time
    ----
    6756222 pairs aligned 0 times concordantly or discordantly; of these:
      13512444 mates make up the pairs; of these:
        12316684 (91.15%) aligned 0 times
        846830 (6.27%) aligned exactly 1 time
        348930 (2.58%) aligned >1 times
72.52% overall alignment rate
[bam_sort_core] merging from 5 files and 8 in-memory blocks...
21204959 reads; of these:
  21204959 (100.00%) were paired; of these:
    5630028 (26.55%) aligned concordantly 0 times
    10575974 (49.88%) aligned concordantly exactly 1 time
    4998957 (23.57%) aligned concordantly >1 times
    ----
    5630028 pairs aligned concordantly 0 times; of these:
      188403 (3.35%) aligned discordantly 1 time
    ----
    5441625 pairs aligned 0 times concordantly or discordantly; of these:
      10883250 mates make up the pairs; of these:
        9800975 (90.06%) aligned 0 times
        819802 (7.53%) aligned exactly 1 time
        262473 (2.41%) aligned >1 times
76.89% overall alignment rate
[bam_sort_core] merging from 4 files and 8 in-memory blocks...
17901354 reads; of these:
  17901354 (100.00%) were unpaired; of these:
    3546701 (19.81%) aligned 0 times
    12504646 (69.85%) aligned exactly 1 time
    1850007 (10.33%) aligned >1 times
80.19% overall alignment rate
```

Alignment rates look great. Time to assemble! `nano assemble.sh`

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
module load uri/main
module load all/StringTie/2.2.1-GCC-11.2.0

echo "Assembling transcripts using stringtie" $(date)

cd /scratch/workspace/jillashey_uri_edu-ashey_scratch/Mcap2023/egg_sperm_trim

array1=($(ls *.bam)) #Make an array of sequences to assemble

for i in ${array1[@]}; do
    stringtie -p 8 -e -B -G /work/pi_hputnam_uri_edu/HI_Genomes/MCapV3/Montipora_capitata_HIv3.genes_fixed.gff3 -A ${i}.gene_abund.tab -o ${i}.gtf ${i}
done

echo "Assembly complete " $(date)
```

Submitted batch job 33316146

Do the rest of the assembly in interactive mode. 

```
cd /scratch/workspace/jillashey_uri_edu-ashey_scratch/Mcap2023/egg_sperm_trim
salloc 
module load uri/main
module load all/StringTie/2.2.1-GCC-11.2.0

# List of gtfs
ls *.gtf > gtf_list.txt

# Merge gtfs into single gtf 
stringtie --merge -e -p 8 -G /work/pi_hputnam_uri_edu/HI_Genomes/MCapV3/Montipora_capitata_HIv3.genes_fixed.gff3 -o Mcap_egg_sperm_merged.gtf gtf_list.txt #Merge GTFs 

module load all/GffCompare/0.12.6-GCC-11.2.0
gffcompare -r /work/pi_hputnam_uri_edu/HI_Genomes/MCapV3/Montipora_capitata_HIv3.genes_fixed.gff3 -G -o merged Mcap_egg_sperm_merged.gtf #Compute the accuracy 
54384 reference transcripts loaded.
  2 duplicate reference transcripts discarded.
  54384 query transfrags loaded.
```

Look at merge stats 

```
# gffcompare v0.12.6 | Command line was:
#gffcompare -r /work/pi_hputnam_uri_edu/HI_Genomes/MCapV3/Montipora_capitata_HIv3.genes.gff3 -G -o merged Mcap_egg_sperm_merged.gtf
#
#= Summary for dataset: Mcap_egg_sperm_merged.gtf 
#     Query mRNAs :   54384 in   54185 loci  (36023 multi-exon transcripts)
#            (141 multi-transcript loci, ~1.0 transcripts per locus)
# Reference mRNAs :   54382 in   54185 loci  (36023 multi-exon)
# Super-loci w/ reference transcripts:    54185
#-----------------| Sensitivity | Precision  |
        Base level:   100.0     |   100.0    |
        Exon level:   100.0     |   100.0    |
      Intron level:   100.0     |   100.0    |
Intron chain level:   100.0     |   100.0    |
  Transcript level:   100.0     |   100.0    |
       Locus level:   100.0     |   100.0    |
     Matching intron chains:   36023
       Matching transcripts:   54375
              Matching loci:   54185
          Missed exons:       0/256029  (  0.0%)
           Novel exons:       0/256029  (  0.0%)
        Missed introns:       0/201643  (  0.0%)
         Novel introns:       0/201643  (  0.0%)
           Missed loci:       0/54185   (  0.0%)
            Novel loci:       0/54185   (  0.0%)
 Total union super-loci across all input datasets: 54185 
54384 out of 54384 consensus transcripts written in merged.annotated.gtf (0 discarded as redundant)
```

Make gtf list text file for count matrix generation

```
for filename in *bam.gtf; do echo $filename $PWD/$filename; done > listGTF.txt
```

Download the prepDE.py script from [here](https://github.com/gpertea/stringtie/blob/master/prepDE.py) and put it in the scripts folder. Load python and compile the gene count matrix

```
salloc 
module load uri/main
module load Python/2.7.18-GCCcore-9.3.0

python /work/pi_hputnam_uri_edu/jillashey/Mcap_2023/egg_sperm/scripts/prepDE.py -g Mcap_egg_sperm_gene_count_matrix.csv -i listGTF.txt
```

Keeps saying `Illegal instruction`. Try running as a job. `nano prepDE.sh`

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

module load uri/main
module load Python/2.7.18-GCCcore-9.3.0

cd /scratch/workspace/jillashey_uri_edu-ashey_scratch/Mcap2023/egg_sperm_trim

echo "Creating count matrices " $(date)

python /work/pi_hputnam_uri_edu/jillashey/Mcap_2023/egg_sperm/scripts/prepDE.py -g Mcap_egg_sperm_gene_count_matrix.csv -i listGTF.txt

echo "Count matrices complete" $(date)
```

Submitted batch job 33316331. Success! Move csvs to work directory

```
mv /scratch/workspace/jillashey_uri_edu-ashey_scratch/Mcap2023/egg_sperm_trim/*csv /work/pi_hputnam_uri_edu/jillashey/Mcap_2023/egg_sperm/output
```

Download csvs to computer. 


