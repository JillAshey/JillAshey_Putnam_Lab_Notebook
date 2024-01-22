---
layout: post
title: Pacuta HI 2022
date: '2024-01-18'
categories: Analysis
tags: [Bioinformatics, Pacuta, mRNA]
projects: Pacuta HI 2022
---

## Pacuta 2022 mRNA analysis

These data came from the Pacuta 2022 experiment in Hawaii, done by Federica and myself. In this experiment, larval and spat *Pocillopora acuta* were subjected to a combination of high pH and temperature treatments. The github for that project is [here](https://github.com/fscucchia/Hawaii2022_pH_Temp_Mcap_Pacu). 

Files were downloaded to this location: `/data/putnamlab/KITT/hputnam/20231127_Scucchia_HI`

Make a new directory in my own folder on Andromeda for this experiment

```
cd /data/putnamlab/jillashey
mkdir Pacuta_HI_2022
cd Pacuta_HI_2022
mkdir data scripts output
cd data 
mkdir raw trim 
cd raw
```

Copy the files into the raw data folder 

```
cp /data/putnamlab/KITT/hputnam/20231127_Scucchia_HI/* .
```

Check md has already been done by Hollie. Let's count how many reads each file has. 

```
zgrep -c "@LH00" *.gz

A_R1_001.fastq.gz:27485525
A_R2_001.fastq.gz:27485525
B_R1_001.fastq.gz:30895421
B_R2_001.fastq.gz:30895421
C_R1_001.fastq.gz:31829448
C_R2_001.fastq.gz:31829448
D_R1_001.fastq.gz:29779596
D_R2_001.fastq.gz:29779596
E_R1_001.fastq.gz:32423337
E_R2_001.fastq.gz:32423337
F_R1_001.fastq.gz:28875199
F_R2_001.fastq.gz:28875199
G_R1_001.fastq.gz:31635945
G_R2_001.fastq.gz:31635945
H_R1_001.fastq.gz:29455896
H_R2_001.fastq.gz:29455896
L_R1_001.fastq.gz:29250631
L_R2_001.fastq.gz:29250631
M_R1_001.fastq.gz:28581558
M_R2_001.fastq.gz:28581558
N_R1_001.fastq.gz:28609255
N_R2_001.fastq.gz:28609255
```

QC the raw data 

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
#SBATCH -D /data/putnamlab/jillashey/Pacuta_HI_2022/scripts            
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

module load FastQC/0.11.9-Java-11
module load MultiQC/1.9-intel-2020a-Python-3.8.2

for file in /data/putnamlab/jillashey/Pacuta_HI_2022/data/raw/*fastq.gz
do 
fastqc $file --outdir /data/putnamlab/jillashey/Pacuta_HI_2022/output/fastqc/raw
done

cd /data/putnamlab/jillashey/Pacuta_HI_2022/output/fastqc/raw

multiqc *
```

Submitted batch job 293007

Check to make sure multiqc ran and look at the multiQC output. The phred scores look really good, as does the duplication rate. There is some substantial adapter content so I will trim. Raw QC report can be found [here](https://github.com/fscucchia/Hawaii2022_pH_Temp_Mcap_Pacu/tree/main/output/QC). 

Trim reads w/ fastp. In scripts folder: `nano fastp_qc.sh`

```
#!/bin/bash
#SBATCH -t 100:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Pacuta_HI_2022/scripts            
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

# Load modules needed 
module load fastp/0.19.7-foss-2018b
module load FastQC/0.11.8-Java-1.8
module load MultiQC/1.9-intel-2020a-Python-3.8.2

# Make an array of sequences to trim in raw data directory 

cd /data/putnamlab/jillashey/Pacuta_HI_2022/data/raw/
array1=($(ls *R1_001.fastq.gz))

echo "Read trimming of adapters started." $(date)

# fastp and fastqc loop 
for i in ${array1[@]}; do
    fastp --in1 ${i} \
        --in2 $(echo ${i}|sed s/_R1/_R2/)\
        --out1 /data/putnamlab/jillashey/Pacuta_HI_2022/data/trim/trim.${i} \
        --out2 /data/putnamlab/jillashey/Pacuta_HI_2022/data/trim/trim.$(echo ${i}|sed s/_R1/_R2/) \
        --detect_adapter_for_pe \
        --qualified_quality_phred 30 \
        --unqualified_percent_limit 10 \
        --length_required 100 \
        --cut_right cut_right_window_size 5 cut_right_mean_quality 20
	 fastqc /data/putnamlab/jillashey/Pacuta_HI_2022/output/fastqc/trim/trim.${i}
    fastqc /data/putnamlab/jillashey/Pacuta_HI_2022/output/fastqc/trim/trim.$(echo ${i}|sed s/_R1/_R2/)
done

echo "Read trimming of adapters and QC complete. Starting multiqc" $(date)

cd /data/putnamlab/jillashey/Pacuta_HI_2022/output/fastqc/trim
multiqc *

echo "MultiQC complete" $(date)
```

Submitted batch job 293017. Took a couple of hours, but looks like the fastqc step didn't work...Run fastqc as its own script. In scripts folder: `nano fastqc_trim.sh`

```
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Pacuta_HI_2022/scripts            
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

module load FastQC/0.11.9-Java-11
module load MultiQC/1.9-intel-2020a-Python-3.8.2

for file in /data/putnamlab/jillashey/Pacuta_HI_2022/data/trim/*fastq.gz
do 
fastqc $file --outdir /data/putnamlab/jillashey/Pacuta_HI_2022/output/fastqc/trim
done

cd /data/putnamlab/jillashey/Pacuta_HI_2022/output/fastqc/trim

multiqc *
```

Submitted batch job 293029


Count the number of reads each trimmed file has: 

```
zgrep -c "@LH00" *.gz

trim.A_R1_001.fastq.gz:22328993
trim.A_R2_001.fastq.gz:22328993
trim.B_R1_001.fastq.gz:24904172
trim.B_R2_001.fastq.gz:24904172
trim.C_R1_001.fastq.gz:25854719
trim.C_R2_001.fastq.gz:25854719
trim.D_R1_001.fastq.gz:24265173
trim.D_R2_001.fastq.gz:24265173
trim.E_R1_001.fastq.gz:26184558
trim.E_R2_001.fastq.gz:26184558
trim.F_R1_001.fastq.gz:23339918
trim.F_R2_001.fastq.gz:23339918
trim.G_R1_001.fastq.gz:25448242
trim.G_R2_001.fastq.gz:25448242
trim.H_R1_001.fastq.gz:23634641
trim.H_R2_001.fastq.gz:23634641
trim.L_R1_001.fastq.gz:23596961
trim.L_R2_001.fastq.gz:23596961
trim.M_R1_001.fastq.gz:23003891
trim.M_R2_001.fastq.gz:23003891
trim.N_R1_001.fastq.gz:22825941
trim.N_R2_001.fastq.gz:22825941
```

Once trimming is done, look at the multiqc plots to make sure adapter content is gone and sample quality is high. QC looks good! Trimmed QC report can be found [here](https://github.com/fscucchia/Hawaii2022_pH_Temp_Mcap_Pacu/tree/main/output/QC).

Align reads to Pacuta genome using hisat2. Download version 2 of the pocillopora genome [here](http://cyanophora.rutgers.edu/Pocillopora_acuta/) to Andromeda and unzip the assembly fasta file. In the scripts folder: `nano align.sh`

```
#!/bin/bash
#SBATCH -t 120:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Pacuta_HI_2022/scripts            
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

# load modules needed
module load HISAT2/2.2.1-foss-2019b #Alignment to reference genome: HISAT2
module load SAMtools/1.9-foss-2018b #Preparation of alignment for assembly: SAMtools

echo "Building genome reference" $(date)

# index the reference genome for Pacuta output index to working directory
hisat2-build -f /data/putnamlab/jillashey/genome/Pacuta/V2/Pocillopora_acuta_HIv2.assembly.fasta Pacuta_ref
echo "Referece genome indexed. Starting alingment" $(date)

# This script exports alignments as bam files
# sorts the bam file because Stringtie takes a sorted file for input (--dta)
# removes the sam file because it is no longer needed

array=($(ls /data/putnamlab/jillashey/Pacuta_HI_2022/data/trim/*fastq.gz)) # call the clean sequences - make an array to align

for i in ${array[@]}; do
        sample_name=`echo $i| awk -F [.] '{print $2}'`
        hisat2 -p 8 --rna-strandness RF --dta -x Pacuta_ref -1 ${i} -2 $(echo ${i}|sed s/_R1/_R2/) -S ${sample_name}.sam
        samtools sort -@ 8 -o ${sample_name}.bam ${sample_name}.sam
                echo "${i} bam-ified!"
        rm ${sample_name}.sam
done

echo "Alignment complete!" $(date)
```

Submitted batch job 293038. Took about 7 hours. Bam files are in the scripts folder, move them to hisat2 output folder. 

```
cd /data/putnamlab/jillashey/Pacuta_HI_2022/scripts
mv *bam ../output/hisat2/
cd ../output/hisat2/
```



View mapping percentages 

```
module load SAMtools/1.9-foss-2018b #Preparation of alignment for assembly: SAMtools

for i in *.bam; do
    echo "${i}" >> mapped_reads_counts_Pacuta
    samtools flagstat ${i} | grep "mapped (" >> mapped_reads_counts_Pacuta
done

42773589 + 0 mapped (83.14% : N/A)
D_R2_001.bam
42679238 + 0 mapped (83.07% : N/A)
E_R1_001.bam
45317858 + 0 mapped (81.71% : N/A)
E_R2_001.bam
45247270 + 0 mapped (81.64% : N/A)
F_R1_001.bam
39470870 + 0 mapped (79.99% : N/A)
F_R2_001.bam
39408754 + 0 mapped (79.91% : N/A)
G_R1_001.bam
41604092 + 0 mapped (77.69% : N/A)
G_R2_001.bam
41533656 + 0 mapped (77.62% : N/A)
H_R1_001.bam
38316350 + 0 mapped (77.30% : N/A)
H_R2_001.bam
38212084 + 0 mapped (77.21% : N/A)
L_R1_001.bam
41921009 + 0 mapped (83.53% : N/A)
L_R2_001.bam
41948865 + 0 mapped (83.50% : N/A)
M_R1_001.bam
38624730 + 0 mapped (79.63% : N/A)
M_R2_001.bam
38570811 + 0 mapped (79.57% : N/A)
N_R1_001.bam
27333684 + 0 mapped (57.24% : N/A)
N_R2_001.bam
27314674 + 0 mapped (57.21% : N/A)
```

Pretty high mapping percentages (>75%) for all samples! Sample N had the lowest % at 57%, but that's still relatively high compared to other coral genome alignments that I've done. 

Now assemble and quantify the gene counts using stringtie. In the scripts folder: `nano assemble.sh`

```
#!/bin/bash
#SBATCH -t 100:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=128GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Pacuta_HI_2022/scripts             
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

module load StringTie/2.2.1-GCC-11.2.0

echo "Assembling transcripts using stringtie" $(date)

cd /data/putnamlab/jillashey/Pacuta_HI_2022/output/hisat2

array1=($(ls *.bam)) #Make an array of sequences to assemble

for i in ${array1[@]}; do
    stringtie -p 8 -e -B -G /data/putnamlab/jillashey/genome/Pacuta/V2/Pocillopora_acuta_HIv2.genes.gff3 -A ${i}.gene_abund.tab -o ${i}.gtf ${i}
done

echo "Assembly for each sample complete " $(date)
```

Submitted batch job 293063. Why are the files not combining? Are they supposed to combine in the stringtie step or the hisat2 step? Took about 40 mins and there is no info in the gene abundance tab file, but there is some info in the gtf file. Also the files are still separated by R1 and R2 even though they should be combined at this point. I'm going to rerun the alignment step because I think that is where things went wrong. Going to double check my alignment script. Ah I think it's because I used this array: `array=($(ls /data/putnamlab/jillashey/Pacuta_HI_2022/data/trim/*fastq.gz))` instead of `array=($(ls /data/putnamlab/jillashey/Pacuta_HI_2022/data/trim/*_R1_001.fastq.gz`. Ie I didn't specify R1 for the array. I'm also commenting out the hisat2 build line in the script because the genome has already been indexed. Submitted batch job 293076. Now the output should have the reads combined in one file. Run the `assemble.sh` script. Submitted batch job 293087

While assembly is running, look at mapping percentages. 

```
module load SAMtools/1.9-foss-2018b #Preparation of alignment for assembly: SAMtools

for i in *.bam; do
     echo "${i}" >> mapped_reads_counts_Pacuta
     samtools flagstat ${i} | grep "mapped (" >> mapped_reads_counts_Pacuta
 done
 
A_R1_001.bam
38453816 + 0 mapped (81.09% : N/A)
B_R1_001.bam
43150815 + 0 mapped (82.01% : N/A)
C_R1_001.bam
45052096 + 0 mapped (82.07% : N/A)
D_R1_001.bam
42773589 + 0 mapped (83.14% : N/A)
E_R1_001.bam
45317858 + 0 mapped (81.71% : N/A)
F_R1_001.bam
39470870 + 0 mapped (79.99% : N/A)
G_R1_001.bam
41604092 + 0 mapped (77.69% : N/A)
H_R1_001.bam
38316350 + 0 mapped (77.30% : N/A)
L_R1_001.bam
41921009 + 0 mapped (83.53% : N/A)
M_R1_001.bam
38624730 + 0 mapped (79.63% : N/A)
N_R1_001.bam
27333684 + 0 mapped (57.24% : N/A)
```

Once assembly is done, move the gtf and tab files to the stringtie output folder 

```
cd /data/putnamlab/jillashey/Pacuta_HI_2022/output/hisat2
mv *gtf ../stringtie/
mv *tab ../stringtie/
cd ../stringtie 
```

Make a list of the gtfs

```
ls *.gtf > gtf_list.txt
```

Merge the gtfs into one with stringtie 

```
module purge
module load StringTie/2.1.4-GCC-9.3.0

stringtie --merge -e -p 8 -G /data/putnamlab/jillashey/genome/Pacuta/V2/Pocillopora_acuta_HIv2.genes.gff3 -o Pacuta_merged.gtf gtf_list.txt #Merge GTFs 
```

Assess the accuracy of the merged assembly with gffcompare

```
module purge
module load GffCompare/0.12.1-GCCcore-8.3.0

gffcompare -r /data/putnamlab/jillashey/genome/Pacuta/V2/Pocillopora_acuta_HIv2.genes.gff3 -G -o merged Pacuta_merged.gtf #Compute the accuracy 

  33730 reference transcripts loaded.
  33730 query transfrags loaded.
```

Check out the merged file. 

```
less merged.stats

# gffcompare v0.12.1 | Command line was:
#gffcompare -r /data/putnamlab/jillashey/genome/Pacuta/V2/Pocillopora_acuta_HIv2.genes.gff3 -G -o merged Pacuta_merged.gtf
#

#= Summary for dataset: Pacuta_merged.gtf 
#     Query mRNAs :   33730 in   33611 loci  (25853 multi-exon transcripts)
#            (85 multi-transcript loci, ~1.0 transcripts per locus)
# Reference mRNAs :   33730 in   33611 loci  (25853 multi-exon)
# Super-loci w/ reference transcripts:    33611
#-----------------| Sensitivity | Precision  |
        Base level:   100.0     |   100.0    |
        Exon level:   100.0     |   100.0    |
      Intron level:   100.0     |   100.0    |
Intron chain level:   100.0     |   100.0    |
  Transcript level:   100.0     |   100.0    |
       Locus level:   100.0     |   100.0    |

     Matching intron chains:   25853
       Matching transcripts:   33717
              Matching loci:   33604

          Missed exons:       0/222629  (  0.0%)
           Novel exons:       0/222622  (  0.0%)
        Missed introns:       0/188898  (  0.0%)
         Novel introns:       0/188898  (  0.0%)
           Missed loci:       0/33611   (  0.0%)
            Novel loci:       0/33611   (  0.0%)

 Total union super-loci across all input datasets: 33611 
33730 out of 33730 consensus transcripts written in merged.annotated.gtf (0 discarded as redundant)
```

Make gtf list text file for gene count matrix creation 

```
for filename in *bam.gtf; do echo $filename $PWD/$filename; done > listGTF.txt
```

Download the prepDE.py script from [here](https://github.com/gpertea/stringtie/blob/master/prepDE.py) and put it in the stringtie output folder. Load python and compile the gene count matrix

```
module purge
module load Python/2.7.18-GCCcore-9.3.0

python prepDE.py -g Pacuta_gene_count_matrix.csv -i listGTF.txt
```

The gene count matrix has many STRG gene ids in it, but the transcript count matrix appears to have the gene id names from the gff. Going to download both and check them out. 

