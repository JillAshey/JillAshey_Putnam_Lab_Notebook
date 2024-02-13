---
layout: post
title: Apulchra genome assembly 
date: '2024-02-06'
categories: Analysis, Genome Assembly
tags: [Bioinformatics, Genome, Assembly]
projects: Apulchra genome
---

## Apulchra genome assembly 

Sperm and tissue from adult *Acropora pulchra* colonies were collected from Moorea, French Polynesia and sequencing with PacBio (long reads) and Illumina (short reads). This post will detail the genome assembly notes. The github for this project is [here](https://github.com/hputnam/Apulchra_genome/tree/main). 

I'm going to write notes and code chronologically so that I can keep track of what I'm doing each day. When assembly is complete, I will compile the workflow in a separate post. 

### 20240206

Met w/ Ross and Hollie today re Apulchra genome assembly. We decided to move forward with the workflow from [Stephens et al. 2022](https://academic.oup.com/gigascience/article/doi/10.1093/gigascience/giac098/6815755) which assembled genomes for 4 Hawaiian corals. 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/Stephens_et_al_2022_workflow.png)

PacBio long reads were received in late Jan/early Feb 2024. According to reps from Genohub, the [PacBio raw output](https://github.com/hputnam/Apulchra_genome/blob/main/DNA_Seq_Info/20240129_Project_6693786_Acropora_pulchra_PacBio_Seq_Summary.pdf) looks good. 

We decided to move forward with [Canu](https://canu.readthedocs.io/en/latest/index.html) to assembly the genome. Canu is specialized to assemble PacBio sequences, operating in three phases: correction, trimming and assembly. According to the Canu website, "The correction phase will improve the accuracy of bases in reads. The trimming phase will trim reads to the portion that appears to be high-quality sequence, removing suspicious regions such as remaining SMRTbell adapter. The assembly phase will order the reads into contigs, generate consensus sequences and create graphs of alternate paths." 

The PacBio files that will be used for assembly are located here on Andromeda: `/data/putnamlab/KITT/hputnam/20240129_Apulchra_Genome_LongRead`. The files in the folder that we will use are: 

```
m84100_240128_024355_s2.hifi_reads.bc1029.bam
m84100_240128_024355_s2.hifi_reads.bc1029.bam.pbi
```

The bam file contains all of the read information in computer language and the pbi file is an index file of the bam. Both are needed in the same folder to run Canu. 

For Canu, input files must be fasta or fastq format. I'm going to use `bam2fastq`from the [PacBio github](https://github.com/pacificbiosciences/pbtk/). This module is not on Andromeda so I will need to install it via conda. 

The PacBio sequencing for the Apul genome were done with HiFi sequencing that are produced with circular consensus sequencing on PacoBio long read systems. Here's how HiFi reads are generated from the [PacBio website](https://www.pacb.com/technology/hifi-sequencing/):

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/hifi_read_generation.png)

Since Hifi sequencing was used, a specific HiCanu flag (`-pacbio-hifi`) must be used. Additionally, in the Canu [tutorial](https://canu.readthedocs.io/en/latest/quick-start.html#), it says that if this flag is used, the program is assuming that the reads are trimmed and corrected already. However, our reads are not. I'm going to try to run the first pass at Canu with the `-raw` flag. 

Before Canu, I will run `bam2fastq`. This is a part of the PacBio BAM toolkit package `pbtk`. I need to create the conda environment and install the package. Load miniconda module: `module load Miniconda3/4.9.2`. Create a conda environment. 

```
conda create --prefix /data/putnamlab/conda/pbtk
```

Install the package. Once the package is installed on the server, anyone in the Putnam lab group can use it. 

```
conda activate /data/putnamlab/conda/pbtk
conda install -c bioconda pbtk
```

In my own directory, make a new directory for genome assembly

```
cd /data/putnamlab/jillashey
mkdir Apul_Genome
cd Apul_Genome
mkdir assembly structural functional
cd assembly
mkdir scripts data output
```

Run code to make the PacBio bam file to a fastq file. In the scripts folder: `nano bam2fastq.sh`

```
#!/bin/bash -i
#SBATCH -t 500:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=500GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Apul_Genome/assembly/scripts
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

conda activate /data/putnamlab/conda/pbtk

echo "Convert PacBio bam file to fastq file" $(date)

bam2fastq -o /data/putnamlab/jillashey/Apul_Genome/assembly/data/m84100_240128_024355_s2.hifi_reads.bc1029.fastq /data/putnamlab/KITT/hputnam/20240129_Apulchra_Genome_LongRead/m84100_240128_024355_s2.hifi_reads.bc1029.bam

echo "Bam to fastq complete!" $(date)

conda deactivate
```

Submitted batch job 294235

### 20240208

Job pended for about a day, then ran in 1.5 hours. I got this error message: `bash: cannot set terminal process group (-1): Function not implemented bash: no job control in this shell`, but not sure why. A fastq.gz file was produced! The file is pretty large (35G).

```
less m84100_240128_024355_s2.hifi_reads.bc1029.fastq.fastq.gz 

@m84100_240128_024355_s2/261887593/ccs
TATAAGTTTTACAGCTGCCTTTTGCTCAGCAAAGAAAGCAGCATTGTTATTGAACAGAAAAAGCCTTTTGGTGATATAAAGGTTTCTAAGGGACCAAAGTTTGATCTAGTATGCTAAGTGTGGTGGGTTAAAACTTTGTTTCACCTTTTTCCCGTGATGTTACAAAATTGGTGCAAATATTCATGACGTCGTACTCAAGTCTGACACTAAGAGCATGCAATTACTTAAACAAACAAGCCATACACCAATAAACTGAAGCTCTGTCAACTAGAAAACCTTTGAGTATTTTCATTGTAAAGTACAAGTGGTATATTGTCACTTGCTTTTACAACTTGAAGAACTCACAGTTAAGTTACTAGATTCACCATAGTGCTTGGCAATGAAGAAGCCAAATCACATAAAGTCGGAGCATGTGGTGTTTAGACCTAATCAAACAAGAACACAATATTTAGTACCTGCATCCTGTCTAAGGAGGAAATTTTAAGCTGCTTTCTTTTAAATTTTTTTTATTAGCATTTCAATGGTTGAGGTCGATTATAGGTGCTAGGCTTTAATTCCGTACTATGAAAGAAGAAAGGTCGTTGTTATTAACCATGTCAAACAGAGAAACACATGGTAAAAAATTGACTTCCTTTTCCTCTCGTTGCCACTTAAGCTTAATGATGGTGTTTGACCTGAAAGATGTTACAATTGTTTTAGATGAAAAGACTGTTCTGCGTAAAACAGTGAAGCCTCCCAACTTATTTTGTTATGTGGATTTTGTTGTCTTGTTAGTAACATGTATTGGACTATCTTTTGTGAGTACATAGCTTTTTTTCCATCAACTGACTATATACGTGGTGTAATTTGAGATCATGCCTCCAAGTGTTAGTCTTTTGTTTGGGGCTAACTCGTAAAAGACAAAGGGAGGGGGGTTGTCTAATTCCTAAGCAAAGCATTAAGTTTAACACAGGAAATTGTTTGCGTTGATATTGCTATCCTTTCAGCCCCAAACAAAAAATTTAATGGTTATTTTATTTTACATCTATTGTAAATATATTTTAACATTAATTTTTATTATTGCACTGTAAATACTTGTACTAATGTTCTGTTTGAATTAATTTTGATTCATTCCTTGTGCTTACAACAACAGGGATACAAAACCGATATGTATAATAATACTATTAGAGATGCTTATTTGCATTTTTAGCCCATACCATGAGTTTTAATAACGCCAGGCCATTGGAGATTTTATGGAGTGAGGATTCATTGTACAAACATGGTTGATTTAATATTAAAGTTGTATCCAAATAATTAATATCTGCTGTGATCAGTGAAAGATTGACCTTTCAGTTGTTTGGTTGCACCTTCATCTTATTGGAAACAACTGAATGGAGCATCTTTCCAGTTTAAAAATGTACCACTGCCCACTTTCATGAAGTTATGCCACATATTAATAATGACTATTAATTGTTGAAAACCCTTCTTCCAAAATGTTTCCATTTATTTGTAATAGCATATGTGGTCCATCAGAACAATAATTTAAATCATTACTATTAATAATTTTCCAATAACTGACTTTCAAACCTAGCCAACAAGCATAAGTCAGTAAGCCACAGAGTCCAGAGATACACTTACACTTTACTTTTCACTTCTGAAACATTTTATAATCTCAGTATGAGCATAGAACTTTTCAGTTGGGCAGCATGGAATAGAACCTTTGGACCCCTCTGTGAATATCAAAAATAGGCAACCACTTCCAGCATACACTCTAGCCTCCTTCATAAAGCAAGCCTTAGTGTTTTAGCTTCTACTAGTTAGATTCATTTTAAAAGAAGTTCAGTATACTTAATCTTATAGAAGCTGATTGTGATATAATTGCATAGGTGGATCTCAGAAAAGTGAGATGTAGCTGTCAAATTAAAAGAAGTCCTTTCCAAGCGTAGCTTCTGATAAACAATGCATTTTAGTTAACATTGGATTATGGTTTCAAAGGACTTGTAAAGCTAAATTCAAGTTTTTATGACAACTTGAAAGCCTTTGCCACAGTCTCCGCTGATTTAAGACTTCCATCAAAGTTAGAGTGGTGTGAATGCATCTCCACATGCAATTAATAAAGGTGAGGCAGACAACACAAAACACCCTGGTGCACCATCAACTCCACGGATCACTTGACTGTAACGCCATCTTATACAGCGACTGCCAATTGGAACTGGAAGATCAGGAATGATCTCTTTCACATGGGAAATGAGCATGGTCTTGATAATGCTTTCATCAGTATCCCAATGTTGAAGACAGAAGGACTTTGTTGAATGTACTAAAAGAGAGGGACCAACATCCACTGGATCTGGAACAGAATGCCAAAGCATATGCAAAATGTTATCATGTATAATAAAGAAAAAACCAAACTCATGAAGTTCAACGGTCACAGGCTTTGGTAAAAGACTGTACATGGTAAAGTACTCAGGGTGAAAGTTTGTATCCAAAAAACGAAAAGCTAACCTGTACCACATCTCTTCCGAGAGTCCACTGACACAAATGCAATATTAGGATTGCCAGTCGCATATTTACACGTCCATGGCACATTGATAACTGCTTCATGATCAAAAAAGTATCCTACTGCAAACCGTGATGAATATTCCACATTTTGTAAGGATTTGATTTCATTTTGTAAGAATGCTTGAATTGAACCTTGTAGCTGAAGAAGTTGTGGTACTGGAATAGTTACTATGACTGA
```

See how many `@m84100` are in the file. I'm not sure what these stand for, maybe contigs? 

```
zgrep -c "@m84100" m84100_240128_024355_s2.hifi_reads.bc1029.fastq.fastq.gz 
5898386
```

More than 5 million, so many not contigs? I guess it represents the number of HiFi reads generated. Now time to run Canu! Canu is already installed on the server, which is nice. 

In the scripts folder: `nano canu.sh`

```
#!/bin/bash 
#SBATCH -t 500:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=500GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Apul_Genome/assembly/scripts
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

module load canu/2.2-GCCcore-11.2.0 

cd /data/putnamlab/jillashey/Apul_Genome/assembly/data

echo "Unzip paco-bio fastq file" $(date)

gunzip m84100_240128_024355_s2.hifi_reads.bc1029.fastq.fastq.gz

echo "Unzip complete, starting assembly" $(date)

canu -p apul -d /data/putnamlab/jillashey/Apul_Genome/assembly/data genomeSize=475m -raw -pacbio-hifi m84100_240128_024355_s2.hifi_reads.bc1029.fastq.fastq

echo "Canu assembly complete" $(date)
```

I'm not sure if the `raw` and `-pacbio-hifi` will be compatible, as the Canu tutorial says that the `-pacbio-hifi` assumes that the input is trimmed and corrected (still not sure what this means). Submitted batch job 294325

### 20240212

The canu script appears to have ran. The script `canu.sh` itself ran in ~40 mins, but it spawned 10000s of other jobs on the server (for parallel processing, I'm guessing). Since I didn't start those jobs, I didn't get emails when they finished, so I just had to check the server every few hours. It took about a day for the rest of the jobs to finish. First, I'm looking at the `slurm-294325.error` output file. There's a lot in this file, but I will break it down. 

It first provides details on the slurm support and associated memory, as well as the number of threads that each portion of canu will need to run.

```
-- Slurm support detected.  Resources available:
--     25 hosts with  36 cores and  124 GB memory.
--      3 hosts with  24 cores and  124 GB memory.
--      1 host  with  36 cores and  123 GB memory.
--      1 host  with  48 cores and  753 GB memory.
--      8 hosts with  36 cores and   61 GB memory.
--      1 host  with  48 cores and  250 GB memory.
--      2 hosts with  36 cores and  250 GB memory.
--      2 hosts with  48 cores and  375 GB memory.
--      4 hosts with  36 cores and  502 GB memory.
--
--                         (tag)Threads
--                (tag)Memory         |
--        (tag)             |         |  algorithm
--        -------  ----------  --------  -----------------------------
-- Grid:  meryl     12.000 GB    6 CPUs  (k-mer counting)
-- Grid:  hap       12.000 GB   12 CPUs  (read-to-haplotype assignment)
-- Grid:  cormhap   13.000 GB   12 CPUs  (overlap detection with mhap)
-- Grid:  obtovl     8.000 GB    6 CPUs  (overlap detection)
-- Grid:  utgovl     8.000 GB    6 CPUs  (overlap detection)
-- Grid:  cor        -.--- GB    4 CPUs  (read correction)
-- Grid:  ovb        4.000 GB    1 CPU   (overlap store bucketizer)
-- Grid:  ovs       16.000 GB    1 CPU   (overlap store sorting)
-- Grid:  red       10.000 GB    6 CPUs  (read error detection)
-- Grid:  oea        8.000 GB    1 CPU   (overlap error adjustment)
-- Grid:  bat       64.000 GB    8 CPUs  (contig construction with bogart)
-- Grid:  cns        -.--- GB    8 CPUs  (consensus)
--
-- Found trimmed raw PacBio HiFi reads in the input files.
```

The file also says that it skipped the correction and trimming steps. This indicates that adding the `-raw` flag didn't work. 

```
--   Stages to run:
--     assemble HiFi reads.
--
--
-- Correction skipped; not enabled.
--
-- Trimming skipped; not enabled.
```

Two histograms are then presented. The first is a histogram of correct reads: 

```
-- In sequence store './apul.seqStore':
--   Found 5897694 reads.
--   Found 79182880880 bases (166.7 times coverage).
--    Histogram of corrected reads:
--    
--    G=79182880880                      sum of  ||               length     num
--    NG         length     index       lengths  ||                range    seqs
--    ----- ------------ --------- ------------  ||  ------------------- -------
--    00010        26133    264284   7918288640  ||       1150-2287         8297|--
--    00020        22541    592551  15836596549  ||       2288-3425        76768|-------------
--    00030        20184    964605  23754880672  ||       3426-4563       316378|---------------------------------------------------
--    00040        18285   1377072  31673163019  ||       4564-5701       397324|---------------------------------------------------------------
--    00050        16551   1832151  39591448183  ||       5702-6839       374507|------------------------------------------------------------
--    00060        14786   2337786  47509730799  ||       6840-7977       351843|--------------------------------------------------------
--    00070        12854   2910882  55428020435  ||       7978-9115       339864|------------------------------------------------------
--    00080        10612   3585762  63346313974  ||       9116-10253      340204|------------------------------------------------------
--    00090         7724   4449769  71264594524  ||      10254-11391      341160|-------------------------------------------------------
--    00100         1150   5897693  79182880880  ||      11392-12529      343170|-------------------------------------------------------
--    001.000x             5897694  79182880880  ||      12530-13667      340043|------------------------------------------------------
--                                               ||      13668-14805      336214|------------------------------------------------------
--                                               ||      14806-15943      328927|-----------------------------------------------------
--                                               ||      15944-17081      316290|---------------------------------------------------
--                                               ||      17082-18219      293393|-----------------------------------------------
--                                               ||      18220-19357      261212|------------------------------------------
--                                               ||      19358-20495      225351|------------------------------------
--                                               ||      20496-21633      188701|------------------------------
--                                               ||      21634-22771      154123|-------------------------
--                                               ||      22772-23909      124842|--------------------
--                                               ||      23910-25047       99659|----------------
--                                               ||      25048-26185       78280|-------------
--                                               ||      26186-27323       61530|----------
--                                               ||      27324-28461       47760|--------
--                                               ||      28462-29599       37668|------
--                                               ||      29600-30737       29339|-----
--                                               ||      30738-31875       22548|----
--                                               ||      31876-33013       17292|---
--                                               ||      33014-34151       13055|---
--                                               ||      34152-35289        9797|--
--                                               ||      35290-36427        7251|--
--                                               ||      36428-37565        5005|-
--                                               ||      37566-38703        3560|-
--                                               ||      38704-39841        2474|-
--                                               ||      39842-40979        1553|-
--                                               ||      40980-42117        1006|-
--                                               ||      42118-43255         604|-
--                                               ||      43256-44393         334|-
--                                               ||      44394-45531         184|-
--                                               ||      45532-46669         102|-
--                                               ||      46670-47807          45|-
--                                               ||      47808-48945          18|-
--                                               ||      48946-50083          11|-
--                                               ||      50084-51221           4|-
--                                               ||      51222-52359           1|-
--                                               ||      52360-53497           2|-
--                                               ||      53498-54635           0|
--                                               ||      54636-55773           0|
--                                               ||      55774-56911           0|
--                                               ||      56912-58049           1|-
```

There's also a histogram of corrected-trimmed reads, but it is the exact same as the histogram above. The histogram represents the length ranges for each sequence and the number of sequences that have that length range. For example, if we look at the top row, there are 8297 sequences that range in length from 1150-2287 bp. 

It looks like canu did run jobs by itself: 

```
--  For 5897694 reads with 79182880880 bases, limit to 791 batches.
--  Will count kmers using 16 jobs, each using 13 GB and 6 threads.
--
-- Finished stage 'merylConfigure', reset canuIteration.
--
-- Running jobs.  First attempt out of 2.
--
-- 'meryl-count.jobSubmit-01.sh' -> job 294326 tasks 1-16.
--
----------------------------------------
-- Starting command on Thu Feb  8 14:32:12 2024 with 8923.722 GB free disk space

    cd /glfs/brick01/gv0/putnamlab/jillashey/Apul_Genome/assembly/data
    sbatch \
      --depend=afterany:294326 \
      --cpus-per-task=1 \
      --mem-per-cpu=4g   \
      -D `pwd` \
      -J 'canu_apul' \
      -o canu-scripts/canu.01.out  canu-scripts/canu.01.sh

-- Finished on Thu Feb  8 14:32:13 2024 (one second) with 8923.722 GB free disk space
```

Let's look at the output files that canu produced in `/data/putnamlab/jillashey/Apul_Genome/assembly/data`. I'll be using the [Canu tutorial output](https://canu.readthedocs.io/en/latest/tutorial.html#outputs) information to understand outputs. 

```
-rwxr-xr-x. 1 jillashey 1.1K Feb  8 14:13 apul.seqStore.sh
-rw-r--r--. 1 jillashey  951 Feb  8 14:31 apul.seqStore.err
drwxr-xr-x. 3 jillashey 4.0K Feb  8 14:32 apul.seqStore
drwxr-xr-x. 9 jillashey 4.0K Feb  8 23:57 unitigging
drwxr-xr-x. 2 jillashey 4.0K Feb  9 01:36 canu-scripts
lrwxrwxrwx. 1 jillashey   24 Feb  9 01:36 canu.out -> canu-scripts/canu.09.out
drwxr-xr-x. 2 jillashey 4.0K Feb  9 01:36 canu-logs
-rw-r--r--. 1 jillashey  23K Feb  9 01:47 apul.report
-rw-r--r--. 1 jillashey 7.0M Feb  9 01:51 apul.contigs.layout.tigInfo
-rw-r--r--. 1 jillashey 155M Feb  9 01:51 apul.contigs.layout.readToTig
-rw-r--r--. 1 jillashey 2.8G Feb  9 01:57 apul.unassembled.fasta
-rw-r--r--. 1 jillashey 943M Feb  9 02:01 apul.contigs.fasta
```

The `apul.report` file will provide information about the analysis during assembly, including histogram of read lengths, the histogram or k-mers in the raw and corrected reads, the summary of corrected data, summary of overlaps, and the summary of contig lengths. The histogram of read lengths is the same as in the error file above. There is also a histogram (?) of the mer information: 

```
--  22-mers                                                                                           Fraction
--    Occurrences   NumMers                                                                         Unique Total
--       1-     1         0                                                                        0.0000 0.0000
--       2-     2   4316301 ****                                                                   0.0150 0.0002
--       3-     4    855268                                                                        0.0172 0.0002
--       5-     7    183637                                                                        0.0183 0.0002
--       8-    11     62578                                                                        0.0187 0.0002
--      12-    16     29616                                                                        0.0188 0.0002
--      17-    22     22278                                                                        0.0189 0.0003
--      23-    29     20222                                                                        0.0190 0.0003
--      30-    37     28236                                                                        0.0190 0.0003
--      38-    46     73803                                                                        0.0192 0.0003
--      47-    56    556391                                                                        0.0194 0.0003
--      57-    67   6166302 ******                                                                 0.0218 0.0010
--      68-    79  35800144 ***************************************                                0.0477 0.0098
--      80-    92  63190280 ********************************************************************** 0.1835 0.0633
--      93-   106  26607016 *****************************                                          0.3988 0.1607
--     107-   121   2495822 **                                                                     0.4799 0.2024
--     122-   137   2376701 **                                                                     0.4871 0.2066
--     138-   154  16117515 *****************                                                      0.4965 0.2131
--     155-   172  47510630 ****************************************************                   0.5575 0.2605
--     173-   191  41907409 **********************************************                         0.7262 0.4060
--     192-   211   9190117 **********                                                             0.8650 0.5375
--     212-   232   1806623 **                                                                     0.8934 0.5669
--     233-   254   4028020 ****                                                                   0.8997 0.5743
--     255-   277   3875382 ****                                                                   0.9140 0.5926
--     278-   301   1456651 *                                                                      0.9271 0.6107
--     302-   326   1905710 **                                                                     0.9319 0.6181
--     327-   352   3003944 ***                                                                    0.9388 0.6293
--     353-   379   1723658 *                                                                      0.9491 0.6477
--     380-   407    968475 *                                                                      0.9549 0.6587
--     408-   436   1249856 *                                                                      0.9583 0.6657
--     437-   466    890757                                                                        0.9626 0.6752
--     467-   497    768714                                                                        0.9656 0.6824
--     498-   529    937920 *                                                                      0.9683 0.6892
--     530-   562    618907                                                                        0.9716 0.6978
--     563-   596    582755                                                                        0.9737 0.7039
--     597-   631    535408                                                                        0.9757 0.7100
--     632-   667    441122                                                                        0.9775 0.7159
--     668-   704    459953                                                                        0.9791 0.7211
--     705-   742    367266                                                                        0.9807 0.7268
--     743-   781    354455                                                                        0.9819 0.7316
--     782-   821    304354                                                                        0.9832 0.7365
```

There are 22-mers. A k-mer are substrings of length k contained in a biological sequence. For example, the term k-mer refers to all of a sequence's subsequences of length k such that the sequence AGAT would have four monomers (A, G, A, and T), three 2-mers (AG, GA, AT), two 3-mers (AGA and GAT) and one 4-mer (AGAT). So if we have 22-mers, we have subsequences of 22 nt? The Canu documentation says that k-mer histograms with more than 1 peak likely indicate a heterozygous genome. I'm not sure if the stars represent peaks or counts but if this is a histogram of k-mer information, it has two peaks, indicating a heterozygous genome.

The Canu documentation states that corrected read reports should be given with information about number of reads, coverage, N50, etc. My log file does not have this information, likely because the trimming and correcting steps were not performed. Instead, I have this information: 

```
--   category            reads     %          read length        feature size or coverage  analysis
--   ----------------  -------  -------  ----------------------  ------------------------  --------------------
--   middle-missing       4114    0.07    10652.45 +- 6328.73       1185.95 +- 2112.86    (bad trimming)
--   middle-hump          4148    0.07    11485.92 +- 4101.92       4734.88 +- 3702.70    (bad trimming)
--   no-5-prime           8533    0.14     9311.86 +- 5288.80       2087.03 +- 3315.79    (bad trimming)
--   no-3-prime          10888    0.18     7663.24 +- 5159.57       1677.51 +- 3090.22    (bad trimming)
--   
--   low-coverage        48831    0.83     6823.01 +- 3725.15         16.35 +- 15.52      (easy to assemble, potential for lower quality consensus)
--   unique            5419159   91.89     9403.44 +- 4761.57        110.87 +- 40.33      (easy to assemble, perfect, yay)
--   repeat-cont         93801    1.59     7730.76 +- 4345.55       1001.60 +- 665.30     (potential for consensus errors, no impact on assembly)
--   repeat-dove           380    0.01    23076.97 +- 3235.06        878.08 +- 510.47     (hard to assemble, likely won't assemble correctly or eve
n at all)
--   
--   span-repeat         64724    1.10    11397.59 +- 5058.57       2335.58 +- 2836.91    (read spans a large repeat, usually easy to assemble)
--   uniq-repeat-cont   182925    3.10     9764.80 +- 4218.80                             (should be uniquely placed, low potential for consensus e
rrors, no impact on assembly)
--   uniq-repeat-dove    14288    0.24    17978.89 +- 4847.85                             (will end contigs, potential to misassemble)
--   uniq-anchor         19659    0.33    11312.36 +- 4510.86       5442.58 +- 3930.94    (repeat read, with unique section, probable bad read)
```

I'm not sure why I'm getting all of this information or what it means. There is a high % of unique reads in the data which is good. In the file, there is also information about edges (not sure what this means), as well as error rates. May discuss further with Hollie. 

The Canu output documentation says that I'm supposed to get a file with corrected and trimmed reads but I don't have those. I do have `apul.unassembled.fasta` and `apul.contigs.fasta`. 

```
head apul.unassembled.fasta
>tig00000838 len=19665 reads=4 class=unassm suggestRepeat=no suggestBubble=no suggestCircular=no trim=0-19665
TAAAAACATTGATTCTTGTTTCAATATGAGACTTGTTTCGGAAGATGTTCGCGCAGGTTACATTTCATAATCCACAAGAAATGCGACATCGCCAACCTTA
CTTTAGTGTTTGCATTTAAGCAAAACAATGATAAAGAAACAAATCTCACATCTCGCAAAAGTATGCATTCTATGAAGAACAATGAAAATTAATGAAAGTG
AATCTTACACCTCCTATTCAAGACGCCGCATTAATTCAACTTGTTGATTTCTCCTAAAACGCTTTCTTTTAGAAGGCTTTCTTAGTCTTTCAATTGTAAA
GAATACATAAAGGACTCATGACCACTTATGTTCTTAAGTGTTACTGCTGCTTTAAAACACGTTACAAACCACATGTGAATATAGTTGCGGCACAGAAGGG
AAAATCGCTGAAATGCTGTCCAAATATACACAATATCATTAAGTAAAGTACGATCGTCCGGGTGAGTGTAGTCCTGAGAAGGACTGTTTGAGATGACATT
GACTGACGTTTCGACAACCTGAGCGGAAGTCATCTTCAGAGTCATCTTCACTTGACTCTGAAGATGACTTCCGCTCAGGTTGTCGAAACGTCAGTCAATG
TCATCTCAAACAGTCCTTCTCAGGACTACACTCACCCGGACGATCGTACTTTACTTAATGATATGACTCCTGGGTTCAAACCATTTACAAATATACACAA
GTTTGAAAGATCATATCGCCTGCCAGTTTTACAACTCGTCTTAGACACAATGGAATACAAAACCCTACCGAACGAATACCTTTGATTGAGATTTATGAAT
GTGAAACAGCGACTTCGAGAGAAAAACGAATTCTTAAAAATGCAGTTCAACTCTATCGTCATTCAACTGAAGCCAGCCGTGGCTCGGCTTCACAAGAAAG

head apul.contigs.fasta 
>tig00000004 len=43693 reads=3 class=contig suggestRepeat=no suggestBubble=yes suggestCircular=no trim=0-43693
TACAATTTTAGAACACGGGACCAGCTTAGCATAATAGCTTCACCTTTCGTCTATCTAACTCTAGGAAGTTTTAATTTTTTCAAGTATTATAAAGGGCTCC
GTCGACTGTCAAAGATTTGCCTTTTCAAGCTCCAATGGAAACTGTAAAAGTTGCGTATTTTTACGAGCTTGAAACGCATCTTGGTATGCCCGATATCAAG
TCGAAATTAATATTGAGATAATCCTTTGGCCTTCTTCTATCAACATCTGAAATAAAAATCTGGGCAGTTTGAACGCGCTCTATCAAACATAACAGATTTG
AAGGTAGGTGATAACTTATTTGCATAATCTACGTTAACAAAAAGTCTATTTATAGAATGACTACTCGGCATATTTCTAACAGTGGTACTTCAGATACGTT
TTGATGGACTTATTATTCTGTCGTTTGTATTGTTTTCTTCAATTATTTAGCCTTAATAATTCCAAATAATAAAGAAATAAGGAAAGTCTTTGGTGTAAGT
CACACTCAAAAGGTGAGTTTCAACAGTTACTGAACACCCTTACGTATTAAACAGTCATTTCAATTTCCAGATTCTAACAGAAAATGTCAAATCGTTGTTT
TATAGTAGAAATCCATCTTCAAAAGTTATTCCCCGCTTATGCAGGCTTGATTCTCGCGGCTCTTTCCAGCTCGGTTTACAATATAAGACACCGGTGCAGA
TACCATTGAACTTGTAAACAATGTCACGCAAATTAAACTGTACTTCAATTTGCAAGCCATACAGCTTTAAGTCAGGTCTTTATTGAACTTTCTAAGTCAA
GGTTGGGGAATATAAAGATATTTTATTACCAGTATATTTTCGGTGAAAATTACAACGGATACATGTTATGGGCCTGTTCTTTAAACTCAGTTACATACAT
```

The unassembled file contains the reads and low-coverage contigs which couldn't be incorporated into the primary assembly. The contigs file contains the full assembly, including unique, repetitive, and bubble elements. What does this mean? Unsure, but the header line provides metadata on each sequence. 

```
len
   Length of the sequence, in bp.

reads
   Number of reads used to form the contig.

class
   Type of sequence.  Unassembled sequences are primarily low-coverage sequences spanned by a single read.

suggestRepeat
   If yes, sequence was detected as a repeat based on graph topology or read overlaps to other sequences.

suggestBubble
   If yes, sequence is likely a bubble based on potential placement within longer sequences.

suggestCircular
   If yes, sequence is likely circular.  The fasta line will have a trim=X-Y to indicate the non-redundant coordinates
```

If we are looking at this header: `>tig00000004 len=43693 reads=3 class=contig suggestRepeat=no suggestBubble=yes suggestCircular=no trim=0-43693`, the length is the contig is 43693 bp, 3 reads were used to form the contig, class is a contig, no repeats used to form the contig, the sequence is likely a bubble based on placement within longer sequences, the sequence is not circular, and the entire contig is non-redundant. Not sure what bubble means...

Because Canu didn't trim or correct the Hifi reads, I may need to use a different assembly tool. I'm going to try [Hifiasm](https://hifiasm.readthedocs.io/en/latest/index.html), which is a fast haplotype-resolved de novo assembler designed for PacBio HiFi reads. According to [Hifiasm github](https://github.com/chhylp123/hifiasm), here are some good reasons to use Hifiasm: 

- Hifiasm delivers high-quality telomere-to-telomere assemblies. It tends to generate longer contigs and resolve more segmental duplications than other assemblers.

- Hifiasm can purge duplications between haplotigs without relying on third-party tools such as purge_dups. Hifiasm does not need polishing tools like pilon or racon, either. This simplifies the assembly pipeline and saves running time.

- Hifiasm is fast. It can assemble a human genome in half a day and assemble a ~30Gb redwood genome in three days. No genome is too large for hifiasm.

- Hifiasm is trivial to install and easy to use. It does not required Python, R or C++11 compilers, and can be compiled into a single executable. The default setting works well with a variety of genomes.

If I use this tool, I may not need to pilon to polish the assembly. This tool isn't on the server, so I'll need to create a conda environment and install the package. 

```
cd /data/putnamlab/conda
mkdir hifiasm
cd hifiasm

module load Miniconda3/4.9.2
conda create --prefix /data/putnamlab/conda/hifiasm
conda activate /data/putnamlab/conda/hifiasm
conda install -c bioconda hifiasm
```

Once this package is installed, run code for hifiasm assembly. In the scripts folder: `nano hifiasm.sh`

```
#!/bin/bash -i
#SBATCH -t 500:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=500GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Apul_Genome/assembly/scripts
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

conda activate /data/putnamlab/conda/hifiasm

cd /data/putnamlab/jillashey/Apul_Genome/assembly/data

echo "Starting assembly with hifiasm" $(date)

hifiasm -o apul.hifiasm m84100_240128_024355_s2.hifi_reads.bc1029.fastq.fastq

echo "Assembly with hifiasm complete!" $(date)

conda deactivate
```

Submitted batch job 300534

### 20240213

Even though the reads were not trimmed or corrected with Canu, I am going to run [Busco](https://busco.ezlab.org/) on the output. This will provide information about how well the genome was assembled and its completeness based on evolutionarily informed expectations of gene content from near-universal single-copy orthologs. Danielle and Kevin have both run BUSCO before and used similar scripts but I think I'll adapt mine a little to fit my needs and personal preferences for code. 

From the [Busco user manual](https://busco.ezlab.org/busco_userguide.html#running-busco), the mandatory parameters are `-i`, which defines the input fasta file and `-m`, which sets the assessment mode (in our case, genome). Some recommended parameters incude `l` (specify busco lineage dataset; in our case, metazoans), `c` (specify number of cores to use), and `-o` (assigns specific label to output). 

In `/data/putnamlab/shared/busco/scripts`, the script `busco_init.sh` has information about the modules to load and in what order. Both Danielle and Kevin sourced this file specifically in their code, but I will probably just copy and paste the modules. In the same folder, they also used `busco-config.ini` as input for the `--config` flag in busco, which provides a config file as an alternative to command line parameters. I am not going to use this config file (yet), as Danielle and Kevin were assembling transcriptomes and I'm not sure what the specifics of the file are (or what they should be for genomes). In `/data/putnamlab/shared/busco/downloads/lineages/metazoa_odb10`, there is information about the metazoan database. 

In `/data/putnamlab/jillashey/Apul_Genome/assembly/scripts`, `nano busco_canu.sh`

```
#!/bin/bash 
#SBATCH -t 500:00:00
#SBATCH --nodes=1 --ntasks-per-node=15
#SBATCH --export=NONE
#SBATCH --mem=250GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Apul_Genome/assembly/scripts
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

module load BUSCO/5.2.2-foss-2020b
module load BLAST+/2.11.0-gompi-2020b
module load AUGUSTUS/3.4.0-foss-2020b
module load SEPP/4.4.0-foss-2020b
module load prodigal/2.6.3-GCCcore-10.2.0
module load HMMER/3.3.2-gompi-2020b

cd /data/putnamlab/jillashey/Apul_Genome/assembly/data

echo "Begin busco on canu-assembled fasta" $(date)

busco -i apul.contigs.fasta -m genome -l /data/putnamlab/shared/busco/downloads/lineages/metazoa_odb10 -c 15 -o apul.busco.canu

echo "busco complete for canu-assembled fasta" $(date)
```

Submitted batch job 301588. This failed and gave me some errors. This one seemed to have been the fatal one: `Message: BatchFatalError(AttributeError("'NoneType' object has no attribute 'remove_tmp_files'"))`. Danielle ran into a similar error in her busco [code](https://github.com/daniellembecker/DanielleBecker_Lab_Notebook/blob/master/_posts/2023-08-31-Acropora-pulchra-denovo-transcriptome.md) so I am going to try to set the `--config` file as `"$EBROOTBUSCO/config/config.ini"`. 

In the script, `nano busco_canu.sh`

```
#!/bin/bash 
#SBATCH -t 500:00:00
#SBATCH --nodes=1 --ntasks-per-node=15
#SBATCH --export=NONE
#SBATCH --mem=250GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Apul_Genome/assembly/scripts
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

#module load BUSCO/5.2.2-foss-2020b
#module load BLAST+/2.11.0-gompi-2020b
#module load AUGUSTUS/3.4.0-foss-2020b
#module load SEPP/4.4.0-foss-2020b
#module load prodigal/2.6.3-GCCcore-10.2.0
#module load HMMER/3.3.2-gompi-2020b

cd /data/putnamlab/jillashey/Apul_Genome/assembly/data

echo "Begin busco on canu-assembled fasta" $(date)

source "/data/putnamlab/shared/busco/scripts/busco_init.sh"  # sets up the modules required for this in the right order

busco --config "$EBROOTBUSCO/config/config.ini" -f -c 15 --long -i apul.contigs.fasta -m genome -l /data/putnamlab/shared/busco/downloads/lineages/metazoa_odb10 -o apul.busco.canu

echo "busco complete for canu-assembled fasta" $(date)
```

Submitted batch job 301594. Failed, same error as before. Going to try copying Kevin and Danielle code directly, even though its a little messy and confusing with paths. 

In the script, `nano busco_canu.sh`

```
#!/bin/bash 
#SBATCH -t 500:00:00
#SBATCH --nodes=1 --ntasks-per-node=15
#SBATCH --export=NONE
#SBATCH --mem=250GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Apul_Genome/assembly/scripts
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

echo "Begin busco on canu-assembled fasta" $(date)

labbase=/data/putnamlab
busco_shared="${labbase}/shared/busco"
[ -z "$query" ] && query="${labbase}/jillashey/Apul_Genome/assembly/data/apul.contigs.fasta" # set this to the query (genome/transcriptome) you are running
[ -z "$db_to_compare" ] && db_to_compare="${busco_shared}/downloads/lineages/metazoa_odb10"

source "${busco_shared}/scripts/busco_init.sh"  # sets up the modules required for this in the right order

# This will generate output under your $HOME/busco_output
cd "${labbase}/${Apul_Genome/assembly/data}"
busco --config "$EBROOTBUSCO/config/config.ini"  -f -c 20 --long -i "${query}" -l metazoa_odb10 -o apul.busco.canu -m genome

echo "busco complete for canu-assembled fasta" $(date)
```

Submitted batch job 301599. This appears to have worked! Took about an hour to run. This is the primary result in the out file: 

```
# BUSCO version is: 5.2.2 
# The lineage dataset is: metazoa_odb10 (Creation date: 2024-01-08, number of genomes: 65, number of BUSCOs: 954)
# Summarized benchmarking in BUSCO notation for file /data/putnamlab/jillashey/Apul_Genome/assembly/data/apul.contigs.fasta
# BUSCO was run in mode: genome
# Gene predictor used: metaeuk

        ***** Results: *****

        C:94.4%[S:9.4%,D:85.0%],F:2.7%,M:2.9%,n:954        
        901     Complete BUSCOs (C)                        
        90      Complete and single-copy BUSCOs (S)        
        811     Complete and duplicated BUSCOs (D)         
        26      Fragmented BUSCOs (F)                      
        27      Missing BUSCOs (M)                         
        954     Total BUSCO groups searched                
```

We have 94.4% completeness with this assembly but 85% complete and duplicated BUSCOs. The busco manual says this on high levels of duplication: "BUSCO completeness results make sense only in the context of the biology of your organism. You have to understand whether missing or duplicated genes are of biological or technical origin. For instance, a high level of duplication may be explained by a recent whole duplication event (biological) or a chimeric assembly of haplotypes (technical). Transcriptomes and protein sets that are not filtered for isoforms will lead to a high proportion of duplicates. Therefore you should filter them before a BUSCO analysis". 

Danielle also got a high number (78.9%) of duplicated BUSCOs in her [de novo transcriptome](https://github.com/daniellembecker/DanielleBecker_Lab_Notebook/blob/master/_posts/2023-08-31-Acropora-pulchra-denovo-transcriptome.md) of Apulchra, but Kevin got much less duplication (6.9%) in his Past [transcriptome assembly](https://github.com/kevinhwong1/KevinHWong_Notebook/blob/master/_posts/2021-01-04-20210104-BUSCO-on-P.-astreoides-transcriptome-assembly.md). I need to ask Danielle if she ended up using her Trinity results (which had a high duplication percentage) for her alignment for Apul. I also need to ask her if she thinks the high duplication percentage is biologically meaningful. 

Might be worth running [HiFiAdapter Filt](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-022-08375-1)
