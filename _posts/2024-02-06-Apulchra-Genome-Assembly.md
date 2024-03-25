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

### 20240215 

Last night, the hifiasm job failed after almost 2 days but the email says PREEMPTED, ExitCode0. Two minutes after the job failed, job 300534 started again on the server and it says its a hifiasm job...I did not start this job myself, not sure what happened. Looking on the server now, hifiasm is running but has only been running for about 18 hours (as of 2pm today). It's the same job number though which is strange. 

### 20240220

Hifiasm job is still running after ~5 days. In the meantime, I'm going to run [HiFiAdapterFilt](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-022-08375-1#Sec6), which is an adapter filtering command for PacBio HiFi data. On the github [page](https://github.com/sheinasim/HiFiAdapterFilt/tree/master), it says that the tool converts .bam to .fastq and removes reads with remnant PacBio adapter sequences. Required dependencies are BamTools and BLAST+; optional dependencies are NCBI FCS Adaptor and pigz. It looks like I'll need to use the original bam file instead of the converted fastq file. 

The github says I should add the script and database to my path using: 

```
export PATH=$PATH:[PATH TO HiFiAdapterFilt]
export PATH=$PATH:[PATH TO HiFiAdapterFilt]/DB
```

I will do this in the script for the adapter filt code that I write myself. In the scripts folder, make a folder for hifi information

```
mkdir HiFiAdapterFilt
cd HiFiAdapterFilt
``` 

I need to make a script for the hifi adapter [code](https://raw.githubusercontent.com/sheinasim/HiFiAdapterFilt/master/hifiadapterfilt.sh). In the `scripts/HiFiAdapterFilt` folder, I copy and pasted the linked code into `hifiadapterfilt.sh`. Make a folder for pacbio databases and copy in the [db information](https://github.com/sheinasim/HiFiAdapterFilt/blob/master/DB/pacbio_vectors_db) from the github. 

```
mkdir DB
cd DB
nano pacbio_vectors_db
>gnl|uv|NGB00972.1:1-45 Pacific Biosciences Blunt Adapter
ATCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGAT
>gnl|uv|NGB00973.1:1-35 Pacific Biosciences C2 Primer
AAAAAAAAAAAAAAAAAATTAACGGAGGAGGAGGA
```

In the `scripts/HiFiAdapterFilt` folder: `nano hifiadapterfilt_JA.sh`

```
#!/bin/bash 
#SBATCH -t 500:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=500GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Apul_Genome/assembly/scripts/HiFiAdapterFilt
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

module load GCCcore/11.3.0 # need this to resolve conflicts between GCCcore/8.3.0 and loaded GCCcore/11.3.0
module load BamTools/2.5.1-GCC-8.3.0 
module load BLAST+/2.9.0-iimpi-2019b 

cd /data/putnamlab/jillashey/Apul_Genome/assembly/scripts/HiFiAdapterFilt

echo "Setting paths" $(date)

export PATH=$PATH:[/data/putnamlab/jillashey/Apul_Genome/assembly/scripts/HiFiAdapterFilt] # path to original script
export PATH=$PATH:[/data/putnamlab/jillashey/Apul_Genome/assembly/scripts/HiFiAdapterFilt]/DB # path to db info 

echo "Paths set, starting adapter filtering" $(date)

bash hifiadapterfilt.sh -p /data/putnamlab/KITT/hputnam/20240129_Apulchra_Genome_LongRead/m84100_240128_024355_s2.hifi_reads.bc1029 -l 44 -m 97 -o /data/putnamlab/jillashey/Apul_Genome/assembly/data

echo "Completing adapter filtering" $(date)
```

The `-l` and `-m` refer to the minimum length of adapter match to remove and the minumum percent match of adapter to remove, respectively. I left them as the default settings for now. Submitted batch job 303636. Giving me this error: 

```
hifiadapterfilt.sh: line 63: /data/putnamlab/KITT/hputnam/20240129_Apulchra_Genome_LongRead/m84100_240128_024355_s2.hifi_reads.bc1029.temp_file_list: Permission denied
cat: /data/putnamlab/KITT/hputnam/20240129_Apulchra_Genome_LongRead/m84100_240128_024355_s2.hifi_reads.bc1029.bam.temp_file_list: No such file or directory
cat: /data/putnamlab/KITT/hputnam/20240129_Apulchra_Genome_LongRead/m84100_240128_024355_s2.hifi_reads.bc1029.bam.temp_file_list: No such file or directory
```

Giving me permission denied to write in the folder? I'll sym link the bam file to my data folder

```
cd /data/putnamlab/jillashey/Apul_Genome/assembly/data
ln -s /data/putnamlab/KITT/hputnam/20240129_Apulchra_Genome_LongRead/m84100_240128_024355_s2.hifi_reads.bc1029.bam
```

Editing script so that the prefix is connecting with the sym linked file. Submitted batch job 303637. Still giving me the same error. Hollie may need to give me permission to write and access files in that specific folder. 

### 20240221

Probably need to run [haplomerger2](https://github.com/mapleforest/HaploMerger2/releases/), which is installed on the server already. 

I'm also looking at the Canu FAQs to see if there is any info about using PacBio HiFi reads. Under the question ["What parameters should I use for my reads?"](https://canu.readthedocs.io/en/latest/faq.html), they have this info:

```
The defaults for -pacbio-hifi should work on this data. There is still some variation in data quality between samples. If you have poor continuity, it may be because the data is lower quality than expected. Canu will try to auto-adjust the error thresholds for this (which will be included in the report). If that still doesn’t give a good assembly, try running the assembly with -untrimmed. You will likely get a genome size larger than you expect, due to separation of alleles. See My genome size and assembly size are different, help! for details on how to remove this duplication.
```

When I look at the question ["My genome size and assembly size are different, help!"](https://canu.readthedocs.io/en/latest/faq.html#my-genome-size-and-assembly-size-are-different-help), it says that this difference could be due to a heterozygous genome where the assembly separated some loci or the previous estimate is incorrect. They recommended running BUSCO to check completeness of the assembly (which I already did) and using [purge_dups](https://github.com/dfguan/purge_dups?tab=readme-ov-file#usg) to remove duplication. I will look into this. 

Next steps 

- Run Canu with `-untrimmed` option 
- Run purge_dups - targets the removal of duplicated sequences to enhance overall quality of assembly 
- Run haplomerger - merges haplotypes to addresws heterozygosity 

In the scripts folder: `nano canu_untrimmed.sh`

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

#echo "Unzip paco-bio fastq file" $(date)

#gunzip m84100_240128_024355_s2.hifi_reads.bc1029.fastq.fastq.gz

echo "Starting assembly w/ untrimmed flag" $(date)

canu -p apul.canu.untrimmed -d /data/putnamlab/jillashey/Apul_Genome/assembly/data genomeSize=475m -raw -pacbio-hifi m84100_240128_024355_s2.hifi_reads.bc1029.fastq.fastq

echo "Canu assembly complete" $(date)
```

Submitted batch job 303660

On the purge_dups [github](https://github.com/dfguan/purge_dups?tab=readme-ov-file#usg), they say to install using the following: 

```
git clone https://github.com/dfguan/purge_dups.git
cd purge_dups/src && make

# only needed if running run_purge_dups.py
git clone https://github.com/dfguan/runner.git
cd runner && python3 setup.py install --user
```

Cloned both into the assembly folder (ie `/data/putnamlab/jillashey/Apul_Genome/assembly/`). 

First, use pd_config.py to generate a configuration file. Here's possible usage: 

```
usage: pd_config.py [-h] [-s SRF] [-l LOCD] [-n FN] [--version] ref pbfofn

generate a configuration file in json format

positional arguments:
  ref                   reference file in fasta/fasta.gz format
  pbfofn                list of pacbio file in fastq/fasta/fastq.gz/fasta.gz format (one absolute file path per line)

optional arguments:
  -h, --help            show this help message and exit
  -s SRF, --srfofn SRF  list of short reads files in fastq/fastq.gz format (one record per line, the
                        record is a tab splitted line of abosulte file path
                        plus trimmed bases, refer to
                        https://github.com/dfguan/KMC) [NONE]
  -l LOCD, --localdir LOCD
                        local directory to keep the reference and lists of the
                        pacbio, short reads files [.]
  -n FN, --name FN      output config file name [config.json]
  --version             show program's version number and exit

# Example 
./scripts/pd_config.py -l iHelSar1.pri -s 10x.fofn -n config.iHelSar1.PB.asm1.json ~/vgp/release/insects/iHelSar1/iHlSar1.PB.asm1/iHelSar1.PB.asm1.fa.gz pb.fofn
```

I need to make a list of the pacbio files that I'll use in the script and put it in a file called pb.fofn. 

```
for filename in *.fastq; do echo $PWD/$filename; done > pb.fofn
```

In my scripts folder: `nano pd_config_apul.py`

```
#!/bin/bash 
#SBATCH -t 100:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=250GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Apul_Genome/assembly/scripts
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

module load Python/3.9.6-GCCcore-11.2.0 # do i need python?

echo "Making config file for purge dups scripts" $(date)

/data/putnamlab/jillashey/Apul_Genome/assembly/purge_dups/scripts/pd_config.py -l /data/putnamlab/jillashey/Apul_Genome/assembly/data -n config.apul.canu.json /data/putnamlab/jillashey/Apul_Genome/assembly/data/m84100_240128_024355_s2.hifi_reads.bc1029.fastq.fastq /data/putnamlab/jillashey/Apul_Genome/assembly/data/pb.fofn

echo "Config file complete" $(date)
```

Submitted batch job 303661. Ran in 1 second. Got this in the error file: 

```
cp: ‘/data/putnamlab/jillashey/Apul_Genome/assembly/data/m84100_240128_024355_s2.hifi_reads.bc1029.fastq.fastq’ and ‘/data/putnamlab/jillashey/Apul_Genome/assembly/data/m84100_240128_024355_s2.hifi_reads.bc1029.fastq.fastq’ are the same file
cp: ‘/data/putnamlab/jillashey/Apul_Genome/assembly/data/pb.fofn’ and ‘/data/putnamlab/jillashey/Apul_Genome/assembly/data/pb.fofn’ are the same file
```

But it did generate a config file, which looks like this: 

```
{
  "cc": {
    "fofn": "/data/putnamlab/jillashey/Apul_Genome/assembly/data/pb.fofn",
    "isdip": 1,
    "core": 12,
    "mem": 20000,
    "queue": "normal",
    "mnmp_opt": "",
    "bwa_opt": "",
    "ispb": 1,
    "skip": 0
  },
  "sa": {
    "core": 12,
    "mem": 10000,
    "queue": "normal"
  },
  "busco": {
    "core": 12,
    "mem": 20000,
    "queue": "long",
    "skip": 0,
    "lineage": "mammalia",
    "prefix": "m84100_240128_024355_s2.hifi_reads.bc1029.fastq_purged",
    "tmpdir": "busco_tmp"
  },
  "pd": {
    "mem": 20000,
    "queue": "normal"
  },
  "gs": {
    "mem": 10000,
    "oe": 1
  },
  "kcp": {
    "core": 12,
    "mem": 30000,
    "fofn": "",
    "prefix": "m84100_240128_024355_s2.hifi_reads.bc1029.fastq_purged_kcm",
    "tmpdir": "kcp_tmp",
    "skip": 1
  },
  "ref": "/data/putnamlab/jillashey/Apul_Genome/assembly/data/m84100_240128_024355_s2.hifi_reads.bc1029.fastq.fastq",
  "out_dir": "m84100_240128_024355_s2.hifi_reads.bc1029.fastq"
}
```

My config file looks basically the same as the example one on the github. Manually edited the config file so that the out_dir was `/data/putnamlab/jillashey/Apul_Genome/assembly/data/`. Now the purging can begin using `run_purge_dups.py`. Here's possible usage: 

```
usage: run_purge_dups.py [-h] [-p PLTFM] [-w WAIT] [-r RETRIES] [--version]
                         config bin_dir spid

purge_dups wrapper

positional arguments:
  config                configuration file
  bin_dir               directory of purge_dups executable files
  spid                  species identifier

optional arguments:
  -h, --help            show this help message and exit
  -p PLTFM, --platform PLTFM
                        workload management platform, input bash if you want to run locally
  -w WAIT, --wait WAIT  <int> seconds sleep intervals
  -r RETRIES, --retries RETRIES
                        maximum number of retries
  --version             show program's version number and exit
  
# Example 
python scripts/run_purge_dups.py config.iHelSar1.json src iHelSar1
```

In my scripts folder: `nano run_purge_dups_apul.py`

```
#!/bin/bash 
#SBATCH -t 100:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=250GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Apul_Genome/assembly/scripts
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

module load Python/3.9.6-GCCcore-11.2.0 # do i need python?

echo "Starting to purge duplications" $(date)

/data/putnamlab/jillashey/Apul_Genome/assembly/purge_dups/scripts/run_purge_dups.py /data/putnamlab/jillashey/Apul_Genome/assembly/scripts/config.apul.canu.json /data/putnamlab/jillashey/Apul_Genome/assembly/purge_dups/src apul

echo "Duplication purge complete" $(date)
```

Submitted batch job 303662. Failed immediately with this error: 

```
Traceback (most recent call last):
  File "/data/putnamlab/jillashey/Apul_Genome/assembly/purge_dups/scripts/run_purge_dups.py", line 3, in <module>
    from runner.manager import manager
ModuleNotFoundError: No module named 'runner'
```

I installed runner but the code is not seeing it...where am I supposed to put it? inside of the purge_dups github? Okay going to move `runner` folder inside of the `purge_dups` folder. 

```
cd /data/putnamlab/jillashey/Apul_Genome/assembly
mv runner/ purge_dups/scripts/
```

Submitting job again, Submitted batch job 303663. Got the same error. In the `run_purge_dups.py` script itself, the first few lines are: 

```
#!/usr/bin/env python3

from runner.manager import manager
from runner.hpc import hpc
from multiprocessing import Process, Pool
import sys, os, json
import argparse
```

So it isn't seeing that the runner module is there. This [issue](https://github.com/dfguan/purge_dups/issues/96) and this [issue](https://github.com/dfguan/purge_dups/issues/83) on the github were reported but never really answered in a clear way. Will have to look into this more. 

In other news, the canu script finished running but looks like it failed. This is the bottom of the error message: 

```
ERROR:
ERROR:  Failed with exit code 139.  (rc=35584)
ERROR:

ABORT:
ABORT: canu 2.2
ABORT: Don't panic, but a mostly harmless error occurred and Canu stopped.
ABORT: Try restarting.  If that doesn't work, ask for help.
ABORT:
ABORT:   failed to configure the overlap store.
ABORT:
ABORT: Disk space available:  8477.134 GB
ABORT:
ABORT: Last 50 lines of the relevant log file (unitigging/apul.canu.untrimmed.ovlStore.config.err):
ABORT:
ABORT:   
ABORT:   Finding number of overlaps per read and per file.
ABORT:   
ABORT:      Moverlaps
ABORT:   ------------ ----------------------------------------
ABORT:   
ABORT:   Failed with 'Segmentation fault'; backtrace (libbacktrace):
ABORT:
```

Unsure what it means...

### 20240301

BIG NEWS!!!!! This week, a paper came out that assembled and annotated the *Orbicella faveolata* genome using PacBio HiFi reads ([Young et al. 2024](https://link.springer.com/article/10.1186/s12864-024-10092-w?utm_source=rct_congratemailt&utm_medium=email&utm_campaign=oa_20240229&utm_content=10.1186/s12864-024-10092-w#Sec12334225451))!!!!!!! The [github](https://github.com/benyoung93/orbicella_faveolata_pacbio_genome_transcriptome/blob/main) for this paper has a detailed pipeline for how the genome was put together. Since I am also using HiFi reads, I will be following their methodology! I am using this [pipeline](https://github.com/benyoung93/orbicella_faveolata_pacbio_genome_transcriptome/blob/main/ofav_genome_pipeline.Rmd) starting at line 260. 

I changed the file from bam to fastq, but now I need to change it to fasta with [`seqtk`](https://github.com/lh3/seqtk). In the scripts folder: `nano seqtk.sh`

```
#!/bin/bash 
#SBATCH -t 100:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --export=NONE
#SBATCH --mem=500GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Apul_Genome/assembly/scripts
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

module load seqtk/1.3-GCC-9.3.0

echo "Convert PacBio fastq file to fasta file" $(date)

cd /data/putnamlab/jillashey/Apul_Genome/assembly/data

seqtk seq -a m84100_240128_024355_s2.hifi_reads.bc1029.fastq.fastq > m84100_240128_024355_s2.hifi_reads.bc1029.fasta

echo "Fastq to fasta complete! Summarize read lengths" $(date)

awk '/^>/{printf("%s\t",substr($0,2));next;} {print length}' m84100_240128_024355_s2.hifi_reads.bc1029.fasta > rr_read_lengths.txt

echo "Read length summary complete" $(date)
```

Submitted batch job 304257. 

In R, I looked at the data to quantify length for each read. See code [here](https://github.com/hputnam/Apulchra_genome/blob/main/scripts/genome_analysis.Rmd). 

```{r, echo=F}read.table(file = "../data/rr_read_lengths.txt",            header = F) %>%   dplyr::rename("hifi_read_name" = 1,          "length" = 2) -> hifi_read_lengthnrow(hifi_read_length) # 5,898,386 total readsmean(hifi_read_length$length) # mean length of reads is 13,424.64sum(hifi_read_length$length) #length sum 79,183,709,778. Will need this for the NCBI submission```

Make histogram for read bins from raw hifi data
```{r, echo = F}ggplot(data = hifi_read_length,        aes(x = length, fill = "blue")) +  geom_histogram(binwidth = 2000) +   labs(x = "Raw Read Length", y = "Count", title = "Histogram of Raw HiFi Read Lengths") +   scale_fill_manual(values = c("blue")) +   scale_y_continuous(labels = function(x) format(x, scientific = FALSE))```

![](https://raw.githubusercontent.com/hputnam/Apulchra_genome/main/output/rr_length_histogram.png)

### 20240303 

My next step is to remove any contaminant reads from the raw hifi reads. From Young et al. 2024: "Raw HiFi reads first underwent a contamination screening, following the methodology in [68], using BLASTn [32, 68] against the assembled mitochondrial O. faveolata genome and the following databases: common eukaryote contaminant sequences (ftp.ncbi.nlm.nih. gov/pub/kitts/contam_in_euks.fa.gz), NCBI viral (ref_ viruses_rep_genomes) and prokaryote (ref_prok_rep_ genomes) representative genome sets". 

I tried to run the `update_blastdb.pl` script (included in the blast program) with the `BLAST+/2.13.0-gompi-2022a` module but I got this error: 

```
Can't locate Archive/Tar.pm in @INC (@INC contains: /usr/local/lib64/perl5 /usr/local/share/perl5 /usr/lib64/perl5/vendor_perl /usr/share/perl5/vendor_perl /usr/lib64/perl5 /usr/share/perl5 .) at /opt/software/BLAST+/2.13.0-gompi-2022a/bin/update_blastdb.pl line 41.
BEGIN failed--compilation aborted at /opt/software/BLAST+/2.13.0-gompi-2022a/bin/update_blastdb.pl line 41.
```

Not sure what this means...will email Kevin Bryan to ask about it, as I don't want to mess with anything on the installed modules. I did download the `contam_in_euks.fa.gz` db to my computer so I'm going to copy it to Andromeda. This file is considerably smaller than the viral or prok dbs. 

Make a database folder in the Apul genome folder. 

```
cd /data/putnamlab/jillashey/Apul_Genome
mkdir dbs 
cd dbs

zgrep -c ">" contam_in_euks.fa 
3554
```

Now I ran run a script that blasts the pacbio fasta against these sequences. In `/data/putnamlab/jillashey/Apul_Genome/assembly/scripts`, `nano blast_contam_euk.sh`

```
#!/bin/bash 
#SBATCH -t 100:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --export=NONE
#SBATCH --mem=500GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Apul_Genome/assembly/scripts
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

module load BLAST+/2.13.0-gompi-2022a

echo "BLASTing hifi fasta against eukaryote contaminant sequences" $(date)

cd /data/putnamlab/jillashey/Apul_Genome/assembly/data


blastn -query m84100_240128_024355_s2.hifi_reads.bc1029.fasta -subject /data/putnamlab/jillashey/Apul_Genome/dbs/contam_in_euks.fa -task megablast -outfmt 6 -evalue 4 -perc_identity 90 -num_threads 15 -out contaminant_hits_euks_rr.txt

echo "BLAST complete, remove contaminant seqs from hifi fasta" $(date)

awk '{ if( ($4 >= 50 && $4 <= 99 && $3 >=98 ) ||
         ($4 >= 100 && $4 <= 199 && $3 >= 94 ) ||
         ($4 >= 200 && $3 >= 90) )  {print $0}
    }' contaminant_hits_euks_rr.txt > contaminants_pass_filter_euks_rr.txt

echo "Contaminant seqs removed from hifi fasta" $(date)
```

Submitted batch job 304389. Finished in about 2.5 hours. Looked at the output in R (code [here](https://github.com/hputnam/Apulchra_genome/blob/main/scripts/genome_analysis.Rmd)). 

### 20240304

Emailed Kevin Bryan this morning and asked if he knew anything about why the `update_blastdb.pl` wasn't working. Still waiting to hear back from him. 

Kevin Bryan also emailed me this morning about my hifiasm job that has been running for 18 days and said: "This job has been running for 18 days. I just took a look at it and it appears you didn’t specify `-t $SLURM_CPUS_ON_NODE` (and also `#SBATCH --exclusive`) to make use of all of the CPU cores on the node. You might want to consider re-submitting this job with those parameters. Because the nodes generally have 36 cores, it should be able to catch up to where it is now in a little over half a day, assuming perfect scaling." 

I need to add those parameters into the hifiasm code, so cancelling the `300534` job. In `/data/putnamlab/jillashey/Apul_Genome/assembly/scripts`, `nano hifiasm.sh`:

```
#!/bin/bash -i
#SBATCH -t 30-00:00:00
#SBATCH --nodes=1 --ntasks-per-node=36
#SBATCH --export=NONE
#SBATCH --mem=500GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH --exclusive
#SBATCH -D /data/putnamlab/jillashey/Apul_Genome/assembly/scripts
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

conda activate /data/putnamlab/conda/hifiasm

cd /data/putnamlab/jillashey/Apul_Genome/assembly/data

echo "Starting assembly with hifiasm" $(date)

hifiasm -o apul.hifiasm m84100_240128_024355_s2.hifi_reads.bc1029.fastq.fastq -t 36

echo "Assembly with hifiasm complete!" $(date)

conda deactivate
```

Submitted batch job 304463

### 20240304

Response from Kevin Bryan about viral and prok blast databases: "Ok, I downloaded those databases, and actually consolidated the rest of them into `/data/shared/ncbi-db/`, under which is a directory for today, and then there will be a new one next Sunday and following Sundays. There’s a file `/data/shared/ncbi-db/.ncbirc` that gets updated to point to the current directory, which the blast* tools will automatically pick up, so you can just do `-db ref_prok_rep_ genomes`, for example.

For other tools that can read the ncbi databases, you can use `blastdb_path -db  ref_viruses_rep_genomes` to get the path, although for some reason with the nr database you need to specify `-dbtype prot`, i.e., `blastdb_path -db  nr -dbtype prot`.

The reason for the extra complication is because otherwise a job that runs while the database is being updated may fail or return strange results. The dated directories should resolve this issue. Note that Unity blast-plus modules work in a similar way with a different path, `/datasets/bio/ncbi-db`."

Amazing! Now I can move forward with the blasting against viral and prok genomes. In the scripts folder: `nano blastn_viral.sh`

```
#!/bin/bash 
#SBATCH -t 30-00:00:00
#SBATCH --nodes=1 --ntasks-per-node=36
#SBATCH --export=NONE
#SBATCH --mem=500GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH --exclusive
#SBATCH -D /data/putnamlab/jillashey/Apul_Genome/assembly/scripts
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

module load BLAST+/2.13.0-gompi-2022a

cd /data/putnamlab/jillashey/Apul_Genome/assembly/data

echo "Blasting hifi reads against viral genomes to look for contaminants" $(date)

blastn -query m84100_240128_024355_s2.hifi_reads.bc1029.fasta -db ref_viruses_rep_genomes -outfmt 6 -evalue 1e-4 -perc_identity 90 -out viral_contaminant_hits_rr.txt

echo "Blast complete!" $(date)
```

Submitted batch job 304500

In the scripts folder: `nano blastn_prok.sh`

```
#!/bin/bash 
#SBATCH -t 30-00:00:00
#SBATCH --nodes=1 --ntasks-per-node=36
#SBATCH --export=NONE
#SBATCH --mem=500GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH --exclusive
#SBATCH -D /data/putnamlab/jillashey/Apul_Genome/assembly/scripts
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

module load BLAST+/2.13.0-gompi-2022a

cd /data/putnamlab/jillashey/Apul_Genome/assembly/data

echo "Blasting hifi reads against prokaryote genomes to look for contaminants" $(date)

blastn -query m84100_240128_024355_s2.hifi_reads.bc1029.fasta -db ref_prok_rep_genomes -outfmt 6 -evalue 1e-4 -perc_identity 90 -out prok_contaminant_hits_rr.txt

echo "Blast complete!" $(date)
```

Submitted batch job 304502

### 20240306

Making list of all software programs that Young et al. 2024 used and if they are on Andromeda 

- blastn 
	- On Andromeda? YES
- Meryl
	- On Andromeda? NO
- Genome-Scope2
	- On Andromeda? NO
- Hifiasm
	- On Andromeda? NO but I added it to `putnamlab` via conda
- Quast
	- On Andromeda? YES
- Busco
	- On Andromeda? YES
- Merqury
	- On Andromeda? NO
- RepeatModeler2
	- On Andromeda? NO
- Repeat-Masker 
	- On Andromeda? YES
- TeloScafs 
	- On Andromeda? NO
- PASA
	- On Andromeda? NO
- funnannotate 
	- On Andromeda? NO
- Augustus
	- On Andromeda? YES
- GeneMark-ES/ET
	- On Andromeda? YES but only GeneMark-ET
- snap 
	- On Andromeda? YES
- glimmerhmm
	- On Andromeda? NO
- Evidence Modeler
	- On Andromeda? NO
- tRNAscan-SE
	- On Andromeda? YES
- Trinity 
	- On Andromeda? Yes
- InterproScan 
	- On Andromeda? YES
	
### 20240311

Hifiasm (with unfiltered reads) finished running over the weekend and the prok blast script preemptively ended and then restarted in the early hours of this morning. I think this might be because I am not making use of all cores on the node (similar to my earlier hifiasm script). I cancelled the prok blast job (`304502`) and edited the script so that it includes the flag `-num_threads 36`. Submitted batch job 305351

It created many files:

```
-rw-r--r--. 1 jillashey  19G Mar  7 23:58 apul.hifiasm.ec.bin
-rw-r--r--. 1 jillashey  47G Mar  8 00:08 apul.hifiasm.ovlp.source.bin
-rw-r--r--. 1 jillashey  17G Mar  8 00:12 apul.hifiasm.ovlp.reverse.bin
-rw-r--r--. 1 jillashey 1.2G Mar  8 01:37 apul.hifiasm.bp.r_utg.gfa
-rw-r--r--. 1 jillashey  21M Mar  8 01:37 apul.hifiasm.bp.r_utg.noseq.gfa
-rw-r--r--. 1 jillashey 8.6M Mar  8 01:41 apul.hifiasm.bp.r_utg.lowQ.bed
-rw-r--r--. 1 jillashey 1.1G Mar  8 01:42 apul.hifiasm.bp.p_utg.gfa
-rw-r--r--. 1 jillashey  21M Mar  8 01:42 apul.hifiasm.bp.p_utg.noseq.gfa
-rw-r--r--. 1 jillashey 8.2M Mar  8 01:46 apul.hifiasm.bp.p_utg.lowQ.bed
-rw-r--r--. 1 jillashey 506M Mar  8 01:47 apul.hifiasm.bp.p_ctg.gfa
-rw-r--r--. 1 jillashey  11M Mar  8 01:47 apul.hifiasm.bp.p_ctg.noseq.gfa
-rw-r--r--. 1 jillashey 2.0M Mar  8 01:49 apul.hifiasm.bp.p_ctg.lowQ.bed
-rw-r--r--. 1 jillashey 469M Mar  8 01:50 apul.hifiasm.bp.hap1.p_ctg.gfa
-rw-r--r--. 1 jillashey 9.9M Mar  8 01:50 apul.hifiasm.bp.hap1.p_ctg.noseq.gfa
-rw-r--r--. 1 jillashey 2.0M Mar  8 01:52 apul.hifiasm.bp.hap1.p_ctg.lowQ.bed
-rw-r--r--. 1 jillashey 468M Mar  8 01:52 apul.hifiasm.bp.hap2.p_ctg.gfa
-rw-r--r--. 1 jillashey 9.9M Mar  8 01:52 apul.hifiasm.bp.hap2.p_ctg.noseq.gfa
-rw-r--r--. 1 jillashey 1.9M Mar  8 01:54 apul.hifiasm.bp.hap2.p_ctg.lowQ.bed
```

This [page](https://hifiasm.readthedocs.io/en/latest/interpreting-output.html#interpreting-output) gives a brief overview of the hifiasm output files, which is super helpful. It generates the assembly graphs in Graphical Fragment Assembly ([GFA](https://github.com/GFA-spec/GFA-spec)) format. 

- prefix.r_utg.gfa: haplotype-resolved raw unitig graph. This graph keeps all haplotype information.
	- A unitig is a portion of a contig. It is a nondisputed and assembled group of fragments. A contiguous sequence of ordered unitigs is a contig, and a single unitig can be in multiple contigs.
- prefix.p_utg.gfa: haplotype-resolved processed unitig graph without small bubbles. Small bubbles might be caused by somatic mutations or noise in data, which are not the real haplotype information. Hifiasm automatically pops such small bubbles based on coverage. The option --hom-cov affects the result. See homozygous coverage setting for more details. In addition, the option -p forcedly pops bubbles.
	- Confused about the bubbles, but it looks like a medium level (what is a "medium" level"?) of heterozygosity will result in bubbles (see image in this [post](https://medium.com/pacbio/hifi-assembler-series-part-1-hifiasm-a-fast-haplotype-resolved-genome-assembler-bd2f30dab571))
	- Homozygous coverage refers to coverage threshold for homozygous reads. Hifiasm prints it as: `[M::purge_dups] homozygous read coverage threshold: X`. If it is not around homozygous coverage, the final assembly might be either too large or too small. 
- prefix.p_ctg.gfa: assembly graph of primary contigs. This graph includes a complete assembly with long stretches of phased blocks.
	- From my understanding based on this [post](https://lh3.github.io/2021/04/17/concepts-in-phased-assemblies) discussing concepts in phased assemblies, a phased assembly identifies different alleles
- prefix.a_ctg.gfa: assembly graph of alternate contigs. This graph consists of all contigs that are discarded in primary contig graph.
	- There were none of these files in my output. Does this mean all contigs created were used in the final assembly? 
- prefix.hap*.p_ctg.gfa: phased contig graph. This graph keeps the phased contigs for haplotype 1 and haplotype 2. 

I believe this file (`apul.hifiasm.bp.p_ctg.gfa`) contains the sequence information for the assembled contigs. `zgrep -c "S" apul.hifiasm.bp.p_ctg.gfa` showed that there were 188 Segments, or continuous sequences, in this assembly, meaning there are 188 contigs (to my understanding). This file (`apul.hifiasm.bp.p_ctg.noseq.gfa`) contains information about the reads used to construct the contigs in a plain text format. 

```
head apul.hifiasm.bp.p_ctg.noseq.gfa
S	ptg000001l	*	LN:i:21642937	rd:i:82
A	ptg000001l	0	+	m84100_240128_024355_s2/165481517/ccs	0	25806	id:i:5084336	HG:A:a
A	ptg000001l	6512	+	m84100_240128_024355_s2/221910786/ccs	0	21626	id:i:5609370	HG:A:a
A	ptg000001l	7287	-	m84100_240128_024355_s2/128778986/ccs	0	24026	id:i:2691352	HG:A:a
A	ptg000001l	7493	+	m84100_240128_024355_s2/262476615/ccs	0	30540	id:i:287120	HG:A:a
A	ptg000001l	15042	+	m84100_240128_024355_s2/191824085/ccs	0	27619	id:i:2783058	HG:A:a
A	ptg000001l	16616	-	m84100_240128_024355_s2/29886096/ccs	0	28664	id:i:2336674	HG:A:a
A	ptg000001l	17527	-	m84100_240128_024355_s2/37553051/ccs	0	31883	id:i:2336554	HG:A:a
A	ptg000001l	21104	-	m84100_240128_024355_s2/242419829/ccs	0	28551	id:i:5120251	HG:A:a
A	ptg000001l	21788	+	m84100_240128_024355_s2/266536351/ccs	0	27994	id:i:5577462	HG:A:a
```

The `S` line is the Segment and it acts as a header for the the contig. `LN:i:` is the segment length (in this case, 21,642,937 bp). The `rd:i:` is the read coverage, calculated by the reads coming from the same contig (in this case, read coverage is 82, which is high). The `A` lines provide information about the sequences that make up the contig. Here's what each column means 

- Column 1: should always be A 
- Column 2: contig name 
- Column 3: contig start coordinate of subregion constructed by read 
- Column 4: read strand (+ or -)
- Column 5: read name 
- Column 6: read start coordinate of subregion which is used to construct contig
- Column 7: read end coordinate of subregion which is used to construct contig
- Column 8: read ID
- Column 9: haplotype status of read. `HG:A:a`, `HG:A:p`, and `HG:A:m` indicate that the read is non-binnable (ie heterozygous), father/hap1 specific, or mother/hap2 specific. 

If I'm interpreting this correctly, it looks like most of the reads in the first contig are heterozygous. 

The error file (`slurm-304463.error`) for this script contains the histogram of the kmers. It has 4 iterations of kmer histograms with some sort of analysis in between the histograms. Here's what the last histogram looks like: 

```
[M::ha_analyze_count] lowest: count[9] = 2315
[M::ha_analyze_count] highest: count[83] = 417609
[M::ha_hist_line]     1: ****************************************************************************************************> 1732549
[M::ha_hist_line]     2: ************* 53793
[M::ha_hist_line]     3: **** 16969
[M::ha_hist_line]     4: ** 9184
[M::ha_hist_line]     5: * 5507
[M::ha_hist_line]     6: * 4361
[M::ha_hist_line]     7: * 3389
[M::ha_hist_line]     8: * 2959
[M::ha_hist_line]     9: * 2315
[M::ha_hist_line]    10: * 2369
[M::ha_hist_line]    11:  2037
[M::ha_hist_line]    12:  1889
[M::ha_hist_line]    13:  1724
[M::ha_hist_line]    14:  1729
[M::ha_hist_line]    15:  1721
[M::ha_hist_line]    16:  1641
[M::ha_hist_line]    17:  1530
[M::ha_hist_line]    18:  1687
[M::ha_hist_line]    19:  1389
[M::ha_hist_line]    20:  1341
[M::ha_hist_line]    21:  1370
[M::ha_hist_line]    22:  1269
[M::ha_hist_line]    23:  1249
[M::ha_hist_line]    24:  1329
[M::ha_hist_line]    25:  1327
[M::ha_hist_line]    26:  1317
[M::ha_hist_line]    27:  1274
[M::ha_hist_line]    28:  1370
[M::ha_hist_line]    29:  1495
[M::ha_hist_line]    30:  1468
[M::ha_hist_line]    31:  1677
[M::ha_hist_line]    32:  1625
[M::ha_hist_line]    33:  1729
[M::ha_hist_line]    34:  1697
[M::ha_hist_line]    35:  1825
[M::ha_hist_line]    36:  1919
[M::ha_hist_line]    37:  1966
[M::ha_hist_line]    38:  2066
[M::ha_hist_line]    39: * 2101
[M::ha_hist_line]    40: * 2195
[M::ha_hist_line]    41: * 2119
[M::ha_hist_line]    42: * 2100
[M::ha_hist_line]    43: * 2325
[M::ha_hist_line]    44: * 2644
[M::ha_hist_line]    45: * 2807
[M::ha_hist_line]    46: * 3080
[M::ha_hist_line]    47: * 3289
[M::ha_hist_line]    48: * 3661
[M::ha_hist_line]    49: * 3984
[M::ha_hist_line]    50: * 4856
[M::ha_hist_line]    51: * 5391
[M::ha_hist_line]    52: ** 6627
[M::ha_hist_line]    53: ** 7648
[M::ha_hist_line]    54: ** 9319
[M::ha_hist_line]    55: *** 11051
[M::ha_hist_line]    56: *** 13316
[M::ha_hist_line]    57: **** 16452
[M::ha_hist_line]    58: ***** 19650
[M::ha_hist_line]    59: ****** 24229
[M::ha_hist_line]    60: ******* 29998
[M::ha_hist_line]    61: ********* 37438
[M::ha_hist_line]    62: *********** 45813
[M::ha_hist_line]    63: ************* 54367
[M::ha_hist_line]    64: **************** 67165
[M::ha_hist_line]    65: ******************* 79086
[M::ha_hist_line]    66: *********************** 95901
[M::ha_hist_line]    67: *************************** 111990
[M::ha_hist_line]    68: ******************************* 128413
[M::ha_hist_line]    69: *********************************** 147679
[M::ha_hist_line]    70: ***************************************** 171224
[M::ha_hist_line]    71: ********************************************** 193638
[M::ha_hist_line]    72: **************************************************** 217984
[M::ha_hist_line]    73: ********************************************************** 242258
[M::ha_hist_line]    74: **************************************************************** 266643
[M::ha_hist_line]    75: ********************************************************************** 291355
[M::ha_hist_line]    76: **************************************************************************** 317522
[M::ha_hist_line]    77: ********************************************************************************* 337960
[M::ha_hist_line]    78: ************************************************************************************** 358656
[M::ha_hist_line]    79: ******************************************************************************************* 378437
[M::ha_hist_line]    80: ********************************************************************************************** 393447
[M::ha_hist_line]    81: ************************************************************************************************* 405399
[M::ha_hist_line]    82: *************************************************************************************************** 413774
[M::ha_hist_line]    83: **************************************************************************************************** 417609
[M::ha_hist_line]    84: **************************************************************************************************** 416365
[M::ha_hist_line]    85: **************************************************************************************************** 417459
[M::ha_hist_line]    86: *************************************************************************************************** 413637
[M::ha_hist_line]    87: ************************************************************************************************* 404341
[M::ha_hist_line]    88: ********************************************************************************************* 387519
[M::ha_hist_line]    89: ***************************************************************************************** 372331
[M::ha_hist_line]    90: ************************************************************************************ 352702
[M::ha_hist_line]    91: ******************************************************************************** 333308
[M::ha_hist_line]    92: ************************************************************************* 305452
[M::ha_hist_line]    93: ******************************************************************* 279706
[M::ha_hist_line]    94: ************************************************************** 257317
[M::ha_hist_line]    95: ******************************************************* 230115
[M::ha_hist_line]    96: ************************************************* 205580
[M::ha_hist_line]    97: ******************************************* 181564
[M::ha_hist_line]    98: ************************************** 159113
[M::ha_hist_line]    99: ********************************* 139005
[M::ha_hist_line]   100: ***************************** 120518
[M::ha_hist_line]   101: ************************* 102686
[M::ha_hist_line]   102: ********************* 86025
[M::ha_hist_line]   103: ****************** 73567
[M::ha_hist_line]   104: *************** 61207
[M::ha_hist_line]   105: ************ 50380
[M::ha_hist_line]   106: ********** 41491
[M::ha_hist_line]   107: ******** 34384
[M::ha_hist_line]   108: ******* 28223
[M::ha_hist_line]   109: ***** 22483
[M::ha_hist_line]   110: **** 18607
[M::ha_hist_line]   111: **** 14975
[M::ha_hist_line]   112: *** 12513
[M::ha_hist_line]   113: ** 10316
[M::ha_hist_line]   114: ** 8237
[M::ha_hist_line]   115: ** 6969
[M::ha_hist_line]   116: * 6015
[M::ha_hist_line]   117: * 5348
[M::ha_hist_line]   118: * 4850
[M::ha_hist_line]   119: * 4508
[M::ha_hist_line]   120: * 4436
[M::ha_hist_line]   121: * 4296
[M::ha_hist_line]   122: * 4655
[M::ha_hist_line]   123: * 4414
[M::ha_hist_line]   124: * 4850
[M::ha_hist_line]   125: * 5053
[M::ha_hist_line]   126: * 5326
[M::ha_hist_line]   127: * 6256
[M::ha_hist_line]   128: ** 6763
[M::ha_hist_line]   129: ** 7359
[M::ha_hist_line]   130: ** 8371
[M::ha_hist_line]   131: ** 9116
[M::ha_hist_line]   132: ** 10114
[M::ha_hist_line]   133: *** 11557
[M::ha_hist_line]   134: *** 12951
[M::ha_hist_line]   135: *** 14573
[M::ha_hist_line]   136: **** 16195
[M::ha_hist_line]   137: **** 17982
[M::ha_hist_line]   138: ***** 19859
[M::ha_hist_line]   139: ***** 22041
[M::ha_hist_line]   140: ****** 24033
[M::ha_hist_line]   141: ****** 26500
[M::ha_hist_line]   142: ******* 30035
[M::ha_hist_line]   143: ******** 32677
[M::ha_hist_line]   144: ********* 36297
[M::ha_hist_line]   145: ********* 39324
[M::ha_hist_line]   146: ********** 43146
[M::ha_hist_line]   147: *********** 47105
[M::ha_hist_line]   148: ************ 52168
[M::ha_hist_line]   149: ************* 55935
[M::ha_hist_line]   150: *************** 61001
[M::ha_hist_line]   151: **************** 65990
[M::ha_hist_line]   152: ***************** 70743
[M::ha_hist_line]   153: ****************** 74363
[M::ha_hist_line]   154: ******************* 79804
[M::ha_hist_line]   155: ******************** 83768
[M::ha_hist_line]   156: ********************* 88057
[M::ha_hist_line]   157: ********************** 92818
[M::ha_hist_line]   158: *********************** 97623
[M::ha_hist_line]   159: ************************* 103918
[M::ha_hist_line]   160: ************************** 107072
[M::ha_hist_line]   161: ************************** 110105
[M::ha_hist_line]   162: *************************** 113902
[M::ha_hist_line]   163: **************************** 117243
[M::ha_hist_line]   164: **************************** 118933
[M::ha_hist_line]   165: ***************************** 122058
[M::ha_hist_line]   166: ****************************** 123371
[M::ha_hist_line]   167: ****************************** 125091
[M::ha_hist_line]   168: ****************************** 125263
[M::ha_hist_line]   169: ****************************** 125254
[M::ha_hist_line]   170: ****************************** 123856
[M::ha_hist_line]   171: ****************************** 123656
[M::ha_hist_line]   172: ***************************** 121423
[M::ha_hist_line]   173: ***************************** 121400
[M::ha_hist_line]   174: **************************** 117980
[M::ha_hist_line]   175: **************************** 115344
[M::ha_hist_line]   176: *************************** 112655
[M::ha_hist_line]   177: ************************** 109202
[M::ha_hist_line]   178: ************************* 105297
[M::ha_hist_line]   179: ************************ 102172
[M::ha_hist_line]   180: *********************** 97507
[M::ha_hist_line]   181: ********************** 93418
[M::ha_hist_line]   182: ********************* 88400
[M::ha_hist_line]   183: ******************** 83674
[M::ha_hist_line]   184: ******************* 77971
[M::ha_hist_line]   185: ***************** 72480
[M::ha_hist_line]   186: **************** 68366
[M::ha_hist_line]   187: *************** 63165
[M::ha_hist_line]   188: ************** 58702
[M::ha_hist_line]   189: ************* 54012
[M::ha_hist_line]   190: ************ 50360
[M::ha_hist_line]   191: *********** 45887
[M::ha_hist_line]   192: ********** 40846
[M::ha_hist_line]   193: ********* 36887
[M::ha_hist_line]   194: ******** 33506
[M::ha_hist_line]   195: ******* 30266
[M::ha_hist_line]   196: ******* 27487
[M::ha_hist_line]   197: ****** 24333
[M::ha_hist_line]   198: ***** 21602
[M::ha_hist_line]   199: ***** 19303
[M::ha_hist_line]   200: **** 17108
[M::ha_hist_line]   201: **** 15156
[M::ha_hist_line]   202: *** 13661
[M::ha_hist_line]   203: *** 12076
[M::ha_hist_line]   204: *** 10526
[M::ha_hist_line]   205: ** 9215
[M::ha_hist_line]   206: ** 8143
[M::ha_hist_line]   207: ** 7395
[M::ha_hist_line]   208: ** 6602
[M::ha_hist_line]   209: * 5949
[M::ha_hist_line]   210: * 5447
[M::ha_hist_line]   211: * 4869
[M::ha_hist_line]   212: * 4270
[M::ha_hist_line]   213: * 3890
[M::ha_hist_line]   214: * 3731
[M::ha_hist_line]   215: * 3629
[M::ha_hist_line]   216: * 3494
[M::ha_hist_line]   217: * 3613
[M::ha_hist_line]   218: * 3512
[M::ha_hist_line]   219: * 3618
[M::ha_hist_line]   220: * 3772
[M::ha_hist_line]   221: * 3774
[M::ha_hist_line]   222: * 3708
[M::ha_hist_line]   223: * 3818
[M::ha_hist_line]   224: * 3986
[M::ha_hist_line]   225: * 4029
[M::ha_hist_line]   226: * 4380
[M::ha_hist_line]   227: * 4386
[M::ha_hist_line]   228: * 4510
[M::ha_hist_line]   229: * 4678
[M::ha_hist_line]   230: * 4797
[M::ha_hist_line]   231: * 5106
[M::ha_hist_line]   232: * 5197
[M::ha_hist_line]   233: * 5242
[M::ha_hist_line]   234: * 5474
[M::ha_hist_line]   235: * 5733
[M::ha_hist_line]   236: * 6021
[M::ha_hist_line]   237: * 6265
[M::ha_hist_line]   238: * 6246
[M::ha_hist_line]   239: ** 6646
[M::ha_hist_line]   240: ** 6722
[M::ha_hist_line]   241: ** 6844
[M::ha_hist_line]   242: ** 6733
[M::ha_hist_line]   243: ** 7254
[M::ha_hist_line]   244: ** 7250
[M::ha_hist_line]   245: ** 7243
[M::ha_hist_line]   246: ** 7275
[M::ha_hist_line]   247: ** 7405
[M::ha_hist_line]   248: ** 7421
[M::ha_hist_line]   249: ** 7571
[M::ha_hist_line]   250: ** 7291
[M::ha_hist_line]   251: ** 7331
[M::ha_hist_line]   252: ** 7309
[M::ha_hist_line]   253: ** 7278
[M::ha_hist_line]   254: ** 7264
[M::ha_hist_line]   255: ** 7092
[M::ha_hist_line]   256: ** 6912
[M::ha_hist_line]   257: ** 6958
[M::ha_hist_line]   258: ** 6689
[M::ha_hist_line]   259: ** 6607
[M::ha_hist_line]   260: ** 6542
[M::ha_hist_line]   261: * 6242
[M::ha_hist_line]   262: * 6185
[M::ha_hist_line]   263: * 5934
[M::ha_hist_line]   264: * 5717
[M::ha_hist_line]   265: * 5388
[M::ha_hist_line]   266: * 5448
[M::ha_hist_line]   267: * 5279
[M::ha_hist_line]   268: * 4944
[M::ha_hist_line]   269: * 4724
[M::ha_hist_line]   270: * 4549
[M::ha_hist_line]   271: * 4404
[M::ha_hist_line]   272: * 4417
[M::ha_hist_line]   273: * 4041
[M::ha_hist_line]   274: * 3888
[M::ha_hist_line]   275: * 3760
[M::ha_hist_line]   276: * 3592
[M::ha_hist_line]   277: * 3490
[M::ha_hist_line]   278: * 3150
[M::ha_hist_line]   279: * 3001
[M::ha_hist_line]   280: * 2926
[M::ha_hist_line]   281: * 3084
[M::ha_hist_line]   282: * 2758
[M::ha_hist_line]   283: * 2672
[M::ha_hist_line]   284: * 2526
[M::ha_hist_line]   285: * 2432
[M::ha_hist_line]   286: * 2313
[M::ha_hist_line]   287: * 2314
[M::ha_hist_line]   288: * 2255
[M::ha_hist_line]   289: * 2266
[M::ha_hist_line]   290: * 2230
[M::ha_hist_line]   291: * 2124
[M::ha_hist_line]   292: * 2109
[M::ha_hist_line]   293: * 2208
[M::ha_hist_line]   294: * 2187
[M::ha_hist_line]   295: * 2154
[M::ha_hist_line]   296: * 2139
[M::ha_hist_line]   297: * 2160
[M::ha_hist_line]   298: * 2211
[M::ha_hist_line]   299: * 2220
[M::ha_hist_line]   300: * 2364
[M::ha_hist_line]   301: * 2340
[M::ha_hist_line]   302: * 2367
[M::ha_hist_line]   303: * 2554
[M::ha_hist_line]   304: * 2611
[M::ha_hist_line]   305: * 2559
[M::ha_hist_line]   306: * 2525
[M::ha_hist_line]   307: * 2573
[M::ha_hist_line]   308: * 2791
[M::ha_hist_line]   309: * 2797
[M::ha_hist_line]   310: * 2819
[M::ha_hist_line]   311: * 2945
[M::ha_hist_line]   312: * 2958
[M::ha_hist_line]   313: * 3034
[M::ha_hist_line]   314: * 3275
[M::ha_hist_line]   315: * 3267
[M::ha_hist_line]   316: * 3351
[M::ha_hist_line]   317: * 3284
[M::ha_hist_line]   318: * 3513
[M::ha_hist_line]   319: * 3510
[M::ha_hist_line]   320: * 3692
[M::ha_hist_line]   321: * 3704
[M::ha_hist_line]   322: * 3755
[M::ha_hist_line]   323: * 3906
[M::ha_hist_line]   324: * 3910
[M::ha_hist_line]   325: * 3888
[M::ha_hist_line]   326: * 4038
[M::ha_hist_line]   327: * 4052
[M::ha_hist_line]   328: * 4249
[M::ha_hist_line]   329: * 4048
[M::ha_hist_line]   330: * 3908
[M::ha_hist_line]   331: * 4098
[M::ha_hist_line]   332: * 3993
[M::ha_hist_line]   333: * 4066
[M::ha_hist_line]   334: * 4055
[M::ha_hist_line]   335: * 4079
[M::ha_hist_line]   336: * 4027
[M::ha_hist_line]   337: * 3871
[M::ha_hist_line]   338: * 3942
[M::ha_hist_line]   339: * 3919
[M::ha_hist_line]   340: * 3896
[M::ha_hist_line]   341: * 3891
[M::ha_hist_line]   342: * 3721
[M::ha_hist_line]   343: * 3773
[M::ha_hist_line]   344: * 3657
[M::ha_hist_line]   345: * 3596
[M::ha_hist_line]   346: * 3377
[M::ha_hist_line]   347: * 3273
[M::ha_hist_line]   348: * 3190
[M::ha_hist_line]   349: * 3224
[M::ha_hist_line]   350: * 3142
[M::ha_hist_line]   351: * 3076
[M::ha_hist_line]   352: * 3058
[M::ha_hist_line]   353: * 2916
[M::ha_hist_line]   354: * 2776
[M::ha_hist_line]   355: * 2764
[M::ha_hist_line]   356: * 2785
[M::ha_hist_line]   357: * 2701
[M::ha_hist_line]   358: * 2490
[M::ha_hist_line]   359: * 2416
[M::ha_hist_line]   360: * 2361
[M::ha_hist_line]   361: * 2298
[M::ha_hist_line]   362: * 2325
[M::ha_hist_line]   363: * 2146
[M::ha_hist_line]   364: * 2136
[M::ha_hist_line]   365: * 2099
[M::ha_hist_line]  rest: *********************************************************************************************** 396842
[M::ha_analyze_count] left: none
[M::ha_analyze_count] right: none
[M::ha_pt_gen] peak_hom: 83; peak_het: -1
[M::ha_ct_shrink::285039.075*35.33] ==> counted 17036513 distinct minimizer k-mers
[M::ha_pt_gen::] counting in normal mode
[M::yak_count] collected 2137320573 minimizers
[M::ha_pt_gen::285482.798*35.31] ==> indexed 2135588024 positions, counted 17036513 distinct minimizer k-mers
[M::ha_assemble::297514.365*35.34@246.283GB] ==> found overlaps for the final round
[M::ha_print_ovlp_stat] # overlaps: 1183659490
[M::ha_print_ovlp_stat] # strong overlaps: 596745657
[M::ha_print_ovlp_stat] # weak overlaps: 586913833
[M::ha_print_ovlp_stat] # exact overlaps: 1149035029
[M::ha_print_ovlp_stat] # inexact overlaps: 34624461
[M::ha_print_ovlp_stat] # overlaps without large indels: 1180991551
[M::ha_print_ovlp_stat] # reverse overlaps: 410771757
Writing reads to disk... 
Reads has been written.
```

Thats a lot of information. The Hifiasm [output](https://hifiasm.readthedocs.io/en/latest/interpreting-output.html#homcov) page says that for heterozygous samples (which mine likely are), there should be 2 peaks in the k-mer plot, where the smaller peak is around the heterozygous read coverage and the larger peak is around the homozygous read coverage. This is true for all k-mer plots produced from this data. In all of my k-mer plots, the homozygous peak is 83 and the heterozygous peak is 168. I'm going to include the information without the k-mer plots below because they are so large: 

```
[M::ha_analyze_count] lowest: count[17] = 11309
[M::ha_analyze_count] highest: count[85] = 9499858

## first k-mer plot

[M::ha_analyze_count] left: none
[M::ha_analyze_count] right: count[168] = 3093394
[M::ha_ft_gen] peak_hom: 168; peak_het: 85
[M::ha_ct_shrink::3427.856*4.32] ==> counted 2684382 distinct minimizer k-mers
[M::ha_ft_gen::3431.212*4.31@20.882GB] ==> filtered out 2684382 k-mers occurring 840 or more times
[M::ha_opt_update_cov] updated max_n_chain to 840
[M::yak_count] collected 2139678804 minimizers
[M::ha_pt_gen::4355.642*5.80] ==> counted 31723566 distinct minimizer k-mers
[M::ha_pt_gen] count[4095] = 0 (for sanity check)
[M::ha_analyze_count] lowest: count[17] = 2230
[M::ha_analyze_count] highest: count[83] = 418590

## second k-mer plot

[M::ha_analyze_count] left: none
[M::ha_analyze_count] right: count[168] = 125423
[M::ha_pt_gen] peak_hom: 168; peak_het: 83
[M::ha_ct_shrink::4355.907*5.81] ==> counted 17720475 distinct minimizer k-mers
[M::ha_pt_gen::] counting in normal mode
[M::yak_count] collected 2139678804 minimizers
[M::ha_pt_gen::4818.141*7.37] ==> indexed 2125675713 positions, counted 17720475 distinct minimizer k-mers
[M::ha_assemble::100006.981*34.52@107.237GB] ==> corrected reads for round 1
[M::ha_assemble] # bases: 79183709778; # corrected bases: 91988177; # recorrected bases: 115414
[M::ha_assemble] size of buffer: 61.709GB
[M::yak_count] collected 2137476150 minimizers
[M::ha_pt_gen::100401.244*34.48] ==> counted 19037497 distinct minimizer k-mers
[M::ha_pt_gen] count[4095] = 0 (for sanity check)
[M::ha_analyze_count] lowest: count[13] = 1719
[M::ha_analyze_count] highest: count[83] = 417869

## third k-mer plot

[M::ha_analyze_count] left: none
[M::ha_analyze_count] right: none
[M::ha_pt_gen] peak_hom: 83; peak_het: -1
[M::ha_ct_shrink::100401.532*34.48] ==> counted 17060420 distinct minimizer k-mers
[M::ha_pt_gen::] counting in normal mode
[M::yak_count] collected 2137476150 minimizers
[M::ha_pt_gen::100851.319*34.43] ==> indexed 2135499073 positions, counted 17060420 distinct minimizer k-mers
[M::ha_assemble::192251.063*35.13@143.827GB] ==> corrected reads for round 2
[M::ha_assemble] # bases: 79186546424; # corrected bases: 1694280; # recorrected bases: 16078
[M::ha_assemble] size of buffer: 60.552GB
[M::yak_count] collected 2137334800 minimizers
[M::ha_pt_gen::192627.741*35.11] ==> counted 18787923 distinct minimizer k-mers
[M::ha_pt_gen] count[4095] = 0 (for sanity check)
[M::ha_analyze_count] lowest: count[9] = 2352
[M::ha_analyze_count] highest: count[83] = 417664

## fourth k-mer plot

[M::ha_analyze_count] left: none
[M::ha_analyze_count] right: none
[M::ha_pt_gen] peak_hom: 83; peak_het: -1
[M::ha_ct_shrink::192628.074*35.11] ==> counted 17040202 distinct minimizer k-mers
[M::ha_pt_gen::] counting in normal mode
[M::yak_count] collected 2137334800 minimizers
[M::ha_pt_gen::193066.706*35.08] ==> indexed 2135587079 positions, counted 17040202 distinct minimizer k-mers
[M::ha_assemble::284661.787*35.35@241.898GB] ==> corrected reads for round 3
[M::ha_assemble] # bases: 79186730278; # corrected bases: 235080; # recorrected bases: 12650
[M::ha_assemble] size of buffer: 60.546GB
[M::yak_count] collected 2137320573 minimizers
[M::ha_pt_gen::285038.830*35.33] ==> counted 18769062 distinct minimizer k-mers
[M::ha_pt_gen] count[4095] = 0 (for sanity check)
[M::ha_analyze_count] lowest: count[9] = 2315
[M::ha_analyze_count] highest: count[83] = 417609

## fifth k-mer plot

[M::ha_analyze_count] left: none
[M::ha_analyze_count] right: none
[M::ha_pt_gen] peak_hom: 83; peak_het: -1
[M::ha_ct_shrink::285039.075*35.33] ==> counted 17036513 distinct minimizer k-mers
[M::ha_pt_gen::] counting in normal mode
[M::yak_count] collected 2137320573 minimizers
[M::ha_pt_gen::285482.798*35.31] ==> indexed 2135588024 positions, counted 17036513 distinct minimizer k-mers
[M::ha_assemble::297514.365*35.34@246.283GB] ==> found overlaps for the final round
[M::ha_print_ovlp_stat] # overlaps: 1183659490
[M::ha_print_ovlp_stat] # strong overlaps: 596745657
[M::ha_print_ovlp_stat] # weak overlaps: 586913833
[M::ha_print_ovlp_stat] # exact overlaps: 1149035029
[M::ha_print_ovlp_stat] # inexact overlaps: 34624461
[M::ha_print_ovlp_stat] # overlaps without large indels: 1180991551
[M::ha_print_ovlp_stat] # reverse overlaps: 410771757
Writing reads to disk... 
Reads has been written.
Writing ma_hit_ts to disk... 
ma_hit_ts has been written.
Writing ma_hit_ts to disk... 
ma_hit_ts has been written.
bin files have been written.
[M::purge_dups] homozygous read coverage threshold: 168
[M::purge_dups] purge duplication coverage threshold: 210
Writing raw unitig GFA to disk... 
Writing processed unitig GFA to disk... 
[M::purge_dups] homozygous read coverage threshold: 168
[M::purge_dups] purge duplication coverage threshold: 210
[M::mc_solve_core::0.284] ==> Partition
[M::adjust_utg_by_primary] primary contig coverage range: [142, infinity]
Writing apul.hifiasm.bp.p_ctg.gfa to disk... 
[M::adjust_utg_by_trio] primary contig coverage range: [142, infinity]
Writing apul.hifiasm.bp.hap1.p_ctg.gfa to disk... 
[M::adjust_utg_by_trio] primary contig coverage range: [142, infinity]
Writing apul.hifiasm.bp.hap2.p_ctg.gfa to disk... 
Inconsistency threshold for low-quality regions in BED files: 70%
[M::main] Version: 0.16.1-r375
[M::main] CMD: hifiasm -o apul.hifiasm -t 36 m84100_240128_024355_s2.hifi_reads.bc1029.fastq.fastq
[M::main] Real time: 304623.035 sec; CPU: 10520661.634 sec; Peak RSS: 246.283 GB
```

Okay so 4 rounds of assembly were done, hypothetically improving the assembly each time. Weirdly, the first 2 rounds had the heterozygous peak at 68, but the last couple of rounds had the heterozygous peak at -1. In all k-mer plots, there is a high peak at the 1 location, but this doesn't seem unusual (based on the example posts [here](https://github.com/chhylp123/hifiasm/issues/49#issue-729106823) and [here](https://github.com/chhylp123/hifiasm/issues/10#issuecomment-616213684) provided by the hifiasm log interpretation section). Overall, I think this file is providing information on the assembly iterations. 

Given all of this information, let's now run BUSCO to assess completeness of assembly. In the scripts folder: `nano busco_unfilt_hifiasm.sh`

```
#!/bin/bash 
#SBATCH -t 100:00:00
#SBATCH --nodes=1 --ntasks-per-node=15
#SBATCH --export=NONE
#SBATCH --mem=250GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Apul_Genome/assembly/scripts
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

echo "Convert from gfa to fasta for downstream use" $(date)

cd /data/putnamlab/jillashey/Apul_Genome/assembly/data/

awk '/^S/{print ">"$2;print $3}' apul.hifiasm.bp.p_ctg.gfa > apul.hifiasm.bp.p_ctg.fa

echo "Begin busco on unfiltered hifiasm-assembled fasta" $(date)

labbase=/data/putnamlab
busco_shared="${labbase}/shared/busco"
[ -z "$query" ] && query="${labbase}/jillashey/Apul_Genome/assembly/data/apul.hifiasm.bp.p_ctg.fa" # set this to the query (genome/transcriptome) you are running
[ -z "$db_to_compare" ] && db_to_compare="${busco_shared}/downloads/lineages/metazoa_odb10"

source "${busco_shared}/scripts/busco_init.sh"  # sets up the modules required for this in the right order

# This will generate output under your $HOME/busco_output
cd "${labbase}/${Apul_Genome/assembly/data}"
busco --config "$EBROOTBUSCO/config/config.ini"  -f -c 20 --long -i "${query}" -l metazoa_odb10 -o apul.busco.canu -m genome

echo "busco complete for unfiltered hifiasm-assembled fasta" $(date)
```

Submitted batch job 305426. I could also run some analysis on the hap assemblies, but I'm going to wait until I am done with the filtering and re-assembly, as that assembly is the one that I will likely be moving forward with. How are they making the distinction between hap1/father and hap2/mother? 

Busco ran in ~45 mins and results look a lot better than they did with canu! This is the primary result from the output file: 

```
2024-03-11 13:37:21 INFO:       

        --------------------------------------------------
        |Results from dataset metazoa_odb10               |
        --------------------------------------------------
        |C:93.3%[S:92.0%,D:1.3%],F:3.1%,M:3.6%,n:954      |
        |890    Complete BUSCOs (C)                       |
        |878    Complete and single-copy BUSCOs (S)       |
        |12     Complete and duplicated BUSCOs (D)        |
        |30     Fragmented BUSCOs (F)                     |
        |34     Missing BUSCOs (M)                        |
        |954    Total BUSCO groups searched               |
        --------------------------------------------------
```

This looks so much better than the canu assembly! The previous canu assembly was 94.4% complete, but had only 9.4% single copy BUSCOs and 85% duplicated BUSCOs. Ideally, the duplication level should be lower. The hifiasm assembly had 93.3% completeness, 92% single copy BUSCOs, and 1.3% duplicated BUSCOs. Hifiasm is definitely the way to go for assembly. Now I just have to wait until the blast prok is done so I can remove any contamination. After I remove contamination, I will re-assemble using hifiasm. 

### 20240313

Blast prok failed after 2 days and then restarted on Andromeda. Maybe I should stop the code after 2 days...idk. Maybe I need to increase the memory for the job? Canceling the job (`305351`) and increasing the memory (`#SBATCH --mem=500GB`). Submitted batch job 308997

### 20240316

Copied the prok data to my local computer. It is still running on the server, but I'm nervous it will restart again. If it does restart, I'll cancel the job and just use this data from today. On my local computer, I combined the prok and viral blast results and then removed any hits whose bit score was <1000. 

```
cd /Users/jillashey/Desktop/PutnamLab/Apulchra_genome
cat viral_contaminant_hits_rr.txt prok_contaminant_hits_rr.txt > all_contaminant_hits_rr.txt
awk '$12 > 1000 {print $0}' all_contaminant_hits_rr.txt > contaminant_hits_pv_passfilter_rr.txt
```

I then looked at the contamination hits in R. See code [here](https://github.com/hputnam/Apulchra_genome/blob/main/scripts/genome_analysis.Rmd). 

As a summary, I first read in the eukaryotic blast hits that passed the contamination threshold. I found that only 2 reads had any euk contamination (`m84100_240128_024355_s2/48759857/ccs` and `m84100_240128_024355_s2/234751852/ccs`). I then read in the prokaryotic and viral blast hits that passed the threshold (only 224 blast hits passed). I calculated the percentage of each hits align length to the contigs so if there was a result that had 100%, that would mean that the whole contig was a contaminant. I looked at a histogram of the % alignments and found that most of the % alignments are on the lower size and there are not many 100% sequences. I summarized the contigs that were to be filtered out and found that 222 contigs had some level of pv contamination. I added the euk + pv contamination reads together (224 total) and calculated the proportion of contamination to raw reads. The contamination ended up being only 0.003797649% of the raw reads, which is pretty amazing! I calculated the mean length of the filtered reads (13,242.1 bp) and the sums of the unfiltered read length (79183709778 total bp) and filtered read length (79181142809 total bp).  Using these lengths, I calculated a rough estimation of sequencing depth and found we have roughly 100x coverage! That's similar to the Young et al. 2024 results as well. Finally, I wrote the list of filtered reads to a text file on my local computer. This information will be used to filter the raw reads on Andromeda prior to assembly. I'm impressed with the low contamination and the high coverage of the PacBio HiFi reads. 

### 20240319 

Prok blast script finally finished running! Took about 5 days. Cat the prok and viral results together and remove anything that has a bit score <1000. 

```
cd /data/putnamlab/jillashey/Apul_Genome/assembly/data

cat viral_contaminant_hits_rr.txt prok_contaminant_hits_rr.txt > all_contaminant_hits_rr.txt
awk '$12 > 1000 {print $0}' all_contaminant_hits_rr.txt > contaminant_hits_pv_passfilter_rr.txt
```

Similarly to what I did above, I looked at the contamination hits in R. See code [here](https://github.com/hputnam/Apulchra_genome/blob/main/scripts/genome_analysis.Rmd). 

As a summary, I first read in the eukaryotic blast hits that passed the contamination threshold. I found that only 2 reads had any euk contamination (`m84100_240128_024355_s2/48759857/ccs` and `m84100_240128_024355_s2/234751852/ccs`). I then read in the prokaryotic and viral blast hits that passed the threshold (only 2 blast hits passed). I calculated the percentage of each hits align length to the contigs so if there was a result that had 100%, that would mean that the whole contig was a contaminant. I looked at a histogram of the % alignments and found that most of the % alignments are on the lower size and there are not many 100% sequences. I summarized the contigs that were to be filtered out and found that 494 contigs had some level of pv contamination. I added the euk + pv contamination reads together (496 contigs total) and calculated the proportion of contamination to raw reads. The contamination ended up being only 0.00840908% of the raw reads, which is pretty amazing! I calculated the mean length of the filtered reads (13,424.84 bp) and the sums of the unfiltered read length (79183709778 total bp) and filtered read length (79178244314 total bp).  Using these lengths, I calculated a rough estimation of sequencing depth and found we have roughly 100x coverage! That's similar to the Young et al. 2024 results as well. Finally, I wrote the list of filtered reads to a text file on my local computer. This information will be used to filter the raw reads on Andromeda prior to assembly. I'm impressed with the low contamination and the high coverage of the PacBio HiFi reads. 

### 20240320

Before filtered the hifi reads, I'm going to clean up the `/data/putnamlab/jillashey/Apul_Genome/assembly/data` folder so that I have more memory for the next steps. Here's whats in there right now: 

```
total 312G
-rw-r--r--. 1 jillashey 148G Feb  8 02:37 m84100_240128_024355_s2.hifi_reads.bc1029.fastq.fastq
-rwxr-xr-x. 1 jillashey 1.1K Feb  8 14:13 apul.seqStore.sh
-rw-r--r--. 1 jillashey  951 Feb  8 14:31 apul.seqStore.err
drwxr-xr-x. 3 jillashey 4.0K Feb  8 14:32 apul.seqStore
-rw-r--r--. 1 jillashey  23K Feb  9 01:47 apul.report
-rw-r--r--. 1 jillashey 7.0M Feb  9 01:51 apul.contigs.layout.tigInfo
-rw-r--r--. 1 jillashey 155M Feb  9 01:51 apul.contigs.layout.readToTig
-rw-r--r--. 1 jillashey 2.8G Feb  9 01:57 apul.unassembled.fasta
-rw-r--r--. 1 jillashey 943M Feb  9 02:01 apul.contigs.fasta
drwxr-xr-x. 2 jillashey 4.0K Feb 13 14:01 busco_output
drwxr-xr-x. 3 jillashey 4.0K Feb 13 14:01 busco_downloads
-rw-r--r--. 1 jillashey 7.2K Feb 13 14:02 busco_96228.log
lrwxrwxrwx. 1 jillashey  108 Feb 20 14:08 m84100_240128_024355_s2.hifi_reads.bc1029.bam -> /data/putnamlab/KITT/hputnam/20240129_Apulchra_Genome_LongRead/m84100_240128_024355_s2.hifi_reads.bc1029.bam
-rw-r--r--. 1 jillashey  106 Feb 21 16:48 pb.fofn
drwxr-xr-x. 9 jillashey 4.0K Feb 21 16:53 unitigging
-rw-r--r--. 1 jillashey  74G Mar  1 16:08 m84100_240128_024355_s2.hifi_reads.bc1029.fasta
-rw-r--r--. 1 jillashey 244M Mar  1 16:36 rr_read_lengths.txt
-rw-r--r--. 1 jillashey 106K Mar  3 16:27 contaminant_hits_euks_rr.txt
-rw-r--r--. 1 jillashey 2.4K Mar  3 16:29 contaminants_pass_filter_euks_rr.txt
-rw-r--r--. 1 jillashey  69M Mar  5 23:20 viral_contaminant_hits_rr.txt
-rw-r--r--. 1 jillashey  19G Mar  7 23:58 apul.hifiasm.ec.bin
-rw-r--r--. 1 jillashey  47G Mar  8 00:08 apul.hifiasm.ovlp.source.bin
-rw-r--r--. 1 jillashey  17G Mar  8 00:12 apul.hifiasm.ovlp.reverse.bin
-rw-r--r--. 1 jillashey 1.2G Mar  8 01:37 apul.hifiasm.bp.r_utg.gfa
-rw-r--r--. 1 jillashey  21M Mar  8 01:37 apul.hifiasm.bp.r_utg.noseq.gfa
-rw-r--r--. 1 jillashey 8.6M Mar  8 01:41 apul.hifiasm.bp.r_utg.lowQ.bed
-rw-r--r--. 1 jillashey 1.1G Mar  8 01:42 apul.hifiasm.bp.p_utg.gfa
-rw-r--r--. 1 jillashey  21M Mar  8 01:42 apul.hifiasm.bp.p_utg.noseq.gfa
-rw-r--r--. 1 jillashey 8.2M Mar  8 01:46 apul.hifiasm.bp.p_utg.lowQ.bed
-rw-r--r--. 1 jillashey 506M Mar  8 01:47 apul.hifiasm.bp.p_ctg.gfa
-rw-r--r--. 1 jillashey  11M Mar  8 01:47 apul.hifiasm.bp.p_ctg.noseq.gfa
-rw-r--r--. 1 jillashey 2.0M Mar  8 01:49 apul.hifiasm.bp.p_ctg.lowQ.bed
-rw-r--r--. 1 jillashey 469M Mar  8 01:50 apul.hifiasm.bp.hap1.p_ctg.gfa
-rw-r--r--. 1 jillashey 9.9M Mar  8 01:50 apul.hifiasm.bp.hap1.p_ctg.noseq.gfa
-rw-r--r--. 1 jillashey 2.0M Mar  8 01:52 apul.hifiasm.bp.hap1.p_ctg.lowQ.bed
-rw-r--r--. 1 jillashey 468M Mar  8 01:52 apul.hifiasm.bp.hap2.p_ctg.gfa
-rw-r--r--. 1 jillashey 9.9M Mar  8 01:52 apul.hifiasm.bp.hap2.p_ctg.noseq.gfa
-rw-r--r--. 1 jillashey 1.9M Mar  8 01:54 apul.hifiasm.bp.hap2.p_ctg.lowQ.bed
-rw-r--r--. 1 jillashey 495M Mar 11 12:52 apul.hifiasm.bp.p_ctg.fa
-rw-r--r--. 1 jillashey 246M Mar 19 14:05 prok_contaminant_hits_rr.txt
-rw-r--r--. 1 jillashey 315M Mar 19 16:08 all_contaminant_hits_rr.txt
-rw-r--r--. 1 jillashey 1.1M Mar 19 16:09 contaminant_hits_pv_passfilter_rr.txt
```

I removed the following: 

```
rm -r apul.seqStore
rm apul*
rm -r busco*
rm pb.fofn 
rm -r unitigging/
```

Now I have more space. Copy the file `all_contam_rem_good_hifi_read_list.txt` that was generated from the [R code](https://github.com/hputnam/Apulchra_genome/blob/main/scripts/genome_analysis.Rmd) mentioned above. This specific file was written starting on line 242. It contains the reads that have passed contamination filtering. I copied this file into `/data/putnamlab/jillashey/Apul_Genome/assembly/data`. 

```
wc -l all_contam_rem_good_hifi_read_list.txt
 5897892 all_contam_rem_good_hifi_read_list.txt
```

The vast majority of the hifi reads are retained after contamination filtering, which is a good sign of high quality sequencing. My next step is to subset the raw hifi fasta file to remove the contaminants identified above. I can do this with the [seqtk subseq](https://github.com/lh3/seqtk) command. In the scripts folder: `nano subseq.sh`

```
#!/bin/bash 
#SBATCH -t 100:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --export=NONE
#SBATCH --mem=500GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Apul_Genome/assembly/scripts
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

module load seqtk/1.3-GCC-9.3.0

echo "Subsetting hifi reads that passed contamination filtering" $(date)

cd /data/putnamlab/jillashey/Apul_Genome/assembly/data

seqtk subseq m84100_240128_024355_s2.hifi_reads.bc1029.fasta all_contam_rem_good_hifi_read_list.txt > hifi_rr_allcontam_rem.fasta

echo "Subsetting complete!" $(date)
```

Submitted batch job 309659. Finished in about 8 mins but the output file looks like this: 

```
>m84100_240128_024355_s2/261887593/ccs:18224-18224
A
>m84100_240128_024355_s2/255530003/ccs:21870-21870
A
>m84100_240128_024355_s2/249237028/ccs:23691-23691
A
>m84100_240128_024355_s2/262606536/ccs:14772-14772
A
>m84100_240128_024355_s2/217322854/ccs:12923-12923
A
>m84100_240128_024355_s2/256512826/ccs:12914-12914
A
>m84100_240128_024355_s2/245632166/ccs:28440-28440
A
>m84100_240128_024355_s2/250548903/ccs:23076-23076
A
>m84100_240128_024355_s2/256054930/ccs:15405-15405
A
>m84100_240128_024355_s2/241242930/ccs:15521-15521
A
>m84100_240128_024355_s2/254348773/ccs:14578-14578
C
>m84100_240128_024355_s2/252319399/ccs:12407-12407
A
>m84100_240128_024355_s2/229183717/ccs:5757-5757
```

Not ideal. It looks like it only took the first letter from each sequence. Maybe I need to remove the length information from the `all_contam_rem_good_hifi_read_list.txt` file? 

```
awk '{$2=""; print $0}' all_contam_rem_good_hifi_read_list.txt > output_file.txt
```

Edit the script so that the list of reads to keep is `output_file.txt` and decrease mem to 250GB. Submitted batch job 309672. This appears to have worked! 

```
zgrep -c ">" hifi_rr_allcontam_rem.fasta
5897892
```

Following Young et al. 2024, I need to use [Merqury](https://github.com/marbl/merqury/tree/master) and [GenomeScope2](https://github.com/tbenavi1/genomescope2.0?tab=readme-ov-file) analysis of the cleaned reads. Merqury and Meryl seem to be related somehow, but not sure. [Meryl](https://github.com/marbl/meryl?tab=readme-ov-file) is a tool for counting and working with sets of k-mers. It is a part of Canu as well. So it seems like Meryl counts the k-mers and Merqury estimates accuracy and completeness? Young et al. 2024 used only meryl here to generate a kmer database, which was then used as input to genomescope2. Merqury was used after hifiasm assembly. 

First, the best value for k needs to be determined for use in Meryl. This can be done with the [`best_k.sh`](https://github.com/marbl/merqury/blob/master/best_k.sh) script from Merqury. I'm going to copy this script into my own scripts folder. Young et al. 2024 used 500mb estimated genome size. This is what I will use too, as previous coral/Acropora genomes are about this size. I will need to install Meryl and Merqury first. I'm following the Meryl and Merqury githubs for installation instructions. 

For Meryl: 

```
wget https://github.com/marbl/meryl/releases/download/v1.4.1/meryl-1.4.1.Linux-amd64.tar.xz
tar -xJf meryl-1.4.1.Linux-amd64.tar.xz
export PATH=/data/putnamlab/jillashey/Apul_Genome/assembly/meryl-1.4.1/bin:$PATH
```

Now that I have exported the PATH variable to include the directory where the Meryl executable files are, I can run Meryl commands (I think), regardless of working directory. For now, lets run meryl with k=18. The k-mer database will be generated using k=18 and the cleaned reads. In the scripts folder: `nano meryl.sh`

```
#!/bin/bash 
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=250GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Apul_Genome/assembly/scripts
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

export PATH=/data/putnamlab/jillashey/Apul_Genome/assembly/meryl-1.4.1/bin:$PATH

echo "Creating meryl k-mer db" $(date)

cd /data/putnamlab/jillashey/Apul_Genome/assembly/data

meryl k=18 \
count hifi_rr_allcontam_rem.fasta \
output meryl_merc

echo "Meryl k-mer db complete!" $(date)
```

Submitted batch job 309682. While this runs, lets try to install Merqury. 

```
git clone https://github.com/marbl/merqury.git
cd merqury
export MERQURY=$PWD:$PATH
```

Hypothetically, Merqury is now installed. Perhaps now we can run the best k script? I am not sure what genome size estimate to use. [Shinzato et al. 2020](https://academic.oup.com/mbe/article/38/1/16/5900672) assessed several Acropora genomes and found the genome size ranged from 384 (Acropora microphthalma) to 447 (Acropora hyacinthus) Mb. I think I will go with 450 Mb. Let's try to run the best k script. 

```
$MERQURY/best_k.sh 450000000
```

Giving me this error: 

```
-bash: /data/putnamlab/jillashey/Apul_Genome/assembly/merqury:/data/putnamlab/jillashey/Apul_Genome/assembly/meryl-1.4.1/bin:/path/to/meryl-1.4.1/bin:/path/to/meryl/…/bin:/opt/software/Miniconda3/4.9.2/bin:/opt/software/Miniconda3/4.9.2/condabin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/home/jillashey/.local/bin:/home/jillashey/bin/best_k.sh: No such file or directory
```

I'm not sure what it means or how to fix it. On the Merqury [github](https://github.com/marbl/merqury), it lists dependencies that Merqury may need to run: 

- gcc 10.2.0 or higher (for installing Meryl)
- Meryl v1.4.1
- Java run time environment (JRE)
- R with argparse, ggplot2, and scales (recommend R 4.0.3+)
- bedtools
- samtools

Maybe I need to load these? Idk I am confused. I may just continue with the initial assembly...because I don't really see the importance of this step. In Young et al. 2024, it says "the kmer profile of cleaned raw HiFi reads was generated with Meryl [34], and used for genome profiling with GenomeScope2 [69] to estimate genome size, repetitiveness, heterozygosity, and ploidy." I will come back to this. I may need to email Kevin Bryan to install meryl, merqury and genomescope2. 

While I wait for his response, I will start running the initial hifiasm assembly with default flags. I have a script for hifiasm on the server, but I'm going to make a new one for the initial assembly. In the scripts folder: `nano initial_hifiasm.sh`

```
#!/bin/bash -i
#SBATCH -t 30-00:00:00
#SBATCH --nodes=1 --ntasks-per-node=36
#SBATCH --export=NONE
#SBATCH --mem=500GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH --exclusive
#SBATCH -D /data/putnamlab/jillashey/Apul_Genome/assembly/scripts
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

conda activate /data/putnamlab/conda/hifiasm

cd /data/putnamlab/jillashey/Apul_Genome/assembly/data

echo "Starting assembly with hifiasm" $(date)

hifiasm -o apul.hifiasm.intial hifi_rr_allcontam_rem.fasta -t 36 2> apul_hifiasm_allcontam_rem_initial.asm.log

echo "Assembly with hifiasm complete!" $(date)

conda deactivate
``` 

Submitted batch job 309689

### 20240325

Initial assembly ran in about 3 days. These are the files that were generated: 

```
cd /data/putnamlab/jillashey/Apul_Genome/assembly/data

-rw-r--r--. 1 jillashey  19G Mar 24 02:10 apul.hifiasm.intial.ec.bin
-rw-r--r--. 1 jillashey  47G Mar 24 02:20 apul.hifiasm.intial.ovlp.source.bin
-rw-r--r--. 1 jillashey  17G Mar 24 02:23 apul.hifiasm.intial.ovlp.reverse.bin
-rw-r--r--. 1 jillashey 1.2G Mar 24 03:40 apul.hifiasm.intial.bp.r_utg.gfa
-rw-r--r--. 1 jillashey  21M Mar 24 03:40 apul.hifiasm.intial.bp.r_utg.noseq.gfa
-rw-r--r--. 1 jillashey 8.6M Mar 24 03:44 apul.hifiasm.intial.bp.r_utg.lowQ.bed
-rw-r--r--. 1 jillashey 1.1G Mar 24 03:45 apul.hifiasm.intial.bp.p_utg.gfa
-rw-r--r--. 1 jillashey  21M Mar 24 03:45 apul.hifiasm.intial.bp.p_utg.noseq.gfa
-rw-r--r--. 1 jillashey 8.2M Mar 24 03:49 apul.hifiasm.intial.bp.p_utg.lowQ.bed
-rw-r--r--. 1 jillashey 506M Mar 24 03:50 apul.hifiasm.intial.bp.p_ctg.gfa
-rw-r--r--. 1 jillashey  11M Mar 24 03:50 apul.hifiasm.intial.bp.p_ctg.noseq.gfa
-rw-r--r--. 1 jillashey 2.0M Mar 24 03:52 apul.hifiasm.intial.bp.p_ctg.lowQ.bed
-rw-r--r--. 1 jillashey 469M Mar 24 03:52 apul.hifiasm.intial.bp.hap1.p_ctg.gfa
-rw-r--r--. 1 jillashey 9.9M Mar 24 03:52 apul.hifiasm.intial.bp.hap1.p_ctg.noseq.gfa
-rw-r--r--. 1 jillashey 2.0M Mar 24 03:54 apul.hifiasm.intial.bp.hap1.p_ctg.lowQ.bed
-rw-r--r--. 1 jillashey 468M Mar 24 03:55 apul.hifiasm.intial.bp.hap2.p_ctg.gfa
-rw-r--r--. 1 jillashey 9.9M Mar 24 03:55 apul.hifiasm.intial.bp.hap2.p_ctg.noseq.gfa
-rw-r--r--. 1 jillashey 1.9M Mar 24 03:56 apul.hifiasm.intial.bp.hap2.p_ctg.lowQ.bed
-rw-r--r--. 1 jillashey  80K Mar 24 03:57 apul_hifiasm_allcontam_rem_initial.asm.log
```

Many files. The output file description can be found above and also on the hifiasm [output](https://hifiasm.readthedocs.io/en/latest/interpreting-output.html#interpreting-output) website. The log output file contains the k-mer histogram (similar to what is posted above), which shows two peaks, indicative of a heterozygous genome assembly. This is what the bottom of the log file looks like: 

```
[M::ha_pt_gen::] counting in normal mode
[M::yak_count] collected 2137171007 minimizers
[M::ha_pt_gen::281249.073*35.34] ==> indexed 2135566893 positions, counted 17030417 distinct minimizer k-mers
[M::ha_assemble::292970.825*35.36@202.994GB] ==> found overlaps for the final round
[M::ha_print_ovlp_stat] # overlaps: 1183659340
[M::ha_print_ovlp_stat] # strong overlaps: 596745652
[M::ha_print_ovlp_stat] # weak overlaps: 586913688
[M::ha_print_ovlp_stat] # exact overlaps: 1149035007
[M::ha_print_ovlp_stat] # inexact overlaps: 34624333
[M::ha_print_ovlp_stat] # overlaps without large indels: 1180991403
[M::ha_print_ovlp_stat] # reverse overlaps: 410771728
Writing reads to disk... 
Reads has been written.
Writing ma_hit_ts to disk... 
ma_hit_ts has been written.
Writing ma_hit_ts to disk... 
ma_hit_ts has been written.
bin files have been written.
[M::purge_dups] homozygous read coverage threshold: 168
[M::purge_dups] purge duplication coverage threshold: 210
Writing raw unitig GFA to disk... 
Writing processed unitig GFA to disk... 
[M::purge_dups] homozygous read coverage threshold: 168
[M::purge_dups] purge duplication coverage threshold: 210
[M::mc_solve_core::0.308] ==> Partition
[M::adjust_utg_by_primary] primary contig coverage range: [142, infinity]
Writing apul.hifiasm.intial.bp.p_ctg.gfa to disk... 
[M::adjust_utg_by_trio] primary contig coverage range: [142, infinity]
[M::adjust_utg_by_trio] primary contig coverage range: [142, infinity]
Writing apul.hifiasm.intial.bp.hap1.p_ctg.gfa to disk... 
[M::adjust_utg_by_trio] primary contig coverage range: [142, infinity]
Writing apul.hifiasm.intial.bp.hap2.p_ctg.gfa to disk... 
Inconsistency threshold for low-quality regions in BED files: 70%
[M::main] Version: 0.16.1-r375
[M::main] CMD: hifiasm -o apul.hifiasm.intial -t 36 hifi_rr_allcontam_rem.fasta
[M::main] Real time: 299500.413 sec; CPU: 10366612.139 sec; Peak RSS: 202.994 GB
```

Does the `M::purge_dups` mean that this is the value that I should be using for the `purge_dups` flag in hifiasm? It gives me two lines for `M::purge_dups`: `homozygous read coverage threshold: 168` and `purge duplication coverage threshold: 210`. I may need to talk with Ross and Hollie more about this. I'm going to QC the assembly and the haplotype assemblies with [busco](https://busco.ezlab.org/busco_userguide.html) and [quast](https://github.com/ablab/quast). In the scripts folder: `nano initial_qc.sh`

```
#!/bin/bash 
#SBATCH -t 100:00:00
#SBATCH --nodes=1 --ntasks-per-node=15
#SBATCH --export=NONE
#SBATCH --mem=250GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Apul_Genome/assembly/scripts
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

echo "Convert from gfa to fasta for downstream use" $(date)

cd /data/putnamlab/jillashey/Apul_Genome/assembly/data/

awk '/^S/{print ">"$2;print $3}' apul.hifiasm.intial.bp.p_ctg.gfa > apul.hifiasm.intial.bp.p_ctg.fa

awk '/^S/{print ">"$2;print $3}' apul.hifiasm.intial.bp.hap1.p_ctg.gfa > apul.hifiasm.intial.bp.hap1.p_ctg.fa

awk '/^S/{print ">"$2;print $3}' apul.hifiasm.intial.bp.hap2.p_ctg.gfa > apul.hifiasm.intial.bp.hap2.p_ctg.fa

echo "Begin busco on filtered hifiasm-assembled fasta (initial run)" $(date)

labbase=/data/putnamlab
busco_shared="${labbase}/shared/busco"
[ -z "$query" ] && query="${labbase}/jillashey/Apul_Genome/assembly/data/apul.hifiasm.intial.bp.p_ctg.fa" # set this to the query (genome/transcriptome) you are running
[ -z "$db_to_compare" ] && db_to_compare="${busco_shared}/downloads/lineages/metazoa_odb10"

source "${busco_shared}/scripts/busco_init.sh"  # sets up the modules required for this in the right order

# This will generate output under your $HOME/busco_output
cd "${labbase}/${Apul_Genome/assembly/data}"
busco --config "$EBROOTBUSCO/config/config.ini" -f -c 20 --long -i "${query}" -l metazoa_odb10 -o apul.initial.busco -m genome

echo "busco complete for unfiltered hifiasm-assembled fasta (initial run)" $(date)
echo "Begin busco on hifiasm-assembled haplotype 1 fasta" $(date)

# Reset query 
[ -z "$query" ] && query="${labbase}/jillashey/Apul_Genome/assembly/data/apul.hifiasm.intial.bp.hap1.p_ctg.fa" # set this to the query (genome/transcriptome) you are running

source "${busco_shared}/scripts/busco_init.sh"  # sets up the modules required for this in the right order

# This will generate output under your $HOME/busco_output
cd "${labbase}/${Apul_Genome/assembly/data}"
busco --config "$EBROOTBUSCO/config/config.ini" -f -c 20 --long -i "${query}" -l metazoa_odb10 -o apul.initial.hap1.busco -m genome

echo "busco complete for hifiasm-assembled haplotype 1 fasta (initial run)" $(date)
echo "Begin busco on hifiasm-assembled haplotype 2 fasta" $(date)

# Reset query 
[ -z "$query" ] && query="${labbase}/jillashey/Apul_Genome/assembly/data/apul.hifiasm.intial.bp.hap2.p_ctg.fa" # set this to the query (genome/transcriptome) you are running

source "${busco_shared}/scripts/busco_init.sh"  # sets up the modules required for this in the right order

# This will generate output under your $HOME/busco_output
cd "${labbase}/${Apul_Genome/assembly/data}"
busco --config "$EBROOTBUSCO/config/config.ini" -f -c 20 --long -i "${query}" -l metazoa_odb10 -o apul.initial.hap2.busco -m genome

echo "busco complete for hifiasm-assembled haplotype 2 fasta (initial run)" $(date)
echo "busco complete all assemblies of interest (initial run)" $(date)
echo "Begin quast of primary and haplotypes (initial run)" $(date)

module load QUAST/5.2.0-foss-2021b 
# there is another version of quast if the one above does not work: QUAST/5.0.2-foss-2020b-Python-2.7.18

quast -t 15 --eukaryote \
apul.hifiasm.intial.bp.p_ctg.fa \
apul.hifiasm.intial.bp.hap1.p_ctg.fa \
apul.hifiasm.intial.bp.hap2.p_ctg.fa \
/data/putnamlab/jillashey/genome/Amil_v2.01/Amil.v2.01.chrs.fasta \
-o /data/putnamlab/jillashey/Apul_Genome/assembly/output/quast

echo "Quast complete (initial run); all QC complete!" $(date)

```

Submitted batch job 309969. Ran in ~2 hours. Busco ran but only for the primariy assembly. For some reason, the hap1 and hap2 busco did not run. Quast appears to have failed with these errors: 

```
foss/2021b(24):ERROR:150: Module 'foss/2021b' conflicts with the currently loaded module(s) 'foss/2020b'
foss/2021b(24):ERROR:102: Tcl command execution failed: conflict foss

Python/3.9.6-GCCcore-11.2.0(61):ERROR:150: Module 'Python/3.9.6-GCCcore-11.2.0' conflicts with the currently loaded module(s) 'Python/3.8.6-GCCcore-10.2.0'
Python/3.9.6-GCCcore-11.2.0(61):ERROR:102: Tcl command execution failed: conflict Python

Perl/5.34.0-GCCcore-11.2.0(133):ERROR:150: Module 'Perl/5.34.0-GCCcore-11.2.0' conflicts with the currently loaded module(s) 'Perl/5.32.0-GCCcore-10.2.0'
Perl/5.34.0-GCCcore-11.2.0(133):ERROR:102: Tcl command execution failed: conflict Perl

foss/2021b(24):ERROR:150: Module 'foss/2021b' conflicts with the currently loaded module(s) 'foss/2020b'
foss/2021b(24):ERROR:102: Tcl command execution failed: conflict foss

Python/3.9.6-GCCcore-11.2.0(61):ERROR:150: Module 'Python/3.9.6-GCCcore-11.2.0' conflicts with the currently loaded module(s) 'Python/3.8.6-GCCcore-10.2.0'
Python/3.9.6-GCCcore-11.2.0(61):ERROR:102: Tcl command execution failed: conflict Python

SciPy-bundle/2021.10-foss-2021b(30):ERROR:150: Module 'SciPy-bundle/2021.10-foss-2021b' conflicts with the currently loaded module(s) 'SciPy-bundle/2020.11-foss-2020b'
SciPy-bundle/2021.10-foss-2021b(30):ERROR:102: Tcl command execution failed: conflict SciPy-bundle

GCCcore/11.2.0(24):ERROR:150: Module 'GCCcore/11.2.0' conflicts with the currently loaded module(s) 'GCCcore/10.2.0'
GCCcore/11.2.0(24):ERROR:102: Tcl command execution failed: conflict GCCcore
```

Basically just a lot of conflicting modules. Below are the results for the Busco code for the primary assembly: 

```
        --------------------------------------------------
        |Results from dataset metazoa_odb10               |
        --------------------------------------------------
        |C:93.3%[S:92.0%,D:1.3%],F:3.1%,M:3.6%,n:954      |
        |890    Complete BUSCOs (C)                       |
        |878    Complete and single-copy BUSCOs (S)       |
        |12     Complete and duplicated BUSCOs (D)        |
        |30     Fragmented BUSCOs (F)                     |
        |34     Missing BUSCOs (M)                        |
        |954    Total BUSCO groups searched               |
        --------------------------------------------------
```

Going to edit the `initial_qc.sh` script so that I am including all things for busco to run properly. I'm also commenting out the lines that ran successfully and switching the module to `QUAST/5.0.2-foss-2020b-Python-2.7.18`. Submitted batch job 309984. It ran, but only the hap1 busco scores were generated, not hap2. Below are the results for the hap1 assembly: 

```
        --------------------------------------------------
        |Results from dataset metazoa_odb10               |
        --------------------------------------------------
        |C:93.3%[S:92.8%,D:0.5%],F:3.0%,M:3.7%,n:954      |
        |890    Complete BUSCOs (C)                       |
        |885    Complete and single-copy BUSCOs (S)       |
        |5      Complete and duplicated BUSCOs (D)        |
        |29     Fragmented BUSCOs (F)                     |
        |35     Missing BUSCOs (M)                        |
        |954    Total BUSCO groups searched               |
        --------------------------------------------------
```

Quast also failed again with this error:

```
Python/2.7.18-GCCcore-10.2.0(58):ERROR:150: Module 'Python/2.7.18-GCCcore-10.2.0' conflicts with the currently loaded module(s) 'Python/3.8.6-GCCcore-10.2.0'
Python/2.7.18-GCCcore-10.2.0(58):ERROR:102: Tcl command execution failed: conflict Python

Python/2.7.18-GCCcore-10.2.0(58):ERROR:150: Module 'Python/2.7.18-GCCcore-10.2.0' conflicts with the currently loaded module(s) 'Python/3.8.6-GCCcore-10.2.0'
Python/2.7.18-GCCcore-10.2.0(58):ERROR:102: Tcl command execution failed: conflict Python
```

I'm going to add `module purge` and `module load Python/2.7.18-GCCcore-10.2.0` prior to loading quast. I'm also going to comment out lines that have already run successfully. Submitted batch job 310034. Ran in ~45 mins. Below are the results for the hap2 assembly: 

```
        --------------------------------------------------
        |Results from dataset metazoa_odb10               |
        --------------------------------------------------
        |C:94.0%[S:92.7%,D:1.3%],F:2.9%,M:3.1%,n:954      |
        |896    Complete BUSCOs (C)                       |
        |884    Complete and single-copy BUSCOs (S)       |
        |12     Complete and duplicated BUSCOs (D)        |
        |28     Fragmented BUSCOs (F)                     |
        |30     Missing BUSCOs (M)                        |
        |954    Total BUSCO groups searched               |
        --------------------------------------------------
```

Once again, quast failed to run and I got this error: 

```
ERROR! File not found (contigs): apul.hifiasm.intial.bp.p_ctg.fa
In case you have troubles running QUAST, you can write to quast.support@cab.spbu.ru
or report an issue on our GitHub repository https://github.com/ablab/quast/issues
Please provide us with quast.log file from the output directory.
```

I'm going to write quast its own script. Quast can be run with or without a reference genome. I'm going to try both, using Amillepora as the reference genome. In the scripts folder: `nano initial_quast.sh`

```
#!/bin/bash 
#SBATCH -t 100:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --export=NONE
#SBATCH --mem=250GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Apul_Genome/assembly/scripts
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

cd /data/putnamlab/jillashey/Apul_Genome/assembly/data

module purge
module load Python/2.7.18-GCCcore-10.2.0
module load QUAST/5.0.2-foss-2020b-Python-2.7.18
# previously used QUAST/5.2.0-foss-2021b but it failed and produced module conflict errors

echo "Begin quast of primary and haplotypes (initial run) w/ reference" $(date)

quast -t 10 --eukaryote \
apul.hifiasm.intial.bp.p_ctg.fa \
apul.hifiasm.intial.bp.hap1.p_ctg.fa \
apul.hifiasm.intial.bp.hap2.p_ctg.fa \
/data/putnamlab/jillashey/genome/Amil_v2.01/Amil.v2.01.chrs.fasta \
-o /data/putnamlab/jillashey/Apul_Genome/assembly/output/quast

echo "Quast complete (initial run); all QC complete!" $(date)
```

Submitted batch job 310038. I need to read more about [quast command line](https://quast.sourceforge.net/docs/manual.html#sec2) options, as there seem to be a lot of options. Also need to look into including a reference vs not including a reference. Ran super fast but success! Quast is super informative!!!!! The most useful information is in the `report.*` files. This is from the `report.txt` file: 

```
All statistics are based on contigs of size >= 500 bp, unless otherwise noted (e.g., "# contigs (>= 0 bp)" and "Total length (>= 0 bp)" include all contigs).

Assembly                    apul.hifiasm.intial.bp.p_ctg  apul.hifiasm.intial.bp.hap1.p_ctg  apul.hifiasm.intial.bp.hap2.p_ctg  Amil.v2.01.chrs
# contigs (>= 0 bp)         188                           275                                162                                854            
# contigs (>= 1000 bp)      188                           275                                162                                851            
# contigs (>= 5000 bp)      188                           273                                162                                748            
# contigs (>= 10000 bp)     186                           271                                162                                672            
# contigs (>= 25000 bp)     166                           246                                153                                545            
# contigs (>= 50000 bp)     98                            163                                124                                445            
Total length (>= 0 bp)      518528298                     481372407                          480341213                          475381253      
Total length (>= 1000 bp)   518528298                     481372407                          480341213                          475378544      
Total length (>= 5000 bp)   518528298                     481363561                          480341213                          475052084      
Total length (>= 10000 bp)  518514885                     481349140                          480341213                          474498957      
Total length (>= 25000 bp)  518188097                     480901871                          480181619                          472383091      
Total length (>= 50000 bp)  515726224                     477880931                          479141568                          468867721      
# contigs                   188                           275                                162                                854            
Largest contig              45111900                      21532546                           22038975                           39361238       
Total length                518528298                     481372407                          480341213                          475381253      
GC (%)                      39.05                         39.03                              39.04                              39.06          
N50                         16268372                      12353884                           13054353                           19840543       
N75                         13007972                      7901416                            8791894                            1469964        
L50                         11                            15                                 15                                 9              
L75                         20                            28                                 26                                 23             
# N's per 100 kbp           0.00                          0.00                               0.00                               7.79
```

Look at all of that info! The initial assembly appears to be the best in terms of all the stats. The total length is longer than the Amillepora and the haplotype assemblies. Additionally, it has the largest contig. The Amillepora genome has a higher N50 but the N50 for the primary assembly still looks good. There were 188 contigs generated in the primary assembly. Haplotype 1 assembly had more contigs (275), while haplotype 2 had less (162). Amillepora has 854 contigs which is so high! The primary assembly contig number is much lower than Atenuis (614), Adigitifera (955), or Amillepora (854). Quast, you have converted me <3. My next steps are to play with the `-s` flag in hifiasm to determine the threshold at which duplicate haplotigs should be purged. The default is 0.55 and Young et al. (2024) ran it with a range of values (0.55, 0.50, 0.45, 0.40, 0.35, 0.30). They found that all worked well to resolve haplotypes, so they stuck with the default of 0.55. They also used the `--primary` flag, which outputs a primary and alternate assembly as opposed to the primary, hap1 and hap2 assemblies. In their [code](https://github.com/benyoung93/orbicella_faveolata_pacbio_genome_transcriptome/blob/main/ofav_genome_pipeline.Rmd), they justified this by saying "running hifiasm using the primary flag as we have no real way of knowing if the haplotypes produced are real or not" (line 826). 
