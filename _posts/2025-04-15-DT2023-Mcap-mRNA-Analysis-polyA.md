---
layout: post
title: Developmental 2023 Timeseries mRNA analysis - polyA
date: '2025-04-15'
categories: Analysis
tags: [Bioinformatics, mRNA, Montipora capitata]
projects: Developmental Timseries 2023
---

## Developmental 2023 Timeseries mRNA analysis - polyA

These data came from my developmental timeseries experiment in 2023 with *Montipora capitata* in Hawaii. In this experiment, *Montipora capitata* embryos and larvae were exposed to ambient and heat stress over 72 hours from embryo to swimming larvae. Samples were collected at 8 time points: 1, 4, 9, 14, 22, 28, 48 and 72 hpf. The github for this project is [here](https://github.com/JillAshey/Hawaii_Developmental_TimeSeries). 

We sequenced these samples using rRNA depletion library methods and we wanted to compare that with libraries prepared using polyA selection. The library prep was done by me using the [Zymo 3' Switchfree](https://www.zymoresearch.com/products/zymo-seq-switchfree-3-mrna-library-kit) library prep kit (see my notebook [posts](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook) for more information on library prep). The sequencing was done through Oklahoma Medical Research Foundation NGS Core, which I found on Genohub. Here is the sequencing information: 

- Instrument: Illumina NovaSeq X Plus - 10B - PE 150 Cycle
- Read length: 2 x 150bp (Paired End)
- Pricing unit: per lane
- Number of samples: 8 libraries
- Guaranteed number of pass filter PE reads/sample: 20M (10M in each direction)
- Deliverables: FastQ files uploaded to Genohub project bucket

Files were downloaded to this location on Unity: XXXX. Goodbye Andromeda, hello Unity. See Unity [documentation](https://docs.unity.rc.umass.edu/documentation/) and Putnam lab Unity [documentation](https://github.com/Putnam-Lab/Lab_Management/blob/master/Bioinformatics_%26_Coding/Unity.md) (thanks Zoe :')) for info on using Unity. Time to analyze!

Check files to ensure transfer is complete and no files were corrupted. 

```
md5sum *fastq.gz > checkmd5.txt
md5sum -c checkmd5.txt 
M10_S136_R1_001.fastq.gz: OK
M10_S136_R2_001.fastq.gz: OK
M11_S137_R1_001.fastq.gz: OK
M11_S137_R2_001.fastq.gz: OK
M13_S138_R1_001.fastq.gz: OK
M13_S138_R2_001.fastq.gz: OK
M14_S139_R1_001.fastq.gz: OK
M14_S139_R2_001.fastq.gz: OK
M23_S140_R1_001.fastq.gz: OK
M23_S140_R2_001.fastq.gz: OK
M24_S141_R1_001.fastq.gz: OK
M24_S141_R2_001.fastq.gz: OK
M26_S142_R1_001.fastq.gz: OK
M26_S142_R2_001.fastq.gz: OK
M28_S143_R1_001.fastq.gz: OK
M28_S143_R2_001.fastq.gz: OK
M35_S144_R1_001.fastq.gz: OK
M35_S144_R2_001.fastq.gz: OK
M36_S145_R1_001.fastq.gz: OK
M36_S145_R2_001.fastq.gz: OK
M37_S146_R1_001.fastq.gz: OK
M37_S146_R2_001.fastq.gz: OK
M39_S147_R1_001.fastq.gz: OK
M39_S147_R2_001.fastq.gz: OK
M47_S148_R1_001.fastq.gz: OK
M47_S148_R2_001.fastq.gz: OK
M48_S149_R1_001.fastq.gz: OK
M48_S149_R2_001.fastq.gz: OK
M51_S150_R1_001.fastq.gz: OK
M51_S150_R2_001.fastq.gz: OK
M52_S151_R1_001.fastq.gz: OK
M52_S151_R2_001.fastq.gz: OK
M60_S152_R1_001.fastq.gz: OK
M60_S152_R2_001.fastq.gz: OK
M61_S153_R1_001.fastq.gz: OK
M61_S153_R2_001.fastq.gz: OK
M62_S154_R1_001.fastq.gz: OK
M62_S154_R2_001.fastq.gz: OK
M63_S155_R1_001.fastq.gz: OK
M63_S155_R2_001.fastq.gz: OK
M6_S132_R1_001.fastq.gz: OK
M6_S132_R2_001.fastq.gz: OK
M72_S156_R1_001.fastq.gz: OK
M72_S156_R2_001.fastq.gz: OK
M73_S157_R1_001.fastq.gz: OK
M73_S157_R2_001.fastq.gz: OK
M74_S158_R1_001.fastq.gz: OK
M74_S158_R2_001.fastq.gz: OK
M75_S159_R1_001.fastq.gz: OK
M75_S159_R2_001.fastq.gz: OK
M7_S133_R1_001.fastq.gz: OK
M7_S133_R2_001.fastq.gz: OK
M85_S160_R1_001.fastq.gz: OK
M85_S160_R2_001.fastq.gz: OK
M86_S161_R1_001.fastq.gz: OK
M86_S161_R2_001.fastq.gz: OK
M87_S162_R1_001.fastq.gz: OK
M87_S162_R2_001.fastq.gz: OK
M88_S163_R1_001.fastq.gz: OK
M88_S163_R2_001.fastq.gz: OK
M8_S134_R1_001.fastq.gz: OK
M8_S134_R2_001.fastq.gz: OK
M9_S135_R1_001.fastq.gz: OK
M9_S135_R2_001.fastq.gz: OK
```

Make directories on Unity. Sym link the raw data in the project directory to the work directory. 

```
cd /work/pi_hputnam_uri_edu/jillashey/Mcap_2023/polyA/data/raw
ln -s /project/pi_hputnam_uri_edu/jillashey/Mcap_2023/polyA/*fastq.gz .
```

Run fastqc on raw reads. `nano raw_fastqc.sh`

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --signal=2
#SBATCH --no-requeue
#SBATCH --mem=200GB
#SBATCH -t 48:00:00
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error
#SBATCH -D /work/pi_hputnam_uri_edu/jillashey/Mcap_2023/polyA/scripts 

# load modules needed
module load parallel/20240822
module load fastqc/0.12.1
module load uri/main
module load all/MultiQC/1.12-foss-2021b

cd /work/pi_hputnam_uri_edu/jillashey/Mcap_2023/polyA/data/raw

# Create an array of fastq files to process
files=($('ls' *001.fastq.gz)) 

# Run fastqc in parallel
echo "Starting fastqc..." $(date)
parallel -j 20 "fastqc {} -o /work/pi_hputnam_uri_edu/jillashey/Mcap_2023/polyA/output/fastqc/raw && echo 'Processed {}'" ::: "${files[@]}"
echo "fastQC done." $(date)

cd /work/pi_hputnam_uri_edu/jillashey/Mcap_2023/polyA/output/fastqc/raw

echo "Starting multiqc..." $(date)
multiqc *

echo "Initial QC of polyA data complete." $(date)
```

Submitted batch job 32423847. Raw QC data is [here](https://github.com/JillAshey/DevelopmentalTimeseries/blob/main/data/Molecular/mRNA_polyA/raw_multiqc_report_polyA.html). Time to trim!

This info is from the Zymo switchfree [protocol](https://www.zymoresearch.com/products/zymo-seq-switchfree-3-mrna-library-kit):

Read 1 provides UMI information and Read 2 provides the transcript sequence. Read 1 will sequence the 8-nucleotide UMI, the 6 non-T random nucleotides, and the oligo dT region sequentially before sequencing the insert region. Therefore, it is recommended to assign no more than 25 cycles to Read 1. This is sufficient to allow for cluster identification, quality metrics calculation, and coverage of UMI information while avoiding more low-quality bases associated with sequencing through the oligo dT region.

It is possible to sequence longer for Read 1 (e.g., 100 bp or 150 bp pairedend sequencing) when other libraries of high complexity are sequenced on the same lane with the Zymo-Seq SwitchFree™ 3′ mRNA libraries. Read 2 will provide the transcript sequences. 

Generally, it is sufficient to trim Read 2 with the Illumina® TruSeq® adapter sequence: “AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT”. When the read length is longer than the average library insert size, Read 2 will possibly sequence into the oligo dT sequence, the random hexamer, and even the UMI sequence. Additionally, poly A tail may also present at the end of Read 2. These additional bases can be trimmed to potentially improve alignment. An example using Cutadapt2 for such 2-step trimming
is shown as below.

```
cutadapt -a A{8}B{6}N{8}AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o
adapterTrimmed.fastq.gz sample.fastq.gz
cutadapt -a “A{100}” -o completeTrimmed.fastq.gz adapterTrimmed.fastq.gz
```

We highly recommend discarding Read 1 after UMI extraction in any bioinformatic analysis because of the poor sequence quality associated with the oligo dT region.

Hollie also went through some trimming iterations on her switchfree data. See her code [here](https://github.com/hputnam/Poc_RAPID/blob/main/RAnalysis/data/Tagseq_zymo_genohub/20240424_Pverrucosa_tagseq_workflow.md). In the scripts folder: `nano trim_R2_cutadapt.sh`

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=15
#SBATCH --partition=uri-cpu
#SBATCH --no-requeue
#SBATCH --mem=300GB
#SBATCH -t 24:00:00
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error
#SBATCH -D /work/pi_hputnam_uri_edu/jillashey/Mcap_2023/polyA/scripts 

# load modules needed
module load parallel/20240822
module load fastqc/0.12.1
module load uri/main
module load all/cutadapt/3.5-GCCcore-11.2.0
module load all/MultiQC/1.12-foss-2021b

cd /work/pi_hputnam_uri_edu/jillashey/Mcap_2023/polyA/data/raw

# Create an array of fastq files to process
array=($(ls *R2_001.fastq.gz)) 

echo "trim using cutadapt" $(date)
for i in "${array[@]}"; do 
    sample_name=$(echo "$i" | awk -F '_' '{print $1}')  
    
    # Use original filename ($i) as input
    cutadapt -a A{8}B{6}N{8}AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o "${sample_name}_adapterTrimmed.fastq.gz" "$i"
    cutadapt -a A{100} -o "${sample_name}_completeTrimmed.fastq.gz" "${sample_name}_adapterTrimmed.fastq.gz"
    
    echo "trimming complete" $(date)
    
    # Move trimmed file
    mv "${sample_name}_completeTrimmed.fastq.gz" /work/pi_hputnam_uri_edu/jillashey/Mcap_2023/polyA/data/trim/
    
    # FastQC
    fastqc "/work/pi_hputnam_uri_edu/jillashey/Mcap_2023/polyA/data/trim/${sample_name}_completeTrimmed.fastq.gz" \
        -o /work/pi_hputnam_uri_edu/jillashey/Mcap_2023/polyA/output/fastqc/trim
done

# Run MultiQC once after all samples are processed
cd /work/pi_hputnam_uri_edu/jillashey/Mcap_2023/polyA/output/fastqc/trim
echo "Starting multiqc..." $(date)
multiqc *
echo "QC of trimmed polyA data complete." $(date)
```

Submitted batch job 32583551. Lots of polyG still left and adapter content. Discussed with Zoe. Based on cutadapt [manual](https://cutadapt.readthedocs.io/en/stable/guide.html#adapter-types), I need to add `--nextseq-trim=20` because: "Quality trimming of reads using two-color chemistry (NextSeq) Some Illumina instruments use a two-color chemistry to encode the four bases. This includes the NextSeq and the NovaSeq. In those instruments, a ‘dark cycle’ (with no detected color) encodes a G. However, dark cycles also occur when sequencing “falls off” the end of the fragment. The read then contains a run of high-quality, but incorrect “G” calls at its 3’ end."

I am also going to do some quality trimming to toss any bad reads (ie <20 phred score). Because cutadapt does the adapter removal first, then the polyG trim, I will need to do two trimming iterations, similar to what I did above. Also going to do both reads this time. In the scripts folder: `nano trim_cutadapt.sh`

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=15
#SBATCH --partition=uri-cpu
#SBATCH --no-requeue
#SBATCH --mem=300GB
#SBATCH -t 24:00:00
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error
#SBATCH -D /work/pi_hputnam_uri_edu/jillashey/Mcap_2023/polyA/scripts 

# load modules needed
module load parallel/20240822
module load fastqc/0.12.1
module load uri/main
module load all/cutadapt/3.5-GCCcore-11.2.0
module load all/MultiQC/1.12-foss-2021b

cd /work/pi_hputnam_uri_edu/jillashey/Mcap_2023/polyA/data/raw

# Create an array of fastq files to process
array=($(ls *001.fastq.gz)) 

echo "trim using cutadapt" $(date)
for i in "${array[@]}"; do 
	sample_name=$(echo "$i" | awk -F'_' '{print $1 "_" $3}')
	cutadapt --nextseq-trim=20 -q 20,20 --minimum-length 20 -o "${sample_name}_PolyGTrimmed.fastq.gz" "$i"
	cutadapt -a A{8}B{6}N{8}AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o "${sample_name}_AdapterPolyGTrimmed.fastq.gz" "${sample_name}_PolyGTrimmed.fastq.gz"
	cutadapt -a A{100} -o "${sample_name}_completeTrimmed.fastq.gz" "${sample_name}_AdapterPolyGTrimmed.fastq.gz"
   
   echo "trimming complete" $(date)
    
    # Move trimmed file
    mv "${sample_name}_completeTrimmed.fastq.gz" /work/pi_hputnam_uri_edu/jillashey/Mcap_2023/polyA/data/trim/
    
    # FastQC
    fastqc "/work/pi_hputnam_uri_edu/jillashey/Mcap_2023/polyA/data/trim/${sample_name}_completeTrimmed.fastq.gz" \
        -o /work/pi_hputnam_uri_edu/jillashey/Mcap_2023/polyA/output/fastqc/trim
done

# Run MultiQC once after all samples are processed

cd /work/pi_hputnam_uri_edu/jillashey/Mcap_2023/polyA/output/fastqc/trim

echo "Starting multiqc..." $(date)
multiqc *
echo "Multiqc complete." $(date)
```

Submitted batch job 32648465. Still a lot of adapter content. Discussed with Zoe. I have paired end reads but only really care about R2. I still need to provide cutadapt with two input files. Since I am only interested in trimming R2 and , I should be using `-A`. 

Also beginning to scratch directories on Unity (see Unity [docuemntation](https://docs.unity.rc.umass.edu/documentation/managing-files/hpc-workspace/) and Putnam lab [documentation](https://github.com/Putnam-Lab/Lab_Management/blob/master/Bioinformatics_%26_Coding/Unity.md)). 

In the scripts folder: `nano trim_round2_cutadapt.sh`

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=15
#SBATCH --partition=uri-cpu
#SBATCH --no-requeue
#SBATCH --mem=300GB
#SBATCH -t 24:00:00
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error
#SBATCH -D /work/pi_hputnam_uri_edu/jillashey/Mcap_2023/polyA/scripts 

# load modules needed
module load parallel/20240822
module load fastqc/0.12.1
module load uri/main
module load all/cutadapt/3.5-GCCcore-11.2.0
module load all/MultiQC/1.12-foss-2021b

cd /work/pi_hputnam_uri_edu/jillashey/Mcap_2023/polyA/data/raw

# Make array of files to process
R1_raw=(*R1*.fastq.gz)
R2_raw=(*R2*.fastq.gz)

for i in "${!R1_raw[@]}"; do
    sample_name=$(basename "${R1_raw[$i]}" | sed -E 's/_S[0-9]+_R1_001.fastq.gz//')

    # Step 1: PolyG trimming
    cutadapt --nextseq-trim=20 -q 20,20 --minimum-length=20 \
        -o "${sample_name}_R1_PolyGTrimmed.fastq.gz" \
        -p "${sample_name}_R2_PolyGTrimmed.fastq.gz" \
        "${R1_raw[$i]}" "${R2_raw[$i]}"

    # Step 2: Adapter trimming
    cutadapt -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -a A{100} -q 20,20 --minimum-length=20 \
        -o "${sample_name}_R1_AdapterPolyGTrimmed.fastq.gz" \
        -p "${sample_name}_R2_AdapterPolyGTrimmed.fastq.gz" \
        "${sample_name}_R1_PolyGTrimmed.fastq.gz" "${sample_name}_R2_PolyGTrimmed.fastq.gz"

    echo "trimming complete for $sample_name" $(date)

    # Move only the correct files
    mv "${sample_name}_R1_AdapterPolyGTrimmed.fastq.gz" "${sample_name}_R2_AdapterPolyGTrimmed.fastq.gz" \
       /scratch/workspace/jillashey_uri_edu-ashey_scratch/Mcap2023/trim_data_round2

    # FastQC on the final adapter+PolyG trimmed files
    fastqc "/scratch/workspace/jillashey_uri_edu-ashey_scratch/Mcap2023/trim_data_round2/${sample_name}_R1_AdapterPolyGTrimmed.fastq.gz" \
           -o /work/pi_hputnam_uri_edu/jillashey/Mcap_2023/polyA/output/fastqc/trim_round2

    fastqc "/scratch/workspace/jillashey_uri_edu-ashey_scratch/Mcap2023/trim_data_round2/${sample_name}_R2_AdapterPolyGTrimmed.fastq.gz" \
           -o /work/pi_hputnam_uri_edu/jillashey/Mcap_2023/polyA/output/fastqc/trim_round2

done

# Run MultiQC once after all samples are processed

cd /work/pi_hputnam_uri_edu/jillashey/Mcap_2023/polyA/output/fastqc/trim_round2

echo "Starting multiqc..." $(date)
multiqc *
echo "Multiqc complete." $(date)
```

Submitted batch job 32744961. The cutadapt version on Unity is 3.5 which does not have the capacity for polyA tail trimming. I installed the updated version via conda on Unity. See cutadapt installation [instructions](https://cutadapt.readthedocs.io/en/stable/installation.html) and Unity information about [conda environments](https://docs.unity.rc.umass.edu/documentation/software/conda/). 

```
cd /work/pi_hputnam_uri_edu/conda/envs
module load conda/latest # need to load before making any conda envs
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
conda create --prefix /work/pi_hputnam_uri_edu/conda/envs/cutadapt cutadapt

conda activate /work/pi_hputnam_uri_edu/conda/envs/cutadapt 
cutadapt --version
5.0                                                                           
```

Successful installation. Okay lets try running this. In the scripts folder: `nano trim_round2_cutadaptv5.sh`

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

# load modules needed
module load parallel/20240822
module load fastqc/0.12.1

# Activate conda env
module load conda/latest 
conda activate /work/pi_hputnam_uri_edu/conda/envs/cutadapt 

cd /work/pi_hputnam_uri_edu/jillashey/Mcap_2023/polyA/data/raw

# Array of target samples (MXX numbers only)
target_samples=("M6" "M7" "M8" "M9" "M10" "M11" "M13" "M14" "M23" "M24" "M26" "M28" "M35" "M36" "M37" "M39" "M47" "M48" "M51" "M52" "M60" "M61" "M62" "M63" "M72" "M73" "M74" "M75" "M85" "M86" "M87" "M88")

# Process each target sample
for sample in "${target_samples[@]}"; do
    # Find corresponding raw files
    R1_raw=(${sample}_S*_R1_001.fastq.gz)
    R2_raw=(${sample}_S*_R2_001.fastq.gz)
    
    # Verify files exist
    if [ ${#R1_raw[@]} -eq 1 ] && [ ${#R2_raw[@]} -eq 1 ]; then
        # Step 1: PolyG trimming
        cutadapt --nextseq-trim=20 -q 20,20 --minimum-length=20 \
            -o "${sample}_R1_PolyGTrimmed.fastq.gz" \
            -p "${sample}_R2_PolyGTrimmed.fastq.gz" \
            "${R1_raw[0]}" "${R2_raw[0]}"

        # Step 2: Adapter trimming
        cutadapt -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --poly-a -q 20,20 --minimum-length=20 \
            -o "${sample}_R1_AdapterPolyGTrimmed.fastq.gz" \
            -p "${sample}_R2_AdapterPolyGTrimmed.fastq.gz" \
            "${sample}_R1_PolyGTrimmed.fastq.gz" "${sample}_R2_PolyGTrimmed.fastq.gz"

        echo "Trimming complete for ${sample}" $(date)

        # Move final files
        mv "${sample}_R1_AdapterPolyGTrimmed.fastq.gz" "${sample}_R2_AdapterPolyGTrimmed.fastq.gz" \
           /scratch/workspace/jillashey_uri_edu-ashey_scratch/Mcap2023/trim_data_round2

        # FastQC analysis
        fastqc "/scratch/workspace/jillashey_uri_edu-ashey_scratch/Mcap2023/trim_data_round2/${sample}_R1_AdapterPolyGTrimmed.fastq.gz" \
               -o /work/pi_hputnam_uri_edu/jillashey/Mcap_2023/polyA/output/fastqc/trim_round2

        fastqc "/scratch/workspace/jillashey_uri_edu-ashey_scratch/Mcap2023/trim_data_round2/${sample}_R2_AdapterPolyGTrimmed.fastq.gz" \
               -o /work/pi_hputnam_uri_edu/jillashey/Mcap_2023/polyA/output/fastqc/trim_round2
    else
        echo "Error: Found ${#R1_raw[@]} R1 and ${#R2_raw[@]} R2 files for ${sample}"
    fi
done

# Deactivate conda env
conda deactivate

# Load MultiQC modules 
module load uri/main
module load all/MultiQC/1.12-foss-2021b

cd /work/pi_hputnam_uri_edu/jillashey/Mcap_2023/polyA/output/fastqc/trim_round2
multiqc *

echo "MultiQC complete!" $(date)
```

Submitted batch job 33026985. Trimming occurred! But I forgot to create the `/scratch/workspace/jillashey_uri_edu-ashey_scratch/Mcap2023/trim_data_round2/` directory so no fastqc was done. Let's do that in an interactive session 

```
salloc

mv /work/pi_hputnam_uri_edu/jillashey/Mcap_2023/polyA/data/raw/*AdapterPolyGTrimmed.fastq.gz /scratch/workspace/jillashey_uri_edu-ashey_scratch/Mcap2023/trim_data_round2
```

Run QC as a job. I also accidently deleted the M85 R1 but that is okay for now. `nano trim_round2_QC.sh`

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

# load modules needed
#module load parallel/20240822
module load fastqc/0.12.1

cd /scratch/workspace/jillashey_uri_edu-ashey_scratch/Mcap2023/trim_data_round2

fastqc *AdapterPolyGTrimmed.fastq.gz -o /work/pi_hputnam_uri_edu/jillashey/Mcap_2023/polyA/output/fastqc/trim_round2 

echo "fastqc complete" $(date)

# Load MultiQC modules 
module load uri/main
module load all/MultiQC/1.12-foss-2021b

cd /work/pi_hputnam_uri_edu/jillashey/Mcap_2023/polyA/output/fastqc/trim_round2
multiqc *

echo "MultiQC complete!" $(date)
```

Submitted batch job 33041843. MultiQC report can be found [here](https://github.com/JillAshey/DevelopmentalTimeseries/blob/main/data/Molecular/mRNA_polyA/trim_round2_multiqc_report.html). There is still polyA tail on the R2 samples...and many of them have a lot of reads ~30bp. wtf????? I'm going to run QC on the intermediate cutadapt files (`*PolyGTrimmed.fastq.gz`) to see what the differences are. 

In the scripts folder: `nano trim_round2_QC_PolyGTrimmed.sh`

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

# load modules needed
#module load parallel/20240822
module load fastqc/0.12.1

cd /work/pi_hputnam_uri_edu/jillashey/Mcap_2023/polyA/data/raw

fastqc *PolyGTrimmed.fastq.gz -o /work/pi_hputnam_uri_edu/jillashey/Mcap_2023/polyA/output/fastqc/trim_round2_PolyGTrimmed 

echo "fastqc complete" $(date)

# Load MultiQC modules 
module load uri/main
module load all/MultiQC/1.12-foss-2021b

cd /work/pi_hputnam_uri_edu/jillashey/Mcap_2023/polyA/output/fastqc/trim_round2_PolyGTrimmed
multiqc *

echo "MultiQC complete!" $(date)
```

Submitted batch job 33173364. Here, we see that the reads are still 150bp long and adapter is still present. 


Let's try trimming `*PolyGTrimmed.fastq.gz` by specifying the length of polyA we want to trim. `nano trim_round3_cutadaptv5.sh`

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

# load modules needed
module load parallel/20240822
module load fastqc/0.12.1

# Activate conda env
module load conda/latest 
conda activate /work/pi_hputnam_uri_edu/conda/envs/cutadapt 

cd /work/pi_hputnam_uri_edu/jillashey/Mcap_2023/polyA/data/raw

# Array of target samples (MXX numbers only)
target_samples=("M6" "M7" "M8" "M9" "M10" "M11" "M13" "M14" "M23" "M24" "M26" "M28" "M35" "M36" "M37" "M39" "M47" "M48" "M51" "M52" "M60" "M61" "M62" "M63" "M72" "M73" "M74" "M75" "M85" "M86" "M87" "M88")

echo "Trimming started at" $(date)

# Process each target sample
for sample in "${target_samples[@]}"; do
    # Find corresponding raw files
    R1_polyg_trim=(${sample}_R1_PolyGTrimmed.fastq.gz)
    R2_polyg_trim=(${sample}_R2_PolyGTrimmed.fastq.gz)
    
    # Run cutadapt
    cutadapt \
        -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
        -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
        -a "A{20}" \
        -A "A{20}" \
        -q 20,20 \
        --minimum-length=20 \
        -o "${sample}_R1_AdapterPolyATrimmed.fastq.gz" \
        -p "${sample}_R2_AdapterPolyATrimmed.fastq.gz" \
        "${R1_polyg_trim[0]}" "${R2_polyg_trim[0]}"
done

echo "Trimming complete, starting QC" $(date)   

fastqc *AdapterPolyGATrimmed.fastq.gz

# Load MultiQC modules 
module load uri/main
module load all/MultiQC/1.12-foss-2021b
multiqc *

echo "MultiQC complete!" $(date)   
```

Submitted batch job 33316139. QC didn't work. `nano trim_round3_QC.sh`

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

# load modules needed
#module load parallel/20240822
module load fastqc/0.12.1

cd /scratch/workspace/jillashey_uri_edu-ashey_scratch/Mcap2023/trim_round3

fastqc *AdapterPolyATrimmed.fastq.gz -o /work/pi_hputnam_uri_edu/jillashey/Mcap_2023/polyA/output/fastqc/trim_round3

echo "fastqc complete" $(date)

# Load MultiQC modules 
module load uri/main
module load all/MultiQC/1.12-foss-2021b

cd /work/pi_hputnam_uri_edu/jillashey/Mcap_2023/polyA/output/fastqc/trim_round3
multiqc *

echo "MultiQC complete!" $(date)
```

Submitted batch job 33320990. Still not removing all of the polyA tail. Use A10 a10 instead of A20. `nano trim_round4_cutadaptv5.sh`

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

# load modules needed
module load parallel/20240822
module load fastqc/0.12.1

# Activate conda env
module load conda/latest 
conda activate /work/pi_hputnam_uri_edu/conda/envs/cutadapt 

cd /work/pi_hputnam_uri_edu/jillashey/Mcap_2023/polyA/data/raw

# Array of target samples (MXX numbers only)
target_samples=("M6" "M7" "M8" "M9" "M10" "M11" "M13" "M14" "M23" "M24" "M26" "M28" "M35" "M36" "M37" "M39" "M47" "M48" "M51" "M52" "M60" "M61" "M62" "M63" "M72" "M73" "M74" "M75" "M85" "M86" "M87" "M88")

echo "Trimming started at" $(date)

# Process each target sample
for sample in "${target_samples[@]}"; do
    # Find corresponding raw files
    R1_polyg_trim=(${sample}_R1_PolyGTrimmed.fastq.gz)
    R2_polyg_trim=(${sample}_R2_PolyGTrimmed.fastq.gz)
    
    # Run cutadapt
    cutadapt \
        -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
        -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
        -a "A{10}" \
        -A "A{10}" \
        -q 20,20 \
        --minimum-length=20 \
        -o "${sample}_R1_AdapterPolyA10Trimmed.fastq.gz" \
        -p "${sample}_R2_AdapterPolyA10Trimmed.fastq.gz" \
        "${R1_polyg_trim[0]}" "${R2_polyg_trim[0]}"
done
```

Submitted batch job 33322136. Run QC. `nano trim_round4_QC.sh`

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

# load modules needed
#module load parallel/20240822
module load fastqc/0.12.1

mv /work/pi_hputnam_uri_edu/jillashey/Mcap_2023/polyA/data/raw/*AdapterPolyA10Trimmed.fastq.gz /scratch/workspace/jillashey_uri_edu-ashey_scratch/Mcap2023/trim_round4

cd /scratch/workspace/jillashey_uri_edu-ashey_scratch/Mcap2023/trim_round4

fastqc *AdapterPolyA10Trimmed.fastq.gz -o /work/pi_hputnam_uri_edu/jillashey/Mcap_2023/polyA/output/fastqc/trim_round4

echo "fastqc complete" $(date)

# Load MultiQC modules 
module load uri/main
module load all/MultiQC/1.12-foss-2021b

cd /work/pi_hputnam_uri_edu/jillashey/Mcap_2023/polyA/output/fastqc/trim_round4
multiqc *

echo "MultiQC complete!" $(date)
```

Submitted batch job 33328319. 

may need to run the adapter and polyA trimming parts separately... `nano trim_round5_cutadaptv5.sh`

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

# load modules needed
module load parallel/20240822
module load fastqc/0.12.1

# Activate conda env
module load conda/latest 
conda activate /work/pi_hputnam_uri_edu/conda/envs/cutadapt 

cd /work/pi_hputnam_uri_edu/jillashey/Mcap_2023/polyA/data/raw

# Array of target samples (MXX numbers only)
target_samples=("M6" "M7" "M8" "M9" "M10" "M11" "M13" "M14" "M23" "M24" "M26" "M28" "M35" "M36" "M37" "M39" "M47" "M48" "M51" "M52" "M60" "M61" "M62" "M63" "M72" "M73" "M74" "M75" "M85" "M86" "M87" "M88")

# Process each target sample
for sample in "${target_samples[@]}"; do
    # Find corresponding raw files
    R1_raw=(${sample}_S*_R1_001.fastq.gz)
    R2_raw=(${sample}_S*_R2_001.fastq.gz)
    
    # Verify files exist
    if [ ${#R1_raw[@]} -eq 1 ] && [ ${#R2_raw[@]} -eq 1 ]; then
        # Step 1: PolyG trimming
        cutadapt --nextseq-trim=20 -q 20,20 --minimum-length=20 \
            -o "${sample}_R1_PolyGTrimmed.fastq.gz" \
            -p "${sample}_R2_PolyGTrimmed.fastq.gz" \
            "${R1_raw[0]}" "${R2_raw[0]}"
          echo "polyG trimming complete for ${sample}" $(date)

		# Step 2: adapter trimming 
    cutadapt \
        -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
        -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
        -q 20,20 \
        --minimum-length=20 \
        -o "${sample}_R1_AdapterPolyGTrimmed.fastq.gz" \
        -p "${sample}_R2_AdapterPolyGTrimmed.fastq.gz" \
           "${sample}_R1_PolyGTrimmed.fastq.gz" "${sample}_R2_PolyGTrimmed.fastq.gz"
     echo "Adapter trimming complete for ${sample}" $(date)
     
     		# Step 3: polyA trimming 
	cutadapt \
        -a "A{10}" \
        -A "A{10}" \
        -q 20,20 \
        --minimum-length=20 \
        -o "${sample}_R1_AdapterPolyGA10Trimmed.fastq.gz" \
        -p "${sample}_R2_AdapterPolyGA10Trimmed.fastq.gz" \
        "${sample}_R1_AdapterPolyGTrimmed.fastq.gz" "${sample}_R2_AdapterPolyGTrimmed.fastq.gz"
      echo "PolyA trimming complete for ${sample}" $(date)
	fi
done

echo "Trimming complete, run fastQC" $(date)

mv *Trimmed.fastq.gz /scratch/workspace/jillashey_uri_edu-ashey_scratch/Mcap2023/trim_round5

cd /scratch/workspace/jillashey_uri_edu-ashey_scratch/Mcap2023/trim_round5

fastqc *AdapterPolyGA10Trimmed.fastq.gz -o /work/pi_hputnam_uri_edu/jillashey/Mcap_2023/polyA/output/fastqc/trim_round5

echo "fastqc complete" $(date)

# Load MultiQC modules 
module load uri/main
module load all/MultiQC/1.12-foss-2021b

cd /work/pi_hputnam_uri_edu/jillashey/Mcap_2023/polyA/output/fastqc/trim_round5
multiqc *

echo "MultiQC complete!" $(date)
```

Submitted batch job 33332147. FINALLY!!!!! SUCCESS!!!!!!! 

TIME TO ALIGN!!!! In the scripts folder: `nano align_R2.sh`

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

## Genome has already been indexed for hisat2

# load modules needed
module load uri/main
module load all/HISAT2/2.2.1-gompi-2022a
module load all/SAMtools/1.18-GCC-12.3.0

cd /scratch/workspace/jillashey_uri_edu-ashey_scratch/Mcap2023/trim_round5

echo "Aligning polyA R2 samples!" $(date)

for READ in *_R2_AdapterPolyGA10Trimmed.fastq.gz; do
    ID=$(echo "$READ" | cut -d'_' -f1)

    # Run alignment
    hisat2 --rna-strandness F -p 8 --dta -x /work/pi_hputnam_uri_edu/HI_Genomes/MCapV3/McapV3_hisat2_ref -U "$READ" -S "${ID}.sam"
    
    # Convert to BAM and sort
    samtools sort -@ 8 -o "${ID}.bam" "${ID}.sam"
    echo "${ID} bam-ified!"
    
    # Remove SAM file
    rm "${ID}.sam"
done

echo "PolyA R2 samples alignment complete!" $(date)
```

Submitted batch job 33409631

Assemble! `nano assemble_R2.sh`

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

# load modules needed
module load uri/main
module load all/StringTie/2.2.1-GCC-11.2.0

echo "Assembling transcripts using stringtie" $(date)

cd /scratch/workspace/jillashey_uri_edu-ashey_scratch/Mcap2023/trim_round5

array1=($(ls *.bam)) #Make an array of sequences to assemble

for i in ${array1[@]}; do
    stringtie -p 8 -e -B -G /work/pi_hputnam_uri_edu/HI_Genomes/MCapV3/Montipora_capitata_HIv3.genes_fixed.gff3 -A ${i}.gene_abund.tab -o ${i}.gtf ${i}
done

echo "Assembly complete " $(date)
```

Submitted batch job 33420828. TPM information can be seen in the gtf files. 

Interactive mode! 

```
cd /scratch/workspace/jillashey_uri_edu-ashey_scratch/Mcap2023/trim_round5
salloc 
module load uri/main
module load all/StringTie/2.2.1-GCC-11.2.0

# List of gtfs
ls *.gtf > gtf_list.txt

# Merge gtfs into single gtf 
stringtie --merge -e -p 8 -G /work/pi_hputnam_uri_edu/HI_Genomes/MCapV3/Montipora_capitata_HIv3.genes_fixed.gff3 -o Mcap_polyA_merged.gtf gtf_list.txt #Merge GTFs 

module load all/GffCompare/0.12.6-GCC-11.2.0
gffcompare -r /work/pi_hputnam_uri_edu/HI_Genomes/MCapV3/Montipora_capitata_HIv3.genes_fixed.gff3 -G -o merged Mcap_polyA_merged.gtf #Compute the accuracy 
54384 reference transcripts loaded.
  2 duplicate reference transcripts discarded.
  54384 query transfrags loaded.
```

Look at merge stats 

```

# gffcompare v0.12.6 | Command line was:
#gffcompare -r /work/pi_hputnam_uri_edu/HI_Genomes/MCapV3/Montipora_capitata_HIv3.genes_fixed.gff3 -G -o merged Mcap_polyA_merged.gtf
#
#= Summary for dataset: Mcap_polyA_merged.gtf 
#     Query mRNAs :   54384 in   54185 loci  (36023 multi-exon transcripts)
#            (141 multi-transcript loci, ~1.0 transcripts per locus)
# Reference mRNAs :   54382 in   54185 loci  (36023 multi-exon)
# Super-loci w/ reference transcripts:    54185
#-----------------| Sensitivity | Precision  |
        Base level:   100.0     |   100.0    |
        Exon level:   100.0     |   100.0    |
      Intron level:   100.0     |   100.0    |
Intron chain level:   100.0     |   100.0    |
  Transcript level:    99.9     |    99.9    |
       Locus level:   100.0     |   100.0    |
     Matching intron chains:   36023
       Matching transcripts:   54349
              Matching loci:   54183
          Missed exons:       0/256029  (  0.0%)
           Novel exons:       0/256028  (  0.0%)
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

Run as a job. `nano prepDE.sh`

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
module load Python/2.7.18-GCCcore-9.3.0

cd /scratch/workspace/jillashey_uri_edu-ashey_scratch/Mcap2023/trim_round5

echo "Creating count matrices " $(date)

python /work/pi_hputnam_uri_edu/jillashey/Mcap_2023/egg_sperm/scripts/prepDE.py -g Mcap_polyA_gene_count_matrix.csv -i listGTF.txt

echo "Count matrices complete" $(date)
```

Submitted batch job 33422229. Move csvs to output folder and download to computer

```
cd /scratch/workspace/jillashey_uri_edu-ashey_scratch/Mcap2023/trim_round5

mv *count_matrix.csv /work/pi_hputnam_uri_edu/jillashey/Mcap_2023/polyA/output/
```




