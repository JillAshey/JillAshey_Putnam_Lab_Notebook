Data is [here](https://owl.fish.washington.edu/nightingales/E5-coral-DevTimeseries2023/30-1155978938/) and `/project/pi_hputnam_uri_edu/raw_sequencing_data/20250624_Mcap_2023_smRNA` on unity. 



In the scripts folder: `nano raw_qc.sh`

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
#SBATCH -D /work/pi_hputnam_uri_edu/jillashey/Mcap_2023/smRNA/scripts 

# load modules needed
module load parallel/20240822
module load fastqc/0.12.1
module load uri/main
module load all/MultiQC/1.12-foss-2021b

cd /project/pi_hputnam_uri_edu/raw_sequencing_data/20250624_Mcap_2023_smRNA

echo "Running md5sum" $(date)
md5sum *fastq.gz > checkmd5.txt
md5sum -c checkmd5.txt 

# Create an array of fastq files to process
files=($('ls' *001.fastq.gz)) 

# Run fastqc in parallel
echo "Starting fastqc..." $(date)
parallel -j 20 "fastqc {} -o /work/pi_hputnam_uri_edu/jillashey/Mcap_2023/smRNA/output/raw && echo 'Processed {}'" ::: "${files[@]}"
echo "fastQC done." $(date)

cd /work/pi_hputnam_uri_edu/jillashey/Mcap_2023/smRNA/output/raw

echo "Starting multiqc..." $(date)
multiqc *

echo "Initial QC of smRNA data complete." $(date)
```

Submitted batch job 38607468

Time to trim using fastp! In the scripts folder: `nano fastp_trim.sh`

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
#SBATCH -D /work/pi_hputnam_uri_edu/jillashey/Mcap_2023/smRNA/scripts 

# load modules needed
module load parallel/20240822
module load fastqc/0.12.1
module load uri/main
module load all/MultiQC/1.12-foss-2021b


```






