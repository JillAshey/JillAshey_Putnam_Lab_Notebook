---
layout: post
title: e5 closest
date: '2024-08-19'
categories: Analysis
tags: [Bioinformatics]
projects: e5 deep dive
---

## e5 deep dive - extending genes to identify 3' UTRs 

For the e5 deep dive project, we want to identify the genes that the miRNAs bind to. Mature miRNAs usually bind to the 3' UTR portion of an mRNA but our coral gtfs/gffs do not have 3' UTRs annotated. To annotate them, I'm going to use [GeneExt](https://github.com/sebepedroslab/GeneExt), which is a program that extends the genes to obtain the 3' UTR information. [Zoe](https://github.com/zdellaert/LaserCoral/blob/main/code/RNA-seq-bioinf.md) the gene extending queen, has gotten gene ext on unity (I was having a lot of trouble installing it on Andromeda), so I am going to run gene ext on [Unity](https://ood.unity.rc.umass.edu/pun/sys/dashboard/batch_connect/sys/bc_desktop/session_contexts/new). 


### Apulchra 

I am going to test gene ext on Apulchra first. I downloaded the [Apul bam files](https://gannet.fish.washington.edu/seashell/bu-github/deep-dive-expression/D-Apul/output/07.2-Apul-Hisat/) from gannet -- SR alined the deep dive reads to the Apul genome. 

Make folders on unity 

```
cd /project/pi_hputnam_uri_edu/
mkdir jillashey
cd jillashey 
mkdir e5_deepdive
cd e5_deepdive
mkdir data scripts output 
```

Download the bam files from gannet to unity. 

```
cd data 
wget https://gannet.fish.washington.edu/seashell/bu-github/deep-dive-expression/D-Apul/output/07.2-Apul-Hisat/RNA-ACR-140.sorted.bam
https://gannet.fish.washington.edu/seashell/bu-github/deep-dive-expression/D-Apul/output/07.2-Apul-Hisat/RNA-ACR-145.sorted.bam
https://gannet.fish.washington.edu/seashell/bu-github/deep-dive-expression/D-Apul/output/07.2-Apul-Hisat/RNA-ACR-150.sorted.bam
https://gannet.fish.washington.edu/seashell/bu-github/deep-dive-expression/D-Apul/output/07.2-Apul-Hisat/RNA-ACR-173.sorted.bam
https://gannet.fish.washington.edu/seashell/bu-github/deep-dive-expression/D-Apul/output/07.2-Apul-Hisat/RNA-ACR-178.sorted.bam
```

In the scripts folder: `nano merge_bam.sh`

```
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --mem=100GB
#SBATCH --export=NONE
#SBATCH --error="%x_error.%j" #write out slurm error reports
#SBATCH --output="%x_output.%j" #write out any program outpus
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH -D /project/pi_hputnam_uri_edu/jillashey/e5_deepdive/scripts #set working directory
#SBATCH --constraint=avx512

#load modules
module load uri/main
module load SAMtools/1.16.1-GCC-11.3.0

cd /project/pi_hputnam_uri_edu/jillashey/e5_deepdive/data

#use samtools merge to merge all the files
samtools merge merge.bam *.bam
```

Submitted batch job 25761547. The resulting file, `merge.bam`, should be as large or larger than the sum of the individual bam files. Once this is done running, run gene ext! 

In the scripts folder: `nano GeneExt.sh`

```
#!/bin/bash
#SBATCH -t 18:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=500GB
#SBATCH --export=NONE
#SBATCH --error=../scripts/outs_errs/“%x_error.%j” #write out slurm error reports
#SBATCH --output=../scripts/outs_errs/“%x_output.%j” #write out any program outpus
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH -D /project/pi_hputnam_uri_edu/jillashey/e5_deepdive/scripts #set working directory

# activate environment
module load conda/latest #load miniconda
conda activate /work/pi_hputnam_uri_edu/conda/envs/GeneExt/geneext

echo "Environment activated, run Gene ext" $(date)

# use --clip_strand both to not allow GeneExt to create overlaps on the same strand
python /work/pi_hputnam_uri_edu/conda/envs/GeneExt/geneext.py --verbose=3 \
    -g /project/pi_hputnam_uri_edu/jillashey/e5_deepdive/data/Apulchra-genome.gtf \
    -b /project/pi_hputnam_uri_edu/jillashey/e5_deepdive/data/merge.bam \
    -o /project/pi_hputnam_uri_edu/jillashey/e5_deepdive/output/Apul_GeneExt.gtf \
    -j $SLURM_CPUS_PER_TASK --clip_strand both --force
    
echo "Gene Ext complete, deactivate env" $(date)
conda deactivate
```

Submitted batch job 25762636. I got output but not what I thought...There is not gtf file in the output folder, only a log file (`Apul_GeneExt.gtf.GeneExt.log`. The log file has this info at the end: 

```
FUN_044369-T1
FUN_044370-T1
1 transcripts found.
FUN_044370-T1
FUN_044371-T1
1 transcripts found.
FUN_044371-T1
Fix done, annotation with gene features: tmp/genome.fixed.gtf
Added gene feature names. New file name: tmp/genome.fixed.gtf
```

Maybe it needed the gene features added before it could extend the genes? In the tmp folder:

```
wc -l genome.fixed.gtf 
298279 genome.fixed.gtf

head genome.fixed.gtf 
ntLink_0        funannotate     gene    1105    7056    .       +       .       gene_id "FUN_000001"
ntLink_0        funannotate     transcript      1105    7056    .       +       .       transcript_id "FUN_000001-T1"; gene_id "FUN_000001";
ntLink_0        funannotate     exon    1105    1188    .       +       .       transcript_id "FUN_000001-T1"; gene_id "FUN_000001";
ntLink_0        funannotate     exon    1861    1941    .       +       .       transcript_id "FUN_000001-T1"; gene_id "FUN_000001";
ntLink_0        funannotate     exon    2762    2839    .       +       .       transcript_id "FUN_000001-T1"; gene_id "FUN_000001";
ntLink_0        funannotate     exon    5044    7056    .       +       .       transcript_id "FUN_000001-T1"; gene_id "FUN_000001";
ntLink_0        funannotate     gene    10215   15286   .       +       .       gene_id "FUN_000002"
ntLink_0        funannotate     transcript      10215   15286   .       +       .       transcript_id "FUN_000002-T1"; gene_id "FUN_000002";
ntLink_0        funannotate     exon    13074   14383   .       +       .       transcript_id "FUN_000002-T1"; gene_id "FUN_000002";
ntLink_0        funannotate     exon    14722   14900   .       +       .       transcript_id "FUN_000002-T1"; gene_id "FUN_000002";
```

Compared to the actual gtf: 

```
wc -l Apulchra-genome.gtf 
455521 Apulchra-genome.gtf

head Apulchra-genome.gtf 
ntLink_0        funannotate     transcript      1105    7056    .       +       .       transcript_id "FUN_000001-T1"; gene_id "FUN_000001"
ntLink_0        funannotate     exon    1105    1188    .       +       .       transcript_id "FUN_000001-T1"; gene_id "FUN_000001";
ntLink_0        funannotate     exon    1861    1941    .       +       .       transcript_id "FUN_000001-T1"; gene_id "FUN_000001";
ntLink_0        funannotate     exon    2762    2839    .       +       .       transcript_id "FUN_000001-T1"; gene_id "FUN_000001";
ntLink_0        funannotate     exon    5044    7056    .       +       .       transcript_id "FUN_000001-T1"; gene_id "FUN_000001";
ntLink_0        funannotate     CDS     1105    1188    .       +       0       transcript_id "FUN_000001-T1"; gene_id "FUN_000001";
ntLink_0        funannotate     CDS     1861    1941    .       +       0       transcript_id "FUN_000001-T1"; gene_id "FUN_000001";
ntLink_0        funannotate     CDS     2762    2839    .       +       0       transcript_id "FUN_000001-T1"; gene_id "FUN_000001";
ntLink_0        funannotate     CDS     5044    7056    .       +       0       transcript_id "FUN_000001-T1"; gene_id "FUN_000001";
ntLink_0        funannotate     transcript      10215   15286   .       +       .       transcript_id "FUN_000002-T1"; gene_id "FUN_000002"
```

It looks like gene ext added gene rows and removed CDSs from the file? Am I supposed to use this gtf to run gene ext? Maybe. Going to move it from the tmp directory into the data directory in case gene ext deletes the tmp directory. 

```
mv genome.fixed.gtf ../../data/
```

Edit the `Gene_Ext.sh` script so that the `-g` option is directed to the fixed gtf. Submitted batch job 25767632. Success! Downloaded the output file (`Apul_GeneExt.gtf`) to my computer. This is what the file looks like: 

```
ntLink_7	funannotate	gene	79	5033	.	+	.	gene_id "FUN_002303";
ntLink_7	GeneExt	transcript	79	5033	.	+	.	transcript_id "GeneExt~FUN_002303-T1"; gene_id "FUN_002303"; three_prime_ext "354";
ntLink_7	GeneExt	exon	79	179	.	+	.	transcript_id "GeneExt~FUN_002303-T1"; gene_id "FUN_002303";
ntLink_7	GeneExt	exon	1098	1312	.	+	.	transcript_id "GeneExt~FUN_002303-T1"; gene_id "FUN_002303";
ntLink_7	GeneExt	exon	2302	2608	.	+	.	transcript_id "GeneExt~FUN_002303-T1"; gene_id "FUN_002303";
ntLink_7	GeneExt	exon	3242	3337	.	+	.	transcript_id "GeneExt~FUN_002303-T1"; gene_id "FUN_002303";
ntLink_7	GeneExt	exon	3545	3678	.	+	.	transcript_id "GeneExt~FUN_002303-T1"; gene_id "FUN_002303";
ntLink_7	GeneExt	exon	4187	5033	.	+	.	transcript_id "GeneExt~FUN_002303-T1"; gene_id "FUN_002303"; three_prime_ext "354";
ntLink_7	funannotate	gene	12385	16904	.	-	.	gene_id "FUN_002304";
ntLink_7	funannotate	transcript	12385	16904	.	-	.	transcript_id "FUN_002304-T1"; gene_id "FUN_002304";
ntLink_7	funannotate	exon	12385	13137	.	-	.	transcript_id "FUN_002304-T1"; gene_id "FUN_002304";
ntLink_7	funannotate	exon	13624	14387	.	-	.	transcript_id "FUN_002304-T1"; gene_id "FUN_002304";
ntLink_7	funannotate	exon	16898	16904	.	-	.	transcript_id "FUN_002304-T1"; gene_id "FUN_002304";
ntLink_7	funannotate	gene	18480	24541	.	+	.	gene_id "FUN_002305";
ntLink_7	GeneExt	transcript	18480	24541	.	+	.	transcript_id "GeneExt~FUN_002305-T1"; gene_id "FUN_002305"; three_prime_ext "354";
ntLink_7	GeneExt	exon	18480	19242	.	+	.	transcript_id "GeneExt~FUN_002305-T1"; gene_id "FUN_002305";
ntLink_7	GeneExt	exon	19586	19686	.	+	.	transcript_id "GeneExt~FUN_002305-T1"; gene_id "FUN_002305";
```

Not all genes got a 3' UTR extension from gene ext. From the log file: 

```
╭───────────╮
│ All done! │
╰───────────╯
Extended 8177/44371 genes
Median extension length: 1010.0 bp
Running:
	Rscript geneext/plot_extensions.R tmp/extensions.tsv /project/pi_hputnam_uri_edu/jillashey/e5_deepdive/output/Apul_GeneExt.gtf.extension_length.pdf
Running:
	Rscript geneext/peak_density.R tmp/genic_peaks.bed tmp/allpeaks_noov.bed /project/pi_hputnam_uri_edu/jillashey/e5_deepdive/output/Apul_GeneExt.gtf.peak_coverage.pdf 25 
Removing tmp/plus.bam
Removing tmp/minus.bam
Removing tmp/_peaks_tmp
Removing tmp/_genes_tmp
Removing tmp/_peaks_tmp_sorted
Removing tmp/_genes_tmp_sorted
Removing tmp/_genes_peaks_closest
```

Median extension length is around 1000bp, which is what I estimated before when I was extending the genes by 1kb. Only 8177 out of 44371 genes were extended (~18%). Look at the difference in one of the genes (FUN_002322) in the `genome.fixed.gtf` and `Apul_GeneExt.gtf`:

```
# genome.fixed.gtf 
ntLink_7	funannotate	gene	155009	160717	.	+	.	gene_id "FUN_002322"
ntLink_7	funannotate	transcript	155009	160717	.	+	.	transcript_id "FUN_002322-T1"; gene_id "FUN_002322";
ntLink_7	funannotate	exon	155009	155771	.	+	.	transcript_id "FUN_002322-T1"; gene_id "FUN_002322";
ntLink_7	funannotate	exon	156115	156215	.	+	.	transcript_id "FUN_002322-T1"; gene_id "FUN_002322";
ntLink_7	funannotate	exon	157136	157350	.	+	.	transcript_id "FUN_002322-T1"; gene_id "FUN_002322";
ntLink_7	funannotate	exon	158340	158646	.	+	.	transcript_id "FUN_002322-T1"; gene_id "FUN_002322";
ntLink_7	funannotate	exon	159280	159375	.	+	.	transcript_id "FUN_002322-T1"; gene_id "FUN_002322";
ntLink_7	funannotate	exon	159583	159716	.	+	.	transcript_id "FUN_002322-T1"; gene_id "FUN_002322";
ntLink_7	funannotate	exon	160225	160717	.	+	.	transcript_id "FUN_002322-T1"; gene_id "FUN_002322";

# Apul_GeneExt.gtf
ntLink_7	funannotate	gene	155009	161071	.	+	.	gene_id "FUN_002322";
ntLink_7	GeneExt	transcript	155009	161071	.	+	.	transcript_id "GeneExt~FUN_002322-T1"; gene_id "FUN_002322"; three_prime_ext "354";
ntLink_7	GeneExt	exon	155009	155771	.	+	.	transcript_id "GeneExt~FUN_002322-T1"; gene_id "FUN_002322";
ntLink_7	GeneExt	exon	156115	156215	.	+	.	transcript_id "GeneExt~FUN_002322-T1"; gene_id "FUN_002322";
ntLink_7	GeneExt	exon	157136	157350	.	+	.	transcript_id "GeneExt~FUN_002322-T1"; gene_id "FUN_002322";
ntLink_7	GeneExt	exon	158340	158646	.	+	.	transcript_id "GeneExt~FUN_002322-T1"; gene_id "FUN_002322";
ntLink_7	GeneExt	exon	159280	159375	.	+	.	transcript_id "GeneExt~FUN_002322-T1"; gene_id "FUN_002322";
ntLink_7	GeneExt	exon	159583	159716	.	+	.	transcript_id "GeneExt~FUN_002322-T1"; gene_id "FUN_002322";
ntLink_7	GeneExt	exon	160225	161071	.	+	.	transcript_id "GeneExt~FUN_002322-T1"; gene_id "FUN_002322"; three_prime_ext "354";
```

In the new gtf, the gene, transcript and the last exon on the + strand are extended by 354. Only 8177 genes got extended but we need to extend to rest of them for miRNA target prediction. Upload `genome.fixed.gtf` and `Apul_GeneExt.gtf` to Andromeda (this is where miranda, the target prediction software lives). I put them here: `/data/putnamlab/jillashey/e5/refs`. 

Filter so that there is only 3'UTR information 

```
interactive 

# Input files
ORIGINAL_GTF="genome.fixed.gtf"
EXTENDED_GTF="Apul_GeneExt.gtf"

# Output file
OUTPUT_GTF="genome_with_UTR.gtf"

# Temporary file
TEMP_FILE="temp_utr.gtf"

# Function to extract UTR information
extract_utr_info() {
    grep 'three_prime_ext' "$EXTENDED_GTF" | while read -r line; do
        gene_id=$(echo "$line" | awk -F 'gene_id "' '{print $2}' | awk -F '"' '{print $1}')
        utr_length=$(echo "$line" | awk -F 'three_prime_ext "' '{print $2}' | awk -F '"' '{print $1}')
        strand=$(echo "$line" | awk '{print $7}')
        echo "$gene_id $utr_length $strand"
    done
}

# Extract UTR information
extract_utr_info > "$TEMP_FILE"

# Process the original GTF and add UTR information
awk -v temp_file="$TEMP_FILE" '
BEGIN {
    while ((getline < temp_file) > 0) {
        split($0, a, " ")
        utr_length[a[1]] = a[2]
        utr_strand[a[1]] = a[3]
    }
    close(temp_file)
}
{
    print $0
    if ($3 == "gene") {
        gene_id = $0
        sub(/.*gene_id "/, "", gene_id)
        sub(/".*/, "", gene_id)
        if (gene_id in utr_length) {
            strand = $7
            if (strand == "+") {
                utr_start = $5 + 1
                utr_end = $5 + utr_length[gene_id]
            } else if (strand == "-") {
                utr_end = $4 - 1
                utr_start = $4 - utr_length[gene_id]
            }
            if (utr_start > 0 && utr_end > 0) {
                printf "%s\tfunannotate\tthree_prime_UTR\t%d\t%d\t.\t%s\t.\ttranscript_id \"%s-T1\"; gene_id \"%s\";\n", 
                       $1, utr_start, utr_end, strand, gene_id, gene_id
            }
        }
    }
}' "$ORIGINAL_GTF" > "$OUTPUT_GTF"

# Clean up
rm "$TEMP_FILE"

echo "Processing complete. Output written to $OUTPUT_GTF"

head genome_with_UTR.gtf 
ntLink_7	funannotate	transcript	79	4679	.	+	.	transcript_id "FUN_002303-T1"; gene_id "FUN_002303";
ntLink_7	funannotate	gene	79	4679	.	+	.	gene_id "FUN_002303"
ntLink_7	funannotate	three_prime_UTR	4680	5033	.	+	.	transcript_id "FUN_002303-T1"; gene_id "FUN_002303";
ntLink_7	funannotate	exon	79	179	.	+	.	transcript_id "FUN_002303-T1"; gene_id "FUN_002303";
ntLink_7	funannotate	exon	1098	1312	.	+	.	transcript_id "FUN_002303-T1"; gene_id "FUN_002303";
ntLink_7	funannotate	exon	2302	2608	.	+	.	transcript_id "FUN_002303-T1"; gene_id "FUN_002303";
ntLink_7	funannotate	exon	3242	3337	.	+	.	transcript_id "FUN_002303-T1"; gene_id "FUN_002303";
ntLink_7	funannotate	exon	3545	3678	.	+	.	transcript_id "FUN_002303-T1"; gene_id "FUN_002303";
ntLink_7	funannotate	exon	4187	4679	.	+	.	transcript_id "FUN_002303-T1"; gene_id "FUN_002303";
ntLink_7	funannotate	gene	12385	16904	.	-	.	gene_id "FUN_002304"
```

Success! Also checked to make sure strandedness was taken into account. Keep only genes that do not have an associated three_prime_UTR. 

```
INPUT_GTF="genome_with_UTR.gtf"
OUTPUT_GTF="genome_without_UTR_genes.gtf"

# Create a temporary file to store gene IDs with 3' UTRs
TEMP_FILE="temp_genes_with_utr.txt"

# Find all gene IDs that have an associated 3' UTR
awk '$3 == "three_prime_UTR" {print $0}' "$INPUT_GTF" | \
    sed -n 's/.*gene_id "\([^"]*\)".*/\1/p' | sort | uniq > "$TEMP_FILE"

# Process the GTF file, excluding genes with 3' UTRs and their associated features
awk -v temp_file="$TEMP_FILE" '
BEGIN {
    while ((getline < temp_file) > 0) {
        genes_with_utr[$0] = 1
    }
    close(temp_file)
}
{
    gene_id = $0
    sub(/.*gene_id "/, "", gene_id)
    sub(/".*/, "", gene_id)
    
    if (!(gene_id in genes_with_utr)) {
        print $0
    }
}' "$INPUT_GTF" > "$OUTPUT_GTF"

# Remove the temporary file
rm "$TEMP_FILE"
```

Success!

```
head genome_without_UTR_genes.gtf
ntLink_7	funannotate	gene	12385	16904	.	-	.	gene_id "FUN_002304"
ntLink_7	funannotate	transcript	12385	16904	.	-	.	transcript_id "FUN_002304-T1"; gene_id "FUN_002304";
ntLink_7	funannotate	exon	12385	13137	.	-	.	transcript_id "FUN_002304-T1"; gene_id "FUN_002304";
ntLink_7	funannotate	exon	13624	14387	.	-	.	transcript_id "FUN_002304-T1"; gene_id "FUN_002304";
ntLink_7	funannotate	exon	16898	16904	.	-	.	transcript_id "FUN_002304-T1"; gene_id "FUN_002304";
ntLink_7	funannotate	gene	31894	36413	.	-	.	gene_id "FUN_002306"
ntLink_7	funannotate	transcript	31894	36413	.	-	.	transcript_id "FUN_002306-T1"; gene_id "FUN_002306";
ntLink_7	funannotate	exon	31894	32646	.	-	.	transcript_id "FUN_002306-T1"; gene_id "FUN_002306";
ntLink_7	funannotate	exon	33133	33896	.	-	.	transcript_id "FUN_002306-T1"; gene_id "FUN_002306";
ntLink_7	funannotate	exon	36407	36413	.	-	.	transcript_id "FUN_002306-T1"; gene_id "FUN_002306";
```

Filter so only genes remain

```
awk '$3 == "gene" { print }' genome_without_UTR_genes.gtf > gene_without_UTRs.gtf

wc -l gene_without_UTRs.gtf 
36193 gene_without_UTRs.gtf

head gene_without_UTRs.gtf 
ntLink_7	funannotate	gene	12385	16904	.	-	.	gene_id "FUN_002304"
ntLink_7	funannotate	gene	31894	36413	.	-	.	gene_id "FUN_002306"
ntLink_7	funannotate	gene	37989	41453	.	+	.	gene_id "FUN_002307"
ntLink_7	funannotate	gene	51376	55895	.	-	.	gene_id "FUN_002308"
ntLink_7	funannotate	gene	70886	75405	.	-	.	gene_id "FUN_002310"
ntLink_7	funannotate	gene	90396	94918	.	-	.	gene_id "FUN_002312"
ntLink_7	funannotate	gene	96494	97288	.	+	.	gene_id "FUN_002313"
ntLink_7	funannotate	gene	108757	110163	.	+	.	gene_id "FUN_002315"
ntLink_7	funannotate	gene	110424	111490	.	-	.	gene_id "FUN_002316"
ntLink_7	funannotate	gene	116079	118064	.	+	.	gene_id "FUN_002317"
```

There are 36193 genes that do not have 3'UTR annotation. 8177 genes were extended by gene ext. These two numbers combined equal the total number of genes in the genome. 

Sort `gene_without_UTRs.gtf` by chromosome. 

```
module load BEDTools/2.30.0-GCC-11.3.0
sortBed -faidx apul.Chromosome_names.txt -i gene_without_UTRs.gtf > gene_without_UTRs_sorted.gtf
```

File was already sorted but just in case. Extract 1kb 3' UTRs with `bedflank`. 

```
bedtools flank -i gene_without_UTRs_sorted.gtf -g apul.Chromosome_lenghts.txt -l 0 -r 1000 -s | \
awk '{gsub("gene","3prime_UTR",$3); print $0}' | \
awk '{if($5-$4 > 3) print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9" "$10" "$11" "$12}' | \
tr ' ' '\t' > genes.3UTR_1kb.gtf

wc -l genes.3UTR_1kb.gff
36193 genes.3UTR_1kb.gff

head genes.3UTR_1kb.gff
ntLink_7	funannotate	3prime_UTR	11385	12384	.	-	.	gene_id	"FUN_002304"		
ntLink_7	funannotate	3prime_UTR	30894	31893	.	-	.	gene_id	"FUN_002306"		
ntLink_7	funannotate	3prime_UTR	41454	42453	.	+	.	gene_id	"FUN_002307"		
ntLink_7	funannotate	3prime_UTR	50376	51375	.	-	.	gene_id	"FUN_002308"		
ntLink_7	funannotate	3prime_UTR	69886	70885	.	-	.	gene_id	"FUN_002310"		
ntLink_7	funannotate	3prime_UTR	89396	90395	.	-	.	gene_id	"FUN_002312"		
ntLink_7	funannotate	3prime_UTR	97289	98288	.	+	.	gene_id	"FUN_002313"		
ntLink_7	funannotate	3prime_UTR	110164	111163	.	+	.	gene_id	"FUN_002315"		
ntLink_7	funannotate	3prime_UTR	109424	110423	.	-	.	gene_id	"FUN_002316"		
ntLink_7	funannotate	3prime_UTR	118065	119064	.	+	.	gene_id	"FUN_002317"
```

Success. Remove extra tabs

```
awk -F'\t' 'OFS="\t" {$1=$1; print}' genes.3UTR_1kb.gtf > genes.3UTR_1kb_clean.gtf
awk -F'\t' 'OFS="\t" {$1=$1; print}' gene_without_UTRs_sorted.gtf > gene_without_UTRs_sorted_clean.gtf
```

Make sure that gtfs have 9 columns, as per gtf format

```
awk -F'\t' 'OFS="\t" {print $1, $2, $3, $4, $5, $6, $7, $8, $9" "$10}' genes.3UTR_1kb_clean.gtf > genes.3UTR_1kb_fixed.gtf
awk -F'\t' '{print NF}' genes.3UTR_1kb_fixed.gtf | sort | uniq -c
  36193 9

awk -F'\t' '{print NF}' gene_without_UTRs_sorted_clean.gtf | sort | uniq -c
  36193 9
```

Next, subtract portions of UTRs that overlap nearby genes 

```
bedtools subtract -a genes.3UTR_1kb_fixed.gtf -b gene_without_UTRs_sorted_clean.gtf > genes.3UTR_1kb_corrected.gtf

wc -l genes.3UTR_1kb_corrected.gtf
44238 genes.3UTR_1kb_corrected.gtf
```

There are more lines in this file because when the 3' UTR region partially overlaps with a gene, bedtools subtract will split the 3' UTR region into multiple parts, keeping only the non-overlapping portions. Will need to look into which part is the 3'UTR for a specific gene. 

From the Gene Ext output, keep only the 3'UTR information 

```
awk '$3 == "three_prime_UTR"' genome_with_UTR.gtf > UTR_only_GeneExt.gtf

head UTR_only_GeneExt.gtf 
ntLink_7	funannotate	three_prime_UTR	4680	5033	.	+	.	transcript_id "FUN_002303-T1"; gene_id "FUN_002303";
ntLink_7	funannotate	three_prime_UTR	24188	24541	.	+	.	transcript_id "FUN_002305-T1"; gene_id "FUN_002305";
ntLink_7	funannotate	three_prime_UTR	63180	63540	.	+	.	transcript_id "FUN_002309-T1"; gene_id "FUN_002309";
ntLink_7	funannotate	three_prime_UTR	82690	83043	.	+	.	transcript_id "FUN_002311-T1"; gene_id "FUN_002311";
ntLink_7	funannotate	three_prime_UTR	102203	102556	.	+	.	transcript_id "FUN_002314-T1"; gene_id "FUN_002314";
ntLink_7	funannotate	three_prime_UTR	121583	122060	.	+	.	transcript_id "FUN_002318-T1"; gene_id "FUN_002318";
ntLink_7	funannotate	three_prime_UTR	141305	141567	.	+	.	transcript_id "FUN_002320-T1"; gene_id "FUN_002320";
ntLink_7	funannotate	three_prime_UTR	160718	161071	.	+	.	transcript_id "FUN_002322-T1"; gene_id "FUN_002322";
ntLink_7	funannotate	three_prime_UTR	180226	180589	.	+	.	transcript_id "FUN_002326-T1"; gene_id "FUN_002326";
ntLink_8	funannotate	three_prime_UTR	59307	60891	.	+	.	transcript_id "FUN_002328-T1"; gene_id "FUN_002328";
```

Join the two 3'UTR files - one that was created by Gene ext (`UTR_only_GeneExt.gtf`) and one that was manually created by adding 1kb right flanks (`genes.3UTR_1kb_corrected.gtf`). 

```
cat UTR_only_GeneExt.gtf genes.3UTR_1kb_corrected.gtf > Apul_all_3UTRs.gtf
```

Extract 3'UTR sequences from genome fasta 

```
bedtools getfasta -fi /data/putnamlab/REFS/Apul/apul.hifiasm.s55_pa.p_ctg.fa.k32.w100.z1000.ntLink.5rounds.masked.fa -bed Apul_all_3UTRs.gtf -fo Apul_all_3UTRs.fasta -fullHeader

head Apul_all_3UTRs.fasta
>ntLink_7:4679-5033
TCTGACATCGAATGCTGCCCAAAGCAGGGAAAAGCTAAGGACAGGAACTCTGAACTCCGTCAAAAGTTGCGCTAGATCGAGTCTAGGATAGGATAGGAATTTCACCGAACTTTCTATTTGTGTATTTGTATTTCGCTCTTCGAGCCTATTCAGTCAAACGGTTCAAAAGTCCATGTTCATGTAAGATGCATTAAATCAGTTTAGATCTCGATTCTCGTTGGATGTAAATGTAGGTTTAATTAAAAAAAAAAAAAGATTGCGAGGAAACGATCAGTGTTGGTCGATGTGTTGGTGCCCCCGGGGAAACGGGTTTTGAAGGAGAAATTCTACGCTCTAAAGAGCACTTTTATACAG
>ntLink_7:24187-24541
TCTGACATCGAATGCTGCCCAAAGCAGGGAAAAGCTAAGGACAGGAACTCTGAACTCCGTCAAAAGTTGCGCTAGATCGAGTCTAGGATAGGATAGGAATTTCACCGAACTTTCTATTTGTGTATTTGTATTTCGCTCTTCGAGCCTATTCAGTCAAACGGTTCAAAAGTCCATGTTCATGTAAGATGCATTAAATCAGTTTAGATCTCGATTCTCGTTGGATGTAAATGTAGGTTTAATTAAAAAAAAAAAAAGATTGCGAGGAAACGATCAGTGTTGGTCGATGTGTTGGTGCCCCCGGGGAAACGGGTTTTGAAGGAGAAATTCTACGCTCTAAAGAGCACTTTTATACAG
>ntLink_7:63179-63540
TCTGACATCGAATGCTGCCCAAAGCAGGGAAAAGCTAAGGACAGGAACTCTGAACTCCGTCAAAAGTTGCGCTAGATCGAGTCTAGGATAGGATAGGAATTTCACCGAACTTTCTATTTGTGTATTTGTATTTCGCTCTTCGAGCCTATTCAGTCAAACGGTTCAAAAGTCCATGTTCATGTAAGATGCATTAAATCAGTTTAGATCTCGATTCTCGTTGGATGTAAATGTAGGTTTAATTAAAAAAAAAAAAAGATTGCGAGGAAACGATCAGTGTTGGTCGATGTGTTGGTGCCCCCGGGGAAACGGGTTTTGAAGGAGAAATTCTACGCTCTAAAGAGCACTTTTATACAGAACAAGA
>ntLink_7:82689-83043
TCTGACATCGAATGCTGCCCAAAGCAGGGAAAAGCTAAGGACAGGAACTCTGAACTCCGTCAAAAGTTGCGCTAGATCGAGTCTAGGATAGGATAGGAATTTCACCGAACTTTCTATTTGTGTATTTGTATTTCGCTCTTCGAGCCTATTCAGTCAAACGGTTCAAAAGTCCATGTTCATGTAAGATGCATTAAATCAGTTTAGATCTCGATTCTCGTTGGATGTAAATGTAGGTTTAATTAAAAAAAAAAAAAGATTGCGAGGAAACGATCAGTGTTGGTCGATGTGTTGGTGCCCCCGGGGAAACGGGTTTTGAAGGAGAAATTCTACGCTCTAAAGAGCACTTTTATACAG
>ntLink_7:102202-102556
TCTGACATCGAATGCTGCCCAAAGCAGGGAAAAGCTAAGGACAGGAACTCTGAACTCCGTCAAAAGTTGCGCTAGATCGAGTCTAGGATAGGATAGGAATTTCACCGAACTTTCTATTTGTGTATTTGTATTTCGCTCTTCGAGCCTATTCAGTCAAACGGTTCAAAAGTCCATGTTCATGTAAGATGCATTAAATCAGTTTAGATCTCGATTCTCGTTGGATGTAAATGTAGGTTTAATTAAAAAAAAAAAAAGATTGCGAGGAAACGATCAGTGTTGGTCGATGTGTTGGTGCCCCCGGGGAAACGGGTTTTGAAGGAGAAATTCTACGCTCTAAAGAGCACTTTTATACAG

zgrep -c ">" Apul_all_3UTRs.fasta
52416
```

For some reason, bedtools fasta always subtracts 1bp from the start site...I attempted to fix it with some grep and awk but no luck. Also it is not the same number of genes because some UTRs were split if they intersected with a gene. Will come back to this to see which UTRs are appropriate to use for each gene. 

Move all of the Apul stuff into its own folder 

```
cd /data/putnamlab/jillashey/e5/refs
mkdir Apul Pmea Peve
mv Apul* Apul/
mv apul* Apul/
mv UTR_only_GeneExt.gtf Apul/
mv gen* Apul/
```

Run miranda. In the scripts folder: `nano miranda_strict_all_3UTRs_Apul.sh`

```
#!/bin/bash -i
#SBATCH -t 72:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --export=NONE
#SBATCH --mem=500GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/e5/scripts
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

echo "Apul starting miranda run with all genes and miRNAs with score cutoff >100, energy cutoff <-10, and strict binding invoked"$(date)

module load Miniconda3/4.9.2
conda activate /data/putnamlab/conda/miranda 

miranda /data/putnamlab/jillashey/e5/mirna_seqs/Apul_results_mature_named.fasta /data/putnamlab/jillashey/e5/refs/Apul/Apul_all_3UTRs.fasta -sc 100 -en -10 -strict -out /data/putnamlab/jillashey/e5/output/miranda/miranda_strict_all_3UTRs_Apul.tab

conda deactivate

echo "miranda run finished!"$(date)
echo "counting number of putative interactions predicted" $(date)

zgrep -c "Performing Scan" /data/putnamlab/jillashey/e5/output/miranda/miranda_strict_all_3UTRs_Apul.tab

echo "Parsing output" $(date)
grep -A 1 "Scores for this hit:" /data/putnamlab/jillashey/e5/output/miranda/miranda_strict_all_3UTRs_Apul.tab | sort | grep '>' > /data/putnamlab/jillashey/e5/output/miranda/miranda_strict_all_3UTRs_parsed_Apul.txt

echo "counting number of putative interactions predicted" $(date)
wc -l /data/putnamlab/jillashey/e5/output/miranda/miranda_strict_all_3UTRs_parsed_Apul.txt

echo "Apul miranda script complete" $(date)
```

Submitted batch job 345245. Ran in about an hour. 

```
counting number of putative interactions predicted Mon Oct 28 13:11:27 EDT 2024
1991808
Parsing output Mon Oct 28 13:11:32 EDT 2024
counting number of putative interactions predicted Mon Oct 28 13:11:34 EDT 2024
99345 /data/putnamlab/jillashey/e5/output/miranda/miranda_strict_all_3UTRs_parsed_Apul.txt
```

I am looking at some of the coral miRNA papers and they all mention a consecutive 8 bp "seed" region of the miRNA that directly binds with 8 bp of the mRNA. Even though I added the `-strict` flag, there only seems to be 7bp consectively if that makes sense. Example: 

```
Read Sequence:ntLink_8:4298823-4300537 (1714 nt)
=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
Performing Scan: apul-mir-novel-2 vs ntLink_8:4298823-4300537
=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

   Forward:     Score: 167.000000  Q:2 to 21  R:105 to 129 Align Len (20) (80.00%) (80.00%)

   Query:    3' ugagAAGUAAAUAGA-UGAGCAAUu 5'
                    || |  ||||| |||||||| 
   Ref:      5' taaaTTAAAGTATCTGACTCGTTAt 3'

   Energy:  -15.730000 kCal/Mol

Scores for this hit:
>apul-mir-novel-2       ntLink_8:4298823-4300537        167.00  -15.73  2 21    105 129 20      80.00%  80.00%

Score for this Scan:
Seq1,Seq2,Tot Score,Tot Energy,Max Score,Max Energy,Strand,Len1,Len2,Positions
>>apul-mir-novel-2      ntLink_8:4298823-4300537        167.00  -15.73  167.00  -15.73  77      24      1714     105
Complete

Read Sequence:ntLink_8:8821127-8822825 (1698 nt)
=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
Performing Scan: apul-mir-novel-2 vs ntLink_8:8821127-8822825
=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

   Forward:     Score: 152.000000  Q:2 to 15  R:10 to 32 Align Len (13) (84.62%) (92.31%)

   Query:    3' ugagaaguaaAUAGAUGAGCAAUu 5'
                          |||:| ||||||| 
   Ref:      5' agattaagagTATTT-CTCGTTAg 3'

   Energy:  -13.480000 kCal/Mol

Scores for this hit:
>apul-mir-novel-2       ntLink_8:8821127-8822825        152.00  -13.48  2 15    10 32   13      84.62%  92.31%

Score for this Scan:
Seq1,Seq2,Tot Score,Tot Energy,Max Score,Max Energy,Strand,Len1,Len2,Positions
>>apul-mir-novel-2      ntLink_8:8821127-8822825        152.00  -13.48  152.00  -13.48  158     24      1698     10
Complete
```

Going to rerun and set the `-scale` flag to 8. In the (paltry) [miranda information](https://www.animalgenome.org/bioinfo/resources/manuals/miranda.html) online, this flag does the following: Set the scaling parameter to scale. This scaling is applied to match / mismatch scores in the  critical 10bp  region  of  the  5' end of the microRNA. Many known examples of miRNA:Target duplexes are  highly complementary in this region. This parameter can be thought of as a contrast function  to  more  effectively detect alignments of this type. The default is 4 and I'm going to set it to 8. Submitted batch job 345472. Took about 10 hrs to run. 

```
counting number of putative interactions predicted Tue Oct 29 04:06:29 EDT 2024
1991808
Parsing output Tue Oct 29 04:06:36 EDT 2024
counting number of putative interactions predicted Tue Oct 29 04:06:37 EDT 2024
99403 /data/putnamlab/jillashey/e5/output/miranda/miranda_strict_all_3UTRs_8seed_parsed_Apul.txt
```

Interesting that the number of putative interactions actually increased from the previous run (99345 vs 99403). Look at the parsed files: 

```
head miranda_strict_all_3UTRs_parsed_Apul.txt 
>apul-mir-100	ntLink_6:10255357-10256357	145.00	-15.40	2 19	264 286	20	65.00%	70.00%
>apul-mir-100	ntLink_6:10563192-10564192	153.00	-20.42	2 14	791 810	12	83.33%	91.67%
>apul-mir-100	ntLink_6:10617301-10618424	145.00	-15.84	2 16	679 697	14	78.57%	78.57%
>apul-mir-100	ntLink_6:10672091-10673508	161.00	-23.97	2 16	951 969	14	92.86%	92.86%
>apul-mir-100	ntLink_6:11090988-11091647	146.00	-15.95	2 11	545 564	9	88.89%	100.00%
>apul-mir-100	ntLink_6:11405289-11406669	149.00	-16.68	2 16	1000 1018	14	78.57%	85.71%
>apul-mir-100	ntLink_6:1142082-1142891	142.00	-16.61	2 11	9 28	9	88.89%	88.89%
>apul-mir-100	ntLink_6:11697492-11698492	154.00	-18.81	2 15	540 559	13	84.62%	84.62%
>apul-mir-100	ntLink_6:12100733-12101733	154.00	-15.72	2 15	563 582	13	84.62%	84.62%
>apul-mir-100	ntLink_6:12100733-12101733	157.00	-23.22	2 16	467 485	14	85.71%	92.86%

head miranda_strict_all_3UTRs_8seed_parsed_Apul.txt
>apul-mir-100	ntLink_6:10255357-10256357	285.00	-15.40	2 19	264 286	20	65.00%	70.00%
>apul-mir-100	ntLink_6:10563192-10564192	293.00	-20.42	2 14	791 810	12	83.33%	91.67%
>apul-mir-100	ntLink_6:10617301-10618424	285.00	-15.84	2 16	679 697	14	78.57%	78.57%
>apul-mir-100	ntLink_6:10672091-10673508	301.00	-23.97	2 16	951 969	14	92.86%	92.86%
>apul-mir-100	ntLink_6:11090988-11091647	286.00	-15.95	2 11	545 564	9	88.89%	100.00%
>apul-mir-100	ntLink_6:11405289-11406669	289.00	-16.68	2 16	1000 1018	14	78.57%	85.71%
>apul-mir-100	ntLink_6:1142082-1142891	282.00	-16.61	2 11	9 28	9	88.89%	88.89%
>apul-mir-100	ntLink_6:11697492-11698492	294.00	-18.81	2 15	540 559	13	84.62%	84.62%
>apul-mir-100	ntLink_6:12100733-12101733	294.00	-15.72	2 15	563 582	13	84.62%	84.62%
>apul-mir-100	ntLink_6:12100733-12101733	297.00	-23.22	2 16	467 485	14	85.71%	92.86%
```

The interactions appear to be the same except the score is much higher, maybe overly inflated? Are the scans the same? 

```
# from original run
Read Sequence:ntLink_8:8646169-8647749 (1580 nt)
=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
Performing Scan: apul-mir-novel-2 vs ntLink_8:8646169-8647749
=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

   Forward:     Score: 152.000000  Q:2 to 15  R:386 to 408 Align Len (13) (84.62%) (92.31%)

   Query:    3' ugagaaguaaAUAGAUGAGCAAUu 5'
                          |||:| ||||||| 
   Ref:      5' agattaagagTATTT-CTCGTTAg 3'

   Energy:  -13.480000 kCal/Mol

Scores for this hit:
>apul-mir-novel-2       ntLink_8:8646169-8647749        152.00  -13.48  2 15    386 408 13      84.62%  92.31%

Score for this Scan:
Seq1,Seq2,Tot Score,Tot Energy,Max Score,Max Energy,Strand,Len1,Len2,Positions
>>apul-mir-novel-2      ntLink_8:8646169-8647749        152.00  -13.48  152.00  -13.48  156     24      1580     386
Complete

# from scale 8 run 
Read Sequence:ntLink_8:8646169-8647749 (1580 nt)
=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
Performing Scan: apul-mir-novel-2 vs ntLink_8:8646169-8647749
=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

   Forward:     Score: 292.000000  Q:2 to 15  R:386 to 408 Align Len (13) (84.62%) (92.31%)

   Query:    3' ugagaaguaaAUAGAUGAGCAAUu 5'
                          |||:| ||||||| 
   Ref:      5' agattaagagTATTT-CTCGTTAg 3'

   Energy:  -13.480000 kCal/Mol

Scores for this hit:
>apul-mir-novel-2       ntLink_8:8646169-8647749        292.00  -13.48  2 15    386 408 13      84.62%  92.31%

Score for this Scan:
Seq1,Seq2,Tot Score,Tot Energy,Max Score,Max Energy,Strand,Len1,Len2,Positions
>>apul-mir-novel-2      ntLink_8:8646169-8647749        292.00  -13.48  292.00  -13.48  156     24      1580     386
Complete
```

The only difference is the score so I think I am inflating the calculations somehow, would need to look at source code. After looking at other papers that have used miranda, they typically employ strict binding, energy <-10, default scale (4), and 140-150 score (see example from human flu hozt response to virus [here](https://www.frontiersin.org/journals/molecular-biosciences/articles/10.3389/fmolb.2022.866072/full)). [Gajigan & Conaco (2017)](https://onlinelibrary.wiley.com/doi/10.1111/mec.14130) used miranda on an Acropora spp as well and they used strict seed binding and energy <-10. They further narrowed targets with exact seed match (is seed match mean 6, 7, 8bp?) and an A in position 1 (of miRNA or mRNA?); they cited a mammalian paper here for these choices specifically (Agarwal, Bell, Nam & Bartel 2015). They also used `fasta` to check the targets for more extensive complementarity and scored the alignments as described in Moran et al. 2014 (classic nematostella miRNA paper). Moran et al. 2014 did the following: "mapped sequences with FASTA v36 (Pearson and Lipman 1988) using the parameters -n -H -Q -f -16 -r +15/-10 -g -10 -w 100 -W 25 -E 100000 -i -U and scored the alignments using a weighted sum of the number of mismatches (scoremismatch = 1, scoreG:U = 0.5) with mismatches for guide RNA nucleotides g2–g13 counting double (scoremismatch = 2, scoreG:U = 1)." Moran et al. got this scoring method from [Fahlgren et al. 2007](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0000219#s3), who did the following: "miRNA targets were computationally predicted as described [34]. Briefly, potential targets from FASTA searches (+15/−10 match/mismatch scoring ratio, -16 gap penalty and a RNA scoring matrix) were scored using a position-dependent, mispair penalty system [34]. Penalties were assessed for mismatches, bulges, and gaps (+1 per position) and G∶U pairs (+0.5 per position). Penalties were doubled if the mismatch, bulge, gap, or G∶U pair occurred at positions 2 to 13 relative to the 5′ end of the miRNA. Only one single-nt bulge or single-nt gap was allowed. Based on a reference set of validated miRNA targets, only predicted targets with scores of four or less were considered reasonable." Need to look more into this scoring. 

I am going to run miranda once more with the following flags: `-en -10 -strict`. The score (140), scale (4), gap-open penalty (-4), and gap-extend pentaly (-9) will be kept as the defaults. Submitted batch job 345844





things to ask e5 ppl: 

- do we want to run gene ext or estimate the 3'UTR as 1000bp right flank? if yes, i need bam files from Peve and Pmea 

redo the 1kb stuff 

I want to check to see how similar the results would be between the 1kb estimate and the gene ext + 1kb estimate. Use this [post](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/_posts/2024-06-15-e5-deepdive-miRNA-TargetPrediction.md) as reference. 

```
cd /data/putnamlab/jillashey/e5/refs/Apul

grep $'\tgene\t' genome.fixed.gtf > Apul_gene.gtf
wc -l Apul_gene.gtf 
44371 Apul_gene.gtf
```

Sort by chromosome 

```
module load BEDTools/2.30.0-GCC-11.3.0

sortBed -faidx apul.Chromosome_names.txt -i Apul_gene.gtf > Apul_gene_sorted.gtf
```

Extract 1kb 3'UTRs and keep associated gene info

```
bedtools flank -i Apul_gene_sorted.gtf -g apul.Chromosome_lenghts.txt -l 0 -r 1000 -s | \
awk 'BEGIN{OFS="\t"} {
    gsub("gene","3prime_UTR",$3);
    if($5-$4 > 3) {
        gene_id = $NF;
        sub(/;$/, "", gene_id);
        print $1, $2, $3, $4, $5, $6, $7, $8, gene_id;
    }
}' > Apul_3UTR_1kb.gtf

```

Subtract overlaps of other genes

```
bedtools subtract -a Apul_3UTR_1kb.gtf -b Apul_gene_sorted.gtf > Apul_3UTR_1kb_corrected.gtf
```

Extract 3'UTR seqs from genome fasta 

```
bedtools getfasta -fi /data/putnamlab/REFS/Apul/apul.hifiasm.s55_pa.p_ctg.fa.k32.w100.z1000.ntLink.5rounds.masked.fa -bed Apul_3UTR_1kb_corrected.gtf -fo Apul_3UTR_1kb.fasta -fullHeader
```

Run miranda - edit `miranda_strict_all_1kb_apul.sh` in scripts folder so that input 3'UTR fasta is `/data/putnamlab/jillashey/e5/refs/Apul/Apul_3UTR_1kb.fasta`. Submitted batch job 345851

