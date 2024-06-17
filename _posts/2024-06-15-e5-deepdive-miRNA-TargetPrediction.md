On my local computer, make fasta file with the given miRNA names. 

```
cd /Users/jillashey/Desktop/PutnamLab/Repositories/deep-dive/DEF-cross-species/output/10-shortRNA-ShortStack-comparison

awk -F, 'NR>1 {print ">" $22 "\n" $11}' Apul_results_mature_named.csv > Apul_results_mature_named.fasta

awk -F, 'NR>1 {print ">" $22 "\n" $11}' Peve_results_mature_named.csv > Peve_results_mature_named.fasta

awk -F, 'NR>1 {print ">" $22 "\n" $11}' Ptuh_results_mature_named.csv > Ptuh_results_mature_named.fasta
```

Copy the fasta files to the server 

```
cd /data/putnamlab/jillashey/e5
mkdir mirna_seqs
cd mirna_seqs
```

### Apul

To run miranda, I need to identify the 3' UTR ends. 

First, identify counts of each feature from gff file 

```
GFF_FILE="/data/putnamlab/jillashey/genome/Amil_v2.01/Amil.all.maker.noseq.gff"
genome="/data/putnamlab/jillashey/genome/Amil_v2.01/Amil.v2.01.chrs.fasta"

grep -v '^#' ${GFF_FILE} | cut -s -f 3 | sort | uniq -c | sort -rn > all_features.txt

1572090 match_part
 558048 match
 201232 exon
 176923 CDS
  88262 expressed_sequence_match
  84724 protein_match
  38247 gene
  28188 mRNA
  19507 three_prime_UTR
  15776 five_prime_UTR
  10059 tRNA
    854 contig
```

There are a decent amount of three prime UTRs annotated in the Amil gff. Calculate the length of the 3' UTRs. 

```
awk '$3 == "three_prime_UTR" {print $0, $5 - $4 + 1}' OFS="\t" Amil.all.maker.noseq.gff > three_prime_UTR_lengths.txt

awk '{sum += $NF} END {if (NR > 0) print sum / NR}' three_prime_UTR_lengths.txt
```

Average length of 3' UTRs is 424.277. Extract feature types and generate individual gffs for each feature. 

```
grep $'\texon\t' ${GFF_FILE} > amil_GFFannotation.exon.gff
grep $'\tCDS\t' ${GFF_FILE} > amil_GFFannotation.CDS.gff
grep $'\t expressed_sequence_match\t' ${GFF_FILE} > amil_GFFannotation.expressed_sequence_match.gff
grep $'\tprotein_match\t' ${GFF_FILE} > amil_GFFannotation.protein_match.gff
grep $'\tgene\t' ${GFF_FILE} > amil_GFFannotation.gene.gff
grep $'\tmRNA\t' ${GFF_FILE} > amil_GFFannotation.mRNA.gff
grep $'\tthree_prime_UTR\t' ${GFF_FILE} > amil_GFFannotation.three_prime_UTR.gff
grep $'\tfive_prime_UTR\t' ${GFF_FILE} > amil_GFFannotation.five_prime_UTR.gff
grep $'\ttRNA\t' ${GFF_FILE} > amil_GFFannotation.tRNA.gff
grep $'\tcontig\t' ${GFF_FILE} > amil_GFFannotation.contig.gff
```

Extract chromosome lengths 

```
cat is ${genome} | awk '$0 ~ ">" {if (NR > 1) {print c;} c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' > amil.Chromosome_lenghts.txt

wc -l amil.Chromosome_lenghts.txt 
854 amil.Chromosome_lenghts.txt
```

Extract scaffold names 

```
awk -F" " '{print $1}' amil.Chromosome_lenghts.txt > amil.Chromosome_names.txt
```

Sort gffs by chromosome name. In the the e5 scripts folder: `nano bed_sort.sh`

```
#!/bin/bash 
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=5
#SBATCH --export=NONE
#SBATCH --mem=250GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/e5/scripts
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

module load BEDTools/2.30.0-GCC-11.3.0

cd /data/putnamlab/jillashey/genome/Amil_v2.01

echo "Sorting gffs by chromosome" $(date)

sortBed -faidx amil.Chromosome_names.txt -i amil_GFFannotation.exon.gff > amil_GFFannotation.exon_sorted.gff
sortBed -faidx amil.Chromosome_names.txt -i amil_GFFannotation.CDS.gff > amil_GFFannotation.CDS_sorted.gff
sortBed -faidx amil.Chromosome_names.txt -i amil_GFFannotation.gene.gff > amil_GFFannotation.gene_sorted.gff
sortBed -faidx amil.Chromosome_names.txt -i amil_GFFannotation.mRNA.gff > amil_GFFannotation.mRNA_sorted.gff
sortBed -faidx amil.Chromosome_names.txt -i amil_GFFannotation.three_prime_UTR.gff > amil_GFFannotation.three_prime_UTR_sorted.gff
sortBed -faidx amil.Chromosome_names.txt -i amil_GFFannotation.five_prime_UTR.gff > amil_GFFannotation.five_prime_UTR_sorted.gff
sortBed -faidx amil.Chromosome_names.txt -i amil_GFFannotation.tRNA.gff > amil_GFFannotation.tRNA_sorted.gff
sortBed -faidx amil.Chromosome_names.txt -i amil_GFFannotation.contig.gff > amil_GFFannotation.contig_sorted.gff

echo "Sorting complete!" $(date)
```

Submitted batch job 323530. Ran super fast! Now use bedtools `flank` and `subtract`. Flank will extract the 3' UTRs and subtract will remove any 3' UTR regions that overlap with known genes in the genome. Javi used 2kb and 3kb as his cutoffs (ie he extracted the flanks that were 2000 and 3000 bp around his gene of interest). In his actual analysis, he used 3kb. I'm going to use 3kb for now, but may have to redo with a different number. In the e5 scripts folder: `nano flank_sub_bed.sh`

```
#!/bin/bash 
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=5
#SBATCH --export=NONE
#SBATCH --mem=250GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/e5/scripts
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

module load BEDTools/2.30.0-GCC-11.3.0

cd /data/putnamlab/jillashey/genome/Amil_v2.01

echo "Extracting 3kb 3' UTRs" $(date)

bedtools flank -i amil_GFFannotation.gene_sorted.gff -g amil.Chromosome_lenghts.txt -l 0 -r 3000 -s | awk '{gsub("gene","3prime_UTR",$3); print $0 }' | awk '{if($5-$4 > 3)print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' | tr ' ' '\t' > amil.GFFannotation.3UTR_3kb.gff

echo "Subtract portions of UTRs that overlap nearby genes" $(date)

bedtools subtract -a amil.GFFannotation.3UTR_3kb.gff -b amil_GFFannotation.gene_sorted.gff > amil.GFFannotation.3UTR_3kb_corrected.gff 

echo "3' UTRs identified!" $(date)
```

Submitted batch job 323531. Again, ran super fast. 

```
wc -l amil.GFFannotation.3UTR_3kb_corrected.gff
48393 amil.GFFannotation.3UTR_3kb_corrected.gff
```

About 48,000 3' UTRs. Using bedtools, obtain a fasta file for the 3' UTRs. In the e5 scripts folder: `nano bed_getfasta_3UTR.sh`

```
#!/bin/bash 
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=5
#SBATCH --export=NONE
#SBATCH --mem=250GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/e5/scripts
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

module load BEDTools/2.30.0-GCC-11.3.0

cd /data/putnamlab/jillashey/genome/Amil_v2.01

echo "Extracting 3' UTR sequences" $(date)

bedtools getfasta -fi Amil.v2.01.chrs.fasta -bed amil.GFFannotation.3UTR_3kb_corrected.gff -fo amil_3UTR_3kb.fasta -fullHeader

echo "Sequence extraction complete!" $(date)
```

Submitted batch job 323533. Ran fast. Now miranda can be run! In the scripts folder: `nano miranda_strict_all_apul.sh`

```
#!/bin/bash -i
#SBATCH -t 200:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --export=NONE
#SBATCH --mem=500GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/e5/scripts
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

echo "APUL starting miranda run with all genes and miRNAs with score cutoff >100, energy cutoff <-10, and strict binding invoked"$(date)

module load Miniconda3/4.9.2
conda activate /data/putnamlab/conda/miranda 

miranda /data/putnamlab/jillashey/e5/mirna_seqs/Apul_results_mature_named.fasta /data/putnamlab/jillashey/genome/Amil_v2.01/amil_3UTR_3kb.fasta -sc 100 -en -10 -strict -out /data/putnamlab/jillashey/e5/output/miranda/miranda_strict_all_apul.tab

conda deactivate

echo "miranda run finished!"$(date)
echo "counting number of putative interactions predicted" $(date)

zgrep -c "Performing Scan" /data/putnamlab/jillashey/e5/output/miranda/miranda_strict_all_apul.tab

echo "Parsing output" $(date)
grep -A 1 "Scores for this hit:" /data/putnamlab/jillashey/e5/output/miranda/miranda_strict_all_apul.tab | sort | grep '>' > /data/putnamlab/jillashey/e5/output/miranda/miranda_strict_all_parsed_apul.txt

echo "counting number of putative interactions predicted" $(date)
wc -l /data/putnamlab/jillashey/e5/output/miranda/miranda_strict_all_parsed_apul.txt

echo "APUL miranda script complete" $(date)
```

Submitted batch job 323536. Ran in about 4 hours. Convert 3' UTR fasta to bed file. 

```
awk -F'[:-]' '/^>/ {chromosome=substr($1, 2); start=$2; end=$3; print chromosome "\t" start "\t" end}' amil_3UTR_3kb.fasta > amil_3UTR_3kb.bed
```

Run `closest` in bedtools to find the closest gene to the 3' UTRs that I created above. In the scripts folder: `nano bed_close_3UTR.sh`

```
#!/bin/bash 
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=250GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/e5/scripts
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

module load BEDTools/2.30.0-GCC-11.3.0

cd /data/putnamlab/jillashey/genome/Amil_v2.01

echo "Finding closest gene to 3'UTR seqs " $(date)

sort -k1,1 -k2,2n -o amil_3UTR_3kb_sorted.bed amil_3UTR_3kb.bed
sort -k1,1 -k4,4n amil_GFFannotation.gene_sorted.gff > sorted_amil_GFFannotation.gene_sorted.gff

bedtools closest -a amil_3UTR_3kb.bed -b amil_GFFannotation.gene_sorted.gff > closest_genes_3UTR_3kb_apul.txt

echo "Complete!" $(date)
```

Submitted batch job 323636. Got an error message but output file still produced. 

```
wc -l closest_genes_3UTR_3kb_apul.txt
73144 closest_genes_3UTR_3kb_apul.txt
```

Make ID column

```
awk 'BEGIN{FS=OFS="\t"} {split($NF, id, ";"); split(id[1], id_value, "="); $NF=id_value[2]; print $0, id_value[1]}' closest_genes_3UTR_3kb_apul.txt > modified_closest_genes_3UTR_3kb_apul.txt

wc -l modified_closest_genes_3UTR_3kb_apul.txt
73144 modified_closest_genes_3UTR_3kb_apul.txt

head modified_closest_genes_3UTR_3kb_apul.txt
chr1	48020	48410	chr1	maker	gene	15503	48020	.	+	.	Amillepora00001	ID
chr1	48020	48410	chr1	maker	gene	48411	49719	.	-	.	Amillepora00002	ID
chr1	48020	48410	chr1	maker	gene	15503	48020	.	+	.	Amillepora00001	ID
chr1	48020	48410	chr1	maker	gene	48411	49719	.	-	.	Amillepora00002	ID
chr1	49719	51020	chr1	maker	gene	48411	49719	.	-	.	Amillepora00002	ID
chr1	52175	53065	chr1	maker	gene	53066	54685	.	+	.	Amillepora00003	ID
chr1	54685	55175	chr1	maker	gene	53066	54685	.	+	.	Amillepora00003	ID
chr1	54685	55175	chr1	maker	gene	55176	69114	.	-	.	Amillepora00004	ID
chr1	54685	55175	chr1	maker	gene	53066	54685	.	+	.	Amillepora00003	ID
chr1	54685	55175	chr1	maker	gene	55176	69114	.	-	.	Amillepora00004	ID
```

Remove duplicate rows

```
awk '!seen[$0]++' modified_closest_genes_3UTR_3kb_apul.txt > uniq_modified_closest_genes_3UTR_3kb_apul.txt

wc -l uniq_modified_closest_genes_3UTR_3kb_apul.txt 
48791 uniq_modified_closest_genes_3UTR_3kb_apul.txt

head uniq_modified_closest_genes_3UTR_3kb_apul.txt 
chr1	48020	48410	chr1	maker	gene	15503	48020	.	+	.	Amillepora00001	ID
chr1	48020	48410	chr1	maker	gene	48411	49719	.	-	.	Amillepora00002	ID
chr1	49719	51020	chr1	maker	gene	48411	49719	.	-	.	Amillepora00002	ID
chr1	52175	53065	chr1	maker	gene	53066	54685	.	+	.	Amillepora00003	ID
chr1	54685	55175	chr1	maker	gene	53066	54685	.	+	.	Amillepora00003	ID
chr1	54685	55175	chr1	maker	gene	55176	69114	.	-	.	Amillepora00004	ID
chr1	73239	76239	chr1	maker	gene	69380	73239	.	+	.	Amillepora00005	ID
chr1	84565	87495	chr1	maker	gene	77282	84565	.	+	.	Amillepora00006	ID
chr1	84565	87495	chr1	maker	gene	87496	87568	.	-	.	Amillepora00007	ID
chr1	92360	95360	chr1	maker	gene	95361	108461	.	-	.	Amillepora00008	ID
```

Make new column so that the new column represents the file headers in apoc_3UTR.fasta. Ie the new column should look like this: chr1:17663-20663. 

```
awk '{print $1":"$2"-"$3, $0}' uniq_modified_closest_genes_3UTR_3kb_apul.txt > uniq_modified_closest_genes_3UTR_3kb_apul_3UTRid.txt

head uniq_modified_closest_genes_3UTR_3kb_apul_3UTRid.txt 
chr1:48020-48410 chr1	48020	48410	chr1	maker	gene	15503	48020	.	+	.	Amillepora00001	ID
chr1:48020-48410 chr1	48020	48410	chr1	maker	gene	48411	49719	.	-	.	Amillepora00002	ID
chr1:49719-51020 chr1	49719	51020	chr1	maker	gene	48411	49719	.	-	.	Amillepora00002	ID
chr1:52175-53065 chr1	52175	53065	chr1	maker	gene	53066	54685	.	+	.	Amillepora00003	ID
chr1:54685-55175 chr1	54685	55175	chr1	maker	gene	53066	54685	.	+	.	Amillepora00003	ID
chr1:54685-55175 chr1	54685	55175	chr1	maker	gene	55176	69114	.	-	.	Amillepora00004	ID
chr1:73239-76239 chr1	73239	76239	chr1	maker	gene	69380	73239	.	+	.	Amillepora00005	ID
chr1:84565-87495 chr1	84565	87495	chr1	maker	gene	77282	84565	.	+	.	Amillepora00006	ID
chr1:84565-87495 chr1	84565	87495	chr1	maker	gene	87496	87568	.	-	.	Amillepora00007	ID
chr1:92360-95360 chr1	92360	95360	chr1	maker	gene	95361	108461	.	-	.	Amillepora00008	ID
```

### Peve

To run miranda, I need to identify the 3' UTR ends. 

First, identify counts of each feature from gff file 

```
GFF_FILE="/data/putnamlab/jillashey/genome/Peve/Porites_evermanni_v1.annot.gff"
genome="/data/putnamlab/jillashey/genome/Peve/Porites_evermanni_v1.fa"

grep -v '^#' ${GFF_FILE} | cut -s -f 3 | sort | uniq -c | sort -rn > all_features_peve.txt

 231320 CDS
  40389 mRNA
  15098 UTR
```

There are a decent amount of UTRs annotated in the Peve gff but they aren't differentiated between 5' and 3'. Calculate the length of the 3' UTRs. 

```
awk '$3 == "UTR" {print $0, $5 - $4 + 1}' OFS="\t" Porites_evermanni_v1.annot.gff > UTR_lengths_Peve.txt

awk '{sum += $NF} END {if (NR > 0) print sum / NR}' UTR_lengths_Peve.txt
```

Average length of UTRs is 296.414. Extract feature types and generate individual gffs for each feature. 

```
grep $'\tCDS\t' ${GFF_FILE} > peve_GFFannotation.CDS.gff
grep $'\tmRNA\t' ${GFF_FILE} > peve_GFFannotation.mRNA.gff
grep $'\tUTR\t' ${GFF_FILE} > peve_GFFannotation.UTR.gff
```

Extract chromosome lengths 

```
cat is ${genome} | awk '$0 ~ ">" {if (NR > 1) {print c;} c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' > peve.Chromosome_lenghts.txt

wc -l peve.Chromosome_lenghts.txt 
8186 peve.Chromosome_lenghts.txt
```

Extract scaffold names 

```
awk -F" " '{print $1}' peve.Chromosome_lenghts.txt > peve.Chromosome_names.txt
```

I'm going to combine all of the bed stuff into one script. In the e5 scripts folder: `nano bed_peve.sh`

```
#!/bin/bash 
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=5
#SBATCH --export=NONE
#SBATCH --mem=250GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/e5/scripts
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

module load BEDTools/2.30.0-GCC-11.3.0

cd /data/putnamlab/jillashey/genome/Peve

echo "Sorting gffs by chromosome" $(date)

sortBed -faidx peve.Chromosome_names.txt -i peve_GFFannotation.CDS.gff > peve_GFFannotation.CDS_sorted.gff
sortBed -faidx peve.Chromosome_names.txt -i peve_GFFannotation.mRNA.gff > peve_GFFannotation.mRNA_sorted.gff
sortBed -faidx peve.Chromosome_names.txt -i peve_GFFannotation.UTR.gff > peve_GFFannotation.UTR_sorted.gff

echo "Sorting complete!" $(date)

echo "Extracting 3kb 3' UTRs" $(date)

bedtools flank -i peve_GFFannotation.mRNA_sorted.gff -g peve.Chromosome_lenghts.txt -l 0 -r 3000 -s | awk '{gsub("gene","3prime_UTR",$3); print $0 }' | awk '{if($5-$4 > 3)print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' | tr ' ' '\t' > peve.GFFannotation.3UTR_3kb.gff

echo "Subtract portions of UTRs that overlap nearby genes" $(date)

bedtools subtract -a peve.GFFannotation.3UTR_3kb.gff -b peve_GFFannotation.mRNA_sorted.gff > peve.GFFannotation.3UTR_3kb_corrected.gff 

echo "3' UTRs identified!" $(date)

echo "Extracting 3' UTR sequences" $(date)

bedtools getfasta -fi Porites_evermanni_v1.fa -bed peve.GFFannotation.3UTR_3kb_corrected.gff -fo peve_3UTR_3kb.fasta -fullHeader

echo "Sequence extraction complete!" $(date)
```

Submitted batch job 323539. Done! Now run miranda for the Peve data. In the scripts folder: `nano miranda_strict_all_peve.sh`

```
#!/bin/bash -i
#SBATCH -t 200:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --export=NONE
#SBATCH --mem=500GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/e5/scripts
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

echo "PEVE starting miranda run with all genes and miRNAs with score cutoff >100, energy cutoff <-10, and strict binding invoked"$(date)

module load Miniconda3/4.9.2
conda activate /data/putnamlab/conda/miranda 

miranda /data/putnamlab/jillashey/e5/mirna_seqs/Peve_results_mature_named.fasta /data/putnamlab/jillashey/genome/Peve/peve_3UTR_3kb.fasta -sc 100 -en -10 -strict -out /data/putnamlab/jillashey/e5/output/miranda/miranda_strict_all_peve.tab

conda deactivate

echo "miranda run finished!"$(date)
echo "counting number of putative interactions predicted" $(date)

zgrep -c "Performing Scan" /data/putnamlab/jillashey/e5/output/miranda/miranda_strict_all_peve.tab

echo "Parsing output" $(date)
grep -A 1 "Scores for this hit:" /data/putnamlab/jillashey/e5/output/miranda/miranda_strict_all_peve.tab | sort | grep '>' > /data/putnamlab/jillashey/e5/output/miranda/miranda_strict_all_parsed_peve.txt

echo "counting number of putative interactions predicted" $(date)
wc -l /data/putnamlab/jillashey/e5/output/miranda/miranda_strict_all_parsed_peve.txt

echo "PEVE miranda script complete" $(date)
```

Submitted batch job 323540


